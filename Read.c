#include "Read.h"

// GetPBC() //{{{
/*
 * Function to get box dimensions from the provided coordinate file.
 */
VECTOR GetPBC(FILE *vcf, char *input_coor) {

  VECTOR BoxLength;

  char line[LINE], line2[LINE];
  while (fgets(line, sizeof(line), vcf)) {
    strcpy(line2, line); // copy line to print in case of error

    char split[30][100], delim[8];
    strcpy(delim, " \t");
    int words = SplitLine(split, line, delim);

    if (strcmp(split[0], "pbc") == 0) {
      BoxLength.x = atof(split[1]);
      BoxLength.y = atof(split[2]);
      BoxLength.z = atof(split[3]);
      break;
    // only certain keywords besides pbc can be present before the first coordinate block
    // 1) t(imestep) ... starting the coordinate block
    // 2) t(imestep) i(ndexed)/o(ordered) ... starting the coordinate block
    // 3) i(ndexed)/o(ordered) ... starting the coordinate block
    // 4) a(tom) ... in case of vtf file
    // 5) b(ond) ... in case of vtf file
    // 6) empty line
    // 7) comment
    } else if (!(split[0][0] == 't' && words == 1) && // 1)
               !(split[0][0] == 't' && words > 1 && (split[1][0] == 'o' || split[1][0] =='i')) && // 2)
               !(split[0][0] == 'o' || split[0][0] == 'i') && // 3)
               split[0][0] != 'a' && // 4)
               split[0][0] != 'b' && // 5)
               split[0][0] != '\0' && // 6)
               split[0][0] != '#') { // 7)
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - unrecognised line\n", input_coor);
      ErrorPrintLine(split, words);
      fprintf(stderr, "\033[1;31m");
      if (IsInteger(split[0]) && words > 2) {
        fprintf(stderr, "       Possibly missing pbc line\n");
      }
      putc('\n', stderr);
      fprintf(stderr, "\033[0m");
      exit(1);
    }
  };

  return BoxLength;
} //}}}

// ReadAggCommand() //{{{
void ReadAggCommand(BEADTYPE *BeadType, COUNTS Counts,
                    char *input_coor, char *input_agg,
                    double *distance, int *contacts) {

  // open input aggregate file
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  // read first line (Aggregate command)
  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t");
  fgets(line, sizeof(line), agg);
  int words = SplitLine(split, line, delim);
  // errof if not enough striings for an Aggregate command //{{{
  if (words < 6) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError - \033[1;33m%s\033[1;31m: first line must contain a valid Aggregates command\n", input_agg);
    fprintf(stderr, "\033[0m");
    ErrorPrintLine(split, words);
    exit(1);
  } //}}}
  // read minimum distance for closeness check (<distance> argument in Aggregates utility) //{{{
  if (!IsPosDouble(split[2])) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError - \033[1;33m%s\033[1;31m: <distance> from Aggregate command must be a non-negative real number\n", input_agg);
    fprintf(stderr, "\033[0m");
    ErrorPrintLine(split, words);
    exit(1);
  }
  *distance = atof(split[2]); //}}}
  // read number of contacts aggregate check (<contacts> argument in Aggregates utility) //{{{
  if (!IsInteger(split[3])) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError - \033[1;33m%s\033[1;31m: <contacts> from Aggregate command must be a non-negative whole number\n", input_agg);
    fprintf(stderr, "\033[0m");
    ErrorPrintLine(split, words);
    exit(1);
  }
  *contacts = atof(split[3]); //}}}
  // warn if a differently named vcf file is used than the one in agg file //{{{
  if (strcmp(split[1], input_coor) != 0) {
    fprintf(stderr, "\033[1;33m");
    fprintf(stderr, "\nWARNING: the coordinate file (\033[1;36m%s\033[1;33m) ", input_coor);
    fprintf(stderr, "is different to the one in the aggregate file (\033[1;36m%s\033[1;33m).\n", split[1]);
    fprintf(stderr, "         Mismatch between beads present in both files can lead to undefined behaviour.\n\n");
    fprintf(stderr, "\033[0m");
  } //}}}
  // read <type names> from Aggregates command //{{{
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType(split[i], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m -", input_agg);
      fprintf(stderr, " - non-existent bead name \033[1;33m%s\033[1;31m in Aggregate command\n", split[i]);
      fprintf(stderr, "\033[0m");
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}
  fclose(agg);
} //}}}

// ReadStructure() //{{{
/** Function reading information about beads and molecules from DL_MESO
 * `FIELD` file, `.vsf` structure file, and `.vcf` coordinate file. Charge
 * and mass of beads is read from `FIELD` file, but all other information
 * is read from `vsf` structure file. If `vcf` coordinate is passed to
 * `ReadStructure()`, bead types not present in the `vcf` file get ignored.
 *
 * Overly complicated, some things done more than once, needs complete
 * overhaul. Thankfully, it works.
 */
bool ReadStructure(char *vsf_file, char *vcf_file, COUNTS *Counts,
                   BEADTYPE **BeadType, BEAD **Bead, int **Index,
                   MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {

  FILE *vsf;

  // zeroize stuff //{{{
  (*Counts).TypesOfBeads = 0;
  (*Counts).TypesOfMolecules = 0;
  (*Counts).Molecules = 0; //}}}

  // initial allocations (realloced later)
  *BeadType = calloc(1,sizeof(struct BeadType));
  *MoleculeType = calloc(1,sizeof(struct MoleculeType));

  // delimiters for SplitLine - whitespace in general, whitespace+colon in bond lines
  char delim[8], bond_delim[8];
  strcpy(delim, " \t");
  strcpy(bond_delim, " \t:");

  // first, read through vsf to find highest bead/mol id and all bead/mol types and identify default bead type (if exist) //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}
  int max_bead = 0, // highest bead id detected
      max_mol = 0, // highes mol id detected
      total_atom_lines = 0, // number of atom lines in vsf (including comments and blanks)
      atom_lines = 0, // number of atom lines in vsf (excluding comments and blanks)
      type_default = -1; // default bead type; -1 remains if no default exists
  char line[LINE];
  while(fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    int words = SplitLine(split, line, delim);

    if (strncmp(split[0], "atom", 1) == 0) {
      total_atom_lines++;
      atom_lines++;
      // error - odd number of values: a(tom) lines are composed of 'keyword <value>' pairs //{{{
      if ((words % 2) == 1) { // because it's split[0] ... split[2n-1]
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - odd number of strings an atom line\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // line begins with 'a(tom) default' or 'a(tom) <int>' //{{{
      bool name = false; // 'name' is necessary
      double charge = 1000, mass = 1000;
      if (strcmp("default", split[1]) == 0) {
        type_default = (*Counts).TypesOfBeads;
        // increment number of bead types
        // TODO: 'default' must be first - so why realloc here?
        (*Counts).TypesOfBeads++;
        // realloc BeadType array
        *BeadType = realloc(*BeadType, (*Counts).TypesOfBeads*sizeof(struct BeadType));
        (*BeadType)[type_default].Number = 0;
        // search the line for bead name
        for (int j = 2; j < words; j += 2) {
          if (strncmp(split[j], "name", 1) == 0) {
            // copy new name to BeadType[].Name
            strcpy((*BeadType)[type_default].Name, split[j+1]);
            (*BeadType)[type_default].Charge = 1000;
            (*BeadType)[type_default].Mass = 1000;
            name = true;
          } else if (strcmp(split[j], "charge") == 0 || strcmp(split[j], "q") == 0) {
            if (!IsDouble(split[j+1])) {
              fprintf(stderr, "\033[1;31m");
              fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom charge must be a real number\n", vsf_file);
              fprintf(stderr, "\033[0m");
              ErrorPrintLine(split, words);
              exit(1);
            }
            charge = atof(split[j+1]);
          } else if (strncmp(split[j], "mass", 1) == 0) {
            if (!IsPosDouble(split[j+1])) {
              fprintf(stderr, "\033[1;31m");
              fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom charge must be a non-negative real number\n", vsf_file);
              fprintf(stderr, "\033[0m");
              ErrorPrintLine(split, words);
              exit(1);
            }
            mass = atof(split[j+1]);
          }
        }
        if (charge != 1000) {
          (*BeadType)[type_default].Charge = charge;
        }
        if (mass != 1000) {
          (*BeadType)[type_default].Mass = mass;
        }
        if (!name) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom line must contain 'name'\n", vsf_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        }
        continue; // skip the rest of the line
      } else if (!IsInteger(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'a(tom)' must be followed by 'default' or a whole number\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } else if (atoi(split[1]) > max_bead) {
        max_bead = atoi(split[1]);
      } //}}}
      // check that line contains either both 'resid' and 'resname' or neither //{{{
      // TODO: only resid is necessary, if not resname, add some autogenerated name, creating a new MoleculeType
      int test = 0;
      for (int i = 2; i < words && test < 2; i+= 2) {
        if (strcmp("resname", split[i]) == 0) {
          test++;
        } else if (strcmp("resid", split[i]) == 0) {
          test++;
          if (!IsInteger(split[i+1]) || atoi(split[i+1]) == 0) {
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'resid' must be followed by a positive whole number\n", vsf_file);
            fprintf(stderr, "\033[0m");
            ErrorPrintLine(split, words);
            exit(1);
          }
        }
      }
      if (test == 1) { // test == 0: unbonded bead; test == 1: only one keyword present; test == 2: both keywords present
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom line must contain both 'resid' and 'resname'\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // go through the rest of the line to find bead name
      // by twos - first is always keyword, second is always value
      int type_qm = -1; // for possible mass/charge in vsf
      for (int i = 2; i < words; i += 2) {
        if (strncmp("name", split[i], 1) == 0) { // bead type name //{{{
          name = true;
          // if the bead name doesn't exist, add it to the structures
          int type = FindBeadType(split[i+1], *Counts, *BeadType);
          if (type == -1) {
            // increment number of bead types
            type = (*Counts).TypesOfBeads++;
            // realloc BeadType array
            *BeadType = realloc(*BeadType, (*Counts).TypesOfBeads*sizeof(struct BeadType));
            (*BeadType)[type].Number = 0;
            // copy new name to BeadType[].Name
            strcpy((*BeadType)[type].Name, split[i+1]);
            // initialize charge and mass by 'impossible' value
            (*BeadType)[type].Charge = 1000;
            (*BeadType)[type].Mass = 1000;
            type_qm = type;
          }
          // increment number of beads of given type
          (*BeadType)[type].Number++; //}}}
        } else if (strcmp("resname", split[i]) == 0) { // molecule type name //{{{
          // if the molecule name doesn't exist, add it to the structures
          int type = FindMoleculeType(split[i+1], *Counts, *MoleculeType);
          if (type == -1) {
            // increment number of molecule types
            type = (*Counts).TypesOfMolecules++;
            // realloc MoleculeType array
            *MoleculeType = realloc(*MoleculeType, (*Counts).TypesOfMolecules*sizeof(struct MoleculeType));
            // initialize stuff
            (*MoleculeType)[type].Number = 0;
            (*MoleculeType)[type].nBonds = 0;
            (*MoleculeType)[type].nBeads = 0;
            (*MoleculeType)[type].nBTypes = 0;
            // copy new name to MoleculeType[].Name
            strcpy((*MoleculeType)[type].Name, split[i+1]);
          } //}}}
        } else if (strcmp("resid", split[i]) == 0) { // molecule id //{{{
          if (atoi(split[i+1]) > max_mol) {
            if (!IsInteger(split[i+1]) || atoi(split[i+1]) == 0) {
              fprintf(stderr, "\033[1;31m");
              fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'resid' must be followed by a positive whole number\n", vsf_file);
              fprintf(stderr, "\033[0m");
              ErrorPrintLine(split, words);
              exit(1);
            }
            max_mol = atoi(split[i+1]);
          }
        //}}}
        } else if (strcmp("charge", split[i]) == 0 || strcmp("q", split[i]) == 0) { //{{{
          if (!IsDouble(split[i+1])) {
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom charge must be a real number\n", vsf_file);
            fprintf(stderr, "\033[0m");
            ErrorPrintLine(split, words);
            exit(1);
          }
          charge = atof(split[i+1]);
        //}}}
        } else if (strncmp("mass", split[i], 1) == 0) { //{{{
          if (!IsPosDouble(split[i+1])) {
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom mass must be a non-negative real number\n", vsf_file);
            fprintf(stderr, "\033[0m");
            ErrorPrintLine(split, words);
            exit(1);
          }
          mass = atof(split[i+1]);
        } //}}}
      }

      if (type_qm != -1 && (*BeadType)[type_qm].Charge == 1000) {
        (*BeadType)[type_qm].Charge = charge;
      }
      if (type_qm != -1 && (*BeadType)[type_qm].Mass == 1000) {
        (*BeadType)[type_qm].Mass = mass;
      }

      // error - no 'name' //{{{
      if (!name) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom line must contain 'name'\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
    } else if (split[0][0] == '#' || split[0][0] == '\0') {
      // count the line if comment or blank
      total_atom_lines++;
    } else {
      break;
    }
  }
  fclose(vsf);
  // test whether there's enough atom lines if there's no 'default' line
  if (type_default == -1 && atom_lines != (max_bead+1)) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError - \033[1;33m%s\033[1;31m: too few atom lines (or 'default' atom missing)\n", vsf_file);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  (*Counts).Molecules = max_mol; // mol ids start from 1 in vsf
  (*Counts).BeadsInVsf = max_bead + 1; // bead ids start from 0 in vsf

  // reverse of Bead[].Index
  *Index = malloc((*Counts).BeadsInVsf*sizeof(int));

  // allocate Bead and Molecule structures
  *Bead = calloc((*Counts).BeadsInVsf, sizeof(struct Bead));
  int mol_alloced = (*Counts).Molecules; // redundant as the earlier max_mol is not used later
  *Molecule = calloc(mol_alloced, sizeof(struct Molecule));

  // assign type -1 to all beads to later identify 'default' beads //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    (*Bead)[i].Type = -1;
  }
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Molecule)[i].Type = -1;
  } //}}}

  // second, read through vsf to find stuff about beads and mols //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}

  // go through atom lines - no error checking,
  // because it's been done in the first read-through
  for (int count = 0; count < total_atom_lines; count++) {
    // read line
    fgets(line, sizeof(line), vsf);
    char split[30][100];
    int words = SplitLine(split, line, delim);
    // go through the line
    if (split[0][0] == 'a' && strcmp("default", split[1]) != 0) { // non-default a(tom) line
      int bead_id = -1, bead_type = -1, mol_id = -1, mol_type = -1;
      bead_id = atoi(split[1]);
      // by twos - first is always keyword, second is always value
      for (int i = 2; i < words; i += 2) {
        if (strncmp("name", split[i], 1) == 0) { // bead name
          bead_type = FindBeadType(split[i+1], *Counts, *BeadType);
        } else if (strcmp("resname", split[i]) == 0) { // molecule name
          mol_type = FindMoleculeType(split[i+1], *Counts, *MoleculeType);
        } else if (strcmp("resid", split[i]) == 0) { // molecule id
          mol_id = atoi(split[i+1]) - 1; // mol ids start with 1 in vsf
        }
      }
      // assign values to Bead & Molecule
      (*Index)[bead_id] = bead_id; // reverse of Bead[].Index
      (*Bead)[bead_id].Type = bead_type;
      (*Bead)[bead_id].Index = bead_id;
      (*Bead)[bead_id].Molecule = mol_id;
      if (mol_id > -1) {
        // TODO: now - mtype according to last molecule of the type; should be - according to first molecule; i.e., add some if stuff
        (*Molecule)[mol_id].Type = mol_type;
        (*MoleculeType)[mol_type].nBeads++;
      }
    }
  } //}}}

  // assign 'type_default' to default beads //{{{
  if (type_default != -1) {
    int count = 0;
    for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
      if ((*Bead)[i].Type == -1) {
        (*Index)[i] = i; // reverse of Bead[].Index
        (*Bead)[i].Type = type_default;
        (*Bead)[i].Index = i;
        (*Bead)[i].Molecule = -1; // default beads aren't in molecules
        count++;
      } else if ((*Bead)[i].Type == type_default) { // default type beads explicitly specified by 'atom' line
        count++;
      }
    }
    (*BeadType)[type_default].Number = count;
  } //}}}

  // calculate number of molecules/beads in each molecule type //{{{
  int *mols = calloc((*Counts).TypesOfMolecules, sizeof(int)); // number of molecules for each molecule type
  for (int i = 0; i < (*Counts).Molecules; i++) {
    mols[(*Molecule)[i].Type]++;
  }
  // assign info to MoleculeType
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Number = mols[i];
    (*MoleculeType)[i].nBeads /= (*MoleculeType)[i].Number;
  }
  free(mols); //}}}

  // third, go through the bonds section of vsf to find the number of bonds in molecules //{{{
  int bonds[(*Counts).TypesOfMolecules]; // helper array //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  } //}}}
  while (fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    int words = SplitLine(split, line, bond_delim);

    if (strncmp(split[0], "bond", 1) == 0) {
      if (!IsInteger(split[1]) || !IsInteger(split[2])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bond' must be followed by two colon-separated numbers\n", vsf_file);
        fprintf(stderr, "\033[0m");
        exit(1);
      }
      int mol_id = (*Bead)[atoi(split[1])].Molecule;
      int mol_type = (*Molecule)[mol_id].Type;
      // is this mol_type already in use in bond numbers calculation?
      if (bonds[mol_type] == -1) {
        bonds[mol_type] = mol_id;
      }
      // increment number of bonds if correct molecule
      if (bonds[mol_type] == mol_id) {
        (*MoleculeType)[mol_type].nBonds++;
      }
    /* Break while loop if at end of bond section
     * Excepting comments or white lines, behind bonds section can be either coordinates (vtf file):
     * timestep starting line can be
     * 1) t(imestep) ...and nothing behind it
     * 2) t(imestep) o(rdered)/i(ndexed)
     * 3) o(rdered)/i(ndexed)
     * or pbc line:
     * 4) pbc <number> <number> <number>
     */
    } else if ((strncmp(split[0], "timestep", 1) == 0 && words == 1) || // 1)
               (strncmp(split[0], "timestep", 1) == 0 && words > 1 && // 2) 1st string: t(imestep)
                (strncmp(split[1], "ordered", 1) == 0 || strncmp(split[1], "indexed", 1) == 0)) || // 2) 2nd string: o/i
               (strncmp(split[0], "ordered", 1) == 0 || strncmp(split[0], "indexed", 1) == 0) || // 3)
               (strcmp("pbc", split[0]) == 0 && words > 3 && // 4)
                IsPosDouble(split[1]) && // 4) 1st <number>
                IsPosDouble(split[2]) && // 4) 2nd <number>
                IsPosDouble(split[2]))) { // 4) 3rd <number>
      break;
    // exit with error if not empty line or comment
    } else if (split[0][0] != '\0' && // empty line
               split[0][0] != '#') { // comment
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError - \033[1;33m%s\033[1;31m: unrecognised line\n", vsf_file);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
  } //}}}

  fclose(vsf); // will be later reopened to read from the beginning again

  // find the highest number of bonds //{{{
  int max_bonds = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    if ((*MoleculeType)[i].nBonds > max_bonds) {
      max_bonds = (*MoleculeType)[i].nBonds;
    }
  } //}}}

  // allocate MoleculeType[].Bond array //{{{
  int moltype_alloced = (*Counts).TypesOfMolecules;
  for (int i = 0; i < moltype_alloced; i++) {
    (*MoleculeType)[i].Bond = calloc(max_bonds, sizeof(int *));
    for (int j = 0; j < max_bonds; j++) {
      (*MoleculeType)[i].Bond[j] = calloc(3, sizeof(int));
      (*MoleculeType)[i].Bond[j][2] = -1; // no bond type assigned
    }
  } //}}}

  // fourth, go through vsf and assign bead ids to molecules //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}

  // go through atom lines - no error checking,
  // because it's been done in the first read-through
  int *beads = calloc((*Counts).BeadsInVsf, sizeof(int)); // helper array
  for (int count = 0; count < total_atom_lines; count++) {
    fgets(line, sizeof(line), vsf);
    char split[30][100];
    int words = SplitLine(split, line, delim);

    // go through the line
    if (strncmp(split[0], "atom", 1) == 0 && strcmp("default", split[1]) != 0) { // non-default a(tom) line
      int bead_id = atoi(split[1]);
      // by twos - first is always keyword, second is always value
      int mol_id = -1, mtype = -1;
      for (int i = 2; i < words; i += 2) {
        if (strcmp("resid", split[i]) == 0) { // molecule id
          mol_id = atoi(split[i+1]) - 1; // mol ids start with 1 in vsf
        } else if (strcmp("resname", split[i]) == 0) { // molecule name (i.e., type)
          mtype = FindMoleculeType(split[i+1], *Counts, *MoleculeType);
        }
      }
      if (mol_id > -1) {
        // allocate Molecule[].Bead array if bead_id is the first bead of the molecule
        if (beads[mol_id] == 0) {
          (*Molecule)[mol_id].Bead = calloc((*MoleculeType)[mtype].nBeads, sizeof(int));
        }
        (*Molecule)[mol_id].Bead[beads[mol_id]] = bead_id;
        beads[mol_id]++;
      }
    }
  }
  free(beads);

  // sort ids in molecule according to ascending id in vsf //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    SortArray(&(*Molecule)[i].Bead, (*MoleculeType)[mtype].nBeads, 0);
  } //}}}
  //}}}

  // fill MoleculeType[x].Bead array according to the first molecule with the name MoleculeType[x].Name //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Bead = calloc((*MoleculeType)[i].nBeads, sizeof(int));
    (*MoleculeType)[i].Bead[0] = -1; // just to test whether this molecule type was used already
  }
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    if ((*MoleculeType)[mtype].Bead[0] == -1) {
      for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
        int bead = (*Molecule)[i].Bead[j];
        int btype = (*Bead)[bead].Type;
        (*MoleculeType)[mtype].Bead[j] = btype;
      }
    }
  } //}}}

  // fifth, go again through the bonds and assign ids of bonded beads //{{{
  // initialize helper arrays //{{{
  int *count_bonds;
  count_bonds = calloc((*Counts).TypesOfMolecules, sizeof(int));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  } //}}}
  // again, no error checking, as it was done earlier already
  while (fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    SplitLine(split, line, bond_delim);
    if (strncmp(split[0], "bond", 1) == 0) {
      int mol_id = (*Bead)[atoi(split[1])].Molecule;
      int mol_type = (*Molecule)[mol_id].Type;
      // is this mol_type already in use?
      if (bonds[mol_type] == -1) {
        bonds[mol_type] = mol_id;
      }
      // increment number of bonds if correct molecule
      if (bonds[mol_type] == mol_id) {
        int id1 = atoi(split[1]);
        int id2 = atoi(split[2]);
        bool test_id1 = false, test_id2 = false;
        for (int i = 0; i < (*MoleculeType)[mol_type].nBeads; i++) {
          if (!test_id1 && id1 == (*Molecule)[mol_id].Bead[i]) {
            (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][0] = i;
            test_id1 = true;
          }
          if (!test_id2 && id2 == (*Molecule)[mol_id].Bead[i]) {
            (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][1] = i;
            test_id2 = true;
          }
          if (test_id1 && test_id2) {
            break;
          }
        }
        count_bonds[mol_type]++;
      }
    } else if (split[0][0] != '\0' && split[0][0] != '#') {
      break;
    }
  }
  free(count_bonds);
  fclose(vsf); //}}}

  // sort the Bond arrays //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    SortBonds((*MoleculeType)[i].Bond, (*MoleculeType)[i].nBonds);
  } //}}}

  // fill BTypes array //{{{
  // allocate & initialize arrays //{{{
  bool used[(*Counts).TypesOfMolecules];
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].BType = calloc((*Counts).TypesOfBeads, sizeof(int));
    (*MoleculeType)[i].nBTypes = 0;
    used[i] = false;
  } //}}}
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mol_type = (*Molecule)[i].Type;
    if (!used[mol_type]) { // go only through the first molecule of the given type
      for (int j = 0; j < (*MoleculeType)[mol_type].nBeads; j++) {
        int id = (*Molecule)[i].Bead[j];
        // test if bead id's type is present in BType array //{{{
        bool in_mol = false;
        for (int k = 0; k < (*MoleculeType)[mol_type].nBTypes; k++) {
          if ((*MoleculeType)[mol_type].BType[k] == (*Bead)[id].Type) {
            in_mol = true;
            break;
          }
        } //}}}
        // if bead id is of a type not yet present in BType array, add it there
        if (!in_mol) {
          (*MoleculeType)[mol_type].BType[(*MoleculeType)[mol_type].nBTypes] = (*Bead)[id].Type;
          (*MoleculeType)[mol_type].nBTypes++;
        }
      }
    }
  } //}}}

  // sixth, go through vcf to find which bead types are there //{{{
  bool indexed = true;
  if (vcf_file[0] != '\0') {
    // open vcf coordinate file //{{{
    FILE *vcf;
    if ((vcf = fopen(vcf_file, "r")) == NULL) {
      ErrorFileOpen(vcf_file, 'r');
      exit(1);
    } //}}}

    // skip initial stuff //{{{
    // skip lines including 't(imestep)' line //{{{
    char line2[LINE], *split[30], str[30];
    str[0] = '\0';
    do {
      fgets(line, sizeof(line), vcf);
      strcpy(line2, line);
      split[0] = strtok(line, " \t");

      // 't(imestep)' line
      if (split[0][0] == 't' || split[0][0] == 'T') {
        split[1] = strtok(NULL, " \t");
        // if only 't(imestep)' present, assume ordered timestep
        if (split[1] == NULL) {
          str[0] = 'o';
          break;
        // otherwise, error if 'i(ndexed)' or 'o(ordered)' not present
        } else if (split[1][0] != 'o' && split[1][0] != 'O' &&
            split[1][0] != 'i' && split[1][0] != 'I') {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", vsf_file);
          fprintf(stderr, " - unrecognised keywords '\033[1;33m%s %s\033[1;31m'", split[0], split[1]);
          fprintf(stderr, "\033[0m");
          exit(1);
        // if no error, find if 'i(ndexed)' or 'o(rdered)' timestep
        } else {
          str[0] = split[1][0];
          break;
        }

      // 'o(rdered)' or 'i(ndexed)' line
      } else if (split[0][0] == 'o' || split[0][0] == 'O' ||
                 split[0][0] == 'i' || split[0][0] == 'I') {
        str[0] = split[0][0];
        break;
      }
    } while (true); //}}}

    // skip the rest until first coordinates //{{{
    fpos_t pos;
    do {
      fgetpos(vcf, &pos); // save vcf file pointer

      fgets(line, sizeof(line), vcf);
      split[0] = strtok(line, " \t");
    } while (split[0][0] < '0' || split[0][0] > '9');
    fsetpos(vcf, &pos); // restore vcf file pointer to before the first coordinate line //}}}
    //}}}

    if (str[0] == 'o' || str[0] == 'O') { // ordered timesteps //{{{
      indexed = false;

      // set all bead types to 'use'
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        (*BeadType)[i].Use = true;
      }
      (*Counts).Beads = (*Counts).BeadsInVsf;
      //}}}
    } else if (str[0] == 'i' || str[0] == 'I') { // indexed timesteps //{{{
      indexed = true;
      // set all bead types to 'do not use' //{{{
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        (*BeadType)[i].Use = false;
      } //}}}

      // read data
      // the first coordinate line
      fgets(line, sizeof(line), vcf);
      strcpy(line, TrimLine(line)); // trim excess whitespace
      // for lines containing only whitespace
      if (strlen(line) == 1) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - blank line instead of a first coordinate line\n\n", vcf_file);
        fprintf(stderr, "\033[0m");
        exit(1);
      }

      // split the line into array //{{{
      split[0] = strtok(line, " \t");
      int i = 0;
      while (i < 29 && split[i] != NULL) {
        i++;
        split[i] = strtok(NULL, " \t");
      } //}}}

      // first split is bead index; read data till there is <int> at the beginning of the line //{{{
      while (split[0][0] >= '0' && split[0][0] <= '9') {
        int type = (*Bead)[atoi(split[0])].Type;
        (*BeadType)[type].Use = true;
        char test;
        if ((test = getc(vcf)) == EOF) {
          break;
        }
        ungetc(test, vcf);
        // read line
        fgets(line, sizeof(line), vcf);
        strcpy(line, TrimLine(line)); // trim excess whitespace
        // split the line into array //{{{
        split[0] = strtok(line, " \t");
        i = 0;
        while (i < 29 && split[i] != NULL) {
          i++;
          split[i] = strtok(NULL, " \t");
        } //}}}
      } //}}}

      // count the number of beads in vcf //{{{
      (*Counts).Beads = 0;
      for (i = 0; i < (*Counts).TypesOfBeads; i++) {
        if ((*BeadType)[i].Use) {
          (*Counts).Beads += (*BeadType)[i].Number;
        }
      } //}}}
    //}}}
    } else { // error //{{{
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'i(ndexed)' or 'o(rdered)' keyword\n\n", vcf_file);
      fprintf(stderr, "\033[0m");
      exit(1);
    } //}}}
    fclose(vcf);

  // use all beads if no vcf file provided
  } else {
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      (*BeadType)[i].Use = true;
    }
    (*Counts).Beads = (*Counts).BeadsInVsf;
  } //}}}

  // determine which molecule types are present in vcf //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Use = false;
  }
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // find what bead types are in vcf
    int test = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int type = (*MoleculeType)[i].BType[j];
      if ((*BeadType)[type].Use) {
        test++;
      }
    }
    if (test != 0) {
      (*MoleculeType)[i].Use = true;
    }
  } //}}}

  // remove unused molecules from Molecule struct  //{{{
  int mols_in_vcf = 0;
  int count = 0;
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int type = (*Molecule)[i].Type;
    if ((*MoleculeType)[type].Use) {
      (*Molecule)[count].Type = (*Molecule)[i].Type;
      for (int j = 0; j < (*MoleculeType)[type].nBeads; j++) {
        (*Molecule)[count].Bead[j] = (*Molecule)[i].Bead[j];
      }
      mols_in_vcf++;
      count++;
    }
  }
  (*Counts).Molecules = mols_in_vcf; //}}}

  // copy molecule type names //{{{
  char **name;
  name = calloc(moltype_alloced, sizeof(char *));
  for (int i = 0; i < moltype_alloced; i++) {
    name[i] = calloc(16, sizeof(char));
  }
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    strcpy(name[i], (*MoleculeType)[i].Name);
  } //}}}

  // remove unused molecule types //{{{
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    if ((*MoleculeType)[i].Use) {
      if (i != count) {
        strcpy((*MoleculeType)[count].Name, (*MoleculeType)[i].Name);
        (*MoleculeType)[count].Number = (*MoleculeType)[i].Number;
        (*MoleculeType)[count].nBeads = (*MoleculeType)[i].nBeads;
        for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
          (*MoleculeType)[count].Bead[j] = (*MoleculeType)[i].Bead[j];
        }
        (*MoleculeType)[count].nBonds = (*MoleculeType)[i].nBonds;
        for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
          (*MoleculeType)[count].Bond[j][0] = (*MoleculeType)[i].Bond[j][0];
          (*MoleculeType)[count].Bond[j][1] = (*MoleculeType)[i].Bond[j][1];
        }
        (*MoleculeType)[count].nBTypes = (*MoleculeType)[i].nBTypes;
        for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
          (*MoleculeType)[count].BType[j] = (*MoleculeType)[i].BType[j];
        }
      }
      count++;
    }
  }
  (*Counts).TypesOfMolecules = count; //}}}

  // correct molecule ids (in case of some molecule type absent from vcf) //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int type = (*Molecule)[i].Type;
    (*Molecule)[i].Type = FindMoleculeType(name[type], *Counts, *MoleculeType);
  } //}}}

  // remove unused beads from bonds in molecules //{{{
  // helper array indicating if a given type was already changed //{{{
  int mol_type[(*Counts).TypesOfMolecules];
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    mol_type[i] = -1;
  } //}}}
  // go through all molecules
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    // is this molecule type done?
    if (mol_type[mtype] == -1) {
      count = 0;
      // go through all bonds and also find initial minimal id //{{{
      int min1 = 1000000;
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        if ((*MoleculeType)[mtype].Bond[j][0] < min1) {
          min1 = (*MoleculeType)[mtype].Bond[j][0];
        }
        int id1 = (*Molecule)[i].Bead[(*MoleculeType)[mtype].Bond[j][0]];
        int btype1 = (*Bead)[id1].Type;
        int id2 = (*Molecule)[i].Bead[(*MoleculeType)[mtype].Bond[j][1]];
        int btype2 = (*Bead)[id2].Type;
        if ((*BeadType)[btype1].Use && (*BeadType)[btype2].Use) {
          (*MoleculeType)[mtype].Bond[count][0] = (*MoleculeType)[mtype].Bond[j][0];
          (*MoleculeType)[mtype].Bond[count][1] = (*MoleculeType)[mtype].Bond[j][1];
          count++;
        }
      } //}}}
      // change number of bonds
      (*MoleculeType)[mtype].nBonds = count;
      // find minimum id in Bond array //{{{
      int min2 = 1000000;
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        if ((*MoleculeType)[mtype].Bond[j][0] < min2) {
          min2 = (*MoleculeType)[mtype].Bond[j][0];
        }
      } //}}}
      // minimize ids in Bond array //{{{
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        (*MoleculeType)[mtype].Bond[j][0] += min1 - min2;
        (*MoleculeType)[mtype].Bond[j][1] += min1 - min2;
      } //}}}
      // go through all beads to correct MoleculeType[].Bead array //{{{
      count = 0;
      for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
        int btype = (*MoleculeType)[mtype].Bead[j];
        if ((*BeadType)[btype].Use) {
          (*MoleculeType)[mtype].Bead[count] = (*MoleculeType)[mtype].Bead[j];
          count++;
        }
      } //}}}
      // mark this molecule type as done
      mol_type[mtype] = i;
    }
  } //}}}

  // remove unused bead types in molecules (Molecule[].Bead arrays) //{{{
  int beads_in_mols[(*Counts).TypesOfMolecules];
  for (int i = 0; i < (*Counts).Molecules; i++) {
    count = 0;
    int mtype = (*Molecule)[i].Type;
    for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
      int type = (*Bead)[(*Molecule)[i].Bead[j]].Type;
      if ((*BeadType)[type].Use) {
        (*Molecule)[i].Bead[count] = (*Molecule)[i].Bead[j];
        (*Bead)[(*Molecule)[i].Bead[j]].Molecule = i;
        count++;
      }
    }
    beads_in_mols[mtype] = count;
  }
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nBeads = beads_in_mols[i];
  } //}}}

  // remove unused bead types in molecules (MoleculeType[].BType arrays) //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    count = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int type = (*MoleculeType)[i].BType[j];
      if ((*BeadType)[type].Use) {
        (*MoleculeType)[i].BType[count] = type;
        count++;
      }
    }
    (*MoleculeType)[i].nBTypes = count;
  } //}}}

  // remove unused beads from Bead struct //{{{
  count = 0;
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int type = (*Bead)[i].Type;
    if ((*BeadType)[type].Use) {
      (*Bead)[count].Type = (*Bead)[i].Type;
      (*Bead)[count].Molecule = (*Bead)[i].Molecule;
      (*Bead)[count].Index = (*Bead)[i].Index;
      (*Index)[(*Bead)[i].Index] = count;
      int id = (*Bead)[count].Molecule;
      if (id != -1) {
        type = (*Molecule)[id].Type;
        for (int j = 0; j < (*MoleculeType)[type].nBeads; j++) {
          if ((*Molecule)[id].Bead[j] == (*Bead)[count].Index) {
            (*Molecule)[id].Bead[j] = count;
          }
        }
      }
      count++;
    }
  } //}}}

  // copy bead type names //{{{
  // first, free name[] ...to prevent memory leak
  for (int i = 0; i < moltype_alloced; i++) {
    free(name[i]);
  }
  // second, realloc name
  int beadtype_alloced = (*Counts).TypesOfMolecules;
  name = realloc(name, beadtype_alloced*sizeof(char *));
  // third, calloc name[]
  for (int i = 0; i < beadtype_alloced; i++) {
    name[i] = calloc(16, sizeof(char));
  }
  // fourth, copy names
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    strcpy(name[i], (*MoleculeType)[i].Name);
  } //}}}

  // remove unused bead types from BeadType struct and molecule bonds //{{{
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if ((*BeadType)[i].Use) {
      if (count != i) {
        strcpy((*BeadType)[count].Name, (*BeadType)[i].Name);
        (*BeadType)[count].Number = (*BeadType)[i].Number;
        (*BeadType)[count].Charge = (*BeadType)[i].Charge;
        (*BeadType)[count].Mass = (*BeadType)[i].Mass;
        (*BeadType)[count].Use = (*BeadType)[i].Use;
        for (int j = 0; j < (*Counts).TypesOfMolecules; j++) {
          for (int k = 0; k < (*MoleculeType)[j].nBTypes; k++) {
            if ((*MoleculeType)[j].BType[k] == i) {
              (*MoleculeType)[j].BType[k] = count;
            }
          }
          for (int k = 0; k < (*MoleculeType)[j].nBeads; k++) {
            if ((*MoleculeType)[j].Bead[k] == i) {
              (*MoleculeType)[j].Bead[k] = count;
            }
          }
        }
        for (int j = 0; j < (*Counts).Beads; j++) {
          if ((*Bead)[j].Type == i) {
            (*Bead)[j].Type = count;
          }
        }
      }
      count++;
    }
  } //}}}

  (*Counts).TypesOfBeads = count;

  // free name array //{{{
  for (int i = 0; i < beadtype_alloced; i++) {
    free(name[i]);
  }
  free(name); //}}}

  // free extraneous arrays to prevent memory leak //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = (*MoleculeType)[i].nBonds; j < max_bonds; j++) {
      free((*MoleculeType)[i].Bond[j]);
    }
  }
  for (int i = (*Counts).TypesOfMolecules; i < moltype_alloced; i++) {
    for (int j = 0; j < max_bonds; j++) {
      free((*MoleculeType)[i].Bond[j]);
    }
    free((*MoleculeType)[i].Bond);
    free((*MoleculeType)[i].BType);
  }
  for (int i = (*Counts).Molecules; i < mol_alloced; i++) {
    free((*Molecule)[i].Bead);
  } //}}}

  // assign default charge (0) and mass (1) to bead types if not read from vsf //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if ((*BeadType)[i].Charge == 1000) {
      (*BeadType)[i].Charge = 0;
    }
    if ((*BeadType)[i].Mass == 1000) {
      (*BeadType)[i].Mass = 1;
    }
  } //}}}

  // calculate molecule mass //{{{
  count = 0; // to store id of the first molecule id of the given type
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int id = (*Molecule)[count].Bead[j];
      int type = (*Bead)[id].Type;
      (*MoleculeType)[i].Mass += (*BeadType)[type].Mass;
    }
    count += (*MoleculeType)[i].Number;
  } //}}}

  // calculate number of (un)bonded beads //{{{
  (*Counts).Unbonded = 0;
  (*Counts).Bonded = 0;
  for (int i = 0; i < (*Counts).Beads; i++) {
    if ((*Bead)[i].Molecule == -1) {
      (*Counts).Unbonded++;
    } else {
      (*Counts).Bonded++;
    }
  } //}}}

  // allocate Bead[].Aggregate array //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(1, sizeof(int));
    (*Bead)[i].Aggregate[0] = -1;
  } //}}}

  // fill Molecule[].Aggregate //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Molecule)[i].Aggregate = -1;
  } //}}}

  // set all molecule & bead types to be unused //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Use = false;
    (*MoleculeType)[i].Write = false;
  }
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i].Use = false;
    (*BeadType)[i].Write = false;
  } //}}}

  // allocate memory for MoleculeType[].Angles (just to free later) //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nAngles = 0;
    (*MoleculeType)[i].Angle = calloc(1, sizeof(int *));
  } //}}}

  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, vsf_file);

  return indexed;
} //}}}

// ReadTimestepPreamble() //{{{
int ReadTimestepPreamble(bool indexed, char *input_coor, FILE *vcf_file, char **stuff) {
  // save pointer position in vcf_file
  fpos_t position;
  fgetpos(vcf_file, &position);

  (*stuff)[0] = '\0'; // empty the array
  int words, count_lines = 0;
  char split[30][100];
  bool timestep = false;
  do {
    char line[LINE];
    fgets(line, sizeof(line), vcf_file);
    if (feof(vcf_file)) {
      return -1;
    }
    words = SplitLine(split, line, " \t");
    // comment line - copy to stuff array
    if (split[0][0] == '#') {
      char temp[LINE];
      strcpy(temp, *stuff);
      sprintf(temp, "%s\n", *stuff);
      strcpy(*stuff, temp);
      for (int i = 0; i < words; i++) {
        sprintf(temp, "%s %s", *stuff, split[i]);
        strcpy(*stuff, temp);
      }
    // error if not timestep or pbc line, blank line, or a double (i.e., start of the coordinate lines)
    } else if (split[0][0] != 't' && split[0][0] != 'i' && split[0][0] != 'o' &&
               words != 0 && strcasecmp(split[0], "pbc") != 0 && !IsDouble(split[0])) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
      fprintf(stderr, " - unrecognised line in the timestep preamble\n");
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    // test for a t(imestep) i(ndexed)/o(rdered) line - error if different than the first timestep //{{{
    if ((words > 1 && split[0][0] == 't' && (split[1][0] == 'o' || split[1][0] == 'i')) ||
        (split[0][0] == 'o' || split[0][0] == 'i')) {
      timestep = true;
      if ((split[0][0] == 'o' || split[1][0] == 'o') && indexed) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
        fprintf(stderr, " - ordered timestep instead of an indexed one\n");
        fprintf(stderr, "\033[0m");
        exit(1);
      } else if ((split[0][0] == 'i' || split[1][0] == 'i') && !indexed) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
        fprintf(stderr, " - indexed timestep instead of an ordered one\n");
        fprintf(stderr, "\033[0m");
        exit(1);
      }
    } //}}}
    count_lines++;
  } while (words == 0 || !IsDouble(split[0]));
  count_lines--; // the last counted line contained the first coordinate line
  // error - missing timestep line //{{{
  if (!timestep) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
    fprintf(stderr, " - missing [t(imestep)] o(rdered)/i(indexed) line\n");
    if (indexed) {
      fprintf(stderr, "       possibly more coordinate lines than in the first timestep\n");
    }
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}
  fsetpos(vcf_file, &position); // restore pointer position
  return count_lines;
} //}}}

// ReadCoordinates() //{{{
/**
 * Function reading coordinates from .vcf file with indexed timesteps (\ref IndexedCoorFile).
 */
void ReadCoordinates(bool indexed, char *input_coor, FILE *vcf_file, COUNTS Counts, int *Index, BEAD **Bead, char **stuff) {
  int count_lines = ReadTimestepPreamble(indexed, input_coor, vcf_file, stuff);
  // skip the preamble lines //{{{
  for (int i = 0; i < count_lines; i++) {
    char line[LINE];
    fgets(line, sizeof(line), vcf_file);
  } //}}}
  if (indexed) { // indexed timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[30][100];
      fgets(line, sizeof(line), vcf_file);
      // error - end of file
      if (feof(vcf_file)) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
        fprintf(stderr, " - unexpected end of file (possibly fewer beads in a timestep than in the first timestep)\n");
        fprintf(stderr, "\033[0m");
        exit(1);
      }
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <int> <double> <double> <double>
      if (words < 4 || !IsInteger(split[0]) || !IsDouble(split[1]) ||
          !IsDouble(split[2]) || !IsDouble(split[3])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
        fprintf(stderr, " - cannot read coordinates\n");
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      int index = atoi(split[0]);
      // bead coordinates
      (*Bead)[Index[index]].Position.x = atof(split[1]);
      (*Bead)[Index[index]].Position.y = atof(split[2]);
      (*Bead)[Index[index]].Position.z = atof(split[3]);
    } //}}}
  } else { // ordered timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[30][100];
      fgets(line, sizeof(line), vcf_file);
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <double> <double> <double>
      if (words < 3 || !IsDouble(split[0]) ||
          !IsDouble(split[1]) || !IsDouble(split[2])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
        fprintf(stderr, " - cannot read coordinates\n");
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      // bead coordinates
      (*Bead)[i].Position.x = atof(split[0]);
      (*Bead)[i].Position.y = atof(split[1]);
      (*Bead)[i].Position.z = atof(split[2]);
    }
  } //}}}
} //}}}

// SkipCoor() //{{{
/**
 * Function to skip one timestep in coordinates file. It works with both
 * indexed and ordered vcf files.
 */
bool SkipCoor(FILE *vcf_file, COUNTS Counts, char **stuff) {

  bool error = false;
  // initial stuff //{{{
  (*stuff)[0] = '\0'; // no comment line
  char line[LINE];
  fpos_t position;
  do {
    fgetpos(vcf_file, &position); // save pointer position
    fgets(line, sizeof(line), vcf_file);

    while (line[0] == '\t' || line[0] == ' ') {
      int i;
      for (i = 1; line[i] != '\n'; i++) {
        line[i-1] = line[i];
      }
      line[i-1] = '\n';
      line[i] = '\0';
    }

    if (line[0] == '#') {
      strcat(*stuff, line);
    }
  } while (line[0] < '0' || line[0] > '9');

  // return file pointer to before the first coordinate line
  fsetpos(vcf_file, &position); //}}}
  for (int i = 0; i < Counts.Beads; i++) {
    int test;
    while ((test = getc(vcf_file)) != '\n' && test != EOF)
      ;
    // premature end of file
    if (test == EOF) {
      error = true;
      break;
    }
  }
  if (!error) {
    getc(vcf_file);
  }
  return error;
} //}}}

// ReadAggregates() //{{{
/**
 * Function reading information about aggregates from `.agg` file (\ref AggregateFile) generated by Aggregates utility.
 */
bool ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts, AGGREGATE **Aggregate,
                    BEADTYPE *BeadType, BEAD **Bead,
                    MOLECULETYPE *MoleculeType, MOLECULE **Molecule, int *Index) {

  bool error = false;

  char line[LINE], split[30][100], delim[8];

  // read (Last) Step line
  fgets(line, sizeof(line), fr);
  strcpy(delim, " \t"); // delimiters for SplitLine()
  int words = SplitLine(split, line, delim);

  // is there a Step? i.e., isn't this the line 'Last Step'?
  if (split[0][0] == 'L') {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file", agg_file);
    fprintf(stderr, " (%s %d)\n\n", split[0], atoi(split[1]));
    fprintf(stderr, "\033[0m");
    exit(1);
  // Step line must be 'Step: <int>'
  } else if (words < 2 || !IsInteger(split[1])) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - wrong 'Step' line\n", agg_file);
    ErrorPrintLine(split, words);
    fprintf(stderr, "\033[0m");
    exit(1);
  } else {
    // initialize array of number of aggregates per bead //{{{
    for (int i = 0; i < (*Counts).Beads; i++) {
      (*Bead)[i].nAggregates = 0;
    } //}}}

    // get number of aggregates //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // number of aggregates must be <int>
    if (words != 0 && !IsInteger(split[0])) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - number of aggregates must be a whole number\n", agg_file);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*Counts).Aggregates = atoi(split[0]); //}}}

    // skip blank line
    fgets(line, sizeof(line), fr);

    // go through all aggregates
    for (int i = 0; i < (*Counts).Aggregates; i++) {
      // read molecules in Aggregate 'i' //{{{
      fscanf(fr, "%d :", &(*Aggregate)[i].nMolecules);
      for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
        int mol;
        fscanf(fr, "%d", &mol);
        mol--; // in agg file the numbers correspond to vmd

        (*Aggregate)[i].Molecule[j] = mol;
        (*Molecule)[mol].Aggregate = i;
      }

      while (getc(fr) != '\n')
       ; //}}}

      // read monomeric beads in Aggregate 'i' //{{{
      int count;
      fscanf(fr, "%d :", &count);
      (*Aggregate)[i].nMonomers = 0; // their number will be counted according to which beads are in vcf
      for (int j = 0; j < count; j++) {
        int vsf_id;
        fscanf(fr, "%d", &vsf_id); // monomer index in vsf file
        int id = Index[vsf_id]; // monomer index in Bead structure
        if (id > -1) {
          int beads = (*Aggregate)[i].nMonomers++;
          (*Aggregate)[i].Monomer[beads] = id;
          (*Bead)[id].nAggregates++;
          (*Bead)[id].Aggregate = realloc((*Bead)[id].Aggregate, (*Bead)[id].nAggregates*sizeof(int));
          (*Bead)[id].Aggregate[(*Bead)[id].nAggregates-1] = i;
        }
      }

      while (getc(fr) != '\n')
       ; //}}}
    }

    // skip blank line at the end of every entry
    while (getc(fr) != '\n')
      ;

    // fill Aggregate[].Bead, Bead[].Aggregate, and Molecule[].Aggregate arrays //{{{
    for (int i = 0; i < (*Counts).Aggregates; i++) {
      (*Aggregate)[i].nBeads = 0;

      for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
        int mol = (*Aggregate)[i].Molecule[j];
        (*Molecule)[mol].Aggregate = i;
        for (int k = 0; k < MoleculeType[(*Molecule)[mol].Type].nBeads; k++) {
          int bead = (*Molecule)[mol].Bead[k];
          (*Aggregate)[i].Bead[(*Aggregate)[i].nBeads+k] = bead;
          (*Bead)[bead].Aggregate[(*Bead)[bead].nAggregates] = i;
          (*Bead)[bead].nAggregates++;
        }

        // increment number of beads in aggregate 'i'
        (*Aggregate)[i].nBeads += MoleculeType[(*Molecule)[mol].Type].nBeads;
      }
    }

    // fill Bead[].Aggregate for monomeric Beads
//  for (int i = 0, i < (*Counts).Bead)
    //}}}

    // calculate aggregates' masses //{{{
    for (int i = 0; i < (*Counts).Aggregates; i++) {
      (*Aggregate)[i].Mass = 0;

      for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
        int type = (*Molecule)[(*Aggregate)[i].Molecule[j]].Type;

        (*Aggregate)[i].Mass += MoleculeType[type].Mass;
      }
    } //}}}
  }

  return error;
} //}}}

// ReadFieldPbc() //{{{
/**
 * Function reading box size from the first line of the FIELD-like file.
 */
void ReadFieldPbc(char *field, VECTOR *BoxLength) {
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}

  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t"); // delimiters for SplitLine()

  // read box size //{{{
  fgets(line, sizeof(line), fr);
  int words = SplitLine(split, line, delim);
  // first line of the FIELD must be: <double> <double> <double> ...whatever
  // Error if:
  // 1) too few strings
  // 2) enough strings, but not numbers
  if (words < 3 || // 1
      !IsPosDouble(split[0]) || // 2) first <double>
      !IsPosDouble(split[1]) || // 2) first <double>
      !IsPosDouble(split[2])) { // 2) first <double>
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - first line must start with box size (i.e., three positive numbers)\n", field);
    fprintf(stderr, "\033[0m");
    ErrorPrintLine(split, words);
    exit(1);
  }
  (*BoxLength).x = atof(split[0]);
  (*BoxLength).y = atof(split[1]);
  (*BoxLength).z = atof(split[2]); //}}}

  fclose(fr);
} //}}}

// ReadFieldBeadType() //{{{
/**
 * Function reading the species section of the FIELD-like file.
 */
void ReadFieldBeadType(char *field, COUNTS *Counts, BEADTYPE **BeadType, BEAD **Bead) {

  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}

  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t"); // delimiters for SplitLine()

  // read number of bead types //{{{
  bool missing = true; // is 'species' keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    int words = SplitLine(split, line, delim);
    if (strcasecmp(split[0], "species") == 0) {
      missing = false;
      // check if the next string is a number
      if (words < 2 ||            // missing next string
          !IsInteger(split[1])) { // next string isn't a number
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing number of species\n", field);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'species' keyword\n\n", field);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  *BeadType = calloc((*Counts).TypesOfBeads,sizeof(struct BeadType));

  // read info about bead types //{{{
  (*Counts).Unbonded = 0;
  (*Counts).Beads = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    fgets(line, sizeof(line), fr);
    int words = SplitLine(split, line, delim);

    // Error on the line: //{{{
    // 1) empty line
    // 2) less then four strings
    // 3) second string isn't a positive double (mass)
    // 4) third string isn't a double (charge)
    // 5) fifth string isn't an integer (unbonded beads)
    if (words == 1 && split[0][0] == '\0') { // 1)
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - no blank spaces permitted in the species section\n\n", field);
      fprintf(stderr, "\033[0m");
      exit(1);
    } else if (words < 4 ||                  // 2)
               !IsPosDouble(split[1]) ||     // 3)
               !IsDouble(split[2]) ||        // 4)
               !IsInteger(split[3])) {       // 5)
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - wrong species line\n", field);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}

    strcpy((*BeadType)[i].Name, split[0]);
    (*BeadType)[i].Mass = atof(split[1]);
    (*BeadType)[i].Charge = atof(split[2]);
    (*BeadType)[i].Number = atoi(split[3]);
    (*Counts).Unbonded += (*BeadType)[i].Number;
  } //}}}

  (*Counts).Beads = (*Counts).Unbonded;

  // allocate & fill Bead array //{{{
  *Bead = calloc((*Counts).Unbonded, sizeof(struct Bead));
  int count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    for (int j = 0; j < (*BeadType)[i].Number; j++) {
      (*Bead)[count].Type = i;
      (*Bead)[count].Molecule = -1;
      (*Bead)[count].Index = count;
      count++;
    }
  } //}}}

  fclose(fr);
} //}}}

// ReadFieldMolecules() //{{{
/**
 * Function reading the molecules section of the FIELD-like file.
 */
void ReadFieldMolecules(char *field, COUNTS *Counts,
                        BEADTYPE **BeadType, BEAD **Bead,
                        MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
                        PARAMS **bond_type, PARAMS **angle_type) {

  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}

  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t"); // delimiters for SplitLine()

  // read number of molecule types //{{{
  bool missing = true; // is molecule keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    int words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "molecule", 8) == 0) {
      missing = false;
      // error - next string isn't a number
      if (!IsInteger(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing number of molecule types\n", field);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfMolecules = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'molecule' line\n\n", field);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // allocate molecule type struct //{{{
  *MoleculeType = calloc((*Counts).TypesOfMolecules,sizeof(struct MoleculeType));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Number = 0;
    (*MoleculeType)[i].nBeads = 0;
    (*MoleculeType)[i].nBonds = 0;
    (*MoleculeType)[i].BType = calloc((*Counts).TypesOfBeads, sizeof(int));
  } //}}}

  // test proper number of 'finish' keywords //{{{
  fpos_t position;
  fgetpos(fr, &position); // save pointer position
  int count = 0;
  while(fgets(line, sizeof(line), fr)) {
    SplitLine(split, line, delim);
    if (strcmp(split[0], "finish") == 0) {
      count++;
    }
  }
  if (count < (*Counts).TypesOfMolecules) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'finish' keyword\n\n", field);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  fsetpos(fr, &position); //}}}

  // read info about molecule types //{{{
  fgetpos(fr, &position); // save pointer position
  (*Counts).TypesOfBonds = 0;
  (*Counts).TypesOfAngles = 0;
  (*bond_type) = calloc(1, sizeof(struct Params));
  (*angle_type) = calloc(1, sizeof(struct Params));
  // stored bond & angle types - temporary
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // read i-th molecule name //{{{
    fgets(line, sizeof(line), fr);
    SplitLine(split, line, delim);
    if (split[0][0] == '\0') {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing molecule name\n\n", field);
      fprintf(stderr, "\033[0m");
      exit(1);
    }
    strcpy((*MoleculeType)[i].Name, split[0]); //}}}
    // read number of the molecules //{{{
    fgets(line, sizeof(line), fr);
    int words = SplitLine(split, line, delim);
    if (strcasecmp(split[0], "nummols") != 0) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'nummols' keyword\n", field);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'nummols' must be followed by a non-negative integer\n", field);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].Number = atoi(split[1]); //}}}
    // read number of beads in the molecule //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "beads", 4) != 0) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'beads' keyword\n", field);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'beads' must be followed by a non-negative integer\n", field);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].nBeads = atoi(split[1]);
    (*MoleculeType)[i].Bead = calloc((*MoleculeType)[i].nBeads, sizeof(int)); //}}}
    // read bead info //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line, delim);
      // bead line must be 'name <double> <double> <double>'
      if (words < 4 ||
          !IsDouble(split[1]) || !IsDouble(split[2]) || !IsDouble(split[3])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - wrong bead line in molecule ", field);
        fprintf(stderr, "\033[1;33m%s\n", (*MoleculeType)[i].Name);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      int type = FindBeadType(split[0], *Counts, *BeadType);
      if (type == -1) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - non-existent bead name \033[1;33m%s\033[1;31m ", field, split[0]);
        fprintf(stderr, "in molecule \033[1;33m%s\n", (*MoleculeType)[i].Name);
        fprintf(stderr, "\033[0m");
        ErrorBeadType(*Counts, *BeadType);
        exit(1);
      }
      (*MoleculeType)[i].Bead[j] = type;
    } //}}}
    // read number of bonds in the molecule //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "bonds", 4) != 0) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'bonds' keyword\n", field);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bonds' must be followed by a non-negative integer\n", field);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].nBonds = atoi(split[1]);
    (*MoleculeType)[i].Bond = calloc((*MoleculeType)[i].nBonds, sizeof(int *));
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      (*MoleculeType)[i].Bond[j] = calloc(3, sizeof(int));
      (*MoleculeType)[i].Bond[j][2] = -1; // no bond type assigned
    } //}}}
    // read bond info //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line, delim);
      // error - bead line must be '<type> <bead id> <bead id> <bond strenthe> <distance>' //{{{
      if (words < 5 ||
          !IsInteger(split[1]) || atoi(split[1]) == 0 || // bead ids must start from 1
          !IsInteger(split[2]) || atoi(split[2]) == 0 ||
          !IsPosDouble(split[3]) || !IsPosDouble(split[4])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - wrong bond line in molecule %s\n", field, (*MoleculeType)[i].Name);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      (*MoleculeType)[i].Bond[j][0] = atoi(split[1]) - 1;
      (*MoleculeType)[i].Bond[j][1] = atoi(split[2]) - 1;
      // find if the bond type is new
      bool known = false;
      for (int k = 0; k < (*Counts).TypesOfBonds; k++) {
        if ((*bond_type)[k].a == atof(split[3]) && (*bond_type)[k].b == atof(split[4])) {
          known = true;
          break;
        }
      }
      if (!known) {
        (*Counts).TypesOfBonds++;
        (*bond_type) = realloc((*bond_type), (*Counts).TypesOfBonds*sizeof(struct Params));
        (*bond_type)[(*Counts).TypesOfBonds-1].a = atof(split[3]);
        (*bond_type)[(*Counts).TypesOfBonds-1].b = atof(split[4]);
      }
    } //}}}
    // read angle info if present //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // are angles present?
    if (strncasecmp(split[0], "angles", 5) == 0) {
      // get number of angles //{{{
      if (words < 2 || !IsInteger(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'angles' must be followed by a non-negative integer\n", field);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*MoleculeType)[i].nAngles = atoi(split[1]);
      (*MoleculeType)[i].Angle = calloc((*MoleculeType)[i].nAngles, sizeof(int *));
      for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
        (*MoleculeType)[i].Angle[j] = calloc(4, sizeof(int));
        (*MoleculeType)[i].Angle[j][3] = -1; // no angle type assigned
      } //}}}
      // get info about angles
      for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // bead line must be '<type> <bead id> <bead id> <bead id> <angle strength> <angle>'
        if (words < 6 ||
            !IsInteger(split[1]) || atoi(split[1]) == 0 || // bead ids must start from 1
            !IsInteger(split[2]) || atoi(split[2]) == 0 ||
            !IsInteger(split[3]) || atoi(split[3]) == 0 ||
            !IsPosDouble(split[4]) || !IsPosDouble(split[5])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", field);
          fprintf(stderr, " - wrong angle line in molecule \033[1;33m%s\n", (*MoleculeType)[i].Name);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        }
        (*MoleculeType)[i].Angle[j][0] = atoi(split[1]) - 1;
        (*MoleculeType)[i].Angle[j][1] = atoi(split[2]) - 1;
        (*MoleculeType)[i].Angle[j][2] = atoi(split[3]) - 1;
        bool known = false;
        for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
          if ((*angle_type)[k].a == atof(split[4]) && (*angle_type)[k].b == atof(split[5])) {
            known = true;
            break;
          }
        }
        if (!known) {
          (*Counts).TypesOfAngles++;
          (*angle_type) = realloc((*angle_type), (*Counts).TypesOfAngles*sizeof(struct Params));
          (*angle_type)[(*Counts).TypesOfAngles-1].a = atof(split[4]);
          (*angle_type)[(*Counts).TypesOfAngles-1].b = atof(split[5]);
        }
      }
    // error - isn't there more bonds, by any chance?
    } else if (strncasecmp(split[0], "harm", 5) == 0) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m ", field);
      fprintf(stderr, "- extra bond line in molecule \033[1;33m%s\n", (*MoleculeType)[i].Name);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}
    (*Counts).Bonded += (*MoleculeType)[i].Number * (*MoleculeType)[i].nBeads;
    (*Counts).Molecules += (*MoleculeType)[i].Number;
    // skip till 'finish' //{{{
    if (strcasecmp(split[0], "finish") == 0) {
      continue;
    }
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  } //}}}

  // return the file pointer to the beginning of the molecules section
  fsetpos(fr, &position);

  // calculate molecule masses //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      (*MoleculeType)[i].Mass += (*BeadType)[btype].Mass;
    }
  } //}}}

  // update the number beads in the system //{{{
  (*Counts).Beads += (*Counts).Bonded;
  (*Counts).BeadsInVsf = (*Counts).Beads;
  *Bead = realloc(*Bead, (*Counts).Beads*sizeof(struct Bead));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      (*BeadType)[btype].Number += (*MoleculeType)[i].Number;
    }
  } //}}}

  // read coordinates of bonded beads & assign bond and angle types //{{{
  // get the coordinates from FIELD to each molecule
  count = (*Counts).Unbonded;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // skip name, nummols, and beads lines //{{{
    fgets(line, sizeof(line), fr);
    fgets(line, sizeof(line), fr);
    fgets(line, sizeof(line), fr); //}}}
    // read bead lines //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      for (int k = 0; k < (*MoleculeType)[i].Number; k++) {
        int id = count + k * (*MoleculeType)[i].nBeads;
        (*Bead)[id].Position.x = atof(split[1]);
        (*Bead)[id].Position.y = atof(split[2]);
        (*Bead)[id].Position.z = atof(split[3]);
      }
      count++;
    } //}}}
    // set count to the first bead of the next molecule type
    count += ((*MoleculeType)[i].Number - 1) * (*MoleculeType)[i].nBeads;
    // skip bonds line
    fgets(line, sizeof(line), fr);
    // read bond types //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      for (int k = 0; k < (*Counts).TypesOfBonds; k++) {
        if ((*bond_type)[k].a == atof(split[3]) && (*bond_type)[k].b == atof(split[4])) {
          (*MoleculeType)[i].Bond[j][2] = k;
          break;
        }
      }
    } //}}}
    // read angle types //{{{
    for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
      // skip 'angles' line at the beginning of the angle section //{{{
      if (j == 0) {
        fgets(line, sizeof(line), fr);
      } //}}}
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
        if ((*angle_type)[k].a == atof(split[4]) && (*angle_type)[k].b == atof(split[5])) {
          (*MoleculeType)[i].Angle[j][3] = k;
          break;
        }
      }
    } //}}}
    // skip till 'finish' //{{{
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  } //}}}

  // return the file pointer to the beginning of the molecules section
  fsetpos(fr, &position);

  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    SortBonds((*MoleculeType)[i].Bond, (*MoleculeType)[i].nBonds);
    SortAngles((*MoleculeType)[i].Angle, (*MoleculeType)[i].nAngles);
  }

  fclose(fr);
} //}}}

// ReadField() //{{{
/**
 * Function reading the FIELD-like file; it completely fills provided structs,
 * overriding any possible data in there. If the FIELD-like file is a source of
 * some additional structure information, new structs must be used, and the
 * data copied from there afterwards.
 */
void ReadField(char *field, VECTOR *BoxLength, COUNTS *Counts,
               BEADTYPE **BeadType, BEAD **Bead, int **Index,
               MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
               PARAMS **bond_type, PARAMS **angle_type) {

  // read pbc only if required
  if (BoxLength != NULL) {
    ReadFieldPbc(field, BoxLength);
  }
  ReadFieldBeadType(field, Counts, BeadType, Bead);
  ReadFieldMolecules(field, Counts, BeadType, Bead, MoleculeType, Molecule, bond_type, angle_type);

  // allocate Bead[].Aggregate array - needed only to free() //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(10, sizeof(int));
  } //}}}

  // fill MoleculeType[].BType array //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nBTypes = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int type = (*MoleculeType)[i].Bead[j];
      bool present = false;
      for (int k = 0; k < (*MoleculeType)[i].nBTypes; k++) {
        if (type == (*MoleculeType)[i].BType[k]) {
          present = true;
          break;
        }
      }
      if (!present) {
        (*MoleculeType)[i].BType[(*MoleculeType)[i].nBTypes] = type;
        (*MoleculeType)[i].nBTypes++;
      }
    }
  } //}}}

  // fill Molecule & Bead structs //{{{
  *Molecule = calloc((*Counts).Molecules,sizeof(struct Molecule));
  int count_mol = 0, count_bead = (*Counts).Unbonded;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].Number; j++) {
      (*Molecule)[count_mol].Type = i;
      (*Molecule)[count_mol].Bead = calloc((*MoleculeType)[i].nBeads, sizeof(int));
      for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
        int btype = (*MoleculeType)[i].Bead[k];
        (*Molecule)[count_mol].Bead[k] = count_bead;
        (*Bead)[count_bead].Molecule = count_mol;
        (*Bead)[count_bead].Type = btype;
        (*Bead)[count_bead].Index = count_bead;
        count_bead++;
      }
      count_mol++;
    }
  } //}}}

  // fill Index - I don't thinks it's used right now //{{{
  *Index = calloc((*Counts).Beads, sizeof(int));
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Index)[i] = i;
  } //}}}

  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, field);
} //}}}

// ReadLmpData() //{{{
void ReadLmpData(char *data_file, int *bonds, PARAMS **bond_type,
                 int *angles, PARAMS **angle_type,
                 VECTOR *BoxLength, VECTOR *box_lo, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {
  // open data file //{{{
  FILE *fr;
  if ((fr = fopen(data_file, "r")) == NULL) {
    ErrorFileOpen(data_file, 'r');
    exit(1);
  } //}}}

  // ignore first line (comment) //{{{
  char line[LINE];
  fgets(line, sizeof(line), fr); //}}}

  // read data file header //{{{
  // data file header lines must start with a number (or '#' for comment),
  // therefore read until something else is encountered
  char split[30][100], delim[8];
  int words;
  strcpy(delim, " \t");
  do {
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // number of atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;33m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'atoms' keyword must be preceded by integer\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).BeadsInVsf = atoi(split[0]);
      (*Counts).Beads = atoi(split[0]);
    } //}}}
    // number of bonds //{{{
    if (words > 1 && strcmp(split[1], "bonds") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bonds' keyword must be preceded by integer\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      *bonds = atoi(split[0]);
    } //}}}
    // number of angles //{{{
    if (words > 1 && strcmp(split[1], "angles") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'angles' keyword must be preceded by integer\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      *angles = atoi(split[0]);
    } //}}}
    // number of bead types //{{{
    if (words > 2 && strcmp(split[1], "atom") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'atom types' keyword must be preceded by integer\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[0]);
    } //}}}
    // number of bond types //{{{
    if (words > 2 && strcmp(split[1], "bond") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bond types' keyword must be preceded by integer\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBonds = atoi(split[0]);
    } //}}}
    // number of angle types //{{{
    if (words > 2 && strcmp(split[1], "angle") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'angle types' keyword must be preceded by integer\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfAngles = atoi(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "xlo") == 0 && strcmp(split[3], "xhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'xlo xhi' keyword must be preceded by two floats\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).x = atof(split[1]) - atof(split[0]);
      (*box_lo).x = atof(split[0]);
    } //}}}
    // box length in y //{{{
    if (words > 3 && strcmp(split[2], "ylo") == 0 && strcmp(split[3], "yhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'ylo yhi' keyword must be preceded by two floats\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).y = atof(split[1]) - atof(split[0]);
      (*box_lo).y = atof(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "zlo") == 0 && strcmp(split[3], "zhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'zlo zhi' keyword must be preceded by two floats\n", data_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).z = atof(split[1]) - atof(split[0]);
      (*box_lo).z = atof(split[0]);
    } //}}}
  } while (words == 0 ||
           split[0][0] == '#' ||
           IsDouble(split[0]) ||
           IsInteger(split[0])); //}}}

  // some error checking //{{{
  if ((*Counts).TypesOfBeads == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'atom types' line (or is 0)\n\n", data_file);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if ((*Counts).BeadsInVsf == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'atoms' line (or is 0)\n\n", data_file);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if ((*BoxLength).x == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'xlo xhi' line (or is 0 0)\n\n", data_file);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if ((*BoxLength).y == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'ylo yhi' line (or is 0 0)\n\n", data_file);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if ((*BoxLength).z == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'zlo zhi' line (or is 0 0)\n\n", data_file);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // fill something in BeadType struct //{{{
  *BeadType = calloc((*Counts).TypesOfBeads, sizeof(struct BeadType));
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    sprintf((*BeadType)[i].Name, "bead%d", i+1);
    (*BeadType)[i].Use = true;
    (*BeadType)[i].Write = true;
  } //}}}

  // bead struct memory allocation //{{{
  *Bead = calloc((*Counts).Beads, sizeof(struct Bead));
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(1, sizeof(int));
  } //}}}

  *Index = calloc((*Counts).Beads, sizeof(int)); // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  *bond_type = calloc((*Counts).TypesOfBonds, sizeof(PARAMS));
  *angle_type = calloc((*Counts).TypesOfAngles, sizeof(PARAMS));

  // read body of data file //{{{
  int test,
      *mols, // number of beads in each molecule
      monomer = 0; // monomer beads designated by mol_ID = 0 in data file
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);
    // atom masses //{{{
    if (words > 0 && strcmp(split[0], "Masses") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 2 || !IsInteger(split[0]) || !IsPosDouble(split[1])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each line in 'Masses' section must start with '<int> <float>'\n", data_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*BeadType)[i].Mass = atof(split[1]);
        // if there's a comment at the end of the line, consider it bead name
        if (words > 2 && split[2][0] == '#') {
          if (strlen(split[2]) == 1 && words > 3) { // comment form '# name'
            strcpy((*BeadType)[i].Name, split[3]);
          } else if (strlen(split[2]) > 1) { // comment form '#name'
            for (int j = 0; j < strlen(split[2]); j++) {
              split[2][j] = split[2][j+1];
            }
            strcpy((*BeadType)[i].Name, split[2]);
          }
        }
      }
    } //}}}
    // bond coefficients //{{{
    if (words > 1 && strcmp(split[0], "Bond") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfBonds; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosDouble(split[1]) || !IsPosDouble(split[2])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each line in 'Bond Coeffs' section must start with '<int> <float> <float>'\n", data_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*bond_type)[i].a = atof(split[1]);
        (*bond_type)[i].b = atof(split[2]);
      }
    } //}}}
    // angle coefficients //{{{
    if (words > 1 && strcmp(split[0], "Angle") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfAngles; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosDouble(split[1]) || !IsPosDouble(split[2])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each line in 'Angle Coeffs' section must start with '<int> <float> <float>'\n", data_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*angle_type)[i].a = atof(split[1]);
        (*angle_type)[i].b = atof(split[2]);
      }
    } //}}}
    // atoms section //{{{
    if (words > 0 && strcmp(split[0], "Atoms") == 0) {
      // array for number of beads in each molecule
      mols = calloc((*Counts).Beads, sizeof(int));
      // array for list of beads in each molecule
      // bead_mols[i] ... molecule's id; bead_mols[][i] ... molecule's beads
      int **bead_mols = calloc((*Counts).Beads, sizeof(int *));
      for (int i = 0; i < (*Counts).Beads; i++) {
        bead_mols[i] = calloc(1, sizeof(int));
      }
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // go through atom section to get basic info //{{{
      fpos_t pos; // set file counter
      fgetpos(fr, &pos); // save file pointer
      for (int i = 0; i < (*Counts).Beads; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <id> <mol_id> <btype> <charge> <x> <y> <z>
        // Error - incorrect format //{{{
        if (words < 7 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) || !IsInteger(split[2]) ||
            !IsDouble(split[3]) || !IsDouble(split[4]) || !IsDouble(split[5]) || !IsDouble(split[6])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each 'Atoms' line must be <id> <mol_id> <bead type> <charge> <x> <y> <z>\n", data_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int id = atoi(split[0]) - 1, // in lammps, these start at 1
            mol_id = atoi(split[1]) - 1, // in lammps, molecules start with 1; unbonded atoms can be 0
            btype = atoi(split[2]) - 1; // in lammps, these start at 1

        if (mol_id == -1) { // corresponds to 0 in the data file
          monomer++;
          (*Bead)[id].Molecule = -1;
        } else { // possibly in a molecule (if more beads share its mol_id)
          mols[mol_id]++;
          bead_mols[mol_id] = realloc(bead_mols[mol_id], mols[mol_id]*sizeof(int));
          bead_mols[mol_id][mols[mol_id]-1] = id;
          (*Bead)[id].Molecule = mol_id;
        }
        (*BeadType)[btype].Charge = atof(split[3]);
        (*BeadType)[btype].Number++;
        (*Bead)[id].Position.x = atof(split[4]) - (*box_lo).x;
        (*Bead)[id].Position.y = atof(split[5]) - (*box_lo).y;
        (*Bead)[id].Position.z = atof(split[6]) - (*box_lo).z;
        (*Bead)[id].Type = btype;
        (*Bead)[id].Index = id;
        (*Index)[id] = id;
      } //}}}
      (*Counts).Unbonded = monomer;
      // go through possible molecules and remove 1-bead molecules //{{{
      int count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < (*Counts).Beads; i++) {
        if (mols[i] == 1) {
          int bead = bead_mols[i][0];
          (*Bead)[bead].Molecule = -1;
          (*Counts).Unbonded++;
        } else if (mols[i] > 1){
          (*Counts).Molecules++;
          (*Counts).Bonded += mols[i];
          for (int j = 0; j < mols[i]; j++) {
            int bead = bead_mols[i][j];
            (*Bead)[bead].Molecule = count;
          }
          count++;
        }
      } //}}}
      // TODO: is the 'remove 1-bead...' necessary? Join with 'allocate Molecule struct...'
      // remove 1-bead molecules from mols array and sort bead_mols arrayy according to ascending bead index //{{{
      count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < (*Counts).Beads; i++) {
        if (mols[i] > 1) {
          mols[count] = mols[i];
          // realloc bead_mols[count] to bead_mols[i] as it realloc'd in the first place
          bead_mols[count] = realloc(bead_mols[count], mols[count]*sizeof(int));
          for (int j = 0; j < mols[i]; j++) {
            bead_mols[count][j] = bead_mols[i][j];
          }
          // sort molecules in bead_mols[count][] according to ascending id
          SortArray(&bead_mols[count], mols[i], 0);
          count++;
        }
      } //}}}
      // zeroize unused part of mols array - just to be on the save side //{{{
      for (int i = (*Counts).Molecules; i < (*Counts).Beads; i++) {
        mols[i] = 0;
      } //}}}
      // allocate Molecule struct and fill Molecule[].Bead array //{{{
      (*Molecule) = calloc((*Counts).Molecules, sizeof(struct Molecule));
      for (int i = 0; i < (*Counts).Molecules; i++) {
        (*Molecule)[i].Bead = calloc(mols[i], sizeof(int));
        for (int j = 0; j < mols[i]; j++) {
          (*Molecule)[i].Bead[j] = bead_mols[i][j];
        }
      } //}}}
      // free helper array  //{{{
      for (int i = 0; i < (*Counts).Beads; i++) {
        free(bead_mols[i]);
      }
      free(bead_mols); //}}}
    } //}}}
    // bonds section //{{{
    if (words > 0 && strcmp(split[0], "Bonds") == 0) {
      // allocate helper arrays to hold bond info //{{{
      // number of bonds in each molecule
      int *bonds_per_mol = calloc((*Counts).Molecules, sizeof(int));
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] ... ids of connected beads in molecule 'i'
      // [i][][2] ... bond type
      int ***mol_bonds = calloc((*Counts).Molecules, sizeof(int **));
      for (int i = 0; i < (*Counts).Molecules; i++) {
         mol_bonds[i] = calloc(1, sizeof(int *));
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // read all bonds //{{{
      for (int i = 0; i < *bonds; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <bond id> <bond type> <bead1> <bead2>
        // Error - incorrect format //{{{
        if (words < 4 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each 'Bonds' line must be <bond id> <bond type> <bead1d> <bead2>\n", data_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        // assign molecule to the bond
        int mol = (*Bead)[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != (*Bead)[bead2].Molecule) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", data_file);
          fprintf(stderr, " - beads in \033[1;33mbond %d\033[1;31m are in different molecules\n\n", atoi(split[0]));
          fprintf(stderr, "\033[0m");
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        bonds_per_mol[mol]++;
        int bond = bonds_per_mol[mol];
        // add bonded beads to the molecule they belong to
        mol_bonds[mol] = realloc(mol_bonds[mol], bond*sizeof(int *));
        mol_bonds[mol][bond-1] = calloc(3, sizeof(int));
        mol_bonds[mol][bond-1][0] = bead1;
        mol_bonds[mol][bond-1][1] = bead2;
        mol_bonds[mol][bond-1][2] = type;
      } //}}}
      // sort bonds according to the id of the first bead in a bond //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        SortBonds(mol_bonds[i], bonds_per_mol[i]);
      } //}}}
      // minimize mol_bonds based on lowest id in each molecule //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int lowest = (*Counts).Beads; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if ((*Molecule)[i].Bead[j] < lowest) {
            lowest = (*Molecule)[i].Bead[j];
          }
        }
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          mol_bonds[i][j][0] -= lowest;
          mol_bonds[i][j][1] -= lowest;
        }
      } //}}}
      // identify molecule type based on bead order and connectivity //{{{
      *MoleculeType = calloc(1, sizeof(struct MoleculeType));
      // number of molecule types differing in connectivity but not in beads - naming purposes
      int *diff_conn = calloc(1, sizeof(int));
      // number of molecule types differing in bead order - naming purposes
      int mtype_bead = 0;
      for (int i = 0; i < (*Counts).Molecules; i++) {
        // is molecule 'i' of known type?
        bool exists = false;
        // do 'i' and given molecule type share connectivity and bead order?
        bool same_conn = true, same_bead = true;
        // molecule type with which 'i' shares bead order - naming purposes
        int type_bead = -1;
        for (int j = 0; j < (*Counts).TypesOfMolecules; j++) { //{{{
          if ((*MoleculeType)[j].nBeads == mols[i] && // same number of molecules?
              (*MoleculeType)[j].nBonds == bonds_per_mol[i]) { // same number of bonds?
            // same connectivity?
            same_conn = true;
            for (int k = 0; k < (*MoleculeType)[j].nBonds; k++) {
              if ((*MoleculeType)[j].Bond[k][0] != mol_bonds[i][k][0] ||
                  (*MoleculeType)[j].Bond[k][1] != mol_bonds[i][k][1] ||
                  (*MoleculeType)[j].Bond[k][2] != mol_bonds[i][k][2]) {
                same_conn = false;
                break;
              }
            }
            // same bead types?
            same_bead = true;
            for (int k = 0; k < (*MoleculeType)[j].nBeads; k++) {
              int btype = (*Bead)[(*Molecule)[i].Bead[k]].Type;
              if ((*MoleculeType)[j].Bead[k] != btype) {
                same_bead = false;
                break;
              }
            }
            if (same_bead && type_bead == -1) {
              type_bead = j;
            }
            // if the molecule has the same connectivity and bead order, it's of known type
            if (same_conn && same_bead) {
              exists = true;
              (*MoleculeType)[j].Number++;
              (*Molecule)[i].Type = j;
              break;
            }
          }
        } //}}}
        // add new type? //{{{
        if (!exists) {
          int mtype = (*Counts).TypesOfMolecules;
          *MoleculeType = realloc(*MoleculeType, (mtype+1)*sizeof(struct MoleculeType));
          diff_conn = realloc(diff_conn, (mtype+1)*sizeof(int));
          diff_conn[mtype] = 0;
          // molecule name
          if (same_bead && !same_conn) { // same beads - mol<type_bead>-b#
            diff_conn[type_bead]++;
            char name[5];
            strcpy(name, (*MoleculeType)[type_bead].Name);
            sprintf((*MoleculeType)[mtype].Name, "%s-b%d", name, diff_conn[type_bead]);
          } else { // same connectivity or both different - mol#
            mtype_bead++;
            sprintf((*MoleculeType)[mtype].Name, "mol%d", mtype_bead);
          }
          (*Molecule)[i].Type = mtype;
          (*MoleculeType)[mtype].Number = 1;
          // copy bead sequence and determine BTypes stuff
          (*MoleculeType)[mtype].nBeads = mols[i];
          (*MoleculeType)[mtype].Bead = calloc((*MoleculeType)[mtype].nBeads, sizeof(int));
          (*MoleculeType)[mtype].nBTypes = 0;
          (*MoleculeType)[mtype].BType = calloc(1, sizeof(int));
          for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
            int btype = (*Bead)[(*Molecule)[i].Bead[j]].Type;
            (*MoleculeType)[mtype].Bead[j] = btype;
            exists = false; // recycling the bool to check if btype is already in BType[]
            for (int k = 0; k < (*MoleculeType)[mtype].nBTypes; k++) {
              if ((*MoleculeType)[mtype].BType[k] == btype) {
                exists = true;
              }
            }
            if (!exists) { // recycled exists
              int types = (*MoleculeType)[mtype].nBTypes;
              (*MoleculeType)[mtype].nBTypes++;
              (*MoleculeType)[mtype].BType = realloc((*MoleculeType)[mtype].BType, (types+1)*sizeof(int));
              (*MoleculeType)[mtype].BType[types] = btype;
            }
          }
          // copy bonds
          (*MoleculeType)[mtype].nBonds = bonds_per_mol[i];
          (*MoleculeType)[mtype].Bond = calloc((*MoleculeType)[mtype].nBonds, sizeof(int *));
          for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
            (*MoleculeType)[mtype].Bond[j] = calloc(3, sizeof(int));
            (*MoleculeType)[mtype].Bond[j][0] = mol_bonds[i][j][0];
            (*MoleculeType)[mtype].Bond[j][1] = mol_bonds[i][j][1];
            (*MoleculeType)[mtype].Bond[j][2] = mol_bonds[i][j][2];
          }
          (*MoleculeType)[mtype].nAngles = 0;
          (*MoleculeType)[mtype].Angle = calloc(1, sizeof(int *));
          (*MoleculeType)[mtype].Write = true;
          (*Counts).TypesOfMolecules++;
        } //}}}
      } //}}}
      // free helper arrays //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          free(mol_bonds[i][j]);
        }
        free(mol_bonds[i]);
      }
      free(mol_bonds);
      free(bonds_per_mol);
      free(diff_conn); //}}}
    } //}}}
    // angles section //{{{
    if (words > 0 && strcmp(split[0], "Angles") == 0) {
      // allocate helper arrays to hold angle info //{{{
      // number of angles in each molecule
      int *angles_per_mol = calloc((*Counts).Molecules, sizeof(int));
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] & [i][][2] ... ids of connected beads in molecule 'i'
      // [i][][3] ... angle type
      int ***mol_angles = calloc((*Counts).Molecules, sizeof(int **));
      for (int i = 0; i < (*Counts).Molecules; i++) {
         mol_angles[i] = calloc(1, sizeof(int *));
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // read all angles //{{{
      for (int i = 0; i < *angles; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <angle id> <angle type> <bead1> <bead2> <bead3>
        // Error - incorrect format //{{{
        if (words < 5 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3]) || !IsInteger(split[4])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each 'Angles' line must be <angle id> <angle type> <bead1d> <bead2> <bead3>\n", data_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        int bead3 = atoi(split[4]) - 1;
        // assign molecule to the bond
        int mol = (*Bead)[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != (*Bead)[bead2].Molecule || mol != (*Bead)[bead3].Molecule) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m ", data_file);
          fprintf(stderr, " - atoms in \033[1;33mangle %d\033[1;31m are in different molecules\n\n", atoi(split[0]));
          fprintf(stderr, "\033[0m");
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        angles_per_mol[mol]++;
        int angle = angles_per_mol[mol];
        // add angle beads to the molecule they belong to
        mol_angles[mol] = realloc(mol_angles[mol], angle*sizeof(int *));
        mol_angles[mol][angle-1] = calloc(4, sizeof(int));
        mol_angles[mol][angle-1][0] = bead1;
        mol_angles[mol][angle-1][1] = bead2;
        mol_angles[mol][angle-1][2] = bead3;
        mol_angles[mol][angle-1][3] = type;
      } //}}}
      // sort angles according to the id of the first bead in a bond //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        SortAngles(mol_angles[i], angles_per_mol[i]);
      } //}}}
      // minimize mol_angles based on lowest id in each molecule //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int lowest = (*Counts).Beads; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if ((*Molecule)[i].Bead[j] < lowest) {
            lowest = (*Molecule)[i].Bead[j];
          }
        }
        for (int j = 0; j < angles_per_mol[i]; j++) {
          mol_angles[i][j][0] -= lowest;
          mol_angles[i][j][1] -= lowest;
          mol_angles[i][j][2] -= lowest;
        }
      } //}}}
      // add angles to existing molecule types (or create a new one differing only in angles)
      // ...'basic' molecule types already generated in bonds section
      int *extra = calloc((*Counts).TypesOfMolecules, sizeof(int)); // number of molecule types differing in angles
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int mtype = (*Molecule)[i].Type;
        if ((*MoleculeType)[mtype].nAngles == 0) { // add angles if there are no angles in the molecule type //{{{
          (*MoleculeType)[mtype].nAngles = angles_per_mol[i];
          (*MoleculeType)[mtype].Angle = realloc((*MoleculeType)[mtype].Angle, (*MoleculeType)[mtype].nAngles*sizeof(int *));
          for (int j = 0; j < (*MoleculeType)[mtype].nAngles; j++) {
            (*MoleculeType)[mtype].Angle[j] = calloc(4, sizeof(int));
            (*MoleculeType)[mtype].Angle[j][0] = mol_angles[i][j][0];
            (*MoleculeType)[mtype].Angle[j][1] = mol_angles[i][j][1];
            (*MoleculeType)[mtype].Angle[j][2] = mol_angles[i][j][2];
            (*MoleculeType)[mtype].Angle[j][3] = mol_angles[i][j][3];
          } //}}}
        } else { // if angles already present, check whether they're the same //{{{
          bool exists = true;
          // same number of angles? //{{{
          if ((*MoleculeType)[mtype].nAngles != angles_per_mol[i]) {
            exists = false;
          } //}}}
          // same angles as in mtype? //{{{
          for (int j = 0; j < angles_per_mol[i] && j < (*MoleculeType)[mtype].nAngles; j++) {
            if ((*MoleculeType)[mtype].Angle[j][0] != mol_angles[i][j][0] ||
                (*MoleculeType)[mtype].Angle[j][1] != mol_angles[i][j][1] ||
                (*MoleculeType)[mtype].Angle[j][2] != mol_angles[i][j][2] ||
                (*MoleculeType)[mtype].Angle[j][3] != mol_angles[i][j][3] ) {
              exists = false;
            }
          } //}}}
          // check against angle-generated molecule types //{{{
          // If its angles aren't the same as in mtype, check against other
          // types (i.e., against newly generated thanks to different angles)
          if (!exists) {
            for (int j = (mtype+1); j < (*Counts).TypesOfMolecules; j++) {
              if ((*MoleculeType)[j].nAngles == angles_per_mol[i]) {
                int count = 0;
                for (int k = 0; k < (*MoleculeType)[j].nAngles; k++) {
                  if ((*MoleculeType)[j].Angle[k][0] == mol_angles[i][k][0] &&
                      (*MoleculeType)[j].Angle[k][1] == mol_angles[i][k][1] &&
                      (*MoleculeType)[j].Angle[k][2] == mol_angles[i][k][2] &&
                      (*MoleculeType)[j].Angle[k][3] == mol_angles[i][k][3] ) {
                    count++;
                  }
                }
                if (count == angles_per_mol[i]) {
                  exists = true;
                  (*MoleculeType)[j].Number++;
                  (*MoleculeType)[mtype].Number--;
                  (*Molecule)[i].Type = j;
                  break;
                }
              }
            }
          } //}}}
          // create a new molecule type //{{{
          if (!exists) {
            int new = (*Counts).TypesOfMolecules;
            extra[mtype]++;
            (*Molecule)[i].Type = new;
            (*Counts).TypesOfMolecules++;
            *MoleculeType = realloc(*MoleculeType, (new+1)*sizeof(struct MoleculeType));
            char name[19];
            sprintf(name, "%s-a%d", (*MoleculeType)[mtype].Name, extra[mtype]);
            strcpy((*MoleculeType)[new].Name, name);
            (*MoleculeType)[new].Number = 1;
            (*MoleculeType)[mtype].Number--;
            (*MoleculeType)[new].nBeads = (*MoleculeType)[mtype].nBeads;
            (*MoleculeType)[new].Bead = calloc((*MoleculeType)[new].nBeads, sizeof(int));
            for (int j = 0; j < (*MoleculeType)[new].nBeads; j++) {
              (*MoleculeType)[new].Bead[j] = (*MoleculeType)[mtype].Bead[j];
            }
            (*MoleculeType)[new].nBonds = (*MoleculeType)[mtype].nBonds;
            (*MoleculeType)[new].Bond = calloc((*MoleculeType)[new].nBonds, sizeof(int *));
            for (int j = 0; j < (*MoleculeType)[new].nBonds; j++) {
              (*MoleculeType)[new].Bond[j] = calloc(3, sizeof(int));
              (*MoleculeType)[new].Bond[j][0] = (*MoleculeType)[mtype].Bond[j][0];
              (*MoleculeType)[new].Bond[j][1] = (*MoleculeType)[mtype].Bond[j][1];
              (*MoleculeType)[new].Bond[j][2] = (*MoleculeType)[mtype].Bond[j][2];
            }
            (*MoleculeType)[new].nAngles = (*MoleculeType)[mtype].nAngles;
            (*MoleculeType)[new].Angle = calloc((*MoleculeType)[new].nAngles, sizeof(int *));
            for (int j = 0; j < (*MoleculeType)[new].nAngles; j++) {
              (*MoleculeType)[new].Angle[j] = calloc(4, sizeof(int));
              (*MoleculeType)[new].Angle[j][0] = mol_angles[i][j][0];
              (*MoleculeType)[new].Angle[j][1] = mol_angles[i][j][1];
              (*MoleculeType)[new].Angle[j][2] = mol_angles[i][j][2];
              (*MoleculeType)[new].Angle[j][3] = mol_angles[i][j][3];
            }
            (*MoleculeType)[new].nBTypes = (*MoleculeType)[mtype].nBTypes;
            (*MoleculeType)[new].BType = calloc((*MoleculeType)[new].nBTypes, sizeof(int));
            for (int j = 0; j < (*MoleculeType)[new].nBTypes; j++) {
              (*MoleculeType)[new].BType[j] = (*MoleculeType)[mtype].BType[j];
            }
            (*MoleculeType)[new].Mass = (*MoleculeType)[mtype].Mass;
            (*MoleculeType)[new].Write = true;
          } //}}}
        } //}}}
      }
      // free helper arrays //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        for (int j = 0; j < angles_per_mol[i]; j++) {
          free(mol_angles[i][j]);
        }
        free(mol_angles[i]);
      }
      free(mol_angles);
      free(angles_per_mol);
      free(extra); //}}}
    } //}}}
    // read and split next line
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
  }
  free(mols); //}}}

  // calculate molecular mass //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      (*MoleculeType)[i].Mass += (*BeadType)[(*MoleculeType)[i].Bead[j]].Mass;
    }
  } //}}}

  WarnElNeutrality(*Counts, *BeadType, data_file);

  fclose(fr);
} //}}}
