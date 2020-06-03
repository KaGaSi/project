#include "Read.h"

// GetPBC() //{{{
/*
 * Function to get box dimensions from the provided coordinate file.
 */
Vector GetPBC(FILE *vcf, char *input_coor) {

  Vector BoxLength;

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
      fprintf(stderr, "\nError - %s: unrecognised line '%s'\n", input_coor, line2);
      if (line2[0] >= '0' && line2[0] <= '9' && words > 2) {
        fprintf(stderr, "        Possibly missing pbc line\n");
      }
      putc('\n', stderr);
      exit(1);
    }
  };

  return BoxLength;
} //}}}

// ReadAggCommand() //{{{
void ReadAggCommand(BeadType *BeadType, Counts Counts,
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

  // errof if not enough striings for an Aggregate command
  if (words < 6) {
    fprintf(stderr, "\nError - %s: first line must contain a valid Aggregates command\n", input_agg);
    ErrorPrintLine(split, words);
    exit(1);
  }

  // read minimum distance for closeness check (<distance> argument in Aggregates utility)
  if (!IsPosDouble(split[2])) {
    fprintf(stderr, "\nError - %s: <distance> from Aggregate command must be a non-negative real number\n", input_agg);
    ErrorPrintLine(split, words);
    exit(1);
  }
  *distance = atof(split[2]);
  // read number of contacts aggregate check (<contacts> argument in Aggregates utility)
  if (!IsInteger(split[3])) {
    fprintf(stderr, "\nError - %s: <contacts> from Aggregate command must be a non-negative whole number\n", input_agg);
    ErrorPrintLine(split, words);
    exit(1);
  }
  *contacts = atof(split[3]);
  // warn if a differently named vcf file is used than the one in agg file
  if (strcmp(split[1], input_coor) != 0) {
    fprintf(stderr, "\nWARNING: the coordinate file (%s) ", input_coor);
    fprintf(stderr, "is different to the one in the aggregate file (%s).\n", split[1]);
    fprintf(stderr, "         Mismatch between beads present in both files can lead to undefined behaviour.\n\n");
  }

  // read <type names> from Aggregates command
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType(split[i], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\nError: %s - non-existent bead name '%s' in Aggregate command\n", input_agg, split[i]);
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    BeadType[type].Use = true;
  }

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
bool ReadStructure(char *vsf_file, char *vcf_file, Counts
    *Counts, BeadType **BeadType, Bead **Bead, int **Index,
    MoleculeType **MoleculeType, Molecule **Molecule) {

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
      atom_lines = 0, // number of atom lines in vsf (including comments and blanks)
      type_default = -1; // default bead type
  char line[LINE];
  while(fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    int words = SplitLine(split, line, delim);

    if (strncmp(split[0], "atom", 1) == 0) {
      atom_lines++;

      // error - odd number of values: a(tom) lines are composed of 'keyword <value>' pairs //{{{
      if ((words % 2) == 1) { // because it's split[0] ... split[2n-1]
        fprintf(stderr, "Error: %s - odd number of strings an atom line\n", vsf_file);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}

      // line begins with 'a(tom) default' or 'a(tom) <int>' //{{{
      bool name = false; // 'name' is necessary
      double charge = 1000, mass = 1000;
      if (strcmp("default", split[1]) == 0) {
        type_default = (*Counts).TypesOfBeads;
        // increment number of bead types
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
              fprintf(stderr, "Error: %s - atom charge must be a real number\n", vsf_file);
              ErrorPrintLine(split, words);
              exit(1);
            }
            charge = atof(split[j+1]);
          } else if (strncmp(split[j], "mass", 1) == 0) {
            if (!IsPosDouble(split[j+1])) {
              fprintf(stderr, "Error: %s - atom charge must be a non-negative real number\n", vsf_file);
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
          fprintf(stderr, "Error: %s - atom line must contain 'name'\n", vsf_file);
          ErrorPrintLine(split, words);
          exit(1);
        }
        continue; // skip the rest of the line
      } else if (!IsInteger(split[1])) {
        fprintf(stderr, "Error: %s - 'a(tom)' must be followed by 'default' or a whole number\n", vsf_file);
        ErrorPrintLine(split, words);
        exit(1);
      } else if (atoi(split[1]) > max_bead) {
        max_bead = atoi(split[1]);
      } //}}}

      // check that line contains either both 'resid' and 'resname' or neither //{{{
      // TODO: only resid is necessary, if not resname, add some autogenerated, creating a new MoleculeType
      int test = 0;
      for (int i = 2; i < words && test < 2; i+= 2) {
        if (strcmp("resname", split[i]) == 0) {
          test++;
        } else if (strcmp("resid", split[i]) == 0) {
          test++;
          if (!IsInteger(split[i+1]) || atoi(split[i+1]) == 0) {
            fprintf(stderr, "Error: %s - 'resid' must be followed by a positive whole number\n", vsf_file);
            ErrorPrintLine(split, words);
            exit(1);
          }
        }
      }
      if (test == 1) {
        fprintf(stderr, "Error: %s - atom line must contain both 'resid' and 'resname'\n", vsf_file);
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
          int type;
          if ((type = FindBeadType(split[i+1], *Counts, *BeadType)) == -1) {
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
            (*MoleculeType)[type].Bead = calloc(1, sizeof(int));
            (*MoleculeType)[type].nBTypes = 0;
            // copy new name to MoleculeType[].Name
            strcpy((*MoleculeType)[type].Name, split[i+1]);
          }
          //}}}
        } else if (strcmp("resid", split[i]) == 0) { // molecule id //{{{
          if (atoi(split[i+1]) > max_mol) {
            if (!IsInteger(split[i+1]) || atoi(split[i+1]) == 0) {
              fprintf(stderr, "Error: %s - 'resid' must be followed by a positive whole number\n", vsf_file);
              ErrorPrintLine(split, words);
              exit(1);
            }
            max_mol = atoi(split[i+1]);
          }
        //}}}
        } else if (strcmp("charge", split[i]) == 0 || strcmp("q", split[i]) == 0) { //{{{
          if (!IsDouble(split[i+1])) {
            fprintf(stderr, "Error: %s - atom charge must be a real number\n", vsf_file);
            ErrorPrintLine(split, words);
            exit(1);
          }
          charge = atof(split[i+1]);
        //}}}
        } else if (strncmp("mass", split[i], 1) == 0) { //{{{
          if (!IsPosDouble(split[i+1])) {
            fprintf(stderr, "Error: %s - atom mass must be a non-negative real number\n", vsf_file);
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
        fprintf(stderr, "Error: %s - atom line must contain 'name'\n", vsf_file);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
    } else if (split[0][0] == '#' || split[0][0] == '\0') {
      // count the line if comment or blank
      atom_lines++;
    } else {
      break;
    }
  }
  fclose(vsf); //}}}

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
  for (int count = 0; count < atom_lines; count++) {
    // read line
    fgets(line, sizeof(line), vsf);
    char split[30][100];
    int words = SplitLine(split, line, delim);

    // go through the line
    int bead_id = -1, bead_type = -1, mol_id = -1, mol_type = -1;
    if (split[0][0] == 'a' && strcmp("default", split[1]) != 0) { // non-default a(tom) line
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
        (*Molecule)[mol_id].Type = mol_type;
        (*MoleculeType)[mol_type].nBeads++;
        (*MoleculeType)[mol_type].Bead = realloc((*MoleculeType)[mol_type].Bead, (*MoleculeType)[mol_type].nBeads*sizeof(int));
        (*MoleculeType)[mol_type].Bead[(*MoleculeType)[mol_type].nBeads-1] = bead_type;
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
  int max = 0; // the highest number of beads in a molecule - for memory allocation later
  // count number of molecules for each molecule type
  int *mols;
  mols = calloc((*Counts).TypesOfMolecules, sizeof(int));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    mols[(*Molecule)[i].Type]++;
  }
  // assign info to MoleculeType
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Number = mols[i];
    (*MoleculeType)[i].nBeads /= (*MoleculeType)[i].Number;
    if ((*MoleculeType)[i].nBeads > max) {
      max = (*MoleculeType)[i].nBeads;
    }
  }
  free(mols); //}}}

  // allocate memory for Molecule[].Bead array //{{{
  // use maximum number of beads in any molecule - no molecule can have more
  for (int i = 0; i < mol_alloced; i++) {
    (*Molecule)[i].Bead = calloc(max, sizeof(int));
  } //}}}

  // save pointer position in the vsf file //{{{
  fpos_t pos;
  fgetpos(vsf, &pos); // }}}

  // third, go through the bonds section of vsf to find the number of bonds in molecules //{{{
  int bonds[(*Counts).TypesOfMolecules]; // helper array //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  } //}}}
  while (fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    int words = SplitLine(split, line, bond_delim);

    if (strncmp(split[0], "bond", 1) == 0) {
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
    } else if (split[0][0] != '\n' && // empty line
               split[0][0] != '#') { // comment
      fprintf(stderr, "\nError - %s: unrecognised line\n", vsf_file);
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
      (*MoleculeType)[i].Bond[j] = calloc(2, sizeof(int));
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
  int *beads; // helper array
  beads = calloc((*Counts).BeadsInVsf, sizeof(int));
  for (int count = 0; count < atom_lines; count++) {
    fgets(line, sizeof(line), vsf);
    char split[30][100];
    int words = SplitLine(split, line, delim);

    // go through the line
    int bead_id = -1, mol_id;
    if (strncmp(split[0], "atom", 1) == 0 && strcmp("default", split[1]) != 0) { // non-default a(tom) line
      bead_id = atoi(split[1]);
      // by twos - first is always keyword, second is always value
      for (int i = 2; i < words; i += 2) {
        if (strcmp("resid", split[i]) == 0) { // molecule id
          mol_id = atoi(split[i+1]) - 1; // mol ids start with 1 in vsf
          (*Molecule)[mol_id].Bead[beads[mol_id]++] = bead_id;
          break;
        }
      }
    }
  }
  free(beads); //}}}

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
    } else if (split[0][0] != '\n' && split[0][0] != '#') {
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
    char line[LINE], line2[LINE], *split[32], str[32];
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
          // remove newline if present at the end of split[1]
          if (split[1][strlen(split[1])-1] == '\n') {
            split[1][strlen(split[1])-1] = '\0';
          }
          fprintf(stderr, "\nError: %s - unrecognised keywords '%s %s'", vsf_file, split[0], split[1]);
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
        fprintf(stderr, "Error: %s - blank line instead of a first coordinate line\n", vcf_file);
        exit(1);
      }

      // split the line into array //{{{
      char *split[30];
      split[0] = strtok(line, " \t");
      int i = 0;
      while (split[i] != NULL && i < 29) {
        split[++i] = strtok(NULL, " \t");
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
        int i = 0;
        while (split[i] != NULL && i < 29) {
          split[++i] = strtok(NULL, " \t");
        } //}}}
      } //}}}

      // count the number of beads in vcf //{{{
      (*Counts).Beads = 0;
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        if ((*BeadType)[i].Use) {
          (*Counts).Beads += (*BeadType)[i].Number;
        }
      } //}}}
    //}}}
    } else { // error //{{{
      fprintf(stderr, "Error: %s - missing 'i(ndexed)' or 'o(rdered)' keyword\n", vcf_file);
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
    int mol_type = (*Molecule)[i].Type;
    for (int j = 0; j < (*MoleculeType)[mol_type].nBeads; j++) {
      int type = (*Bead)[(*Molecule)[i].Bead[j]].Type;
      if ((*BeadType)[type].Use) {
        (*Molecule)[i].Bead[count] = (*Molecule)[i].Bead[j];
        (*Bead)[(*Molecule)[i].Bead[j]].Molecule = i;
        count++;
      }
    }
    beads_in_mols[mol_type] = count;
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
        int type = (*Molecule)[id].Type;
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
  int beadtype_alloced = (*Counts).TypesOfBeads;
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

  return indexed;
} //}}}

// ReadCoordinates() //{{{
/**
 * Function reading coordinates from .vcf file with indexed timesteps (\ref IndexedCoorFile).
 */
int ReadCoordinates(bool indexed, FILE *vcf_file, Counts Counts, int *Index, Bead **Bead, char **stuff) {

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

  if (indexed) { // indexed timestep
    // read data //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      fgets(line, sizeof(line), vcf_file);

      // split the line into array //{{{
      char *split[5];
      split[0] = strtok(line, " \t");
      int j = 0;
      while (split[j] != NULL && j < 4) {
        split[++j] = strtok(NULL, " \t");
      } //}}}

      // error - less then four whitespace-separated strings //{{{
      if (j < 4) {
        return i+1;
      } //}}}

      // test if split[0] is integer //{{{
      for (int j = 0; j < strlen(split[0]); j++) {
        if (split[0][j] < '0' || split[0][j] > '9') {
          return i+1;
        }
      } //}}}
      int index = atoi(split[0]);

      // test if split[1-3] are doubles //{{{
      // first two coordinates
      for (int j = 1; j < 4; j++) {
        // first character can be '-' (but must be longer) or a number
        if ((split[j][0] < '0' || split[j][0] > '9') &&
            split[j][0] != '-') {
          return i+1;
        } else if (split[j][0] == '-' && strlen(split[j]) == 1) {
          return i+1;
        }
        // other characters can be numbers or decimal point or newline (last character of 4th split)
        // they can also be 'e' or '-' in case of number format 1.0e-1
        for (int k = 1; k < strlen(split[j]); k++) {
          if ((split[j][k] < '0' || split[j][k] > '9') &&
              split[j][k] != '.' &&
              split[j][k] != '\n' &&
              split[j][k] != 'e' &&
              split[j][k] != '-' ) {
            putc('\n', stderr);
            return i+1;
          }
        }
      } //}}}

      // bead coordinates
      (*Bead)[Index[index]].Position.x = atof(split[1]);
      (*Bead)[Index[index]].Position.y = atof(split[2]);
      (*Bead)[Index[index]].Position.z = atof(split[3]);
    } //}}}
  } else { // ordered timestep
    // read data //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      fgets(line, sizeof(line), vcf_file);

      // split the line into array //{{{
      char *split[4];
      split[0] = strtok(line, " \t");
      int j = 0;
      while (split[j] != NULL && j < 3) {
        split[++j] = strtok(NULL, " \t");
      } //}}}

      // error - less than 3 whitespace-separated strings //{{{
      if (j < 3) {
        return i+1;
      } //}}}

      // test if split[0-2] are doubles //{{{
      for (int j = 0; j < 3; j++) {
        // first character can be '-' (but must be longer) or a number
        if ((split[j][0] < '0' || split[j][0] > '9') &&
            split[j][0] != '-') {
          return i+1;
        } else if (split[j][0] == '-' && strlen(split[j]) == 1) {
          return i+1;
        }
        // other characters can be numbers or decimal point or newline (last character of 4th split)
        for (int k = 1; k < strlen(split[j]); k++) {
          if ((split[j][k] < '0' || split[j][k] > '9') &&
              split[j][k] != '.' &&
              split[j][k] != '\n') { // last split ends with newline
            return i+1;
          }
        }
      } //}}}

      // bead coordinates
      (*Bead)[i].Position.x = atof(split[0]);
      (*Bead)[i].Position.y = atof(split[1]);
      (*Bead)[i].Position.z = atof(split[2]);
    } //}}}
  }

  return 0;
} //}}}

// SkipCoor() //{{{
/**
 * Function to skip one timestep in coordinates file. It works with both
 * indexed and ordered vcf files.
 */
bool SkipCoor(FILE *vcf_file, Counts Counts, char **stuff) {

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
bool ReadAggregates(FILE *agg_file, Counts *Counts, Aggregate **Aggregate,
                    BeadType *BeadType, Bead **Bead,
                    MoleculeType *MoleculeType, Molecule **Molecule, int *Index) {

  bool error = false;

  // is there a Step? i.e., isn't this the line 'Last Step'?
  if (getc(agg_file) == 'L') {
    error = true;
  } else {
    // initialize array of number of aggregates per bead //{{{
    for (int i = 0; i < (*Counts).Beads; i++) {
      (*Bead)[i].nAggregates = 0;
    } //}}}

    // skip 'Step: #' line //{{{
    while (getc(agg_file) != '\n')
      ; //}}}

    fscanf(agg_file, "%d", &(*Counts).Aggregates);

    // skip rest of the line and blank line //{{{
    while (getc(agg_file) != '\n')
      ;
    while (getc(agg_file) != '\n')
      ; //}}}

    // go through all aggregates
    for (int i = 0; i < (*Counts).Aggregates; i++) {
      // read molecules in Aggregate 'i' //{{{
      fscanf(agg_file, "%d :", &(*Aggregate)[i].nMolecules);
      for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
        int mol;
        fscanf(agg_file, "%d", &mol);
        mol--; // in agg file the numbers correspond to vmd

        (*Aggregate)[i].Molecule[j] = mol;
        (*Molecule)[mol].Aggregate = i;
      }

      while (getc(agg_file) != '\n')
       ; //}}}

      // read monomeric beads in Aggregate 'i' //{{{
      int count;
      fscanf(agg_file, "%d :", &count);
      (*Aggregate)[i].nMonomers = 0; // their number will be counted according to which beads are in vcf
      for (int j = 0; j < count; j++) {
        int vsf_id;
        fscanf(agg_file, "%d", &vsf_id); // monomer index in vsf file
        int id = Index[vsf_id]; // monomer index in Bead structure
        if (id > -1) {
          int beads = (*Aggregate)[i].nMonomers++;
          (*Aggregate)[i].Monomer[beads] = id;
          (*Bead)[id].nAggregates++;
          (*Bead)[id].Aggregate = realloc((*Bead)[id].Aggregate, (*Bead)[id].nAggregates*sizeof(int));
          (*Bead)[id].Aggregate[(*Bead)[id].nAggregates-1] = i;
        }
      }

      while (getc(agg_file) != '\n')
       ; //}}}
    }

    // skip blank line at the end of every entry //{{{
    while (getc(agg_file) != '\n')
      ; //}}}

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
void ReadFieldPbc(char *field, Vector *BoxLength) {
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
    fprintf(stderr, "\nError: %s - first line must start with box size (i.e., three positive numbers)\n", field);
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
void ReadFieldBeadType(char *field, Counts *Counts, BeadType **BeadType, Bead **Bead) {

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
        fprintf(stderr, "\nError: %s - missing number of species\n", field);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    fprintf(stderr, "\nError: %s - missing 'species' keyword\n\n", field);
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
      fprintf(stderr, "\nError: %s - no blank spaces permitted in the species section\n\n", field);
      exit(1);
    } else if (words < 4 ||                  // 2)
               !IsPosDouble(split[1]) ||     // 3)
               !IsDouble(split[2]) ||        // 4)
               !IsInteger(split[3])) {       // 5)
      fprintf(stderr, "\nError: %s - wrong species line\n", field);
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
void ReadFieldMolecules(char *field, Counts *Counts,
                        BeadType **BeadType, Bead **Bead,
                        MoleculeType **MoleculeType, Molecule **Molecule) {

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
        fprintf(stderr, "\nError: %s - missing number of molecule types\n", field);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfMolecules = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    fprintf(stderr, "\nError: %s - missing 'molecule' line\n\n", field);
    exit(1);
  } //}}}

  *MoleculeType = calloc((*Counts).TypesOfMolecules,sizeof(struct MoleculeType));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Number = 0;
    (*MoleculeType)[i].nBeads = 0;
    (*MoleculeType)[i].nBonds = 0;
    (*MoleculeType)[i].BType = calloc((*Counts).TypesOfBeads, sizeof(int));
  }

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
    fprintf(stderr, "\nError: %s - missing 'finish' keyword\n\n", field);
    exit(1);
  }
  fsetpos(fr, &position); //}}}

  // read info about molecule types //{{{
  fgetpos(fr, &position); // save pointer position
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // read i-th molecule name //{{{
    fgets(line, sizeof(line), fr);
    SplitLine(split, line, delim);
    if (split[0][0] == '\0') {
      fprintf(stderr, "\nError: %s - missing molecule name\n\n", field);
      exit(1);
    }
    strcpy((*MoleculeType)[i].Name, split[0]); //}}}
    // read number of the molecules //{{{
    fgets(line, sizeof(line), fr);
    int words = SplitLine(split, line, delim);
    if (strcasecmp(split[0], "nummols") != 0) {
      fprintf(stderr, "\nError: %s - missing 'nummols' keyword\n", field);
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      fprintf(stderr, "\nError: %s - 'nummols' must be followed by a non-negative integer\n", field);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].Number = atoi(split[1]); //}}}
    // read number of beads in the molecule //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "beads", 4) != 0) {
      fprintf(stderr, "\nError: %s - missing 'beads' keyword\n", field);
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      fprintf(stderr, "\nError: %s - 'beads' must be followed by a non-negative integer\n", field);
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
        fprintf(stderr, "\nError: %s - wrong bead line in molecule %s\n", field, (*MoleculeType)[i].Name);
        ErrorPrintLine(split, words);
        exit(1);
      }
      int type = FindBeadType(split[0], *Counts, *BeadType);
      if (type == -1) {
        fprintf(stderr, "\nError: %s - non-existent bead name '%s' in molecule %s\n", field, split[0], (*MoleculeType)[i].Name);
        ErrorBeadType(*Counts, *BeadType);
        exit(1);
      }
      (*MoleculeType)[i].Bead[j] = type;
    } //}}}
    // read number of bonds in the molecule //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "bonds", 4) != 0) {
      fprintf(stderr, "\nError: %s - missing 'bonds' keyword\n", field);
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      fprintf(stderr, "\nError: %s - 'bonds' must be followed by a non-negative integer\n", field);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].nBonds = atoi(split[1]);
    (*MoleculeType)[i].Bond = calloc((*MoleculeType)[i].nBonds, sizeof(int *));
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      (*MoleculeType)[i].Bond[j] = calloc(2, sizeof(int));
    } //}}}
    // read bond info //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line, delim);
      // bead line must be '<type> <bead id> <bead id> <parameters - ignored>'
      if (words < 3 ||
          !IsInteger(split[1]) || atoi(split[1]) == 0 || // bead ids must start from 1
          !IsInteger(split[2]) || atoi(split[2]) == 0) {
        fprintf(stderr, "\nError: %s - wrong bond line in molecule %s\n", field, (*MoleculeType)[i].Name);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*MoleculeType)[i].Bond[j][0] = atoi(split[1]) - 1;
      (*MoleculeType)[i].Bond[j][1] = atoi(split[2]) - 1;
    } //}}}
    (*Counts).Bonded += (*MoleculeType)[i].Number * (*MoleculeType)[i].nBeads;
    (*Counts).Molecules += (*MoleculeType)[i].Number;
    // skip till 'finish' //{{{
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  }
  // return the file pointer to the beginning of the molecules section
  fsetpos(fr, &position); //}}}

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

  // read coordinates of bonded beads //{{{
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
    // skip till 'finish' //{{{
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  } //}}}

  fclose(fr);
} //}}}

// ReadField() //{{{
/**
 * Function reading the FIELD-like file; it completely fills provided structs,
 * overriding any possible data in there. If the FIELD-like file is a source of
 * some additional structure information, new structs must be used, and the
 * data copied from there afterwards.
 */
void ReadField(char *field, Vector *BoxLength, Counts *Counts,
               BeadType **BeadType, Bead **Bead, int **Index,
               MoleculeType **MoleculeType, Molecule **Molecule) {

  // read pbc only if required
  if (BoxLength != NULL) {
    ReadFieldPbc(field, &(*BoxLength));
  }
  ReadFieldBeadType(field, &(*Counts), &(*BeadType), &(*Bead));
  ReadFieldMolecules(field, &(*Counts), &(*BeadType), &(*Bead), &(*MoleculeType), &(*Molecule));

  // allocate Bead[].Aggregate array - needed only to free() //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(10, sizeof(int));
  } //}}}

  // fill MoleculeType[].BType array //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
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
} //}}}
