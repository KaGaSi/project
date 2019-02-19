#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "AnalysisTools.h"
#include "Errors.h"

// ReadFIELD() - auxiliary for ReadStructure() //{{{
/*
 * Function reading information about bead types (mass and charge) from
 * 'species' lines in FIELD. It only reads information about bead types with
 * Use=true flag.
 */
void ReadFIELD(char *bonds_file, Counts *Counts,
               BeadType **BeadType, MoleculeType **MoleculeType) {

  // open FIELD file //{{{
  FILE *field;
  if ((field = fopen("FIELD", "r")) == NULL) {
    ErrorFileOpen("FIELD", 'r');
    exit(1);
  } //}}}

  // skip lines in FIELD till 'species' //{{{
  char *split;
  do {
    // get whole line - max 1000 chars
    char line[1024];
    fgets(line, 1024, field);

    // first string of the line
    split = strtok(line, " \t");

  } while (strcmp(split, "species") != 0 &&
           strcmp(split, "SPECIES") != 0 &&
           strcmp(split, "Species") != 0); //}}}

  // read number of bead types //{{{
  split = strtok(NULL, " \n");
  int bead_types = atoi(split); //}}}

  // read all bead types from FIELD //{{{
  for (int i = 0; i < bead_types; i++) {

    // get whole line - max 1000 chars //{{{
    char line[1024];
    fgets(line, 1024, field); //}}}

    // read bead type name and test it is against data from vsf and vcf //{{{
    split = strtok(line, " \t");
    int type;
    if ((type = FindBeadType(split, *Counts, *BeadType)) == -1) {
      continue;
    } //}}}

    // read bead type mass //{{{
    if ((split = strtok(NULL, " \t")) == NULL) {
      fprintf(stderr, "Error: cannot read bead mass from FIELD\n");
      exit(1);
    }
    (*BeadType)[type].Mass = atof(split); //}}}

    // read bead type charge //{{{
    if ((split = strtok(NULL, " \t")) == NULL) {
      fprintf(stderr, "Error: cannot read bead charge from FIELD\n");
      exit(1);
    }
    (*BeadType)[type].Charge = atof(split); //}}}
  } //}}}
  fclose(field);
} //}}}

// CommonHelp() //{{{
/**
 * Function to print help for common options, either for `-h` help option
 * or program error.
 */
void CommonHelp(bool error) {
  if (error) {
    fprintf(stderr, "   <standard options>\n");
    fprintf(stderr, "      -i <name>      use input .vsf file different from traject.vsf\n");
//  fprintf(stderr, "      -b <name>      file containing bond alternatives to FIELD\n");
    fprintf(stderr, "      -v             verbose output\n");
    fprintf(stderr, "      -V             more verbose output\n");
    fprintf(stderr, "      -s             no output (overrides verbose options)\n");
    fprintf(stderr, "      -h             print this help and exit\n");
    fprintf(stderr, "      --script       do not reprint line (useful when output goes to file)\n");
  } else {
    fprintf(stdout, "   <standard options>\n");
    fprintf(stdout, "      -i <name>      use input .vsf file different from traject.vsf\n");
//  fprintf(stdout, "      -b <name>      file containing bond alternatives to FIELD\n");
    fprintf(stdout, "      -v             verbose output\n");
    fprintf(stdout, "      -V             more verbose output\n");
    fprintf(stdout, "      -s             no output (overrides verbose options)\n");
    fprintf(stdout, "      -h             print this help and exit\n");
    fprintf(stdout, "      --script       do not reprint line (useful when output goes to file)\n");
  }
} //}}}

// CommonOptions() //{{{
/**
 * Function for options common to most of the utilities.
 */
bool CommonOptions(int argc, char **argv, char **vsf_file,char **bonds_file,
                   bool *verbose, bool *verbose2, bool *silent, bool *script) {

  bool error = false;

  // -i <name> option - filename of input structure file //{{{
  (*vsf_file)[0] = '\0'; // check if -i option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");

        error = true;
        break;
      }

      // check if .vsf ending is present
      char *vsf = strrchr(argv[i+1], '.');
      if (!vsf || strcmp(vsf, ".vsf")) {
        fprintf(stderr, "'-i' arguments does not have .vsf ending!\n");
        error = true;
      }

      strcpy(*vsf_file, argv[i+1]);
    }
  }

  // -i option is not used
  if ((*vsf_file)[0] == '\0') {
    strcpy(*vsf_file, "traject.vsf");
  } //}}}

  // -b <name> option - filename of input bond file //{{{
  (*bonds_file)[0] = '\0'; // check if -b option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-b") == 0) {

      // wrong argument to -b option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-b' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n\n");

        error = true;
        break;
      }

      strcpy(*bonds_file, argv[i+1]);
    }
  } //}}}

  // -v option - verbose output //{{{
  *verbose = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      *verbose = true;

      break;
    }
  } //}}}

  // -V option - verbose output with comments from input .vcf file //{{{
  *verbose2 = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-V") == 0) {
      *verbose = true;
      *verbose2 = true;

      break;
    }
  } //}}}

  // -s option - silent mode //{{{
  *silent = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-s") == 0) {
      *verbose = false;
      *verbose2 = false;
      *silent = true;

      break;
    }
  } //}}}

  // --script  option - meant for when output is routed to file, so don't use flush & \r //{{{
  *script = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--script") == 0) {
      *script = true;

      break;
    }
  } //}}}

  return error;
} //}}}

// VerboseOutput() //{{{
/**
 * Function providing standard verbose output (for cases when verbose
 * option is used). It prints most of the information about used system.
 */
void VerboseOutput(bool Verbose2, char *input_vcf, char *bonds_file, Counts Counts,
                   BeadType *BeadType, Bead *Bead,
                   MoleculeType *MoleculeType, Molecule *Molecule) {

//if (input_vcf[0] != '\0')
  printf("Counts.{");
  printf("TypesOfBeads =%3d, ", Counts.TypesOfBeads);
  printf("Bonded =%7d, ", Counts.Bonded);
  printf("Unbonded =%7d, ", Counts.Unbonded);
  printf("Beads = %7d, ", Counts.Beads);
  printf("BeadsInVsf = %7d, ", Counts.BeadsInVsf);
  printf("TypesOfMolecules =%3d, ", Counts.TypesOfMolecules);
  printf("Molecules =%4d}\n", Counts.Molecules);
  printf("\ntotal number of beads: %d\n\n", Counts.Bonded+Counts.Unbonded);

  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    printf("BeadType[%2d].{", i);
    printf("Name =%10s, ", BeadType[i].Name);
    printf("Number =%7d, ", BeadType[i].Number);
    printf("Charge =%6.2f, ", BeadType[i].Charge);
    printf("Mass =%5.2f, ", BeadType[i].Mass);
    printf("Use = %3s, ", BeadType[i].Use? "Yes":"No");
    printf("Write = %3s}\n", BeadType[i].Write? "Yes":"No");
  }
  putchar('\n');

  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf("MoleculeType[%d].{", i);
    printf("Name =%10s", MoleculeType[i].Name);
    printf(", Number =%4d", MoleculeType[i].Number);
    printf(", nBeads =%3d", MoleculeType[i].nBeads);
    printf(", nBonds =%3d", MoleculeType[i].nBonds);

    printf(", nBTypes =%2d, BType{", MoleculeType[i].nBTypes);
    printf("%8s", BeadType[MoleculeType[i].BType[0]].Name);
    for (int j = 1; j < MoleculeType[i].nBTypes; j++) {
      printf(",%8s", BeadType[MoleculeType[i].BType[j]].Name);
    }

    printf("}, Mass =%7.2f", MoleculeType[i].Mass);

    if (bonds_file[0] == '\0') { // all bonds taken from vsf file
//    printf(", Bonds from vsf,");
    } else {
      // go through bond file to find out if molecule type 'i' is there
      FILE *bond;
      if ((bond = fopen(bonds_file, "r")) == NULL) {
        ErrorFileOpen(bonds_file, 'r');
        exit(1);
      }

      int test;
      char str[32];
      while ((test = getc(bond)) != EOF) {
        ungetc(test, bond);

        if ((fscanf(bond, "%s %d", str, &test)) != 2) {
          fprintf(stderr, "Cannot read string or number of bonds from %s with '-v' option!\n", bonds_file);
          exit(1);
        }

        if (strcmp(str, MoleculeType[i].Name) == 0) {
          printf(", Bonds from '%s',", bonds_file);
          break;
        }

        while (getc(bond) != '\n')
          ;
      }

      // if not in bonds_file, then bonds taken from FIELD
      if (test == EOF) {
        printf(", Bonds from 'FIELD',");
      }

      fclose(bond);
    }
    printf(" Use = %3s,", MoleculeType[i].Use? "Yes":"No");
    printf(" Write = %3s}\n", MoleculeType[i].Write? "Yes":"No");
  }

  printf("Bonds:\n");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf("MoleculeType[%d].nBonds = %d:", i, MoleculeType[i].nBonds);
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      printf(" %d-%d", MoleculeType[i].Bond[j][0], MoleculeType[i].Bond[j][1]);
      if (j < (MoleculeType[i].nBonds-1)) {
        putchar(',');
      }
    }
    putchar('\n');
  }

  // print bead ids of all molecules if '-V' option is used
  for (int i = 0; Verbose2 && i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    printf("Molecule %3d (%s):\n", i+1, MoleculeType[type].Name);
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      printf(" %d (%d)", Molecule[i].Bead[j], Bead[Molecule[i].Bead[j]].Index);
    }
    printf("\n");
  }
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
bool ReadStructure(char *vsf_file, char *vcf_file, char *bonds_file, Counts
    *Counts, BeadType **BeadType, Bead **Bead, MoleculeType **MoleculeType,
    Molecule **Molecule) {

  FILE *vsf;

  // zeroize stuff //{{{
  (*Counts).TypesOfBeads = 0;
  (*Counts).TypesOfMolecules = 0;
  (*Counts).Molecules = 0; //}}}

  // initial allocations (realloced later)
  *BeadType = calloc(1,sizeof(struct BeadType));
  *MoleculeType = calloc(1,sizeof(struct MoleculeType));

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
  char line[1024];
  while(fgets(line, sizeof(line), vsf)) {
    // trim trailing whitespace in line //{{{
    int length = strlen(line);
    // last string character needs to be '\0'
    while (length > 1 &&
           (line[length-1] == ' ' ||
            line[length-1] == '\n' ||
            line[length-1] == '\t')) {
      line[length-1] = '\0';
      length--;
    }
    // for lines containing only whitespace
    if (length == 1) {
      atom_lines++;
      continue;
    } //}}}

    // split the line into array //{{{
    char *split[30];
    split[0] = strtok(line, " \t");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t");
    }
    int strings = i; //}}}
    if (split[0][0] == 'a') {
      atom_lines++;

      // error - odd number of values: a(tom) lines are composed of 'keyword <value>' pairs //{{{
      if ((strings % 2) == 1) { // because it's split[0] ... split[2n-1]
        fprintf(stderr, "Error: vsf - odd number of strings (%d) in the following a(tom) line:\n|", strings);
        for (int i = 0; i < strings; i++) {
          fprintf(stderr, "%s ", split[i]);
        }
        fprintf(stderr, "|\n");
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
        for (int j = 2; j < strings; j += 2) {
          if (strncmp(split[j], "name", 1) == 0) {
            // copy new name to BeadType[].Name
            strcpy((*BeadType)[type_default].Name, split[j+1]);
            (*BeadType)[type_default].Charge = 1000;
            (*BeadType)[type_default].Mass = 1000;
            name = true;
          } else if (strcmp(split[j], "charge") == 0 || strcmp(split[j], "q") == 0) {
            charge = atof(split[j+1]);
          } else if (strncmp(split[j], "mass", 1) == 0) {
            mass = atof(split[j+1]);
          }
        }
        if (charge != 1000) {
          (*BeadType)[type_default].Charge = charge;
        }
        if (mass != 1000) {
          (*BeadType)[type_default].Mass = mass;
        }
        continue; // skip the rest of the line
      } else if (atoi(split[1]) > max_bead) {
        max_bead = atoi(split[1]);
      } //}}}

      // check that line contains either both 'resid' and 'resname' or neither //{{{
      int test = 0;
      for (int i = 2; i < strings && test < 2; i+= 2) {
        if (strcmp("resname", split[i]) == 0) {
          test++;
        } else if (strcmp("resid", split[i]) == 0) {
          test++;
        }
      }
      if (test == 1) {
        fprintf(stderr, "Error: vsf - the following a(tom) line contains only one of 'resid' and 'resname' keywords:\n");
        for (int i = 0; i < strings; i++) {
          fprintf(stderr, "%s ", split[i]);
        }
        fprintf(stderr, "\n");
        exit(1);
      } //}}}

      // go through the rest of the line to find bead name
      // by twos - first is always keyword, second is always value
      int type_qm = -1; // for possible mass/charge in vsf
      for (int i = 2; i < strings; i += 2) {
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
          (*BeadType)[type].Number++;
          //}}}
        } else if (strcmp("resname", split[i]) == 0) { // molecule type name //{{{
          // if the molecule name doesn't exist, add it to the structures
          int type;
          if ((type = FindMoleculeType(split[i+1], *Counts, *MoleculeType)) == -1) {
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
          }
          //}}}
        } else if (strcmp("resid", split[i]) == 0) { // molecule id //{{{
          if (atoi(split[i+1]) > max_mol) {
            max_mol = atoi(split[i+1]);
          }
        //}}}
        } else if (strcmp("charge", split[i]) == 0 || strcmp("q", split[i]) == 0) { //{{{
          charge = atof(split[i+1]);
        //}}}
        } else if (strncmp("mass", split[i], 1) == 0) { //{{{
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
        fprintf(stderr, "\nError: vsf - the following a(tom) line does not contain 'name':\n");
        for (int i = 0; i < strings; i++) {
          fprintf(stderr, "%s ", split[i]);
        }
        fprintf(stderr, "\n");
        exit(1);
      } //}}}
    } else if (split[0][0] == '#') {
      // count the line if comment or blank
      atom_lines++;
    } else {
      break;
    }
  }
  fclose(vsf); //}}}

  (*Counts).Molecules = max_mol; // mol ids start from 1 in vsf
  (*Counts).BeadsInVsf = max_bead + 1; // bead ids start from 0 in vsf

  // allocate Bead and Molecule structures
  *Bead = calloc((*Counts).BeadsInVsf, sizeof(struct Bead));
  int mol_alloced = (*Counts).Molecules;
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

  // go through atom lines
  for (int count = 0; count < atom_lines; count++) {
    // read line
    fgets(line, sizeof(line), vsf);
    // trim trailing whitespace in line //{{{
    int length = strlen(line);
    // last string character needs to be '\0'
    while (length > 1 &&
           (line[length-1] == ' ' ||
            line[length-1] == '\n' ||
            line[length-1] == '\t')) {
      line[length-1] = '\0';
      length--;
    }
    // for lines containing only whitespace
    if (length == 1) {
      continue;
    }//}}}

    // split the line into array //{{{
    char *split[30];
    split[0] = strtok(line, " \t");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t");
    }
    int strings = i; //}}}

    // go through the line //{{{
    int bead_id = -1, bead_type = -1, mol_id = -1, mol_type = -1;
    if (split[0][0] == 'a' && strcmp("default", split[1]) != 0) { // non-default a(tom) line
      bead_id = atoi(split[1]);
      // by twos - first is always keyword, second is always value
      for (int i = 2; i < strings; i += 2) {
        if (strncmp("name", split[i], 1) == 0) { // bead name
          bead_type = FindBeadType(split[i+1], *Counts, *BeadType);
        } else if (strcmp("resname", split[i]) == 0) { // molecule name
          mol_type = FindMoleculeType(split[i+1], *Counts, *MoleculeType);
        } else if (strcmp("resid", split[i]) == 0) { // molecule id
          mol_id = atoi(split[i+1]) - 1; // mol ids start with 1 in vsf
        }
      }

      // assign values to Bead & Molecule
      (*Bead)[bead_id].Type = bead_type;
      (*Bead)[bead_id].Index = bead_id;
      (*Bead)[bead_id].Molecule = mol_id;
      if (mol_id > -1) {
        (*Molecule)[mol_id].Type = mol_type;
        (*MoleculeType)[mol_type].nBeads++;
      }
    } //}}}
  } //}}}

  // assign 'type_default' to default beads //{{{
  int count = 0;
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if ((*Bead)[i].Type == -1) {
      (*Bead)[i].Type = type_default;
      (*Bead)[i].Index = i;
      (*Bead)[i].Molecule = -1; // default beads aren't in molecules
      count++;
    } else if ((*Bead)[i].Type == type_default) { // default type beads explicitly specified by 'atom' line
      count++;
    }
  }
  (*BeadType)[type_default].Number = count; //}}}

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
    // split the line into array - no need to trim trailing white space (the number of columns is known) //{{{
    char *split[30];
    split[0] = strtok(line, " \t:");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t:");
    } //}}}

    // b(ond) lines in the form 'bond <int>: <int>'
    if (split[0][0] == 'b') {
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
    }
  } //}}}

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

  // restore pointer position in vsf file
  fsetpos(vsf, &pos);

  // fourth, go again through the bonds and assign ids of bonded beads //{{{
  // initialize helper arrays //{{{
  int *count_bonds;
  count_bonds = calloc((*Counts).TypesOfMolecules, sizeof(int));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  } //}}}
  while (fgets(line, sizeof(line), vsf)) {
    // split the line into array //{{{
    char *split[30];
    split[0] = strtok(line, " \t:");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t:");
    } //}}}
    if (split[0][0] == 'b') { // b(ond) lines
      int mol_id = (*Bead)[atoi(split[1])].Molecule;
      int mol_type = (*Molecule)[mol_id].Type;
      // is this mol_type already in use?
      if (bonds[mol_type] == -1) {
        bonds[mol_type] = mol_id;
      }
      // increment number of bonds if correct molecule
      if (bonds[mol_type] == mol_id) {
        (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][0] = atoi(split[1]);
        (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][1] = atoi(split[2]);
        count_bonds[mol_type]++;
      }
    }
  }
  free(count_bonds);
  fclose(vsf); //}}}

  // minimize bead ids in bonds //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    int min = 1000000;
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      // swap beads if first bond id is higher than the second
      if ((*MoleculeType)[i].Bond[j][0] > (*MoleculeType)[i].Bond[j][1]) {
        int swap = (*MoleculeType)[i].Bond[j][0];
        (*MoleculeType)[i].Bond[j][0] = (*MoleculeType)[i].Bond[j][1];
        (*MoleculeType)[i].Bond[j][1] = swap;
      }
      // find minimum id
      if ((*MoleculeType)[i].Bond[j][0] < min) {
        min = (*MoleculeType)[i].Bond[j][0];
      }
      if ((*MoleculeType)[i].Bond[j][1] < min) {
        min = (*MoleculeType)[i].Bond[j][1];
      }
    }
    // subtract minimum id from all bonded beads' ids
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      (*MoleculeType)[i].Bond[j][0] -= min;
      (*MoleculeType)[i].Bond[j][1] -= min;
    }
  } //}}}

  // bubble sort the bond arrays //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // go through all bonds in molecule type 'i'
    for (int j = 0; j < ((*MoleculeType)[i].nBonds-1); j++) {
      int swap = 0;
      for (int k = 0; k < ((*MoleculeType)[i].nBonds-j-1); k++) {
        // swap bonds if the first beads are in wrong order
        if ((*MoleculeType)[i].Bond[k][0] > (*MoleculeType)[i].Bond[k+1][0]) {
          swap = (*MoleculeType)[i].Bond[k][0];
          (*MoleculeType)[i].Bond[k][0] = (*MoleculeType)[i].Bond[k+1][0];
          (*MoleculeType)[i].Bond[k+1][0] = swap;

          swap  = (*MoleculeType)[i].Bond[k][1];
          (*MoleculeType)[i].Bond[k][1] = (*MoleculeType)[i].Bond[k+1][1];
          (*MoleculeType)[i].Bond[k+1][1] = swap;
        // swap bonds if the first beads are the same,
        // but second ones are in wrong order
        } else if ((*MoleculeType)[i].Bond[k][0] == (*MoleculeType)[i].Bond[k+1][0] &&
                   (*MoleculeType)[i].Bond[k][1] > (*MoleculeType)[i].Bond[k+1][1]) {
          swap = (*MoleculeType)[i].Bond[k][0];
          (*MoleculeType)[i].Bond[k][0] = (*MoleculeType)[i].Bond[k+1][0];
          (*MoleculeType)[i].Bond[k+1][0] = swap;

          swap  = (*MoleculeType)[i].Bond[k][1];
          (*MoleculeType)[i].Bond[k][1] = (*MoleculeType)[i].Bond[k+1][1];
          (*MoleculeType)[i].Bond[k+1][1] = swap;
        }
      }

      // if no swap was made, the array is sorted
      if (swap == 0)
        break;
    }
  } //}}}

  // fifth, go through vsf and assign bead ids to molecules //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}

  // go through atom lines
  int *beads; // helper array
  beads = calloc((*Counts).BeadsInVsf, sizeof(int));
  for (int count = 0; count < atom_lines; count++) {
    // read line
    fgets(line, sizeof(line), vsf);
    // trim trailing whitespace in line //{{{
    int length = strlen(line);
    // last string character needs to be '\0'
    while (length > 1 &&
           (line[length-1] == ' ' ||
            line[length-1] == '\n' ||
            line[length-1] == '\t')) {
      line[length-1] = '\0';
      length--;
    }
    // for lines containing only whitespace
    if (length == 1) {
      continue;
    }//}}}

    // split the line into array //{{{
    char *split[30];
    split[0] = strtok(line, " \t");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t");
    }
    int strings = i; //}}}

    // go through the line //{{{
    int bead_id = -1, mol_id;
    if (split[0][0] == 'a' && strcmp("default", split[1]) != 0) { // non-default a(tom) line
      bead_id = atoi(split[1]);
      // by twos - first is always keyword, second is always value
      for (int i = 2; i < strings; i += 2) {
        if (strcmp("resid", split[i]) == 0) { // molecule id
          mol_id = atoi(split[i+1]) - 1; // mol ids start with 1 in vsf
          (*Molecule)[mol_id].Bead[beads[mol_id]++] = bead_id;
          break;
        }
      }
    } //}}}
  }
  free(beads);
  fclose(vsf); //}}}

  // bubble sort the bead arrays //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    // go through all besds in the molecule
    int type = (*Molecule)[i].Type;
    for (int j = 0; j < ((*MoleculeType)[type].nBeads-1); j++) {
      int swap = -1;
      for (int k = 0; k < ((*MoleculeType)[type].nBeads-j-1); k++) {
        if ((*Molecule)[i].Bead[k] > (*Molecule)[i].Bead[k+1]) {
          swap = (*Molecule)[i].Bead[k];
          (*Molecule)[i].Bead[k] = (*Molecule)[i].Bead[k+1];
          (*Molecule)[i].Bead[k+1] = swap;
        }
      }

      // if no swap was made, the array is sorted
      if (swap == -1)
        break;
    }
  } //}}}

  // fill BTypes array //{{{
  // allocate arrays //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].BType = calloc((*Counts).TypesOfBeads, sizeof(int));
  } //}}}
  int tmp = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nBTypes = 0;
    // go through all beads of the first molecule of the given molecule type
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int id = (*Molecule)[tmp].Bead[j];
      // test if bead id's type is present in BType array //{{{
      bool in_mol = false;
      for (int k = 0; k < (*MoleculeType)[i].nBTypes; k++) {
        if ((*MoleculeType)[i].BType[k] == (*Bead)[id].Type) {
          in_mol = true;
          break;
        }
      } //}}}
      // if bead id is of a type not yet present in BType array, add it there
      if (!in_mol) {
        (*MoleculeType)[i].BType[(*MoleculeType)[i].nBTypes] = (*Bead)[id].Type;

        (*MoleculeType)[i].nBTypes++;
      }
    }
    // count total number of molecules
    tmp += (*MoleculeType)[i].Number;
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
    char line[1024], line2[1024], *split[32], str[32];
    str[0] = '\0';
    do {
      fgets(line, sizeof(line), vcf);
      strcpy(line2, line);
      split[0] = strtok(line, " \t");

      // 't(imestep)' line
      if (split[0][0] == 't' || split[0][0] == 'T') {
        split[1] = strtok(NULL, " \t");
        if (split[1][0] != 'o' && split[1][0] != 'O' &&
            split[1][0] != 'i' && split[1][0] != 'I') {
          fprintf(stderr, "\nError: %s - no 'i(ndexed)' or 'o(rdered)' keyword after 't(imestep)' keyword\n", vcf_file);
          exit(1);
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

      // ignore the rest of the line //{{{
      while(getc(vcf) != '\n')
        ; //}}}

      // read data //{{{
      // the first coordinate line
      fgets(line, sizeof(line), vcf);
      // trim trailing whitespace in line //{{{
      int length = strlen(line);
      // last string character needs to be '\0'
      while (length > 1 &&
             (line[length-1] == ' ' ||
              line[length-1] == '\n' ||
              line[length-1] == '\t')) {
        line[length-1] = '\0';
        length--;
      }
      // for lines containing only whitespace
      if (length == 1) {
        fprintf(stderr, "Error: %s - blank line instead of a first coordinate line\n", vcf_file);
        exit(1);
      }//}}}
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
        // trim trailing whitespace in line //{{{
        int length = strlen(line);
        // last string character needs to be '\0'
        while (length > 1 &&
               (line[length-1] == ' ' ||
                line[length-1] == '\n' ||
                line[length-1] == '\t')) {
          line[length-1] = '\0';
          length--;
        } //}}}
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
  count = 0;
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

  // free name array //{{{
  for (int i = 0; i < moltype_alloced; i++) {
    free(name[i]);
  }
  free(name); //}}}

  // remove unused beads from bonds in molecules //{{{
  int mol_type[(*Counts).TypesOfMolecules]; // helper array //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    mol_type[i] = -1;
  } //}}}
  // go through all molecules
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    // is this molecule type done?
    if (mol_type[mtype] == -1) {
      count = 0;
      // go through all bonds //{{{
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
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
      // find minimum id in Bond array //{{{
      int min = 1000000;
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        if ((*MoleculeType)[mtype].Bond[j][0] < min) {
          min = (*MoleculeType)[mtype].Bond[j][0];
        }
      } //}}}
      // minimize ids in Bond array //{{{
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        (*MoleculeType)[mtype].Bond[j][0] -= min;
        (*MoleculeType)[mtype].Bond[j][1] -= min;
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

  int *index;
  index = malloc((*Counts).BeadsInVsf*sizeof(int));
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    index[i] = -1;
  }

  // remove unused beads from Bead struct //{{{
  count = 0;
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int type = (*Bead)[i].Type;
    if ((*BeadType)[type].Use) {
      index[i] = count;
      (*Bead)[count].Type = (*Bead)[i].Type;
      (*Bead)[count].Molecule = (*Bead)[i].Molecule;
      (*Bead)[count].Index = (*Bead)[i].Index;
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
  int beadtype_alloced = (*Counts).TypesOfBeads;
  name = calloc(beadtype_alloced, sizeof(char *));
  for (int i = 0; i < beadtype_alloced; i++) {
    name[i] = calloc(16, sizeof(char));
  }
  for (int i = 0; i < beadtype_alloced; i++) {
    strcpy(name[i], (*BeadType)[i].Name);
  } //}}}

  // remove unused bead types from BeadType struct and molecule bonds //{{{
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if ((*BeadType)[i].Use) {
      if (count != i) {
        strcpy((*BeadType)[count].Name, (*BeadType)[i].Name);
        (*BeadType)[count].Number = (*BeadType)[i].Number;
        (*BeadType)[count].Use = (*BeadType)[i].Use;
        for (int j = 0; j < (*Counts).TypesOfMolecules; j++) {
          for (int k = 0; k < (*MoleculeType)[j].nBTypes; k++) {
            if ((*MoleculeType)[j].BType[k] == i) {
              (*MoleculeType)[j].BType[k] = count;
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

  // read bonds - again, unfortunately //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}

  // skip atom lines //{{{
  for (int count = 0; count < atom_lines; count++) {
    fgets(line, sizeof(line), vsf);
  } //}}}

  // save pointer position in the vsf file
  fgetpos(vsf, &pos); //

  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nBonds = 0;
  }
  // go through the bonds section of vsf to find the number of bonds in molecules //{{{
  // clean helper array //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  } //}}}
  while (fgets(line, sizeof(line), vsf)) {
    // split the line into array //{{{
    char *split[30];
    split[0] = strtok(line, " \t:");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t:");
    } //}}}
    if (split[0][0] == 'b') { // b(ond) lines
      // find ids of both beads
      int bead1 = index[atoi(split[1])];
      int bead2 = index[atoi(split[2])];

      // are both beads in a vcf? //{{{
      if (bead1 == -1 || bead2 == -1) {
//      printf("%d:%d %d %d\n", atoi(split[1]), atoi(split[2]), bead1, bead2);
        continue;
      } //}}}

      int mol_id = (*Bead)[bead1].Molecule;
      int mol_type = (*Molecule)[mol_id].Type;
//    printf("%5d(%5d):%5d(%5d) %5d %5d\n", bead1, (*Bead)[bead1].Index, bead2, (*Bead)[bead2].Index, mol_type, mol_id);
      // is this mol_type already in use?
      if (bonds[mol_type] == -1) {
        bonds[mol_type] = mol_id;
      }
//    printf("mol_type=%d\n", mol_type);
      // increment number of bonds if correct molecule
      if (bonds[mol_type] == mol_id) {
        (*MoleculeType)[mol_type].nBonds++;
      }
    }
  } //}}}

  // restore pointer position in vsf file
  fsetpos(vsf, &pos);

  // go again through the bond section and assign ids of bonded beads //{{{
  // initialize helper arrays //{{{
  count_bonds = calloc((*Counts).TypesOfMolecules, sizeof(int));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  } //}}}
  while (fgets(line, sizeof(line), vsf)) {
    // split the line into array //{{{
    char *split[30];
    split[0] = strtok(line, " \t:");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t:");
    } //}}}
    if (split[0][0] == 'b') { // b(ond) lines
      // find ids of both beads
      int bead1 = index[atoi(split[1])];
      int bead2 = index[atoi(split[2])];

      // are both beads in a vcf? //{{{
      if (bead1 == -1 || bead2 == -1) {
        continue;
      } //}}}

      int mol_id = (*Bead)[bead1].Molecule;
      int mol_type = (*Molecule)[mol_id].Type;
      // is this mol_type already in use?
      if (bonds[mol_type] == -1) {
        bonds[mol_type] = mol_id;
      }
      // increment number of bonds if correct molecule
      if (bonds[mol_type] == mol_id) {
        (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][0] = atoi(split[1]);
        (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][1] = atoi(split[2]);
        count_bonds[mol_type]++;
      }
    }
  }
  free(count_bonds);
  fclose(vsf); //}}}

  // minimize bead ids in bonds //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    int min = 1000000;
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      // swap beads if first bond id is higher than the second
      if ((*MoleculeType)[i].Bond[j][0] > (*MoleculeType)[i].Bond[j][1]) {
        int swap = (*MoleculeType)[i].Bond[j][0];
        (*MoleculeType)[i].Bond[j][0] = (*MoleculeType)[i].Bond[j][1];
        (*MoleculeType)[i].Bond[j][1] = swap;
      }
      // find minimum id
      if ((*MoleculeType)[i].Bond[j][0] < min) {
        min = (*MoleculeType)[i].Bond[j][0];
      }
      if ((*MoleculeType)[i].Bond[j][1] < min) {
        min = (*MoleculeType)[i].Bond[j][1];
      }
    }
    // subtract minimum id from all bonded beads' ids
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      (*MoleculeType)[i].Bond[j][0] -= min;
      (*MoleculeType)[i].Bond[j][1] -= min;
    }
  } //}}}

  // bubble sort the bond arrays //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // go through all bonds in molecule type 'i'
    for (int j = 0; j < ((*MoleculeType)[i].nBonds-1); j++) {
      int swap = 0;
      for (int k = 0; k < ((*MoleculeType)[i].nBonds-j-1); k++) {
        // swap bonds if the first beads are in wrong order
        if ((*MoleculeType)[i].Bond[k][0] > (*MoleculeType)[i].Bond[k+1][0]) {
          swap = (*MoleculeType)[i].Bond[k][0];
          (*MoleculeType)[i].Bond[k][0] = (*MoleculeType)[i].Bond[k+1][0];
          (*MoleculeType)[i].Bond[k+1][0] = swap;

          swap  = (*MoleculeType)[i].Bond[k][1];
          (*MoleculeType)[i].Bond[k][1] = (*MoleculeType)[i].Bond[k+1][1];
          (*MoleculeType)[i].Bond[k+1][1] = swap;
        // swap bonds if the first beads are the same,
        // but second ones are in wrong order
        } else if ((*MoleculeType)[i].Bond[k][0] == (*MoleculeType)[i].Bond[k+1][0] &&
                   (*MoleculeType)[i].Bond[k][1] > (*MoleculeType)[i].Bond[k+1][1]) {
          swap = (*MoleculeType)[i].Bond[k][0];
          (*MoleculeType)[i].Bond[k][0] = (*MoleculeType)[i].Bond[k+1][0];
          (*MoleculeType)[i].Bond[k+1][0] = swap;

          swap  = (*MoleculeType)[i].Bond[k][1];
          (*MoleculeType)[i].Bond[k][1] = (*MoleculeType)[i].Bond[k+1][1];
          (*MoleculeType)[i].Bond[k+1][1] = swap;
        }
      }

      // if no swap was made, the array is sorted
      if (swap == 0)
        break;
    }
  } //}}}
  //}}}

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

  // read bead type mass and charge from FIELD file if not in vsf file //{{{
  bool charge_mass = true;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if ((*BeadType)[i].Charge == 1000 || (*BeadType)[i].Mass == 1000) {
      charge_mass = false;
    }
  }
  if (!charge_mass) {
    ReadFIELD(bonds_file, Counts, BeadType, MoleculeType);
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
    // quite a lot memory to be on the save side
    (*Bead)[i].Aggregate = calloc(20, sizeof(int));
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

//// prints for debugging //{{{
//// number of a(tom) lines (including comments and blanks) in vsf
//printf("atom_lines=%d\n", atom_lines);
//// Counts struct
//printf("Counts.{TypesOfMolecules = %d,", (*Counts).TypesOfMolecules);
//printf(" Molecules = %d,", (*Counts).Molecules);
//printf(" TypesOfBeads = %d,", (*Counts).TypesOfBeads);
//printf(" BeadsInVsf = %d,", (*Counts).BeadsInVsf);
//printf(" Beads = %d,", (*Counts).Beads);
//printf(" Unbonded = %d,", (*Counts).Unbonded);
//printf(" Bonded = %d}\n", (*Counts).Bonded);
//// MoleculeType struct
//for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
//  printf("MoleculeType[%d].{Name = %s, ", i, (*MoleculeType)[i].Name);
//  printf("Number = %d, ", (*MoleculeType)[i].Number);
//  printf("Mass = %lf, ", (*MoleculeType)[i].Mass);
//  printf("nBeads = %d, ", (*MoleculeType)[i].nBeads);
//  printf("nBTypes = %d, BTypes = {", (*MoleculeType)[i].nBTypes);
//  for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
//    if (j == 0 && (*MoleculeType)[i].nBTypes == 1) {
//      printf("%d} ", (*MoleculeType)[i].BType[j]);
//    } else if (j == 0) {
//      printf("%d,", (*MoleculeType)[i].BType[j]);
//    } else if (j < ((*MoleculeType)[i].nBTypes-1)) {
//      printf(" %d,", (*MoleculeType)[i].BType[j]);
//    } else {
//      printf(" %d}, ", (*MoleculeType)[i].BType[j]);
//    }
//  }
//  printf("nBonds = %d, Bond = {", (*MoleculeType)[i].nBonds);
//  for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
//    if (j == 0 && (*MoleculeType)[i].nBonds == 1) {
//      printf("%d-%d}", (*MoleculeType)[i].Bond[j][0], (*MoleculeType)[i].Bond[j][1]);
//    } else if (j == 0) {
//      printf("%d-%d,", (*MoleculeType)[i].Bond[j][0], (*MoleculeType)[i].Bond[j][1]);
//    } else if (j < ((*MoleculeType)[i].nBonds-1)) {
//      printf(" %d-%d,", (*MoleculeType)[i].Bond[j][0], (*MoleculeType)[i].Bond[j][1]);
//    } else {
//      printf(" %d-%d}", (*MoleculeType)[i].Bond[j][0], (*MoleculeType)[i].Bond[j][1]);
//    }
//  }
//  printf("}\n");
//}
//// BeadType struct
//for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
//  printf("BeadType[%d].{Name = %s, ", i, (*BeadType)[i].Name);
//  printf("Number = %d, ", (*BeadType)[i].Number);
//  printf("Use = %s, ", (*BeadType)[i].Use? "Yes":"No");
//  printf("Charge = %f, ", (*BeadType)[i].Charge);
//  printf("Mass = %f}\n", (*BeadType)[i].Mass);
//}
//// Molecule struct
//for (int i = 0; i < (*Counts).Molecules; i++) {
//  int type = (*Molecule)[i].Type;
//  printf("Molecule[%d].{Type = %d", i, (*Molecule)[i].Type);
//  printf(" (%s),", (*MoleculeType)[type].Name);
//  printf(" Bead (%d) = {", (*MoleculeType)[type].nBeads);
//  for (int j = 0; j < (*MoleculeType)[type].nBeads; j++) {
//    if (j == 0 && j == (*MoleculeType)[type].nBeads) {
//      printf("%d}", (*Molecule)[i].Bead[j]);
//    } else if (j == 0) {
//      printf("%d,", (*Molecule)[i].Bead[j]);
//    } else if (j < ((*MoleculeType)[type].nBeads-1)) {
//      printf(" %d,", (*Molecule)[i].Bead[j]);
//    } else {
//      printf(" %d}", (*Molecule)[i].Bead[j]);
//    }
//  }
//  printf("}\n");
//}
//// Bead struct
//for (int i = 0; i < (*Counts).Beads; i++) {
//  printf("Bead[%d].{Type = %d", i, (*Bead)[i].Type);
//  printf(" (%s),", (*BeadType)[(*Bead)[i].Type].Name);
//  printf(" Molecule = %d,", (*Bead)[i].Molecule);
//  printf(" Index = %d}\n", (*Bead)[i].Index);
//} //}}}

  free(index);

  return indexed;
} //}}}

// MoveCOMMolecules() - auxiliary //{{{
/*
 * Function to move molecules so that their centre of mass is inside the
 * simulation box. Assumes the moleces are joined.
 */
void MoveCOMMolecules(Counts Counts, Vector BoxLength,
                      BeadType *BeadType, Bead **Bead,
                      MoleculeType *MoleculeType, Molecule *Molecule) {

  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    Vector com; // coordinate vector for centre of mass //{{{
    com.x = 0;
    com.y = 0;
    com.z = 0; //}}}

    // calculate centre of mass
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int id = Molecule[i].Bead[j];

      com.x += (*Bead)[id].Position.x * BeadType[(*Bead)[id].Type].Mass;
      com.y += (*Bead)[id].Position.y * BeadType[(*Bead)[id].Type].Mass;
      com.z += (*Bead)[id].Position.z * BeadType[(*Bead)[id].Type].Mass;
    }

    com.x /= MoleculeType[i].Mass;
    com.y /= MoleculeType[i].Mass;
    com.z /= MoleculeType[i].Mass;

    // move molecule 'j' if its com is outside box
    // x-direction //{{{
    if (com.x >= BoxLength.x) {
      for (int j = 0; j < MoleculeType[type].nBeads; j++) {
        (*Bead)[Molecule[i].Bead[j]].Position.x -= BoxLength.x;
      }
    } else if (com.x < 0) {
      for (int j = 0; j < MoleculeType[type].nBeads; j++) {
        (*Bead)[Molecule[i].Bead[j]].Position.x += BoxLength.x;
      }
    } //}}}
    // y-direction //{{{
    if (com.y >= BoxLength.y) {
      for (int j = 0; j < MoleculeType[type].nBeads; j++) {
        (*Bead)[Molecule[i].Bead[j]].Position.y -= BoxLength.y;
      }
    } else if (com.y < 0) {
      for (int j = 0; j < MoleculeType[type].nBeads; j++) {
        (*Bead)[Molecule[i].Bead[j]].Position.y += BoxLength.y;
      }
    } //}}}
    // z-direction //{{{
    if (com.z >= BoxLength.z) {
      for (int j = 0; j < MoleculeType[type].nBeads; j++) {
        (*Bead)[Molecule[i].Bead[j]].Position.z -= BoxLength.z;
      }
    } else if (com.z < 0) {
      for (int j = 0; j < MoleculeType[type].nBeads; j++) {
        (*Bead)[Molecule[i].Bead[j]].Position.z += BoxLength.z;
      }
    } //}}}
  }
} //}}}

// ReadCoordinates() //{{{
/**
 * Function reading coordinates from .vcf file with indexed timesteps (\ref IndexedCoorFile).
 */
int ReadCoordinates(bool indexed, FILE *vcf_file, Counts Counts, Bead **Bead, char **stuff) {

  // initial stuff //{{{
  (*stuff)[0] = '\0'; // no comment line
  char line[1024];
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
    // allocate helper array of coordinates //{{{
    double imp = 1000000; // a value impossible for bead's coordinate
    Vector *pos = malloc(Counts.BeadsInVsf*sizeof(Vector));
    for (int i = 0; i < Counts.BeadsInVsf; i++) {
      pos[i].x = imp;
    } //}}}

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
      pos[index].x = atof(split[1]);
      pos[index].y = atof(split[2]);
      pos[index].z = atof(split[3]);
    } //}}}

    // copy coordinates to Bead struct //{{{
    int count = 0;
    for (int i = 0; i < Counts.BeadsInVsf; i++) {
      if (pos[i].x != imp) { // i.e., if bead i is present in the timestep
        (*Bead)[count].Position.x = pos[i].x;
        (*Bead)[count].Position.y = pos[i].y;
        (*Bead)[count].Position.z = pos[i].z;
        if ((++count) == Counts.Beads) {
          break;
        }
      }
    } //}}}
    free(pos);
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
  char line[1024];
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
                    MoleculeType *MoleculeType, Molecule **Molecule) {

  bool error = false;

  // is there a Step? I.e., isn't this the line 'Last Step'?
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
        fscanf(agg_file, "%d", &(*Aggregate)[i].Molecule[j]);

        (*Aggregate)[i].Molecule[j]--; // in agg file the numbers correspond to vmd
      }

      while (getc(agg_file) != '\n')
       ; //}}}

      // read monomeric beads in Aggregate 'i' //{{{
      fscanf(agg_file, "%d :", &(*Aggregate)[i].nMonomers);
      for (int j = 0; j < (*Aggregate)[i].nMonomers; j++) {
        fscanf(agg_file, "%d", &(*Aggregate)[i].Monomer[j]);
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
    } //}}}

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

// WriteCoorIndexed() //{{{
/**
 * Function writing coordinates to a `.vcf` file. According to the Use flag
 * in BeadType structure only certain bead types will be saved into the
 * indexed timestep in .vcf file (\ref IndexedCoorFile).
 */
void WriteCoorIndexed(FILE *vcf_file, Counts Counts,
                      BeadType *BeadType, Bead *Bead,
                      MoleculeType *MoleculeType, Molecule *Molecule,
                      char *stuff) {

  // print blank line
  putc('\n', vcf_file);
  // print comment at the beginning of a timestep if present in initial vcf file
  fprintf(vcf_file, "%s", stuff);
  // print 'indexed' on the next
  fprintf(vcf_file, "indexed\n");

  for (int i = 0; i < Counts.Beads; i++) {
    int type_b = Bead[i].Type;
    if (BeadType[type_b].Write) {

      if (Bead[i].Molecule != -1) { // bead in a molecule

        int mol_type = Molecule[Bead[i].Molecule].Type;
        if (MoleculeType[mol_type].Write) {
          fprintf(vcf_file, "%6d %7.3f %7.3f %7.3f\n", Bead[i].Index,
                                                       Bead[i].Position.x,
                                                       Bead[i].Position.y,
                                                       Bead[i].Position.z);
        }
      } else { // monomer bead
        fprintf(vcf_file, "%6d %7.3f %7.3f %7.3f\n", Bead[i].Index,
                                                     Bead[i].Position.x,
                                                     Bead[i].Position.y,
                                                     Bead[i].Position.z);
      }
    }
  }
} //}}}

// WriteCoorXYZ() //{{{
void WriteCoorXYZ(FILE *xyz_file, Counts Counts,
                  BeadType *BeadType, Bead *Bead) {

  // count beads to write
  int count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Write) {
      count += BeadType[i].Number;
    }
  }

  // print number of beads to file
  fprintf(xyz_file, "%d\n\n", count);

  // print coordinates
  for (int i = 0; i < Counts.Beads; i++) {
    int type = Bead[i].Type;
    if (BeadType[type].Write) {
      fprintf(xyz_file, "%8s %7.3f %7.3f %7.3f\n", BeadType[type].Name, Bead[i].Position.x, Bead[i].Position.y, Bead[i].Position.z);
    }
  }
} //}}}

// FindBeadType() //{{{
int FindBeadType(char *name, Counts Counts, BeadType *BeadType) {
  int type;

  // compare give 'name' with all known bead types & return bead type id
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (strcmp(name, BeadType[i].Name) == 0) {
      type = i;
      return type;
    }
  }

  // name isn't in BeadType struct
  return -1;
} //}}}

// FindMoleculeType() //{{{
int FindMoleculeType(char *name, Counts Counts, MoleculeType *MoleculeType) {
  int type;

  // compare give 'name' with all known bead types & return bead type id
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (strcmp(name, MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }

  // name isn't in MoleculeType struct
  return(-1);
} //}}}

// Distance() //{{{
/**
 * Function calculating distance vector between two beads. It removes
 * periodic boundary conditions and returns x, y, and z distances in the
 * range <0, BoxLength/2).
 */
Vector Distance(Vector id1, Vector id2, Vector BoxLength) {

  Vector rij;

  // distance vector
  rij.x = id1.x - id2.x;
  rij.y = id1.y - id2.y;
  rij.z = id1.z - id2.z;

  // remove periodic boundary conditions in x-direction
  while (rij.x >= (BoxLength.x/2))
    rij.x = rij.x - BoxLength.x;
  while (rij.x < -(BoxLength.x/2))
    rij.x = rij.x + BoxLength.x;
  // in y-direction
  while (rij.y >= (BoxLength.y/2))
    rij.y = rij.y - BoxLength.y;
  while (rij.y < -(BoxLength.y/2))
    rij.y = rij.y + BoxLength.y;
  // in z-direction
  while (rij.z >= (BoxLength.z/2))
    rij.z = rij.z - BoxLength.z;
  while (rij.z < -(BoxLength.z/2))
    rij.z = rij.z + BoxLength.z;

  return rij;
} //}}}

// RemovePBCMolecules() //{{{
/**
 * Function to remove periodic boundary conditions from all individual
 * molecules, thus joining them
 */
void RemovePBCMolecules(Counts Counts, Vector BoxLength,
                        BeadType *BeadType, Bead **Bead,
                        MoleculeType *MoleculeType, Molecule *Molecule) {

  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBonds; j++) {
      int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
      int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];

      Vector dist = Distance((*Bead)[id1].Position, (*Bead)[id2].Position, BoxLength);

      (*Bead)[id2].Position.x = (*Bead)[id1].Position.x - dist.x;
      (*Bead)[id2].Position.y = (*Bead)[id1].Position.y - dist.y;
      (*Bead)[id2].Position.z = (*Bead)[id1].Position.z - dist.z;
    }
  }
} //}}}

// RemovePBCAggregates() //{{{
/**
 * Function to remove periodic boundary conditions from all aggregates,
 * thus joining them.
 */
void RemovePBCAggregates(double distance, Aggregate *Aggregate, Counts Counts,
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule) {

  bool *moved = malloc(Counts.Molecules*sizeof(bool));

  // go through all aggregates larger than unimers
  for (int i = 0; i < Counts.Aggregates; i++) {

    // negate moved array, while first molecule is not to move //{{{
    for (int j = 1; j < Counts.Molecules; j++) {
      moved[j] = false;
    }
    moved[0] = true; //}}}

    bool done = false;
    while (!done) {

      // go through all molecule pairs
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        for (int k = 0; k < Aggregate[i].nMolecules; k++) {

          // use only moved molecule 'mol1' and unmoved molecule 'mol2'
          if (moved[j] && !moved[k]) { // automatically follows that j != k
            int mol1 = Aggregate[i].Molecule[j];
            int mol2 = Aggregate[i].Molecule[k];

            // go through all bead pairs in the two molecules
            for (int l = 0; l < MoleculeType[Molecule[mol1].Type].nBeads; l++) {
              for (int m = 0; m < MoleculeType[Molecule[mol2].Type].nBeads; m++) {
                int bead1 = Molecule[mol1].Bead[l];
                int bead2 = Molecule[mol2].Bead[m];

                // use only bead types that were used to assign molecules to aggregates
                if (BeadType[(*Bead)[bead1].Type].Use &&
                    BeadType[(*Bead)[bead2].Type].Use) {

                  // calculate distance between 'bead1' and 'bead2'
                  Vector dist = Distance((*Bead)[bead1].Position, (*Bead)[bead2].Position, BoxLength);
                  dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

                  // move 'mol2' (or 'k') if 'bead1' and 'bead2' are in contact
                  if (dist.x < distance) {

                    // distance vector between 'bead1' and 'bead2' //{{{
                    dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z; //}}}

                    // if 'bead1' and 'bead2' are too far in x-direction, move 'mol2' in x-direction //{{{
                    while (dist.x > (BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x += BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    }
                    while (dist.x <= -(BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x -= BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    } //}}}

                    // if 'bead1' and 'bead2' are too far in y-direction, move 'mol2' in y-direction //{{{
                    while (dist.y > (BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y += BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    }
                    while (dist.y <= -(BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y -= BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    } //}}}

                    // if 'bead1' and 'bead2' are too far in z-direction, move 'mol2' in x-direction //{{{
                    while (dist.z > (BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z += BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    }
                    while (dist.z <= -(BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z -= BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    } //}}}

                    moved[k] = true;

                    // skip remainder of 'mol2' (or 'k')
                    break;
                  }
                }
              }
              // if molekule 'k' (or 'mol2') has been moved, skip also remainder of molecules 'mol1' //{{{
              if (moved[k]) {
                break;
              } //}}}
            }
          }
        }
      }

      // check if all molecules have moved //{{{
      done = true;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (!moved[j]) {
          done = false;
          break;
        }
      } //}}}
    }
  }

  // put aggregates' centre of mass into the simulation box //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    Vector com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead,
                              *Bead, BeadType);

    // by how many BoxLength's should com by moved?
    // for distant aggregates - it shouldn't happen, but better safe than sorry
    IntVector move;
    move.x = com.x / BoxLength.x;
    move.y = com.y / BoxLength.y;
    move.z = com.z / BoxLength.z;
    if (com.x < 0) {
      move.x--;
    }
    if (com.y < 0) {
      move.y--;
    }
    if (com.z < 0) {
      move.z--;
    }

    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      (*Bead)[bead].Position.x -= move.x * BoxLength.x;
      (*Bead)[bead].Position.y -= move.y * BoxLength.y;
      (*Bead)[bead].Position.z -= move.z * BoxLength.z;
    }
  } //}}}

  free(moved);
} //}}}

// RestorePBC() //{{{
/**
 * Function to restore removed periodic boundary conditions. Used in case
 * of cell linked list, because it needs coordinates <0, BoxLength>.
 */
void RestorePBC(Counts Counts, Vector BoxLength, Bead **Bead) {

  for (int i = 0; i < (Counts.Unbonded+Counts.Bonded); i++) {
    // x direction
    while ((*Bead)[i].Position.x >= BoxLength.x) {
      (*Bead)[i].Position.x -= BoxLength.x;
    }
    while ((*Bead)[i].Position.x < 0) {
      (*Bead)[i].Position.x += BoxLength.x;
    }
    // y direction
    while ((*Bead)[i].Position.y >= BoxLength.y) {
      (*Bead)[i].Position.y -= BoxLength.y;
    }
    while ((*Bead)[i].Position.y < 0) {
      (*Bead)[i].Position.y += BoxLength.y;
    }
    // z direction
    while ((*Bead)[i].Position.z >= BoxLength.z) {
      (*Bead)[i].Position.z -= BoxLength.z;
    }
    while ((*Bead)[i].Position.z < 0) {
      (*Bead)[i].Position.z += BoxLength.z;
    }
  }
} //}}}

// CentreOfMass() //{{{
/**
 * Function to calculate centre of mass for a given list of beads.
 */
Vector CentreOfMass(int n, int *list, Bead *Bead, BeadType *BeadType) {

  Vector com;
  com.x = 0;
  com.y = 0;
  com.z = 0;

  for (int i = 0; i < n; i++) {
    com.x += Bead[list[i]].Position.x * BeadType[Bead[list[i]].Type].Mass;
    com.y += Bead[list[i]].Position.y * BeadType[Bead[list[i]].Type].Mass;
    com.z += Bead[list[i]].Position.z * BeadType[Bead[list[i]].Type].Mass;
  }
  com.x /= n;
  com.y /= n;
  com.z /= n;

  return com;
} //}}}

// Gyration() //{{{
/**
 * Function to calculate the principle moments of the gyration tensor.
 */
Vector Gyration(int n, int *list, Counts Counts, Vector BoxLength, BeadType *BeadType, Bead **Bead) {
  // gyration tensor (3x3 array) //{{{
  // use long double to ensure precision -- previous problem with truncation in short chains
  struct Tensor {
    LongVector x, y, z;
  } GyrationTensor;

  GyrationTensor.x.x = 0;
  GyrationTensor.x.y = 0;
  GyrationTensor.x.z = 0;
  GyrationTensor.y.x = 0;
  GyrationTensor.y.y = 0;
  GyrationTensor.y.z = 0;
  GyrationTensor.z.x = 0;
  GyrationTensor.z.y = 0;
  GyrationTensor.z.z = 0; //}}}

// test print of given coordinates -- uncomment if need be //{{{
//for (int i = 0; i < n; i++) {
//  fprintf(stderr, " %10.5f %10.5f %10.5f \n", (*Bead)[list[i]].Position.x,
//                                              (*Bead)[list[i]].Position.y,
//                                              (*Bead)[list[i]].Position.z);
//} //}}}

  Vector com = CentreOfMass(n, list, *Bead, BeadType);
//fprintf(stderr, "%lf %lf %lf\n", com.x, com.y, com.z);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= com.x;
    (*Bead)[list[i]].Position.y -= com.y;
    (*Bead)[list[i]].Position.z -= com.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    GyrationTensor.x.x += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.x;
    GyrationTensor.x.y += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.y;
    GyrationTensor.x.z += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.z;
    GyrationTensor.y.y += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.y;
    GyrationTensor.y.z += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.z;
    GyrationTensor.z.z += (*Bead)[list[i]].Position.z * (*Bead)[list[i]].Position.z;
  }
  GyrationTensor.x.x /= n;
  GyrationTensor.x.y /= n;
  GyrationTensor.x.z /= n;
  GyrationTensor.y.y /= n;
  GyrationTensor.y.z /= n;
  GyrationTensor.z.z /= n; //}}}
//fprintf(stderr, "Tensor: (%lf, %lf, %lf)\n", GyrationTensor.x.x, GyrationTensor.x.y, GyrationTensor.x.z);
//fprintf(stderr, "        (%lf, %lf, %lf)\n", GyrationTensor.x.y, GyrationTensor.y.y, GyrationTensor.y.z);
//fprintf(stderr, "        (%lf, %lf, %lf)\n", GyrationTensor.x.z, GyrationTensor.y.z, GyrationTensor.z.z);

  // char polynomial: a_cube * x^3 + b_cube * x^2 + c_cube * x + d_cube = 0 //{{{
  long double a_cube = -1;
  long double b_cube = GyrationTensor.x.x + GyrationTensor.y.y + GyrationTensor.z.z;
  long double c_cube = - GyrationTensor.x.x * GyrationTensor.y.y
                  - GyrationTensor.x.x * GyrationTensor.z.z
                  - GyrationTensor.y.y * GyrationTensor.z.z
                  + SQR(GyrationTensor.y.z)
                  + SQR(GyrationTensor.x.y)
                  + SQR(GyrationTensor.x.z);
  long double d_cube = + GyrationTensor.x.x * GyrationTensor.y.y * GyrationTensor.z.z
                  + 2 * GyrationTensor.x.y * GyrationTensor.y.z * GyrationTensor.x.z
                  - SQR(GyrationTensor.x.z) * GyrationTensor.y.y
                  - SQR(GyrationTensor.x.y) * GyrationTensor.z.z
                  - SQR(GyrationTensor.y.z) * GyrationTensor.x.x; //}}}
//fprintf(stderr, "character: %lfx^3 + %lfx^2 + %lfx^1 + %lfx^0;\n", a_cube, b_cube, c_cube, d_cube);

  // first root: either 0 or Newton's iterative method to get it //{{{
  long double root0 = 0;
  if (fabs(d_cube) > 0.0000000001L) {
    // derivative of char. polynomial: a_deriv * x^2 + b_deriv * x + c_deriv
    long double a_deriv = 3 * a_cube;
    long double b_deriv = 2 * b_cube;
    long double c_deriv = c_cube;

    long double root1 = 1;

    while (fabs(root0-root1) > 0.0000000001L) {
      long double f_root0 = (a_cube * CUBE(root0) + b_cube * SQR(root0) + c_cube * root0 + d_cube);
      long double f_deriv_root0 = (a_deriv * SQR(root0) + b_deriv * root0 + c_deriv);
      root1 = root0 - f_root0 / f_deriv_root0;

      // swap root0 and root1 for the next iteration
      long double tmp = root0;
      root0 = root1;
      root1 = tmp;
    }
  } //}}}
//fprintf(stderr, "root0=%lf; ", root0);

  // determine paremeters of quadratic equation a_quad * x^2 + b_quad * x + c_quad = 0 //{{{
  // derived by division: (x^3 + (b_cube/a_cube) * x^2 + (c_cube/a_cube) * x + (d_cube/a_cube)):(x - root0)
  long double a_quad = 1;
  long double b_quad = b_cube / a_cube + root0;
  long double c_quad = SQR(root0) + b_cube / a_cube * root0 + c_cube/a_cube; //}}}
//fprintf(stderr, "quad: %lfx^2 + %lfx + %lf; ", a_quad, b_quad, c_quad);

  // calculate & sort eigenvalues //{{{
  LongVector eigen;
  eigen.x = root0; // found out by Newton's method
  // roots of the quadratic equation
  eigen.y = (-b_quad + sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);
  eigen.z = (-b_quad - sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);

  Vector eigen2;
  eigen2.x = eigen.x;
  eigen2.y = eigen.y;
  eigen2.z = eigen.z;
  eigen2 = Sort3(eigen2); //}}}
//fprintf(stderr, "eigen=(%lf, %lf, %lf)\n", eigen.x, eigen.y, eigen.z);

  return eigen2;
} //}}}

// Min3() //{{{
/**
 * Function returning the lowest number from three floats.
 */
double Min3(double x, double y, double z) {

  double min;
  if (x > y) {
    if (y > z) {
      min = z;
    } else {
      min = y;
    }
  } else if (x > z) {
    min = z;
  } else {
    min = x;
  }

  return min;
} //}}}

// Sort3() //{{{
/**
 * Function returning sorted numbers x < y < z.
 */
Vector Sort3(Vector in) {

  Vector out;

  if (in.x < in.y) {
    if (in.y < in.z) {
      out.x = in.x;
      out.y = in.y;
      out.z = in.z;
    } else if (in.x < in.z) {
      out.x = in.x;
      out.y = in.z;
      out.z = in.y;
    } else {
      out.x = in.z;
      out.y = in.x;
      out.z = in.y;
    }
  } else {
    if (in.x < in.z) {
      out.x = in.y;
      out.y = in.x;
      out.z = in.z;
    } else if (in.y < in.z) {
      out.x = in.y;
      out.y = in.z;
      out.z = in.x;
    } else {
      out.x = in.z;
      out.y = in.y;
      out.z = in.x;
    }
  }

  return out;
} //}}}

// FreeBead() //{{{
/**
 * Free memory allocated for Bead struct array. This function makes it
 * easier to add other arrays to the Bead struct in the future
 */
void FreeBead(Counts Counts, Bead **Bead) {
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    free((*Bead)[i].Aggregate);
  }
  free(*Bead);
} //}}}

// FreeMolecule() //{{{
/**
 * Free memory allocated for Molecule struct array. This function makes it
 * easier other arrays to the Molecule struct in the future
 */
void FreeMolecule(Counts Counts, Molecule **Molecule) {
  for (int i = 0; i < Counts.Molecules; i++) {
    free((*Molecule)[i].Bead);
  }
  free(*Molecule);
} //}}}

// FreeMoleculeType() //{{{
/**
 * Free memory allocated for MoleculeType struct array. This function makes
 * it easier other arrays to the MoleculeType struct in the future
 */
void FreeMoleculeType(Counts Counts, MoleculeType **MoleculeType) {
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      free((*MoleculeType)[i].Bond[j]);
    }
    free((*MoleculeType)[i].Bond);
    free((*MoleculeType)[i].BType);
  }
  free(*MoleculeType);
} //}}}

// FreeAggregate() //{{{
/**
 * Free memory allocated for Aggregate struct array. This function makes it
 * easier other arrays to the Aggregate struct in the future
 */
void FreeAggregate(Counts Counts, Aggregate **Aggregate) {
  for (int i = 0; i < Counts.Molecules; i++) {
    free((*Aggregate)[i].Molecule);
    free((*Aggregate)[i].Bead);
    free((*Aggregate)[i].Monomer);
  }
  free(*Aggregate);
} //}}}
