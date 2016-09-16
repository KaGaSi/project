#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "AnalysisTools.h"

// ReadFIELD() - auxiliary for ReadStructure() //{{{
/*
 * Function reading information about all bead types (name, mass, charge) from
 * 'species' lines in FIELD. Information about molecule types from 'molecule'
 * lines are read: number of molecule types, their names, number of molecules
 * of every type, number and type of beads in every molecule type and
 * information about bonds in every molecule type. Total number of beads as
 * well as that of molecules is determined. Only information about bead
 * types contained in a .vcf file (coordinate file to be used in
 * calculations in any analysis utility) are read.
 */
bool ReadFIELD(char *bonds_file, char *vcf_file, Counts *Counts,
               BeadType **BeadType, MoleculeType **MoleculeType) {

  // zeroize all Counts structure //{{{
  (*Counts).TypesOfBeads = 0;
  (*Counts).TypesOfMolecules = 0;
  (*Counts).Unbonded = 0;
  (*Counts).Bonded = 0;
  (*Counts).Molecules = 0; //}}}

  // open FIELD file //{{{
  FILE *field;
  if ((field = fopen("FIELD", "r")) == NULL) {
    fprintf(stderr, "Cannot open FIELD!\n");
    exit(1);
  } //}}}

  // skip lines in FIELD till 'species' //{{{
  char str[16];
  do {
    if (fscanf(field, "%s", str) != 1) {
      fprintf(stderr, "Cannot read a string from FIELD!\n");
      exit(1);
    }
  } while (strcmp(str, "species") != 0 &&
           strcmp(str, "SPECIES") != 0 &&
           strcmp(str, "Species") != 0); //}}}

  // read number of bead types in FIELD //{{{
  int types;
  if (fscanf(field, "%d", &types) != 1) {
    fprintf(stderr, "Cannot read number of species types from FIELD!\n");
    exit(1);
  }

  // skip the rest of the line
  while (getc(field) != '\n')
    ; //}}}

  // allocate BeadType structure
  *BeadType = calloc(types,sizeof(**BeadType));

  // open vcf file (if exists) //{{{
  bool no_vcf = false;
  FILE *vcf;
  if (vcf_file != '\0') {
    if ((vcf = fopen(vcf_file, "r")) == NULL) {
      fprintf(stderr, "Cannot open %s for reading!\n", vcf_file);
      exit(1);
    }
  } else {
   no_vcf = true;
  } //}}}

  // if first character in vcf file is '#' => read info for indexed timesteps //{{{
  bool indexed = false;
  int test;
  while (!no_vcf && (test = getc(vcf)) == '#') {
    indexed = true;

    // read bead type name from vcf //{{{
    char str[16];
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read bead type name from %s!\n", vcf_file);
      exit(1);
    }

    // skip the rest of the line in vcf
    while (getc(vcf) != '\n')
      ; //}}}

    // save pointer position in FIELD file //{{{
    fpos_t pos;
    fgetpos(field, &pos); //}}}

    // go through all bead types in FIELD //{{{
    for (int i = 0; i < types; i++) {

      // read name, mass, charge of a bead type //{{{
      char name[32];
      double m, q;
      int num;
      if (fscanf(field, "%s %lf %lf %d", name,
                                         &m,
                                         &q,
                                         &num) != 4) {
        fprintf(stderr, "Cannot read species name, mass or charge from FIELD!\n");
        exit(1);
      } //}}}

      // copy info if correct bead type found in FIELD //{{{
      if (strcmp(name, str) == 0) {
        strcpy((*BeadType)[(*Counts).TypesOfBeads].Name, name);
        (*BeadType)[(*Counts).TypesOfBeads].Mass = m;
        (*BeadType)[(*Counts).TypesOfBeads].Charge = q;
        (*BeadType)[(*Counts).TypesOfBeads].Number = num;

        // only unbonded beads are listed in 'species' lines
        (*Counts).Unbonded += (*BeadType)[(*Counts).TypesOfBeads].Number;

        // skip the remaining bead types
        break;
      } //}}}

      // Error - bead type name from vcf file not found in FIELD //{{{
      if (i == types) {
        fprintf(stderr, "Cannot find bead type '%s' from vcf file %s in FIELD!\n", str, vcf_file);
        exit(0);
      } //}}}

      while (getc(field) != '\n')
        ;
    } //}}}

    // restore pointer position in FIELD file
    fsetpos(field, &pos);

    // increment total number of bead types
    (*Counts).TypesOfBeads++;
  }
  if (vcf_file != '\0') {
    fclose(vcf);
  } //}}}

  // if vcf file contains ordered timesteps, read info about bead types from FIELD //{{{
  if (!indexed) {
    (*Counts).TypesOfBeads = types;

    // read name, mass, charge of every bead type
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {

      if (fscanf(field, "%s %lf %lf %d", (*BeadType)[i].Name,
                                        &(*BeadType)[i].Mass,
                                        &(*BeadType)[i].Charge,
                                        &(*BeadType)[i].Number) != 4) {
        fprintf(stderr, "Cannot read species name, mass or charge from FIELD!\n");
        exit(1);
      }

      // only unbonded beads are listed in 'species' lines
      (*Counts).Unbonded += (*BeadType)[i].Number;

      while (getc(field) != '\n')
        ;
    } //}}}
  // or if vcf file contains indexed timesteps, skip 'species' section in FIELD //{{{
  } else {
    for (int i = 0; i < types; i++) {
      while (getc(field) != '\n')
        ;
    }
  } //}}}

  // read next string - 'molecule' or 'interactions' //{{{
  if (fscanf(field, "%s", str) != 1) {
    fprintf(stderr, "Cannot read a string from FIELD!\n");
    exit(1);
  } //}}}

  // read info about molecules if there are any //{{{
  if (strncmp(str, "molecule", 3) == 0 ||
      strncmp(str, "Molecule", 3) == 0 ||
      strncmp(str, "MOLECULE", 3) == 0) {

    // read number of types of molecules //{{{
    if (fscanf(field, "%d", &(*Counts).TypesOfMolecules) != 1) {
      fprintf(stderr, "Cannot read number of molecule types from FIELD!\n");
      exit(1);
    } //}}}

    // allocate MoleculeType structure
    *MoleculeType = calloc((*Counts).TypesOfMolecules,sizeof(**MoleculeType));

    // read info about molecule types //{{{
    int mols = 0; // types of molecules
    for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {

      // skip strings till 'finish' (except for first molecule type) //{{{
      while (i != 0 && strcmp(str, "finish") != 0 &&
                       strcmp(str, "Finish") != 0 &&
                       strcmp(str, "FINISH") != 0) {

        if (fscanf(field, "%s", str) != 1) {
          fprintf(stderr, "Cannot read a string from FIELD!\n");
          exit(1);
        }
      } //}}}

      // read name of molecule type 'i' //{{{
      if (fscanf(field, "%s", (*MoleculeType)[mols].Name) != 1) {
        fprintf(stderr, "Cannot read a string from FIELD!\n");
        exit(1);
      } //}}}

      // read 'nummols' keyword and number of molecules of type 'i' //{{{
      if (fscanf(field, "%s %d", str, &(*MoleculeType)[mols].Number) != 2) {
        fprintf(stderr, "Cannot read molecule type name from FIELD!\n");
        exit(1);
      } //}}}

      // read 'beads' keyword and number of beads in molecule type 'i' //{{{
      if (fscanf(field, "%s %d", str, &(*MoleculeType)[mols].nBeads) != 2) {
        fprintf(stderr, "Cannot read number of molecules from FIELD!\n");
        exit(1);
      } //}}}

      // read bead types in molecule //{{{
      bool in_vcf = false;

      // the number of beads in molecule type can be lower if some bead types not in vcf
      int beads = (*MoleculeType)[mols].nBeads;
      (*MoleculeType)[mols].nBeads = 0;

      for (int j = 0; j < beads; j++) {
        if (fscanf(field, "%s", str) != 1) {
          fprintf(stderr, "Cannot read bead type name from 'Molecule' %s in FIELD", (*MoleculeType)[mols].Name);
          exit(1);
        }

        // determine type of the bead
        int type = FindBeadType(str, *Counts, *BeadType);

        // increment total number of beads of given type if the type is in vcf
        if (type != -1) {
          (*BeadType)[type].Number += (*MoleculeType)[mols].Number;
          (*Counts).Bonded += (*MoleculeType)[mols].Number;

          // to get number of beads in molecules (taking into account missing bead types in vcf)
          (*MoleculeType)[mols].nBeads++;

          in_vcf = true;
        }

        while (getc(field) != '\n')
          ;
      } //}}}

      (*MoleculeType)[mols].nBonds = -1;

      // if molecule type is in vcf search for bond info in bonds_file //{{{
      if (in_vcf && bonds_file[0] != '\0') {
        // open bonds_file //{{{
        FILE *bond;
        if ((bond = fopen(bonds_file, "r")) == NULL) {
          fprintf(stderr, "Cannot open bond file %s!\n", bonds_file);
          exit(1);
        } //}}}

        // search for the name //{{{
        do {
          if (fscanf(bond, "%s", str) != 1) {
            fprintf(stderr, "Cannot read a string from bond file %s!\n", bonds_file);
            exit(1);
          }

          if (getc(bond) != '\n')
            ;

          // test end of file
          int test;
          if ((test = getc(bond)) == EOF)
            break;
          else
            ungetc(test, bond);
        } while (strcmp(str, (*MoleculeType)[mols].Name) != 0); //}}}

        // read bonds if molecule type in bond file
        if (strcmp(str, (*MoleculeType)[mols].Name) == 0) {

          // read number of bonds //{{{
          if (fscanf(bond, "%d", &(*MoleculeType)[mols].nBonds) != 1) {
            fprintf(stderr, "Cannot read number of bonds from %s file\n!", bonds_file);
          } //}}}

          // allocate memory for Bond array //{{{
          (*MoleculeType)[mols].Bond = calloc((*MoleculeType)[mols].nBonds,sizeof(int*));
          for (int j = 0; j < (*MoleculeType)[mols].nBonds; j++) {
            (*MoleculeType)[mols].Bond[j] = calloc(2,sizeof(int));
          } //}}}

          // read bead numbers to Bond array //{{{
          for (int j = 0; j < (*MoleculeType)[mols].nBonds; j++) {
            if (fscanf(bond, "%d %d", &(*MoleculeType)[mols].Bond[j][0],
                                      &(*MoleculeType)[mols].Bond[j][1]) != 2) {
              fprintf(stderr, "Cannot read bond data from %s file!\n", bonds_file);
              exit(1);
            }

            // skip remainder of line
            while (getc(bond) != '\n')
              ;

            // decrement bead numbers, because in bonds_file they start from 1
            (*MoleculeType)[mols].Bond[j][0]--;
            (*MoleculeType)[mols].Bond[j][1]--;

            // make sure the first bead id is lower then the second //{{{
            if ((*MoleculeType)[mols].Bond[j][0] > (*MoleculeType)[mols].Bond[j][1]) {
              int swap = (*MoleculeType)[mols].Bond[j][0];
              (*MoleculeType)[mols].Bond[j][0] = (*MoleculeType)[mols].Bond[j][1];
              (*MoleculeType)[mols].Bond[j][1] = swap;
            } //}}}

          } //}}}

          // skip bonds part in FIELD //{{{
          // read 'bonds' keyword and number of bonds in molecule type 'i'
          int bin;
          if (fscanf(field, "%s %d", str, &bin) != 2) {
            fprintf(stderr, "Cannot read number of bonds from FIELD!\n");
            exit(1);
          }

          // skip 'bin' lines
          for (int j = 0; j < bin; j++) {
            while (getc(field) != '\n')
             ;
          }
        } //}}}

        fclose(bond);
      } //}}}

      // if molecule type is in vcf read bond info from FIELD, if not already read from bonds_file //{{{
      if (in_vcf && (*MoleculeType)[mols].nBonds == -1) {
        // read 'bonds' keyword and number of bonds in molecule type 'i' //{{{
        if (fscanf(field, "%s %d", str, &(*MoleculeType)[mols].nBonds) != 2) {
          fprintf(stderr, "Cannot read number of bonds from FIELD!\n");
          exit(1);
        } //}}}

        // allocate memory for Bond array //{{{
        (*MoleculeType)[mols].Bond = calloc((*MoleculeType)[mols].nBonds,sizeof(int*));
        for (int j = 0; j < (*MoleculeType)[mols].nBonds; j++) {
          (*MoleculeType)[mols].Bond[j] = calloc(2,sizeof(int));
        } //}}}

        // read bead numbers to Bond array //{{{
        for (int j = 0; j < (*MoleculeType)[mols].nBonds; j++) {
          if (fscanf(field, "%s %d %d", str, &(*MoleculeType)[mols].Bond[j][0],
                                             &(*MoleculeType)[mols].Bond[j][1]) != 3) {
            fprintf(stderr, "Cannot read bond data from FIELD!\n");
            exit(1);
          }

          // decrement bead numbers, because in FIELD they start from 1
          (*MoleculeType)[mols].Bond[j][0]--;
          (*MoleculeType)[mols].Bond[j][1]--;

          // make sure the first bead id is lower then the second //{{{
          if ((*MoleculeType)[mols].Bond[j][0] > (*MoleculeType)[mols].Bond[j][1]) {
            int swap = (*MoleculeType)[i].Bond[j][0];
            (*MoleculeType)[mols].Bond[j][0] = (*MoleculeType)[mols].Bond[j][1];
            (*MoleculeType)[mols].Bond[j][1] = swap;
          } //}}}

          while (getc(field) != '\n')
           ;
        } //}}}
      } //}}}

      // increment total numbers of molecules and molecule types //{{{
      if (in_vcf) {
        (*Counts).Molecules += (*MoleculeType)[mols].Number;

        mols++;
      } //}}}
    } //}}}

    (*Counts).TypesOfMolecules = mols;
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

  fclose(field);

  return (indexed);
} //}}}

// ReadVsf() - auxiliary for ReadStructure() //{{{
/*
 * Function reading bead id numbers from provided vsf file. It pairs the
 * bead ids with their bead types.
 */
void ReadVsf(char *vsf_file, Counts *Counts, BeadType *BeadType, Bead **Bead) {

  char str[20];

  FILE *vsf;

  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", vsf_file);
    exit(1);
  } //}}}

  // read first line -- 'a(tom) default' //{{{
  if (fscanf(vsf, "%s %s", str, str) != 2) {
    fprintf(stderr, "Cannot read strings from %s file!\n", vsf_file);
    exit(1);
  } //}}}

  // skip in line till 'n(ame)' keyword //{{{
  while (strncmp(str, "name", 1) != 0) {
    if (fscanf(vsf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read 'n(ame)' keywoerd from %s file!\n", vsf_file);
      exit(1);
    }
  } //}}}

  // read name of the first (default) bead //{{{
  if (fscanf(vsf, "%s", str) != 1) {
    fprintf(stderr, "Cannot read bead name from %s file!\n", vsf_file);
    exit(1);
  } //}}}

  // define default bead type
  int type_def = FindBeadType(str, *Counts, BeadType);

  // read the first string of the next line //{{{
  while (getc(vsf) != '\n')
    ;
  if (fscanf(vsf, "%s", str) != 1) {
    fprintf(stderr, "Cannot read a string from %s file!\n", vsf_file);
    exit(1);
  } //}}}

  // every atom line begins with 'a' or 'atom'
  int id_old = -1, id = 0;
  while (strncmp(str, "atom", 1) == 0) {

    // read bead id //{{{
    int id_new;
    if (fscanf(vsf, "%d", &id_new) != 1) {
      fprintf(stderr, "Cannot read bead <id> from %s file!\n", vsf_file);
      exit(1);
    } //}}}

    // skip in line till 'n(ame)' keyword //{{{
    while (strncmp(str, "name", 1) != 0) {
      if (fscanf(vsf, "%s", str) != 1) {
        fprintf(stderr, "Cannot read 'n(ame)' keyword from %s file!\n", vsf_file);
        exit(1);
      }
    } //}}}

    // read name of bead 'id_new' //{{{
    if (fscanf(vsf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read bead name from %s file!\n", vsf_file);
      exit(1);
    } //}}}

    // if default bead type in vcf, assign it to beads between 'id_old' and 'id_new' //{{{
    for (int i = (id_old+1); type_def != -1 && i < id_new; i++) {
      (*Bead)[id].Type = type_def;

      (*Bead)[id].Index = i;

      id++;
    } //}}}

    // determine type of bead 'id_new'
    int type = FindBeadType(str, *Counts, BeadType);

    // if type in vcf file, assign 'type' to bead 'id' //{{{
    if (type != -1) {
      (*Bead)[id].Type = type;

      (*Bead)[id].Index = id_new;

      id++;
    } //}}}

    // save id_new for next atom line
    id_old = id_new;

    (*Counts).BeadsInVsf = id_new + 1;

    // read the first string of the next line //{{{
    while (getc(vsf) != '\n')
      ;
    if (fscanf(vsf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read a string from %s file!\n", vsf_file);
      exit(1);
    } //}}}
  }

  fclose(vsf);
} //}}}

// CommonHelp() //{{{
/**
 * Function to print help for common options, either for `-h` help option
 * or program error.
 */
void CommonHelp(bool error) {
  if (error) {
    fprintf(stderr, "   <common options>\n");
    fprintf(stderr, "      -i <name>      use input .vsf file different from dl_meso.vsf\n");
    fprintf(stderr, "      -b <name>      file containing bond alternatives to FIELD\n");
    fprintf(stderr, "      -v             verbose output\n");
    fprintf(stderr, "      -V             verbose output with comments from input .vcf file\n");
    fprintf(stderr, "      -s             no output (overrides verbose options)\n");
    fprintf(stderr, "      -h             print this help and exit\n");
    fprintf(stderr, "      --script       do not reprint line (useful when output goes to file)\n");
  } else {
    printf("   <common options>\n");
    printf("      -i <name>      use input .vsf file different from dl_meso.vsf\n");
    printf("      -b <name>      file containing bond alternatives to FIELD\n");
    printf("      -v             verbose output\n");
    printf("      -V             verbose output with comments from input .vcf file\n");
    printf("      -s             no output (overrides verbose options)\n");
    printf("      -h             print this help and exit\n");
    printf("      --script       do not reprint line (useful when output goes to file)\n");
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
    strcpy(*vsf_file, "dl_meso.vsf");
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

  return (error);
} //}}}

// VerboseOutput() //{{{
/**
 * Function providing standard verbose output (for cases when verbose
 * option is used). It prints most of the information about used system.
 */
void VerboseOutput(bool Verbose2, char *input_vcf, char *bonds_file, Counts Counts,
                   BeadType *BeadType, Bead *Bead,
                   MoleculeType *MoleculeType, Molecule *Molecule) {

  if (input_vcf[0] != '\0')
    printf("\n   Read from %s file\n", input_vcf);
  else
    printf("\n   Read from FIELD file\n");
  printf("Counts.{");
  printf("TypesOfBeads =%3d, ", Counts.TypesOfBeads);
  printf("Bonded =%7d, ", Counts.Bonded);
  printf("Unbonded =%7d, ", Counts.Unbonded);
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

    if (bonds_file[0] == '\0') { // all bonds taken from FIELD
      printf(", Bonds from 'FIELD',");
    } else {
      // go through bond file to find out if molecule type 'i' is there
      FILE *bond;
      if ((bond = fopen(bonds_file, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s with '-v' option!\n", bonds_file);
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
    printf(" Use = %3s}\n", MoleculeType[i].Use? "Yes":"No");

  }

  // print bead ids of all molecules if '-V' option is used
  for (int i = 0; Verbose2 && i < Counts.Molecules; i++) {
    printf("Molecule %3d:\n", i+1);
    for (int j = 0; j < MoleculeType[Molecule[i].Type].nBeads; j++) {
      printf(" %d", Molecule[i].Bead[j]);
    }
    printf("\n");
  }
} //}}}

// ReadStructure() //{{{
/**
 * Function reading information about beads and molecules from DL_MESO `FIELD`
 * file and a `.vsf` structure file.  Name, mass and charge of every bead type is
 * read from `species` lines in `FIELD`. The number of molecule types are read
 * from `molecule` section.  For each molecule type its name, the number of
 * molecules, the number of beads and bonds in each molecule and the bonds
 * themselves are read. Input structure file provides information about what
 * bead is of which type. Optional file with bond declarations provides
 * an alternative for bonds of any molecule type in `FIELD`. If optional
 * bond file is not used, an empty string is passed to this function.
 */
bool ReadStructure(char *vsf_file, char *vcf_file, char *bonds_file, Counts *Counts,
                   BeadType **BeadType, Bead **Bead,
                   MoleculeType **MoleculeType, Molecule **Molecule) {

  // Counts is actually *Counts - so no &Counts
  bool indexed = ReadFIELD(bonds_file, vcf_file, Counts, BeadType, MoleculeType);

  // no bead types are used initially - to be adjusted in individual utilities //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i].Use = false;
    (*BeadType)[i].Write = false;
  } //}}}

  // allocate memory for Molecule struct
  *Molecule = calloc((*Counts).Molecules,sizeof(**Molecule));

  // fill array of Molecule structs //{{{
  int count = 0,
      bead = (*Counts).Unbonded; // because Counts.whatever shouldn't change

  // go through all types of molecules
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {

    // go through all molecules of type 'i'
    for (int j = 0; j < (*MoleculeType)[i].Number; j++) {
      (*Molecule)[count].Type = i;

      // allocate memory for beads in molecule 'count'
      (*Molecule)[count].Bead = calloc((*MoleculeType)[i].nBeads,sizeof(int));

      for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
        (*Molecule)[count].Bead[k] = bead++;
      }

      count++;
    }
  } //}}}

  // allocate memory for Bead struct
  *Bead = calloc(((*Counts).Bonded+(*Counts).Unbonded),sizeof(**Bead));

  // Counts is pointer, so *Counts to pass by value
  // Bead is in reality **Bead, so no &Bead
  ReadVsf(vsf_file, Counts, *BeadType, Bead);

  // assign 'in no molecule' and 'in no aggregate' status to all beads //{{{
  for (int i = 0; i < ((*Counts).Unbonded+(*Counts).Bonded); i++) {
    (*Bead)[i].Molecule = -1;
    (*Bead)[i].nAggregates = 0;
  } //}}}

  // allocate memory for indices of aggregates individual beads are in //{{{
  for (int i = 0; i < (*Counts).Unbonded; i++) {
    (*Bead)[i].Aggregate = calloc(20,sizeof(int)); // so much memory to be on the save side
  }
  for (int i = (*Counts).Unbonded; i < ((*Counts).Unbonded+(*Counts).Bonded); i++) {
    (*Bead)[i].Aggregate = calloc(1,sizeof(int)); // so much memory to be on the save side
  }
  //}}}

  // assign correct molecule id for beads in molecules //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    for (int j = 0; j < (*MoleculeType)[(*Molecule)[i].Type].nBeads; j++) {
      (*Bead)[(*Molecule)[i].Bead[j]].Molecule = i;
    }
  } //}}}

  // calculate mass of molecules //{{{
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;

    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      (*MoleculeType)[i].Mass += (*BeadType)[(*Bead)[(*Molecule)[count].Bead[j]].Type].Mass;
    }

    count += (*MoleculeType)[i].Number;
  } //}}}

  // allocate MoleculeType[].BType arrays //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].BType = malloc((*Counts).TypesOfBeads*sizeof(int));

    // initialize array
    for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
      (*MoleculeType)[i].BType[j] = -1;
    }
  } //}}}

  // determine beadtypes in molecules //{{{
  int mols = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nBTypes = 0;

    // go through all beads of the first molecule of the given molecule type
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int id = (*Molecule)[mols].Bead[j];

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
    mols += (*MoleculeType)[i].Number;
  } //}}}

  return (indexed);
} //}}}

// MoveCOMMolecules() - auxiliary //{{{
/*
 * Function to move molecules so that their center of mass is inside the
 * simulation box. Assumes the moleces are joined.
 */
void MoveCOMMolecules(Counts Counts, Vector BoxLength,
                      BeadType *BeadType, Bead **Bead,
                      MoleculeType *MoleculeType, Molecule *Molecule) {

  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    Vector com; // coordinate vector for center of mass //{{{
    com.x = 0;
    com.y = 0;
    com.z = 0; //}}}

    // calculate center of mass
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

// ReadCoorOrdered() //{{{
/**
 * Function reading coordinates from .vcf file with ordered timesteps (\ref OrderedCoorFile).
 */
int ReadCoorOrdered(FILE *vcf_file, Counts Counts, Bead **Bead, char **stuff) {

  // save the first line containing '# <number>' //{{{
  int i = 0;
  while (((*stuff)[i++] = getc(vcf_file)) != '\n')
    ;
  // skip the second line containing 't(imestep)'
  while (getc(vcf_file) != '\n')
    ; //}}}

  for (i = 0; i < (Counts.Unbonded+Counts.Bonded); i++) {
    if (fscanf(vcf_file, "%lf %lf %lf\n", &(*Bead)[i].Position.x,
                                          &(*Bead)[i].Position.y,
                                          &(*Bead)[i].Position.z) != 3) {
      return (i+1); // don't want to return 0, since that generally means no error
    }
  }

  return 0;
} //}}}

// ReadCoorIndexed() //{{{
/**
 * Function reading coordinates from .vcf file with indexed timesteps (\ref IndexedCoorFile).
 */
int ReadCoorIndexed(FILE *vcf_file, Counts Counts, Bead **Bead, char **stuff) {

  // save the first line containing '# <number>' //{{{
  int i = 0;
  while (((*stuff)[i++] = getc(vcf_file)) != '\n')
    ;
  // skip the second line containing 't(imestep)'
  while (getc(vcf_file) != '\n')
    ; //}}}

  int total = Counts.Unbonded + Counts.Bonded;

  // allocate helper array of coordinates //{{{
  Vector *pos = malloc(Counts.BeadsInVsf*sizeof(Vector));
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    pos[i].x = 1000; // a value impossible for bead's coordinate
  } //}}}

  // read data //{{{
  for (i = 0; i < total; i++) {
    // bead index
    int index;
    if (fscanf(vcf_file, "%d", &index) != 1) {
      return (i+1); // don't want to return 0, since that generally means no error
    }

    // bead coordinates
    if (fscanf(vcf_file, "%lf %lf %lf\n", &pos[index].x,
                                          &pos[index].y,
                                          &pos[index].z) != 3) {
      return (i+1); // don't want to return 0, since that generally means no error
    }
  } //}}}

  // copy coordinates to Bead struct //{{{
  int count = 0;
  for (i = 0; i < Counts.BeadsInVsf; i++) {
    if (pos[i].x != 1000) { // a value impossible for bead coordinates
      (*Bead)[count].Position.x = pos[i].x;
      (*Bead)[count].Position.y = pos[i].y;
      (*Bead)[count].Position.z = pos[i].z;

      count++;

      if (count == total)
        break;
    }
  } //}}}

  free(pos);

  return 0;
} //}}}

// SkipCoor() //{{{
/**
 * Function to skip one timestep in coordinates file. It works with both
 * indexed and ordered vcf files.
 */
int SkipCoor(FILE *vcf_file, Counts Counts, char **stuff) {

  // save the first line containing '# <number>' //{{{
  int i = 0;
  while (((*stuff)[i++] = getc(vcf_file)) != '\n')
    ;
  // skip the second line containing 't(imestep)'
  while (getc(vcf_file) != '\n')
    ;
    //}}}

  for (int i = 0; i < (Counts.Unbonded+Counts.Bonded); i++) {
    int test;
    while ((test = getc(vcf_file)) != '\n' && test != EOF)
      ;

    // premature end of file
    if (test == EOF)
      return(1);
  }

  getc(vcf_file);

  return(0);
} //}}}

// ReadAggregates() //{{{
/**
 * Function reading information about aggregates from `.agg` file (\ref AggregateFile) generated by Aggregates utility.
 */
void ReadAggregates(FILE *agg_file, Counts *Counts, Aggregate **Aggregate,
                    MoleculeType *MoleculeType, Molecule *Molecule) {

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

  // fill Aggregate[].Bead arrays //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    (*Aggregate)[i].nBeads = 0;

    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];
      for (int k = 0; k < MoleculeType[Molecule[mol].Type].nBeads; k++) {
        (*Aggregate)[i].Bead[(*Aggregate)[i].nBeads+k] = Molecule[mol].Bead[k];
      }

      // increment number of beads in aggregate 'i'
      (*Aggregate)[i].nBeads += MoleculeType[Molecule[mol].Type].nBeads;
    }
  } //}}}

  // calculate aggregates' masses //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    (*Aggregate)[i].Mass = 0;

    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int type = Molecule[(*Aggregate)[i].Molecule[j]].Type;

      (*Aggregate)[i].Mass += MoleculeType[type].Mass;
    }
  } //}}}
} //}}}

// WriteCoorIndexed() //{{{
/**
 * Function writing coordinates to a `.vcf` file. According to the Use flag
 * in BeadType structure only certain bead types will be saved into the
 * indexed timestep in .vcf file (\ref IndexedCoorFile).
 */
void WriteCoorIndexed(FILE *vcf_file, Counts Counts, BeadType *BeadType, Bead *Bead, char *stuff) {

  // print comment at the beginning of a timestep and 'indexed' on second line
  fprintf(vcf_file, "\n%sindexed\n", stuff);

  for (int i = 0; i < (Counts.Bonded+Counts.Unbonded); i++) {
    if (BeadType[Bead[i].Type].Write) {
      fprintf(vcf_file, "%6d %7.3f %7.3f %7.3f\n", Bead[i].Index,
                                                   Bead[i].Position.x,
                                                   Bead[i].Position.y,
                                                   Bead[i].Position.z);
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
      return (type);
    }
  }

  // name isn't in BeadType struct
  return (-1);
} //}}}

// FindMoleculeType() //{{{
int FindMoleculeType(char *name, Counts Counts, MoleculeType *MoleculeType) {
  int type;

  // compare give 'name' with all known bead types & return bead type id
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (strcmp(name, MoleculeType[i].Name) == 0) {
      type = i;
      return (type);
    }
  }

  // name isn't in MoleculeType struct
  return (-1);
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

  return (rij);
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

// CenterOfMass() //{{{
/**
 * Function to calculate center of mass for a given list of beads.
 */
Vector CenterOfMass(int n, int *list, Bead *Bead, BeadType *BeadType) {

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

  return (com);
} //}}}

// Min3() //{{{
/**
 * Function returning the lowest number from three floats.
 */
double Min3(double x, double y, double z) {

  if (x > y) {
    if (y > z) {
      return (z);
    } else {
      return (y);
    }
  } else if (x > z) {
    return (z);
  } else {
    return (x);
  }
} //}}}

// FreeBead() //{{{
/**
 * Free memory allocated for Bead struct array. This function makes it
 * easier to add other arrays to the Bead struct in the future
 */
void FreeBead(Counts Counts, Bead **Bead) {
  for (int i = 0; i < (Counts.Unbonded+Counts.Bonded); i++) {
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
