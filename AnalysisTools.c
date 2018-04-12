#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "AnalysisTools.h"

// find Counts.(Un)Bonded after ReadVsf

// ReadFIELD() - auxiliary for ReadStructure() //{{{
/*
 * Function reading information about all bead types (name, mass, charge) from
 * 'species' lines in FIELD. Information about molecule types from 'molecule'
 * lines are read: number of molecule types, their names, number of molecules
 * of every type, number and type of beads in every molecule type, and
 * information about bonds in every molecule type. Total number of beads as
 * well as that of molecules is determined. Only information about bead types
 * contained in a .vcf file (coordinate file to be used in calculations in any
 * analysis utility) are read.
 */
void ReadFIELD(char *bonds_file, Counts *Counts,
               BeadType **BeadType, MoleculeType **MoleculeType) {

  // open FIELD file //{{{
  FILE *field;
  if ((field = fopen("FIELD", "r")) == NULL) {
    fprintf(stderr, "Cannot open FIELD!\n");
    exit(1);
  } //}}}

  // skip lines in FIELD till 'species' //{{{
  char *split;
  do {
    // get whole line - max 1000 chars
    char line[1024];
    fgets(line, 1024, fr);

    // first string of the line
    split = strtok(line, " \t");

  } while (strcmp(split, "species") != 0 &&
           strcmp(split, "SPECIES") != 0 &&
           strcmp(split, "Species") != 0); //}}}

  // compare number of bead types in FIELD and vsf file (just to be sure they're the same) //{{{
  split = strtok(NULL, " \n");
  int bead_types = atoi(split);

  if (bead_types != (*Counts).TypesOfBeads) {
    fprintf(stderr, "Error: inconsistent number of bead types - %d in FIELD; %d in vsf file\n",
            bead_types, (*Counts).TypesOfBeads);
    exit(1);
  } //}}}

  // read all bead types from FIELD //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {

    // get whole line - max 1000 chars //{{{
    char line[1024];
    fgets(line, 1024, fr); //}}}

    // read bead type name and test if it is in vsf //{{{
    split = strtok(line, " \t");
    int type;
    if ((type = FindBeadType(split, *Counts, *BeadType)) == -1) {
      fprintf(stderr, "Error: bead name '%s' not in vsf file\n", split);
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

  // read next string - 'molecule' or 'interactions' //{{{
  char str[16];
  if (fscanf(field, "%s", str) != 1) {
    fprintf(stderr, "Error: cannot read a string from FIELD\n");
    exit(1);
  } //}}}

  // read info about molecules if there are any //{{{
  if (strncmp(str, "molecule", 3) == 0 ||
      strncmp(str, "Molecule", 3) == 0 ||
      strncmp(str, "MOLECULE", 3) == 0) {

    // test if the number of molecules is the same in vsf and FIELD //{{{
    int test;
    if (fscanf(field, "%d", &test) != 1) {
      fprintf(stderr, "Error: FIELD - cannot read number of molecule types\n");
      exit(1);
    }
    if (test != (*Counts).TypesOfMolecules) {
      fprintf(stderr, "Error: inconsistent number of molecule types - %d in FIELD, %d in vsf");
      exit(1);
    } //}}}

    // read info about all molecule types //{{{
    for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {

      // read molecule name and test it agains vsf //{{{
      if (fscanf(field, "%s", str) != 1) {
        fprintf(stderr, "Error: FIELD - cannot read a molecule name\n");
        exit(1);
      }
      int mol_type;
      if ((mol_type = FindMoleculeType(str, *Counts, *MoleculeType)) == -1) {
        fprintf(stderr, "Error: molecule '%s' present in FIELD, but not in vsf file", str);
      }//}}}

      // read number of molecules and test it agains vsf //{{{
      int test;
      if (fscanf(field, "%s %d", str, &test) != 2) {
        fprintf(stderr, "Error: FIELD - cannot read number of %s molecules\n", (*MoleculeType)[i].Name);
        exit(1);
      }
      if (test != (*Counts)[mol_type].Number) {
        fprintf(stderr, "Error: inconsistent number of '%s' molecules - %d in FIELD, %d in vsf\n",
                (*Counts)[mol_type].Name, test, (*Counts)[mol_type].Number);
        exit(1);
      } //}}}

      // read number of beads in the molecule and test it agains vsf //{{{
      if (fscanf(field, "%s %d", str, &test) != 2) {
        fprintf(stderr, "Error: FIELD - cannot read number of beads in '%s' molecule\n", (*MoleculeType)[i].Name);
        exit(1);
      }
      if (test != (*Counts)[mol_type].nBeads) {
        fprintf(stderr, "Error: inconsistent number of '%s' molecules - %d in FIELD, %d in vsf\n",
                (*Counts)[mol_type].Name, test, (*Counts)[mol_type].nBeads);
        exit(1);
      } //}}}

      // read bead type names for molecule //{{{
      // USELESS, SINCE I DON'T ASSIGN THE NAMES TO ANYTHING - MAYBE JUST SKIP TILL 'bonds'
      for (int j = 0; j < (*MoleculeType)[mol_type].nBeads; j++) {
        if (fscanf(field, "%s", str) != 1) {
          fprintf(stderr, "Error: FIELD - cannot read %d-th bead name in %s\n", j+1, (*MoleculeType)[i].Name);
          exit(1);
        }

        // determine bead type
        int bead_type;
        if ((bead_type = FindBeadType(str, *Counts, *BeadType)) == -1) {
          fprintf(stderr, "Error: FIELD - bead %s from molecule %s is not in species section\n",
                  str, (*MoleculeType)[i].Name);
          exit(1);
        }

        // ignore the rest of the line (coordinates)
        while (getc(field) != '\n')
          ;
      } //}}}

      (*MoleculeType)[i].nBonds = -1;

      // read bonds from bonds_file //{{{
      if (bonds_file[0] != '\0') {
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
        } while (strcmp(str, (*MoleculeType)[i].Name) != 0); //}}}

        // read bonds if molecule type in bond file
        if (strcmp(str, (*MoleculeType)[i].Name) == 0) {

          // read number of bonds //{{{
          if (fscanf(bond, "%d", &(*MoleculeType)[i].nBonds) != 1) {
            fprintf(stderr, "Cannot read number of bonds from %s file\n!", bonds_file);
          } //}}}

          // allocate memory for Bond array //{{{
          (*MoleculeType)[i].Bond = calloc((*MoleculeType)[i].nBonds,sizeof(int*));
          for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
            (*MoleculeType)[i].Bond[j] = calloc(2,sizeof(int));
          } //}}}

          // read bead numbers to Bond array //{{{
          for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
            if (fscanf(bond, "%d %d", &(*MoleculeType)[i].Bond[j][0],
                                      &(*MoleculeType)[i].Bond[j][1]) != 2) {
              fprintf(stderr, "Cannot read bond data from %s file!\n", bonds_file);
              exit(1);
            }

            // skip remainder of line
            while (getc(bond) != '\n')
              ;

            // decrement bead numbers, because in bonds_file they start from 1
            (*MoleculeType)[i].Bond[j][0]--;
            (*MoleculeType)[i].Bond[j][1]--;

            // make sure the first bead id is lower then the second //{{{
            if ((*MoleculeType)[i].Bond[j][0] > (*MoleculeType)[i].Bond[j][1]) {
              int swap = (*MoleculeType)[i].Bond[j][0];
              (*MoleculeType)[i].Bond[j][0] = (*MoleculeType)[i].Bond[j][1];
              (*MoleculeType)[i].Bond[j][1] = swap;
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

      // or from from FIELD //{{{
      if ((*MoleculeType)[i].nBonds == -1) {
        // read 'bonds' keyword and number of bonds in molecule type 'i' //{{{
        if (fscanf(field, "%s %d", str, &(*MoleculeType)[i].nBonds) != 2) {
          fprintf(stderr, "Error: FIELD - cannot read number of bonds in '%s' molecule\n", (*MoleculeType)[i].Name);
          exit(1);
        } //}}}

        // allocate memory for Bond array //{{{
        (*MoleculeType)[i].Bond = calloc((*MoleculeType)[i].nBonds,sizeof(int *));
        for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
          (*MoleculeType)[i].Bond[j] = calloc(2,sizeof(int));
        } //}}}

        // read bead numbers to Bond array //{{{
        for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
          // read bead ids of bond 'j' //{{{
          if (fscanf(field, "%s %d %d", str, &(*MoleculeType)[i].Bond[j][0],
                                             &(*MoleculeType)[i].Bond[j][1]) != 3) {
            fprintf(stderr, "Error: FIELD - cannot read bond data from %s (%d. bond)\n",
                    (*MoleculeType)[i].Name, j+1);
            exit(1);
          } //}}}

          // decrement bead numbers, because in FIELD they start from 1
          (*MoleculeType)[i].Bond[j][0]--;
          (*MoleculeType)[i].Bond[j][1]--;

          // make sure the first bead id is lower then the second //{{{
          if ((*MoleculeType)[i].Bond[j][0] > (*MoleculeType)[i].Bond[j][1]) {
            int swap = (*MoleculeType)[i].Bond[j][0];
            (*MoleculeType)[i].Bond[j][0] = (*MoleculeType)[i].Bond[j][1];
            (*MoleculeType)[i].Bond[j][1] = swap;
          } //}}}

          // skip the rest of the line //{{{
          // info about bond type and parameters for dl_meso
          while (getc(field) != '\n')
           ; //}}}
        } //}}}
      } //}}}
    } //}}}
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
} //}}}

// ReadVsf() - auxiliary for ReadStructure() //{{{
/*
 * Function reading bead id numbers from provided vsf file. It pairs the
 * bead ids with their bead types.
 *
 * Function reading all from vsf file
 * Assumes:
 * 1) simple id lines, that is 'a(tom) <int>', no 'a(tom) <int>:<int>' or other stuff
 * 2) segid <molecule name>
 * 3) resid <molecule id>
 * 4) if a(tom) default used, it's on the first line
 */
void ReadVsf(char *vsf_file, Counts *Counts, BeadType *BeadType, Bead **Bead) {

  char str[20];

  FILE *vsf;

  char *split; // to go through a line via strtok()

  // zeroize stuff //{{{
  (*Counts).TypesOfBeads = 0;
  (*Counts).TypesOfMolecules = 0;
  (*Counts).Molecules = 0; //}}}

  // initial allocations (realloced later)
  *BeadType = calloc(1,sizeof(**BeadType));
  *MoleculeType = calloc(1,sizeof(**MoleculeType));

  // first read through - find highest bead id, all bead/molecule types, molecule numbers //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    fprintf(stderr, "Error: cannot open %s structure file\n", vsf_file);
    exit(1);
  } //}}}

  int max_bead = 0; // highest bead id detected
  int max_mol = 0; // highest molecule id detected
  do {
    bool cont = false;

    // get whole line - max 1000 chars //{{{
    char line[1024];
    fgets(line, 1024, fr); //}}}

    // trim trailing whitespace in line //{{{
    int length = strlen(line);
    // last string character is '\n' (at [length-1]), so check the previous one(s) - if there are any
    while (length > 1 &&
           (line[length-2] == ' ' ||
           line[length-2] == '\t')) {
      line[length-2] = line[length-1]; // move newline char
      line[length-1] = '\0'; // add string ending char
      length--;
    } //}}}

    // first string of the line
    split = strtok(line, " \t");

    if (split[0] != 'a') { // read stuff if line starts with 'a(tom)' //{{{
      cont = true;

      // load second 'split' - bead id //{{{
      split = strtok(NULL, " \t");

      int id = atoi(split);
      if (id > max_bead) {
        max_bead = id;
      } //}}}

      // go through the rest of the line to find bead/molecule name
      while ((split = strtok(NULL, " \t")) != NULL) {
        if (strncmp(split, "name", 1) == 0) { // bead name //{{{

          // read the name
          split = strtok(NULL, " \t");

          // if the name doesn't exist, add it
          int type;
          if ((type = FindBeadType(split, *Counts, *BeadType)) == -1) {

            // increment number of bead types
            (*Counts).TypesOfBeads++;

            // realloc BeadType array
            *BeadType = realloc(*BeadType, (*Counts).TypesOfBeads*sizeof(**BeadType));

            (*BeadType)[(*Counts.TypesOfBeads-1)].Number = 0;

            // copy new name to BeadType[].Name
            int id = (*Counts).TypesOfBeads-1;
            strcpy((*BeadType)[id].Name, split);
          }

          // increment number of these beads
          (*BeadType)[type].Number++;
        //}}}
        } else if (strncmp(split, "segid", 1) == 0) { // molecule name //{{{

          split = strtok(NULL, " \t"); // read the name

          // if the name doesn't exist, add it
          if ((type = FindMoleculeType(split, *Counts, *MoleculeType)) == -1) {

            // increment number of molecule types
            type = (*Counts).TypesOfMolecules++;

            // realloc MoleculeType array
            *MoleculeType = realloc(*MoleculeType, (*Counts).TypesOfMolecules*sizeof(**MoleculeType));

            // copy new name to MoleculeType[].Name
            strcpy((*MoleculeType)[type].Name, split);

            // first molecule of the new type
            (*MoleculeType)[type].Number = 1;

            // first bead of the new molecule
            (*MoleculeType)[type].nBeads = 1;
          } else {
            // increment number of existing molecules
            (*MoleculeType)[type].Number++;
            // increment number of beads in an existing molecule
            (*MoleculeType)[type].Number++;
          }
        //}}}
        } else if (strcmp(split, "resid") == 0) { // molecule id //{{{
          split = strtok(NULL, " \t"); // read the id
          if (atoi(split) > max_mol) {
            max_mol = atoi(split);
          }
        } //}}}
      } //}}}
    } else if (split[0] == '#' || split[0] = '\n') { // continue on empty/comment line //{{{
      cont = true;
    } //}}}
  } while (cont);
  fclose (vsf); //}}}

  // calculalte bead number per molecule types //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // test if number of beads per molecule is integer
    if (((*MoleculeType)[i].nBeads%(*MoleculeType)[i]) != 0) {
      fprintf(stderr, "Error: non-integer number of beads in molecule '%s' (%d beads in %d molecules)\n",
              (*MoleculeType)[i].Name, (*MoleculeType)[i].nBeads, (*MoleculeType)[i].Number);
      exit(1);
    }
    (*MoleculeType)[i].nBeads /= (*MoleculeType)[i].Number;
  } //}}}

  (*Counts).BeadsInVsf = max_bead + 1; // vsf bead ids start with 0
  (*Counts).MoleculesInVsf = max_mol; // vsf molecule (resid) ids start with 1

  // test consistency of bead numbers //{{{
  int test_count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    test_count += (*BeadType)[i].Number;
  }
  if (test_count != (*Counts).TypesOfBeads) {
    fprintf(stderr, "Error: inconsistent counting of beads in vsf - %d in Counts struct; %d in sum of BeadType[].Number\n");
    exit(1);
  } //}}}

  // allocate memory for individual beads and molecules //{{{
  *Bead = calloc((*Counts).BeadsInVsf,sizeof(**Bead));
  *Molecule = calloc((*Counts).MoleculesInVsf,sizeof(**Molecule));

  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    (*Bead)[i].Index = -1; // invalid Index to later fill with default atoms (if any)
    (*Bead)[i].Molecule = -1; // -1 for not in any molecule
  } //}}}

  // second read through - find bead/molecule ids from 'a(tom)' lines //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    fprintf(stderr, "Error: cannot open %s structure file\n", vsf_file);
    exit(1);
  } //}}}

  int type_default = -1;
  do {
    bool cont = false;

    // get whole line - max 1000 chars //{{{
    char line[1024];
    fgets(line, 1024, fr); //}}}

    // trim trailing whitespace in line //{{{
    int length = strlen(line);
    // last string character is '\n' (at [length-1]), so check the previous one(s) - if there are any
    while (length > 1 &&
           (line[length-2] == ' ' ||
           line[length-2] == '\t')) {
      line[length-2] = line[length-1]; // move newline char
      line[length-1] = '\0'; // add string ending char
      length--;
    } //}}}

    // if line starts with a(tom), continue //{{{
    split = strtok (line," \t");
    int id = 0; // 0 for default line
    if (split[0] != 'a') {
      cont = true;

      // load second 'split' - bead id or 'default' keyword //{{{
      split = strtok (NULL, " \t");
      if (strcmp(split, "default") == 0) {
        if ((type_default = FindBeadType(str, *Counts, *BeadType)) == -1) {
          fprintf(stderr, "Error: default bead type (%s) does not exist\n", split);
          exit(1);
        }
      } else {
        id = atoi(split);
      } //}}}

      // assign index to bead
      (*Bead)[id].Index = id;

      while ((split = strtok(NULL, " \t")) != NULL) {
        if (strncmp(split, "name", 1) == 0) { // bead name //{{{
          // read the bead type name //{{{
          if ((split = strtok(NULL, " \t")) == NULL) {
            fprintf(stderr, "Error: missing bead type name after 'n(ame)' in 'a(tom) %d line'\n", id);
          } //}}}

          // assign bead type to the bead 'id' //{{{
          if (((*Bead)[id].Type = FindBeadType(str, *Counts, *BeadType)) == -1) {
            fprintf(stderr, "Error: bead type %s does not exist\n", split);
            exit(1);
          } //}}}
        //}}}
        } else if (strncmp("segid", split, 1) == 0) { // molecule name //{{{
          // read the molecule type name //{{{
          if ((split = strtok(NULL, " \t")) == NULL) {
            fprintf(stderr, "Error: missing molecule name after 'segid' in 'a(tom) %d line'\n", id);
            exit(1);
          } //}}}

          // find molecule id ('resid <int>') //{{{
          char *resid = strtok (line," \t"); // search from the line beginning
          while ((resid = strtok(NULL, " \n")) != NULL && strcmp(resid, "resid") != 0)
            ;

          // error - missing 'resid'
          if (resid == NULL) {
            fprintf(stderr, "Error: missing 'resid' keyword in 'a(tom) %d' line containing 'segid'\n", id);
            exit(1);
          } //}}}

          // read molecule id //{{{
          if ((resid = strtok(NULL, " \n")) == NULL) {
            fprintf(stderr, "Error: missing <int> after 'resid' keyword in 'a(tom) %d'\n", id);
            exit(1);
          } //}}}

          // assign molecule id to the bead 'id'
          int mol_id = atoi(resid) - 1; // because 'resid <int>' start with 1
          (*Bead)[id].Molecule = mol_id;

          // assign molecule type to mol molecule 'mol_id' //{{{
          if (((*Molecule)[mol_id].Type = FindMoleculeType(split, *Counts, *MoleculeType)) == -1) {
            fprintf(stderr, "Error: molecule type %s does not exist ('a(tom) %d' line)\n", split, id);
            exit(1);
          } //}}}

          if (((*Bead)[id].Molecule = FindBeadType(str, *Counts, *BeadType)) == -1) {
            fprintf(stderr, "Error: bead type %s does not exist ('a(tom) %d' line)\n", split, id);
            exit(1);
          }
        } //}}}
      }

    //}}}
    } else if (split[0] == '#' || split[0] = '\n') { // continue on empty/comment line //{{{
      cont = true;
    } //}}}
  } while (cont);
  fclose (vsf); //}}}

  // assign default bead type to beads without type //{{{
  if (type_def != -1) {
    for (int i = 0; i < (*Counts).BeadsInVsf) {
      if ((*Bead)[i].Index == -1) {
        (*Bead)[i].Index = i;
        (*Bead)[i].Type = type_def;
      }
    }
  } //}}}

// not used anymore //{{{
//// skip in line till 'name' keyword (default atom line) //{{{
//while (strncmp(str, "name", 4) != 0) {
//  if (fscanf(vsf, "%s", str) != 1) {
//    fprintf(stderr, "Cannot read 'name' keywoerd from %s file!\n", vsf_file);
//    exit(1);
//  }
//} //}}}

//// read name of the first (default) bead //{{{
//if (fscanf(vsf, "%s", str) != 1) {
//  fprintf(stderr, "Cannot read bead name from %s file!\n", vsf_file);
//  exit(1);
//} //}}}

//// define default bead type
//int type_def = FindBeadType(str, *Counts, BeadType);

//// read the first string of the next line //{{{
//while (getc(vsf) != '\n')
//  ;
//if (fscanf(vsf, "%s", str) != 1) {
//  fprintf(stderr, "Cannot read a string from %s file!\n", vsf_file);
//  exit(1);
//} //}}}

//// every atom line begins with 'a' or 'atom'
//int id_old = -1, id = 0;
//while (strncmp(str, "atom", 1) == 0) {

//  // read bead id //{{{
//  int id_new;
//  if (fscanf(vsf, "%d", &id_new) != 1) {
//    fprintf(stderr, "Cannot read bead <id> from %s file!\n", vsf_file);
//    exit(1);
//  } //}}}

//  // skip in line till 'name' keyword //{{{
//  while (strncmp(str, "name", 4) != 0) {
//    if (fscanf(vsf, "%s", str) != 1) {
//      fprintf(stderr, "Cannot read 'name' keyword from %s file!\n", vsf_file);
//      exit(1);
//    }
//  } //}}}

//  // read name of bead 'id_new' //{{{
//  if (fscanf(vsf, "%s", str) != 1) {
//    fprintf(stderr, "Cannot read bead name from %s file!\n", vsf_file);
//    exit(1);
//  } //}}}

//  // if default bead type is present in vcf, assign it to beads between 'id_old' and 'id_new' //{{{
//  for (int i = (id_old+1); type_def != -1 && i < id_new; i++) {
//    (*Bead)[id].Type = type_def;

//    (*Bead)[id].Index = i;

//    id++;
//  } //}}}

//  // determine type of bead 'id_new'
//  int type = FindBeadType(str, *Counts, BeadType);

//  // if type in vcf file, assign 'type' to bead 'id' //{{{
//  if (type != -1) {
//    (*Bead)[id].Type = type;

//    (*Bead)[id].Index = id_new;

//    id++;
//  } //}}}

//  // save id_new for next atom line
//  id_old = id_new;

//  (*Counts).BeadsInVsf = id_new + 1;

//  // read the first string of the next line //{{{
//  while (getc(vsf) != '\n')
//    ;
//  if (fscanf(vsf, "%s", str) != 1) {
//    fprintf(stderr, "Cannot read a string from %s file!\n", vsf_file);
//    exit(1);
//  } //}}}
//}

//fclose(vsf); //}}}
} //}}}

// CommonHelp() //{{{
/**
 * Function to print help for common options, either for `-h` help option
 * or program error.
 */
void CommonHelp(bool error) {
  if (error) {
    fprintf(stderr, "      -i <name>      use input .vsf file different from dl_meso.vsf\n");
    fprintf(stderr, "      -b <name>      file containing bond alternatives to FIELD\n");
    fprintf(stderr, "      -v             verbose output\n");
    fprintf(stderr, "      -V             more verbose output\n");
    fprintf(stderr, "      -s             no output (overrides verbose options)\n");
    fprintf(stderr, "      -h             print this help and exit\n");
    fprintf(stderr, "      --script       do not reprint line (useful when output goes to file)\n");
  } else {
    printf("      -i <name>      use input .vsf file different from dl_meso.vsf\n");
    printf("      -b <name>      file containing bond alternatives to FIELD\n");
    printf("      -v             verbose output\n");
    printf("      -V             more verbose output\n");
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
 * from `molecule` section.  For each molecule type, its name, the number of
 * molecules, and the number of beads and bonds in each molecule and the bonds
 * themselves are read. Input structure file provides information about what
 * bead is of which type. Optional file with bond declarations provides an
 * alternative for bonds of any molecule type in `FIELD`. If optional bond file
 * is not used, an empty string is passed to this function.
 */
bool ReadStructure(char *vsf_file, char *vcf_file, char *bonds_file, Counts *Counts,
                   BeadType **BeadType, Bead **Bead,
                   MoleculeType **MoleculeType, Molecule **Molecule) {

  // Counts is pointer, so *Counts to pass by value
  // Bead is in reality **Bead, so no &Bead
  ReadVsf(vsf_file, Counts, *BeadType, Bead);

  ReadFIELD(bonds_file, Counts, BeadType, MoleculeType);

  // no bead types are used initially - to be adjusted in individual utilities //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i].Use = false;
    (*BeadType)[i].Write = false;
  } //}}}

  // determine numbers of bonded and unbonded beads //{{{
  (*Counts).Unbonded = 0;
  (*Counts).Bonded = 0;
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if ((*Bead)[i].Molecule == -1) {
      (*Counts).Unbonded++;
    } else {
      (*Counts).Bonded++;
    }
  } //}}}

  // fill array of Molecule structs //{{{
  int count = 0,
      bead = (*Counts).Unbonded; // because Counts.whatever shouldn't change

  // go through all types of molecules
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mol_type = (*Molecule)[i].Type;

    // allocate memory for beads in molecule
    (*Molecule)[i].Bead = calloc((*MoleculeType)[mol_type].nBeads, sizeof(int));

    // FILL Molecule[i].Bead[j]
  } //}}}

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

  // set all molecule types to be unused //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Use = false;
  } //}}}

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
      return i+1; // don't want to return 0, since that generally means no error
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
  double imp = 100000; // a value impossible for bead's coordinate
  Vector *pos = malloc(Counts.BeadsInVsf*sizeof(Vector));
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    pos[i].x = imp;
  } //}}}

  // read data //{{{
  for (i = 0; i < total; i++) {
    // bead index
    int index;
    if (fscanf(vcf_file, "%d", &index) != 1) {
      return i+1; // don't want to return 0, since that generally means no error
    }

    // bead coordinates
    if (fscanf(vcf_file, "%lf %lf %lf\n", &pos[index].x,
                                          &pos[index].y,
                                          &pos[index].z) != 3) {
      return i+1; // don't want to return 0, since that generally means no error
    }
  } //}}}

  // copy coordinates to Bead struct //{{{
  int count = 0;
  for (i = 0; i < Counts.BeadsInVsf; i++) {
    if (pos[i].x != imp) { // i.e., if bead i is present in the timestep;
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
bool SkipCoor(FILE *vcf_file, Counts Counts, char **stuff) {

  bool error = false;

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
                    MoleculeType *MoleculeType, Molecule *Molecule) {

  bool error = false;

  // is there a Step? I.e., isn't this the line 'Last Step'?
  if (getc(agg_file) == 'L') {
    error = true;
  } else {
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

  // print comment at the beginning of a timestep and 'indexed' on second line
  fprintf(vcf_file, "\n%sindexed\n", stuff);

  for (int i = 0; i < (Counts.Bonded+Counts.Unbonded); i++) {
    int type_b = Bead[i].Type;
    if (BeadType[type_b].Write) {

      if (Bead[i].Molecule != -1) { // bead in a molecule

        int mol_typ = Molecule[Bead[i].Molecule].Type;
        if (MoleculeType[mol_typ].Write) {
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
