#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "AnalysisTools.h"

// ReadFIELD() - auxiliary //{{{
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

  // first character in vcf file is '#' => read info for indexed timesteps //{{{
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

  // read info for case of vcf file with ordered timesteps //{{{
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
  // or skip 'species' section if vcf file with indexed timesteps //{{{
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
          (*MoleculeType)[i].Bond[j][0]--;
          (*MoleculeType)[i].Bond[j][1]--;

          // make sure the first bead id is lower then the second //{{{
          if ((*MoleculeType)[i].Bond[j][0] > (*MoleculeType)[i].Bond[j][1]) {
            int swap = (*MoleculeType)[i].Bond[j][0];
            (*MoleculeType)[i].Bond[j][0] = (*MoleculeType)[i].Bond[j][1];
            (*MoleculeType)[i].Bond[j][1] = swap;
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

// ReadVsf() - auxiliary //{{{
/*
 * Function reading bead id numbers from provided vsf file. It pairs the
 * bead ids with their bead types.
 */
void ReadVsf(char *vsf_file, Counts Counts, BeadType *BeadType, Bead **Bead) {

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
  int type_def = FindBeadType(str, Counts, BeadType);

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
      (*Bead)[id++].Type = type_def;
    } //}}}

    // determine type of bead 'id_new'
    int type = FindBeadType(str, Counts, BeadType);

    // if type in vcf file, assign 'type' to bead 'id' //{{{
    if (type != -1) {
      (*Bead)[id++].Type = type;
    } //}}}

    // save id_new for next atom line
    id_old = id_new;

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
 *
 * \todo Assigning Bead[].Index should be in auxiliary ReadVsf().
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

  // for indexed timestep find out bead types //{{{
  if (indexed) {
    // open vcf file (if exists) //{{{
    FILE *vcf;
    if ((vcf = fopen(vcf_file, "r")) == NULL) {
      fprintf(stderr, "Cannot open %s for reading!\n", vcf_file);
      exit(1);
    } //}}}

    // skip in vcf till first coordinate line //{{{
    char str[16];
    // skipt till 'pbc'
    do {
      // read first string from the line
      if (fscanf(vcf, "%s", str) != 1) {
        fprintf(stderr, "Cannot read a string from %s file while getting bead indices!\n", vcf_file);
        exit(1);
      }

      // skip the rest of the line
      while(getc(vcf) != '\n')
        ;

    } while (strcmp(str, "pbc") != 0);

    // skip remaining three lines
    for (int i = 0; i < 3; i++) {
      while(getc(vcf) != '\n')
        ;
    } //}}}

    // open vsf structure file //{{{
    FILE *vsf;
    if ((vsf = fopen(vsf_file, "r")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", vsf_file);
      exit(1);
    } //}}}

    // default type from vsf //{{{
    // skip in first line till 'n(ame)'
    do {
      if (fscanf(vsf, "%s", str) != 1) {
        fprintf(stderr, "Cannot read a string from first line of %s!\n", vsf_file);
        exit(1);
      }
    } while (strncmp(str, "name", 1) != 0);

    // read bead type name
    if (fscanf(vsf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read a string from first line of %s!\n", vsf_file);
      exit(1);
    }

    int def_type = FindBeadType(str, *Counts, *BeadType);

    // skip the remainder of first line
    while (getc(vsf) != '\n')
      ;
    //}}}

    // assign types for all beads in indexed vcf file //{{{
    int id_vsf = -1; // index number from vsf file
    for (int i = 0; i < ((*Counts).Bonded+(*Counts).Unbonded); i++) {

      // read bead index number from vcf file //{{{
      int id_vcf;
      if (fscanf(vcf, "%d", &id_vcf) != 1) {
        fprintf(stderr, "Cannot read bead id from %s!\n", vcf_file);
        exit(1);
      }
      (*Bead)[i].Index = id_vcf;

      // skp remainder of line
      while(getc(vcf) != '\n')
        ;
      //}}}

      // read index from vsf file till it's lower then that from vcf file //{{{
      while (id_vsf < id_vcf) {

        // skip 'a(tom)' //{{{
        if (fscanf(vsf, "%s", str) != 1) {
          fprintf(stderr, "Cannot read 'a(tom)' string from %s!\n", vsf_file);
          exit(1);
        } //}}}

        // read bead index number //{{{
        if (fscanf(vsf, "%d", &id_vsf) != 1) {
          fprintf(stderr, "Cannot read bead id from %s!\n", vsf_file);
          exit(1);
        } //}}}

        // skip in line till 'n(ame)' //{{{
        do {
          if (fscanf(vsf, "%s", str) != 1) {
            fprintf(stderr, "Cannot read a string from %s!\n", vsf_file);
            exit(1);
          }
        } while (strncmp(str, "name", 1) != 0); //}}}

        // read bead type name //{{{
        if (fscanf(vsf, "%s", str) != 1) {
          fprintf(stderr, "Cannot read a string from first line of %s!\n", vsf_file);
          exit(1);
        } //}}}

        // skip rest of line //{{{
        while (getc(vsf) != '\n')
          ; //}}}
      } //}}}

      // assign default type to bead 'i' //{{{
      if (id_vcf > id_vsf) {
        (*Bead)[i].Type = def_type; //}}}
      // or assign another type //{{{
      } else {
        (*Bead)[i].Type = FindBeadType(str, *Counts, *BeadType);
      } //}}}
    } //}}}

    fclose(vcf);
    fclose(vsf); //}}}
  // or for ordered timestep assign numbers to Bead[].Index //{{{
  } else {
    for (int i = 0; i < ((*Counts).Bonded+(*Counts).Unbonded); i++) {
      (*Bead)[i].Index = i;
    }
  } //}}}

  // Counts is pointer, so *Counts to pass by value
  // Bead is in reality **Bead, so no &Bead
  ReadVsf(vsf_file, *Counts, *BeadType, Bead);

  // assign 'in no molecule' status to all beads //{{{
  for (int i = 0; i < ((*Counts).Unbonded+(*Counts).Bonded); i++) {
    (*Bead)[i].Molecule = -1;
  } //}}}

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

//for (int i = 0; i < ((*Counts).Unbonded+(*Counts).Bonded); i++) {
//  printf("%5d %5d %s %4d\n", i, (*Bead)[i].Index, (*BeadType)[(*Bead)[i].Type].Name, (*Bead)[i].Molecule);
//}

  return (indexed);
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

  for (i = 0; i < (Counts.Unbonded+Counts.Bonded); i++) {

    // read bead coordinates (and unused index)
    int index;
    if (fscanf(vcf_file, "%d %lf %lf %lf\n", &index, &(*Bead)[i].Position.x,
                                                     &(*Bead)[i].Position.y,
                                                     &(*Bead)[i].Position.z) != 4) {
      return (i+1); // don't want to return 0, since that generally means no error
    }
  }

  return 0;
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

// VerboseOutput() //{{{
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

// DistanceBetweenBeads() //{{{
/**
 * Function calculating distance vector between two beads. It removes
 * periodic boundary conditions and returns x, y, and z distances in the
 * range <0, BoxLength/2).
 */
Vector DistanceBetweenBeads(int id1, int id2, Bead *Bead, Vector BoxLength) {

  Vector rij;

  // distance vector
  rij.x = Bead[id1].Position.x - Bead[id2].Position.x;
  rij.y = Bead[id1].Position.y - Bead[id2].Position.y;
  rij.z = Bead[id1].Position.z - Bead[id2].Position.z;

  // remove periodic boundary conditions in x-direction
  while (rij.x > (BoxLength.x/2))
    rij.x = rij.x - BoxLength.x;
  while (rij.x <= -(BoxLength.x/2))
    rij.x = rij.x + BoxLength.x;
  // in y-direction
  while (rij.y > (BoxLength.y/2))
    rij.y = rij.y - BoxLength.y;
  while (rij.y <= -(BoxLength.y/2))
    rij.y = rij.y + BoxLength.y;
  // in z-direction
  while (rij.z > (BoxLength.z/2))
    rij.z = rij.z - BoxLength.z;
  while (rij.z <= -(BoxLength.z/2))
    rij.z = rij.z + BoxLength.z;

  return (rij);
} //}}}

// FillAggregateBeads() //{{{
/**
 * Function to assign bead ids to aggregates accoding to molecules in the
 * aggregates. It essentially duplicates the information for the
 * convenience of easy access to all beads in the aggregates.
 */
void FillAggregateBeads(Aggregate **Aggregate, Counts Counts,
                        MoleculeType *MoleculeType, Molecule *Molecule) {

  // go through all aggregates
  for (int i = 0; i < Counts.Aggregates; i++) {

    // go through all molecules in aggregate 'i'
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];

      // copy all bead in molecule 'mol' to Aggregate struct
      for (int k = 0; k < MoleculeType[Molecule[mol].Type].nBeads; k++) {
        (*Aggregate)[i].Bead[(*Aggregate)[i].nBeads] = Molecule[mol].Bead[k];
        (*Aggregate)[i].nBeads++;
      }
    }
  }
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

      Vector dist = DistanceBetweenBeads(id1, id2, *Bead, BoxLength);

      (*Bead)[id2].Position.x = (*Bead)[id1].Position.x - dist.x;
      (*Bead)[id2].Position.y = (*Bead)[id1].Position.y - dist.y;
      (*Bead)[id2].Position.z = (*Bead)[id1].Position.z - dist.z;
    }
  }
} //}}}

// RemovePBCAggregates() //{{{
/**
 * Function to remove periodic boundary conditions from all aggregates,
 * thus joining them
 */
void RemovePBCAggregates(double distance, Aggregate *Aggregate, Counts Counts,
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule) {

  bool *moved = malloc(Counts.Molecules*sizeof(bool));

  // go through all aggregates larger than unimers
  for (int i = 0; i < Counts.Aggregates; i++) {
//  printf("i = %d\n", i);
//  for (int j = 0; j < Aggregate[i].nMolecules; j++) {
//    printf("%4d", Aggregate[i].Molecule[j]);
//  }
//  putchar('\n');

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
                  Vector dist = DistanceBetweenBeads(bead1, bead2, *Bead, BoxLength);
                  dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

                  // move 'mol2' (or 'k') if 'bead1' and 'bead2' are in contact
                  if (dist.x < distance) {

                    if (i == 1 && mol2 == 2) {
                      printf("%d %d\n", mol1, mol2);
                      printf("%d: %lf %lf %lf\n", (*Bead)[bead1].Index, (*Bead)[bead1].Position.x, (*Bead)[bead1].Position.y, (*Bead)[bead1].Position.z);
                      printf("%d: %lf %lf %lf\n", (*Bead)[bead2].Index, (*Bead)[bead2].Position.x, (*Bead)[bead2].Position.y, (*Bead)[bead2].Position.z);
                      printf("dist: %lf\n", dist.x);
                    }

                    // distance vector between 'bead1' and 'bead2' //{{{
                    dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z; //}}}

                    // if 'bead1' and 'bead2' are too far in x-direction, move 'mol2' in x-direction //{{{
                    if (dist.x > (BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x += BoxLength.x;
                      }
                    } else if (dist.x <= -(BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x -= BoxLength.x;
                      }
                    } //}}}

                    // if 'bead1' and 'bead2' are too far in y-direction, move 'mol2' in y-direction //{{{
                    if (dist.y > (BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y += BoxLength.y;
                      }
                    } else if (dist.y <= -(BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y -= BoxLength.y;
                      }
                    } //}}}

                    // if 'bead1' and 'bead2' are too far in z-direction, move 'mol2' in x-direction //{{{
                    if (dist.z > (BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z += BoxLength.z;
                      }
                    } else if (dist.z <= -(BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z -= BoxLength.z;
                      }
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
//        printf("not moved: %d %d\n", j+1, Aggregate[i].Molecule[j]+1);
          done = false;
          break;
        }
      } //}}}
    }
  }

  free(moved);
} //}}}
