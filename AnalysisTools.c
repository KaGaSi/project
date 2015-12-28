#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "AnalysisTools.h"

// ReadFIELD() - auxiliary //{{{
/*
 * Function reading information about all bead types (name, mass, charge) from
 * 'species' lines in FIELD. Information about molecule types from 'molecule'
 * lines are read: number of molecule types, their names, number of molecules
 * of every type, number and type of beads in every molecule type and
 * information about bonds in every molecule type. Total number of beads as
 * well as that of molecules is determined.
 */
void ReadFIELD(char *bonds_file, Counts *Counts, BeadType **BeadType,
    MoleculeType **MoleculeType) {

  // zeroize all Counts structure //{{{
  (*Counts).TypesOfBeads = 0;
  (*Counts).TypesOfMolecules = 0;
  (*Counts).Unbonded = 0;
  (*Counts).Bonded = 0;
  (*Counts).Molecules = 0; //}}}

  char str[16];

  FILE *fr;

  // open FIELD file //{{{
  if ((fr = fopen("FIELD", "r")) == NULL) {
    fprintf(stderr, "Cannot open FIELD!\n");
    exit(1);
  } //}}}

  // skip lines till 'species' //{{{
  do {
    if (fscanf(fr, "%s", str) != 1) {
      fprintf(stderr, "Cannot read a string from FIELD!\n");
      exit(1);
    }
  } while (strcmp(str, "species") != 0 &&
           strcmp(str, "SPECIES") != 0 &&
           strcmp(str, "Species") != 0); //}}}

  // read number of bead types //{{{
  if (fscanf(fr, "%d", &(*Counts).TypesOfBeads) != 1) {
    fprintf(stderr, "Cannot read number of species types from FIELD!\n");
    exit(1);
  }
  while (getc(fr) != '\n')
    ; //}}}

  // allocate BeadType structure
  *BeadType = malloc((*Counts).TypesOfBeads*sizeof(**BeadType));

  // read name, mass, charge of every bead type and sum their numbers //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {

    if (fscanf(fr, "%s %lf %lf %d", (*BeadType)[i].Name,
                                   &(*BeadType)[i].Mass,
                                   &(*BeadType)[i].Charge,
                                   &(*BeadType)[i].Number) != 4) {
      fprintf(stderr, "Cannot read species name, mass or charge from FIELD!\n");
      exit(1);
    }

    // only unbonded beads are listed in 'species' lines
    (*Counts).Unbonded += (*BeadType)[i].Number;

    while (getc(fr) != '\n')
      ;
  } //}}}

  // read next string - 'molecule' or 'interactions' //{{{
  if (fscanf(fr, "%s", str) != 1) {
    fprintf(stderr, "Cannot read a string from FIELD!\n");
    exit(1);
  } //}}}

  // read info about molecules if there are any //{{{
  if (strncmp(str, "molecule", 3) == 0 ||
      strncmp(str, "Molecule", 3) == 0 ||
      strncmp(str, "MOLECULE", 3) == 0) {

    // read number of types of molecules //{{{
    if (fscanf(fr, "%d", &(*Counts).TypesOfMolecules) != 1) {
      fprintf(stderr, "Cannot read number of molecule types from FIELD!\n");
      exit(1);
    } //}}}

    // allocate MoleculeType structure
    *MoleculeType = malloc((*Counts).TypesOfMolecules*sizeof(**MoleculeType));

    // read info about molecule types //{{{
    for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {

      // skip strings till 'finish' (except for first molecule type) //{{{
      while (i != 0 && strcmp(str, "finish") != 0 &&
                       strcmp(str, "Finish") != 0 &&
                       strcmp(str, "FINISH") != 0) {

        if (fscanf(fr, "%s", str) != 1) {
          fprintf(stderr, "Cannot read a string from FIELD!\n");
          exit(1);
        }
      } //}}}

      // read name of molecule type 'i' //{{{
      if (fscanf(fr, "%s", (*MoleculeType)[i].Name) != 1) {
        fprintf(stderr, "Cannot read a string from FIELD!\n");
        exit(1);
      } //}}}

      // read 'nummols' keyword and number of molecules of type 'i' //{{{
      if (fscanf(fr, "%s %d", str, &(*MoleculeType)[i].Number) != 2) {
        fprintf(stderr, "Cannot read molecule type name from FIELD!\n");
        exit(1);
      } //}}}

      // read 'beads' keyword and number of beads in molecule type 'i' //{{{
      if (fscanf(fr, "%s %d", str, &(*MoleculeType)[i].nBeads) != 2) {
        fprintf(stderr, "Cannot read number of molecules from FIELD!\n");
        exit(1);
      } //}}}

      // read bead types in molecule //{{{
      for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
        if (fscanf(fr, "%s", str) != 1) {
          fprintf(stderr, "Cannot read bead type name from 'Molecule' %s in FIELD", (*MoleculeType)[i].Name);
          exit(1);
        }

        // determine type of the bead
        int type = FindType(str, *Counts, *BeadType);

        // increment total number of beads of given type
        (*BeadType)[type].Number += (*MoleculeType)[i].Number;

        while (getc(fr) != '\n')
          ;
      } //}}}

      (*MoleculeType)[i].nBonds = -1;

      // search for bond info in bonds_file //{{{
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

        // read bonds
        if (strcmp(str, (*MoleculeType)[i].Name) == 0) {

          // read number of bonds //{{{
          if (fscanf(bond, "%d", &(*MoleculeType)[i].nBonds) != 1) {
            fprintf(stderr, "Cannot read number of bonds from %s file\n!", bonds_file);
          } //}}}

          // allocate memory for Bond array //{{{
          (*MoleculeType)[i].Bond = malloc((*MoleculeType)[i].nBonds*sizeof(int*));
          for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
            (*MoleculeType)[i].Bond[j] = malloc(2*sizeof(int));
          } //}}}

          // read bead numbers to Bond array //{{{
          for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
            if (fscanf(bond, "%d %d", &(*MoleculeType)[i].Bond[j][0],
                                      &(*MoleculeType)[i].Bond[j][1]) != 2) {
              fprintf(stderr, "Cannot read bond data from %s file!\n", bonds_file);
              exit(1);
            }

            // decrement bead numbers, because in bonds_file they start from 1
            (*MoleculeType)[i].Bond[j][0]--;
            (*MoleculeType)[i].Bond[j][1]--;

            // make sure the first bead id is lower then the second //{{{
            if ((*MoleculeType)[i].Bond[j][0] > (*MoleculeType)[i].Bond[j][1]) {
              int swap = (*MoleculeType)[i].Bond[j][0];
              (*MoleculeType)[i].Bond[j][0] = (*MoleculeType)[i].Bond[j][1];
              (*MoleculeType)[i].Bond[j][1] = swap;
            } //}}}

            while (getc(fr) != '\n')
             ;
          } //}}}
        }

        fclose(bond);
      } //}}}

      // read bond info from FIELD, if not already read from bonds_file //{{{
      if ((*MoleculeType)[i].nBonds == -1) {
        // read 'bonds' keyword and number of bonds in molecule type 'i' //{{{
        if (fscanf(fr, "%s %d", str, &(*MoleculeType)[i].nBonds) != 2) {
          fprintf(stderr, "Cannot read number of bonds from FIELD!\n");
          exit(1);
        } //}}}

        // allocate memory for Bond array //{{{
        (*MoleculeType)[i].Bond = malloc((*MoleculeType)[i].nBonds*sizeof(int*));
        for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
          (*MoleculeType)[i].Bond[j] = malloc(2*sizeof(int));
        } //}}}

        // read bead numbers to Bond array //{{{
        for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
          if (fscanf(fr, "%s %d %d", str, &(*MoleculeType)[i].Bond[j][0],
                                          &(*MoleculeType)[i].Bond[j][1]) != 3) {
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

          while (getc(fr) != '\n')
           ;
        } //}}}
      } //}}}

      // bubble sort the bond array //{{{
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
      } //}}}

      // increment total numbers of molecules and beads
      (*Counts).Molecules += (*MoleculeType)[i].Number;
      (*Counts).Bonded += (*MoleculeType)[i].Number * (*MoleculeType)[i].nBeads;
    } //}}}
  } //}}}

  fclose(fr);
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
  int type_def = FindType(str, Counts, BeadType);

  // assign default type to the first bead
  (*Bead)[0].Type = type_def;

  // read the first string of the next line //{{{
  while (getc(vsf) != '\n')
    ;
  if (fscanf(vsf, "%s", str) != 1) {
    fprintf(stderr, "Cannot read a string from %s file!\n", vsf_file);
    exit(1);
  } //}}}

  // every atom line begins with 'a' or 'atom'
  int id_old = 0;
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

    // determine type of bead 'id_new'
    int type = FindType(str, Counts, BeadType);

    /* assign default type & "in no molecule" status to beads between 'id_old' and 'id_new'*/ //{{{
    for (int i = (id_old+1); i < id_new; i++) {
      (*Bead)[i].Type = type_def;
    } //}}}

    // assign type 'type' to bead 'id_new'
    (*Bead)[id_new].Type = type;

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
void ReadStructure(char *vsf_file, char *bonds_file, Counts *Counts,
                   BeadType **BeadType, Bead **Bead,
                   MoleculeType **MoleculeType, Molecule **Molecule) {

  // Counts is actually *Counts - so no &Counts
  ReadFIELD(bonds_file, Counts, BeadType, MoleculeType);

  // no bead types are used initially - to be adjusted in individual utilities //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i].Use = 0;
  } //}}}

  // allocate memory for Molecule struct
  *Molecule = malloc((*Counts).Molecules*sizeof(**Molecule));

  // fill array of Molecule structs //{{{
  int count = 0,
      bead = (*Counts).Unbonded; // because Counts.whatever shouldn't change

  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].Number; j++) {
      (*Molecule)[count].Type = i;

      // allocate memory for beads in molecule 'count'
      (*Molecule)[count].Bead = malloc((*MoleculeType)[i].nBeads*sizeof(int));

      for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
        (*Molecule)[count].Bead[k] = bead++;
      }

      count++;
    }
  } //}}}

  // allocate memory for Bead struct
  *Bead = malloc(((*Counts).Bonded+(*Counts).Unbonded)*sizeof(**Bead));

  // Counts is pointer, so *Counts to pass by value
  // Bead is in reality **Bead, so no &Bead
  ReadVsf(vsf_file, *Counts, *BeadType, Bead);
} //}}}

// ReadCoorOrdered() //{{{
/**
 * Function reading coordinates from .vcf file with ordered timesteps (\ref OrderedCoorFile).
 */
int ReadCoorOrdered(FILE *vcf_file, Counts Counts, Bead **Bead, char
    **stuff) {

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
int ReadCoorIndexed(FILE *vcf_file, int beadcount, Bead **Bead, char **stuff) {

  // save the first line containing '# <number>' //{{{
  int i = 0;
  while (((*stuff)[i++] = getc(vcf_file)) != '\n')
    ;
  // skip the second line containing 't(imestep)'
  while (getc(vcf_file) != '\n')
    ; //}}}

  for (i = 0; i < beadcount; i++) {

    // read bead index
    int index;
    if (fscanf(vcf_file, "%d", &index) != 1) {
      return (i+1); // don't want to return 0, since that generally means no error
    }

    // read bead coordinates
    if (fscanf(vcf_file, "%lf %lf %lf\n", &(*Bead)[index].Position.x,
                                          &(*Bead)[index].Position.y,
                                          &(*Bead)[index].Position.z) != 3) {
      return (i+1); // don't want to return 0, since that generally means no error
    }
  }

  return 0;
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
    if (BeadType[Bead[i].Type].Use > 0) {
      fprintf(vcf_file, "%6d %7.3f %7.3f %7.3f\n", i, Bead[i].Position.x,
                                                      Bead[i].Position.y,
                                                      Bead[i].Position.z);
    }
  }
} //}}}

// FindType() //{{{
int FindType(char *name, Counts Counts, BeadType *BeadType) {
  int type;

  // compare give 'name' with all known bead types & return bead type id
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (strcmp(name, BeadType[i].Name) == 0) {
      type = i;
      return (type);
    }
  }

  // name isn't in BeadType struct
  fprintf(stderr, "Bead type %s doesn't exist!\n", name);
  exit(1);
} //}}}

// DistanceBetweenBeads() //{{{
/**
 * Function calculating distance vector between two beads. It removes
 * periodic boundary conditions and returns x, y, and z distances in the
 * range <0, BoxLength/2).
 */
Vector DistanceBetweenBeads(int id1, int id2, Bead *Bead, Vector BoxLength)
{

  Vector rij;

  // distance vector
  rij.x = Bead[id1].Position.x - Bead[id2].Position.x;
  rij.y = Bead[id1].Position.y - Bead[id2].Position.y;
  rij.z = Bead[id1].Position.z - Bead[id2].Position.z;

  // remove periodic boundary conditions in x-direction
  if (rij.x > (BoxLength.x/2))
    rij.x = rij.x - BoxLength.x;
  else if (rij.x <= -(BoxLength.x/2))
    rij.x = rij.x + BoxLength.x;
  // in y-direction
  if (rij.y > BoxLength.y/2)
    rij.y = rij.y - BoxLength.y;
  else if (rij.y <= -(BoxLength.y/2))
    rij.y = rij.y + BoxLength.y;
  // in z-direction
  if (rij.z > BoxLength.z/2)
    rij.z = rij.z - BoxLength.z;
  else if (rij.z <= -(BoxLength.z/2))
    rij.z = rij.z + BoxLength.z;

  return (rij);
} //}}}