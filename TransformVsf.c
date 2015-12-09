#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CStructs.h"
#include "Structure.h"

int main(int argc, char *argv[]) {

  // initial stuff -- printing command, argument checks, etc. //{{{
  // print command
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  putchar('\n');

  // sanity check for provided arguments (if any)
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-h") != 0) {
      fprintf(stderr, "Incorrect argument: %s!\n\n", argv[i]);
      fprintf(stderr, "Possible arguments:\n");
      fprintf(stderr,"   -v   verbose output\n");
      fprintf(stderr,"   -h   print this help and exit\n\n");
      exit(1);
    }
  }

  // print program description if option '-h' is used
  if ((argc > 1 && strcmp(argv[1], "-h") == 0) ||
      (argc > 2 && strcmp(argv[2], "-h") == 0)) {
    printf("TransformVsf reads information from FIELD and\n");
    printf("dl_meso.vsf files and creates input.vsf file used\n");
    printf("for further analysis and visualisation of trajectory\n");
    printf("(vcf files) via VMD visualisation tool.\n\n");

    printf("Optional arguments:\n");
    printf("   -v   verbose output\n");
    printf("   -h   print this help and exit\n\n");
    exit(0);
  } //}}}

  // variables - structures
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc.

  // read information from FIELD
  ReadFIELD(&Counts, &BeadType, &MoleculeType);

  // allocate memory for Molecule struct
  Molecule = malloc(Counts.Molecules*sizeof(*Molecule));

  // fill array of Molecule structs //{{{
  int count = 0,
      bead = Counts.Unbonded; // because Counts.whatever shouldn't change
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < MoleculeType[i].Number; j++) {
      Molecule[count].Type = i;

      // allocate memory for beads in molecule 'count'
      Molecule[count].Bead = malloc(MoleculeType[i].nBeads*sizeof(int));

      for (int k = 0; k < MoleculeType[i].nBeads; k++) {
        Molecule[count].Bead[k] = bead++;
      }

      count++;
    }
  } //}}}

  // allocate memory for Bead struct
  Bead = malloc((Counts.Bonded+Counts.Unbonded)*sizeof(*Bead));

  // print information read from FIELD if option '-v' is used //{{{
  if (argc > 1 && strcmp(argv[1], "-v") == 0) {
    printf("\n   Read from FIELD\n\n");
    printf("Counts.{");
    printf("TypesOfBeads =%3d, ", Counts.TypesOfBeads);
    printf("Bonded =%7d, ", Counts.Bonded);
    printf("Unboded =%7d, ", Counts.Unbonded);
    printf("TypesOfMolecules =%3d, ", Counts.TypesOfMolecules);
    printf("Molecules =%4d}\n", Counts.Molecules);
    printf("total number of beads: %d\n\n", Counts.Bonded+Counts.Unbonded);

    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      printf("BeadType[%2d].{", i);
      printf("Name =%10s, ", BeadType[i].Name);
      printf("Number =%3d, ", BeadType[i].Number);
      printf("Charge =%6.2f, ", BeadType[i].Charge);
      printf("Mass =%5.2f}\n", BeadType[i].Mass);
    }
    putchar('\n');

    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      printf("MoleculeType[%d].{", i);
      printf("Name =%10s", MoleculeType[i].Name);
      printf(", Number =%4d", MoleculeType[i].Number);
      printf(", nBeads =%3d", MoleculeType[i].nBeads);
      printf(", nBonds =%3d", MoleculeType[i].nBonds);
      printf(", Bonds = {");
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        // print incremented ids, because in FIELD they start from 1
        printf("%3d %3d", MoleculeType[i].Bond[j][0]+1, MoleculeType[i].Bond[j][1]+1);
        if (j != (MoleculeType[i].nBonds-1))
          printf(", ");
      }
      printf("}}\n");
    }
  } //}}}

  // vsf file name
  char file[16];
  strcpy(file, "dl_meso.vsf");

  // read bead ids from 'file' and assign them proper types
  ReadVsf(file, Counts, BeadType, &Bead);

  if (argc > 1 && strcmp(argv[1], "-v") == 0) { //{{{
    printf("\n");
    for (int i = 0; i < (Counts.Bonded+Counts.Unbonded); i++) {
      printf("Bead[%6d].Type =%3d\n", i, Bead[i].Type);
    }

    for (int i = 0; i < Counts.Molecules; i++) {
      printf("Molecule[%4d].{Type =%3d, Bead = {", i, Molecule[i].Type);

      for (int j = 0; j < MoleculeType[Molecule[i].Type].nBeads; j++) {
        printf("%6d", Molecule[i].Bead[j]);

        if (j != (MoleculeType[Molecule[i].Type].nBeads-1)) {
          printf(", ");
        }
      }
      printf("}}\n");
    }
  } //}}}

  // create & fill input.vsf file
  strcpy(file, "input.vsf");
  WriteVsf(file, Counts, BeadType, Bead, MoleculeType, Molecule);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      free(MoleculeType[i].Bond[j]);
    }
    free(MoleculeType[i].Bond);
  }
  free(MoleculeType);
  free(Bead);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(Molecule[i].Bead);
  }
  free(Molecule); //}}}

  return 0;
}
