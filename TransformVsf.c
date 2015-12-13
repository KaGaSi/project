#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CStructs.h"
#include "Structure.h"

void ErrorHelp(char cmd[50]) { //{{{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "   TransformVsf <output.vsf> <options>\n\n");
    fprintf(stderr, "   <output.vsf>    output structure file (*.vsf)\n");
    fprintf(stderr, "   <options>\n");
    fprintf(stderr, "      -v   verbose output\n");
    fprintf(stderr, "      -h   print this help and exit\n");
    fprintf(stderr, "      -i   use input .vsf file different from dl_meso.vsf\n\n");
} //}}}

int main(int argc, char *argv[]) {

  // check if correct number of arguments //{{{
  if (argc < 2) {
    fprintf(stderr, "Too little arguments!\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("TransformVsf reads information from FIELD and\n");
      printf("dl_meso.vsf files and creates .vsf structure file\n");
      printf("used for visualisation of trajectory (.vcf files)\n");
      printf("via VMD visualisation tool.\n\n");

      printf("Usage:\n");
      printf("   TransformVsf <output.vsf>\n\n");
      printf("   <output.vsf>    output structure file (*.vsf)\n");
      printf("   <options>\n");
      printf("      -v   verbose output\n");
      printf("      -h   print this help and exit\n");
      printf("      -i   use input .vsf file different from dl_meso.vsf\n\n");
      exit(0);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  putchar('\n'); //}}}

  // <output.vsf> - file name of output structure file (must end with .vsf) //{{{
  char output[32];
  strcpy(output, argv[1]);

  // test if <output.vsf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(output, '.');
  if (!dot || strcmp(dot, ".vsf")) {
    fprintf(stderr, "<output.vsf> does not have .vsf ending!\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -i option - file name of input structure file //{{{
  char input[32];
  input[0] = '\0'; // check if -i option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");
        exit(1);
      }

      // check if .vsf ending is present
      char *vsf = strrchr(argv[i+1], '.');
      if (!vsf || strcmp(vsf, ".vsf")) {
        fprintf(stderr, "'-i' arguments does not have .vsf ending!\n");
        ErrorHelp(argv[0]);
        exit(1);
      } //}}}

      strcpy(input, argv[i+1]);
    }
  }

  // -i option is not used
  if (input[0] == '\0') {
    strcpy(input, "dl_meso.vsf");
  } //}}}

  // variables - structures
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc.

  // vsf input file name
  ReadStructure(input, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

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

  // create & fill output vsf file
  WriteVsf(output, Counts, BeadType, Bead, MoleculeType, Molecule);

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
