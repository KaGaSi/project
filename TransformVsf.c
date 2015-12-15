#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "CStructs.h"
#include "Structure.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   TransformVsf <output.vsf> <options>\n\n");
  fprintf(stderr, "   <output.vsf>    output structure file (*.vsf)\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -i <name>  use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>  file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v         verbose output\n");
  fprintf(stderr, "      -h         print this help and exit\n");
} //}}}

int main(int argc, char *argv[]) {

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
      printf("      -i   use input .vsf file different from dl_meso.vsf\n");
      printf("      -b   file containing bond alternatives to FIELD\n");
      printf("      -v   verbose output\n");
      printf("      -h   print this help and exit\n");
      exit(0);
    }
  } //}}}

  // check if correct number of arguments //{{{
  if (argc < 2) {
    fprintf(stderr, "Too little arguments!\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  putchar('\n'); //}}}

  // -i <name> option - filename of input structure file //{{{
  char vsf_file[32];
  vsf_file[0] = '\0'; // check if -i option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");
        exit(1);
      }

      // check if .vsf ending is present
      char *vsf = strrchr(argv[i+1], '.');
      if (!vsf || strcmp(vsf, ".vsf")) {
        fprintf(stderr, "'-i' arguments does not have .vsf ending!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(vsf_file, argv[i+1]);
    }
  }

  // -i option is not used
  if (vsf_file[0] == '\0') {
    strcpy(vsf_file, "dl_meso.vsf");
  } //}}}

  // -b <name> option - filename of input bond file //{{{
  char bonds_file[32];
  bonds_file[0] = '\0'; // check if -b option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-b") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-b' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(bonds_file, argv[i+1]);
    }
  } //}}}

  // -v option - verbose output //{{{
  bool verbose = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      verbose = true;

      break;
    }
  } //}}}

  int count = 0; // count arguments

  // <output.vsf> - file name of output structure file (must end with .vsf) //{{{
  char output[32];
  strcpy(output, argv[++count]);

  // test if <output.vsf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(output, '.');
  if (!dot || strcmp(dot, ".vsf")) {
    fprintf(stderr, "<output.vsf> %s does not have .vsf ending!\n", output);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // variables - structures
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc.

  // read system information
  ReadStructure(vsf_file, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // print information read from FIELD if option '-v' is used //{{{
  if (verbose) {
    printf("\n   Read from FIELD\n\n");
    printf("Counts.{");
    printf("TypesOfBeads =%3d, ", Counts.TypesOfBeads);
    printf("Bonded =%7d, ", Counts.Bonded);
    printf("Unboded =%7d, ", Counts.Unbonded);
    printf("TypesOfMolecules =%3d, ", Counts.TypesOfMolecules);
    printf("Molecules =%4d}\n", Counts.Molecules);
    printf("\ntotal number of beads: %d\n\n", Counts.Bonded+Counts.Unbonded);

    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      printf("BeadType[%2d].{", i);
      printf("Name =%10s, ", BeadType[i].Name);
      printf("Number =%7d, ", BeadType[i].Number);
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
      if (bonds_file[0] == '\0') {
        printf(", Bonds in '%s'}\n", vsf_file);
      } else {
        printf(", Bonds in '%s'}\n", bonds_file);
      }
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
