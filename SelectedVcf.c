#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "CStructs.h"
#include "Structure.h"
#include "Aux.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   SelectedVcf <input.vcf> <start> <skip> ");
  fprintf(stderr, "<output.vcf> <type names> <options>\n\n");

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <start>           number of timestep to start from\n");
  fprintf(stderr, "   <skip>            leave out every 'skip' steps\n");
  fprintf(stderr, "   <output.vcf>      output filename (vcf format)\n");
  fprintf(stderr, "   <#> <type names>  number of bead types to save followed by their names\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -i <name>      use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>      file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v             verbose output\n");
  fprintf(stderr, "      -h             print this help and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("SelectedVcf creates new <output.vcf> file from <input.vcf>\n");
      printf("containing only selected bead types. Also <start> timesteps\n");
      printf("can be omitted and every <skip> timesteps can be left out.\n");
      printf("The program uses dl_meso.vsf (or other input structure file)\n");
      printf("and FIELD (along with optional bond file) files to determine\n");
      printf("all information about the system.\n\n");

      printf("Usage:\n");
      printf("   SelectedVcf <input.vcf> <start> <skip> ");
      printf("<output.vcf> <type names> <options>\n\n");

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <start>           number of timestep to start from\n");
      printf("   <skip>            leave out every 'skip' steps\n");
      printf("   <output.vcf>      output filename (vcf format)\n");
      printf("   <#> <type names>  number of bead types to save followed by their names\n");
      printf("   <options>\n");
      printf("      -i <name>      use input .vsf file different from dl_meso.vsf\n");
      printf("      -b <name>      file containing bond alternatives to FIELD\n");
      printf("      -v             verbose output\n");
      printf("      -h             print this help and exit\n");
      exit(0);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  printf("\n\n"); //}}}

  // check if correct number of arguments //{{{
  if (argc < 6) {
    fprintf(stderr, "Too little arguments!\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

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

  // <input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf[32];
  strcpy(input_vcf, argv[++count]);

  // test if <input.vcf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(input_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<input.vcf> '%s' does not have .vcf ending!\n", input_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <start> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <start>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int start = atoi(argv[count]); //}}}

  // <skip> - number of steps to skip per one used //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <skip>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int skip = atoi(argv[count]); //}}}

  // <output.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[32];
  strcpy(output_vcf, argv[++count]);

  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  dot = strrchr(output_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<output.vcf> '%s' does not have .vcf ending!\n", output_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <#> - number of bead types to save //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <skip>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int number_of_saved = atoi(argv[count]); //}}}

  // allocate memory for array holding bead type ids to save //{{{
  int *saved_types;
  saved_types = malloc(number_of_saved*sizeof(int));

  // initialize array
  for (int i = 0; i < number_of_saved; i++) {
    saved_types[i] = -1;
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  ReadStructure(vsf_file, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // <type names> - names of bead types to save //{{{
  { int i = 0;
    while (++count < argc && argv[count][0] != '-') {
      saved_types[i++] = FindType(argv[count], Counts, BeadType);
    }

    // Error - incorrect number of bead type names provided
    if (i != number_of_saved) {
      fprintf(stderr, "Incorrect number of bead type names (should be %d)!\n", number_of_saved);
      ErrorHelp(argv[0]);
      exit(1);
    }

    // Error - does not make sense to use all bead types
    if (number_of_saved == Counts.TypesOfBeads) {
      fprintf(stderr, "All beadtypes are selected which would only copy All.vcf,\n");
      fprintf(stderr, "therefore this is not allowed!\n\n");
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // print information option '-v' is used //{{{
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

    printf("\n   Selected bead types\n");
    printf("number of selected bead types: %d\n", number_of_saved);
    for (int i = 0; i < number_of_saved; i++) {
      printf("BeadType[%2d].{", saved_types[i]);
      printf("Name =%10s, ", BeadType[saved_types[i]].Name);
      printf("Number =%7d, ", BeadType[saved_types[i]].Number);
      printf("Charge =%6.2f, ", BeadType[saved_types[i]].Charge);
      printf("Mass =%5.2f}\n", BeadType[saved_types[i]].Mass);
    }

    printf("\n   starting from %d. timestep\n", start);
    printf("   every %d. timestep used\n", skip+1);
  } //}}}

  // open All.vcf coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen("All.vcf", "r")) == NULL) {
//  fprintf(stderr, "Cannot open file %s!\n", vsf_file);
    fprintf(stderr, "Cannot open file All.vcf!\n");
    exit(1);
  } //}}}

  // read pbc from coordinate file //{{{
  char str[32];
  Vector box_length;
  if (fscanf(vcf, "%s %lf %lf %lf", str, &box_length.x, &box_length.y, &box_length.z) != 4 ||
      strcmp(str, "pbc") != 0) {
    fprintf(stderr, "Cannot read pbc from All.vcf (should be first line)!\n");
  }

  while (getc(vcf) != '\n')
    ;
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    printf("   box size: %lf x %lf x %lf\n\n", box_length.x, box_length.y, box_length.z);
  } //}}}

  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    fflush(stdout);
    printf("\rStep: %6d", ++count); //}}}

    char *stuff; // for the first two lines of a timestep ('# <number>\n t(imestep')
    stuff = malloc(128*sizeof(int));

    // read ordered timestep from All.vcf
    if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
      fprintf(stderr, "Cannot read coordinates from All.vcf! (%d. step; %d. bead)\n", count, test);
      exit(1);
    }
  }

  fflush(stdout);
  printf("\rLast Step: %6d\n", count);

  fclose(vcf);

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
  free(Molecule);
  free(saved_types);
  //}}}

  return 0;
}
