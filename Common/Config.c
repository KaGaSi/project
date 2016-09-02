#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -n <step>      timestep for creating CONFIG (default: last)\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("Config utility generates CONFIG file from given step of a vcf file.         \n");
      printf("Coordinate file needs to contain all beads in the simulation for to work    \n");
      printf("with dl_meso.                                                               \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <options>\n\n", argv[0]);

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <options>\n");
      printf("      -n <step>      timestep for creating CONFIG (default: last)\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 1; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < options) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // standard options //{{{
  char *vsf_file = calloc(32,sizeof(char *));
  char *bonds_file = calloc(32,sizeof(char *));
  bool verbose, verbose2, silent;
  bool error = CommonOptions(argc, argv, &vsf_file, &bonds_file, &verbose, &verbose2, &silent);

  // was there error during CommonOptions()?
  if (error) {
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -n <step> - choose timestep to create CONFIG file from //{{{
  int timestep = -1;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) {

      // Error - non-numeric argument
      if (argv[i+1][0] <= '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric (or non-positive numeric) argement for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      timestep = atoi(argv[i+1]);
    }
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      printf(" %s", argv[i]);
    printf("\n\n");
  } //}}}

  count = 0; // count mandatory arguments

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

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file);

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf);
    exit(1);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read string from '%s' file!\n", input_vcf);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_vcf);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;
  // skip blank line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    printf("   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = calloc(128,sizeof(int)); //}}}

  // main loop //{{{
  fpos_t pos; // for saving pointer position in vcf file
  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF && count != timestep) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      fflush(stdout);
      printf("\rStep: %6d", count);
    }

    // save pointer position in file
    fgetpos(vcf, &pos);

    SkipCoor(vcf, Counts, &stuff);

    // if -V option used, print comment from the beginning of a timestep
    if (verbose2)
      printf("\n%s", stuff);
  }

  // restore pointer position in FIELD file
  fsetpos(vcf, &pos);

  // read indexed timestep from input .vcf file //{{{
  if (indexed) {
    if ((test = ReadCoorIndexed(vcf, Counts, &Bead, &stuff)) != 0) {
      fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
      exit(1);
    } //}}}
  // or read ordered timestep from input .vcf file //{{{
  } else {
    if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
      fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
      exit(1);
    }
  } //}}}

  if (!silent) {
    fflush(stdout);
    printf("\rConfig Step: %6d\n", count);
  }

  fclose(vcf); //}}}

  // create CONFIG //{{{
  // open output CONFIG file for writing //{{{
  FILE *out;
  if ((out = fopen("CONFIG", "w")) == NULL) {
    fprintf(stderr, "Cannot open file CONFIG!\n");
    exit(1);
  } //}}}

  // print CONFIG file initial stuff //{{{
  fprintf(out, "NAME\n       0       1\n");
  fprintf(out, "%lf 0.000000 0.000000\n", BoxLength.x);
  fprintf(out, "0.000000 %lf 0.000000\n", BoxLength.y);
  fprintf(out, "0.000000 0.000000 %lf\n", BoxLength.z); //}}}

  // coordinates of unbonded beads //{{{
  for (int i = 0; i < Counts.Unbonded; i++) {
    fprintf(out, "%s %d\n", BeadType[Bead[i].Type].Name, i+1);
    fprintf(out, "%lf %lf %lf\n", Bead[i].Position.x,
                                  Bead[i].Position.y,
                                  Bead[i].Position.z);
  } //}}}

  // coordinates of bonded beads //{{{
  for (int i = Counts.Unbonded; i < (Counts.Unbonded+Counts.Bonded); i++) {
    fprintf(out, "%s %d\n", BeadType[Bead[i].Type].Name, i+1);
    fprintf(out, "%lf %lf %lf\n", Bead[i].Position.x,
                                  Bead[i].Position.y,
                                  Bead[i].Position.z);
  } //}}}

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  //}}}

  return 0;
}
