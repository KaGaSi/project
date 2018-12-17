#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <output file> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <output file>     name of output file with end-to-end distances\n");
  fprintf(stderr, "   <molecule names>  names of molecule type(s) to use for calculation\n");
  fprintf(stderr, "   <options>\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("\
EndToEnd utility calculates end to end distance for linear chains (no check \
whether the molecules are linear is performed). It calculates distance \
between first and last bead in a molecule.\n\n");

      printf("\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n"); \

      printf("Usage:\n");
      printf("   %s <input.vcf> <output file> <molecule names> <options>\n\n", argv[0]);

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <output file>     name of output file with end-to-end distances\n");
      printf("   <molecule names>  names of molecule type(s) to use for calculation\n");
      printf("   <options>\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count][0] != '-'; i++) {
    count++;
  }

  if (argc < req_args) {
    ErrorArgNumber(count, req_args);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
//      strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than dl_meso.vsf? //{{{
  char *input_vsf = calloc(32,sizeof(char *));
  if (VsfFileOption(argc, argv, &input_vsf)) {
    exit(1);
  } //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(32,sizeof(char *));
  if (BondsFileOption(argc, argv, &bonds_file)) {
    exit(0);
  } //}}}

  // output verbosity //{{{
  bool verbose2, silent;
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  VerboseLongOption(argc, argv, &verbose, &verbose2); // more verbose output
  SilentOption(argc, argv, &verbose, &verbose2, &silent); // no output
  bool script = BoolOption(argc, argv, "--script"); // do not use \r & co.
  // }}}
  //}}}

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
    ErrorExtension(input_vcf, ".vcf");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <output> - file name with end-to-end distances //{{{
  char output[32];
  strcpy(output, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <molecule names> - names of molecule types to use //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (type == -1) {
      fprintf(stderr, "\nError: molecule type '%s' is not in %s file\n\n", argv[count], input_vcf);
      exit(1);
    }

    MoleculeType[type].Use = true;
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open output file and print molecule names //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    fprintf(stderr, "\nError: cannot open file %s for writing\n\n", output);
    exit(1);
  }

  fprintf(out, "# timestep ");

  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, "%10s", MoleculeType[i].Name);
    }
  }
  putc('\n', out);

  fclose(out); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "\nError: cannot open file %s for reading\n\n", input_vcf);
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_vcf);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "\nError: cannot read pbc from %s\n\n", input_vcf);
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
  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      if (script) {
        printf("Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rStep: %6d", count);
      }
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // join all molecules
    RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

    // open output file for appending //{{{
    if ((out = fopen(output, "a")) == NULL) {
      fprintf(stderr, "\nError: cannot open file %s for appending\n\n", output);
      exit(1);
    } //}}}

    // calculate & write end-to-end distance //{{{
    double Re[Counts.TypesOfMolecules];
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      Re[i] = 0;
    }

    // go through all molecules
    for (int i = 0; i < Counts.Molecules; i++) {
      if (MoleculeType[Molecule[i].Type].Use) {
        int id1 = Molecule[i].Bead[0];
        int id2 = Molecule[i].Bead[MoleculeType[Molecule[i].Type].nBeads-1];

        // calculate distance between first and last bead in molecule 'i'
        Vector dist;
        dist.x = Bead[id1].Position.x - Bead[id2].Position.x;
        dist.y = Bead[id1].Position.y - Bead[id2].Position.y;
        dist.z = Bead[id1].Position.z - Bead[id2].Position.z;

        dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

        Re[Molecule[i].Type] += dist.x;
      }
    }

    // write to output file and add to global averages
    fprintf(out, "%6d", count);
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, "%10f", Re[i]/MoleculeType[i].Number);
      }
    }
    putc('\n', out); //}}}

    fclose(out);

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      printf("\n%s", stuff);
  }

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
    }
  }

  fclose(vcf); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(bonds_file);
  //}}}

  return 0;
}
