#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

#define ABS(x) (x>0?x:-x)

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <width> <output.rho> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <width>           width of a single bin\n");
  fprintf(stderr, "   <output.rho>      output density file (automatic ending 'molecule_name.rho' added)\n");
  fprintf(stderr, "   <molecule names>  molecule names to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -j             specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -n <average>   number of bins to average\n");
  fprintf(stderr, "      -c <int>       use <int>-th molecule bead instead of centre of mass\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("Just to calculate proper distances between beads in a nanoparticle.\n");
//    printf("MolDensity utility calculates number beads density for all bead"
//        "types from the centre of mass (or specified bead number in a"
//        "molecule) of specified molecules. Care must be taken with beadtype"
//        "names in various molecules types, because if one beadtype appears"
//        "in more molecule types, the resulting density for that beadtype"
//        "will be averaged without regard for the various types of molecule"
//        "it appears in.\n\n");

//    printf("The utility uses dl_meso.vsf (or other input structure file) and"
//        "FIELD (along with optional bond file) files to determine all"
//        "information about the system.\n\n");

//    printf("Usage:\n");
//    printf("   %s <input.vcf> <width> <output.rho> <molecule names> <options>\n\n", argv[0]);

//    printf("   <input.vcf>       input filename (vcf format)\n");
//    printf("   <width>           width of a single bin\n");
//    printf("   <output.rho>      output density file (automatic ending 'molecule_name.rho' added)\n");
//    printf("   <molecule names>  molecule names to calculate density for\n");
//    printf("   <options>\n");
//    printf("      -j             specify that aggregates with joined coordinates are used\n");
//    printf("      -n <average>   number of bins to average\n");
//    printf("      -c <int>       use <int>-th molecule bead instead of centre of mass\n");
//    CommonHelp(0);
//    exit(0);
    }
  }

  int options = 1; //}}}

//// test if options are given correctly //{{{
//for (int i = 1; i < argc; i++) {
//  if (argv[i][0] == '-' &&
//      strcmp(argv[i], "-i") != 0 &&
//      strcmp(argv[i], "-b") != 0 &&
//      strcmp(argv[i], "-v") != 0 &&
//      strcmp(argv[i], "-V") != 0 &&
//      strcmp(argv[i], "-s") != 0 &&
//      strcmp(argv[i], "-h") != 0 &&
//      strcmp(argv[i], "--script") != 0 &&
//      strcmp(argv[i], "-j") != 0 &&
//      strcmp(argv[i], "-n") != 0 &&
//      strcmp(argv[i], "-c") != 0) {

//    fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
//    ErrorHelp(argv[0]);
//    exit(1);
//  }
//} //}}}

  int count = 0;
  // check if correct number of arguments //{{{
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < options) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of at least %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // standard options //{{{
  char *vsf_file = calloc(32,sizeof(char *));
  char *bonds_file = calloc(32,sizeof(char *));
  bool verbose, verbose2, silent, script;
  bool error = CommonOptions(argc, argv, &vsf_file, &bonds_file, &verbose, &verbose2, &silent, &script);

  // was there error during CommonOptions()?
  if (error) {
    ErrorHelp(argv[0]);
    exit(1);
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

//// <width> - number of starting timestep //{{{
//// Error - non-numeric argument
//if (argv[++count][0] < '0' || argv[count][0] > '9') {
//  fprintf(stderr, "Non-numeric argement for <width>!\n");
//  ErrorHelp(argv[0]);
//  exit(1);
//}
//double width = atof(argv[count]); //}}}

//// <output.rho> - filename with bead densities //{{{
//char output_rho[16];
//strcpy(output_rho, argv[++count]); //}}}

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

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf);
    exit(1);
  } //}}}

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
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    printf("Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        printf(" %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int test;
  Vector coor[105];
  for (int i = 0; i < 105; i++ ) {
    coor[i].x = 10;
    coor[i].y = 10;
    coor[i].z = 10;
  }

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

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf, Counts, &Bead, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      }
    } //}}}

    // only one nanoparticle
    for (int i = 0; i < MoleculeType[0].nBeads; i++) {
      coor[i].x += Bead[Molecule[0].Bead[i]].Position.x - Bead[Molecule[0].Bead[0]].Position.x;
      coor[i].y += Bead[Molecule[0].Bead[i]].Position.y - Bead[Molecule[0].Bead[0]].Position.y;
      coor[i].z += Bead[Molecule[0].Bead[i]].Position.z - Bead[Molecule[0].Bead[0]].Position.z;
    }
    coor[0].x = 0;
    coor[0].y = 0;
    coor[0].z = 0;

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
    }
  } //}}}

  for (int i = 0; i < MoleculeType[0].nBeads; i++) {
    coor[i].x /= count;
    coor[i].y /= count;
    coor[i].z /= count;
  }

  // create bonds //{{{
  count = 0;
  double BondDist[MoleculeType[0].nBeads*(MoleculeType[0].nBeads-1)];
  FILE *write;
  write = fopen("stub", "w");
  for (int i = 1; i < MoleculeType[0].nBeads; i++) {
    for (int j = 0; j < i; j++) {

      BondDist[count] = sqrt(SQR(Bead[Molecule[0].Bead[i]].Position.x - Bead[Molecule[0].Bead[j]].Position.x) +
                             SQR(Bead[Molecule[0].Bead[i]].Position.y - Bead[Molecule[0].Bead[j]].Position.y) +
                             SQR(Bead[Molecule[0].Bead[i]].Position.z - Bead[Molecule[0].Bead[j]].Position.z));
      fprintf(write, "harm %3d %3d 5.000 %.3f\n", j+1, i+1, BondDist[count]);
      count++;
    }

//  double dist = sqrt(SQR(All[i][0]) + SQR(All[i][1]) + SQR(All[i][2]));
//  printf("%lf\n", dist);
  } //}}}

  fclose (write);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff); //}}}

  return 0;
}
