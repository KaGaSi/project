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
  fprintf(stderr, "   %s <input.vcf> <squares> <molecule(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <squares>         number of squares to add\n");
  fprintf(stderr, "   <output.vcf>      output vcf file with added squares\n");
  fprintf(stderr, "   <molecule(s)>     molecule names for squares to away from\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -st <int>      starting timestep for calculation\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
PlaceSquares adds square molecules into the box far away from specified molecules. \
The utility needs a vcf file with all beads in the simulation (such as All.vcf) \
and generates dl_meso CONFIG file. Since dl_meso requires molecules to be after \
unbonded beads, the utility substitutes the beads before the first molecule present \
in the initial system for the beads of the squares. If any of the substituted beads \
is charged, it exchanges this bead with the first uncharged bead in the system. \
The utility assumes neutral squares consisting of four 'P' beads -- if not so \
either change the output CONFIG file or change the square parameters inside \
this utility.\n\n");

      fprintf(stdout, "\
To use this utility, use a vcf file containing all beads from the simulations \
(i.e. All.vcf) along with vsf and FIELD files for that system. The utility will then \
generate the CONFIG file. FIELD file must be adapted manually (the square molecules) \
should be the first molecule type).\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input.vcf> <squares> <molecule(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input.vcf>       input filename (vcf format)\n");
      fprintf(stdout, "   <squares>         number of squares to add\n");
      fprintf(stdout, "   <output.vcf>      output vcf file with added squares\n");
      fprintf(stdout, "   <molecule(s)>     molecule names for squares to away from\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -st <int>      starting timestep for calculation\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 3; //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
//      strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      fprintf(stderr, "\nError: non-existent option '%s'\n\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < options) {
    fprintf(stderr, "\nError: too few mandatory arguments (%d instead of at least %d)\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
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

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf[32];
  strcpy(input_vcf, argv[++count]);

  // test if <input.vcf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(input_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "\nError: <input.vcf> '%s' does not have .vcf ending\n\n", input_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <squares> - number of squares to add //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "\nError: non-numeric argument for <squares>\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int number_of_squares = atoi(argv[count]); //}}}

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

  // <molecule names> - types of molecules for squares to avoid //{{{
  while (++count < argc && argv[count][0] != '-') {

    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);
    fprintf(stdout, "%s\n", argv[count]);

    if (mol_type == -1) {
      fprintf(stderr, "\nError: molecule '%s' does not exist in FIELD\n\n", argv[count]);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "\nError: cannot open file %s for reading\n\n", input_vcf);
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[128];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_vcf);
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
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(stdout, " %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // skip first start-1 steps //{{{
  count = 0;
  int test;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent) {
      if (script) {
        fprintf(stdout, "Discarding step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rDiscarding step: %6d", count);
      }
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "Error: premature end of %s file\n\n", input_vcf);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Discarded steps: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rDiscarded steps: %6d\n", count);
    }
  } //}}}
  //}}}

  // Square definition //{{{
//double bond = 0.7;
  // }}}

  // main loop //{{{
  count = 0; // count timesteps
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
      }
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
      exit(1);
    } //}}}


    // place squares
    for (int i = 0; i < number_of_squares; i++) {
      double closest;
      Vector position[16], distance[16];
      position[0].x = 13.692; //{{{
      position[0].y = 0.628;
      position[0].z = 1.270;
      position[1].x = 13.680;
      position[1].y = 0.831;
      position[1].z = 1.850;
      position[2].x = 14.356;
      position[2].y = 1.034;
      position[2].z = 1.947;
      position[3].x = 14.383;
      position[3].y = 0.810;
      position[3].z = 1.338;
      position[4].x = 13.257;
      position[4].y = 0.484;
      position[4].z = 0.710;
      position[5].x = 13.176;
      position[5].y = 0.907;
      position[5].z = 2.220;
      position[6].x = 14.955;
      position[6].y = 1.316;
      position[6].z = 2.265;
      position[7].x = 14.862;
      position[7].y = 0.885;
      position[7].z = 0.839;
      position[8].x = 14.423;
      position[8].y = 1.070;
      position[8].z = 1.386;
      position[9].x = 13.824;
      position[9].y = 1.235;
      position[9].z = 1.806;
      position[10].x = 13.502;
      position[10].y = 0.765;
      position[10].z = 1.562;
      position[11].x = 13.991;
      position[11].y = 0.578;
      position[11].z = 1.171;
      position[12].x = 15.014;
      position[12].y = 1.350;
      position[12].z = 1.072;
      position[13].x = 13.675;
      position[13].y = 1.807;
      position[13].z = 2.241;
      position[14].x = 12.994;
      position[14].y = 0.498;
      position[14].z = 1.628;
      position[15].x = 14.210;
      position[15].y = 0.060;
      position[15].z = 0.783; //}}}
      for (int j = 0; j < 16; j++) {
        distance[j].x = position[j].x - position[0].x;
        distance[j].y = position[j].y - position[0].y;
        distance[j].z = position[j].z - position[0].z;
      }
      do {
        closest = 1000;
        // generate random coordinate in a box
        Vector random;
        random.x = (double)rand() / (double)RAND_MAX * BoxLength.x;
        random.y = (double)rand() / (double)RAND_MAX * BoxLength.y;
        random.z = (double)rand() / (double)RAND_MAX * BoxLength.z;
        for (int j = 0; j < 16; j++) {
          position[j].x = random.x + distance[j].x;
          position[j].y = random.y + distance[j].y;
          position[j].z = random.z + distance[j].z;
        }
//      position[0].x = random.x;
//      position[0].y = random.y;
//      position[0].z = random.z;
//      position[1].x = position[0].x + bond;
//      position[1].y = position[0].y;
//      position[1].z = position[0].z;
//      position[2].x = position[0].x + bond;
//      position[2].y = position[0].y + bond;
//      position[2].z = position[0].z;
//      position[3].x = position[0].x;
//      position[3].y = position[0].y + bond;
//      position[3].z = position[0].z;
        // find if it is close to a specified molecule
        for (int j = 0; j < 16; j++) {
        printf("\nx%d %lf\n", i, closest);
          for (int k = 0; k < Counts.Molecules; k++) {
            int mol_type = Molecule[k].Type;
            if (MoleculeType[mol_type].Use) {
              for (int l = 0; l < MoleculeType[mol_type].nBeads; l++) {
                Vector dist = Distance(Bead[Molecule[k].Bead[l]].Position, position[j], BoxLength);
                dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                if (dist.x < closest) {
                  closest = dist.x;
                }
              }
            }
          }
        }
      } while (closest < 0.1);
      printf("closest=%lf\n", closest);

      int l = 0;
      for (int j = (Counts.Unbonded-1-16*i); j > (Counts.Unbonded-1-(16*i+16)); j--) {
        // change j's bead type to uncharged and change one uncharged unbonded bead into charged
        if (BeadType[Bead[j].Type].Charge != 0) {
          for (int k = 0; k < Counts.Unbonded; k++) {
            if (BeadType[Bead[k].Type].Charge == 0) {
              Bead[k].Type = Bead[j].Type;
              break;
            }
          }
        }
        if (((j-(Counts.Unbonded-1-16*i))%8) <= 3) {
          Bead[j].Type = 99; // some random high number to later recognise the bead as an unkown type
        } else {
          Bead[j].Type = 100; // some random high number to later recognise the bead as an unkown type

        }
        Bead[j].Position.x = position[15-l].x;
        Bead[j].Position.y = position[15-l].y;
        Bead[j].Position.z = position[15-l].z;
        l++;
        printf("%d %d %d\n", i, j, l);
      }
    }

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      fprintf(stdout, "\n%s", stuff);
    } //}}}
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  } //}}}

  // open output CONFIG file //{{{
  if ((vcf = fopen("CONFIG", "w")) == NULL) {
    fprintf(stderr, "\nError: cannot open file %s for reading\n\n", input_vcf);
    exit(1);
  } //}}}

  fprintf(vcf, "Generated via %s\n0 1\n", argv[0]);
  fprintf(vcf, "%lf 0 0\n", BoxLength.x);
  fprintf(vcf, "0 %lf 0\n", BoxLength.y);
  fprintf(vcf, "0 0 %lf\n", BoxLength.z);

  for (int i = 0; i < Counts.Beads; i++) {
    if (Bead[i].Type < 99) {
      fprintf(vcf, "%s %d\n", BeadType[Bead[i].Type].Name, i+1);
    } else if (Bead[i].Type == 99) {
      fprintf(vcf, "P %d\n", i+1);
    } else {
      fprintf(vcf, "Pq %d\n", i+1);
    }
    fprintf(vcf, "%7.3f %7.3f %7.3f\n", Bead[i].Position.x, Bead[i].Position.y, Bead[i].Position.z);
  }
  fclose(vcf);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff); //}}}

  return 0;
}
