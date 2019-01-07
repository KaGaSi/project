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
  fprintf(stderr, "   %s <input> <width> <output> <bead type(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input>          input coordinate file (either vcf or vtf format)\n");
  fprintf(stderr, "   <width>          width of a single bin\n");
  fprintf(stderr, "   <output>         output file with pair correlation function(s)\n");
  fprintf(stderr, "   <bead type(s)>   bead type name(s) for pcf calculation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -n <int>      number of bins to average\n");
  fprintf(stderr, "      -st <int>     starting timestep for calculation\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, " \
PairCorrel utility calculates pair correlation function for specified \
bead types. All pairs of bead types (including same pair) are calculated - \
given A and B types, pcf between A-A, A-B and B-B are calculated.\n\n");

/*      fprintf(stdout, " \
The utility uses traject.vsf (or other input structure file) and FIELD \
(along with optional bond file) files to determine all information about \
the system.\n\n");
*/

      fprintf(stdout, "   %s <input> <width> <output> <bead type(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input>           input coordinate file (either vcf or vtf format)\n");
      fprintf(stdout, "   <width>           width of a single bin\n");
      fprintf(stdout, "   <output>          output file with pair correlation function(s)\n");
      fprintf(stdout, "   <bead type(s)>    bead type name(s) for pcf calculation \n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -n <int>       number of bins to average\n");
      fprintf(stdout, "      -st <int>      starting timestep for calculation\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
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
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than traject.vsf? //{{{
  char *input_vsf = calloc(32,sizeof(char *));
  if (FileOption(argc, argv, "-i", &input_vsf)) {
    exit(1);
  }
  if (input_vsf[0] == '\0') {
    strcpy(input_vsf, "traject.vsf");
  }

  // test if structure file ends with '.vsf'
  int ext = 2;
  char **extension;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_vsf, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(32,sizeof(char *));
  if (FileOption(argc, argv, "-b", &bonds_file)) {
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

  // number of bins to average //{{{
  int avg = 1;
  if (IntegerOption(argc, argv, "-n", &avg)) {
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

  // <input> - input coordinate file //{{{
  char input_coor[32];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  ext = 2;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(8*sizeof(char));
  }
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_coor, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double bin_width = atof(argv[count]); //}}}

  // <output> - filename with pcf(s) //{{{
  char output_pcf[32];
  strcpy(output_pcf, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <type names> - names of bead types to use for closeness calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\nError: bead type '%s' not present in %s file\n\n", argv[count], input_coor);
      exit(1);
    }

    BeadType[type].Use = true;
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "\nError: cannot read pbc from %s file\n\n", input_coor);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // write initial stuff to output pcf file //{{{
  FILE *out;
  if ((out = fopen(output_pcf, "w")) == NULL) {
    ErrorFileOpen(output_pcf, 'w');
    exit(1);
  }

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print bead type names to output file //{{{
  putc('#', out);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = i; j < Counts.TypesOfBeads; j++) {
      if (BeadType[i].Use && BeadType[j].Use) {
        fprintf(out, " %s-%s", BeadType[i].Name, BeadType[j].Name);
      }
    }
  }
  putc('\n', out); //}}}

  fclose(out); //}}}

  // number of bins - maximum distance is taken as half of the shortes BoxLength //{{{
  double max_dist = Min3(BoxLength.x, BoxLength.y, BoxLength.z) / 2;
  int bins = ceil(max_dist / bin_width); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

  // allocate memory for pcf arrays //{{{
  int *counter = calloc(Counts.TypesOfBeads,sizeof(int *)); // to count number of pairs
  double ***pcf = malloc(Counts.TypesOfBeads*sizeof(double **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    pcf[i] = malloc(Counts.TypesOfBeads*sizeof(double *));
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      pcf[i][j] = calloc(bins,sizeof(double));
    }
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_coor, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
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
      fprintf(stderr, "\nError: premature end of %s file\n\n", input_coor);
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

  // main loop //{{{
  count = 0; // count timesteps
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
      }
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // calculate pair correlation function //{{{
    for (int j = 0; j < (Counts.Bonded+Counts.Unbonded); j++) {
      if (BeadType[Bead[j].Type].Use) {

        for (int k = (j+1); k < (Counts.Bonded+Counts.Unbonded); k++) {
          if (BeadType[Bead[k].Type].Use) {

            int bead1 = j;
            int bead2 = k;

            int type1 = Bead[bead1].Type;
            int type2 = Bead[bead2].Type;

            // type1 shouldn't be larger then type2 //{{{
            if (type1 > type2) {
              int temp = type1;
              type1 = type2;
              type2 = temp;

              temp = bead1;
              bead1 = bead2;
              bead2 = temp;
            } //}}}

            counter[type2]++;

            // distance between bead1 and bead2
            Vector rij = Distance(Bead[bead1].Position, Bead[bead2].Position, BoxLength);
            rij.x = sqrt(SQR(rij.x) + SQR(rij.y) + SQR(rij.z));

            // count only distances up to half of the shortest box length
            if (rij.x < max_dist) {
              int l = rij.x / bin_width;
              pcf[type1][type2][l]++;
            }
          }
        }
      }
    } //}}}

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

  // write data to output file(s) //{{{
  if ((out = fopen(output_pcf, "a")) == NULL) {
    ErrorFileOpen(output_pcf, 'a');
    exit(1);
  }

  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    counter[0] = 0;
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Use) {
      counter[0] += BeadType[i].Number;
    }
  }

  // calculate pcf
  for (int j = 1; j < (bins-avg); j++) {

    // calculate volume of every shell that will be averaged
    double shell[avg];
    for (int k = 0; k < avg; k++) {
      shell[k] = 4.0 / 3 * PI * CUBE(bin_width) * (CUBE(j+k+1) - CUBE(j+k));
    }

    fprintf(out, "%8.5f", bin_width*(j+0.5*avg));

    double volume = BoxLength.x * BoxLength.y * BoxLength.z;
    for (int k = 0; k < Counts.TypesOfBeads; k++) {
      for (int l = k; l < Counts.TypesOfBeads; l++) {
        if (BeadType[k].Use && BeadType[l].Use) {

          double temp = 0; // for normalisation

          // sump pcfs from all shells to be averaged
          for (int m = 0; m < avg; m++) {
            double pairs;
            if (k == l) {
              pairs = ((SQR(BeadType[k].Number) - BeadType[k].Number)) / 2;
            } else {
              pairs = BeadType[k].Number * BeadType[l].Number;
            }
            // for normalisation - WRONG
            double pair_den = volume / pairs;
            double norm_factor = pair_den / shell[m];
            temp += pcf[k][l][j+m] * norm_factor;
          }

          // print average value to output file
          fprintf(out, " %10f", temp/avg);
        }
      }
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      free(pcf[i][j]);
    }
    free(pcf[i]);
  }
  free(pcf);
  free(counter); //}}}

  return 0;
}
