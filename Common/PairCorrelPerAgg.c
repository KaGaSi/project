#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input.agg> <width> <output.pcf> <bead type(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>    input filename (vcf format)\n");
  fprintf(stderr, "   <input.agg>    input filename with information about aggregates (agg format)\n");
  fprintf(stderr, "   <width>        width of a single bin\n");
  fprintf(stderr, "   <output.pcf>   output file with pair correlation function(s)\n");
  fprintf(stderr, "   <bead type(s)> bead type name(s) for pcf calculation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -n <int>    number of bins to average\n");
  fprintf(stderr, "      -st <int>   starting timestep for calculation\n");
  CommonHelp(1);
} //}}}
//TODO: somehow implement Join flag for bead types (throughout all utilities)

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("PairCorrelPerAgg utility calculates pair correlation function for specified \n");
      printf("bead types. The calculation is done per aggregates - that is only beads in  \n");
      printf("the same aggregate are used. If aggregate size(s) is not specified, average \n");
      printf("pcf is calculated (that is, regardless of aggregate size). All pairs of bead\n");
      printf("types (including same pair) are calculated - given A and B types, pcf       \n");
      printf("between A-A, A-B and B-B are calculated.                                    \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");
      printf("   %s <input.vcf> <input.agg> <width> <output.pcf> <bead type(s)> <options>\n\n", argv[0]);

      printf("   <input.vcf>     input filename (vcf format)\n");
      printf("   <input.agg>     input filename with information about aggregates (agg format)\n");
      printf("   <width>         width of a single bin\n");
      printf("   <output.pcf>    output file with pair correlation function(s)\n");
      printf("   <bead type(s)>  bead type name(s) for pcf calculation \n");
      printf("   <options>\n");
      printf("      -n <int>     number of bins to average\n");
      printf("      -st <int>    starting timestep for calculation\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int mandatory = 5; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < mandatory) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of at least %d)!\n\n", count, mandatory);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than dl_meso.vsf? //{{{
  char *vsf_file = calloc(32,sizeof(char *));
  if (VsfFileOption(argc, argv, &vsf_file)) {
    exit(1);
  } //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(32,sizeof(char *));
  if (BondsFileOption(argc, argv, &bonds_file)) {
    exit(0);
  } //}}}

  // output verbosity //{{{
  bool verbose, verbose2, silent;
  SilentOption(argc, argv, &verbose, &verbose2, &silent); // no output
  VerboseShortOption(argc, argv, &verbose); // verbose output
  VerboseLongOption(argc, argv, &verbose, &verbose2); // more verbose output
  bool script = BoolOption(argc, argv, "--script"); // do not use \r & co.
  // }}}
  //}}}

  // -n <int> option - number of bins to average //{{{
  int avg = 1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) {

      // Error - missing or non-numeric argument //{{{
      if ((i+1) >= argc) {
        fprintf(stderr, "Missing numeric argument for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }
      if (argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      } //}}}

      avg = atoi(argv[i+1]);
    }
  } //}}}

  // -st <int> option - number of starting timestep //{{{
  int start = 1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-st") == 0) {

      // Error - non-numeric argument //{{{
      if ((i+1) >= argc) {
        fprintf(stderr, "Missing numeric argument for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }
      if (argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      } //}}}

      start = atoi(argv[i+1]);
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

  // <input.agg> - filename of input file with aggregate information //{{{
  char input_agg[32];
  strcpy(input_agg, argv[++count]); //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <width>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.pcf> - filename with pcf(s) //{{{
  char output_pcf[32];
  strcpy(output_pcf, argv[++count]); //}}}

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

  // <type names> - names of bead types to use for closeness calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", argv[count], input_vcf);
      exit(1);
    }

    BeadType[type].Use = true;
  } //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  // open for the first time to read the distance and type names
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  // read minimum distance for closeness check (<distance> argument in Aggregates utility)
//double distance;
//fscanf(agg, "%*s %*s %lf", &distance);

//// skip <contacts> and <output.agg> in Aggregates command
//fscanf(agg, "%*s %*s");

//// read <type names> from Aggregates command //{{{
//int test;
//// reading ends if next argument (beginning with '-') or the following empty line is read
//while ((test = getc(agg)) != '-' && test != '\n') {
//  ungetc(test, agg);

//  char name[10];
//  fscanf(agg, "%s", name);
//  int type = FindBeadType(name, Counts, BeadType);

//  // Error - specified bead type name not in vcf input file
//  if (type == -1) {
//    fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", name, input_vcf);
//    exit(1);
//  }

//  BeadType[type].Use = true;

//  while ((test = getc(agg)) == ' ')
//    ;
//  ungetc(test, agg);
//} //}}}
  fclose(agg);

  // open again for production run - to ensure the pointer position in file is correct (at first 'Step')
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }
  // skip line with command to produce the agg file
  while (getc(agg) != '\n')
    ;
  // skip empty line
  while (getc(agg) != '\n')
    ; //}}}

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

  // write initial stuff to output pcf file //{{{
  FILE *out;
  if ((out = fopen(output_pcf, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_pcf);
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

  // number of bins - maximum distance is taken as three times box's body diagonal //{{{
  double max_dist = 3 * pow(SQR(BoxLength.x)+SQR(BoxLength.y)+SQR(BoxLength.z), 1.0/3.0);
  int bins = ceil(max_dist / width); //}}}

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
        printf("Discarding step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rDiscarding step: %6d", count);
      }
    } //}}}

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);
    if (SkipCoor(vcf, Counts, &stuff) == 1) {
      fprintf(stderr, "Premature end of %s file!\n", input_vcf);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      printf("Discarded steps: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rDiscarded steps: %6d\n", count);
    }
  } //}}}
  //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int aggs = 0; // number of aggregates
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent) {
      if (script) {
        printf("Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rStep: %6d", count);
      }
    } //}}}

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

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

    // calculate pair correlation function //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {
      aggs++;

      for (int j = 0; j < Aggregate[i].nBeads; j++) {
        if (BeadType[Bead[Aggregate[i].Bead[j]].Type].Use) {

          for (int k = (j+1); k < Aggregate[i].nBeads; k++) {
            if (BeadType[Bead[Aggregate[i].Bead[j]].Type].Use) {

              int bead1 = Aggregate[i].Bead[j];
              int bead2 = Aggregate[i].Bead[k];

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

              int l = rij.x / width;

              pcf[type1][type2][l]++;
            }
          }
        }
      }
    } //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}
  }
  fclose(vcf);
  fclose(agg);

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
    }
  } //}}}

  // write data to output file(s) //{{{
  if ((out = fopen(output_pcf, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", str);
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
  // normalisation factor - just sum up the pcf (for resulting pcf = 1) %{{{
  double **norm = malloc(Counts.TypesOfBeads*sizeof(double *));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    norm[i] = calloc(Counts.TypesOfBeads,sizeof(double));
  }
  for (int j = 0; j < bins; j++) {
    double shell = 4 * PI * CUBE(width) *(CUBE(j+1) - CUBE(j)) / 3;
    for (int k = 0; k < Counts.TypesOfBeads; k++) {
      for (int l = k; l < Counts.TypesOfBeads; l++) {
        if (BeadType[k].Use && BeadType[l].Use) {
          norm[k][l] += pcf[k][l][j] / shell;
        }
      }
    }
  } //}}}

  for (int j = 0; j < bins; j++) {

    // calculate volume of every shell that will be averaged
    double shell[avg];
    for (int k = 0; k < avg; k++) {
      shell[k] = 4 * PI * CUBE(width) *(CUBE(j+k+1) - CUBE(j+k)) / 3;
    }

    fprintf(out, "%8.5f", width*(j+0.5*avg));

    for (int k = 0; k < Counts.TypesOfBeads; k++) {
      for (int l = k; l < Counts.TypesOfBeads; l++) {
        if (BeadType[k].Use && BeadType[l].Use) {

          double temp = 0;

          // sump pcfs from all shells to be averaged
          for (int m = 0; m < avg; m++) {
            temp += pcf[k][l][j+m] / (shell[m] * norm[k][l]);
          }

          // print average value to output file
          fprintf(out, " %10f", temp/avg);
        }
      }
    }
    putc('\n',out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
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
