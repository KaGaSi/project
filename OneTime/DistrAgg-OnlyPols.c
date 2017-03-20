#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input> <output distr file> <output avg file> <options>\n\n", cmd);

  fprintf(stderr, "   <input>              input filename (agg format)\n");
  fprintf(stderr, "   <output distr file>  filename with weight and number distributions\n");
  fprintf(stderr, "   <output avg file>    filename with weight and number averages throughout simulation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -n <int>          start weight distribution calculation with <int>-th step\n");
  fprintf(stderr, "      --no-unimers      do not count unimers into averages\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("DistrAgg calculates weight and number average aggregation numbers during    \n");
      printf("the simulation run as well as overall weight and number distributions and   \n");
      printf("volume fractions of aggregates.                                           \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("This version calculates distributions of aggregate that don't contain       \n");
      printf("nanoparticles (they contain only polymer chains - done for A4B30 chains).   \n");
      printf("It also prints ids of all molecules in aggregates with a single 'nano'.     \n\n");

      printf("Usage:\n");
      printf("   %s <input> <output distr file> <output avg file> <options>\n\n", argv[0]);

      printf("   <input>              input filename (agg format)\n");
      printf("   <output distr file>  filename with weight and number distributions\n");
      printf("   <output avg file>    filename with weight and number averages throughout simulation\n");
      printf("   <options>\n");
      printf("      -n <int>          start weight distribution calculation with <int>-th step\n");
      printf("      --no-unimers      do not count unimers into averages\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 3; //}}}

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
        strcmp(argv[i], "--no-unimers") != 0) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
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
    fprintf(stderr, "Too little mandatory arguments (%d instead of %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -n <int> option - number of starting timestep //{{{
  int start = 1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) {

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

  // --no-unimers - do not count unimers towards averages //{{{
  bool no_uni = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--no-unimers") == 0) {
      no_uni = true;
    }
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

  // <input> - filename of input agg file //{{{
  char input_agg[32];
  strcpy(input_agg, argv[++count]); //}}}

  // open input file and skip the first two lines //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // <output distr file> - filename with weight and number distributions //{{{
  char output_distr[32];
  strcpy(output_distr, argv[++count]); //}}}

  // <output avg file> - filename with weight and number average aggregation numbers //{{{
  char output_avg[32];
  strcpy(output_avg, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  ReadStructure(vsf_file, '\0', bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file);

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));

  for (int i = 0; i < Counts.Molecules; i++) {
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
  } //}}}

  // open output files and print first line //{{{
  FILE *out;
  if ((out = fopen(output_distr, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  fprintf(out, "# A_s  wdistr  ndistr voldistr\n");
  fclose(out);

  if ((out = fopen(output_avg, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_avg);
    exit(1);
  }

  fprintf(out, "# step  w-avg  n-avg\n");
  fclose(out); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    printf("Since no coordinates are used, no structure information is available and therefore the data is for the whole simulated system!\n\n");
    char null[1] = {'\0'};
    VerboseOutput(verbose2, null, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // arrays for distribution //{{{
  int wdistr[Counts.Molecules];
  int ndistr[Counts.Molecules];
  int voldistr[Counts.Molecules];

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    wdistr[i] = 0;
    ndistr[i] = 0;
    voldistr[i] = 0;
  } //}}}

  // main loop //{{{
  int test;
  count = 0;
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    count++;
    if (!silent) {
      if (script) {
        printf("Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rStep: %6d", count);
      }
    }

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

    // print info about aggregates if '-V' is used //{{{
    if (verbose2) {
      for (int i = 0; i < Counts.Aggregates; i++) {
        printf("\nAggregate[%3d].{Mass = %6.2f,\nnMolecules = %3d:", i+1, Aggregate[i].Mass, Aggregate[i].nMolecules);
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          printf(" %d", Aggregate[i].Molecule[j]+1);
        }
        printf(",\n nBeads = %4d:", Aggregate[i].nBeads);
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          printf(" %d", Aggregate[i].Bead[j]);
        }
        printf(",\n nMonomers = %4d:", Aggregate[i].nMonomers);
        for (int j = 0; j < Aggregate[i].nMonomers; j++) {
          printf(" %d", Aggregate[i].Monomer[j]);
        }
        printf("}\n");
      }
      putchar('\n');
    } //}}}

    // go through all aggregates
    double avg_n = 0,
           avg_w = 0;
    int aggs = 0, // number of aggregates (w/o unimers if --no-unimers)
        mols = 0; // number of molecules (w/o those in unimers if --no-unimers)
    for (int i = 0; i < Counts.Aggregates; i++) {
      // check if aggregate contains Nano molecule type
      bool nano = false;
      int nanos = 0;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (strcmp(MoleculeType[mol_type].Name, "Nano") == 0) {
          nano = true;
          nanos++;
//        break;
        }
      }

//    // print all molecule ids for aggregates containing 1 nanoparticle
//    if (nanos == 1 && count >= start) {
//      putchar('\n');
//      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
//        printf("%d ", Aggregate[i].Molecule[j]+1);
//      }
//      putchar('\n');
//    }

      // distribution if no Nano in aggregate
      if (!nano && count >= start) {
        ndistr[Aggregate[i].nMolecules-1]++;
        wdistr[Aggregate[i].nMolecules-1] += Aggregate[i].nMolecules;
        voldistr[Aggregate[i].nMolecules-1] += Aggregate[i].nBeads;
      }

      // average aggregation number if no Nano in aggregate
      if (!nano && (!no_uni || Aggregate[i].nMolecules != 1)) {
        aggs++;
        mols += Aggregate[i].nMolecules;

        avg_n += Aggregate[i].nMolecules;
        avg_w += SQR(Aggregate[i].nMolecules);
      }
    }

    // print averages to output file //{{{
    if ((out = fopen(output_avg, "a")) == NULL) {
      // print newline to stdout if Step... doesn't end with one
      if (!script && !silent) {
        putchar('\n');
      }
      fprintf(stderr, "Cannot open file %s!\n", output_avg);
      exit(1);
    }

    fprintf(out, "%5d %lf %lf\n", count, avg_w/mols, avg_n/aggs);
    fclose(out); //}}}
  }
  fclose(agg);

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
    }
  } //}}}
count -= (start - 1);
printf("count = %d\n", count);

  // total number of Aggregates in the simulation //{{{
  int sum_agg = 0, sum_agg_no_uni = 0, molecules_num = 0, molecules_vol = 0;
  for (int i = 0; i < Counts.Molecules; i++) {
//  if (i < 30) {
//    printf("ndistr[%d] = %3d ", i, ndistr[i]);
//    printf("wdistr[%d] = %d\n", i, wdistr[i]);
//  }
    sum_agg += ndistr[i];
    if (i > 0) {
      sum_agg_no_uni += ndistr[i];
    }
    molecules_num += wdistr[i];
    molecules_vol += voldistr[i];
  } //}}}
//  printf("%d %d %d", sum_agg, sum_agg_no_uni, count);

  // number of species in agg file //{{{
  int mols = 0, beads = 0, mons = 0;
  for (int i = 0; i < Counts.Aggregates; i++) {
    mols += Aggregate[i].nMolecules;
    beads += Aggregate[i].nBeads;
    mons += Aggregate[i].nMonomers;
  }
  if (verbose) {
    printf("\nNumber of species in provided .agg file: \n");
    printf("%10d molecules\n", mols);
    printf("%10d beads in molecules\n", beads);
    printf("%10d monomeric beads\n", mons);
    printf("Average number of chains not in contact with 'Nano' molecule: %lf\n", (double)(molecules_num)/count);
    printf("Average number of aggregates without 'Nano' molecule: %lf\n", (double)(sum_agg)/count);
    printf("Average number of polymer unimers: %lf\n", (double)(sum_agg - sum_agg_no_uni)/count);
  } //}}}

  // print distributions to output file //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  count -= start -1;
  for (int i = 0; i < mols; i++) {
    fprintf(out, "%4d %lf %lf %lf\n", i+1, (double)(wdistr[i])/molecules_num, (double)(ndistr[i])/sum_agg, (double)(voldistr[i])/molecules_vol);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead); //}}}

  return 0;
}
