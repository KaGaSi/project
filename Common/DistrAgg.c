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
  fprintf(stderr, "      -n <int>          start distribution calculation with <int>-th step\n");
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

      printf("Usage:\n");
      printf("   %s <input> <output distr file> <output avg file> <options>\n\n", argv[0]);

      printf("   <input>              input filename (agg format)\n");
      printf("   <output distr file>  filename with weight and number distributions\n");
      printf("   <output avg file>    filename with weight and number averages throughout simulation\n");
      printf("   <options>\n");
      printf("      -n <int>          start distribution calculation with <int>-th step\n");
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

  // read system information //{{{
  char vcf[1];
  vcf[0] = '\0';
  ReadStructure(vsf_file, vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file); //}}}

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
  double wdistr[Counts.Molecules];
  int ndistr[Counts.Molecules];
  int voldistr[Counts.Molecules];

  // molecule typs in aggregates: [agg size][mol type][number or SQR(number)]
  int ***molecules = malloc(Counts.Molecules*sizeof(int **));
  for (int i = 0; i < Counts.Molecules; i++) {
    molecules[i] = malloc(Counts.TypesOfMolecules*sizeof(int *));
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules[i][j] = calloc(2,sizeof(int));
    }
  }

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    wdistr[i] = 0;
    ndistr[i] = 0;
    voldistr[i] = 0;
  } //}}}

  // main loop //{{{
  int test, size_sqr = 0;
  count = 0;
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    // print (or not) step number{{{
    count++;
    if (!silent) {
      if (script) {
        printf("Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rStep: %6d", count);
      }
    } //}}}

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
      if (count >= start) { // start calculation of averages from specified 'start' timestep
        // distribution //{{{
        ndistr[Aggregate[i].nMolecules-1]++;
        wdistr[Aggregate[i].nMolecules-1] += Aggregate[i].nMolecules;
        voldistr[Aggregate[i].nMolecules-1] += Aggregate[i].nBeads;
        size_sqr += SQR(Aggregate[i].nMolecules); //}}}

        // number of various species in the aggregate //{{{
        int *mol_count = calloc(Counts.TypesOfMolecules,sizeof(int));
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          mol_count[Molecule[Aggregate[i].Molecule[j]].Type]++;
        }
        for (int j = 0; j < Counts.TypesOfMolecules; j++) {
          molecules[Aggregate[i].nMolecules-1][j][0] += mol_count[j];
          molecules[Aggregate[i].nMolecules-1][j][1] += SQR(mol_count[j]);
        } //}}}

        free(mol_count);
      }

      // average aggregation number //{{{
      if (!no_uni || Aggregate[i].nMolecules != 1) {
        aggs++;
        mols += Aggregate[i].nMolecules;

        avg_n += Aggregate[i].nMolecules;
        avg_w += SQR(Aggregate[i].nMolecules);
      } //}}}
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

  // number of species in agg file //{{{
  int mols = 0, beads = 0, mons = 0;
  for (int i = 0; i < Counts.Aggregates; i++) {
    mols += Aggregate[i].nMolecules;
    beads += Aggregate[i].nBeads;
    mons += Aggregate[i].nMonomers;
  }
  if (verbose) {
    printf("Number of species in provided .agg file: \n");
    printf("%10d molecules\n", mols);
    printf("%10d beads in molecules\n", beads);
    printf("%10d monomeric beads\n", mons);
  } //}}}

  // total number of Aggregates in the simulation //{{{
  int sum_agg = 0;
  for (int i = 0; i < Counts.Molecules; i++) {
    sum_agg += ndistr[i];
  } //}}}

  // print distributions to output file //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  // mass of all molecules
  double molecules_vol = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    molecules_vol += MoleculeType[i].Number * MoleculeType[i].nBeads;
  }

  count -= start -1;
  for (int i = 0; i < mols; i++) {
    fprintf(out, "%4d %lf %lf %lf\n", i+1, (double)(wdistr[i])/(mols*count), (double)(ndistr[i])/sum_agg, voldistr[i]/(count*molecules_vol));
  }
  fclose(out); //}}}

  // print numbers of various molecules types in the aggregates //{{{
  printf("1:sum");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf(" %d:%s", i+1, MoleculeType[i].Name);
  }
  putchar('\n');
  for (int i = 0; i < Counts.Molecules; i++) {
    if (ndistr[i] > 0) { // are there any aggregates of size 'i'?
      printf("%3d", i+1);
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        printf(" %7.3f", (double)(molecules[i][j][0])/ndistr[i]);
      }
      putchar('\n');
    }
  } //}}}

  // print overall averages //{{{
  // sum up total number (and SQR(number)) of each molecular species
  for (int i = 1; i < Counts.Molecules; i++) {
    ndistr[0] += ndistr[i]; // total number of aggregates
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules[0][j][0] += molecules[i][j][0];
      molecules[0][j][1] += molecules[i][j][1];
    }
  }
  // sum up total number (and SQR(number)) of all molecules
  for (int i = 1; i < Counts.Molecules; i++) {
    wdistr[0] += wdistr[i];
  }

  // open file with time evolution //{{{
  if ((out = fopen(output_avg, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_avg);
    exit(1);
  } //}}}

  // print legend (with column numbers)
  printf("1:<A_s>_w");
  fprintf(out, "# 1:<A_s>_w");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf(" %d:<%s>_w", i+1, MoleculeType[i].Name);
    fprintf(out, " %d:<%s>_w", i+1, MoleculeType[i].Name);
  }
  printf(" %d:<A_s>_n", Counts.TypesOfMolecules+1);
  fprintf(out, " %d:<A_s>_n", Counts.TypesOfMolecules+1);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf(" %d:<%s>_n", Counts.TypesOfMolecules+i+2, MoleculeType[i].Name);
    fprintf(out, " %d:<%s>_n", Counts.TypesOfMolecules+i+2, MoleculeType[i].Name);
  }
  putchar('\n');
  putc('\n', out);

  // print the averages
  printf("%7.3f", (double)(size_sqr)/wdistr[0]); // <A_s>_w
  fprintf(out, "# %7.3f", (double)(size_sqr)/wdistr[0]); // <A_s>_w
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf(" %7.3f", (double)(molecules[0][i][1])/molecules[0][i][0]); // <species>_w
    fprintf(out, " %7.3f", (double)(molecules[0][i][1])/molecules[0][i][0]); // <species>_w
  }
  printf("%7.3f", (double)(wdistr[0])/ndistr[0]); // <A_s>_n
  fprintf(out, "%7.3f", (double)(wdistr[0])/ndistr[0]); // <A_s>_n
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf(" %7.3f", (double)(molecules[0][i][0])/ndistr[0]); // <species>_n
    fprintf(out, " %7.3f", (double)(molecules[0][i][0])/ndistr[0]); // <species>_n
  }
  putc('\n', out);
  putchar('\n'); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.Molecules; i++) {
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      free(molecules[i][j]);
    }
    free(molecules[i]);
  }
  free(molecules); //}}}

  return 0;
}
