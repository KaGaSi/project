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
  fprintf(stderr, "      -i <name>         use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>         file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v                verbose output\n");
  fprintf(stderr, "      -V                verbose output with more information\n");
  fprintf(stderr, "      -h                print this help and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("DistrAgg calculates weight and number average aggregation numbers during    \n");
      printf("the simulation run as well as overall weight and number distributions.      \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input> <output distr file> <output avg file> <options>\n\n", argv[0]);

      printf("   <input>              input filename (agg format)\n");
      printf("   <output distr file>  filename with weight and number distributions\n");
      printf("   <output avg file>    filename with weight and number averages throughout simulation\n");
      printf("   <options>\n");
      printf("      -i <name>         use input .vsf file different from dl_meso.vsf\n");
      printf("      -b <name>         file containing bond alternatives to FIELD\n");
      printf("      -v                verbose output\n");
      printf("      -V                verbose output with more information\n");
      printf("      -h                print this help and exit\n");
      exit(0);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  printf("\n\n"); //}}}

  // check if correct number of arguments //{{{
  if (argc < 4) {
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

  // -V option - verbose output with comments from input .vcf file //{{{
  bool verbose2 = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-V") == 0) {
      verbose = true;
      verbose2 = true;

      break;
    }
  } //}}}

  int count = 0; // count mandatory arguments

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

  fprintf(out, "# step  wdistr  ndistr\n");
  fclose(out);

  if ((out = fopen(output_avg, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_avg);
    exit(1);
  }

  fprintf(out, "# step  w-avg  n-avg\n");
  fclose(out); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    char null[1] = {'\0'};
    VerboseOutput(verbose2, null, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // arrays for distribution //{{{
  double wdistr[Counts.Molecules];
  int ndistr[Counts.Molecules];

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    wdistr[i] = 0;
    ndistr[i] = 0;
  } //}}}

  // main loop //{{{
  int test;
  count = 0;
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    fflush(stdout);
    printf("\rStep: %6d", ++count);

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
    for (int i = 0; i < Counts.Aggregates; i++) {
      // distribution
      ndistr[Aggregate[i].nMolecules-1]++;
      wdistr[Aggregate[i].nMolecules-1] += Aggregate[i].nMolecules;

      // average aggregation number
      avg_n += Aggregate[i].nMolecules;
      avg_w += SQR(Aggregate[i].nMolecules);
    }

    // print averages to output file //{{{
    if ((out = fopen(output_avg, "a")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", output_avg);
      exit(1);
    }

    fprintf(out, "%5d %lf %lf\n", count, avg_w/Counts.Molecules, avg_n/Counts.Aggregates);
    fclose(out); //}}}
  }
  fclose(agg);

  fflush(stdout);
  printf("\rLast Step: %6d\n", count); //}}}

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

  for (int i = 0; i < Counts.Molecules; i++) {
    fprintf(out, "%4d %lf %lf\n", i+1, wdistr[i]/(Counts.Molecules*count), (double)(ndistr[i])/sum_agg);
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
