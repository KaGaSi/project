#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input> <output distr file> <options>\n\n", cmd);

  fprintf(stderr, "   <input>              input filename (agg format)\n");
  fprintf(stderr, "   <output distr file>  filename with weight and number distributions\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -n <int>          start with <int>-th step\n");
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
      printf("   %s <input> <output distr file> <options>\n\n", argv[0]);

      printf("   <input>              input filename (agg format)\n");
      printf("   <output distr file>  filename with weight and number distributions\n");
      printf("      -n <int>          start with <int>-th step\n");
      printf("   <options>\n");
      CommonHelp(0);
      exit(0);
    }
  } //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 0; i < argc && argv[count][0] != '-'; i++) {
    count++;
  }

  if (argc < 3) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of 4)!\n\n", count);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -n <int> option - number of starting timestep //{{{
  int start = 1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) {

      // Error - non-numeric argument
      if (argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      start = atoi(argv[i+1]);
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

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  char vcf[1];
  vcf[0] = '\0';
  ReadStructure(vsf_file, vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));

  for (int i = 0; i < Counts.Molecules; i++) {
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
  } //}}}

  // open output file and print first line //{{{
  FILE *out;
  if ((out = fopen(output_distr, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  fprintf(out, "# A_s  Absorbed D_n chains\n");
  fclose(out); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    printf("Since no coordinates are used, no structure information is available and therefore the data is for the whole simulated system!\n\n");
    char null[1] = {'\0'};
    VerboseOutput(verbose2, null, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // allocate & zeroize arrays //{{{
  int sum[Counts.Molecules], AbsorbedHomopolymers[Counts.Molecules];

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    sum[i] = 0;
    AbsorbedHomopolymers[i] = 0;
  } //}}}

  // main loop //{{{
  int test;
  count = 0;
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    count++;
    if (!silent) {
      if (script) {
        printf("Last Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rLast Step: %6d", count);
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

    if (count >= start) {
      // go through all aggregates
      for (int i = 0; i < Counts.Aggregates; i++) {

        // calculate number of charged polymers in aggregate //{{{
        int neutral = 0;
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          if (strcmp("Neutral", MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name) == 0) {
            neutral++;
          }
        }
        int charged = Aggregate[i].nMolecules - neutral; //}}}

        sum[charged]++;

        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          if (strcmp("Neutral", MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name) == 0) {
            AbsorbedHomopolymers[charged]++;
          }
        }
      }
    }
  }
  fclose(agg);

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d", count);
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

  // print results to output file //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  count -= start -1;
  for (int i = 0; i < mols; i++) {
    fprintf(out, "%4d %lf %d\n", i, (double)(AbsorbedHomopolymers[i])/sum[i], sum[i]);
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
