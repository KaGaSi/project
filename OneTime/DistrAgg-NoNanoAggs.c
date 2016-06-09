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
      CommonHelp(0);
      exit(0);
    }
  } //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 0; i < argc && argv[count][0] != '-'; i++) {
    count++;
  }

  if (argc < 4) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of 4)!\n\n", count);
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
  } //}}}

  // arrays for distribution //{{{
  double wdistr[Counts.Molecules];
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
  int NoNano = 0, NoNanoAgg = 0;
  double NoNanoMass = 0;
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    if (!silent) {
      fflush(stdout);
      printf("\rStep: %6d", ++count);
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
    for (int i = 0; i < Counts.Aggregates; i++) {
      bool Nano = false;

      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (strcmp("Nano", MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name) == 0) {
          Nano = true;
          break;
        }
      }

      if (!Nano) {
        NoNano += Aggregate[i].nMolecules;
        NoNanoAgg++;
        NoNanoMass += Aggregate[i].Mass;

        // distribution
        ndistr[Aggregate[i].nMolecules-1]++;
        wdistr[Aggregate[i].nMolecules-1] += Aggregate[i].nMolecules;
        voldistr[Aggregate[i].nMolecules-1] += Aggregate[i].Mass;

        // average aggregation number
        avg_n += Aggregate[i].nMolecules;
        avg_w += SQR(Aggregate[i].nMolecules);
      }
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

  if (!silent) {
    fflush(stdout);
    printf("\rLast Step: %6d\n", count);
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
    printf("%10d A4B2 molecules without contact with Nano molecules\n", NoNano);
    printf("           (%.2f per step)\n", (double)(NoNano)/count);
    printf("%10d number of aggregates with only A4B2 molecules\n", NoNanoAgg);
    printf("%lf mass of A4B2 molecules without contact with Nano molecules\n", NoNanoMass);
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

  double molecules_mass = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    molecules_mass += MoleculeType[i].Mass * MoleculeType[i].Number;
  }

  for (int i = 0; i < mols; i++) {
    fprintf(out, "%4d %lf %lf %lf\n", i+1, (double)(wdistr[i])/NoNano, (double)(ndistr[i])/NoNanoAgg, voldistr[i]/NoNanoMass);
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
