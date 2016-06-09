#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input.agg> <output.rho> <agg sizes> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <input.agg>         input filename with information about aggregates (agg format)\n");
  fprintf(stderr, "   <output.rho>        output density file (automatic ending 'agg#.rho' added)\n");
  fprintf(stderr, "   <agg sizes>         aggregate sizes to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -j               specify that aggregates with joined coordinates are used\n");
  CommonHelp(1);
} //}}}


Vector Gyration(int n, int *list, Counts Counts, Vector BoxLength,
                BeadType *BeadType, Bead *Bead) {

  // gyration tensor (3x3 array) //{{{
  struct Tensor {
    Vector x, y, z;
  } GyrationTensor;

  GyrationTensor.x.x = 0;
  GyrationTensor.x.y = 0;
  GyrationTensor.x.z = 0;
  GyrationTensor.y.x = 0;
  GyrationTensor.y.y = 0;
  GyrationTensor.y.z = 0;
  GyrationTensor.z.x = 0;
  GyrationTensor.z.y = 0;
  GyrationTensor.z.z = 0; //}}}

  Vector com = CenterOfMass(n, list, Bead, BeadType);

  // move center of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    Bead[list[i]].Position.x =- com.x;
    Bead[list[i]].Position.y =- com.y;
    Bead[list[i]].Position.z =- com.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 1; i < n; i++) {
    GyrationTensor.x.x += Bead[list[i]].Position.x * Bead[list[i]].Position.x;
    GyrationTensor.x.y += Bead[list[i]].Position.x * Bead[list[i]].Position.y;
    GyrationTensor.x.z += Bead[list[i]].Position.x * Bead[list[i]].Position.z;
    GyrationTensor.y.x += Bead[list[i]].Position.y * Bead[list[i]].Position.x;
    GyrationTensor.y.y += Bead[list[i]].Position.y * Bead[list[i]].Position.y;
    GyrationTensor.y.z += Bead[list[i]].Position.y * Bead[list[i]].Position.z;
    GyrationTensor.z.x += Bead[list[i]].Position.z * Bead[list[i]].Position.x;
    GyrationTensor.z.y += Bead[list[i]].Position.z * Bead[list[i]].Position.y;
    GyrationTensor.z.z += Bead[list[i]].Position.z * Bead[list[i]].Position.z;
  }
  GyrationTensor.x.x /= n;
  GyrationTensor.x.y /= n;
  GyrationTensor.x.z /= n;
  GyrationTensor.y.x /= n;
  GyrationTensor.y.y /= n;
  GyrationTensor.y.z /= n;
  GyrationTensor.z.x /= n;
  GyrationTensor.z.y /= n;
  GyrationTensor.z.z /= n; //}}}

  // calculate tensor's eigenvalues //{{{
  double a = GyrationTensor.x.x + GyrationTensor.y.y + GyrationTensor.z.z;
  double b = - GyrationTensor.x.x * GyrationTensor.y.y
             - GyrationTensor.x.x * GyrationTensor.z.z
             - GyrationTensor.y.y * GyrationTensor.z.z
             + GyrationTensor.x.z * GyrationTensor.z.x
             + GyrationTensor.x.y * GyrationTensor.y.x
             + GyrationTensor.y.z * GyrationTensor.z.y;
  double c = GyrationTensor.x.x + GyrationTensor.y.y + GyrationTensor.z.z
           + GyrationTensor.y.x + GyrationTensor.z.y + GyrationTensor.x.z
           + GyrationTensor.z.x + GyrationTensor.x.y + GyrationTensor.y.z
           - GyrationTensor.z.x + GyrationTensor.x.z + GyrationTensor.y.y
           - GyrationTensor.x.y + GyrationTensor.y.x + GyrationTensor.z.z
           - GyrationTensor.z.y + GyrationTensor.y.z + GyrationTensor.x.x; //}}}

  Vector eigen;
  eigen.x = a;
  eigen.y = b;
  eigen.z = c;

  return (eigen);
}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("                                                                            \n");
      printf("                                                                            \n");
      printf("                                                                            \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <input.agg> <output> <agg sizes> <options>\n\n", argv[0]);

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <input.agg>         input filename with information about aggregates (agg format)\n");
      printf("   <output>            output file with radii of gyration\n");
      printf("   <agg sizes>         aggregate sizes to calculate radius of gyration for\n");
      printf("   <options>\n");
      printf("      -j               specify that aggregates with joined coordinates are used\n");
      CommonHelp(0);
      exit(0);
    }
  } //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 0; i < argc && argv[count][0] != '-'; i++) {
    count++;
  }

  if (argc < 5) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of at least 5)!\n\n", count);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -j option - coordinates are joined //{{{
  bool joined = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      joined = true;
    }
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

  // <output> - filename with radii of gyration //{{{
  char output_rho[16];
  strcpy(output_rho, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
  }

  int aggs = 0;

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (argv[count][0] < '1' || argv[count][0] > '9') {
      fprintf(stderr, "Non-numeric option in <agg sizes>!\n");
      exit(1);
    } //}}}

    agg_sizes[aggs][0] = atoi(argv[count]);

    // write initial stuff to output density file //{{{
    FILE *out;
    char str[32];

    sprintf(str, "%s%d.rho", output_rho, agg_sizes[aggs][0]);
    if ((out = fopen(str, "w")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", str);
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
      fprintf(out, " %s", BeadType[i].Name);
    }
    putc('\n', out); //}}}

    fclose(out); //}}}

    aggs++; // number of aggregate sizes
  } //}}}

  // write initial stuff to output file //{{{
  FILE *out;
  if ((out = fopen(output_rho, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_rho);
    exit(1);
  }

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print agg sizes to output file //{{{
  putc('#', out);
  for (int i = 0; i < aggs; i++) {
    fprintf(out, " %10d", agg_sizes[i][0]);
  }
  putc('\n', out); //}}}

  fclose(out); //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  // read minimum distance for closeness check (<distance> argument in Aggregates utility)
  double distance;
  fscanf(agg, "%*s %*s %lf", &distance);

  // skip <contacts> and <output.agg> in Aggregates command
  fscanf(agg, "%*s %*s");

  // read <type names> from Aggregates command //{{{
  int test;
  while ((test = getc(agg)) != '-') {
    ungetc(test, agg);

    char name[10];
    fscanf(agg, "%s ", name);
    int type = FindBeadType(name, Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", name, input_vcf);
      exit(1);
    }

    BeadType[type].Use = true;
  } //}}}

  while (getc(agg) != '\n')
    ;
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

    printf("Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      printf(" %d", agg_sizes[i][0]);
    }
    putchar('\n');
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    if (!silent) {
      fflush(stdout);
      printf("\rStep: %6d", ++count);
    }

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf, Counts, &Bead, &stuff)) != 0) {
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      }
    } //}}}

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate arrays for the timestep //{{{
    int *agg_counts = calloc(aggs,sizeof(int));
    double *Rg = calloc(aggs,sizeof(int)); //}}}

    for (int i = 0; i < Counts.Aggregates; i++) {

      // test if aggregate is of correct size //{{{
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j][0] == Aggregate[i].nMolecules) {
          correct_size = j;
        }
      } //}}}

      if (correct_size != -1) {
        agg_counts[correct_size]++;

        Vector eigen = Gyration(Aggregate[i].nBeads, Aggregate[i].Bead, Counts, BoxLength, BeadType, Bead);

        Rg[correct_size] += sqrt(eigen.x + eigen.y + eigen.z);

        agg_sizes[correct_size][1]++;
      }
    }

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}

    free(agg_counts);
    free(Rg);
  }
  fclose(vcf);
  fclose(agg);

  if (!silent) {
    fflush(stdout);
    printf("\rLast Step: %6d\n", count);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(agg_sizes);
  free(stuff); //}}}

  return 0;
}
