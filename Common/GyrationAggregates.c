#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input.agg> <output> <agg sizes> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <input.agg>         input filename with information about aggregates (agg format)\n");
  fprintf(stderr, "   <output>            output file with data during simulation run\n");
  fprintf(stderr, "   <agg sizes>         aggregate sizes to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined         specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -bt              specify bead types to be used for calculation (default is all)\n");
  fprintf(stderr, "      -m <name>        agg size means number of <name> molecule types in an aggregate\n");
  fprintf(stderr, "      --no-unimers      do not count unimers into averages\n");
  CommonHelp(1);
} //}}}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void jacobi(double **a, int n, double d[], double **v, int *nrot) { //{{{

  double b[3], z[3];

  for (int ip = 0; ip < n; ip++) {
    for (int iq = 0; iq < n; iq++) {
      v[ip][iq] = 0;
    }
    v[ip][ip] = 1;
  }

  for (int ip = 0; ip < n; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0;
  }

  *nrot = 0;

  for (int i = 0; i < 50; i++) {
    double sm = 0;
    for (int ip = 0; ip < (n-1); ip++) {
      for (int iq = 0; iq < (n-1); iq++) {
        sm += fabs(a[ip][iq]);
      }
    }

    if (sm == 0) {
      return;
    }
    double tresh;
    if (i < 4) {
      tresh = 0.2 * sm / (SQR(n));
    } else {
      tresh = 0;
    }

    for (int ip = 0; ip < (n-1); ip++) {
      for (int iq = 0; iq < (n-1); iq++) {
        double g = 100 * fabs(a[ip][iq]);

        if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
            && (double)(fabs(d[iq])+g) == (double)fabs(d[iq])) {
          a[ip][iq] = 0;
        } else if (fabs(a[ip][iq]) > tresh) {
          double h = d[iq] - d[ip];
          double t, theta;

          if ((double)(fabs(h)+g) == (double)fabs(h)) {
            t = a[ip][iq] / h;
          } else {
            theta = 0.5 * h / a[ip][iq];
            t = 1 / (fabs(theta) + sqrt(1 + SQR(theta)));
            if (theta < 0) {
              t = -t;
            }
          }

          double c = 1 / sqrt(1 + SQR(t));
          double s = t * c;
          double tau = s / (1 + c);
          h = t * a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq] = 0;

          for (int j = 0; j < (iq-1); j++) {
            ROTATE(a, j, ip, j, iq);
          }
          for (int j = (ip+1); j < (iq-1); j++) {
            ROTATE(a, ip, j, j, iq);
          }
          for (int j = (iq+1); j < n; j++) {
            ROTATE(a, ip, j, iq, j);
          }
          for (int j = 0; j < n; j++) {
            ROTATE(v, j, ip, j, iq);
          }
          ++(*nrot);
        }
      }
    }

    for (int ip = 0; ip < n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0;
    }
  }
} //}}}

Vector Gyration(int n, int *list, Counts Counts, Vector BoxLength, BeadType *BeadType, Bead **Bead) { //{{{
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

  Vector com = CentreOfMass(n, list, *Bead, BeadType);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= com.x;
    (*Bead)[list[i]].Position.y -= com.y;
    (*Bead)[list[i]].Position.z -= com.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    GyrationTensor.x.x += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.x;
    GyrationTensor.x.y += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.y;
    GyrationTensor.x.z += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.z;
    GyrationTensor.y.y += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.y;
    GyrationTensor.y.z += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.z;
    GyrationTensor.z.z += (*Bead)[list[i]].Position.z * (*Bead)[list[i]].Position.z;
  }
  GyrationTensor.x.x /= n;
  GyrationTensor.x.y /= n;
  GyrationTensor.x.z /= n;
  GyrationTensor.y.y /= n;
  GyrationTensor.y.z /= n;
  GyrationTensor.z.z /= n;

  // just pro forma
  GyrationTensor.y.x = GyrationTensor.x.y;
  GyrationTensor.z.x = GyrationTensor.x.z;
  GyrationTensor.z.y = GyrationTensor.y.z;
  //}}}

  // create variables and arrays for jacobi() //{{{
  double **a;
  a = malloc(3*sizeof(double *));
  for (int i = 0; i < 3; i++) {
    a[i] = malloc(3*sizeof(double));
  }

  a[0][0] = GyrationTensor.x.x;
  a[0][1] = GyrationTensor.y.x;
  a[0][2] = GyrationTensor.z.x;
  a[1][0] = GyrationTensor.x.y;
  a[1][1] = GyrationTensor.y.y;
  a[1][2] = GyrationTensor.z.y;
  a[2][0] = GyrationTensor.x.z;
  a[2][1] = GyrationTensor.y.z;
  a[2][2] = GyrationTensor.z.z;

  double *d = malloc(3*sizeof(double));
  double **v = malloc(3*sizeof(double *));
  for (int i = 0; i < 3; i++) {
    v[i] = malloc(3*sizeof(double));
  }
  int nrot, size = 3; //}}}
  jacobi(a, size, d, v, &nrot);

  Vector eigen;

  eigen.x = d[0];
  eigen.y = d[1];
  eigen.z = d[2];

  eigen = Sort3(eigen);

  free(d);
  for (int i = 0; i < 3; i++) {
    free(v[i]);
    free(a[i]);
  }
  free(v);
  free(a);

  return (eigen);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("\
GyrationAggregates calculates radii of gyration, acylindricities, \
asphericities and relative shape anisotropies during the simulation for all \
aggregates izes. The shape descriptors are calculated from eigenvalues of \
gyration tensor. It also prints averages to the stdout. Instead of aggregate \
size, a number of specified molecular species in an aggregate can be used and \
only specified bead types can be used for all calculations.\n\n");

      printf("\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <input.agg> <output> <agg sizes> <options>\n\n", argv[0]);

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <input.agg>         input filename with information about aggregates (agg format)\n");
      printf("   <output>            output file with data during simulation run\n");
      printf("   <agg sizes>         aggregate sizes to calculate radius of gyration for\n");
      printf("   <options>\n");
      printf("      --joined         specify that aggregates with joined coordinates are used\n");
      printf("      -bt              specify bead types to be used for calculation (default is all)\n");
      printf("      -m <name>        agg size means number of <name> molecule types in an aggregate\n");
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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
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
    fprintf(stderr, "Too little mandatory arguments (%d instead of at least %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
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
  bool verbose2, silent;
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  VerboseLongOption(argc, argv, &verbose, &verbose2); // more verbose output
  SilentOption(argc, argv, &verbose, &verbose2, &silent); // no output
  bool script = BoolOption(argc, argv, "--script"); // do not use \r & co.
  // }}}

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

  // do not count unimers towards averages //{{{
  bool no_uni = BoolOption(argc, argv, "--no-unimers"); //}}}
  //}}}

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

  // <output> - filename with data during simulation run //{{{
  char output[32];
  strcpy(output, argv[++count]); //}}}

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

  // -m <name> option - specify MoleculeType that is used for determining agg sizes //{{{
  int specific_molecule = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-m") == 0) {

      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument for '-m' option (or molecule name beginning with a dash)!\n");
        ErrorHelp(argv[0]);
        exit(1);
      } //}}}

      specific_molecule = FindMoleculeType(argv[i+1], Counts, MoleculeType);
      if (specific_molecule == -1) {
        fprintf(stderr, "Molecule '%s' does not exist in FIELD ('-m' option)!\n", argv[i+1]);
        exit(1);
      }
    }
  } //}}}

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
  }

  int aggs = 0;

  for (int i = 0; i < Counts.Molecules; i++) {

    agg_sizes[aggs][0] = i + 1;

    aggs++; // number of aggregate sizes
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, Counts, &BeadType)) {
    exit(0);
  } //}}}

  // write initial stuff to output file //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output);
    exit(1);
  }

  // print command to output file
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out);

  // print legend line to output file
  fprintf(out, "# 1:dt 2:<Rg>_n 3:<Rg>_w 4:<Rg>_z ");
  fprintf(out, "5:<Anis>_n 6:<Acyl>_n 7:<Aspher>_n\n");

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
  // reading ends if next argument (beginning with '-') or the following empty line is read
  while ((test = getc(agg)) != '-' && test != '\n') {
    ungetc(test, agg);

    char name[10];
    fscanf(agg, "%s", name);
    int type = FindBeadType(name, Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", name, input_vcf);
      exit(1);
    }

    while ((test = getc(agg)) == ' ')
      ;
    ungetc(test, agg);
  } //}}}
  fclose(agg);

  // open again for production run - to ensure the pointer position in file is correct (at first 'Step')
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

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
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // allocate memory for sum of various things //{{{
  double sum_mass = 0, // total mass of aggregates throughout simulation
         sum_sqr_mass = 0; // sum of squared mass of aggregates throughout simulation
  double **Rg_sum = malloc(aggs*sizeof(double *));
  double *Anis_sum = calloc(aggs,sizeof(double));
  double *Acyl_sum = calloc(aggs,sizeof(double));
  double *Aspher_sum = calloc(aggs,sizeof(double));
  int **Size_sum = malloc(aggs*sizeof(int *));
  int **Molecules_sum = malloc(aggs*sizeof(int *));
  for (int i = 0; i < aggs; i++) {
    Rg_sum[i] = calloc(3,sizeof(double));
    Size_sum[i] = calloc(2,sizeof(int));
    Molecules_sum[i] = calloc((Counts.TypesOfMolecules*2),sizeof(int));
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
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

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate arrays for the timestep //{{{
    int *agg_counts = calloc(aggs,sizeof(int));
    double **Rg = malloc(aggs*sizeof(double *));
    double *Anis = calloc(aggs,sizeof(double));
    double *Acyl = calloc(aggs,sizeof(double));
    double *Aspher = calloc(aggs,sizeof(double));
    for (int i = 0; i < aggs; i++) {
      Rg[i] = calloc(3,sizeof(double));
    } //}}}

    // calculate shape descriptors //{{{
    double step_mass = 0, // total mass of aggregates
           step_sqr_mass = 0; // sum of squares of aggregate masses
    for (int i = 0; i < Counts.Aggregates; i++) {
      if (!no_uni || Aggregate[i].nMolecules != 1) {

        // test if aggregate 'i' should be used //{{{
        int mols = 0; // agg size
        if (specific_molecule != -1) { // agg size = number of molecules of type 'specific_molecule'
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int id = Aggregate[i].Molecule[j];
            if (specific_molecule == Molecule[id].Type) {
              mols++;
            }
          }
        } else { // agg size = total number of all molecules
          mols = Aggregate[i].nMolecules;
        }
        // is 'mols' agg size in provided list?
        int correct_size = -1;
        for (int j = 0; j < aggs; j++) {
          if (agg_sizes[j][0] == mols) {
            correct_size = j;
          }
        } //}}}

        if (correct_size != -1) {
          agg_counts[correct_size]++;
          agg_sizes[correct_size][1]++;

          // copy bead ids to a separate array //{{{
          int *list = malloc(Aggregate[i].nBeads*sizeof(int));
          int n = 0;
          for (int j = 0; j < Aggregate[i].nBeads; j++) {
            int id = Aggregate[i].Bead[j];
            if (BeadType[Bead[id].Type].Use) {
              list[n] = id;
              n++;
            }
          } //}}}

          Vector eigen = Gyration(n, list, Counts, BoxLength, BeadType, &Bead);

  //      // calcule Rg the 'usual way' -- for testing purposes //{{{
  //      double Rg2 = 0;
  //      Vector com = CentreOfMass(n, list, Bead, BeadType);
  //      for (int j = 0; j < n; j++) {
  //        Vector rij = Distance(Bead[list[j]].Position, com, BoxLength);
  //        Rg2 += SQR(rij.x) + SQR(rij.y) + SQR(rij.z);
  //      }
  //      Rg2 /= n; //}}}

          free(list); // free array of bead ids for gyration calculation

          double Rgi = sqrt(eigen.x + eigen.y + eigen.z); 

          // agg masses
          step_mass += Aggregate[i].Mass; // for this timestep
          step_sqr_mass += SQR(Aggregate[i].Mass); // for this timestep
          sum_mass += Aggregate[i].Mass; // for all timesteps
          sum_sqr_mass += SQR(Aggregate[i].Mass); // for all timesteps
          // Radius of gyration
          Rg[correct_size][0] += Rgi; // number avg
          Rg[correct_size][1] += Rgi * Aggregate[i].Mass; // weight average
          Rg[correct_size][2] += Rgi * SQR(Aggregate[i].Mass); // z-average
          // relative shape anisotropy
          Anis[correct_size] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) / SQR(eigen.x + eigen.y + eigen.z) - 0.5;
          // acylindricity
          Acyl[correct_size] += eigen.y - eigen.x;
          // asphericity
          Aspher[correct_size] += eigen.z - 0.5 * (eigen.x + eigen.y);
          // aggregate size
          Size_sum[correct_size][0] +=     Aggregate[i].nMolecules;
          Size_sum[correct_size][1] += SQR(Aggregate[i].nMolecules);

          // count number of Diblocks and Surfacts in aggrate
          // TODO: generalise
          int Diblock = 0, Surfact = 0, Fluores = 0;
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name, "Diblock") == 0) {
              Diblock++;
            } else if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name, "Surfact") == 0) {
              Surfact++;
            } else {
              Fluores++;
            }
          }

          Molecules_sum[correct_size][0] +=     Diblock;
          Molecules_sum[correct_size][1] += SQR(Diblock);
          Molecules_sum[correct_size][2] +=     Surfact;
          Molecules_sum[correct_size][3] += SQR(Surfact);
          Molecules_sum[correct_size][4] +=     Fluores;
          Molecules_sum[correct_size][5] += SQR(Fluores);
        }
      }
    } //}}}

    // add values to sums //{{{
    for (int i = 0; i < aggs; i++) {
      Rg_sum[i][0] += Rg[i][0];
      Rg_sum[i][1] += Rg[i][1];
      Rg_sum[i][2] += Rg[i][2];
      Anis_sum[i] += Anis[i];
      Acyl_sum[i] += Acyl[i];
      Aspher_sum[i] += Aspher[i];
    } //}}}

    // print data to output file //{{{
    FILE *out;
    if ((out = fopen(output, "a")) == NULL) { // out file opened fine? //{{{
      // print newline to stdout if Step... doesn't end with one
      if (!script && !silent) {
        putchar('\n');
      }
      fprintf(stderr, "Cannot open file %s!\n", output);
      exit(1);
    } //}}}

    fprintf(out, "%5d", count); // timestep
    // radius of gyration
    for (int i = 1; i < aggs; i++) {
      Rg[0][0] += Rg[i][0];
      Rg[0][1] += Rg[i][1];
      Rg[0][2] += Rg[i][2];
      Anis[0] += Anis[i];
      Acyl[0] += Acyl[i];
      Aspher[0] += Aspher[i];

      agg_counts[0] += agg_counts[i];
    }
    fprintf(out, " %8.5f %8.5f %8.5f", Rg[0][0]/agg_counts[0], Rg[0][1]/step_mass, Rg[0][2]/step_sqr_mass);
    fprintf(out, " %8.5f", Anis[0]/agg_counts[0]);
    fprintf(out, " %8.5f", Acyl[0]/agg_counts[0]);
    fprintf(out, " %8.5f", Aspher[0]/agg_counts[0]);
    putc('\n', out);

    fclose(out); //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}

    // free memory //{{{
    free(agg_counts);
    for (int i = 0; i < aggs; i++) {
      free(Rg[i]);
    }
    free(Rg);
    free(Anis);
    free(Acyl);
    free(Aspher); //}}}
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

  // calculate simple averages //{{{
  printf("1:Size 2:<Rg> 3:<Acyl> 4:<Ashper> 5:<Anis> 6:<Diblocks> 7:<Surfact>\n");
  for (int i = 0; i < aggs; i++) {
    if (agg_sizes[i][1] > 0) {
      printf("%8d", agg_sizes[i][0]);
      printf(" %7.3f", Rg_sum[i][0]/agg_sizes[i][1]);
      printf(" %7.3f", Anis_sum[i]/agg_sizes[i][1]);
      printf(" %7.3f", Acyl_sum[i]/agg_sizes[i][1]);
      printf(" %7.3f", Aspher_sum[i]/agg_sizes[i][1]);
      printf(" %7.3f", (double)(Molecules_sum[i][0])/agg_sizes[i][1]);
      printf(" %7.3f", (double)(Molecules_sum[i][2])/agg_sizes[i][1]);
      putchar('\n');
    }
  } //}}}

  // total averages //{{{
  for (int i = 1; i < aggs; i++) {
    Rg_sum[0][0] += Rg_sum[i][0];
    Rg_sum[0][1] += Rg_sum[i][1];
    Rg_sum[0][2] += Rg_sum[i][2];
    Anis_sum[0] += Anis_sum[i];
    Acyl_sum[0] += Acyl_sum[i];
    Aspher_sum[0] += Aspher_sum[i];

    agg_sizes[0][1] += agg_sizes[i][1];

    Size_sum[0][0] += Size_sum[i][0];
    Size_sum[0][1] += Size_sum[i][1];

    Molecules_sum[0][0] += Molecules_sum[i][0];
    Molecules_sum[0][1] += Molecules_sum[i][1];
    Molecules_sum[0][2] += Molecules_sum[i][2];
    Molecules_sum[0][3] += Molecules_sum[i][3];
    Molecules_sum[0][4] += Molecules_sum[i][4];
    Molecules_sum[0][5] += Molecules_sum[i][5];
  }

  // print to stdout
  printf("1:<A_s>_w 2:<A_s>_n 3:<R_G>_n 4:_w 5:_z 6:<Anis> 7:<Acyl> 8:<Aspher>\n");
  printf("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
    (double)(Size_sum[0][1])/Size_sum[0][0], //<A_s>_w
    (double)(Size_sum[0][0])/agg_sizes[0][1], //<A_s>_n
    Rg_sum[0][0]/agg_sizes[0][1],
    Rg_sum[0][1]/sum_mass,
    Rg_sum[0][2]/sum_sqr_mass,
    Anis_sum[0]/agg_sizes[0][1],
    Acyl_sum[0]/agg_sizes[0][1],
    Aspher_sum[0]/agg_sizes[0][1]);

  if ((out = fopen(output, "a")) == NULL) { //{{{
    fprintf(stderr, "Cannot open file %s!\n", output);
    exit(1);
  } //}}}

  fprintf(out, "# 1:<A_s>_w 2:<A_s>_n 3:<R_G>_n 4:_w 5:_z 6:<Anis> 7:<Acyl> 8:<Aspher>\n");
  fprintf(out, "# %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
    (double)(Size_sum[0][1])/Size_sum[0][0], //<A_s>_w
    (double)(Size_sum[0][0])/agg_sizes[0][1], //<A_s>_n
    Rg_sum[0][0]/agg_sizes[0][1],
    Rg_sum[0][1]/sum_mass,
    Rg_sum[0][2]/sum_sqr_mass,
    Anis_sum[0]/agg_sizes[0][1],
    Acyl_sum[0]/agg_sizes[0][1],
    Aspher_sum[0]/agg_sizes[0][1]);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < aggs; i++) {
    free(agg_sizes[i]);
    free(Rg_sum[i]);
    free(Size_sum[i]);
    free(Molecules_sum[i]);
  }
  free(Molecules_sum);
  free(Size_sum);
  free(agg_sizes);
  free(Rg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(stuff); //}}}

  return 0;
}
