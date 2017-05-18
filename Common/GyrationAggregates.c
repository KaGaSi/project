#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input.agg> <output> <agg sizes> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <input.agg>         input filename with information about aggregates (agg format)\n");
  fprintf(stderr, "   <output>            output file with data during simulation run\n");
  fprintf(stderr, "   <agg sizes>         aggregate sizes to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -j               specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -t               specify bead types to be used for calculation (default is all)\n");
  fprintf(stderr, "      -m <name>        agg size means number of <name> molecule types in an aggregate\n");
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

Vector Gyration(int n, int *list, Counts Counts, Vector BoxLength, BeadType *BeadType, Bead **Bead) { //{{{ // gyration tensor (3x3 array) //{{{
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

  Vector com = CenterOfMass(n, list, *Bead, BeadType);

  // move center of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= com.x;
    (*Bead)[list[i]].Position.y -= com.y;
    (*Bead)[list[i]].Position.z -= com.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    GyrationTensor.x.x += SQR((*Bead)[list[i]].Position.x);
    GyrationTensor.x.y += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.y;
    GyrationTensor.x.z += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.z;
    GyrationTensor.y.y += SQR((*Bead)[list[i]].Position.y);
    GyrationTensor.y.z += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.z;
    GyrationTensor.z.z += SQR((*Bead)[list[i]].Position.z);
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
  double **a = malloc(3*sizeof(double *));
  a[0] = malloc(3*sizeof(double));
  a[1] = malloc(3*sizeof(double));
  a[2] = malloc(3*sizeof(double));

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
  v[0] = malloc(3*sizeof(double));
  v[1] = malloc(3*sizeof(double));
  v[2] = malloc(3*sizeof(double));
  int nrot, size = 3; //}}}
  jacobi(a, size, d, v, &nrot);

  Vector eigen;

  eigen.x = d[0];
  eigen.y = d[1];
  eigen.z = d[2];

  eigen = Sort3(eigen);

  free(d);
  free(v[0]);
  free(v[1]);
  free(v[2]);
  free(v);

  return (eigen);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("GyrationAggregates calculates radii of gyration, acylindricities,           \n");
      printf("asphericities and relative shape anisotropies during the simulation for     \n");
      printf("aggregates of given size(s). The shape descriptors are calculated from      \n");
      printf("eigenvalues of gyration tensor. It also prints simplie averages to the      \n");
      printf("screen. Instead of aggregate size, a number of specified molecular species  \n");
      printf("in an aggregate can be used and only specified bead types can be used for   \n");
      printf("all calculations.                                                           \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <input.agg> <output> <agg sizes> <options>\n\n", argv[0]);

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <input.agg>         input filename with information about aggregates (agg format)\n");
      printf("   <output>            output file with data during simulation run\n");
      printf("   <agg sizes>         aggregate sizes to calculate radius of gyration for\n");
      printf("   <options>\n");
      printf("      -j               specify that aggregates with joined coordinates are used\n");
      printf("      -t               specify bead types to be used for calculation (default is all)\n");
      printf("      -m <name>        agg size means number of <name> molecule types in an aggregate\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 4; //}}}

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
        strcmp(argv[i], "-j") != 0 &&
        strcmp(argv[i], "-t") != 0 &&
        strcmp(argv[i], "-m") != 0) {

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

  // -j option - coordinates are joined //{{{
  bool joined = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      joined = true;
    }
  } //}}}

  // -t option - bead types to be used //{{{
  int types = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-t") == 0) {
      types = i; // positon of the '-t' argument in command
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

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (argv[count][0] < '1' || argv[count][0] > '9') {
      fprintf(stderr, "Non-numeric option in <agg sizes>!\n");
      exit(1);
    } //}}}

    agg_sizes[aggs][0] = atoi(argv[count]);

    aggs++; // number of aggregate sizes
  } //}}}

  // specify what bead types to use - either specified by '-t' option or all //{{{
  if (types != -1) { // '-t' option is present
    // <type names> - names of bead types to save
    while (++types < argc && argv[types][0] != '-') {
      int type = FindBeadType(argv[types], Counts, BeadType);

      BeadType[type].Use = true;
      printf("%s\n", BeadType[type].Name);
    }
  } else {
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      BeadType[i].Use = true;
    }
  } //}}}

  // write initial stuff to output file //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output);
    exit(1);
  }

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print agg sizes to output file //{{{
  fprintf(out, "# dt");
  for (int i = 0; i < aggs; i++) {
    fprintf(out, " |%4d: Rg  Anis     Acyl     Aspher   eigen_x  eigen_y  eigen_z", agg_sizes[i][0]);
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

    printf("Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      printf(" %d", agg_sizes[i][0]);
    }
    putchar('\n');
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // allocate memory for sum of radii of gyration //{{{
  double *Rg_sum = calloc(aggs,sizeof(double));
  double *Anis_sum = calloc(aggs,sizeof(double));
  double *Acyl_sum = calloc(aggs,sizeof(double));
  double *Aspher_sum = calloc(aggs,sizeof(double));
  double *eigen_x_sum = calloc(aggs,sizeof(double));
  double *eigen_y_sum = calloc(aggs,sizeof(double));
  double *eigen_z_sum = calloc(aggs,sizeof(double)); //}}}

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
    double *Rg = calloc(aggs,sizeof(double));
    double *Anis = calloc(aggs,sizeof(double));
    double *Acyl = calloc(aggs,sizeof(double));
    double *Aspher = calloc(aggs,sizeof(double));
    double *eigen_x = calloc(aggs,sizeof(double));
    double *eigen_y = calloc(aggs,sizeof(double));
    double *eigen_z = calloc(aggs,sizeof(double)); //}}}

    // calculate radii of gyration //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {

//    // OLD - test if aggregate is of correct size //{{{
//    int correct_size = -1;
//    for (int j = 0; j < aggs; j++) {
//      if (agg_sizes[j][0] == Aggregate[i].nMolecules) {
//        correct_size = j;
//      }
//    } //}}}
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

        int *list = malloc(Aggregate[i].nBeads*sizeof(int));
        int n = 0;

        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int id = Aggregate[i].Bead[j];
          if (BeadType[Bead[id].Type].Use) {
            list[n] = id;
            n++;
          }
        }

        Vector eigen = Gyration(n, list, Counts, BoxLength, BeadType, &Bead);

        eigen_x[correct_size] += eigen.x;
        eigen_y[correct_size] += eigen.y;
        eigen_z[correct_size] += eigen.z;

        Rg[correct_size] += sqrt(eigen.x + eigen.y + eigen.z);
        Anis[correct_size] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) / SQR(eigen.x + eigen.y + eigen.z) - 0.5;
        Acyl[correct_size] += eigen.y - eigen.x;
        Aspher[correct_size] += eigen.z - 0.5 * (eigen.x + eigen.y);
      }
    } //}}}

    // add values to sums //{{{
    for (int i = 0; i < aggs; i++) {
      eigen_x_sum[i] += eigen_x[i];
      eigen_y_sum[i] += eigen_y[i];
      eigen_z_sum[i] += eigen_z[i];

      Rg_sum[i] += Rg[i];
      Anis_sum[i] += Anis[i];
      Acyl_sum[i] += Acyl[i];
      Aspher_sum[i] += Aspher[i];
    } //}}}

    // print radii of gyration to output file //{{{
    FILE *out;
    if ((out = fopen(output, "a")) == NULL) {
      // print newline to stdout if Step... doesn't end with one
      if (!script && !silent) {
        putchar('\n');
      }
      fprintf(stderr, "Cannot open file %s!\n", output);
      exit(1);
    }

    fprintf(out, "%5d", count);
    for (int i = 0; i < aggs; i++) {
      fprintf(out, " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f", Rg[i]/agg_counts[i], Anis[i]/agg_counts[i], Acyl[i]/agg_counts[i], Aspher[i]/agg_counts[i], eigen_x[i]/agg_counts[i], eigen_y[i]/agg_counts[i], eigen_z[i]/agg_counts[i]);
    }
    putc('\n', out);

    fclose(out); //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}

    free(agg_counts);
    free(Rg);
    free(Anis);
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
  fprintf(out, "Size     Rg     Acyl     Ashper   Anis     tensor_x  tensor_y  tensor_z");
  for (int i = 0; i < aggs; i++) {
    printf("%8d %6f %6f %6f %6f %6f %6f %6f\n", agg_sizes[i][0], Rg_sum[i]/agg_sizes[i][1], Acyl_sum[i]/agg_sizes[i][1], Aspher_sum[i]/agg_sizes[i][1], Anis_sum[i]/agg_sizes[i][1], eigen_x_sum[i]/agg_sizes[i][1], eigen_y_sum[i]/agg_sizes[i][1], eigen_z_sum[i]/agg_sizes[i][1]);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(agg_sizes);
  free(Rg_sum);
  free(stuff); //}}}

  return 0;
}
