#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <output> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <output>            output file with radii of gyration\n");
  fprintf(stderr, "   <molecule names>    molecule types to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -j               specify that joined coordinates are used\n");
  CommonHelp(1);
} //}}}

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

  double **a = malloc(3*sizeof(double *)); //{{{
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
  a[2][2] = GyrationTensor.z.z; //}}}

  double *d = malloc(3*sizeof(double));
  double **v = malloc(3*sizeof(double *));
  v[0] = malloc(3*sizeof(double));
  v[1] = malloc(3*sizeof(double));
  v[2] = malloc(3*sizeof(double));
  int nrot, size = 3;
  jacobi(a, size, d, v, &nrot);

  Vector eigen;
  eigen.x = d[0];
  eigen.y = d[1];
  eigen.z = d[2];

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
      printf("GyrationAggregates calculates radii of gyration during the simulation for   \n");
      printf("given molecule type(s). The radius of gyration is calculated from           \n");
      printf("eigenvalues of gyration tensor. It also prints average radii of gyration to \n");
      printf("the screen. Currently, it uses all beads present in the molecules.          \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <output> <molecule names> <options>\n\n", argv[0]);

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <output>            output file with radii of gyration\n");
      printf("   <molecule names>    molecule types to calculate radius of gyration for\n");
      printf("   <options>\n");
      printf("      -j               specify that joined coordinates are used\n");
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
        strcmp(argv[i], "-j") != 0) {

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

  // <output> - filename with radii of gyration //{{{
  char output[16];
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

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {

    bool test = false;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (strcmp(argv[count], MoleculeType[i].Name) == 0) {
        MoleculeType[i].Use = true;

        test = true;

        break;
      }
    }

    // wrong molecule name //{{{
    if (!test) {
      fprintf(stderr, "Non-existent molecule name: %s!\n", argv[count]);
      exit(1);
    } //}}}
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

  // print molecule names to output file //{{{
  fprintf(out, "# timestep");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, " %10s", MoleculeType[i].Name);
    }
  }
  putc('\n', out); //}}}

  fclose(out); //}}}

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

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // allocate memory for sum of radii of gyration
  double *Rg_sum = calloc(Counts.TypesOfMolecules,sizeof(double));

  // main loop //{{{
  count = 0; // count timesteps
  int test;
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

    // join molecules if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate arrays for the timestep
    double *Rg = calloc(Counts.TypesOfMolecules,sizeof(double));

    // calculate radii of gyration //{{{
    for (int i = 0; i < Counts.Molecules; i++) {

      if (MoleculeType[Molecule[i].Type].Use) {
        Vector eigen = Gyration(MoleculeType[Molecule[i].Type].nBeads, Molecule[i].Bead, Counts, BoxLength, BeadType, &Bead);

        Rg[Molecule[i].Type] += sqrt(eigen.x + eigen.y + eigen.z);
      }
    } //}}}

    // add radii to sum //{{{
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      Rg_sum[i] += Rg[i];
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

    fprintf(out, "%d", count);
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %lf", Rg[i]/MoleculeType[i].Number);
      }
    }
    putc('\n', out);

    fclose(out); //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}

    free(Rg);
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
    }
  } //}}}

  // calculate simple averages //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf("%10s %lf\n", MoleculeType[i].Name, Rg_sum[i]/(count*MoleculeType[i].Number));
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(Rg_sum);
  free(stuff); //}}}

  return 0;
}
