#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <output> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <output>            output file with shape descriptors (automatic ending '-<name>.text')\n");
  fprintf(stderr, "   <molecule names>    molecule types to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined         specify that joined coordinates are used\n");
  fprintf(stderr, "      -bt              specify bead types to be used for calculation (default is all)\n");
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

  Vector com = CentreOfMass(n, list, *Bead, BeadType);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= com.x;
    (*Bead)[list[i]].Position.y -= com.y;
    (*Bead)[list[i]].Position.z -= com.z;
  } //}}}

  double **a = malloc(3*sizeof(double *)); //{{{
  a[0] = calloc(3,sizeof(double));
  a[1] = calloc(3,sizeof(double));
  a[2] = calloc(3,sizeof(double));

//a[0][0] = GyrationTensor.x.x;
//a[0][1] = GyrationTensor.y.x;
//a[0][2] = GyrationTensor.z.x;
//a[1][0] = GyrationTensor.x.y;
//a[1][1] = GyrationTensor.y.y;
//a[1][2] = GyrationTensor.z.y;
//a[2][0] = GyrationTensor.x.z;
//a[2][1] = GyrationTensor.y.z;
//a[2][2] = GyrationTensor.z.z; //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    a[0][0] += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.x;
    a[0][1] += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.y;
    a[0][2] += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.z;
    a[1][1] += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.y;
    a[1][2] += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.z;
    a[2][2] += (*Bead)[list[i]].Position.z * (*Bead)[list[i]].Position.z;
  }
  a[0][0] /= n;
  a[0][1] /= n;
  a[0][2] /= n;
  a[1][1] /= n;
  a[2][1] /= n;
  a[2][2] /= n;

  // just pro forma
  a[1][0] = a[0][1];
  a[2][0] = a[0][2];
  a[2][1] = a[1][2];
  //}}}

  putchar('\n');
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      printf(" %lf", a[i][j]);
    }
    putchar('\n');
  }
  putchar('\n');

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

  eigen = Sort3(eigen);

  printf("d=(%lf, %lf, %lf)\n", d[0], d[1], d[2]);
  putchar('\n');

  // free memory //{{{
  free(d);
  free(v[0]);
  free(v[1]);
  free(v[2]);
  free(v);
  free(a[0]);
  free(a[1]);
  free(a[2]);
  free(a); //}}}

  return (eigen);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("\
GyrationAggregates calculates radii of gyration during the simulation for \
given molecule type(s). The radius of gyration is calculated from eigenvalues \
of gyration tensor. It also prints average radii of gyration to the screen. \
Currently, it uses all beads present in the molecules.\n\n");

      printf("\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <output> <molecule names> <options>\n\n", argv[0]);

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <output>            output file with shape descriptors (automatic ending '-<name>.text')\n");
      printf("   <molecule names>    molecule types to calculate radius of gyration for\n");
      printf("   <options>\n");
      printf("      --joined         specify that joined coordinates are used\n");
      printf("      -bt              specify bead types to be used for calculation (default is all)\n");
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
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "--joined") != 0) {

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

  // <output> - filename with shape descriptors //{{{
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

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, Counts, &BeadType)) {
    exit(0);
  } //}}}

  // write initial stuff to output file //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      char str[32];
      sprintf(str, "%s-%s.txt", output, MoleculeType[i].Name);

      FILE *out;
      if ((out = fopen(str, "w")) == NULL) {
        fprintf(stderr, "Cannot open file %s!\n", str);
        exit(1);
      }

      // print command to output file //{{{
      putc('#', out);
      for (int i = 0; i < argc; i++)
        fprintf(out, " %s", argv[i]);
      putc('\n', out); //}}}

      fprintf(out, "# %s\n", MoleculeType[i].Name);
      fprintf(out, "# 1:dt 2:<Rg>_n 3:<Rg>_w 4:<Rg>_z ");
      fprintf(out, "5:<Anis>_n 6:<Acyl>_n 7:<Aspher>_n\n");
      putc('\n', out);

      fclose(out);
    }
  } //}}}

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

  // allocate memory for sums of shape descriptors //{{{
  double **Rg_sum = malloc(Counts.TypesOfMolecules*sizeof(double *));
  double *Anis_sum = calloc(Counts.TypesOfMolecules,sizeof(double));
  double *Acyl_sum = calloc(Counts.TypesOfMolecules,sizeof(double));
  double *Aspher_sum = calloc(Counts.TypesOfMolecules,sizeof(double));
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    Rg_sum[i] = calloc(3,sizeof(double));
  } //}}}

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

    // allocate memory for shape descriptors //{{{
    double **Rg = malloc(Counts.TypesOfMolecules*sizeof(double *));
    double *Anis = calloc(Counts.TypesOfMolecules,sizeof(double));
    double *Acyl = calloc(Counts.TypesOfMolecules,sizeof(double));
    double *Aspher = calloc(Counts.TypesOfMolecules,sizeof(double));
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      Rg[i] = calloc(3,sizeof(double));
    } //}}}

    // calculate shape descriptors //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type = Molecule[i].Type;

      if (MoleculeType[mol_type].Use) {

        // copy bead ids to a separate array //{{{
        int *list = malloc(MoleculeType[mol_type].nBeads*sizeof(int));
        int n = 0;
        for (int j = 0; j < MoleculeType[mol_type].nBeads; j++) {
          int bead_id = Molecule[i].Bead[j];
          if (BeadType[Bead[bead_id].Type].Use) {
            list[n] = Molecule[i].Bead[j];
            n++;
          }
        } //}}}

        Vector eigen = Gyration(n, list, Counts, BoxLength, BeadType, &Bead);

        free(list); // free array of bead ids for gyration calculation

        // Radius of gyration
        Rg[Molecule[i].Type][0] +=      sqrt(eigen.x + eigen.y + eigen.z);
        Rg[Molecule[i].Type][1] +=           eigen.x + eigen.y + eigen.z;
        Rg[Molecule[i].Type][2] += CUBE(sqrt(eigen.x + eigen.y + eigen.z));
        // relative shape anisotropy
        Anis[Molecule[i].Type] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) / SQR(eigen.x + eigen.y + eigen.z) - 0.5;
        // acylindricity
        Acyl[Molecule[i].Type] += eigen.y - eigen.x;
        // asphericity
        Aspher[Molecule[i].Type] += eigen.z - 0.5 * (eigen.x + eigen.y);
        printf("%lf %lf %lf\n", eigen.x, eigen.y, eigen.z);
      }
    } //}}}

    // add shape descriptors to sum //{{{
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      Rg_sum[i][0] += Rg[i][0];
      Rg_sum[i][1] += Rg[i][1];
      Rg_sum[i][2] += Rg[i][2];
      Anis_sum[i] += Anis[i];
      Acyl_sum[i] += Acyl[i];
      Aspher_sum[i] += Aspher[i];
    } //}}}

    // print shape descriptors to output file(s) //{{{
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        char str[32];
        sprintf(str, "%s-%s.txt", output, MoleculeType[i].Name);

        FILE *out;
        if ((out = fopen(str, "a")) == NULL) {
          fprintf(stderr, "Cannot open file %s!\n", str);
          exit(1);
        }

        fprintf(out, "%5d", count);
        fprintf(out, " %8.5f %8.5f %8.5f", Rg[i][0]/MoleculeType[i].Number, Rg[i][1]/Rg[i][0], Rg[i][2]/Rg[i][1]);
        fprintf(out, " %8.5f", Anis[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Acyl[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Aspher[i]/MoleculeType[i].Number);
        putc('\n', out);

        fclose(out);
      }
    } //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}

    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      free(Rg[i]);
    }
    free(Rg);
    free(Anis);
    free(Acyl);
    free(Aspher);
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
    printf("%10s %lf\n", MoleculeType[i].Name, Rg_sum[i][0]/(count*MoleculeType[i].Number));
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    free(Rg_sum[i]);
  }
  free(Rg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(stuff); //}}}

  return 0;
}
