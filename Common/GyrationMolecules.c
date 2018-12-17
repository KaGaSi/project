#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <output> <molecule(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <output>          output file with shape descriptors (automatic ending '-<name>.txt')\n");
  fprintf(stderr, "   <molecule(s)>     molecule types to calculate shape descriptors for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined       specify that joined coordinates are used\n");
  fprintf(stderr, "      -bt            specify bead types to be used for calculation (default is all)\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
GyrationAggregates calculates radii of gyration during the simulation for \
given molecule(s). The radius of gyration is calculated from eigenvalues \
of gyration tensor. It also prints average radii of gyration to the screen. \
Bead types to be used for calculation can be specified.\n\n");

/*      fprintf(stdout, "\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");
*/

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input.vcf> <output> <molecule(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input.vcf>       input filename (vcf format)\n");
      fprintf(stdout, "   <output>          output file with shape descriptors (automatic ending '-<name>.txt')\n");
      fprintf(stdout, "   <molecule(s)>     molecule types to calculate shape descriptors for\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      --joined       specify that joined coordinates are used\n");
      fprintf(stdout, "      -bt            specify bead types to be used for calculation (default is all)\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
//      strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "--joined") != 0) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than dl_meso.vsf? //{{{
  char *input_vsf = calloc(32,sizeof(char *));
  if (VsfFileOption(argc, argv, &input_vsf)) {
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
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf[32];
  strcpy(input_vcf, argv[++count]);

  // test if <input.vcf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(input_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    ErrorExtension(input_vcf, ".vcf");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <output> - filename with shape descriptors //{{{
  char output[32];
  strcpy(output, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

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
      fprintf(stderr, "\nError: non-existent molecule name: %s\n\n", argv[count]);
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
      char str[128];
      sprintf(str, "%s-%s.txt", output, MoleculeType[i].Name);

      FILE *out;
      if ((out = fopen(str, "w")) == NULL) {
        ErrorFileOpen(str, 'w');
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
    ErrorFileOpen(input_vcf, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[128];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_vcf);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "\nError: cannot read pbc from %s\n\n", input_vcf);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
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
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
      }
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
      exit(1);
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
        sprintf(str, "%s-%s.txt", output, MoleculeType[i].Name);

        FILE *out;
        if ((out = fopen(str, "a")) == NULL) {
          ErrorFileOpen(str, 'a');
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
      fprintf(stdout, "\n%s", stuff);
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
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  } //}}}

//// calculate simple averages //{{{
//for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//  fprintf(stdout, "%10s %lf\n", MoleculeType[i].Name, Rg_sum[i][0]/(count*MoleculeType[i].Number));
//} //}}}

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
