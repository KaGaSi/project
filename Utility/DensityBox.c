#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
DensityBox utility calculates number \
density for all bead types in the direction of specified axis (x, y, or z).\
\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output.rho> <axis> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename (either vcf or vtf format)\n");
  fprintf(ptr, "   <width>           width of a single bin\n");
  fprintf(ptr, "   <output.rho>      output density file (automatic ending '<axis>.rho' added)\n");
  fprintf(ptr, "   <axis>            calculate along x, y, or z axis\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -n <int>       number of bins to average\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       number of timestep to end with\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-x") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  char *input_vsf = calloc(1024,sizeof(char));
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // ending timestep //{{{
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  } //}}}

  // number of bins to average //{{{
  int avg = 1;
  if (IntegerOption(argc, argv, "-n", &avg)) {
    exit(1);
  } //}}}

  // error if ending step is lower than starging step //{{{
  if (end != -1 && start > end) {
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n", start, end);
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[1024];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <width> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.rho> - filename with bead densities //{{{
  char output_rho[1024];
  strcpy(output_rho, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // <axis> - x, y, or z //{{{
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {

    axis = argv[count][0];

    if (axis != 'x' && axis != 'y' && axis != 'z') {
      fprintf(stderr, "\nError: <axis> must be 'x', 'y', or 'z'\n\n");
      exit(1);
    }
  } //}}}

  // write initial stuff to output density file //{{{
  FILE *out;
  char str[1050];

  sprintf(str, "%s%c.rho", output_rho, axis);
  strcpy(output_rho, str);
  if ((out = fopen(output_rho, "w")) == NULL) {
    ErrorFileOpen(output_rho, 'w');
    exit(1);
  }

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print bead type names to output file //{{{
  fprintf(out, "# columns: (1) distance;");
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(out, " (%d) %s", i+2, BeadType[i].Name);
    if (i != (Counts.TypesOfBeads-1)) {
      putc(';', out);
    }
  }
  putc('\n', out);
//for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  fprintf(out, " %d: %s", 4*i+2, BeadType[i].Name);
//}
//fprintf(out, "\n# for each molecule type: rdp | stderr | rnp | stderr\n"); //}}}

  fclose(out); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_coor);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // number of bins //{{{
  double max_dist;
  if (axis == 'x') {
    max_dist = BoxLength.x;
  } else if (axis == 'y') {
    max_dist = BoxLength.y;
  } else {
    max_dist = BoxLength.z;
  }
  int bins = ceil(max_dist / width); //}}}

  // allocate memory for density arrays //{{{
  long int **rho = malloc(Counts.TypesOfBeads*sizeof(long int *));
  long int **rho_2 = malloc(Counts.TypesOfBeads*sizeof(long int *));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    rho[i] = calloc(bins,sizeof(long int));
    rho_2[i] = calloc(bins,sizeof(long int));
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(1024*sizeof(int));

  // initialize the array
  for (int i = 0; i < 1024; i++) {
    stuff[i] = '\0';
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(stdout, " %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
  } //}}}

  // skip first start-1 steps //{{{
  count = 0;
  int test;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %6d", count);
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "Error: premature end of %s file\n\n", input_coor);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step: %6d\n", start);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rStarting step: %6d   \n", start);
    }
  } //}}}
  //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int count_vcf = start - 1;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // write step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %6d", count_vcf);
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count_vcf, stuff, input_vsf);
      exit(1);
    } //}}}

    // add pbc //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      while (Bead[i].Position.x >= BoxLength.x) {
        Bead[i].Position.x -= BoxLength.x;
      }
      while (Bead[i].Position.x < 0) {
        Bead[i].Position.x += BoxLength.x;
      }
      while (Bead[i].Position.y >= BoxLength.y) {
        Bead[i].Position.y -= BoxLength.y;
      }
      while (Bead[i].Position.y < 0) {
        Bead[i].Position.y += BoxLength.y;
      }
      while (Bead[i].Position.z >= BoxLength.z) {
        Bead[i].Position.z -= BoxLength.z;
      }
      while (Bead[i].Position.z < 0) {
        Bead[i].Position.z += BoxLength.z;
      }
    } //}}}

    // allocate memory for temporary density arrays //{{{
    int **temp_rho = malloc(Counts.TypesOfBeads*sizeof(int *));
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      temp_rho[i] = calloc(bins,sizeof(int));
    } //}}}

    // calculate densities //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      bool use = true;
      int mol = Bead[i].Molecule;
      if (mol != -1) {
        int mtype = Molecule[mol].Type;
        use = MoleculeType[mtype].Use;
      }
      if (use) {
        if (axis == 'x') {
          int j = Bead[i].Position.x / width;
          temp_rho[Bead[i].Type][j]++;
        } else if (axis == 'y') {
          int j = Bead[i].Position.y / width;
          temp_rho[Bead[i].Type][j]++;
        } else {
          int j = Bead[i].Position.z / width;
          temp_rho[Bead[i].Type][j]++;
        }
      }
    } //}}}

    // add from temporary density array to global density arrays //{{{
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < bins; k++) {
        rho[j][k] += temp_rho[j][k];
        rho_2[j][k] += SQR(temp_rho[j][k]);
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      free(temp_rho[i]);
    }
    free(temp_rho); //}}}

    if (end == count_vcf)
      break;
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count_vcf);
    }
  } //}}}

  // write densities to output file(s) //{{{
  if ((out = fopen(output_rho, "a")) == NULL) {
    ErrorFileOpen(output_rho, 'a');
    exit(1);
  }

  // calculate rdf
  double volume = width * avg;
  if (axis == 'x') {
    volume *= BoxLength.y * BoxLength.z;
  } else if (axis == 'y') {
    volume *= BoxLength.x * BoxLength.z;
  } else {
    volume *= BoxLength.x * BoxLength.y;
  }
  for (int i = 0; i < (bins-avg); i++) {

    fprintf(out, "%7.3f", width*(i+0.5*avg));

    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      double temp_rho = 0, temp_number = 0,
             temp_rho_err = 0, temp_number_err = 0;

      // sum densities to be averaged
      for (int k = 0; k < avg; k++) {
        temp_rho += rho[j][i+k] / (volume * count);
//      temp_rho_err += rho_2[j][i+k] / (volume * BeadType[j].Number * count);
//      temp_number += rho[j][i+k] / (BeadType[j].Number * count);
//      temp_number_err += rho_2[j][i+k] / (BeadType[j].Number * count);
      }

      temp_rho_err = sqrt(temp_rho_err - temp_rho);
      temp_number_err = sqrt(temp_number_err - temp_number);

      // print average value to output file
      fprintf(out, " %10f", temp_rho);
//    fprintf(out, " %10f %10f", temp_rho, temp_rho_err);
//    fprintf(out, " %10f %10f", temp_number, temp_number_err);
    }
    putc('\n',out);
  }

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
