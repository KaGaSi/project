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
Surface utility determines the first bead (either from the box's centre or \
edges) in each square prism defined by the given <width> parameter, thus \
defining a surface of, e.g., polymer brush or lipid bilayer. The <width> \
slices the box into square prisms along the chosen axis (i.e., if z is the \
chosen axis, the xy plane is chopped into squares, creating \
<width>*<width>*<box length in z> prisms). In each such prism, two beads \
are found corresponding to the two surfaces (e.g., polymer brush on both box \
edges or the two surfaces of a lipid bilayer inside the box).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <axis> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename (either vcf or vtf format)\n");
  fprintf(ptr, "   <width>           width of a single bin\n");
  fprintf(ptr, "   <output>          output file\n");
  fprintf(ptr, "   <axis>            calculate along x, y, or z axis\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -in            start from the box's centre instead of edges\n");
  fprintf(ptr, "      -m <name(s)>   molecule type(s) to use\n");
  fprintf(ptr, "      -bt <name(s)>  bead type(s) to use\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <int>       number of timestep to end with\n");
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
        strcmp(argv[i], "-in") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  char *input_vsf = calloc(LINE,sizeof(char));
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

  // error if ending step is lower than starging step //{{{
  if (end != -1 && start > end) {
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n", start, end);
    exit(1);
  } //}}}

  // go from box edges to the box centre instead of from centre to edges //{{{
  bool in = BoolOption(argc, argv, "-in"); //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
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

  // <output> - output filename //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

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

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", true, Counts, &BeadType)) {
    exit(1);
  } //}}}

  // -m <name(s)> - specify what molecule types to use //{{{
  int *use = malloc(Counts.TypesOfMolecules*sizeof(int *));
  // if -m not present, use all
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    use[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", &use, Counts, &MoleculeType)) {
    exit(1);
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Use = use[i];
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  // skip till 'pbc' keyword
  char str[1050];
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

  // set maximum/minimum as half a box length in the given direction //{{{
  double range[2];
  switch(axis) {
    case 'x':
      range[0] = 0;
      range[1] = BoxLength.x;
      break;
    case 'y':
      range[0] = 0;
      range[1] = BoxLength.y;
      break;
    case 'z':
      range[0] = 0;
      range[1] = BoxLength.z;
      break;
  } //}}}

  // number of bins //{{{
  int bins[2];
  switch(axis) {
    case 'x':
      bins[0] = BoxLength.y / width + 1;
      bins[1] = BoxLength.z / width + 1;
      break;
    case 'y':
      bins[0] = BoxLength.x / width + 1;
      bins[1] = BoxLength.z / width + 1;
      break;
    case 'z':
      bins[0] = BoxLength.x / width + 1;
      bins[1] = BoxLength.y / width + 1;
      break;
  }
  //}}}

  // allocate memory for density arrays //{{{
  double ***surf = calloc(bins[0], sizeof(double **)); // sum
  int ***values = calloc(bins[0], sizeof(int **)); // number of values
  for (int i = 0; i < bins[0]; i++) {
    surf[i] = calloc(bins[1], sizeof(double *));
    values[i] = calloc(bins[0], sizeof(int *));
    for (int j = 0; j < bins[1]; j++) {
      surf[i][j] = calloc(2, sizeof(double));
      values[i][j] = calloc(2, sizeof(int));
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
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
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "Error: premature end of %s file\n\n", input_coor);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step: %d\n", start);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rStarting step: %d\n", start);
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
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count_vcf, stuff, input_vsf);
      exit(1);
    } //}}}

    RestorePBC(Counts, BoxLength, &Bead);

    // allocate temporary memory //{{{
    double ***temp = calloc(bins[0], sizeof(double **));
    for (int i = 0; i < bins[0]; i++) {
      temp[i] = calloc(bins[1], sizeof(double *));
      for (int j = 0; j < bins[1]; j++) {
        temp[i][j] = calloc(2, sizeof(double));
        if (!in) {
          temp[i][j][0] = 0;
          switch(axis) {
            case 'x':
              temp[i][j][1] = BoxLength.x;
              break;
            case 'y':
              temp[i][j][1] = BoxLength.y;
              break;
            case 'z':
              temp[i][j][1] = BoxLength.z;
              break;
          }
        } else {
          switch(axis) {
            case 'x':
              temp[i][j][0] = BoxLength.x;
              break;
            case 'y':
              temp[i][j][0] = BoxLength.y;
              break;
            case 'z':
              temp[i][j][0] = BoxLength.z;
              break;
          }
          temp[i][j][1] = 0;
        }
      }
    } //}}}

    // calculate surface //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      int btype = Bead[i].Type,
          mol = Bead[i].Molecule;
      if (mol == -1) { // consider only beads in molecules
        continue;
      }
      int mtype = Molecule[mol].Type;
      if (Bead[i].Molecule != -1 && MoleculeType[mtype].Use && BeadType[btype].Use) {
        double coor[3];
        switch(axis) {
          case 'x':
            coor[0] = Bead[i].Position.y;
            coor[1] = Bead[i].Position.z;
            coor[2] = Bead[i].Position.x;
            break;
          case 'y':
            coor[0] = Bead[i].Position.x;
            coor[1] = Bead[i].Position.z;
            coor[2] = Bead[i].Position.y;
            break;
          case 'z':
            coor[0] = Bead[i].Position.x;
            coor[1] = Bead[i].Position.y;
            coor[2] = Bead[i].Position.z;
            break;
        }
        int bin[2];
        bin[0] = coor[0] / width;
        bin[1] = coor[1] / width;
        if (!in) { // go from box centre to edges
          if (coor[2] >= temp[bin[0]][bin[1]][0] && coor[2] >= range[0] && coor[2] <= ((range[0]+range[1])/2)) {
            temp[bin[0]][bin[1]][0] = coor[2];
          }
          if (coor[2] <= temp[bin[0]][bin[1]][1] && coor[2] >= ((range[0]+range[1])/2) && coor[2] <= range[1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
          }
        } else if (in && coor[2] >= range[0] && coor[2] <= range[1]) { // go from box edges to centre
          if (coor[2] <= temp[bin[0]][bin[1]][0]) {
            temp[bin[0]][bin[1]][0] = coor[2];
          }
          if (coor[2] >= temp[bin[0]][bin[1]][1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
          }
        }
      }
    } //}}}

    // add to sums //{{{
    for (int i = 0; i < bins[0]; i++) {
      for (int j = 0; j < bins[1]; j++) {
        if (!in) {
          if (temp[i][j][0] > 0) {
            surf[i][j][0] += temp[i][j][0];
            values[i][j][0]++;
          }
          switch(axis) {
            case 'x':
              if (temp[i][j][1] < BoxLength.x) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'y':
              if (temp[i][j][1] < BoxLength.y) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'z':
              if (temp[i][j][1] < BoxLength.z) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
          }
        } else {
          switch(axis) {
            case 'x':
              if (temp[i][j][0] != (BoxLength.x/2)) {
                surf[i][j][0] += temp[i][j][0];
                values[i][j][0]++;
              }
              if (temp[i][j][1] != (BoxLength.x/2)) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'y':
              if (temp[i][j][0] != (BoxLength.y/2)) {
                surf[i][j][0] += temp[i][j][0];
                values[i][j][0]++;
              }
              if (temp[i][j][1] != (BoxLength.y/2)) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
            case 'z':
              if (temp[i][j][0] < BoxLength.z) {
                surf[i][j][0] += temp[i][j][0];
                values[i][j][0]++;
              }
              if (temp[i][j][1] > 0) {
                surf[i][j][1] += temp[i][j][1];
                values[i][j][1]++;
              }
              break;
          }
        }
      }
    } //}}}

    // free the temporary array //{{{
    for (int i = 0; i < bins[0]; i++) {
      for (int j = 0; j < bins[1]; j++) {
        free(temp[i][j]);
      }
      free(temp[i]);
    }
    free(temp); //}}}

    if (end == count_vcf)
      break;
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rLast Step: %d\n", count_vcf);
    }
  } //}}}

  // write surface to output file //{{{
  // open output file //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    ErrorFileOpen(output, 'w');
    exit(1);
  } //}}}

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print legend
  switch(axis) {
    case 'x':
      fprintf(out, "# (1) y coordinate; (2) z coordinate;");
      break;
    case 'y':
      fprintf(out, "# (1) x coordinate; (2) z coordinate;");
      break;
    case 'z':
      fprintf(out, "# (1) x coordinate; (2) y coordinate;");
      break;
  }
  fprintf(out, " (3) surface 1; (4) surface 2\n");

  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2);
      fprintf(out, " %7.4f", width*(2*j+1)/2);
      if (values[i][j][0] > 0) {
        fprintf(out, " %7.4f", surf[i][j][0]/values[i][j][0]);
      } else {
        fprintf(out, "       ?");
      }
      if (values[i][j][1] > 0) {
        fprintf(out, " %7.4f", surf[i][j][1]/values[i][j][1]);
      } else {
        fprintf(out, "       ?");
      }
      fprintf(out, "\n");
    }
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      free(surf[i][j]);
      free(values[i][j]);
    }
    free(surf[i]);
    free(values[i]);
  }
  free(surf);
  free(values); //}}}

  return 0;
}
