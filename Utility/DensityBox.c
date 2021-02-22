#include "../AnalysisTools.h"

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
  fprintf(ptr, "   %s <input> <width> <output> <axis> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename (either vcf or vtf format)\n");
  fprintf(ptr, "   <width>           width of a single bin\n");
  fprintf(ptr, "   <output>          output density file (automatic ending '<axis>.rho' added)\n");
  fprintf(ptr, "   <axis>            calculate along x, y, or z axis\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -n <int>       number of bins to average\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       number of timestep to end with\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
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
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-x") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
  char *input_vsf = calloc(LINE,sizeof(char));
  strcpy(input_coor, argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <width> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (!IsPosDouble(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.rho> - filename with bead densities //{{{
  char output_rho[LINE];
  strcpy(output_rho, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // starting & ending timesteps //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  }
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  }
  ErrorStartEnd(start, end); //}}}

  // number of bins to average //{{{
  int avg = 1;
  if (IntegerOption(argc, argv, "-n", &avg)) {
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  VECTOR BoxLength; // couboid box dimensions
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &BoxLength, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule);
  free(input_vsf); //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // <axis> - x, y, or z //{{{
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {

    axis = argv[count][0];

    if (axis != 'x' && axis != 'y' && axis != 'z') {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m<axis>\033[1;31m must be 'x', 'y', or 'z'\n\n");
      fprintf(stderr, "\033[0m");
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

  // print command to output file
  putc('#', out);
  PrintCommand(out, argc, argv);

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

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(stdout, " %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  count = SkipCoorSteps(vcf, input_coor, Counts, start, silent);

  // main loop //{{{
  count = 0; // count timesteps in the main loop
  int count_vcf = start - 1; // count timesteps from the beginning
  while (true) {
    count++;
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    ReadVcfCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

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
    // if there's no additional timestep, exit the while loop
    bool rubbish; // not used
    if (ReadTimestepPreamble(&rubbish, input_coor, vcf, &stuff, false) == -1) {
      break;
    }
  }
  fclose(vcf);

  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
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
      double temp_rho = 0;
      // sum densities to be averaged
      for (int k = 0; k < avg; k++) {
        temp_rho += rho[j][i+k] / (volume * count);
      }
      // print average value to output file
      fprintf(out, " %10f", temp_rho);
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
