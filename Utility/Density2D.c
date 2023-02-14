#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Density2D utility calculates number \
density for all bead types for a plane perpendicular to the specified axis \
direction (x, y, or z).\
\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <axis> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename \
(either vcf or vtf format)\n");
  fprintf(ptr, "   <width>           width of a single bin \
(in the two planar directions)\n");
  fprintf(ptr, "   <output>          output density file \
(automatic ending '<axis>.rho' added)\n");
  fprintf(ptr, "   <axis>            calculate in plane perpendicular to x, y, \
or z axis\n");
  fprintf(ptr, "   <options>\n");
//fprintf(ptr, "      -n <int>       number of bins to average\n");
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
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
//      strcmp(argv[i], "-n") != 0 &&
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
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  // if vtf, copy to input_vsf
  char *input_vsf = calloc(LINE,sizeof(char));
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
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
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n", start, end);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

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
  if (axis == 'x') {
    fprintf(out, "# columns: (1) y; (2) z;");
  } else if (axis == 'y') {
    fprintf(out, "# columns: (1) x; (2) z;");
  } else {
    fprintf(out, "# columns: (1) x; (2) y;");
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(out, " (%d) %s", i+3, BeadType[i].Name);
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

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // number of bins //{{{
  double max_dist[2];
  if (axis == 'x') {
    max_dist[0] = BoxLength.y;
    max_dist[1] = BoxLength.z;
  } else if (axis == 'y') {
    max_dist[0] = BoxLength.x;
    max_dist[1] = BoxLength.z;
  } else {
    max_dist[0] = BoxLength.x;
    max_dist[1] = BoxLength.y;
  }
  int bins[2];
  bins[0] = ceil(max_dist[0] / width);
  bins[1] = ceil(max_dist[1] / width); //}}}

  // allocate memory for density arrays //{{{
  long int ***rho = malloc(Counts.TypesOfBeads*sizeof(long int **));
  long int ***rho_2 = malloc(Counts.TypesOfBeads*sizeof(long int **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    rho[i] = calloc(bins[0],sizeof(long int *));
    rho_2[i] = calloc(bins[0],sizeof(long int *));
    for (int j = 0; j < bins[0]; j++) {
      rho[i][j] = calloc(bins[1],sizeof(long int));
      rho_2[i][j] = calloc(bins[1],sizeof(long int));
    }
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
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "Error: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
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
  // is the vcf file continuing?
  if (ErrorDiscard(start, count, input_coor, vcf)) {
    exit(1);
  }
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
      ErrorCoorRead(input_coor, test, count_vcf, stuff);
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
    int ***temp_rho = malloc(Counts.TypesOfBeads*sizeof(int **));
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      temp_rho[i] = calloc(bins[0],sizeof(int *));
      for (int j = 0; j < bins[0]; j++) {
        temp_rho[i][j] = calloc(bins[1],sizeof(int));
      }
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
          int j = Bead[i].Position.y / width;
          int k = Bead[i].Position.z / width;
          temp_rho[Bead[i].Type][j][k]++;
        } else if (axis == 'y') {
          int j = Bead[i].Position.x / width;
          int k = Bead[i].Position.z / width;
          temp_rho[Bead[i].Type][j][k]++;
        } else {
          int j = Bead[i].Position.x / width;
          int k = Bead[i].Position.y / width;
          temp_rho[Bead[i].Type][j][k]++;
        }
      }
    } //}}}

    // add from temporary density array to global density arrays //{{{
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < bins[0]; k++) {
        for (int l = 0; l < bins[1]; l++) {
          rho[j][k][l] += temp_rho[j][k][l];
          rho_2[j][k][l] += SQR(temp_rho[j][k][l]);
        }
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      for (int j = 0; j < bins[0]; j++) {
        free(temp_rho[i][j]);
      }
      free(temp_rho[i]);
    }
    free(temp_rho); //}}}

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

  // write densities to output file(s) //{{{
  if ((out = fopen(output_rho, "a")) == NULL) {
    ErrorFileOpen(output_rho, 'a');
    exit(1);
  }

  // calculate rdf
  double volume = SQR(width);
  if (axis == 'x') {
    volume *= BoxLength.x;
  } else if (axis == 'y') {
    volume *= BoxLength.y;
  } else {
    volume *= BoxLength.z;
  }
  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      fprintf(out, "%7.3f", width*(i+0.5)); // TODO: 2D, check Surface
      fprintf(out, "%7.3f", width*(j+0.5)); // TODO: 2D, check Surface
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        double temp_rho = 0;
        temp_rho += rho[k][i][j] / (volume * count);
        // print average value to output file
        fprintf(out, " %10f", temp_rho);
      }
      putc('\n',out);
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
    for (int j = 0; j < bins[0]; j++) {
      free(rho[i][j]);
      free(rho_2[i][j]);
    }
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
