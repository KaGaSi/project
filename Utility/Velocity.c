#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
TBA\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <bead name(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>                 input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <width>                 width of a single distribution bin\n");
  fprintf(ptr, "   <output>                output file with distribution of bond lengths\n");
  fprintf(ptr, "   <bead name(s)>          bead name(s) to calculate velocities for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>            starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>             ending timestep for calculation\n");
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

  if (argc < req_args) {
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
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0) {

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

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (!IsPosDouble(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - file name with bond length distribution //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

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

  // number of steps to skip per one used //{{{
  int skip = 0;
  if (IntegerOption(argc, argv, "-sk", &skip)) {
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

  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Use = false;
  }
  // <bead names> - names of bead types to use //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);
    if (type == -1) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - non-existent bead name \033[1;33m%s\033[1;31m\n", input_coor, argv[count]);
      fprintf(stderr, "\033[0m");
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}
  int total_beads = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Use) {
      total_beads += BeadType[i].Number;
    }
  }

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vtf, vcf, struct_lines); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // number of bins
  VECTOR bins;
  bins.x = ceil(BoxLength.x / width);
  bins.y = ceil(BoxLength.y / width);
  bins.z = ceil(BoxLength.z / width);

  // arrays for distributions //{{{
  VECTOR *velocity_x = calloc(bins.x, sizeof(VECTOR)),
         *velocity_y = calloc(bins.y, sizeof(VECTOR)),
         *velocity_z = calloc(bins.z, sizeof(VECTOR));
  long int *count_x = calloc(bins.x, sizeof(long int)),
           *count_y = calloc(bins.y, sizeof(long int)),
           *count_z = calloc(bins.z, sizeof(long int));
  //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  count = SkipCoorSteps(vcf, input_coor, Counts, start, silent);

  FILE *out;
  if ((out = fopen("temperature.txt", "w")) == NULL) {
    ErrorFileOpen("temperature.txt", 'w');
    exit(1);
  }
  fprintf(out, "# (1) step, (2) temperature\n");
  fclose(out);
  if ((out = fopen("temperature_correct.txt", "w")) == NULL) {
    ErrorFileOpen("temperature_correct.txt", 'w');
    exit(1);
  }
  fprintf(out, "# (1) step, (2) temperature\n");
  fclose(out);

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

    double KE = 0;
    VECTOR temp_step = {0, 0, 0},
           *velocity_x_step = calloc(bins.x, sizeof(VECTOR)),
           *velocity_y_step = calloc(bins.y, sizeof(VECTOR)),
           *velocity_z_step = calloc(bins.z, sizeof(VECTOR));
    int *count_z_step = calloc(bins.z, sizeof(INTVECTOR));
    for (int i = 0; i < Counts.Beads; i++) {
      int btype = Bead[i].Type;
      temp_step.x += 0.5 * BeadType[btype].Mass * SQR(Bead[i].Velocity.x);
      temp_step.y += 0.5 * BeadType[btype].Mass * SQR(Bead[i].Velocity.y);
      temp_step.z += 0.5 * BeadType[btype].Mass * SQR(Bead[i].Velocity.z);
      KE += 0.5 * BeadType[btype].Mass * (SQR(Bead[i].Velocity.x) + SQR(Bead[i].Velocity.y) + SQR(Bead[i].Velocity.z));

      if (BeadType[btype].Use) {
        if (Bead[i].Position.x == BoxLength.x) {
          Bead[i].Position.x -= BoxLength.x;
        }
        int j = Bead[i].Position.x / width;
        if (j < bins.x) {
          velocity_x_step[j].x += Bead[i].Velocity.x;
          velocity_x_step[j].y += Bead[i].Velocity.y;
          velocity_x_step[j].z += Bead[i].Velocity.z;
          count_x[j]++;
        } else {
          printf("\nx: %d %lf %lf\n", j, bins.x, Bead[i].Position.x);
        }

        if (Bead[i].Position.y == BoxLength.y) {
          Bead[i].Position.y -= BoxLength.y;
        }
        j = Bead[i].Position.y / width;
        if (j < bins.y) {
          velocity_y_step[j].x += Bead[i].Velocity.x;
          velocity_y_step[j].y += Bead[i].Velocity.y;
          velocity_y_step[j].z += Bead[i].Velocity.z;
          count_y[j]++;
        } else {
          printf("\ny: %d %lf %lf\n", j, bins.y, Bead[i].Position.y);
        }

        if (Bead[i].Position.z == BoxLength.z) {
          Bead[i].Position.z -= BoxLength.z;
        }
        j = Bead[i].Position.z / width;
        if (j < bins.z) {
          velocity_z_step[j].x += Bead[i].Velocity.x;
          velocity_z_step[j].y += Bead[i].Velocity.y;
          velocity_z_step[j].z += Bead[i].Velocity.z;
          count_z[j]++;
          count_z_step[j]++;
        } else {
          printf("\nz: %d %lf %lf\n", j, bins.z, Bead[i].Position.z);
        }
      }
    }

    for (int i = 0; i < bins.x; i++) {
      velocity_x[i].x += velocity_x_step[i].x;
      velocity_x[i].y += velocity_x_step[i].y;
      velocity_x[i].z += velocity_x_step[i].z;
    }
    for (int i = 0; i < bins.y; i++) {
      velocity_y[i].x += velocity_y_step[i].x;
      velocity_y[i].y += velocity_y_step[i].y;
      velocity_y[i].z += velocity_y_step[i].z;
    }
    for (int i = 0; i < bins.z; i++) {
      velocity_z[i].x += velocity_z_step[i].x;
      velocity_z[i].y += velocity_z_step[i].y;
      velocity_z[i].z += velocity_z_step[i].z;
    }

    double KE_correct = 0;
    VECTOR temp_step_correct = {0, 0, 0};
    for (int i = 0; i < Counts.Beads; i++) {
      if (Bead[i].Position.z == BoxLength.z) {
        Bead[i].Position.z -= BoxLength.z;
      }
      int j = Bead[i].Position.z / width;

      int btype = Bead[i].Type;

      if (count_z_step[j] > 0) {
        temp_step_correct.x += 0.5 * BeadType[btype].Mass * SQR(Bead[i].Velocity.x - velocity_z_step[j].x / count_z_step[j]);
        temp_step_correct.y += 0.5 * BeadType[btype].Mass * SQR(Bead[i].Velocity.y);
        temp_step_correct.z += 0.5 * BeadType[btype].Mass * SQR(Bead[i].Velocity.z);
        KE_correct += 0.5 * BeadType[btype].Mass * (SQR(Bead[i].Velocity.x - velocity_z_step[j].x / count_z_step[j])
                                                  + SQR(Bead[i].Velocity.y)
                                                  + SQR(Bead[i].Velocity.z));
      }
    }

    temp_step.x = 2 * temp_step.x / total_beads;
    temp_step.y = 2 * temp_step.y / total_beads;
    temp_step.z = 2 * temp_step.z / total_beads;
    temp_step_correct.x = 2 * temp_step_correct.x / total_beads;
    temp_step_correct.y = 2 * temp_step_correct.y / total_beads;
    temp_step_correct.z = 2 * temp_step_correct.z / total_beads;
    double temperature = 2 * KE / (3 * total_beads);
    double temperature_correct = 2 * KE_correct / (3 * total_beads);
    FILE *out;
    if ((out = fopen("temperature.txt", "a")) == NULL) {
      ErrorFileOpen("temperature.txt", 'a');
      exit(1);
    }
    fprintf(out, "%d %lf %lf %lf %lf\n", count_vcf, temperature, temp_step.x, temp_step.y, temp_step.z);
    fclose(out);
    if ((out = fopen("temperature_correct.txt", "a")) == NULL) {
      ErrorFileOpen("temperature_correct.txt", 'a');
      exit(1);
    }
    fprintf(out, "%d %lf %lf %lf %lf\n", count_vcf, temperature_correct, temp_step_correct.x, temp_step_correct.y, temp_step_correct.z);
    fclose(out);

    free(velocity_x_step);
    free(velocity_y_step);
    free(velocity_z_step);
    free(count_z_step);

    // -e option - exit main loop if last step is done
    if (end == count_vcf) {
      break;
    }
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

  // write distributions to output files
  // open files //{{{
  // x-slicing distribution
  FILE *out_x;
  char str[1050];
  sprintf(str, "%s-x.txt", output);
  if ((out_x = fopen(str, "w")) == NULL) {
    ErrorFileOpen(str, 'w');
    exit(1);
  }
  // y-slicing distribution
  FILE *out_y;
  sprintf(str, "%s-y.txt", output);
  if ((out_y = fopen(str, "w")) == NULL) {
    ErrorFileOpen(str, 'w');
    exit(1);
  }
  // z-slicing distribution
  FILE *out_z;
  sprintf(str, "%s-z.txt", output);
  if ((out_z = fopen(str, "w")) == NULL) {
    ErrorFileOpen(str, 'w');
    exit(1);
  } //}}}
  fprintf(out_x, "# (1) x-coordinate; (2) <v(x)>; (3) <v(y)>; (4) <v(z)>\n");
  fprintf(out_y, "# (1) y-coordinate; (2) <v(x)>; (3) <v(y)>; (4) <v(z)>\n");
  fprintf(out_z, "# (1) z-coordinate; (2) <v(x)>; (3) <v(y)>; (4) <v(z)>\n");
  for (int i = 0; i < bins.x; i++) {
    if (count_x[i] > 0) {
      fprintf(out_x, "%7.4f %lf %lf %lf\n", width*(2*i+1)/2, velocity_x[i].x/count_x[i],
                                                             velocity_x[i].y/count_x[i],
                                                             velocity_x[i].z/count_x[i]);
    }
  }
  for (int i = 0; i < bins.y; i++) {
    if (count_y[i] > 0) {
      fprintf(out_y, "%7.4f %lf %lf %lf\n", width*(2*i+1)/2, velocity_y[i].x/count_y[i],
                                                             velocity_y[i].y/count_y[i],
                                                             velocity_y[i].z/count_y[i]);
    }
  }
  for (int i = 0; i < bins.z; i++) {
    if (count_z[i] > 0) {
      fprintf(out_z, "%7.4f %lf %lf %lf\n", width*(2*i+1)/2, velocity_z[i].x/count_z[i],
                                                             velocity_z[i].y/count_z[i],
                                                             velocity_z[i].z/count_z[i]);
    }
  }
  fclose(out_x);
  fclose(out_y);
  fclose(out_z);

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  free(velocity_x);
  free(velocity_y);
  free(velocity_z);
  free(count_x);
  free(count_y);
  free(count_z);
  //}}}

  return 0;
}
