#include "../AnalysisTools.h"

// const double wall_dist1 = 0;
// const double wall_dist2 = 17.241-9.055;
double wall_dist1 = 17.241 - (9.85 - 9.055) / 2;
double wall_dist2 = 17.241;
int axis = 2; // wall orientation: 0..x, 1..y, 2..z

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Apart from classic Selected, this 'Extra' version restricts coordinates \
of the saved beads in a specified axis direction to two regions symmetrically \
aligned according to the two provided coordinates - e.g., with z axis box \
size 10 and provided coordinates 1 2, only beads with z coordinates <1,2> \
or <8,9> are saved. \n\n\
Selected creates new coordinate file in the extension-specified format \
that contains specified bead and/or molecule types. Periodic boundary \
conditions can be either stripped away or applied (which happens first if both \
'--join' and '--wrap' options are used).\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <output> <axis> 2×<coor> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output coordinate file\n");
  fprintf(ptr, "<axis>              x, y, or z\n");
  fprintf(ptr, "2×<coor>            interval for constraining coordinates\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -bt <bead type>   bead types to exclude\n");
  fprintf(ptr, "  -mt <mol type>    molecule types to exclude\n");
  fprintf(ptr, "  --reverse         reverse <bead name(s)>, i.e., save the"
          " specified bead types (-bt and/or -mt is needed)\n");
  fprintf(ptr, "  --join            join molecules (remove pbc)\n");
  fprintf(ptr, "  --wrap            wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "  -n <int(s)>       save only specified timesteps"
          "(--last overrides this option)\n");
  fprintf(ptr, "  --last            use only the last step"
          "(-st/-e/-n options are ignored)\n");
  fprintf(ptr, "  -sc <float>       divide all coordinates by given value\n");
  fprintf(ptr, "  -m <3xfloat>      move all coordinates by given vector "
          "(-sc option is applied first)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool bt, mt,               // -bt/-mt
       reverse,              // --reverse
       join, wrap, last;     // --join --wrap --last
  int n_save[100], n_number; // -n
  double scale,              // -sc
         move[3];            // -m
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

static void ScaleCoordinates(SYSTEM *System, const double scale) { //{{{
  if (scale > 1) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[id].Position[dd] /= 50;
      }
    }
  }
} //}}}
static void MoveCoordinates(SYSTEM *System, const double move[3]) { //{{{
  if (move[0] != 0 || move[1] != 0 || move[2] != 0) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[id].Position[dd] += move[dd];
      }
    }
  }
} //}}}

static void WriteWall(SYSTEM System, bool *write, bool *write_wall) { //{{{
  double dist1[2], dist2[2];
  dist1[0] = wall_dist1;
  dist1[1] = System.Box.Length[axis] - wall_dist1;
  dist2[0] = wall_dist2;
  dist2[1] = System.Box.Length[axis] - wall_dist2;
  // printf("\ndist1[0] %lf dist1[1] %lf\n", dist1[0], dist1[1]);
  // printf("dist2[0] %lf dist2[1] %lf\n", dist2[0], dist2[1]);
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    write_wall[id] = false;
    BEAD *b = &System.Bead[id];
    if (write[id]) {
      if (b->Position[axis] >= dist1[0] && b->Position[axis] <= dist2[0]) {
        write_wall[id] = true;
      }
      if (b->Position[axis] <= dist1[1] && b->Position[axis] >= dist2[1]) {
        write_wall[id] = true;
      }
    }
  }
}
static void InvertWriteWall(SYSTEM System, bool *write_wall) {
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    write_wall[id] = !write_wall[id];
  }
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 8, count = 0,
      req_arg = 5;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version", "-bt", "-mt", "--reverse", "--join", "--wrap", "-n",
              "--last", "-sc", "-m");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <output> - output coordinate file
  FILE_TYPE fout;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name);

  if (argv[++count][0] == 'x') {
    axis = 0;
  } else if (argv[count][0] == 'y') {
    axis = 1;
  } else if (argv[count][0] == 'z') {
    axis = 2;
  } else {
    err_msg("<axis> must by 'x', 'y', or 'z'");
    PrintError();
    PrintCommand(stderr, argc, argv);
    Help(argv[0], true, common, option);
    exit(1);
  }

  if (!IsRealNumber(argv[++count], &wall_dist1)) {
    err_msg("<coor> must by a real number");
    PrintError();
    PrintCommand(stderr, argc, argv);
    Help(argv[0], true, common, option);
    exit(1);
  }
  if (!IsRealNumber(argv[++count], &wall_dist2)) {
    err_msg("<coor> must by a real number");
    PrintError();
    PrintCommand(stderr, argc, argv);
    Help(argv[0], true, common, option);
    exit(1);
  }


  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  opt->reverse = BoolOption(argc, argv, "--reverse");
  opt->join = BoolOption(argc, argv, "--join");
  opt->wrap = BoolOption(argc, argv, "--wrap");
  opt->last = BoolOption(argc, argv, "--last");
  if (!OneNumberOption(argc, argv, "-sc", &opt->scale, 'd')) {
    opt->scale = 1;
  }
  if (!ThreeNumbersOption(argc, argv, "-m", opt->move, 'd')) {
    for (int dd = 0; dd < 3; dd++) {
      opt->move[dd] = 0;
    }
  }
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  if (opt->join && Count->Molecule == 0) {
    err_msg("no molecules to join");
    PrintWarning();
  }

  // specify beads to save (possibly using -bt and/or -mt options) //{{{
  /*
   * reverse=true ... save only the specified species
   * reverse=false ... exclude the specified species
   *
   * First setting all bead/molecule types to !reverse and then adjusting this
   * if -bt/-mt options are present correctly specifies which
   * bead/molecule types to save
   */
  // auxiliary arrays holding which bead/molecule types to save
  bool *write_bt = malloc(Count->BeadType * sizeof *write_bt),
       *write_mt = malloc(Count->MoleculeType * sizeof *write_mt);
  // first assume all bead types are saved/excluded based on --reverse option...
  InitBoolArray(write_bt, Count->BeadType, !opt->reverse);
  // ... then adjust if -bt option is present
  opt->bt = BeadTypeOption(argc, argv, "-bt", opt->reverse, write_bt, System);
  // first assume all molecule types are saved/excluded...
  InitBoolArray(write_mt, Count->MoleculeType, !opt->reverse);
  // ... then adjust if -mt option is present
  opt->mt = MoleculeTypeOption(argc, argv, "-mt", opt->reverse,
                               write_mt, System);
  // array for holding which beads to save/exclude
  bool *write = malloc(Count->Bead * sizeof *write);
  // first assume all are saved/excluded...
  InitBoolArray(write, Count->Bead, !opt->reverse);
  // then check possible -mt option...
  if (opt->mt) {
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < mt->Number; j++) {
        int mol = mt->Index[j];
        for (int k = 0; k < mt->nBeads; k++) {
          int id = System.Molecule[mol].Bead[k];
          if (write_mt[i] == opt->reverse) {
            write[id] = opt->reverse; // save/exclude based on --reverse
          }
        }
      }
    }
  }
  // ... and possible -bt option
  if (opt->bt) {
    for (int i = 0; i < Count->Bead; i++) {
      int type = System.Bead[i].Type;
      if (write_bt[type] == opt->reverse) {
        write[i] = opt->reverse;
      }
    }
  }
  // free the auxiliary arrays
  free(write_mt);
  free(write_bt); //}}}

  // '-n' option - specify timestep ids //{{{
  opt->n_number = -1;
  InitIntArray(opt->n_save, 100, 0);
  NumbersOption(argc, argv, 100, "-n", &opt->n_number, opt->n_save, 'i');
  // ignore -st/-e/-sk when -n is used
  if (opt->n_number != -1) {
    opt->c.start = 1;
    opt->c.end = -1;
    opt->c.skip = 1;
  }
  // SortArrayInt(opt->n_save, opt->n_number, 0);
  SortArray(opt->n_save, opt->n_number, 0, 'i'); //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // print initial stuff to output coordinate file //{{{
  if (fout.type == VCF_FILE) {
    PrintByline(fout.name, argc, argv);
  } else if (fout.type == VTF_FILE) {
    WriteStructure(fout, System, -1, false, argc, argv);
  } else { // ensure it's a new file
    FILE *out = OpenFile(fout.name, "w");
    fclose(out);
  } //}}}

  bool *write_wall = malloc(Count->Bead * sizeof *write_wall);
  InitBoolArray(write_wall, Count->Bead, true);

  FILE *fr = OpenFile(in.coor.name, "r");
  // main loop //{{{
  // file pointers for finding the last valid step
  fpos_t *position = calloc(1, sizeof *position);
  // save line count at every fgetpos()
  int *bkp_line_count = calloc(1, sizeof *bkp_line_count);
  int n_opt_count = 0, // count saved steps if -n option is used
      count_coor = 0,  // count steps in the vcf file
      count_saved = 0, // count steps in output file
      line_count = 0;  // count lines in the coor file
  while (true) {
    if (opt->last) {
      if (!opt->c.silent && isatty(STDOUT_FILENO)) {
        count_coor++;
        fprintf(stdout, "\rDiscarding step: %d", count_coor);
      }
    } else {
      PrintStep(&count_coor, opt->c.start, opt->c.silent);
    }
    position = s_realloc(position, count_coor * sizeof *position);
    fgetpos(fr, &position[count_coor-1]);
    bkp_line_count = s_realloc(bkp_line_count, count_coor *
                               sizeof *bkp_line_count);
    bkp_line_count[count_coor-1] = line_count;
    // decide whether this timestep is to be saved //{{{
    bool use = false;
    // no -n option - possibly use only if --last not present
    if (opt->n_number == -1 && !opt->last && UseStep(opt->c, count_coor)) {
      use = true;
    // -n option is used - save the timestep if it's in the list
    } else if (n_opt_count < opt->n_number &&
               opt->n_save[n_opt_count] == count_coor) {
      use = true;
      n_opt_count++;
    } //}}}
    if (use) { // read and write the timestep, if it should be saved //{{{
      if (fout.type == LDATA_FILE && count_saved == 1) {
        err_msg("only one timestep can be saved to lammps data file");
        PrintWarnFile(fout.name, "\0", "\0");
        count_coor--;
        break;
      }
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_saved++;
      ScaleCoordinates(&System, opt->scale);
      MoveCoordinates(&System, opt->move);
      WrapJoinCoordinates(&System, opt->wrap, opt->join);
      WriteWall(System, write, write_wall);
      WriteTimestep(fout, System, count_coor, write_wall, argc, argv);
      //}}}
    } else { // skip the timestep, if it shouldn't be saved //{{{
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // decide whether to exit the main loop //{{{
    /* break the loop if
     *    1) all timesteps in the -n option are saved (and --last isn't used)
     *    or
     *    2) end timestep was reached (-e option)
     */
    if (!opt->last && // never break when --last is used
        (n_opt_count == opt->n_number || // 1)
         count_coor == opt->c.end)) {    // 2)
      break;
    } //}}}
  } //}}}
  // if --last option is used, read & save the last timestep //{{{
  if (opt->last) {
    /* To through all saved file positions (last to first) and save a the first
     * valid step encountered.
     * Start at count_coor as the saved position is at the beginning of the last
     * timestep to be skipped, and count_coor-- is used beore quitting the while
     * loop.
     */
    count_coor--; // decrement as the last step should be eof
    for (int i = (count_coor); i >= 0; i--) {
      fsetpos(fr, &position[i]);
      line_count = bkp_line_count[i];
      if (ReadTimestep(in, fr, &System, &line_count)) {
        ScaleCoordinates(&System, opt->scale);
        WrapJoinCoordinates(&System, opt->wrap, opt->join);
        WriteWall(System, write, write_wall);
        WriteTimestep(fout, System, count_coor, write_wall, argc, argv);
        if (!opt->c.silent) {
          if (isatty(STDOUT_FILENO)) {
            fprintf(stdout, "\r                          \r");
          }
          fprintf(stdout, "Saved Step: %d\n", i+1);
          fflush(stdout);
        }
        break;
      } else {
        snprintf(ERROR_MSG, LINE, "disregarding step %s%d%s", ErrYellow(),
                 i + 1, ErrCyan());
        PrintWarnFile(in.coor.name, "\0", "\0");
      }
    } //}}}
  } else if (count_coor == 0) { // error - input file without a valid timestep //{{{
    remove(fout.name);
    err_msg("no valid timestep found");
    PrintErrorFile(in.coor.name, "\0", "\0"); //}}}
  } else if (opt->c.start > count_coor) { // warn if no timesteps were written //{{{
    remove(fout.name);
    err_msg("no coordinates written (starting timestep is higher "
            "than the total number of timesteps)");
    PrintWarning(); //}}}
  } else if (!opt->c.silent) { // print last step count? //{{{
    if (isatty(STDOUT_FILENO)) {
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (saved %d)\n", count_coor, count_saved);
    fflush(stdout);
  } //}}}
  fclose(fr);

  // free memory
  FreeSystem(&System);
  free(write);
  free(write_wall);
  free(position);
  free(bkp_line_count);
  free(opt);

  return 0;
}
