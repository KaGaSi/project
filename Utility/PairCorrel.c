#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
PairCorrel utility calculates pair correlation function for specified \
bead types. All pairs of bead types (including same type pairs) are \
calculated - given A and B types, pcf between A-A, A-B and B-B are \
calculated.\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <width> <output> <bead(s)> ", cmd);
  fprintf(ptr, "[options]\n\n");

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a distribution bin\n");
  fprintf(ptr, "<output>            output file with pair correlation "
          "function(s)\n");
  fprintf(ptr, "<bead(s)>           bead name(s) for calculation "
          "(optional and ignored if '--all' is used)\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  --all             use all bead types "
          "(overwrites <bead(s)>)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool all; // --all
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 1, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, false, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "--all");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width = -1;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <output> - filename with pcf(s)
  char fout_pcf[LINE] = "";
  s_strcpy(fout_pcf, argv[++count], LINE);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  opt->all = BoolOption(argc, argv, "--all");
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  double *box = System.Box.Length;

  // <bead(s)> - names of bead types to use //{{{
  if (opt->all) {
    for (int i = 0; i < Count->BeadType; i++) {
      System.BeadType[i].Flag = true;
    }
  } else {
    for (int i = 0; i < Count->BeadType; i++) {
      System.BeadType[i].Flag = false;
    }
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType(argv[count], System);
      if (type == -1) {
        ErrorBeadType(argv[count], System);
        exit(1);
      }
      if (System.BeadType[type].Flag) {
        snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
                 ErrYellow(), argv[count], ErrCyan());
        PrintWarning();
      }
      System.BeadType[type].Flag = true;
    }
  } //}}}

  // write initial stuff to output pcf file //{{{
  PrintByline(fout_pcf, argc, argv);
  FILE *out = OpenFile(fout_pcf, "a");
  fprintf(out, "# (1) distance");
  // print bead type names to output file //{{{
  count = 1;
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = i; j < Count->BeadType; j++) {
      if (System.BeadType[i].Flag && System.BeadType[j].Flag) {
        count++;
        fprintf(out, " (%d) %s-%s", count, System.BeadType[i].Name,
                                           System.BeadType[j].Name);
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out); //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  int bins = Max3(box[0], box[1], box[2]) / width;
  double max_dist = 0.5 * Min3(box[0], box[1], box[2]);

  // allocate memory //{{{
  // array counting number of pairs
  long int **counter = calloc(Count->BeadType, sizeof *counter);
  // pair correlation function
  int ***pcf = malloc(Count->BeadType * sizeof **pcf);
  for (int i = 0; i < Count->BeadType; i++) {
    counter[i] = calloc(Count->BeadType, sizeof *counter[i]);
    pcf[i] = malloc(Count->BeadType * sizeof *pcf[i]);
    for (int j = 0; j < Count->BeadType; j++) {
      pcf[i][j] = calloc(bins, sizeof *pcf[i][j]);
    }
  } //}}}

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count timesteps from the beginning
      count_used = 0, // count steps used for calculation
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      for (int i = 0; i < Count->BeadCoor; i++) { //{{{
        int id_i = System.BeadCoor[i];
        BEAD *b_i = &System.Bead[id_i];
        if (System.BeadType[b_i->Type].Flag) {
          for (int j = (i + 1); j < Count->BeadCoor; j++) {
            int id_j = System.BeadCoor[j];
            BEAD *b_j = &System.Bead[id_j];
            if (System.BeadType[b_j->Type].Flag) {
              int btype_i = b_i->Type;
              int btype_j = b_j->Type;
              if (btype_i > btype_j) {
                SwapInt(&btype_i, &btype_j);
              }
              counter[btype_i][btype_j]++;
              double dist[3];
              Distance(b_i->Position, b_j->Position, box, dist);
              dist[0] = VectLength(dist);
              if (dist[0] < max_dist) {
                int l = dist[0] / width;
                pcf[btype_i][btype_j][l]++;
              }
            }
          }
        }
      } //}}}
      //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // write data to output file(s) //{{{
  out = OpenFile(fout_pcf, "a");

  // calculate pcf
  for (int j = 0; j < bins; j++) {

    // calculate volume of every shell that will be averaged
    double shell;
    shell = 4.0 / 3 * PI * Cube(width) * (Cube(j + 1) - Cube(j));
    fprintf(out, "%8.5f", width * (2 * j + 1)/2);

    for (int k = 0; k < Count->BeadType; k++) {
      for (int l = k; l < Count->BeadType; l++) {
        if (System.BeadType[k].Flag && System.BeadType[l].Flag) {
          double norm_factor = System.Box.Volume / counter[k][l] /shell;
          fprintf(out, " %10f", pcf[k][l][j] * norm_factor);
        }
      }
    }
    putc('\n',out);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      free(pcf[i][j]);
    }
    free(pcf[i]);
    free(counter[i]);
  }
  free(pcf);
  free(counter);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
