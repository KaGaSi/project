#include "Options.h"

// STATIC DECLARATIONs
static void SilentOption(const int argc, char **argv,
                         bool *verbose, bool *silent);
static bool VersionOption(const int argc, char **argv);
// some repeated warnings/errors
static void MissingFilenameError(const int argc, char **argv,
                                 const char *opt, const int i);
static bool TooManyArgsWarn(const int max, const int n,
                            const char *opt, int *count);
static void ArgumentNumberErr(const int count, const int n, const char *opt);
static void ArgumentMissingErr(const int n, const char *opt);

// THE VISIBLE FUNCTIONS
// version/help printing and initial check of provided options //{{{
int OptionCheck(const int argc, char **argv, const int req, const int common,
                const int all, const bool check_extra,
                char opt[all][OPT_LENGTH], ...) {
  snprintf(argv[0], LINE, "%s", BareCommand(argv[0]));
  // copy options to an array
  va_list args;
  va_start(args, opt);
  // Loop through the variable arguments, copy each string
  for (int i = 0; i < all; i++) {
    const char *src = va_arg(args, const char *);
    s_strcpy(opt[i], src, OPT_LENGTH);
  }
  va_end(args);
  // --version option?
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  // --help option?
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--help") == 0) {
      Help(argv[0], false, common, opt);
      exit(0);
    }
  }
  // correct number of mandatory options?
  int count = 0;
  while ((count + 1) < argc &&
         // there may be '-' as a mandatory argument
         (argv[count+1][0] != '-' || strnlen(argv[count+1], LINE) == 1)) {
    count++;
  }
  if (count < req) {
    ErrorArgNumber(count, req);
    PrintCommand(stderr, argc, argv);
    Help(argv[0], true, common, opt);
    exit(1);
  }
  // all options exist?
  for (int i = (count+1); i < argc; i++) {
    bool valid = false;
    for (int j = 0; j < all; j++) {
      double value;
      if (argv[i][0] != '-' || // assumes an argument to some option
          IsRealNumber(argv[i], &value) || // assumes negative numeric argument
          strcmp(argv[i], opt[j]) == 0) { // an option
        valid = true;
        break;
      }
    }
    if (!valid) {
      ErrorOption(argv[i]);
      PrintCommand(stderr, argc, argv);
      Help(argv[0], true, common, opt);
      exit(1);
    }
  }
  // warn if extra arguments (between required ones and options)
  if (check_extra && req != count) {
    char extra[LINE] = "\0";
    for (int i = (req + 1); i <= count; i++) {
      char cpy[LINE];
      s_strcpy(cpy, extra, LINE);
      if (snprintf(extra, LINE, "%s %s", cpy, argv[i]) < 0) {
        ErrorSnprintf();
      }
    }
    if (snprintf(ERROR_MSG, LINE, "command line arguments%s%s%s have no effect",
                 ErrYellow(), extra, ErrCyan()) < 0) {
      ErrorSnprintf();
    }
    PrintWarning();
  }
  return count;
} //}}}
// print help for common options //{{{
void CommonHelp(const bool error, const int n,
                const char option[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
  }
  for (int i = 0; i < n; i++) {
    if (strcmp(option[i], "-i") == 0) {
      fprintf(ptr, "  -i <name>         input structure file if different "
                   "than the coordinate file\n");
    } else if (strcmp(option[i], "-st") == 0) {
      fprintf(ptr, "  -st <int>         starting timestep for calculation\n");
    } else if (strcmp(option[i], "-e") == 0) {
      fprintf(ptr, "  -e <end>          ending timestep for calculation\n");
    } else if (strcmp(option[i], "-sk") == 0) {
      fprintf(ptr, "  -sk <int>         leave out every 'skip' steps\n");
    } else if (strcmp(option[i], "--variable") == 0) {
      fprintf(ptr, "  --variable        vtf coordinate file with indexed "
              "timesteps with varying number of beads\n");
    } else if (strcmp(option[i], "-ltrj") == 0) {
      fprintf(ptr, "  -ltrj <int>       does lammpstrj ids go from 0 or 1?\n");
    } else if (strcmp(option[i], "--verbose") == 0) {
      fprintf(ptr, "  --verbose         verbose output\n");
    } else if (strcmp(option[i], "--silent") == 0) {
      fprintf(ptr, "  --silent          no output (overrides --verbose)\n");
    } else if (strcmp(option[i], "--help") == 0) {
      fprintf(ptr, "  --help            print this help and exit\n");
    } else if (strcmp(option[i], "--version") == 0) {
      fprintf(ptr, "  --version         print version number and exit\n");
    } else {
      snprintf(ERROR_MSG, LINE, "unknown common option %s%s%s!", ErrYellow(),
               option[i], ErrRed());
      PrintError();
      exit(1);
    }
  }
} //}}}
// detect options common for most utilities //{{{
COMMON_OPT CommonOptions(const int argc, char **argv, const SYS_FILES f) {
  COMMON_OPT opt;
  opt.start = -1;
  opt.end = -1;
  opt.skip = 0;
  // -v option - verbose output
  opt.verbose = BoolOption(argc, argv, "--verbose");
  // --silent option - silent mode
  SilentOption(argc, argv, &opt.verbose, &opt.silent);
  // starting/ending timestep
  if (IntegerOption1(argc, argv, "-st", &opt.start)) {
    if (opt.start <= 0) {
      s_strcpy(ERROR_MSG, "positive number required", LINE);
      PrintErrorOption("-st");
      exit(1);
    }
  } else {
    opt.start = 1;
  }
  if (IntegerOption1(argc, argv, "-e", &opt.end) && opt.end <= 0) {
    s_strcpy(ERROR_MSG, "positive number required", LINE);
    PrintErrorOption("-e");
    exit(1);
  }
  ErrorStartEnd(opt.start, opt.end);
  // number of timesteps to skip per one used
  if (IntegerOption1(argc, argv, "-sk", &opt.skip) && opt.skip <= 0) {
    s_strcpy(ERROR_MSG, "positive number required", LINE);
    PrintErrorOption("-sk");
    exit(1);
  }
  opt.skip++; // 'skip' steps are skipped, so every 'skip+1'-th step is used
  if (f.coor.type == LDATA_FILE) {
    opt.start = 1;
    opt.skip = 1;
    opt.end = 1;
  }
  return opt;
} //}}}
// exclude specified molecule names (-x <mol name(s)>) //{{{
bool ExcludeOption(const int argc, char **argv, SYSTEM *System) {
  // set all molecules to use
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    System->MoleculeType[i].Flag = true;
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-x") == 0) {
      // wrong argument to -x option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG, "missing an argument "
                 "(or molecule name beginning with a dash)", LINE);
        PrintErrorOption("-x");
        exit(1);
      } //}}}
      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {
        int type = FindMoleculeName(argv[i+1+j], *System);
        if (type == -1) { // is it in vsf?
          snprintf(ERROR_MSG, LINE, "non-existent molecule %s%s",
                   ErrYellow(), argv[i+1+j]);
          PrintErrorOption("-x");
          return true;
        } else {
          // exclude that molecule
          System->MoleculeType[type].Flag = false;
        }
        j++;
      }
    }
  }
  return false;
} //}}}
// tag bead types with true/false (if missing, set all to opposite) //{{{
bool BeadTypeOption(const int argc, char **argv, const char opt[],
                    const bool use, bool *flag, const SYSTEM System) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int pos = i;
      while (++pos < argc && argv[pos][0] != '-') {
        int type = FindBeadType(argv[pos], System);
        if (type == -1) {
          if (snprintf(ERROR_MSG, LINE, "non-existent bead name %s%s",
                       ErrYellow(), argv[pos]) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption(opt);
          ErrorBeadType(argv[pos], System);
          exit(1);
        }
        flag[type] = use;
      }
      if (pos == (i + 1)) {
        s_strcpy(ERROR_MSG, "at least one bead type is required", LINE);
        PrintErrorOption(opt);
        exit(1);
      }
      return true; // option is present
    }
  }
  return false; // option is not present
} //}}}
// tag molecule types with true/false (if missing, set all to opposite) //{{{
bool MoleculeTypeOption(const int argc, char **argv, const char *opt,
                        const bool use, bool *flag, const SYSTEM System) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int pos = i;
      while (++pos < argc && argv[pos][0] != '-') {
        int type = FindMoleculeName(argv[pos], System);
        if (type == -1) {
          if (snprintf(ERROR_MSG, LINE, "non-existent molecule name %s%s",
                       ErrYellow(), argv[pos]) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption(opt);
          ErrorMoleculeType(argv[pos], System);
          exit(1);
        }
        flag[type] = use;
      }
      if (pos == (i + 1)) {
        s_strcpy(ERROR_MSG, "at least one molecule type is required", LINE);
        PrintErrorOption(opt);
        exit(1);
      }
      return true; // option is present
    }
  }
  return false; // option is not present
} // }}}
// general boolean option //{{{
bool BoolOption(const int argc, char **argv, const char *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      return true;
    }
  }
  return false;
} // }}}
// general option with multiple integer arguments (up to 'max') //{{{
bool IntegerOption(const int argc, char **argv, const int max,
                   const char *opt, int *count, int *values) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0, // number of arguments
          arg = i+1+n; // argument ids in the command
      long val;
      while (arg < argc && IsIntegerNumber(argv[arg], &val)) {
        values[n] = val;
        n++;
        arg = i+1+n;
        if (TooManyArgsWarn(max, n, opt, count)) {
          return true;
        }
      }
      ArgumentMissingErr(n, opt);
      *count = n;
      return true;
    }
  }
  return false;
}
bool IntegerOption1(const int argc, char **argv, const char *opt, int *value) {
  int count = 0;
  if (IntegerOption(argc, argv, 1, opt, &count, value)) {
    ArgumentNumberErr(count, 1, opt);
    return true; // option present
  }
  return false; // option not present
}
bool IntegerOption2(const int argc, char **argv,
                    const char *opt, int value[2]) {
  int count = 0;
  if (IntegerOption(argc, argv, 2, opt, &count, value)) {
    ArgumentNumberErr(count, 2, opt);
    return true; // option present
  }
  return false; // option not present
} //}}}
// general option with multiple double arguments (up to 'max') //{{{
bool DoubleOption(const int argc, char **argv, const int max,
                  const char *opt, int *count, double values[max]) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      int n = 0; // number of arguments
      // read integers
      int arg = i+1+n;
      double val;
      while (arg < argc && IsRealNumber(argv[arg], &val)) {
        values[n] = val;
        n++;
        arg = i+1+n;
        if (TooManyArgsWarn(max, n, opt, count)) {
          return true;
        }
      }
      ArgumentMissingErr(n, opt);
      *count = n;
      return true;
    }
  }
  return false;
}
bool DoubleOption1(const int argc, char **argv,
                   const char *opt, double *value) {
  int count = 0;
  if (DoubleOption(argc, argv, 1, opt, &count, value)) {
    ArgumentNumberErr(count, 1, opt);
    return true; // option present
  }
  return false; // option not present
}
bool DoubleOption2(const int argc, char **argv,
                   const char *opt, double value[2]) {
  int count = 0;
  if (DoubleOption(argc, argv, 2, opt, &count, value)) {
    ArgumentNumberErr(count, 2, opt);
    return true; // option present
  }
  return false; // option not present
}
bool DoubleOption3(const int argc, char **argv,
                   const char *opt, double value[3]) {
  int count = 0;
  if (DoubleOption(argc, argv, 3, opt, &count, value)) {
    if (count != 3) {
      s_strcpy(ERROR_MSG, "three numeric arguments required", LINE);
      PrintErrorOption(opt);
      exit(1);
    } else {
      return true; // option present
    }
  }
  return false; // option not present
} //}}}
// general option with filename and integer(s) arguments //{{{
bool FileIntegerOption(const int argc, char **argv, const int min,
                       const int max, const char *opt, int *values,
                       int *count, char *file) {
  int n = 0;
  *count = 0;
  file[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      MissingFilenameError(argc, argv, opt, i);
      snprintf(file, LINE, "%s", argv[i+1]);
      // read integers
      if (max == 0) {
        return true;
      } else {
        while ((i+2+n) < argc && argv[i+2+n][0] != '-') {
          // Error - non-numeric or missing argument
          long val;
          if (!IsIntegerNumber(argv[i+2+n], &val)) {
            s_strcpy(ERROR_MSG, "arguments must be non-negative whole numbers",
                     LINE);
            PrintErrorOption(opt);
            exit(1);
          }
          values[n] = val;
          n++;
          if (TooManyArgsWarn(max, n, opt, count)) {
            return true;
          }
        }
        if (n < min) {
          s_strcpy(ERROR_MSG, "not enough numeric arguments", LINE);
          PrintErrorOption(opt);
          exit(1);
        }
      }
      *count = n;
      return true; // option present
    }
  }
  return false; // option not present
}
bool FileOption(const int argc, char **argv, const char *opt, char *file) {
  int trash;
  if (FileIntegerOption(argc, argv, 0, 0, opt, &trash, &trash, file)) {
    return true;
  }
  return false;
} //}}}
// general option with filename and double(s) arguments //{{{
bool FileDoubleOption(const int argc, char **argv, const int max,
                      const char *opt, double *values, int *count, char *file) {
  int n = 0;
  *count = 0;
  file[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      MissingFilenameError(argc, argv, opt, i);
      snprintf(file, LINE, "%s", argv[i+1]);
      // read doubles
      if (max == 0) {
        return true;
      } else {
        while ((i+2+n) < argc && argv[i+2+n][0] != '-') {
          // Error - non-numeric or missing argument
          double val;
          if (!IsRealNumber(argv[i+2+n], &val)) {
            s_strcpy(ERROR_MSG, "arguments must be non-negative numbers",
                     LINE);
            PrintErrorOption(opt);
            exit(1);
          }
          values[n] = val;
          n++;
          if (TooManyArgsWarn(max, n, opt, count)) {
            return true;
          }
        }
        ArgumentMissingErr(n, opt);
      }
      *count = n;
      return true; // option present
    }
  }
  return false; // option not present
} //}}}

// STATIC IMPLEMENTATIONS
// option for output verbosity (--silent) //{{{
static void SilentOption(const int argc, char **argv,
                         bool *verbose, bool *silent) {
  *silent = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--silent") == 0) {
      *verbose = false;
      *silent = true;
      break;
    }
  }
} //}}}
// print AnalysisTools version number (--version) //{{{
static bool VersionOption(const int argc, char *argv[]) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--version") == 0) {
      fprintf(stdout, "AnalysisTools by Karel Šindelka (KaGaSi), version %s"
              " (released %s)\n", VERSION, DATE);
      fprintf(stdout, "Download at https://github.com/KaGaSi/AnalysisTools/"
              "releases/tag/v%s\n", VERSION);
      return true;
    }
  }
  return false;
} //}}}
// exit if file name is missing in an option where it's required //{{{
static void MissingFilenameError(const int argc, char **argv,
                                 const char *opt, const int i) {
  // Error - no output file name
  if ((i+1) >= argc || argv[i+1][0] == '-') {
    s_strcpy(ERROR_MSG, "missing file name "
             "(or the file name begins with a dash)", LINE);
    PrintErrorOption(opt);
    exit(1);
  }
} //}}}
// warn if too many arguments //{{{
static bool TooManyArgsWarn(const int max, const int n,
                            const char *opt, int *count) {
  if (n > max) {
    snprintf(ERROR_MSG, LINE, "too many arguments; only the first %d "
             "used", max);
    PrintErrorOption(opt);
    *count = max;
    return true;
  } else {
    return false;
  }
} //}}}
// exit if wrong number of arguments  //{{{
static void ArgumentNumberErr(const int count, const int n, const char *opt) {
  if (count != n) {
    snprintf(ERROR_MSG, LINE, "%d numeric argument(s) required", n);
    PrintErrorOption(opt);
    exit(1);
  }
} //}}}
// exit if missing arguments //{{{
static void ArgumentMissingErr(const int n, const char *opt) {
  if (n == 0) {
    s_strcpy(ERROR_MSG, "missing numeric argument(s)", LINE);
    PrintErrorOption(opt);
    exit(1);
  }
} //}}}
