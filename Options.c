#include "Options.h"

// CommonHelp() //{{{
/**
 * Function to print help for common options, either for `-h` help option
 * or program error.
 */
void CommonHelp(bool error) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
  }

  fprintf(ptr, "   <standard options>\n");
  fprintf(ptr, "      -i <name>      use input .vsf file different from traject.vsf\n");
//fprintf(ptr, "      -b <name>      file containing bond alternatives to FIELD\n");
  fprintf(ptr, "      -v             verbose output\n");
  fprintf(ptr, "      --silent       no output (overrides verbose option)\n");
  fprintf(ptr, "      -h             print this help and exit\n");
  fprintf(ptr, "      --script       no progress output (useful if output is routed to a file)\n");
  fprintf(ptr, "      --version      print version number and exit\n");
} //}}}

// CommonOptions() //{{{
/**
 * Function for options common to most of the utilities.
 */
// TODO: copy vsf_file to old_name; then if vsf_file[0] == '\0', copy old_name back to vsf_file
void CommonOptions(int argc, char **argv, char **vsf_file,
                   bool *verbose, bool *silent, bool *script) {

  // -i <name> option - filename of input structure file //{{{
  // save provided file name
  char old_name[LINE];
  strcpy(old_name, *vsf_file);
  // test if '-i' option is there
  if (FileOption(argc, argv, "-i", vsf_file)) {
    exit(1);
  }
  // copy back the old name if '-i' option is not present
  if (*vsf_file[0] == '\0') {
    strcpy(*vsf_file, old_name);
  }
  // test if structure file ends with '.vsf' or '.vtf'
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(*vsf_file, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}
  // -v option - verbose output
  *verbose = BoolOption(argc, argv, "-v");
  // --silent option - silent mode
  SilentOption(argc, argv, verbose, silent);
  // --script - meant for when output is routed to file, so don't use flush & \r
  *script = BoolOption(argc, argv, "--script"); // do not use \r & co.
} //}}}

// SilentOption() //{{{
/**
 * Option to not print anything to stdout (or at least no system
 * definitions and no Step: #). Overrides verbose option. Argument: `--silent`
 */
void SilentOption(int argc, char **argv, bool *verbose, bool *silent) {

  *silent = false;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--silent") == 0) {
      *verbose = false;
      *silent = true;

      break;
    }
  }
} //}}}

// VersionOption() //{{{
/**
 * Option to print version number of the program suite.
 */
bool VersionOption(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--version") == 0) {
      fprintf(stdout, "AnalysisTools by Karel Šindelka (KaGaSi), version %s", VERSION);
      fprintf(stdout, " (released %s)\n", DATE);
      fprintf(stdout, "Download at https://github.com/KaGaSi/AnalysisTools/releases/tag/v%s\n", VERSION);
      return true;
    }
  }

  return false;
} //}}}

// ExcludeOption() //{{{
/**
 * Option to exclude specified molecule types from calculations. Gives
 * specified molecule types `Use = false` and the rest `Use = true`.
 * Arguments: `-x <name(s)>`
 */
bool ExcludeOption(int argc, char **argv, COUNTS Counts,
                   MOLECULETYPE **MoleculeType) {

  // set all molecules to use //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    (*MoleculeType)[i].Use = true;
  } //}}}

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-x") == 0) {

      // wrong argument to -x option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing first argument to '-x' option ");
        fprintf(stderr, "(or molecule name beginning with a dash)!\n");
        exit(1);
      } //}}}

      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {

        int type = FindMoleculeType(argv[i+1+j], Counts, *MoleculeType);

        if (type == -1) { // is it in vsf?
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "Error: non-existent molecule name (%s)\n", argv[i+1+j]);
          fprintf(stderr, "Molecule names in structure file:\n");
          for (int k = 0; k < Counts.TypesOfMolecules; k++) {
            fprintf(stderr, "%3d %s\n", k, (*MoleculeType)[k].Name);
          }
          putc('\n', stderr);
          fprintf(stderr, "\033[0m");
          return(true);
        } else {
          // exclude that molecule
          (*MoleculeType)[type].Use = false;
        }

        j++;
      }
    }
  }

  return(false);
} //}}}

// JoinCoorOption() //{{{
/**
 * Option whether to join aggregates and save joined coordinates into a
 * specified file. Arguments: `-j <joined.vcf>`
 */
bool JoinCoorOption(int argc, char **argv, char *joined_vcf) {

  joined_vcf[0] = '\0'; // no -j option
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {

      // wrong argument to -j option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument to '-j' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");

        return(true);
      } //}}}

      strcpy(joined_vcf, argv[i+1]);

      // test if <joined.vcf> filename ends with '.vcf' (required by VMD) //{{{
      char *dot = strrchr(joined_vcf, '.');
      if (!dot || (strcmp(dot, ".vcf") && strcmp(dot, ".vtf"))) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: '%s' does not have .vcf ending!\n", joined_vcf);
        fprintf(stderr, "\033[0m");
        return(true);
      } //}}}
    }
  }

  return(false);
} //}}}

// BeadTypeOption() //{{{
/**
 * Option to choose which bead types to use for calculation. If the option
 * is absent, all bead types are switched to `Use = true`. Argument: `-bt
 * <name(s)>`
 */
bool BeadTypeOption(int argc, char **argv, char *opt, bool use,
                    COUNTS Counts, BEADTYPE **BeadType) {

  // specify what bead types to use - either specified by '-bt' option or all
  int types = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      types = i; // positon of the '-bt' argument in command
      // <type names> - names of bead types to save
      while (++types < argc && argv[types][0] != '-') {
        int type = FindBeadType(argv[types], Counts, *BeadType);
        if (type == -1) {
          fprintf(stderr, "Bead type '%s' does not exist in structure file ('%s' option)!\n", argv[types], opt);
          return(true);
        }

        (*BeadType)[type].Use = true;
      }
    }
  }
  if (types == -1) {
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      (*BeadType)[i].Use = use;
    }
  }

  return(false);
} // }}}

// BoolOption() //{{{
/**
 * Function for any boolean option (i.e. without argument). The option
 * (e.g. `--script`) is an argument of this function.
 */
bool BoolOption(int argc, char **argv, char *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {
      return(true);
    }
  }
  return(false);
} // }}}

// IntegerOption() //{{{
/**
 * Function for any option with integer argument. The option (e.g. `-n`)
 * is an argument of this function.
 */
bool IntegerOption(int argc, char **argv, char *opt, int *value) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // Error - non-numeric or missing argument //{{{
      if ((i+1) >= argc) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Missing numeric argument for '%s' option!\n", opt);
        fprintf(stderr, "\033[0m");
        return(true);
      }
      if ((i+1) >= argc || argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Non-numeric argement for '%s' option!\n", opt);
        fprintf(stderr, "\033[0m");
        return(true);
      } //}}}

      *value = atoi(argv[i+1]);
    }
  }

  return(false);
} //}}}

// DoubleOption() //{{{
/**
 * Function for any option with double argument. The option (e.g. `-n`)
 * is an argument of this function.
 */
bool DoubleOption(int argc, char **argv, char *opt, double *value) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // Error - non-numeric or missing argument //{{{
      if ((i+1) >= argc) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Missing numeric argument for '%s' option!\n", opt);
        fprintf(stderr, "\033[0m");
        return(true);
      }
      if ((i+1) >= argc || argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Non-numeric argement for '%s' option!\n", opt);
        fprintf(stderr, "\033[0m");
        return(true);
      } //}}}

      *value = atof(argv[i+1]);
    }
  }

  return(false);
} //}}}

// MultiIntegerOption() //{{{
/**
 * Function for any option with two integer arguments. The option (e.g. `-n`)
 * is an argument of this function.
 */
bool MultiIntegerOption(int argc, char **argv, char *opt, int *count, int *values) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      int n = 0; // number of arguments

      // read integers
      int arg = i+1+n;
      while ((arg) < argc && argv[arg][0] != '-') {

        // Error - non-numeric or missing argument //{{{
        if (argv[arg][0] < '0' || argv[i+1+n][0] > '9') {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "Error: option '%s' - non-numeric argement(s)\n", opt);
          fprintf(stderr, "\033[0m");
          return true;
        } //}}}

        values[n] = atoi(argv[arg]);

        n++;
        arg = i+1+n;

        // warning - too many numeric arguments //{{{
        if (n == 100) {
          fprintf(stderr, "Warning: Option '%s' - too many numberic arguments; ", opt);
          fprintf(stderr, "only first 100 used\n");
          *count = n;
          return true;
        } //}}}
      }

      *count = n;
    }
  }

  return false;
} //}}}

// MultiDoubleOption() //{{{
/**
 * Function for any option with two double arguments. The option (e.g. `-n`)
 * is an argument of this function.
 */
bool MultiDoubleOption(int argc, char **argv, char *opt, int *count, double *values) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      int n = 0; // number of arguments

      // read integers
      int arg = i+1+n;
      while ((arg) < argc && argv[arg][0] != '-') {

        // Error - non-numeric or missing argument //{{{
        if (argv[arg][0] < '0' || argv[i+1+n][0] > '9') {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "Error: option '%s' - non-numeric argement(s)\n", opt);
          fprintf(stderr, "\033[0m");
          return true;
        } //}}}

        values[n] = atof(argv[arg]);

        n++;
        arg = i+1+n;

        // warning - too many numeric arguments //{{{
        if (n == 100) {
          fprintf(stderr, "Warning: Option '%s' - too many numberic arguments; ", opt);
          fprintf(stderr, "only first 100 used\n");
          *count = n;
          return true;
        } //}}}
      }

      *count = n;
    }
  }

  return false;
} //}}}

// FileIntsOptions() //{{{
/**
 * Function for any option with a file name and up to 100 integer
 * arguments. The option is an argument of this function.
 */
bool FileIntsOption(int argc, char **argv, char *opt, int *values, int *count, char *file) {

  int n = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // Error - no output file name //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-' ||
          (argv[i+1][0] >= '0' && argv[i+1][0] <= '9')) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Error: option '%s' - missing output file name", opt);
        fprintf(stderr, " (or the file name begins with '-' or a number)\n");
        fprintf(stderr, "\033[0m");
        return false;
      } //}}}
      strcpy(file, argv[i+1+n]);

      // read integers
      while ((i+2+n) < argc && argv[i+2+n][0] != '-') {

        // Error - non-numeric or missing argument //{{{
        if ((i+2+n) >= argc || argv[i+2+n][0] < '0' || argv[i+2+n][0] > '9') {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "Error: option '%s' - non-numeric argement(s)\n", opt);
          fprintf(stderr, "\033[0m");
          return true;
        } //}}}

        values[n] = atoi(argv[i+2+n]);

        n++;
        // warning - too many numeric arguments //{{{
        if (n == 100) {
          fprintf(stderr, "Warning: Option '%s' - too many numberic arguments; ", opt);
          fprintf(stderr, "only first 100 used\n");
          *count = n;
          return true;
        } //}}}
      }
    }
  }

  *count = n;

  return false;
} //}}}

// FileOption() //{{{
/**
 * Generic option for file name. The option is an argument of this function.
 */
bool FileOption(int argc, char **argv, char *opt, char **name) {

  (*name)[0] = '\0';

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // wrong argument to the option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '%s' option ", opt);
        fprintf(stderr, "(or filename beginning with a dash)!\n\n");

        return(true);
      }

      strcpy(*name, argv[i+1]);
    }
  }

  return(false);
} //}}}

// MoleculeTypeOption() //{{{
/**
 * Generic option for molecule type that can take one
 * argument. The option is an argument of this function.
 */
bool MoleculeTypeOption(int argc, char **argv, char *opt, int *moltype,
                        COUNTS Counts, MOLECULETYPE **MoleculeType) {

  *moltype = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Missing argument for '%s' option (or molecule name beginning with a dash)!\n", opt);
        fprintf(stderr, "\033[0m");
        return(true);
      } //}}}

      *moltype = FindMoleculeType(argv[i+1], Counts, *MoleculeType);
      if (*moltype == -1) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Error: molecule '%s' does not exist in structure file ('%s' option)!\n\n", argv[i+1], opt);
        fprintf(stderr, "   Present molecule types:\n");
        for (int j = 0; j < Counts.TypesOfMolecules; j++) {
          fprintf(stderr, "%s\n", (*MoleculeType)[j].Name);
        }
        fprintf(stderr, "\033[0m");
        return(true);
      }
    }
  }

  return(false);
} //}}}

// MoleculeTypeOption2() //{{{
/**
 * Generic option for molecule types. The option is an argument of this
 * function.
 */
bool MoleculeTypeOption2(int argc, char **argv, char *opt, int **moltype,
                         COUNTS Counts, MOLECULETYPE **MoleculeType) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // set all moltypes not to be used
      for (int i = 0; i < Counts.TypesOfMolecules; i++) {
        (*moltype)[i] = 0;
      }

      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "Missing argument for '%s' option (or molecule name beginning with a dash)!\n", opt);
        fprintf(stderr, "\033[0m");
        return(true);
      } //}}}

      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {

        int type = FindMoleculeType(argv[i+1+j], Counts, *MoleculeType);
        if (type == -1) { // is argv[i+1+j] in vsf?
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "Error: molecule '%s' does not exist in the provided coordinate file ('%s' option)!\n\n", argv[i+1+j], opt);
          fprintf(stderr, "   Present molecule types:\n");
          for (int k = 0; k < Counts.TypesOfMolecules; k++) {
            fprintf(stderr, "%s\n", (*MoleculeType)[k].Name);
          }
          fprintf(stderr, "\033[0m");
          return(true);
        }
        (*moltype)[type] = 1;

        j++;
      }
    }
  }

  return(false);
} //}}}
