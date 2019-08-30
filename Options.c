#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "Options.h"
#include "AnalysisTools.h"
#include "Errors.h"

// VerboseLongOption() //{{{
/**
 * Option whether to print detailed data to stdout. Data are printed via
 * VerboseOutput() function (and possibly some in-program code). Argument:
 * `-V`
 */
void VerboseLongOption(int argc, char **argv, bool *verbose, bool *verbose2) {

  *verbose2 = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-V") == 0) {
      *verbose = true;
      *verbose2 = true;

      break;
    }
  }
} //}}}

// SilentOption() //{{{
/**
 * Option to not print anything to stdout (or at least no system
 * definitions and no Step: #). Overrides VerboseShortOption and
 * VerboseLongOption. Argument: `-s`
 */
void SilentOption(int argc, char **argv, bool *verbose, bool *verbose2,
                  bool *silent) {

  *silent = false;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-s") == 0) {
      *verbose = false;
      *verbose2 = false;
      *silent = true;

      break;
    }
  }
} //}}}

// ExcludeOption() //{{{
/**
 * Option to exclude specified molecule types from calculations. Gives
 * specified molecule types `Use = false` and the rest `Use = true`.
 * Arguments: `-x <name(s)>`
 */
bool ExcludeOption(int argc, char **argv, Counts Counts,
                   MoleculeType **MoleculeType) {

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
          fprintf(stderr, "Error: non-existent molecule name (%s)\n", argv[i+1+j]);
          fprintf(stderr, "Molecule names in structure file:\n");
          for (int k = 0; k < Counts.TypesOfMolecules; k++) {
            fprintf(stderr, "%3d %s\n", k, (*MoleculeType)[k].Name);
          }
          putc('\n', stderr);
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
        fprintf(stderr, "\nError: '%s' does not have .vcf ending!\n", joined_vcf);

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
bool BeadTypeOption(int argc, char **argv, bool use,
                    Counts Counts, BeadType **BeadType) {

  // specify what bead types to use - either specified by '-bt' option or all
  int types = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-bt") == 0) {
      types = i; // positon of the '-bt' argument in command
      // <type names> - names of bead types to save
      while (++types < argc && argv[types][0] != '-') {
        int type = FindBeadType(argv[types], Counts, *BeadType);
        if (type == -1) {
          fprintf(stderr, "Bead type '%s' does not exist in structure file ('-bt' option)!\n", argv[types]);
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
        fprintf(stderr, "Missing numeric argument for '%s' option!\n", opt);
        return(true);
      }
      if ((i+1) >= argc || argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '%s' option!\n", opt);
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
        fprintf(stderr, "Missing numeric argument for '%s' option!\n", opt);
        return(true);
      }
      if ((i+1) >= argc || argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '%s' option!\n", opt);
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
          fprintf(stderr, "Error: option '%s' - non-numeric argement(s)\n", opt);
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
          fprintf(stderr, "Error: option '%s' - non-numeric argement(s)\n", opt);
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
        fprintf(stderr, "Error: option '%s' - missing output file name", opt);
        fprintf(stderr, " (or the file name begins with '-' or a number)\n");
        return false;
      } //}}}
      strcpy(file, argv[i+1+n]);

      // read integers
      while ((i+2+n) < argc && argv[i+2+n][0] != '-') {

        // Error - non-numeric or missing argument //{{{
        if ((i+2+n) >= argc || argv[i+2+n][0] < '0' || argv[i+2+n][0] > '9') {
          fprintf(stderr, "Error: option '%s' - non-numeric argement(s)\n", opt);
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
bool MoleculeTypeOption(int argc, char **argv, char *opt, int *moltype, Counts
    Counts, MoleculeType **MoleculeType) {

  *moltype = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument for '%s' option (or molecule name beginning with a dash)!\n", opt);
        return(true);
      } //}}}

      *moltype = FindMoleculeType(argv[i+1], Counts, *MoleculeType);
      if (*moltype == -1) {
        fprintf(stderr, "Error: molecule '%s' does not exist in structure file ('%s' option)!\n\n", argv[i+1], opt);
        fprintf(stderr, "   Present molecule types:\n");
        for (int j = 0; j < Counts.TypesOfMolecules; j++) {
          fprintf(stderr, "%s\n", (*MoleculeType)[j].Name);
        }
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
bool MoleculeTypeOption2(int argc, char **argv, char *opt, int **moltype, Counts
    Counts, MoleculeType **MoleculeType) {

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], opt) == 0) {

      // set all moltypes not to be used
      for (int i = 0; i < Counts.TypesOfMolecules; i++) {
        (*moltype)[i] = 0;
      }

      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument for '%s' option (or molecule name beginning with a dash)!\n", opt);
        return(true);
      } //}}}

      // read molecule(s) names
      int j = 0;
      while ((i+1+j) < argc && argv[i+1+j][0] != '-') {

        int type = FindMoleculeType(argv[i+1+j], Counts, *MoleculeType);
        if (type == -1) { // is argv[i+1+j] in vsf?
          fprintf(stderr, "Molecule '%s' does not exist in the provided coordinate file ('%s' option)!\n\n", argv[i+1+j], opt);
          fprintf(stderr, "   Present molecule types:\n");
          for (int k = 0; k < Counts.TypesOfMolecules; k++) {
            (*MoleculeType)[k].Name[0] = 'X';
            fprintf(stderr, "%s\n", (*MoleculeType)[k].Name);
          }
          return(true);
        }
        (*moltype)[type] = 1;

        j++;
      }
    }
  }

  return(false);
} //}}}
