#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "AnalysisTools.h"
#include "Options.h"

// CommonOptions() //{{{
/**
 * Function for options common to most of the utilities.
 */
bool CommonOptions(int argc, char **argv, char **vsf_file,char **bonds_file,
                   bool *verbose, bool *verbose2, bool *silent, bool *script) {

//// -i <name> option - filename of input structure file //{{{
//(*vsf_file)[0] = '\0'; // check if -i option is used
//for (int i = 1; i < argc; i++) {
//  if (strcmp(argv[i], "-i") == 0) {

//    // wrong argument to -i option
//    if ((i+1) >= argc || argv[i+1][0] == '-') {
//      fprintf(stderr, "\nMissing argument to '-i' option ");
//      fprintf(stderr, "(or filename beginning with a dash)!\n");

//      return(true);
//    }

//    // check if .vsf ending is present
//    char *vsf = strrchr(argv[i+1], '.');
//    if (!vsf || strcmp(vsf, ".vsf")) {
//      fprintf(stderr, "'-i' arguments does not have .vsf ending!\n");
//      return(true);
//    }

//    strcpy(*vsf_file, argv[i+1]);
//  }
//}

//// -i option is not used
//if ((*vsf_file)[0] == '\0') {
//  strcpy(*vsf_file, "dl_meso.vsf");
//} //}}}

//// -b <name> option - filename of input bond file //{{{
//(*bonds_file)[0] = '\0'; // check if -b option is used
//for (int i = 1; i < argc; i++) {
//  if (strcmp(argv[i], "-b") == 0) {

//    // wrong argument to -b option
//    if ((i+1) >= argc || argv[i+1][0] == '-') {
//      fprintf(stderr, "\nMissing argument to '-b' option ");
//      fprintf(stderr, "(or filename beginning with a dash)!\n\n");

//      return(true);
//    }

//    strcpy(*bonds_file, argv[i+1]);
//  }
//} //}}}

//// -v option - verbose output //{{{
//*verbose = false;
//for (int i = 1; i < argc; i++) {
//  if (strcmp(argv[i], "-v") == 0) {
//    *verbose = true;

//    break;
//  }
//} //}}}

//// -V option - verbose output with comments from input .vcf file //{{{
//*verbose2 = false;
//for (int i = 1; i < argc; i++) {
//  if (strcmp(argv[i], "-V") == 0) {
//    *verbose = true;
//    *verbose2 = true;

//    break;
//  }
//} //}}}

//// -s option - silent mode //{{{
//*silent = false;
//for (int i = 1; i < argc; i++) {
//  if (strcmp(argv[i], "-s") == 0) {
//    *verbose = false;
//    *verbose2 = false;
//    *silent = true;

//    break;
//  }
//} //}}}

//// --script  option - meant for when output is routed to file, so don't use flush & \r //{{{
//*script = false;
//for (int i = 1; i < argc; i++) {
//  if (strcmp(argv[i], "--script") == 0) {
//    *script = true;

//    break;
//  }
//} //}}}

  return(false);
} //}}}

// VsfFileOption() //{{{
bool VsfFileOption(int argc, char **argv, char **vsf_file) {

  (*vsf_file)[0] = '\0'; // check if -i option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");

        return(true);
      } //}}}

      // check if .vsf ending is present (required by VMD) //{{{
      char *vsf = strrchr(argv[i+1], '.');
      if (!vsf || strcmp(vsf, ".vsf")) {
        fprintf(stderr, "'-i' argument does not have .vsf ending!\n");

        return (true);
      } //}}}

      strcpy(*vsf_file, argv[i+1]);
    }
  }

  // -i option is not used{{{
  if ((*vsf_file)[0] == '\0') {
    strcpy(*vsf_file, "dl_meso.vsf");
  } //}}}

  return (false);
} //}}}

// BondsFileOption() //{{{
bool BondsFileOption(int argc, char **argv, char **bonds_file) {
  (*bonds_file)[0] = '\0'; // check if -b option is used

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-b") == 0) {

      // wrong argument to -b option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-b' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n\n");

        return(true);
      }

      strcpy(*bonds_file, argv[i+1]);
    }
  }

  return(false);
} //}}}

// VerboseShortOption() //{{{
bool VerboseShortOption(int argc, char **argv, bool *verbose) {
  *verbose = false;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      *verbose = true;

      break;
    }
  }
} //}}}

// VerboseLongOption() //{{{
bool VerboseLongOption(int argc, char **argv, bool *verbose, bool *verbose2) {

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
bool SilentOption(int argc, char **argv, bool *verbose, bool *verbose2,
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

// ScriptOption() //{{{
bool ScriptOption(int argc, char **argv, bool *script) {
  *script = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--script") == 0) {
      *script = true;

      break;
    }
  }
} // }}}

// ExcludeOption() //{{{
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

        int type = FindMoleculeType(argv[i+1+j], Counts, (*MoleculeType));

        if (type == -1) { // is it in FIELD?
          fprintf(stderr, "Non-existent molecule name (%s)!\n", argv[i+1+j]);
          fprintf(stderr, "Molecule names in FIELD:\n");
          for (int k = 0; k < Counts.TypesOfMolecules; k++) {
            fprintf(stderr, "%3d %s\n", k, (*MoleculeType)[k].Name);
          }

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

// JoinedCoorOption() //{{{
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
      if (!dot || strcmp(dot, ".vcf")) {
        fprintf(stderr, "<joined.vcf> '%s' does not have .vcf ending!\n", joined_vcf);

        return(true);
      } //}}}
    }
  }

  return(false);
} //}}}
