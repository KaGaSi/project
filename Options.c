#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "AnalysisTools.h"
#include "Options.h"

// ExcludeOption() //{{{
bool ExcludeOption(int argc, char **argv, Counts Counts,
                   MoleculeType **MoleculeType) {

  bool error = false;

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
          error = true;
          break;
        } else {
          // exclude that molecule
          (*MoleculeType)[type].Use = false;
        }

        j++;
      }
    }

    if (error) {
      break;
    }
  }

  return (false);
} //}}}

///*
// JoinedCoorOption() //{{{
bool JoinCoorOption(int argc, char **argv, char *joined_vcf) {

  bool error = false;

  joined_vcf[0] = '\0'; // no -j option
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {

      // wrong argument to -j option //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument to '-j' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");

        error = true;
        break;
      } //}}}

      strcpy(joined_vcf, argv[i+1]);

      // test if <joined.vcf> filename ends with '.vcf' (required by VMD) //{{{
      char *dot = strrchr(joined_vcf, '.');
      if (!dot || strcmp(dot, ".vcf")) {
        fprintf(stderr, "<joined.vcf> '%s' does not have .vcf ending!\n", joined_vcf);
        ErrorHelp(argv[0]);

        error = true;
        break;
      } //}}}
    }
  }

  return (false);
} //}}}
//*/
