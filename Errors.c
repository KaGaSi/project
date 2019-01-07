#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "Errors.h"

// ErrorCoorRead() //{{{
/** Error when reading vcf file
 */
void ErrorCoorRead(char *input_vcf, int bead, int step, char *stuff, char *input_vsf) {
  fprintf(stderr, "\nError: %s - cannot read coordinates (bead %d; step %d - '%s')\n\n", input_vcf, bead, step, stuff);
  fprintf(stderr, "Possibilities\n");
  fprintf(stderr, "  1) error in '%s', e.g., wrong number of lines per timestep or columns per line, \
incorrect blank lines\n", input_vcf);
  fprintf(stderr, "  2) error in '%s', e.g., wrong number of beads which can cause a complete chaos in \
internal bead numbering\n", input_vsf);
  fprintf(stderr, "  3) used '-x' option in generating %s by SelectedVcf, but reading of systems \
configuration works only if all beads of a given type are in a vcf file \n", input_vcf);
  putchar('\n');
} //}}}

// ErrorArgNumber() //{{{
/** Error when insufficient number of arguments
 */
void ErrorArgNumber(int count, int need) {
  fprintf(stderr, "\nError: too few mandatory arguments (%d instead of %d)\n", count, need);
  putchar('\n');
} //}}}

// ErrorOption() //{{{
/** Error when unknown option specified as argument
 */
void ErrorOption(char *option) {
  fprintf(stderr, "\nError: non-existent option '%s'\n", option);
  putchar('\n');
} //}}}

// ErrorNaN() //{{{
/** Error when unknown non-numeric argument is present instead of a number
 */
void ErrorNaN(char *option) {
  fprintf(stderr, "\nError: non-numeric argument for '%s'\n", option);
  putchar('\n');
} //}}}

// ErrorExtension() //{{{
/** Error when missing or incorrect file extension
 */
bool ErrorExtension(char *file, int number, char **extension) {

  char *dot = strrchr(file, '.');

  for (int i = 0; i < number; i++) {
    if (dot && strcmp(dot, extension[i]) == 0) {
      return true;
    }
  }

  fprintf(stderr, "Error: '%s' does not have a correct extension (", file);
  for (int i = 0; i < number; i++) {
    if (i < (number-1)) {
      fprintf(stderr, "'%s', ", extension[i]);
    } else {
      fprintf(stderr, "'%s')\n", extension[i]);
    }
  }
  return false;
} //}}}

// ErrorFileOpen() //{{{
/** Error when open file
 */
void ErrorFileOpen(char *file, char mode) {
  fprintf(stderr, "\nError: cannot open '%s' for ", file);
  if (mode == 'r') {
    fprintf(stderr, "reading\n");
  } else if (mode == 'w') {
    fprintf(stderr, "writing\n");
  } else if (mode == 'a') {
    fprintf(stderr, "appending\n");
  } else {
    fprintf(stderr, "...well it seems you found a completely new thing to do with a file!\n");
  }
  putchar('\n');
} //}}}
