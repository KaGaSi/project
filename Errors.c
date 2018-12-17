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
  fprintf(stderr, "Possibilities: error in '%s', e.g., wrong number of lines per timestep or columns per line, \
incorrect blank lines\n", input_vcf);
  fprintf(stderr, "               error in '%s', e.g., wrong number of beads which can cause a complete chaos in internal bead numbering\n",
 input_vsf);
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
void ErrorExtension(char *file, char *extension) {
    fprintf(stderr, "Error: '%s' does not have '%s' extension\n", file, extension);
  putchar('\n');
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
