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
  putchar('\n');
} //}}}

// ErrorArgNumber() //{{{
/** Error when insufficient number of arguments
 */
void ErrorArgNumber(int count, int need) {
  fprintf(stderr, "\nError: too few mandatory arguments (%d instead of %d)\n", count, need);
  putchar('\n');
} //}}}

// ErrorDiscard()  //{{{
/** Error when number of starting step is higher then the total number of steps
 * in a coordinate file
 */
bool ErrorDiscard(int start, int step, char *file, FILE *coor) {
  int test;
  if ((test = getc(coor)) == EOF) {
    fflush(stdout);
    fprintf(stderr, "\nError: %s - starting timestep (%d) is higher than the total number of steps (%d)\n\n", file, start, step);
    return true;
  } else {
    ungetc(test,coor);
  }
  return false;
} //}}}

// ErrorExtension() //{{{
/** Error when missing or incorrect file extension
 */
bool ErrorExtension(char *file, int number, char extension[][5]) {

  char *dot = strrchr(file, '.');

  for (int i = 0; i < number; i++) {
    if (dot && strcmp(dot, extension[i]) == 0) {
      return false;
    }
  }

  fprintf(stderr, "\nError: '%s' does not have a correct extension (", file);
  for (int i = 0; i < number; i++) {
    if (i < (number-1)) {
      fprintf(stderr, "'%s', ", extension[i]);
    } else {
      fprintf(stderr, "'%s')\n", extension[i]);
    }
  }
  return true;
} //}}}

// ErrorFileOpen() //{{{
/** Error when open file
 */
void ErrorFileOpen(char *file, char mode) {
  fprintf(stderr, "\nError: cannot open '%s' for ", file);
  switch(mode) {
    case 'r':
      fprintf(stderr, "reading\n");
      break;
    case 'w':
      fprintf(stderr, "writing\n");
      break;
    case 'a':
      fprintf(stderr, "appending\n");
      break;
    default :
      fprintf(stderr, "...well, it seems you found something new to do with a file!\n");
      fprintf(stderr, "Use r(ead), w(rite), or a(ppend).\n\n");
  }
  putchar('\n');
} //}}}

// ErrorNaN() //{{{
/** Error when unknown non-numeric argument is present instead of a number
 */
void ErrorNaN(char *option) {
  fprintf(stderr, "\nError: non-numeric argument for '%s'\n", option);
  putchar('\n');
} //}}}

// ErrorOption() //{{{
/** Error when unknown option specified as argument
 */
void ErrorOption(char *option) {
  fprintf(stderr, "\nError: non-existent option '%s'\n", option);
  putchar('\n');
} //}}}
