#include "Errors.h"

// ErrorCoorRead() //{{{
/**
 * Error when reading vcf file
 */
void ErrorCoorRead(char *input_vcf, int bead, int step, char *stuff, char *input_vsf) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "\nError: %s - cannot read coordinates (bead %d; step %d - '%s')\n\n", input_vcf, bead, step, stuff);
  fprintf(stderr, "\033[0m");
  putchar('\n');
} //}}}

// ErrorArgNumber() //{{{
/**
 * Error when insufficient number of arguments
 */
void ErrorArgNumber(int count, int need) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "\nError: too few mandatory arguments (%d instead of %d)\n\n", count, need);
  fprintf(stderr, "\033[0m");
} //}}}

// ErrorDiscard()  //{{{
/**
 * Error when number of starting step is higher then the total number of steps
 * in a coordinate file
 */
bool ErrorDiscard(int start, int step, char *file, FILE *coor) {
  int test;
  if ((test = getc(coor)) == EOF) {
    fflush(stdout);
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: %s - starting timestep (%d) is higher than the total number of steps (%d)\n\n", file, start, step);
    return true;
  } else {
    ungetc(test,coor);
  }
  fprintf(stderr, "\033[0m");
  return false;
} //}}}

// ErrorExtension() //{{{
/**
 * Error when missing or incorrect file extension
 */
bool ErrorExtension(char *file, int number, char extension[][5]) {

  char *dot = strrchr(file, '.');

  for (int i = 0; i < number; i++) {
    if (dot && strcmp(dot, extension[i]) == 0) {
      return false;
    }
  }

  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "\nError: '%s' does not have a correct extension (", file);
  for (int i = 0; i < number; i++) {
    if (i < (number-1)) {
      fprintf(stderr, "'%s', ", extension[i]);
    } else {
      fprintf(stderr, "'%s')\n\n", extension[i]);
    }
  }
  fprintf(stderr, "\033[0m");
  return true;
} //}}}

// ErrorFileOpen() //{{{
/**
 * Error when open file
 */
void ErrorFileOpen(char *file, char mode) {
  fprintf(stderr, "\033[1;31m");
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
  fprintf(stderr, "\033[0m");
} //}}}

// ErrorNaN() //{{{
/**
 * Error when unknown non-numeric argument is present instead of a number
 */
void ErrorNaN(char *option) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "\nError: non-numeric argument for '%s'\n\n", option);
  fprintf(stderr, "\033[0m");
} //}}}

// ErrorOption() //{{{
/**
 * Error when unknown option specified as argument
 */
void ErrorOption(char *option) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "\nError: non-existent option '%s'\n\n", option);
  fprintf(stderr, "\033[0m");
} //}}}

// ErrorBeadType() //{{{
/**
 * Error when non-existent bead is used.
 */
void ErrorBeadType(COUNTS Counts, BEADTYPE *BeadType) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "       Possible bead names: %s\n", BeadType[0].Name);
  for (int i = 1; i < Counts.TypesOfBeads; i++) {
    fprintf(stderr, "                            %s\n", BeadType[i].Name);
  }
  putc('\n', stderr);
  fprintf(stderr, "\033[0m");
} //}}}

// ErrorPrintLine() //{{{
/**
 * Print provided strings (array of strings generally created using
 * SplitLine()) to error output.
 */
void ErrorPrintLine(char split[30][100], int words) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "        Wrong line:|");
  fprintf(stderr, "\033[0m");
  for (int i = 0; i < words; i++) {
    fprintf(stderr, " %s", split[i]);
  }
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, " |\n\n");
  fprintf(stderr, "\033[0m");
} //}}}

// WarnElNeutrality() //{{{
/**
 * Function to warn if the system is not electrically neutral.
 */
void WarnElNeutrality(COUNTS Counts, BEADTYPE *BeadType, char *file) {
  double charge = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    charge += BeadType[i].Charge * BeadType[i].Number;
  }
  if (charge != 0) {
    fprintf(stderr, "\033[1;33m");
    fprintf(stderr, "\nWarning: system in %s has net electric charge (q = %lf)!\n\n", file, charge);
    fprintf(stderr, "\033[0m");
  }
} //}}}
