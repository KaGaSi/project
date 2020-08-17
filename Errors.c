#include "Errors.h"

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
    fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", file);
    fprintf(stderr, " - starting timestep (\033[1;33m%d\033[1;31m) is higher", step);
    fprintf(stderr, " than the total number of steps (\033[1;33m%d\033[1;31m)\n\n", step);
    fprintf(stderr, "\033[0m");
    return true;
  } else {
    ungetc(test,coor);
    return false;
  }
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
  fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m does not have a correct extension (", file);
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
  fprintf(stderr, "\nError: cannot open \033[1;33m%s\033[1;31m for ", file);
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
 * Error when non-numeric argument is present instead of a number
 */
void ErrorNaN(char *option) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "\nError: non-numeric argument for \033[1;33m%s\033[1;31m\n\n", option);
  fprintf(stderr, "\033[0m");
} //}}}

// ErrorOption() //{{{
/**
 * Error when unknown option specified as argument
 */
void ErrorOption(char *option) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "\nError: non-existent option \033[1;33m%s\033[1;31m\n\n", option);
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

// ErrorMoleculeType() //{{{
/**
 * Error when non-existent molecule is used.
 */
void ErrorMoleculeType(COUNTS Counts, MOLECULETYPE *MoleculeType) {
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "       Possible molecule names: %s\n", MoleculeType[0].Name);
  for (int i = 1; i < Counts.TypesOfMolecules; i++) {
    fprintf(stderr, "                            %s\n", MoleculeType[i].Name);
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
  fprintf(stderr, "       Wrong line: |");
  fprintf(stderr, "\033[1;33m");
  for (int i = 0; i < words; i++) {
    if (i != 0) {
      putc(' ', stderr);
    }
    fprintf(stderr, "%s", split[i]);
  }
  fprintf(stderr, "\033[1;31m");
  fprintf(stderr, "|\n\n");
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
    fprintf(stderr, "\nWarning: system in \033[1;36m%s\033[1;33m", file);
    fprintf(stderr, "has net electric charge (q = \033[1;36m%lf\033[1;33m)!\n\n", charge);
    fprintf(stderr, "\033[0m");
  }
} //}}}
