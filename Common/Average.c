#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input> <column> <discard> <n_blocks>\n\n", cmd);

  fprintf(stderr, "   <input>           input filename\n");
  fprintf(stderr, "   <column>          column number in the file to analyse\n");
  fprintf(stderr, "   <discard>         number of rows discard (from the file beginning)\n");
  fprintf(stderr, "   <n_blocks>        number of blocks for binning\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -h             print this help and exit\n");
} //}}}

int main ( int argc, char** argv ) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
Average uses binning method to analyse data stored in a supplied file. It \
prints average, statistical error and estimate of integrated \
autocorrelation time (tau). Empty lines and lines beginning with '#' are \
skipped. The program prints to standart output 4 values: <n_blocks> \
<simple average> <statistical error> <estimate of tau>.\n\n");

      fprintf(stdout, "\
A way to estimate \
the 'real' value of tau is to use a wide range of <n_blocks> and then plot \
<tau> as a function of <n_blocks>. Since the number of data points in a \
block has to be larger than tau (say, 10 times larger), plotting \
<number of data lines>/10/<n_blocks> vs. <n_blocks> will produce an \
exponential function that will intersect the plotted <tau>. A value of tau \
near the intersection (but to the left where the exponential is above <tau>) \
can be considered a good estimate for tau.\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input> <column> <discard> <n_blocks>\n\n", argv[0]);

      fprintf(stdout, "   <input>           input filename\n");
      fprintf(stdout, "   <column>          column number in the file to analyse\n");
      fprintf(stdout, "   <discard>         number of rows discard (from the file beginning)\n");
      fprintf(stdout, "   <n_blocks>        number of blocks for binning\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -h             print this help and exit\n");
      exit(0);
    }
  }

  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-h") != 0 ) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - filename of input file //{{{
  char input[32];
  strcpy(input, argv[++count]); //}}}

  // <column> - column number to analyze //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<column>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int column = atoi(argv[count]); //}}}

  // <discard> - number of lines to discard from the file beginning //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<discard>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int discard = atoi(argv[count]); //}}}

  // <n_blocks> - number of blocks for binning //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<n_blocks>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int n_blocks = atoi(argv[count]); //}}}

  double *data = calloc(1,sizeof(double));

  // read data from <input> file //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  }

//fprintf(stdout, "input=\"%s\", ", input);
//fprintf(stdout, "column=%d, ", column);
//fprintf(stdout, "discard=%d, ", discard);
//fprintf(stdout, "n_blocks=%d\n", n_blocks);

  int test, lines = 0, all_lines = 0;
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);

    // count data lines, comment, and blank lines in case of error
    all_lines++;

    // get whole line - max 1000 chars //{{{
    char line[1024];
    fgets(line, 1024, fr); //}}}

    // trim trailing whitespace in line //{{{
    int length = strlen(line);
    // last string character is '\n' (at [length-1]), so check the previous one(s) - if there are any
    while (length > 1 &&
           (line[length-2] == ' ' ||
           line[length-2] == '\t')) {
      line[length-2] = line[length-1]; // move newline char
      line[length-1] = '\0'; // add string ending char
      length--;
    } //}}}

    // if not empty line or comment continue //{{{
    char *split = strtok(line," \t");
    if (split[0] != '#' &&
        split[0] != '\n') {

      // number of data lines
      lines++;

      // read correct column //{{{
      int col = 1; // number of columns in the line - first is already loaded in split
      for (col = 1; col < column; col++) {
        // load another 'split'
        split = strtok(NULL, " \t");

        // error - insufficient number of columns //{{{
        if (split == NULL) {
          fprintf(stderr, "Only %d columns in %s on line %d (data line %d)!\n", col, input, all_lines, lines);
          exit(1);
        } //}}}
      } //}}}

      // save the value //{{{
      if (discard < lines) {
        int count = lines - discard - 1;
        data = realloc(data, (count+1)*sizeof(double));
        data[count] = atof(split);
      } //}}}
    } //}}}
  }
  fclose(fr); //}}}

  // error - <discard> is too large //{{{
  if (discard >= lines) {
    fprintf(stderr, "Error: <discard> parameter is too large (%d) - ", discard);
    fprintf(stderr, "there are only %d data lines in %s\n\n", lines, input);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // variables //{{{
  // number of data points must be dividable by 'n_blocks'
  int remainder = (lines - discard) % n_blocks;
  // total number of data points to consider
  count = lines - discard - remainder;
  // number of data points per block
  int data_per_block = count / n_blocks;
  // allocate array for block averages
  long double *avg_block = calloc(n_blocks,sizeof(long double));
  // overall averages
  long double avg_all[2] = {0}; // [0] for avg, [1] for avg of squares //}}}

  // calculate averages //{{{
  int k = remainder; // first datapoint to consider
  for (int i = 0; i < n_blocks; i++) {
    for (int j = 0; j < data_per_block; j++) {
      avg_all[0] += data[k];
      avg_all[1] += SQR(data[k]);
      avg_block[i] += data[k++];
    }
  }

  avg_all[0] /= count;
  avg_all[1] /= count;
  for (int i = 0; i < n_blocks; i++) {
    avg_block[i] /= data_per_block;
  } //}}}

  // standard deviation for block averages
  double block_stdev = 0;
  for (int i = 0; i < n_blocks; i++) {
    block_stdev += (avg_block[i] - avg_all[0])*(avg_block[i] - avg_all[0]);
  }
  block_stdev /= n_blocks - 1;
  // statistical error
  double error = sqrt(block_stdev / n_blocks);
  // approximate integrated autocorrelation time
  double tau_int = 0.5 * data_per_block * block_stdev / (avg_all[1] - SQR(avg_all[0]));

  // print number of blocks, average, statistical error, and estimate of tau
  fprintf(stdout, "%4d %Lf %lf %lf\n", n_blocks, avg_all[0], error, tau_int);

  free(data);
  free(avg_block);

  return 0;
}
