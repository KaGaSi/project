#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Config_from_xyz generates dl_meso CONFIG file from given step of a xyz file. \
If the given timestep is larger than the number of steps the coordinate file, \
the last step is used. Note that box dimensions must be added manually to the \
resulting CONFIG file, as xyz file does not contain box size.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input.xyz> <options>\n\n", cmd);

  fprintf(ptr, "   <input.xyz>       input filename (xyz format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <step>     timestep for creating CONFIG (default: last)\n");
  fprintf(ptr, "      -s             no output\n");
  fprintf(ptr, "      -h             print this help and exit\n");
  fprintf(ptr, "      --script       do not reprint line (useful when output goes to file)\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 1; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  bool verbose, silent;
  SilentOption(argc, argv, &verbose, &silent);
  bool script = BoolOption(argc, argv, "--script");

  // timestep to create CONFIG file from //{{{
  int timestep = -1;
  if (IntegerOption(argc, argv, "-st", &timestep)) {
    exit(1);
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input.xyz> - filename of input xyz file (must end with .xyz) //{{{
  char input_xyz[LINE];
  strcpy(input_xyz, argv[++count]);

  // test if <input> filename ends with '.xyz' (required by VMD)
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".xyz");
  if (!ErrorExtension(input_xyz, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // get number of beads from xyz file //{{{
  // open input coordinate file
  FILE *xyz;
  if ((xyz = fopen(input_xyz, "r")) == NULL) {
    ErrorFileOpen(input_xyz, 'r');
    exit(1);
  }

  int beads;
  if (fscanf(xyz, "%d", &beads) != 1) {
    fprintf(stderr, "Error: cannot read number of beads from %s\n\n", input_xyz);
    exit(1);
  }
  fclose(xyz); //}}}

  // main loop - read timesteps till the chosen one //{{{
  // open input coordinate file //{{{
  if ((xyz = fopen(input_xyz, "r")) == NULL) {
    ErrorFileOpen(input_xyz, 'r');
    exit(1);
  } //}}}

  fpos_t pos; // for saving pointer position in xyz file
  int test;
  count = 0;

  while ((test = getc(xyz)) != EOF && count != timestep) {
    ungetc(test, xyz);

    count++;
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %6d", count);
    }

    // skip remainder of number-of-beads line
    while (getc(xyz) != '\n')
      ;
    // skip comment line
    while (getc(xyz) != '\n')
      ;

    // save pointer position in file
    fgetpos(xyz, &pos);

    // skip timestep
    for (int i = 0; i < beads; i++) {
      while (getc(xyz) != '\n')
        ;
    }
  }

  // restore pointer position in xyz file
  fsetpos(xyz, &pos); //}}}

  // read the coordinates from the chosen timestep and create CONFIG
  // open CONFIG //{{{
  FILE *config;
  if ((config = fopen("CONFIG", "w")) == NULL) {
    ErrorFileOpen("CONFIG", 'w');
    exit(1);
  } //}}}

  // print first stuff to CONFIG
  fprintf(config, "CONFIG from %s\n0 1\n", input_xyz);
  fprintf(config, "x 0 0\n");
  fprintf(config, "0 y 0\n");
  fprintf(config, "0 0 z\n");

  for (int i = 0; i < beads; i++) {
    char line[LINE];
    fgets(line, sizeof(line), xyz);

    // split the line into array //{{{
    char *split[5];
    split[0] = strtok(line, " \t");
    int j = 0;
    while (split[j] != NULL && j < 4) {
      split[++j] = strtok(NULL, " \t");
    } //}}}

    // error - less then four whitespace-separated strings //{{{
    if (j < 4) {
      fprintf(stderr, "\nError: not enough columns in %s in %d. timestep (%d. bead)\n\n", input_xyz, count, i);
    } //}}}

    // test if split[1-3] are doubles //{{{
    for (int j = 1; j < 4; j++) {
      // first character can be '-' (but must be longer) or a number
      if ((split[j][0] < '0' || split[j][0] > '9') &&
          split[j][0] != '-') {
        fprintf(stderr, "\nError: wrong %d. coordinate in %d. timestep for %d. bead (%s)\n\n", j, count, i, split[j]);
        exit(1);
      } else if (split[j][0] == '-' && strlen(split[j]) == 1) {
        fprintf(stderr, "\nError: wrong %d. coordinate in %d. timestep for %d. bead (%s)\n\n", j, count, i, split[j]);
        exit(1);
      }
      // other characters can be numbers or decimal point or newline (last character of 4th split)
      // they can also be 'e' or '-' in case of number format 1.0e-1
      for (int k = 1; k < strlen(split[j]); k++) {
        if ((split[j][k] < '0' || split[j][k] > '9') &&
            split[j][k] != '.' &&
            split[j][k] != '\n' &&
            split[j][k] != 'e' &&
            split[j][k] != '-' ) {
          fprintf(stderr, "\nError: wrong %d. coordinate in %d. timestep for %d. bead (%s)\n\n", j, count, i, split[j]);
          exit(1);
        }
      }
    } //}}}

    // print bead to CONFIG
    fprintf(config, "%s %d\n", split[0], i+1);
    fprintf(config, "%lf %lf %lf\n", atof(split[1]), atof(split[2]), atof(split[3]));
  }
  fclose(config);
  fclose(xyz);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\nConfig Step: %6d\n", count);
    }
  } //}}}

  return 0;
}
