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
    fprintf(stdout, "\
Info simply prints information about a system in the provided structure file. \
The verbose option prints detailed information about every molecule as \
well.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <options>\n\n", cmd);
  fprintf(ptr, "   <input>       input structure file (either vsf or vtf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -c         input coordinate file in either vcf or vtf format (default: none)\n");
  fprintf(ptr, "      -v         verbose output\n");
  fprintf(ptr, "      -h         print this help and exit\n");
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
  for (int i = 1; i < argc && argv[count][0] != '-'; i++) {
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
        strcmp(argv[i], "-c") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    fprintf(stdout, " %s", argv[i]);
  putchar('\n'); //}}}

  // options before reading system data //{{{
  // -c option - use a coordinate file //{{{
  char *input_coor = calloc(LINE,sizeof(char *));
  if (FileOption(argc, argv, "-c", &input_coor)) {
    exit(1);
  }

  // test if coordante file ends with '.vcf' or '.vtf'
  int ext;
  char extension[2][5];
  if (input_coor[0] != '\0') {
    ext = 2;
    strcpy(extension[0], ".vcf");
    strcpy(extension[1], ".vtf");
    if (!ErrorExtension(input_coor, ext, extension)) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // output verbosity //{{{
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  // }}}
  //}}}

  count = 0; // count arguments

  // <input> - input structure file (must end with .vsf or .vtf) //{{{
  char input[LINE];
  strcpy(input, argv[++count]);

  // test if <input> filename ends with '.vsf' or '.vtf' (required by VMD)
  ext = 2;
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  ReadStructure(input, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  if (verbose) {
    fprintf(stdout, "\nInformation about every bead:\n");
    PrintBead(Counts, Index, BeadType, Bead);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecule(Counts, Index, MoleculeType, Molecule, Bead, BeadType);
  }

  // get box dimensions if -c is used //{{{
  Vector BoxLength;
  BoxLength.x = -1;
  if (input_coor[0] != '\0') {
    // open input coordinate file //{{{
    FILE *vcf;
    if ((vcf = fopen(input_coor, "r")) == NULL) {
      ErrorFileOpen(input_coor, 'r');
      exit(1);
    } //}}}

    char str[LINE];
    // skip till 'pbc' keyword //{{{
    do {
      if (fscanf(vcf, "%s", str) != 1) {
        fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
        exit(1);
      }
    } while (strcmp(str, "pbc") != 0); //}}}

    // read pbc
    char line[LINE];
    fgets(line, sizeof(line), vcf);
    // split the line into array
    char *split[30];
    split[0] = strtok(line, " \t");
    int i = 0;
    while (split[i] != NULL && i < 29) {
      split[++i] = strtok(NULL, " \t");
    }
    BoxLength.x = atof(split[0]);
    BoxLength.y = atof(split[1]);
    BoxLength.z = atof(split[2]);

    fclose(vcf);
  }
  free(input_coor); //}}}

  // print information - verbose option //{{{
  VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead); //}}}

  return 0;
}
