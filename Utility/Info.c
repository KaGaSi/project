#include "../AnalysisTools.h"

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
  fprintf(ptr, "   <input>       vtf input structure file\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -c         vtf input coordinate file (default: none)\n");
  fprintf(ptr, "      --detailed differentiate bead types not just by names\n");
  fprintf(ptr, "      -v         verbose output\n");
  fprintf(ptr, "      -h         print this help and exit\n");
  fprintf(ptr, "      --version  print version number and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(EXIT_SUCCESS);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(EXIT_SUCCESS);
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
    exit(EXIT_FAILURE);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-c") != 0 &&
        strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(EXIT_FAILURE);
    }
  } //}}}

  count = 0; // count arguments

  // <input> - input structure file (must end with .vsf or .vtf) //{{{
  char *input_vsf = calloc(LINE,sizeof(char));
  strcpy(input_vsf, argv[++count]);
  // test if <input> filename ends with '.vsf' or '.vtf'
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // options before reading system data //{{{
  // -c option - use a coordinate file //{{{
  char *input_coor = calloc(LINE,sizeof(char));
  if (FileOption(argc, argv, "-c", &input_coor)) {
    exit(1);
  }
  bool vtf = false;
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (input_coor[0] != '\0') {
    int test;
    if ((test=ErrorExtension(input_coor, ext, extension)) == -1) {
      Help(argv[0], true);
      exit(1);
    } else if (test == 1) {
      vtf = true;
    }
  } //}}}

  bool verbose = BoolOption(argc, argv, "-v"); // verbose output?
  bool detailed = BoolOption(argc, argv, "--detailed"); // verbose output?
  //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  VECTOR BoxLength = {-1, -1, -1}; // couboid box dimensions
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, detailed, vtf, &indexed, &struct_lines,
              &BoxLength, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule);
  free(input_vsf); //}}}

  // print information /{{{
  VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule); //}}}

  if (verbose) { //{{{
    fprintf(stdout, "\nInformation about every bead:\n");
    PrintBead(Counts, Index, BeadType, Bead);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecule(Counts.Molecules, MoleculeType, Molecule, BeadType, Bead);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(input_coor); //}}}

  return 0;
}
