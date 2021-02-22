#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Config utility generates CONFIG file from given step of a vcf file. If the \
given timestep is larger than the number of steps the coordinate file, the \
last step is used. Coordinate file needs to contain all beads in the \
simulation for it to work with dl_meso (no such check is made).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <step>     timestep for creating CONFIG (default: last)\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
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
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
  char *input_vsf = calloc(LINE,sizeof(char));
  strcpy(input_coor, argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // timestep to create CONFIG file from //{{{
  int timestep = -1;
  if (IntegerOption(argc, argv, "-st", &timestep)) {
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  VECTOR BoxLength; // couboid box dimensions
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &BoxLength, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // warn if not all beads //{{{
  if (Counts.Beads != Counts.BeadsInVsf) {
    YellowText(STDOUT_FILENO);
    fprintf(stdout, "\nWarning: ");
    CyanText(STDOUT_FILENO);
    fprintf(stdout, "%s", input_coor);
    YellowText(STDOUT_FILENO);
    fprintf(stdout, " does not contain all beads from ");
    CyanText(STDOUT_FILENO);
    fprintf(stdout, "%s", input_vsf);
    YellowText(STDOUT_FILENO);
    fprintf(stdout, "\n\n");
    ResetColour(STDOUT_FILENO);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  // skip structure part of a vtf
  if (vtf) {
    SkipStructVtf(vcf, input_coor);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // main loop //{{{
  fpos_t pos; // for saving pointer position in vcf file
  count = 0;
  while (true) {
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}
    // save pointer position in file
    fgetpos(vcf, &pos);

    SkipVcfCoor(vcf, input_coor, Counts, &stuff);

    if (timestep == count) {
      break;
    }

    // if there's no additional timestep, exit the while loop
    bool rubbish; // not used
    if (ReadTimestepPreamble(&rubbish, input_coor, vcf, &stuff, false) == -1) {
      break;
    }
  }

  // restore pointer position in vcf file & read the coordinates
  fsetpos(vcf, &pos);
  ReadVcfCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);
  fclose(vcf);

  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count);
  } //}}}

  // create CONFIG //{{{
  // open output CONFIG file for writing
  FILE *out;
  if ((out = fopen("CONFIG", "w")) == NULL) {
    ErrorFileOpen("CONFIG", 'w');
    exit(1);
  }

  // print CONFIG file initial stuff
  fprintf(out, "NAME\n       0       1\n");
  fprintf(out, "%lf 0.000000 0.000000\n", BoxLength.x);
  fprintf(out, "0.000000 %lf 0.000000\n", BoxLength.y);
  fprintf(out, "0.000000 0.000000 %lf\n", BoxLength.z);

  // bead coordinates
  // unbonded beads must be first (dl_meso requirement)
  for (int i = 0; i < Counts.Beads; i++) {
    fprintf(out, "%s %d\n", BeadType[Bead[i].Type].Name, i+1);
    fprintf(out, "%lf %lf %lf\n", Bead[i].Position.x-BoxLength.x/2,
                                  Bead[i].Position.y-BoxLength.y/2,
                                  Bead[i].Position.z-BoxLength.z/2);
  }
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(input_vsf);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  //}}}

  return 0;
}
