#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Test reading/writing coordinates and structure.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> ", cmd);
  fprintf(ptr, "<output.vcf> <output.vsf> <output.xyz> <options>\n\n");

  fprintf(ptr, "   <input>           input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <output.vcf>      output coordinate file (vcf format)\n");
  fprintf(ptr, "   <output.vsf>      output structure file (vsf format)\n");
  fprintf(ptr, "   <output.xyz>      output coordinate file (xyz format)\n");
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
  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  // reverse bead type selection? ...do now to check correct number of arguments
  bool reverse = BoolOption(argc, argv, "-r");

  // possible to omit <type name(s)> if '-r' is used
  if (count < (req_args-1) || (count == (req_args-1) && !reverse)) {
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
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0) {

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

  // <output.vcf> - filename of output vcf file //{{{
  char output_vcf[LINE];
  strcpy(output_vcf, argv[++count]);
  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output.vsf> - filename of output vsf file //{{{
  char output_vsf[LINE];
  strcpy(output_vsf, argv[++count]);

  // test if <output.vsf> filename ends with '.vsf'
  ext = 1;
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output.xyz> - filename of output xyz file //{{{
  char output_xyz[LINE];
  strcpy(output_xyz, argv[++count]);

  // test if <output.xyz> filename ends with '.xyz'
  ext = 1;
  strcpy(extension[0], ".xyz");
  if (ErrorExtension(output_xyz, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  }

  // open & close xyz file - later, there's only appending
  FILE *out;
  if ((out = fopen(output_xyz, "w")) == NULL) {
    ErrorFileOpen(output_xyz, 'w');
    exit(1);
  }
  fclose(out); //}}}

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read structure & first timestep  //{{{
  // get box size
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  VECTOR BoxLength = GetPBC(vcf, input_coor);
  fclose(vcf);
  // read the whole structure section
  ReadVtfStructure(input_vsf, true, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);
  // count structure lines if vtf as the coordinate file (-1 otherwise)
  int struct_lines = CountVtfStructLines(vtf, input_coor);
  // determine timestep type & what coordinate file contains from the first timestep
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vtf, vcf, struct_lines);
  bool indexed = CheckVtfTimestep(vcf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);
  fclose(vcf);
  free(input_vsf); //}}}

  // write all bead & molecule types
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = true;
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = true;
  }

  // open input coordinate file //{{{
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  SkipVtfStructure(vtf, vcf, struct_lines);

  // print information - verbose output //{{{
  VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  if (verbose) {
    PrintBead2(Counts.Beads, Index, BeadType, Bead);
    PrintMolecule(Counts.Molecules, MoleculeType, Molecule, BeadType, Bead);
  } //}}}

  // print pbc to output .vcf file //{{{
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  fclose(out); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // main loop //{{{
  while (true) {

    ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

    // write vcf coordinates
    if ((out = fopen(output_vcf, "a")) == NULL) {
      ErrorFileOpen(output_vcf, 'a');
      exit(1);
    }
    WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
    fclose(out);

    // write xyz coordinates
    if ((out = fopen(output_xyz, "a")) == NULL) {
      ErrorFileOpen(output_xyz, 'a');
      exit(1);
    }
    WriteCoorXYZ(out, Counts, BeadType, Bead);
    fclose(out);

    // if there's no additional timestep, exit the while loop
    bool test; // not used
    if (ReadTimestepPreamble(&test, input_coor, vcf, &stuff) == -1) {
      break;
    }
  }
  fclose(vcf); //}}}

  WriteVsf(output_vsf, Counts, BeadType, Bead, MoleculeType, Molecule, false);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead2(Counts.Beads, &Bead);
  free(stuff);
  //}}}

  return 0;
}
