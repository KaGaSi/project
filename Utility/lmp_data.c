#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
lmp_data utility generates lammps data file from FIELD and coordinate file. \
It assumes molecules have bonds and can also have angles, but no dihedrals.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <out.data> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <out.data>        output lammps data file\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -f <name>      FIELD-like file (default: FIELD)\n");
  fprintf(ptr, "      --srp          add one more bead type for srp\n");
  fprintf(ptr, "      -st <step>     timestep for creating the output file (default: last)\n");
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
  int req_args = 2; //}}}

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
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "--srp") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  // if vtf, copy to input_vsf
  char *input_vsf = calloc(LINE,sizeof(char));
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
  } //}}}

  // <out.data> - output lammps data file //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // custom FIELD file //{{{
  char *input = calloc(LINE, sizeof(char));
  if (FileOption(argc, argv, "-f", &input)) {
    exit(1);
  }
  if (input[0] == '\0') {
    strcpy(input, "FIELD");
  } //}}}

  // add srp bead type //{{{
  bool srp = BoolOption(argc, argv, "--srp"); //}}}

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

  // variables - structures for info from vsf //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // read system information from vsf
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // variables - structures for info from FIELD //{{{
  BEADTYPE *BeadType_new; // structure with info about all bead types
  MOLECULETYPE *MoleculeType_new; // structure with info about all molecule types
  BEAD *Bead_new; // structure with info about every bead
  int *Index_new; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule_new; // structure with info about every molecule
  COUNTS Counts_new = InitCounts; // structure with number of beads, molecules, etc.
  PARAMS *bond_type; // information about bond types
  PARAMS *angle_type; // information about angle types //}}}

  // read system information from FIELD (mainly bond & angle types)
  ReadField(input, NULL, &Counts_new, &BeadType_new, &Bead_new, &Index_new, &MoleculeType_new, &Molecule_new, &bond_type, &angle_type);

  // add bond types & angles from field to all molecules from vsf //{{{
  Counts.TypesOfBonds = Counts_new.TypesOfBonds;
  Counts.TypesOfAngles = Counts_new.TypesOfAngles;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // identify molecule from vsf to the one from FIELD (via Name)
    int mol_field = -1;
    for (int j = 0; j < Counts_new.TypesOfMolecules; j++) {
      if (strcmp(MoleculeType[i].Name, MoleculeType_new[j].Name) == 0) {
        mol_field = j;
      }
    }
    // if the molecule is in FIELD, add bond types & angles
    if (mol_field != -1) {
      // copy bond types - assumes sorted bond arrays in both molecules which
      // should have been done while reading vsf and FIELD files
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        MoleculeType[i].Bond[j][2] = MoleculeType_new[mol_field].Bond[j][2];
      }
      // copy angles - vsf contains no angles
      MoleculeType[i].nAngles = MoleculeType_new[mol_field].nAngles;
      MoleculeType[i].Angle = calloc(MoleculeType[i].nAngles, sizeof(int *));
      for (int j = 0; j < MoleculeType[i].nAngles; j++) {
        MoleculeType[i].Angle[j] = calloc(4, sizeof(int));
        MoleculeType[i].Angle[j][0] = MoleculeType_new[mol_field].Angle[j][0];
        MoleculeType[i].Angle[j][1] = MoleculeType_new[mol_field].Angle[j][1];
        MoleculeType[i].Angle[j][2] = MoleculeType_new[mol_field].Angle[j][2];
        MoleculeType[i].Angle[j][3] = MoleculeType_new[mol_field].Angle[j][3];
      }
    // warn if the molecule type isn't in FIELD
    } else {
      fprintf(stderr, "\033[1;33m");
      fprintf(stderr, "\nWarning: molecule %s is not in %s.\n", MoleculeType[i].Name, input);
      fprintf(stderr, "         No bonds or angles will be printed to %s for molecule(s) of this type.\n\n", output);
      fprintf(stderr, "\033[0m");
    }
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    PrintBondTypes(Counts, bond_type);
    PrintAngleTypes(Counts, angle_type);
  } //}}}

  // warn if not all beads //{{{
  if (Counts.Beads != Counts.BeadsInVsf) {
    fprintf(stderr, "\033[1;33m");
    fprintf(stdout, "\nWarning: '%s' does not contain all beads from '%s'\n\n", input_coor, input_vsf);
    fprintf(stderr, "\033[0m");
  } //}}}

  // vsf file is not needed anymore
  free(input_vsf);

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // main loop //{{{
  fpos_t pos; // for saving pointer position in vcf file
  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF && count != timestep) {
    ungetc(test, vcf);

    count++;
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    }

    // save pointer position in file
    fpos_t pos_old = pos;
    fgetpos(vcf, &pos);

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
      pos = pos_old;
      count--;
    }
  }

  // restore pointer position in vcf file
  fsetpos(vcf, &pos);

  ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Config Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\nConfig Step: %6d\n", count);
    }
  }

  fclose(vcf); //}}}

  // count total number of bonds & angles //{{{
  int count_bonds = 0,
      count_angles = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    count_bonds += MoleculeType[i].Number * MoleculeType[i].nBonds;
    count_angles += MoleculeType[i].Number * MoleculeType[i].nAngles;
  } //}}}

  // create lammps data file //{{{
  // open output file for writing //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    ErrorFileOpen(output, 'w');
    exit(1);
  } //}}}

  // print first line that's ignored by lammps
  fprintf(out, "LAMMPS data file via lmp_data (by KaGaSi - https://github.com/KaGaSi/AnalysisTools)\n\n");

  // print number of beads, bonds, etc. //{{{
  fprintf(out, "%7d atoms\n", Counts.BeadsInVsf);
  fprintf(out, "%7d bonds\n", count_bonds);
  fprintf(out, "%7d angles\n", count_angles);
  if (srp) { // an extra srp bead
    fprintf(out, "%7d atom types\n", Counts.TypesOfBeads+1);
  } else {
    fprintf(out, "%7d atom types\n", Counts.TypesOfBeads);
  }
  fprintf(out, "%7d bond types\n", Counts.TypesOfBonds);
  fprintf(out, "%7d angle types\n", Counts.TypesOfAngles);
  putc('\n', out); //}}}

  // print box size //{{{
  fprintf(out, "0.0 %lf xlo xhi\n", BoxLength.x);
  fprintf(out, "0.0 %lf ylo yhi\n", BoxLength.y);
  fprintf(out, "0.0 %lf zlo zhi\n", BoxLength.z);
  putc('\n', out); //}}}

  // print masses of bead types //{{{
  fprintf(out, "Masses\n\n");
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(out, "%2d %lf # %s\n", i+1, BeadType[i].Mass, BeadType[i].Name);
  }
  if (srp) {
    fprintf(out, "%2d %lf # for srp\n", Counts.TypesOfBeads+1, 1.0);
  }
  putc('\n', out); //}}}

  // print bond coefficients //{{{
  if (Counts.TypesOfBonds > 0) {
    fprintf(out, "Bond Coeffs\n\n");
    for (int i = 0; i < Counts.TypesOfBonds; i++) {
      // spring strength divided by 2, because dl_meso (FIELD) uses k/2, but lammps uses k
      fprintf(out, "%2d %lf %lf\n", i+1, bond_type[i].a/2, bond_type[i].b);
    }
    putc('\n', out);
  } //}}}

  // print angle coefficients //{{{
  if (Counts.TypesOfAngles > 0) {
    fprintf(out, "Angle Coeffs\n\n");
    for (int i = 0; i < Counts.TypesOfAngles; i++) {
      // spring strength divided by 2, because dl_meso (FIELD) uses k/2, but lammps uses k
      fprintf(out, "%2d %lf %lf\n", i+1, angle_type[i].a/2, angle_type[i].b);
    }
    putc('\n', out);
  } //}}}

  // print bead coordinates //{{{
  fprintf(out, "Atoms\n\n");
  // coordinates of unbonded beads
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    int type = Bead[i].Type;
    fprintf(out, "%7d", i+1);
    if (Bead[i].Molecule == -1) { // bead not in molecule
      fprintf(out, "%7d", 0);
    } else { // bead in a molecule
      fprintf(out, "%7d", Bead[i].Molecule+1);
    }
    fprintf(out, "%4d   %lf  ", type+1, BeadType[type].Charge);
    if (Bead[i].Position.x == BoxLength.x) {
      Bead[i].Position.x -= 0.00001;
    }
    if (Bead[i].Position.y == BoxLength.y) {
      Bead[i].Position.y -= 0.00001;
    }
    if (Bead[i].Position.z == BoxLength.z) {
      Bead[i].Position.z -= 0.00001;
    }
    fprintf(out, "%lf %lf %lf\n", Bead[i].Position.x,
                                  Bead[i].Position.y,
                                  Bead[i].Position.z);
  }
  putc('\n', out); //}}}

  // print bond information //{{{
  if (count_bonds > 0) {
    fprintf(out, "Bonds\n\n");
    count = 0;
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      for (int j = 0; j < MoleculeType[mtype].nBonds; j++) {
        count++;
        int type = MoleculeType[mtype].Bond[j][2];
        int id1 = Molecule[i].Bead[MoleculeType[mtype].Bond[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[mtype].Bond[j][1]];
        fprintf(out, "%7d %3d %7d %7d\n", count, type+1, id1+1, id2+1);
      }
    }
    putc('\n', out);
  } //}}}

  // print angle information //{{{
  if (count_angles > 0) {
    fprintf(out, "Angles\n\n");
    count = 0;
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      for (int j = 0; j < MoleculeType[mtype].nAngles; j++) {
        count++;
        int type = MoleculeType[mtype].Angle[j][3];
        int id1 = Molecule[i].Bead[MoleculeType[mtype].Angle[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[mtype].Angle[j][1]];
        int id3 = Molecule[i].Bead[MoleculeType[mtype].Angle[j][2]];
        fprintf(out, "%7d %3d %7d %7d %7d\n", count, type+1, id1+1, id2+1, id3+1);
      }
    }
    putc('\n', out);
  } //}}}

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);

  free(BeadType_new);
  free(Index_new);
  FreeMoleculeType(Counts_new, &MoleculeType_new);
  FreeMolecule(Counts_new, &Molecule_new);
  FreeBead(Counts_new, &Bead_new);

  free(stuff);
  free(input);
  free(angle_type);
  free(bond_type);
  //}}}

  return 0;
}
