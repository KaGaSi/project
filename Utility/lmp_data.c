#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
lmp_data utility generates lammps data file from the FIELD and coordinate \
files. It assumes molecules have bonds and can also have angles and \
dihedrals. For all cases, the a harmonic potential is assumed (the \
dihedrals are written to the lmp data file as impropers).\n\n");
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
  char *input_vsf = calloc(LINE,sizeof(char));
  strcpy(input_coor, argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
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

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType;
  MOLECULETYPE *MoleculeType;
  BEAD *Bead;
  int *Index;
  MOLECULE *Molecule;
  COUNTS Counts = InitCounts;
  VECTOR BoxLength;
  bool indexed;
  int struct_lines;
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &BoxLength, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule);
  free(input_vsf); //}}}

  // read FIELD to get bond, angle, and dihedral information //{{{
  BEADTYPE *bt_new;
  MOLECULETYPE *mt_new;
  BEAD *Bead_new;
  int *Index_new;
  MOLECULE *Molecule_new;
  COUNTS Counts_new = InitCounts;
  PARAMS *bond_type;
  PARAMS *angle_type;
  PARAMS *dihedral_type;
  ReadField(input, NULL, &Counts_new, &bt_new, &Bead_new, &Index_new,
            &mt_new, &Molecule_new,
            &bond_type, &angle_type, &dihedral_type); //}}}

  // add bond types & angles from field to all molecules from vsf //{{{
  Counts.TypesOfBonds = Counts_new.TypesOfBonds;
  Counts.TypesOfAngles = Counts_new.TypesOfAngles;
  Counts.TypesOfDihedrals = Counts_new.TypesOfDihedrals;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // identify molecule from vsf to the one from FIELD (via Name)
    int mol_field = -1;
    for (int j = 0; j < Counts_new.TypesOfMolecules; j++) {
      if (strcmp(MoleculeType[i].Name, mt_new[j].Name) == 0) {
        mol_field = j;
      }
    }
    // if the molecule is in FIELD, add bond types & angles
    if (mol_field != -1) {
      // copy bond types - assumes sorted bond arrays in both molecules which
      // should have been done while reading vsf and FIELD files
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        MoleculeType[i].Bond[j][2] = mt_new[mol_field].Bond[j][2];
      }
      // copy angles - vsf contains no angles
      MoleculeType[i].nAngles = mt_new[mol_field].nAngles;
      MoleculeType[i].Angle = calloc(MoleculeType[i].nAngles, sizeof(int *));
      for (int j = 0; j < MoleculeType[i].nAngles; j++) {
        MoleculeType[i].Angle[j] = calloc(4, sizeof(int));
        MoleculeType[i].Angle[j][0] = mt_new[mol_field].Angle[j][0];
        MoleculeType[i].Angle[j][1] = mt_new[mol_field].Angle[j][1];
        MoleculeType[i].Angle[j][2] = mt_new[mol_field].Angle[j][2];
        MoleculeType[i].Angle[j][3] = mt_new[mol_field].Angle[j][3];
      }
      // copy dihedrals - vsf contains no dihedrals
      MoleculeType[i].nDihedrals = mt_new[mol_field].nDihedrals;
      MoleculeType[i].Dihedral = calloc(MoleculeType[i].nDihedrals,
                                        sizeof(int *));
      for (int j = 0; j < MoleculeType[i].nDihedrals; j++) {
        MoleculeType[i].Dihedral[j] = calloc(5, sizeof(int));
        MoleculeType[i].Dihedral[j][0] = mt_new[mol_field].Dihedral[j][0];
        MoleculeType[i].Dihedral[j][1] = mt_new[mol_field].Dihedral[j][1];
        MoleculeType[i].Dihedral[j][2] = mt_new[mol_field].Dihedral[j][2];
        MoleculeType[i].Dihedral[j][3] = mt_new[mol_field].Dihedral[j][3];
        MoleculeType[i].Dihedral[j][4] = mt_new[mol_field].Dihedral[j][4];
      }
    // warn if the molecule type isn't in FIELD
    } else {
      fprintf(stderr, "\033[1;33m");
      fprintf(stderr, "\nWarning: molecule %s is not in %s.\n", MoleculeType[i].Name, input);
      fprintf(stderr, "         No bonds or angles will be printed to %s for molecule(s) of this type.\n\n", output);
      fprintf(stderr, "\033[0m");
    }
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    PrintBondTypes2(Counts.TypesOfBonds, bond_type);
    PrintAngleTypes2(Counts.TypesOfAngles, angle_type);
    PrintDihedralTypes2(Counts.TypesOfDihedrals, dihedral_type);
  } //}}}

  // warn if not all beads //{{{
  if (Counts.Beads != Counts.BeadsInVsf) {
    fprintf(stderr, "\033[1;33m");
    fprintf(stdout, "\nWarning: '%s' does not contain all beads from '%s'\n\n", input_coor, input_vsf);
    fprintf(stderr, "\033[0m");
  } //}}}

  // vsf file is not needed anymore

  char *stuff = calloc(LINE, sizeof(char)); // timestep preamble

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vtf, vcf, struct_lines); //}}}

  // main loop //{{{
  count = 0;
  while (true) {
    count++;
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    }

    ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);
    // if there's no additional timestep, exit the while loop
    bool rubbish; // not used
    if (ReadTimestepPreamble(&rubbish, input_coor, vcf, &stuff, false) == -1) {
      break;
    }
  }

  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count);
  }

  fclose(vcf); //}}}

  // count total number of bonds & angles //{{{
  int count_bonds = 0,
      count_angles = 0,
      count_dihedrals = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    count_bonds += MoleculeType[i].Number * MoleculeType[i].nBonds;
    count_angles += MoleculeType[i].Number * MoleculeType[i].nAngles;
    count_dihedrals += MoleculeType[i].Number * MoleculeType[i].nDihedrals;
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
  fprintf(out, "%7d impropers\n", count_dihedrals);
  if (srp) { // an extra srp bead
    fprintf(out, "%7d atom types\n", Counts.TypesOfBeads+1);
  } else {
    fprintf(out, "%7d atom types\n", Counts.TypesOfBeads);
  }
  fprintf(out, "%7d bond types\n", Counts.TypesOfBonds);
  fprintf(out, "%7d angle types\n", Counts.TypesOfAngles);
  fprintf(out, "%7d improper types\n", Counts.TypesOfDihedrals);
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

  // print dihedral (impropers in lammps speak) coefficients //{{{
  if (Counts.TypesOfDihedrals > 0) {
    fprintf(out, "Improper Coeffs\n\n");
    for (int i = 0; i < Counts.TypesOfDihedrals; i++) {
      // spring strength divided by 2, because dl_meso (FIELD) uses k/2, but lammps uses k
      fprintf(out, "%2d %lf %lf\n", i+1, dihedral_type[i].a/2, dihedral_type[i].b);
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

  // print dihedral (impropers in lammps speak) information //{{{
  if (Counts.TypesOfDihedrals > 0) {
    fprintf(out, "Impropers\n\n");
    count = 0;
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      for (int j = 0; j < MoleculeType[mtype].nDihedrals; j++) {
        count++;
        int type = MoleculeType[mtype].Dihedral[j][4];
        int id1 = Molecule[i].Bead[MoleculeType[mtype].Dihedral[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[mtype].Dihedral[j][1]];
        int id3 = Molecule[i].Bead[MoleculeType[mtype].Dihedral[j][2]];
        int id4 = Molecule[i].Bead[MoleculeType[mtype].Dihedral[j][3]];
        fprintf(out, "%7d %3d %7d %7d %7d %7d\n", count, type+1, id1+1, id2+1, id3+1, id4+1);
      }
    }
    putc('\n', out);
  } //}}}
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  FreeSystemInfo(Counts_new, &mt_new, &Molecule_new, &bt_new, &Bead_new, &Index_new);
  free(stuff);
  free(input);
  free(bond_type);
  free(angle_type);
  free(dihedral_type);
  //}}}

  return 0;
}
