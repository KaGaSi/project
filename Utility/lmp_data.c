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
lmp_data utility generates lammps data file from FIELD and coordinate file. \
It assumes molecules have bonds and can also have angles, but no dihedrals.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <out.data> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <out.data>        output lammps data file\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -f <name>      FIELD-like file (default: FIELD)\n");
  fprintf(ptr, "      --srp          add one more bead type for srp");
  fprintf(ptr, "      -st <step>     timestep for creating the output file (default: last)\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
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
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--srp") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  char *input_vsf = calloc(LINE,sizeof(char));
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

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}
  //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.data> - output lammps data file //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[LINE];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_coor);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // warn if not all beads //{{{
  if (Counts.Beads != Counts.BeadsInVsf) {
    fprintf(stdout, "\nWarning: '%s' does not contain all beads from '%s'\n\n", input_coor, input_vsf);
  } //}}}

  // vsf file is not needed anymore
  free(input_vsf);

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // main loop //{{{
  fpos_t pos, pos_old; // for saving pointer position in vcf file
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
    pos_old = pos;
    fgetpos(vcf, &pos);

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\nError: premature end of %s file\n\n", input_coor);
      pos = pos_old;
      count--;
    }
  }

  // restore pointer position in FIELD file
  fsetpos(vcf, &pos);

  // read coordinates //{{{
  if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
    // print newline to stdout if Step... doesn't end with one
    ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
    exit(1);
  } //}}}

  if (!silent) {
    if (script) {
      fprintf(stdout, "Config Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\nConfig Step: %6d\n", count);
    }
  }

  fclose(vcf); //}}}

  // count number of bonds //{{{
  int count_bonds = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    count_bonds += MoleculeType[i].Number * MoleculeType[i].nBonds;
  } //}}}

  // read stuff from FIELD //{{{
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // read till molecule keyword //{{{
  char line[LINE];
  while(fgets(line, sizeof(line), fr)) {
    char *split;
    split = strtok(line, " \t ");
    if (strncmp(split, "molecule", 8) == 0 ||
        strncmp(split, "Molecule", 8) == 0 ||
        strncmp(split, "MOLECULE", 8) == 0 ) {
      break;
    }
  } //}}}

  // read data for each molecule type
  int count_bond_types = 0, count_angle_types = 0, count_angles = 0;
  // info about bond types (to be realloc'd)
  double **bond_type = malloc(1*sizeof(double *));
  bond_type[0] = malloc(2*sizeof(double));
  bond_type[0][0] = -1; // k (harm spring strength)
  bond_type[0][1] = -1; // r_0 (equilibrium distance)
  // info about angle types (to be realloc'd)
  double **angle_type = malloc(1*sizeof(double *));
  angle_type[0] = malloc(2*sizeof(double));
  angle_type[0][0] = -1; // k (harm spring strength)
  angle_type[0][1] = -1; // angle_0 (equilibrium angle)
  // info about bond and angle types for each molecule
  int *mol_bond_types[Counts.TypesOfMolecules];
  int *mol_angle_types[Counts.TypesOfMolecules];
  // beads in individual angles
  int angle_beads_n[Counts.TypesOfMolecules];
  int **angle_beads[Counts.TypesOfMolecules];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    mol_bond_types[i] = calloc(MoleculeType[i].nBonds, sizeof(int));
    mol_angle_types[i] = calloc(MoleculeType[i].nBonds, sizeof(int));

    // read till bond keyword //{{{
    char *split[3];
    while(fgets(line, sizeof(line), fr)) {
      split[0] = strtok(line, " \t ");
      if (strncmp(split[0], "bond", 4) == 0 ||
          strncmp(split[0], "Bond", 4) == 0 ||
          strncmp(split[0], "BOND", 4) == 0 ) {
        split[0] = strtok(NULL, " \t"); // number of bonds for the molecule
        break;
      }
    } //}}}

    // count bond types //{{{
    int exist = -1;
    count = atoi(split[0]);
    for (int j = 0; j < count; j++) {
      fgets(line, sizeof(line), fr);
      split[0] = strtok(line, " \t"); // 'harm' keyword
      split[0] = strtok(NULL, " \t"); // bead ids
      split[0] = strtok(NULL, " \t");
      split[0] = strtok(NULL, " \t"); // spring strength
      split[1] = strtok(NULL, " \t"); // bond length
      exist = -1;
      for (int k = 0; k < count_bond_types; k++) {
        if (bond_type[k][0] == atof(split[0]) && bond_type[k][1] == atof(split[1])){
          exist = k;
        }
      }
      if (exist == -1) {
        count_bond_types++;
        bond_type = realloc(bond_type, count_bond_types*sizeof(double *));
        bond_type[count_bond_types-1] = malloc(2*sizeof(double));
        bond_type[count_bond_types-1][0] = atof(split[0]);
        bond_type[count_bond_types-1][1] = atof(split[1]);
        exist = count_bond_types - 1;
      }

      mol_bond_types[i][j] = exist;
    } //}}}

    // read angles (if present) //{{{
    angle_beads_n[i] = -1;
    exist = -1;
    while(fgets(line, sizeof(line), fr)) {
      split[0] = strtok(line, " \t ");
      if (strncmp(split[0], "angle", 5) == 0 ||
          strncmp(split[0], "Angle", 5) == 0 ||
          strncmp(split[0], "ANGLE", 5) == 0 ) {
        split[0] = strtok(NULL, " \t"); // number of angles for the molecule
        exist = 1;
        break;
      } else if (strncmp(split[0], "finish", 6) == 0 ||
                 strncmp(split[0], "Finish", 6) == 0 ||
                 strncmp(split[0], "FINISH", 6) == 0 ) {
        break;
      }
    } //}}}

    // count angle types & angles (if any) //{{{
    if (exist == 1) {
      angle_beads_n[i] = atoi(split[0]);
      angle_beads[i] = malloc(angle_beads_n[i]*sizeof(int *));
      count_angles += angle_beads_n[i] * MoleculeType[i].Number;
      for (int j = 0; j < angle_beads_n[i]; j++) {
        angle_beads[i][j] = calloc(3, sizeof(int));
        fgets(line, sizeof(line), fr);
        split[0] = strtok(line, " \t"); // 'harm' keyword
        split[0] = strtok(NULL, " \t"); // bead ids
        angle_beads[i][j][0] = atoi(split[0]) - 1;
        split[0] = strtok(NULL, " \t");
        angle_beads[i][j][1] = atoi(split[0]) - 1;
        split[0] = strtok(NULL, " \t");
        angle_beads[i][j][2] = atoi(split[0]) - 1;
        split[0] = strtok(NULL, " \t"); // spring strength
        split[1] = strtok(NULL, " \t"); // bond length
        // find if the bond type alread exists
        exist = -1;
        for (int k = 0; k < count_angle_types; k++) {
          if (angle_type[k][0] == atof(split[0]) && angle_type[k][1] == atof(split[1])) {
            exist = k;
          }
        }
        if (exist == -1) { // add a new bond type if this doesn't exist
          count_angle_types++;
          angle_type = realloc(angle_type, count_angle_types*sizeof(double *));
          angle_type[count_angle_types-1] = malloc(2*sizeof(double));
          angle_type[count_angle_types-1][0] = atof(split[0]);
          angle_type[count_angle_types-1][1] = atof(split[1]);
          exist = count_angle_types - 1;
        }

        mol_angle_types[i][j] = exist;
      }
    } //}}}
  }
  fclose(fr); //}}}

  // print information about bonds and angles (verbouse output) //{{{
  if (verbose) {
    fprintf(stdout, "Bond types:\n");
    for (int i = 0; i < count_bond_types; i++) {
      fprintf(stdout, " %lf %lf\n", bond_type[i][0], bond_type[i][1]);
    }
    fprintf(stdout, "Bonds in molecules (<type> <bead 1> <bead 2>):\n");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(stdout, " %s (%d bonds)\n", MoleculeType[i].Name, MoleculeType[i].nBonds);
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        fprintf(stdout, "  %2d: %3d %3d %3d\n", j+1, mol_bond_types[i][j]+1, MoleculeType[i].Bond[j][0]+1, MoleculeType[i].Bond[j][1]+1);
      }
    }
    fprintf(stdout, "Angle types:\n");
    for (int i = 0; i < count_angle_types; i++) {
      fprintf(stdout, " %lf %lf\n", angle_type[i][0], angle_type[i][1]);
    }
    fprintf(stdout, "Angles in molecules (<type> <bead 1> <bead 2> <bead 3>):\n");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(stdout, " %s ", MoleculeType[i].Name);
      if (angle_beads_n[i] == -1) {
        fprintf(stdout, "(no angles)\n");
      } else {
        fprintf(stdout, "(%d angles)\n", angle_beads_n[i]);
      }
      for (int j = 0; j < angle_beads_n[i]; j++) {
        int id1 = angle_beads[i][j][0];
        int id2 = angle_beads[i][j][1];
        int id3 = angle_beads[i][j][2];
        fprintf(stdout, "  %2d: %3d %3d %3d %3d\n", j+1, mol_angle_types[i][j]+1, id1+1, id2+1, id3+1);
      }
    }
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
  fprintf(out, "%7d bond types\n", count_bond_types);
  fprintf(out, "%7d angle types\n", count_angle_types);
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
  if (count_bond_types > 0) {
    fprintf(out, "Bond Coeffs\n\n");
    for (int i = 0; i < count_bond_types; i++) {
      // spring strength divided by 2, because dl_meso (FIELD) uses k/2, but lammps uses k
      fprintf(out, "%2d %lf %lf\n", i+1, bond_type[i][0]/2, bond_type[i][1]);
    }
    putc('\n', out);
  } //}}}

  // print angle coefficients //{{{
  if (count_angle_types > 0) {
    fprintf(out, "Angle Coeffs\n\n");
    for (int i = 0; i < count_angle_types; i++) {
      // spring strength divided by 2, because dl_meso (FIELD) uses k/2, but lammps uses k
      fprintf(out, "%2d %lf %lf\n", i+1, angle_type[i][0]/2, angle_type[i][1]);
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
    fprintf(out, "%lf %lf %lf\n", Bead[i].Position.x,
                                  Bead[i].Position.y,
                                  Bead[i].Position.z);
  }
  putc('\n', out); //}}}

  // print bond information //{{{
  if (count_bond_types != 0) {
    fprintf(out, "Bonds\n\n");
    count = 0;
    for (int i = 0; i < Counts.Molecules; i++) {
      int type = Molecule[i].Type;
      for (int j = 0; j < MoleculeType[type].nBonds; j++) {
        count++;
        int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];
        fprintf(out, "%7d %3d %7d %7d\n", count, mol_bond_types[type][j]+1, id1+1, id2+1);
      }
    }
    putc('\n', out);
  } //}}}

  // print angle information //{{{
  if (count_angle_types != 0) {
    fprintf(out, "Angles\n\n");
    count = 0;
    for (int i = 0; i < Counts.Molecules; i++) {

      int type = Molecule[i].Type;
      for (int j = 0; j < angle_beads_n[type]; j++) {
        count++;
        int id1 = Molecule[i].Bead[angle_beads[type][j][0]];
        int id2 = Molecule[i].Bead[angle_beads[type][j][1]];
        int id3 = Molecule[i].Bead[angle_beads[type][j][2]];
        fprintf(out, "%7d %3d %7d %7d %7d\n", count, mol_angle_types[type][j]+1, id1+1, id2+1, id3+1);
      }
    }
  }//}}}

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(input);
  //}}}

  return 0;
}
