#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
lmp_vsf takes lammps data file and transforms it into vsf and vcf files. \
Different bead types are named bead1 to beadN (where N is the number of bead \
types) or as specified in comment in the Masses section of the data file. \
Molecule types are named mol1 to molN (where N is the number of \
molecule types). Molecule types are determined according to bead order as \
well as according to molecule connectivity. Angles are disregarded.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <out.vsf> <out.vcf> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input lammps data file\n");
  fprintf(ptr, "   <out.vsf>         output vsf structure file\n");
  fprintf(ptr, "   <out.vcf>         output vcf coordinate file\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -v             verbose output\n");
  fprintf(ptr, "      -h             print this help and exit\n");
  fprintf(ptr, "      --version      print version number and exit\n");
} //}}}

// TODO: make data file reading into a separate function in Read.c

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
  int req_args = 3; //}}}

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
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // verbose output (-v option)
  bool verbose = BoolOption(argc, argv, "-v");

  // print command to stdout //{{{
  bool silent = false;
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input[LINE];
  strcpy(input, argv[++count]); //}}}

  // <out.vsf> - output structure file //{{{
  char output_vsf[LINE];
  strcpy(output_vsf, argv[++count]);

  // test if <out.vsf> filename ends with '.vsf' (required by VMD)
  int ext = 1;
  char extension[1][5];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output_vsf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output structure file //{{{
  char output_vcf[LINE];
  strcpy(output_vcf, argv[++count]);

  // test if <out.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // some variables //{{{
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  MOLECULE *Molecule;
  MOLECULETYPE *MoleculeType;
  VECTOR BoxLength;
  VECTOR box_lo; // {x,y,z}lo from data file to place beads in (0, BoxLength>
  int bonds = 0; // total number of bonds
  int angles = 0; // total number of angles
  PARAMS *bond_type; // information about bond types
  PARAMS *angle_type; // information about angle types //}}}

  // open <input> //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // ignore first line (comment) //{{{
  char line[LINE];
  fgets(line, sizeof(line), fr); //}}}

  // read data file header //{{{
  // data file header lines must start with a number (or '#' for comment),
  // therefore read until something else is encountered
  char split[30][100], delim[8];
  int words;
  strcpy(delim, " \t");
  do {
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // number of atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;33m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'atoms' keyword must be preceded by integer\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      Counts.BeadsInVsf = atoi(split[0]);
      Counts.Beads = atoi(split[0]);
    } //}}}
    // number of bonds //{{{
    if (words > 1 && strcmp(split[1], "bonds") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bonds' keyword must be preceded by integer\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      bonds = atoi(split[0]);
    } //}}}
    // number of angles //{{{
    if (words > 1 && strcmp(split[1], "angles") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'angles' keyword must be preceded by integer\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      angles = atoi(split[0]);
      if (angles);
    } //}}}
    // number of bead types //{{{
    if (words > 2 && strcmp(split[1], "atom") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'atom types' keyword must be preceded by integer\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      Counts.TypesOfBeads = atoi(split[0]);
    } //}}}
    // number of bond types //{{{
    if (words > 2 && strcmp(split[1], "bond") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bond types' keyword must be preceded by integer\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      Counts.TypesOfBonds = atoi(split[0]);
    } //}}}
    // number of angle types //{{{
    if (words > 2 && strcmp(split[1], "angle") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'angle types' keyword must be preceded by integer\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      Counts.TypesOfAngles = atoi(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "xlo") == 0 && strcmp(split[3], "xhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'xlo xhi' keyword must be preceded by two floats\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      BoxLength.x = atof(split[1]) - atof(split[0]);
      box_lo.x = atof(split[0]);
    } //}}}
    // box length in y //{{{
    if (words > 3 && strcmp(split[2], "ylo") == 0 && strcmp(split[3], "yhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'ylo yhi' keyword must be preceded by two floats\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      BoxLength.y = atof(split[1]) - atof(split[0]);
      box_lo.y = atof(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "zlo") == 0 && strcmp(split[3], "zhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'zlo zhi' keyword must be preceded by two floats\n", input);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      }
      BoxLength.z = atof(split[1]) - atof(split[0]);
      box_lo.z = atof(split[0]);
    } //}}}
  } while (words == 0 ||
           split[0][0] == '#' ||
           IsDouble(split[0]) ||
           IsInteger(split[0])); //}}}

  // some error checking //{{{
  if (Counts.TypesOfBeads == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'atom types' line (or is 0)\n\n", input);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if (Counts.BeadsInVsf == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'atoms' line (or is 0)\n\n", input);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if (BoxLength.x == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'xlo xhi' line (or is 0 0)\n\n", input);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if (BoxLength.y == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'ylo yhi' line (or is 0 0)\n\n", input);
    fprintf(stderr, "\033[0m");
    exit(1);
  }
  if (BoxLength.z == 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError \033[1;33m%s\033[1;31m - missing 'zlo zhi' line (or is 0 0)\n\n", input);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // fill something in BeadType struct //{{{
  BEADTYPE *BeadType = calloc(Counts.TypesOfBeads, sizeof(struct BeadType));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    sprintf(BeadType[i].Name, "bead%d", i+1);
    BeadType[i].Use = true;
    BeadType[i].Write = true;
  } //}}}

  // bead struct memory allocation //{{{
  BEAD *Bead = calloc(Counts.Beads, sizeof(struct Bead));
  for (int i = 0; i < Counts.Beads; i++) {
    Bead[i].Aggregate = calloc(1, sizeof(int));
  } //}}}

  int *Index = calloc(Counts.Beads, sizeof(int)); // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  bond_type = calloc(Counts.TypesOfBonds, sizeof(PARAMS));
  angle_type = calloc(Counts.TypesOfAngles, sizeof(PARAMS));

  // read body of data file //{{{
  int test,
      *mols, // number of beads in each molecule
      monomer = 0; // monomer beads designated by mol_ID = 0 in data file
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);
    // atom masses //{{{
    if (words > 0 && strcmp(split[0], "Masses") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < Counts.TypesOfBeads; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 2 || !IsInteger(split[0]) || !IsPosDouble(split[1])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each line in 'Masses' section must start with '<int> <float>'\n", input);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        BeadType[i].Mass = atof(split[1]);
        // if there's a comment at the end of the line, consider it bead name
        if (words > 2 && split[2][0] == '#') {
          if (strlen(split[2]) == 1 && words > 3) { // comment form '# name'
            strcpy(BeadType[i].Name, split[3]);
          } else if (strlen(split[2]) > 1) { // comment form '#name'
            for (int j = 0; j < strlen(split[2]); j++) {
              split[2][j] = split[2][j+1];
            }
            strcpy(BeadType[i].Name, split[2]);
          }
        }
      }
    } //}}}
    // bond coefficients //{{{
    if (words > 1 && strcmp(split[0], "Bond") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < Counts.TypesOfBonds; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosDouble(split[1]) || !IsPosDouble(split[2])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each line in 'Bond Coeffs' section must start with '<int> <float> <float>'\n", input);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        bond_type[i].a = atof(split[1]);
        bond_type[i].b = atof(split[2]);
      }
    } //}}}
    // angle coefficients //{{{
    if (words > 1 && strcmp(split[0], "Angle") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < Counts.TypesOfAngles; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosDouble(split[1]) || !IsPosDouble(split[2])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each line in 'Angle Coeffs' section must start with '<int> <float> <float>'\n", input);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        angle_type[i].a = atof(split[1]);
        angle_type[i].b = atof(split[2]);
      }
    } //}}}
    // atoms section //{{{
    if (words > 0 && strcmp(split[0], "Atoms") == 0) {
      // array for number of beads in each molecule
      mols = calloc(Counts.Beads, sizeof(int));
      // array for list of beads in each molecule
      // bead_mols[i] ... molecule's id; bead_mols[][i] ... molecule's beads
      int **bead_mols = calloc(Counts.Beads, sizeof(int *));
      for (int i = 0; i < Counts.Beads; i++) {
        bead_mols[i] = calloc(1, sizeof(int));
      }
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // go through atom section to get basic info //{{{
      fpos_t pos; // set file counter
      fgetpos(fr, &pos); // save file pointer
      for (int i = 0; i < Counts.Beads; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <id> <mol_id> <btype> <charge> <x> <y> <z>
        // Error - incorrect format //{{{
        if (words < 7 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) || !IsInteger(split[2]) ||
            !IsDouble(split[3]) || !IsDouble(split[4]) || !IsDouble(split[5]) || !IsDouble(split[6])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each 'Atoms' line must be <id> <mol_id> <bead type> <charge> <x> <y> <z>\n", input);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int id = atoi(split[0]) - 1, // in lammps, these start at 1
            mol_id = atoi(split[1]) - 1, // in lammps, molecules start with 1; unbonded atoms can be 0
            btype = atoi(split[2]) - 1; // in lammps, these start at 1

        if (mol_id == -1) { // corresponds to 0 in the data file
          monomer++;
          Bead[id].Molecule = -1;
        } else { // possibly in a molecule (if more beads share its mol_id)
          mols[mol_id]++;
          bead_mols[mol_id] = realloc(bead_mols[mol_id], mols[mol_id]*sizeof(int));
          bead_mols[mol_id][mols[mol_id]-1] = id;
          Bead[id].Molecule = mol_id;
        }
        BeadType[btype].Charge = atof(split[3]);
        BeadType[btype].Number++;
        Bead[id].Position.x = atof(split[4]) - box_lo.x;
        Bead[id].Position.y = atof(split[5]) - box_lo.y;
        Bead[id].Position.z = atof(split[6]) - box_lo.z;
        Bead[id].Type = btype;
        Bead[id].Index = id;
        Index[id] = id;
      } //}}}
      Counts.Unbonded = monomer;
      // go through possible molecules and remove 1-bead molecules //{{{
      count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < Counts.Beads; i++) {
        if (mols[i] == 1) {
          int bead = bead_mols[i][0];
          Bead[bead].Molecule = -1;
          Counts.Unbonded++;
        } else if (mols[i] > 1){
          Counts.Molecules++;
          Counts.Bonded += mols[i];
          for (int j = 0; j < mols[i]; j++) {
            int bead = bead_mols[i][j];
            Bead[bead].Molecule = count;
          }
          count++;
        }
      } //}}}
      // TODO: is the 'remove 1-bead...' necessary? Join with 'allocate Molecule struct...'
      // remove 1-bead molecules from mols array and sort bead_mols arrayy according to ascending bead index //{{{
      count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < Counts.Beads; i++) {
        if (mols[i] > 1) {
          mols[count] = mols[i];
          // realloc bead_mols[count] to bead_mols[i] as it realloc'd in the first place
          bead_mols[count] = realloc(bead_mols[count], mols[count]*sizeof(int));
          for (int j = 0; j < mols[i]; j++) {
            bead_mols[count][j] = bead_mols[i][j];
          }
          // sort molecules in bead_mols[count][] according to ascending id
          SortArray(&bead_mols[count], mols[i], 0);
          count++;
        }
      } //}}}
      // zeroize unused part of mols array - just to be on the save side //{{{
      for (int i = Counts.Molecules; i < Counts.Beads; i++) {
        mols[i] = 0;
      } //}}}
      // allocate Molecule struct and fill Molecule[].Bead array //{{{
      Molecule = calloc(Counts.Molecules, sizeof(struct Molecule));
      for (int i = 0; i < Counts.Molecules; i++) {
        Molecule[i].Bead = calloc(mols[i], sizeof(int));
        for (int j = 0; j < mols[i]; j++) {
          Molecule[i].Bead[j] = bead_mols[i][j];
        }
      } //}}}
      // free helper array  //{{{
      for (int i = 0; i < Counts.Beads; i++) {
        free(bead_mols[i]);
      }
      free(bead_mols); //}}}
    } //}}}
    // bonds section //{{{
    if (words > 0 && strcmp(split[0], "Bonds") == 0) {
      // allocate helper arrays to hold bond info //{{{
      // number of bonds in each molecule
      int *bonds_per_mol = calloc(Counts.Molecules, sizeof(int));
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] ... ids of connected beads in molecule 'i'
      // [i][][2] ... bond type
      int ***mol_bonds = calloc(Counts.Molecules, sizeof(int **));
      for (int i = 0; i < Counts.Molecules; i++) {
         mol_bonds[i] = calloc(1, sizeof(int *));
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // read all bonds //{{{
      for (int i = 0; i < bonds; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <bond id> <bond type> <bead1> <bead2>
        // Error - incorrect format //{{{
        if (words < 4 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each 'Bonds' line must be <bond id> <bond type> <bead1d> <bead2>\n", input);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        // assign molecule to the bond
        int mol = Bead[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != Bead[bead2].Molecule) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input);
          fprintf(stderr, " - beads in \033[1;33mbond %d\033[1;31m are in different molecules\n\n", atoi(split[0]));
          fprintf(stderr, "\033[0m");
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        bonds_per_mol[mol]++;
        int bond = bonds_per_mol[mol];
        // add bonded beads to the molecule they belong to
        mol_bonds[mol] = realloc(mol_bonds[mol], bond*sizeof(int *));
        mol_bonds[mol][bond-1] = calloc(3, sizeof(int));
        mol_bonds[mol][bond-1][0] = bead1;
        mol_bonds[mol][bond-1][1] = bead2;
        mol_bonds[mol][bond-1][2] = type;
      } //}}}
      // sort bonds according to the id of the first bead in a bond //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        SortBonds(mol_bonds[i], bonds_per_mol[i]);
      } //}}}
      // minimize mol_bonds based on lowest id in each molecule //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        int lowest = Counts.Beads; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if (Molecule[i].Bead[j] < lowest) {
            lowest = Molecule[i].Bead[j];
          }
        }
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          mol_bonds[i][j][0] -= lowest;
          mol_bonds[i][j][1] -= lowest;
        }
      } //}}}
      // identify molecule type based on bead order and connectivity //{{{
      MoleculeType = calloc(1, sizeof(struct MoleculeType));
      // number of molecule types differing in connectivity but not in beads - naming purposes
      int *diff_conn = calloc(1, sizeof(int));
      // number of molecule types differing in bead order - naming purposes
      int mtype_bead = 0;
      for (int i = 0; i < Counts.Molecules; i++) {
        // is molecule 'i' of known type?
        bool exists = false;
        // do 'i' and given molecule type share connectivity and bead order?
        bool same_conn = true, same_bead = true;
        // molecule type with which 'i' shares bead order - naming purposes
        int type_bead = -1;
        for (int j = 0; j < Counts.TypesOfMolecules; j++) { //{{{
          if (MoleculeType[j].nBeads == mols[i] && // same number of molecules?
              MoleculeType[j].nBonds == bonds_per_mol[i]) { // same number of bonds?
            // same connectivity?
            same_conn = true;
            for (int k = 0; k < MoleculeType[j].nBonds; k++) {
              if (MoleculeType[j].Bond[k][0] != mol_bonds[i][k][0] ||
                  MoleculeType[j].Bond[k][1] != mol_bonds[i][k][1] ||
                  MoleculeType[j].Bond[k][2] != mol_bonds[i][k][2]) {
                same_conn = false;
                break;
              }
            }
            // same bead types?
            same_bead = true;
            for (int k = 0; k < MoleculeType[j].nBeads; k++) {
              int btype = Bead[Molecule[i].Bead[k]].Type;
              if (MoleculeType[j].Bead[k] != btype) {
                same_bead = false;
                break;
              }
            }
            if (same_bead && type_bead == -1) {
              type_bead = j;
            }
            // if the molecule has the same connectivity and bead order, it's of known type
            if (same_conn && same_bead) {
              exists = true;
              MoleculeType[j].Number++;
              Molecule[i].Type = j;
              break;
            }
          }
        } //}}}
        // add new type? //{{{
        if (!exists) {
          int mtype = Counts.TypesOfMolecules;
          MoleculeType = realloc(MoleculeType, (mtype+1)*sizeof(struct MoleculeType));
          diff_conn = realloc(diff_conn, (mtype+1)*sizeof(int));
          diff_conn[mtype] = 0;
          // molecule name
          if (same_bead && !same_conn) { // same beads - mol<type_bead>-b#
            diff_conn[type_bead]++;
            char name[5];
            strcpy(name, MoleculeType[type_bead].Name);
            sprintf(MoleculeType[mtype].Name, "%s-b%d", name, diff_conn[type_bead]);
          } else { // same connectivity or both different - mol#
            mtype_bead++;
            sprintf(MoleculeType[mtype].Name, "mol%d", mtype_bead);
          }
          Molecule[i].Type = mtype;
          MoleculeType[mtype].Number = 1;
          // copy bead sequence and determine BTypes stuff
          MoleculeType[mtype].nBeads = mols[i];
          MoleculeType[mtype].Bead = calloc(MoleculeType[mtype].nBeads, sizeof(int));
          MoleculeType[mtype].nBTypes = 0;
          MoleculeType[mtype].BType = calloc(1, sizeof(int));
          for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
            int btype = Bead[Molecule[i].Bead[j]].Type;
            MoleculeType[mtype].Bead[j] = btype;
            exists = false; // recycling the bool to check if btype is already in BType[]
            for (int k = 0; k < MoleculeType[mtype].nBTypes; k++) {
              if (MoleculeType[mtype].BType[k] == btype) {
                exists = true;
              }
            }
            if (!exists) { // recycled exists
              int types = MoleculeType[mtype].nBTypes;
              MoleculeType[mtype].nBTypes++;
              MoleculeType[mtype].BType = realloc(MoleculeType[mtype].BType, (types+1)*sizeof(int));
              MoleculeType[mtype].BType[types] = btype;
            }
          }
          // copy bonds
          MoleculeType[mtype].nBonds = bonds_per_mol[i];
          MoleculeType[mtype].Bond = calloc(MoleculeType[mtype].nBonds, sizeof(int *));
          for (int j = 0; j < MoleculeType[mtype].nBonds; j++) {
            MoleculeType[mtype].Bond[j] = calloc(3, sizeof(int));
            MoleculeType[mtype].Bond[j][0] = mol_bonds[i][j][0];
            MoleculeType[mtype].Bond[j][1] = mol_bonds[i][j][1];
            MoleculeType[mtype].Bond[j][2] = mol_bonds[i][j][2];
          }
          MoleculeType[mtype].nAngles = 0;
          MoleculeType[mtype].Angle = calloc(1, sizeof(int *));
          MoleculeType[mtype].Write = true;
          Counts.TypesOfMolecules++;
        } //}}}
      } //}}}
      // free helper arrays //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          free(mol_bonds[i][j]);
        }
        free(mol_bonds[i]);
      }
      free(mol_bonds);
      free(bonds_per_mol);
      free(diff_conn); //}}}
    } //}}}
    // angles section //{{{
    if (words > 0 && strcmp(split[0], "Angles") == 0) {
      // allocate helper arrays to hold angle info //{{{
      // number of angles in each molecule
      int *angles_per_mol = calloc(Counts.Molecules, sizeof(int));
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] & [i][][2] ... ids of connected beads in molecule 'i'
      // [i][][3] ... angle type
      int ***mol_angles = calloc(Counts.Molecules, sizeof(int **));
      for (int i = 0; i < Counts.Molecules; i++) {
         mol_angles[i] = calloc(1, sizeof(int *));
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // read all angles //{{{
      for (int i = 0; i < angles; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <angle id> <angle type> <bead1> <bead2> <bead3>
        // Error - incorrect format //{{{
        if (words < 5 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3]) || !IsInteger(split[4])) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - each 'Angles' line must be <angle id> <angle type> <bead1d> <bead2> <bead3>\n", input);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        int bead3 = atoi(split[4]) - 1;
        // assign molecule to the bond
        int mol = Bead[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != Bead[bead2].Molecule || mol != Bead[bead3].Molecule) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m ", input);
          fprintf(stderr, " - atoms in \033[1;33mangle %d\033[1;31m are in different molecules\n\n", atoi(split[0]));
          fprintf(stderr, "\033[0m");
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        angles_per_mol[mol]++;
        int angle = angles_per_mol[mol];
        // add angle beads to the molecule they belong to
        mol_angles[mol] = realloc(mol_angles[mol], angle*sizeof(int *));
        mol_angles[mol][angle-1] = calloc(4, sizeof(int));
        mol_angles[mol][angle-1][0] = bead1;
        mol_angles[mol][angle-1][1] = bead2;
        mol_angles[mol][angle-1][2] = bead3;
        mol_angles[mol][angle-1][3] = type;
      } //}}}
      // sort angles according to the id of the first bead in a bond //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        SortAngles(mol_angles[i], angles_per_mol[i]);
      } //}}}
      // minimize mol_angles based on lowest id in each molecule //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        int lowest = Counts.Beads; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if (Molecule[i].Bead[j] < lowest) {
            lowest = Molecule[i].Bead[j];
          }
        }
        for (int j = 0; j < angles_per_mol[i]; j++) {
          mol_angles[i][j][0] -= lowest;
          mol_angles[i][j][1] -= lowest;
          mol_angles[i][j][2] -= lowest;
        }
      } //}}}
      // add angles to existing molecule types (or create a new one differing only in angles)
      // ...'basic' molecule types already generated in bonds section
      int *extra = calloc(Counts.TypesOfMolecules, sizeof(int)); // number of molecule types differing in angles
      for (int i = 0; i < Counts.Molecules; i++) {
        int mtype = Molecule[i].Type;
        if (MoleculeType[mtype].nAngles == 0) { // add angles if there are no angles in the molecule type //{{{
          MoleculeType[mtype].nAngles = angles_per_mol[i];
          MoleculeType[mtype].Angle = realloc(MoleculeType[mtype].Angle, MoleculeType[mtype].nAngles*sizeof(int *));
          for (int j = 0; j < MoleculeType[mtype].nAngles; j++) {
            MoleculeType[mtype].Angle[j] = calloc(4, sizeof(int));
            MoleculeType[mtype].Angle[j][0] = mol_angles[i][j][0];
            MoleculeType[mtype].Angle[j][1] = mol_angles[i][j][1];
            MoleculeType[mtype].Angle[j][2] = mol_angles[i][j][2];
            MoleculeType[mtype].Angle[j][3] = mol_angles[i][j][3];
          } //}}}
        } else { // if angles already present, check whether they're the same //{{{
          bool exists = true;
          // same number of angles? //{{{
          if (MoleculeType[mtype].nAngles != angles_per_mol[i]) {
            exists = false;
          } //}}}
          // same angles as in mtype? //{{{
          for (int j = 0; j < angles_per_mol[i] && j < MoleculeType[mtype].nAngles; j++) {
            if (MoleculeType[mtype].Angle[j][0] != mol_angles[i][j][0] ||
                MoleculeType[mtype].Angle[j][1] != mol_angles[i][j][1] ||
                MoleculeType[mtype].Angle[j][2] != mol_angles[i][j][2] ||
                MoleculeType[mtype].Angle[j][3] != mol_angles[i][j][3] ) {
              exists = false;
            }
          } //}}}
          // check against angle-generated molecule types //{{{
          // If its angles aren't the same as in mtype, check against other
          // types (i.e., against newly generated thanks to different angles)
          if (!exists) {
            for (int j = (mtype+1); j < Counts.TypesOfMolecules; j++) {
              if (MoleculeType[j].nAngles == angles_per_mol[i]) {
                int count = 0;
                for (int k = 0; k < MoleculeType[j].nAngles; k++) {
                  if (MoleculeType[j].Angle[k][0] == mol_angles[i][k][0] &&
                      MoleculeType[j].Angle[k][1] == mol_angles[i][k][1] &&
                      MoleculeType[j].Angle[k][2] == mol_angles[i][k][2] &&
                      MoleculeType[j].Angle[k][3] == mol_angles[i][k][3] ) {
                    count++;
                  }
                }
                if (count == angles_per_mol[i]) {
                  exists = true;
                  MoleculeType[j].Number++;
                  MoleculeType[mtype].Number--;
                  Molecule[i].Type = j;
                  break;
                }
              }
            }
          } //}}}
          // create a new molecule type //{{{
          if (!exists) {
            int new = Counts.TypesOfMolecules;
            extra[mtype]++;
            Molecule[i].Type = new;
            Counts.TypesOfMolecules++;
            MoleculeType = realloc(MoleculeType, (new+1)*sizeof(struct MoleculeType));
            char name[19];
            sprintf(name, "%s-a%d", MoleculeType[mtype].Name, extra[mtype]);
            strcpy(MoleculeType[new].Name, name);
            MoleculeType[new].Number = 1;
            MoleculeType[mtype].Number--;
            MoleculeType[new].nBeads = MoleculeType[mtype].nBeads;
            MoleculeType[new].Bead = calloc(MoleculeType[new].nBeads, sizeof(int));
            for (int j = 0; j < MoleculeType[new].nBeads; j++) {
              MoleculeType[new].Bead[j] = MoleculeType[mtype].Bead[j];
            }
            MoleculeType[new].nBonds = MoleculeType[mtype].nBonds;
            MoleculeType[new].Bond = calloc(MoleculeType[new].nBonds, sizeof(int *));
            for (int j = 0; j < MoleculeType[new].nBonds; j++) {
              MoleculeType[new].Bond[j] = calloc(3, sizeof(int));
              MoleculeType[new].Bond[j][0] = MoleculeType[mtype].Bond[j][0];
              MoleculeType[new].Bond[j][1] = MoleculeType[mtype].Bond[j][1];
              MoleculeType[new].Bond[j][2] = MoleculeType[mtype].Bond[j][2];
            }
            MoleculeType[new].nAngles = MoleculeType[mtype].nAngles;
            MoleculeType[new].Angle = calloc(MoleculeType[new].nAngles, sizeof(int *));
            for (int j = 0; j < MoleculeType[new].nAngles; j++) {
              MoleculeType[new].Angle[j] = calloc(4, sizeof(int));
              MoleculeType[new].Angle[j][0] = mol_angles[i][j][0];
              MoleculeType[new].Angle[j][1] = mol_angles[i][j][1];
              MoleculeType[new].Angle[j][2] = mol_angles[i][j][2];
              MoleculeType[new].Angle[j][3] = mol_angles[i][j][3];
            }
            MoleculeType[new].nBTypes = MoleculeType[mtype].nBTypes;
            MoleculeType[new].BType = calloc(MoleculeType[new].nBTypes, sizeof(int));
            for (int j = 0; j < MoleculeType[new].nBTypes; j++) {
              MoleculeType[new].BType[j] = MoleculeType[mtype].BType[j];
            }
            MoleculeType[new].Mass = MoleculeType[mtype].Mass;
            MoleculeType[new].Write = true;
          } //}}}
        } //}}}
      }
      // free helper arrays //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        for (int j = 0; j < angles_per_mol[i]; j++) {
          free(mol_angles[i][j]);
        }
        free(mol_angles[i]);
      }
      free(mol_angles);
      free(angles_per_mol);
      free(extra); //}}}
    } //}}}
    // read and split next line
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
  }
  free(mols);
  fclose(fr); //}}}

  // calculate molecular mass //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Mass = 0;
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      MoleculeType[i].Mass += BeadType[MoleculeType[i].Bead[j]].Mass;
    }
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput("\0", Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    PrintBondTypes(Counts, bond_type);
    PrintAngleTypes(Counts, angle_type);
  } //}}}

  // write out.vcf file //{{{
  FILE *fw;
  if ((fw = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  }

  // print command to output vcf file
  fprintf(fw, "# Created by: lmp_vtf");
  for (int i = 1; i < argc; i++) {
    fprintf(fw, " %s", argv[i]);
  }
  fprintf(fw, " (AnalysisTools version %s; https://github.com/KaGaSi/AnalysisTools/releases)\n", VERSION);

  fprintf(fw, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  char stuff[LINE];
  strcpy(stuff, "\0");
  WriteCoorIndexed(fw, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

  fclose(fw); //}}}

  // create & fill output vsf file
  WriteVsf(output_vsf, Counts, BeadType, Bead, MoleculeType, Molecule, false);

  // free memory (to make valgrind happy) //{{{
  free(BeadType);
  FreeBead(Counts, &Bead);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  free(Index);
  free(angle_type);
  free(bond_type); //}}}
}
