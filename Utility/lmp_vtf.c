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

// TODO: error checking for input file

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
  Counts Counts = ZeroCounts; // structure with number of beads, molecules, etc.
  Molecule *Molecule;
  MoleculeType *MoleculeType;
  Vector BoxLength;
  Vector box_lo; // {x,y,z}lo from data file to place beads in (0, BoxLength>
  int bonds = 0; // total number of bonds
//int angles = 0; // total number of angles //}}}

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
  line[0] = '\0';
  char split[30][100], delim[8];
  strcpy(delim, " \t");
  // TODO: redo with do while - check split[0], not line[]
  while (line[0] == '\0' || // empty line
         line[0] == '#' || // comment line
         line[0] == '-' ||  // negative number
         line[0] == '\n' ||  // just white space
         (line[0] >= '0' && line[0] <= '9')) { // positive number
    // read one line
    fgets(line, sizeof(line), fr);
    int words = SplitLine(split, line, delim);

    // read header data //{{{
    for (int i = 0; i < words; i++) {
      if (strcmp(split[i], "atoms") == 0) {
        Counts.BeadsInVsf = atoi(split[0]);
        Counts.Beads = atoi(split[0]);
      }
      if (strcmp(split[i], "bonds") == 0) {
        bonds = atoi(split[0]);
      }
//    if (strcmp(split[i], "angles") == 0) {
//      angles = atoi(split[0]);
//    }
      if (strcmp(split[i], "atom") == 0 &&
          (i+1) < words && strcmp(split[i+1], "types") == 0) {
        Counts.TypesOfBeads = atoi(split[0]);
      }
      if (strcmp(split[i], "xlo") == 0 &&
          (i+1) < words && strcmp(split[i+1], "xhi") == 0) {
        BoxLength.x = atof(split[1]) - atof(split[0]);
        box_lo.x = atof(split[0]);
      }
      if (strcmp(split[i], "ylo") == 0 &&
          (i+1) < words && strcmp(split[i+1], "yhi") == 0) {
        BoxLength.y = atof(split[1]) - atof(split[0]);
        box_lo.y = atof(split[0]);
      }
      if (strcmp(split[i], "zlo") == 0 &&
          (i+1) < words && strcmp(split[i+1], "zhi") == 0) {
        BoxLength.z = atof(split[1]) - atof(split[0]);
        box_lo.z = atof(split[0]);
      }
    } //}}}
  } //}}}

  // some error checking //{{{
  if (Counts.TypesOfBeads == 0) {
    fprintf(stderr, "\nError - missing 'atom types' line (or is 0) in %s\n\n", input);
    exit(1);
  }
  if (Counts.BeadsInVsf == 0) {
    fprintf(stderr, "\nError - missing 'atoms' line (or is 0) in %s\n\n", input);
    exit(1);
  }
  if (BoxLength.x == 0) {
    fprintf(stderr, "\nError - missing 'xlo xhi' line (or is 0 0) in %s\n\n", input);
    exit(1);
  }
  if (BoxLength.y == 0) {
    fprintf(stderr, "\nError - missing 'ylo yhi' line (or is 0 0) in %s\n\n", input);
    exit(1);
  }
  if (BoxLength.z == 0) {
    fprintf(stderr, "\nError - missing 'zlo zhi' line (or is 0 0) in %s\n\n", input);
    exit(1);
  } //}}}

  // fill something in BeadType struct //{{{
  BeadType *BeadType = calloc(Counts.TypesOfBeads, sizeof(struct BeadType));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    sprintf(BeadType[i].Name, "bead%d", i+1);
    BeadType[i].Use = true;
    BeadType[i].Write = true;
  } //}}}

  // bead struct memory allocation //{{{
  Bead *Bead = calloc(Counts.Beads, sizeof(struct Bead));
  for (int i = 0; i < Counts.Beads; i++) {
    Bead[i].Aggregate = calloc(1, sizeof(int));
  } //}}}

  int *Index = calloc(Counts.Beads, sizeof(int)); // link between indices in vsf and in program (i.e., opposite of Bead[].Index)

  // read body of data file //{{{
  int test,
      *mols, // number of beads in each molecule
      monomer = 0; // monomer beads designated by mol_ID = 0 in data file
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);

    // atom masses //{{{
    if (strcmp(split[0], "Masses") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      for (int i = 0; i < Counts.TypesOfBeads; i++) {
        fgets(line, sizeof(line), fr);
        strcpy(line, TrimLine(line)); // trim excess whitespace
        int words = SplitLine(split, line, delim);
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

    // atoms section //{{{
    if (strcmp(split[0], "Atoms") == 0) {
      mols = calloc(Counts.Beads, sizeof(int));
      // array for list of beads in each molecule
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
        SplitLine(split, line, delim);

        int id = atoi(split[0]) - 1, // in lammps, these start at 1
            mol_id = atoi(split[1]) - 1, // in lammps, molecules start with 1; unbonded atoms can be 0
            btype = atoi(split[2]) - 1; // in lammps, these start at 1

        if (mol_id == -1) { // corresponds to 0 in the data file
          monomer++;
          Bead[id].Molecule = -1;
        } else {
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
      }
      // zeroize unused part of mols array - just to be on the save side
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
    if (strcmp(split[0], "Bonds") == 0) {
      // allocate helper arrays to hold bond info //{{{
      // number of bonds in each molecule
      int *bonds_per_mol = calloc(Counts.Molecules, sizeof(int));
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] ... ids of connected beads in molecule 'i'
      int ***mol_bonds = calloc(Counts.Molecules, sizeof(int **));
      for (int i = 0; i < Counts.Molecules; i++) {
         mol_bonds[i] = calloc(1, sizeof(int *));
      } //}}}

      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);

      // read all bonds //{{{
      for (int i = 0; i < bonds; i++) {
        fgets(line, sizeof(line), fr);
        SplitLine(split, line, delim);

        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;

        // assign molecule to the bond
        int mol = Bead[bead1].Molecule;

        // error when the second bead is in different molecule //{{{
        if (mol != Bead[bead2].Molecule) {
          fprintf(stderr, "\nError in bond #%d in %s: ", atoi(split[0]), input);
          fprintf(stderr, "atoms %d and %d are in different molecules ", bead1+1, bead2+1);
          fprintf(stderr, "(%d and %d)\n\n", mol+1, Bead[bead2].Molecule+1);
          exit(1);
        } //}}}

        // increment number of bonds in the molecule
        bonds_per_mol[mol]++;
        int bond = bonds_per_mol[mol];
        // add bonded beads to the molecule they belong to
        mol_bonds[mol] = realloc(mol_bonds[mol], bond*sizeof(int *));
        mol_bonds[mol][bond-1] = calloc(2, sizeof(int));
        mol_bonds[mol][bond-1][0] = bead1;
        mol_bonds[mol][bond-1][1] = bead2;
      } //}}}

      // sort bonds according to the id of the first bead in a bond //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        SortBonds(mol_bonds[i], bonds_per_mol[i]);
      } //}}}

      // minimize mol_bonds based on lowest id in each molecule //{{{
      for (int i = 0; i < Counts.Molecules; i++) {
        int lowest = 1000000;
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
      for (int i = 0; i < Counts.Molecules; i++) {
        bool exists = false; // is molecule 'i' of known type?
        for (int j = 0; j < Counts.TypesOfMolecules; j++) { //{{{
          if (MoleculeType[j].nBeads == mols[i] && // same number of molecules?
              MoleculeType[j].nBonds == bonds_per_mol[i]) { // same number of bonds?
            // same connectivity?
            bool same_conn = true;
            for (int k = 0; k < MoleculeType[j].nBonds; k++) {
              if (MoleculeType[j].Bond[k][0] != mol_bonds[i][k][0] ||
                  MoleculeType[j].Bond[k][1] != mol_bonds[i][k][1]) {
                same_conn = false;
                break;
              }
            }
            // same bead types?
            bool same_bead = true;
            for (int k = 0; k < MoleculeType[j].nBeads; k++) {
              int btype = Bead[Molecule[i].Bead[k]].Type;
              if (MoleculeType[j].Bead[k] != btype) {
                same_bead = false;
                break;
              }
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
          // molecule name is 'mol<number>'
          sprintf(MoleculeType[mtype].Name, "mol%d", mtype+1);
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
            exists = false; // recycling the bool
            for (int k = 0; k < MoleculeType[mtype].nBTypes; k++) {
              if (MoleculeType[mtype].BType[k] == btype) {
                exists = true;
              }
            }
            if (!exists) { // recycled exists
              int types = MoleculeType[mtype].nBTypes;
              MoleculeType[mtype].BType = realloc(MoleculeType[mtype].BType, (types+1)*sizeof(int));
              MoleculeType[mtype].BType[types] = btype;
              MoleculeType[mtype].nBTypes++;
            }
          }
          // copy bonds
          MoleculeType[mtype].nBonds = bonds_per_mol[i];
          MoleculeType[mtype].Bond = calloc(MoleculeType[mtype].nBonds, sizeof(int *));
          for (int j = 0; j < MoleculeType[mtype].nBonds; j++) {
            MoleculeType[mtype].Bond[j] = calloc(2, sizeof(int));
            MoleculeType[mtype].Bond[j][0] = mol_bonds[i][j][0];
            MoleculeType[mtype].Bond[j][1] = mol_bonds[i][j][1];
          }
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
      free(bonds_per_mol); //}}}
    } //}}}

//  // ignore angles //{{{
//  if (strcmp(split[0], "Angles") == 0) {
//    fgets(line, sizeof(line), fr);
//    for (int i = 0; i < angles; i++) {
//      fgets(line, sizeof(line), fr);
//    }
//  } //}}}

    // read and split next line
    fgets(line, sizeof(line), fr);
    SplitLine(split, line, delim);
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
  } //}}}

  // write out.vcf file //{{{
  FILE *fw;
  if ((fw = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  }

  // print command
  fprintf(fw, "#Created by: lmp_vtf");
    for (int i = 1; i < argc; i++) {
      fprintf(fw, " %s", argv[i]);
    }
  fprintf(fw, "\n");

  fprintf(fw, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  char stuff[LINE];
  strcpy(stuff, "\0");
  WriteCoorIndexed(fw, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

  fclose(fw); //}}}

  // create & fill output vsf file
  WriteVsf(output_vsf, Counts, BeadType, Bead, MoleculeType, Molecule);

  // free memory (to make valgrind happy) //{{{
  free(BeadType);
  FreeBead(Counts, &Bead);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  free(Index); //}}}
}
