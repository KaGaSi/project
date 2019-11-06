#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
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
lmp_vsf takes lammps data file and transforms it into vcf and vsf files. \n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <out.vsf> <out.vcf> <options>\n\n", cmd);

//fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
//fprintf(ptr, "   <out.data>        output lammps data file\n");
//fprintf(ptr, "   <options>\n");
//fprintf(ptr, "      -f <name>      FIELD file (default: FIELD)\n");
//fprintf(ptr, "      --srp          add one more bead type for srp");
//fprintf(ptr, "      -st <step>     timestep for creating CONFIG (default: last)\n");
//CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
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

//// test if options are given correctly //{{{
//for (int i = 1; i < argc; i++) {
//  if (argv[i][0] == '-' &&
//      strcmp(argv[i], "-i") != 0 &&
//      strcmp(argv[i], "-v") != 0 &&
//      strcmp(argv[i], "-V") != 0 &&
//      strcmp(argv[i], "-s") != 0 &&
//      strcmp(argv[i], "-h") != 0 &&
//      strcmp(argv[i], "--script") != 0 &&
//      strcmp(argv[i], "--srp") != 0 &&
//      strcmp(argv[i], "-f") != 0 &&
//      strcmp(argv[i], "-st") != 0) {

//    ErrorOption(argv[i]);
//    Help(argv[0], true);
//    exit(1);
//  }
//} //}}}

  // print command to stdout //{{{
  bool silent = false;
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input[1024];
  strcpy(input, argv[++count]); //}}}

  // <out.vsf> - output structure file //{{{
  char output_vsf[1024];
  strcpy(output_vsf, argv[++count]);

  // test if <out.vsf> filename ends with '.vsf' (required by VMD)
  int ext = 1;
  char extension[1][5];
  strcpy(extension[0], ".vsf");
  if (!ErrorExtension(output_vsf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output structure file //{{{
  char output_vcf[1024];
  strcpy(output_vcf, argv[++count]);

  // test if <out.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (!ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // some variables //{{{
  Counts Counts; // structure with number of beads, molecules, etc.
  Molecule *Molecule;
  MoleculeType *MoleculeType;
  Vector BoxLength;
  Vector box_lo; // {x,y,z}lo from data file to place beads in (0, BoxLength>
  ZeroCounts(&Counts);
  int bonds = 0; // total number of bonds
//int bond_types = 0; // number of bond types
  int angles = 0; // total number of angles //}}}

  // open <input> //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // ignore first line (comment) //{{{
  char line[1024];
  fgets(line, sizeof(line), fr); //}}}

  // data file header lines must start with a number (or '#' for comment), //{{{
  // therefore read until something else is encountered
  line[0] = '\0';
  char *split[30]; // to hold individual strings from the line
  while (line[0] == '\0' || // empty line
         line[0] == '#' || // comment line
         line[0] == '-' ||  // negative number
         (line[0] >= '0' && line[0] <= '9')) { // positive number
    // read one line
    fgets(line, sizeof(line), fr);
    // trim whitespace in line //{{{
    // 1) trailing whitespace
    int length = strlen(line);
    while (length >= 1 &&
           (line[length-1] == ' ' ||
            line[length-1] == '\n' ||
            line[length-1] == '\t')) {
      // last string character needs to be '\0'
      line[length-1] = '\0';
      length--;
    }
    // 2) preceding whitespace
    while (length > 1 &&
           (line[0] == ' ' ||
            line[0] == '\n' ||
            line[0] == '\t')) {
      for (int i = 0; i < length; i++) { // line[length] contains '\0'
        line[i] = line[i+1];
      }
      length--;
    } //}}}

    count = 0; // number of strings in a line (except for trailing comment)
    // split line //{{{
    split[0] = strtok(line, " \t");
    while (split[count] != NULL && count < 29 && split[count][0] != '#') {
      split[++count] = strtok(NULL, " \t");
    } //}}}

    // read header data //{{{
    for (int i = 0; i < count; i++) {
      if (strcmp(split[i], "atoms") == 0) {
        Counts.BeadsInVsf = atoi(split[0]);
        Counts.Beads = atoi(split[0]);
      }
      if (strcmp(split[i], "bonds") == 0) {
        bonds = atoi(split[0]);
      }
      if (strcmp(split[i], "angles") == 0) {
        angles = atoi(split[0]);
      }
      if (strcmp(split[i], "atom") == 0 &&
          (i+1) < count && strcmp(split[i+1], "types") == 0) {
        Counts.TypesOfBeads = atoi(split[0]);
      }
//    if (strcmp(split[i], "bond") == 0 &&
//        (i+1) < count && strcmp(split[i+1], "types") == 0) {
//      bond_types = atoi(split[0]);
//    }
      if (strcmp(split[i], "xlo") == 0 &&
          (i+1) < count && strcmp(split[i+1], "xhi") == 0) {
        BoxLength.x = atof(split[1]) - atof(split[0]);
        box_lo.x = atof(split[0]);
      }
      if (strcmp(split[i], "ylo") == 0 &&
          (i+1) < count && strcmp(split[i+1], "yhi") == 0) {
        BoxLength.y = atof(split[1]) - atof(split[0]);
        box_lo.y = atof(split[0]);
      }
      if (strcmp(split[i], "zlo") == 0 &&
          (i+1) < count && strcmp(split[i+1], "zhi") == 0) {
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

  // read body of data file //{{{
  int test;
  int *mols; // to determine number of beads in each molecule
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);

    // atom masses //{{{
    if (strcmp(split[0], "Masses") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      for (int i = 0; i < Counts.TypesOfBeads; i++) {
        fgets(line, sizeof(line), fr);
        // trim whitespace in line //{{{
        // 1) trailing whitespace
        int length = strlen(line);
        while (length >= 1 &&
               (line[length-1] == ' ' ||
                line[length-1] == '\n' ||
                line[length-1] == '\t')) {
          // last string character needs to be '\0'
          line[length-1] = '\0';
          length--;
        }
        // 2) preceding whitespace
        while (length > 1 &&
               (line[0] == ' ' ||
                line[0] == '\n' ||
                line[0] == '\t')) {
          for (int i = 0; i < length; i++) { // line[length] contains '\0'
            line[i] = line[i+1];
          }
          length--;
        } //}}}
        // split line //{{{
        split[0] = strtok(line, " \t");
        count = 0; // number of strings in a line (except for trailing comment)
        while (split[count] != NULL && count < 29) {
          split[++count] = strtok(NULL, " \t");
        } //}}}
        BeadType[i].Mass = atof(split[1]);
        // if there's a comment at the end of the line, consider it bead name
        if (split[2][0] == '#') {
          if (strlen(split[2]) == 1 && count > 3) { // comment form '# name'
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
      // 0 for unbonded; 1..N for bonded/unbonded (depends on the data writer's style)
      mols = calloc((Counts.Beads+1), sizeof(int));

      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);

      // go through atom section to get basic info //{{{
      fpos_t pos; // set file counter
      fgetpos(fr, &pos); // save file pointer
      for (int i = 0; i < Counts.Beads; i++) {
        fgets(line, sizeof(line), fr);
        // trim whitespace in line //{{{
        // 1) trailing whitespace
        int length = strlen(line);
        while (length >= 1 &&
               (line[length-1] == ' ' ||
                line[length-1] == '\n' ||
                line[length-1] == '\t')) {
          // last string character needs to be '\0'
          line[length-1] = '\0';
          length--;
        }
        // 2) preceding whitespace
        while (length > 1 &&
               (line[0] == ' ' ||
                line[0] == '\n' ||
                line[0] == '\t')) {
          for (int i = 0; i < length; i++) { // line[length] contains '\0'
            line[i] = line[i+1];
          }
          length--;
        } //}}}
        // split line //{{{
        split[0] = strtok(line, " \t");
        count = 0; // number of strings in a line (except for trailing comment)
        while (split[count] != NULL && count < 29 && split[count][0] != '#') {
          split[++count] = strtok(NULL, " \t");
        } //}}}

        int id = atoi(split[0]) - 1, // in lammps, these start at 1
            mol_id = atoi(split[1]), // in lammps, molecules start with 1, unbonded atoms can be 0
            btype = atoi(split[2]) - 1; // in lammps, these start at 1

        mols[mol_id]++;
        BeadType[btype].Charge = atof(split[3]);
        BeadType[btype].Number++;
        Bead[id].Position.x = atof(split[4]) - box_lo.x;
        Bead[id].Position.y = atof(split[5]) - box_lo.y;
        Bead[id].Position.z = atof(split[6]) - box_lo.z;
        Bead[id].Type = btype;
        Bead[id].Index = id;
      } //}}}

      // count numbers of (un)bonded beads and molecules //{{{
      int max_beads = 0; // largest molecule
      // mols[0] contains beads with molecule id 0, i.e., those that lammps considers unbonded atoms
      Counts.Unbonded = mols[0];
      for (int i = 0; i < mols[0]; i++) {
        Bead[i].Molecule = -1;
      }
      // go through atoms with molecule id > 0 (can be unbonded if there's only 1 atom in the molecule)
      for (int i = 1; i <= Counts.Beads; i++) {
        if (mols[i] != 0) {
          if (mols[i] == 1) {
            Counts.Unbonded++;
            Bead[i-1].Molecule = -1;
          } else {
            // assign molecule id to all beads in molecule 'i'
            int bead_count = Counts.Unbonded + Counts.Bonded;
            for (int j = bead_count; j < (bead_count+mols[i]); j++) {
              Bead[j].Molecule = Counts.Molecules;
            }
            Counts.Molecules++;
            Counts.Bonded += mols[i];
            if (mols[i] > max_beads) {
              max_beads = mols[i];
            }
          }
        }
      }

      // make mols[] array continuous array of only molecules
      count = 0; // molecule id
      for (int i = 0; i < Counts.Beads && count < Counts.Molecules; i++) {
        if (mols[i] > 1) {
          mols[count] = mols[i];
          count++;
        }
      } //}}}

      // allocate Molecule struct //{{{
      Molecule = calloc(Counts.Molecules, sizeof(struct Molecule));
      for (int i = 0; i < Counts.Molecules; i++) {
        Molecule[i].Bead = calloc(max_beads, sizeof(int));
        mols[i] = 0;
      } //}}}

      // go through atom section again to get bead type sequence of each molecule //{{{
      fsetpos(fr, &pos); // restore file pointer
      for (int i = 0; i < Counts.Beads; i++) {
        fgets(line, sizeof(line), fr);
        // trim whitespace in line //{{{
        // 1) trailing whitespace
        int length = strlen(line);
        while (length >= 1 &&
               (line[length-1] == ' ' ||
                line[length-1] == '\n' ||
                line[length-1] == '\t')) {
          // last string character needs to be '\0'
          line[length-1] = '\0';
          length--;
        }
        // 2) preceding whitespace
        while (length > 1 &&
               (line[0] == ' ' ||
                line[0] == '\n' ||
                line[0] == '\t')) {
          for (int i = 0; i < length; i++) { // line[length] contains '\0'
            line[i] = line[i+1];
          }
          length--;
        } //}}}
        // split line //{{{
        split[0] = strtok(line, " \t");
        count = 0; // number of strings in a line (except for trailing comment)
        while (split[count] != NULL && count < 29 && split[count][0] != '#') {
          split[++count] = strtok(NULL, " \t");
        } //}}}

        int id = atoi(split[0]) - 1, // in lammps, these start at 1
            btype = atoi(split[2]) - 1; // in lammps, these start at 1

        if (btype)
          ;
        int mol_id = Bead[id].Molecule;
        if (mol_id != -1) {
          Molecule[mol_id].Bead[mols[mol_id]] = id;
          mols[mol_id]++;
        }
      } //}}}

      // assume that the sequence of beads is the distinguishing feature of a molecule type //{{{
      // TODO: use bonds to discern proper molecules
      Counts.TypesOfMolecules = 0;
      MoleculeType = calloc(1, sizeof(struct MoleculeType));

      for (int i = 0; i < Counts.Molecules; i++) {
        bool exists = true;
        if (i == 0) { // for first molecule, there is no type yet
          exists = false;
        } else {
          Molecule[i].Type = -1;
          for (int j = 0; j < Counts.TypesOfMolecules; j++) {
            bool this_type = true;
            for (int k = 0; k < mols[i]; k++) {
              if (Bead[Molecule[i].Bead[k]].Type != MoleculeType[j].Bead[k]) {
                this_type = false;
                break;
              }
            }
            if (this_type) {
              Molecule[i].Type = j;
              MoleculeType[j].Number++;
            }
          }
          if (Molecule[i].Type == -1) {
            exists = false;
          }
        }

        // add new molecule type //{{{
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
          MoleculeType[mtype].BType = calloc(Counts.TypesOfBeads, sizeof(int));
          for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
            int btype = Bead[Molecule[i].Bead[j]].Type;
            MoleculeType[mtype].Bead[j] = btype;
            exists = false;
            for (int k = 0; k < MoleculeType[mtype].nBTypes; k++) {
              if (MoleculeType[mtype].BType[k] == btype) {
                exists = true;
              }
            }
            if (!exists) { // different exists than than the one 20 or so lines up
              MoleculeType[mtype].BType[MoleculeType[mtype].nBTypes] = btype;
              MoleculeType[mtype].nBTypes++;
            }
          }
          MoleculeType[mtype].Write = true;
          Counts.TypesOfMolecules++;
        } //}}}
      } //}}}
    } //}}}

    // bonds //{{{
    if (strcmp(split[0], "Bonds") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      for (int i = 0; i < bonds; i++) {
        fgets(line, sizeof(line), fr);
        // trim whitespace in line //{{{
        // 1) trailing whitespace
        int length = strlen(line);
        while (length >= 1 &&
               (line[length-1] == ' ' ||
                line[length-1] == '\n' ||
                line[length-1] == '\t')) {
          // last string character needs to be '\0'
          line[length-1] = '\0';
          length--;
        }
        // 2) preceding whitespace
        while (length > 1 &&
               (line[0] == ' ' ||
                line[0] == '\n' ||
                line[0] == '\t')) {
          for (int i = 0; i < length; i++) { // line[length] contains '\0'
            line[i] = line[i+1];
          }
          length--;
        } //}}}
        // split line //{{{
        split[0] = strtok(line, " \t");
        count = 0; // number of strings in a line (except for trailing comment)
        while (split[count] != NULL && count < 29 && split[count][0] != '#') {
          split[++count] = strtok(NULL, " \t");
        } //}}}
      }
    } //}}}

    // ignore angles //{{{
    if (strcmp(split[0], "Angles") == 0) {
      fgets(line, sizeof(line), fr);
      for (int i = 0; i < angles; i++) {
        fgets(line, sizeof(line), fr);
      }
    } //}}}

    // read and split next line //{{{
    fgets(line, sizeof(line), fr);
    // trim whitespace in line //{{{
    // 1) trailing whitespace
    int length = strlen(line);
    while (length >= 1 &&
           (line[length-1] == ' ' ||
            line[length-1] == '\n' ||
            line[length-1] == '\t')) {
      // last string character needs to be '\0'
      line[length-1] = '\0';
      length--;
    }
    // 2) preceding whitespace
    while (length > 1 &&
           (line[0] == ' ' ||
            line[0] == '\n' ||
            line[0] == '\t')) {
      for (int i = 0; i < length; i++) { // line[length] contains '\0'
        line[i] = line[i+1];
      }
      length--;
    } //}}}
    char *split[30];
    count = 0; // number of strings in a line (except for trailing comment)
    split[0] = strtok(line, " \t");
    while (split[count] != NULL && count < 29 && split[count][0] != '#') {
      split[++count] = strtok(NULL, " \t");
    } //}}}
  }
  fclose(fr); //}}}
  free(mols);

  // bonds - for now, assume linear molecules with bead order according to data file //{{{
  // TODO: generalize by readinng Bonds sections
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].nBonds = MoleculeType[i].nBeads - 1;
    MoleculeType[i].Bond = calloc(MoleculeType[i].nBonds, sizeof(int *));
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      MoleculeType[i].Bond[j] = calloc(2, sizeof(int));
      MoleculeType[i].Bond[j][0] = j;
      MoleculeType[i].Bond[j][1] = j + 1;
    }
  } //}}}

  // calculate molecular mass //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Mass = 0;
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      MoleculeType[i].Mass += BeadType[MoleculeType[i].Bead[j]].Mass;
    }
  } //}}}

  // print information - verbose output //{{{
  bool verbose = true;
  bool verbose2 = false;
  if (verbose) {
    fprintf(stdout, "Box = (%lf, %lf, %lf)\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
    VerboseOutput(verbose2, "\0", Counts, BeadType, Bead, MoleculeType, Molecule);
    fprintf(stdout, "\nMolecular prototypes:\n");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(stdout, "   %s\n", MoleculeType[i].Name);
      for (int j = 0; j < Counts.Molecules; j++) {
        if (Molecule[j].Type == i) {
          int id0 = Molecule[j].Bead[0];
          fprintf(stdout, "     0.000000 0.000000 0.000000\n");
          for (int k = 1; k < MoleculeType[i].nBeads; k++) {
            int id = Molecule[j].Bead[k];
            fprintf(stdout, "     %lf %lf %lf\n", Bead[id].Position.x-Bead[id0].Position.x,
                                                  Bead[id].Position.y-Bead[id0].Position.y,
                                                  Bead[id].Position.z-Bead[id0].Position.z);
          }
          putchar('\n');
          break;
        }
      }
    }
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

  char stuff[1024];
  strcpy(stuff, "\0");
  WriteCoorIndexed(fw, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

  fclose(fw); //}}}

  // create & fill output vsf file
  WriteVsf(output_vsf, Counts, BeadType, Bead, MoleculeType, Molecule);

  // free memory (to make valgrind happy) //{{{
  free(BeadType);
  FreeBead(Counts, &Bead);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule); //}}}
}
