#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
GenBrushWall reads information from a FIELD-like file and creates \
vsf structure file and generates coordinates for all beads. \
It creates a simulation box with molecules attached wall in xy planes \
(i.e., with z coordinates of the first bead of each molecule equal to \
0 or the box size in the z direction). The molecules are generated \
on a square grid with specified distance between anchoring points.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <out.vsf> <out.vcf> <options>\n\n", cmd);
  fprintf(ptr, "   <out.vsf>              output structure file (*.vsf)\n");
  fprintf(ptr, "   <out.vcf>              output coordinate file (*.vcf)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -s <float> <float>  spacing in x and y directions (default: 1 1)\n");
  fprintf(ptr, "      -f <name>           FIELD-like file (default: FIELD)\n");
  fprintf(ptr, "      -v                  verbose output\n");
  fprintf(ptr, "      -h                  print this help and exit\n");
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
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
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
  // '-s' option - spacing in x and y directions //{{{
  int test = 2;
  double spacing[2];
  spacing[0] = 1;
  spacing[1] = 1;
  if (MultiDoubleOption(argc, argv, "-s", &test, spacing)) {
    exit(1);
  }
  if (test != 2) {
    fprintf(stderr, "\nError: option '-s' requires two numeric arguments\n\n");
    exit(1);
  } //}}}

  // output verbosity //{{{
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  // }}}

  // FIELD-like file //{{{
  char *input = calloc(1024, sizeof(char *));
  if (FileOption(argc, argv, "-f", &input)) {
    exit(1);
  }
  if (input[0] == '\0') {
    strcpy(input, "FIELD");
  } //}}}
  //}}}

  count = 0; // count arguments

  // <out.vsf> - output structure file (must end with .vsf) //{{{
  char output[1024];
  strcpy(output, argv[++count]);

  // test if <output.vsf> filename ends with '.vsf' or '.vtf' (required by VMD)
  int ext = 1;
  char **extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vsf");
  if (!ErrorExtension(output, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <out.vcf> - output vcf file //{{{
  char *output_vcf = calloc(1024, sizeof(char));
  char stuff[1024];
  strcpy(output_vcf, argv[++count]);

  // test if outpuf_vcf has '.vcf' extension - required by vmd //{{{
  ext = 1;
  extension = malloc(ext*sizeof(char *));
  extension[0] = malloc(5*sizeof(char));
  strcpy(extension[0], ".vcf");
  if (!ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}
  //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information from FIELD-like //{{{
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // read box size //{{{
  char line[1024], *box[3];
  fgets(line, sizeof(line), fr);
  box[0] = strtok(line, " \t");
  box[1] = strtok(NULL, " \t");
  box[2] = strtok(NULL, " \t");
  Vector BoxLength;
  BoxLength.x = atof(box[0]);
  BoxLength.y = atof(box[1]);
  BoxLength.z = atof(box[2]); //}}}

  Counts.BeadsInVsf = 3 * BoxLength.z * BoxLength.y * BoxLength.z;
  Counts.Beads = Counts.BeadsInVsf;
  Bead = malloc(Counts.Beads*sizeof(struct Bead));
  // allocate Bead Aggregate array - needed only to free() //{{{
  for (int i = 0; i < Counts.Beads; i++) {
    Bead[i].Aggregate = calloc(1,sizeof(double));
    Bead[i].nAggregates = 0;
  } //}}}

  int mols[2]; // number of molecules in x ([0]) and y ([1]) directions
  mols[0] = round(BoxLength.x / spacing[0]);
  mols[1] = round(BoxLength.y / spacing[1]);
  printf("Grid of %d x %d mols\n", mols[0], mols[1]);

  Counts.Molecules = mols[0] * mols[1] * 2;
  Molecule = calloc(Counts.Molecules, sizeof(struct Molecule));

  // read number of bead types //{{{
  while(fgets(line, sizeof(line), fr)) {
    char *split;
    split = strtok(line, " \t ");
    if (strcmp(split, "species") == 0 ||
        strcmp(split, "Species") == 0 ||
        strcmp(split, "SPECIES") == 0 ) {
      // solvent isn't in the FIELD-like, but it must be in the resulting vsf
      Counts.TypesOfBeads = atoi(strtok(NULL, " \t")) + 1;
      break;
    }
  } //}}}

  BeadType = calloc(Counts.TypesOfBeads,sizeof(struct BeadType));

  // read info about bead types //{{{
  // number of beads is ignored, because GenBrush only generates the brush;
  // other beads/molecules are added using AddToSystem
  strcpy(BeadType[0].Name, "None"); // BeadType[0] is a placeholder for AddToSystem
  BeadType[0].Mass = 1;
  BeadType[0].Charge = 0;
  BeadType[0].Number = 0;
  for (int i = 1; i < Counts.TypesOfBeads; i++) {
    fgets(line, sizeof(line), fr);

    // split the line into array
    char *split[4];
    split[0] = strtok(line, " \t");
    for (int j = 1; j < 4; j++) {
      split[j] = strtok(NULL, " \t");
    }

    strcpy(BeadType[i].Name, split[0]);
    BeadType[i].Mass = atof(split[1]);
    BeadType[i].Charge = atof(split[2]);
    BeadType[i].Number = 0;
  } //}}}

  // TODO: generalize - for now, only the first molecule type is used
  // read number of molecule types //{{{
  while(fgets(line, sizeof(line), fr)) {
    char *split;
    split = strtok(line, " \t ");
    if (strncmp(split, "molecule", 8) == 0 ||
        strncmp(split, "Molecule", 8) == 0 ||
        strncmp(split, "MOLECULE", 8) == 0 ) {
      Counts.TypesOfMolecules = atoi(strtok(NULL, " \t"));
      break;
    }
  } //}}}

  MoleculeType = calloc(Counts.TypesOfMolecules,sizeof(struct MoleculeType));

  // read info about molecule types (and calculate numbers of beads and stuff) //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // name //{{{
    fgets(line, sizeof(line), fr);
    // trim trailing whitespace in line
    int length = strlen(line);
    // last string character needs to be '\0'
    while (length > 1 &&
           (line[length-1] == ' ' ||
            line[length-1] == '\n' ||
            line[length-1] == '\t')) {
      line[length-1] = '\0';
      length--;
    }
    strcpy(MoleculeType[i].Name, strtok(line, " \t")); //}}}
    // number of molecules - irrelevant, decided by spacing //{{{
    fgets(line, sizeof(line), fr);
    strtok(line, " \t ");
//  MoleculeType[i].Number = atoi(strtok(NULL, " \t")); //}}}
    // number of beads //{{{
    fgets(line, sizeof(line), fr);
    strtok(line, " \t ");
    MoleculeType[i].nBeads = atoi(strtok(NULL, " \t")); //}}}
    // number & ids of beadtypes & mass //{{{
    MoleculeType[i].nBTypes = 0;
    MoleculeType[i].BType = calloc(1,sizeof(int));
    MoleculeType[i].Mass = 0;
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      bool test = false;
      int type = FindBeadType(strtok(line, " \t "), Counts, BeadType);
      for (int k = 0; k < MoleculeType[i].nBTypes; k++) {
        if (MoleculeType[i].BType[k] == type) {
          test = true;
          break;
        }
      }
      if (!test) {
        MoleculeType[i].nBTypes++;
        MoleculeType[i].BType = realloc(MoleculeType[i].BType,MoleculeType[i].nBTypes*sizeof(int));
        MoleculeType[i].BType[MoleculeType[i].nBTypes-1] = type;
      }
      MoleculeType[i].Mass += BeadType[type].Mass;
    } //}}}
    // number of bonds //{{{
    fgets(line, sizeof(line), fr);
    strtok(line, " \t ");
    MoleculeType[i].nBonds = atoi(strtok(NULL, " \t")); //}}}
    // connectivity //{{{
    MoleculeType[i].Bond = malloc(MoleculeType[i].nBonds*sizeof(int *));
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      MoleculeType[i].Bond[j] = calloc(2, sizeof(int));
      fgets(line, sizeof(line), fr);
      strtok(line, " \t ");
      MoleculeType[i].Bond[j][0] = atoi(strtok(NULL, " \t")) - 1;
      MoleculeType[i].Bond[j][1] = atoi(strtok(NULL, " \t")) - 1;
    } //}}}

    // skip till 'finish' //{{{
    while(fgets(line, sizeof(line), fr)) {
      char *split;
      split = strtok(line, " \t\n");
      if (strcmp(split, "finish") == 0 ||
          strcmp(split, "Finish") == 0 ||
          strcmp(split, "FINISH") == 0 ) {
        break;
      }
    } //}}}
  } //}}}
  //}}}

//TODO now: fill Molecule[].Bead
  // fill arrays - based on a single molecule types
  Counts.Bonded = 0;
  int i = 0; // at some point, there will be more molecule types
  MoleculeType[i].Number = Counts.Molecules;
  Counts.Bonded += MoleculeType[i].Number * MoleculeType[i].nBeads;
  Counts.Unbonded = Counts.Beads - Counts.Bonded;
  count = Counts.Unbonded;
  for (int j = 0; j < Counts.Molecules; j++) {
    Molecule[j].Bead = calloc(MoleculeType[i].nBeads, sizeof(int));
    Molecule[j].Type = i;
    for (int k = 0; k < MoleculeType[i].nBeads; k++) {
      Molecule[j].Bead[k] = count;
      count++;
    }
  }
  // fill Bead array
  for (int i = 0; i < Counts.Unbonded; i++) {
    Bead[i].Type = 0;
    Bead[i].Molecule = -1;
    Bead[i].Index = i;
  }
  for (int i = 0; i < Counts.Molecules; i++) {
    int mol_type = Molecule[i].Type; // for now, always 0
    for (int j = 0; j < MoleculeType[mol_type].nBeads; j++) {
    //TODO now: fill Molecule[].Bead - here? Or when reading the FIELD?
    //when reading the FIELD before, prototypes were created - add bead types to that
    }
  }
//// realloc Bead and Molecule arrays //{{{
//Molecule = realloc(Molecule, (Counts.Molecules+MoleculeType[i].Number)*sizeof(struct Molecule));
//for (int j = 0; j < MoleculeType[i].Number; j++) {
//  Molecule[Counts.Molecules+j].Bead = calloc(MoleculeType[i].nBeads, (sizeof(int)));
//  Molecule[Counts.Molecules+j].Type = i;
//} //}}}
//    BeadType[type].Number += MoleculeType[i].Number;
//    // fill Bead & Molecule arrays
//    for (int k = 0; k < MoleculeType[i].Number; k++) {
//      int id = Counts.Beads + MoleculeType[i].nBeads * k + j;
//      Bead[id].Type = type;
//      Bead[id].Molecule = k;
//      Bead[id].Index = id;
//      Molecule[Counts.Molecules+k].Bead[j] = id;
//    }

//// calculate total number of beads //{{{
//Counts.Beads = 0;
//for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  Counts.Beads += BeadType[i].Number;
//}
//Counts.Bonded = Counts.Beads - Counts.Unbonded;
//Counts.BeadsInVsf = Counts.Beads; //}}}

  // calculate total number of molecules //{{{
  Counts.Molecules = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    Counts.Molecules += MoleculeType[i].Number;
  } //}}}

  // write all molecules & beads //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = true;
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = true;
  } //}}}
  fclose(fr); //}}}

  // create & fill output vsf file
  WriteVsf(output, Counts, BeadType, Bead, MoleculeType, Molecule);

  // open FIELD-like file //{{{
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // skip to molecule section //{{{
  while(fgets(line, sizeof(line), fr)) {
    char *split;
    split = strtok(line, " \t ");
    if (strncmp(split, "molecule", 8) == 0 ||
        strncmp(split, "Molecule", 8) == 0 ||
        strncmp(split, "MOLECULE", 8) == 0 ) {
      Counts.TypesOfMolecules = atoi(strtok(NULL, " \t"));
      break;
    }
  } //}}}

  // read bond lenghts for all molecule types and create prototype molecules //{{{
  double *prototype_z[Counts.TypesOfMolecules]; // prototype molecule z coordinate
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // allocate array for bond lengths
    prototype_z[i] = calloc(MoleculeType[i].nBeads,sizeof(double));
    prototype_z[i][0] = 0;

    // skip to bonds section //{{{
    while(fgets(line, sizeof(line), fr)) {
      char *split;
      split = strtok(line, " \t");
      if (strcmp(split, "bonds") == 0 ||
          strcmp(split, "Bonds") == 0 ||
          strcmp(split, "BONDS") == 0 ) {
        break;
      }
    } //}}}

    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      fgets(line, sizeof(line), fr);
      char *split;
      split = strtok(line, " \t"); // bond type
      split = strtok(NULL, " \t"); // first id
      split = strtok(NULL, " \t"); // second id
      split = strtok(NULL, " \t"); // spring constant
      split = strtok(NULL, " \t"); // length - what is needed
      double length = atof(split);
      if (fabs(length) < 0.01) {
        length = 0.7;
      }
      prototype_z[i][j+1] = prototype_z[i][j] + length;
    }

    printf("\nbox = (%lf, %lf, %lf)\n\nPrototype %s molecule (z coordinate):\n", BoxLength.x, BoxLength.y, BoxLength.z, MoleculeType[i].Name);
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      printf("%2d: %lf\n", j, prototype_z[i][j]);
    }
  } //}}}

  for (int i = 0; i < mols[0]; i++) {
    for (int j = 0; j < mols[1]; j++) {

    }
  }

  // open output .vcf file for writing //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  } //}}}

  strcpy(stuff, "# by GenBrush\n");

  // write pbc
  fprintf(out, "pbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);
  // write coordinates
  WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
  fclose(out);

  // print information - verbose option //{{{
  if (verbose) {
    char null[1] = {'\0'};
    putchar('\n');
    putchar('\n');
    VerboseOutput(false, null, null, Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    free(prototype_z[i]);
  }
  free(input); //}}}

  return 0;
}
