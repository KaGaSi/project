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
GenLayers reads information from a FIELD-like file and generates two \
mirror monolayers in xy planes of the simulation box specified distance \
from edges of the box (in z direction). Only the first molecule type in \
the FIELD-like file is used to construct the layers. The molecules are \
generated on a square grid with distance between anchoring points specified \
either explicitly or implicitly (by specifying the number of molecules per \
layer).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <out.vsf> <out.vcf> <options>\n\n", cmd);
  fprintf(ptr, "   <out.vsf>              output structure file (vsf format)\n");
  fprintf(ptr, "   <out.vcf>              output coordinate file (vcf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -s <float> <float>  spacing in x and y directions (default: 1 1)\n");
  fprintf(ptr, "      -n <int>            total number of beads (default: 3*volume_box)\n");
  fprintf(ptr, "      -nm <int>           total number molecules per side (rewrites spacing to fit)\n");
  fprintf(ptr, "      -g <float>          gap between walls and the molecules\n");
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
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-nm") != 0 &&
        strcmp(argv[i], "-g") != 0 &&
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

  // '-n' option - custom total number of beads //{{{
  int number_of_beads = 0;
  if (IntegerOption(argc, argv, "-n", &number_of_beads)) {
    exit(1);
  } //}}}

  // '-nm' option - custom total number of chains //{{{
  int number_of_mols = 0;
  if (IntegerOption(argc, argv, "-nm", &number_of_mols)) {
    exit(1);
  } //}}}

  // '-g' option - gap from both sides of the box in z direction //{{{
  double gap = 0;
  if (DoubleOption(argc, argv, "-g", &gap)) {
    exit(1);
  } //}}}

  // output verbosity
  bool verbose = BoolOption(argc, argv, "-v");

  // FIELD-like file //{{{
  char *input = calloc(LINE, sizeof(char *));
  if (FileOption(argc, argv, "-f", &input)) {
    exit(1);
  }
  if (input[0] == '\0') {
    strcpy(input, "FIELD");
  } //}}}
  //}}}

  count = 0; // count arguments

  // <out.vsf> - output structure file (must end with .vsf) //{{{
  char output[LINE];
  strcpy(output, argv[++count]);

  // test if <output.vsf> filename ends with '.vsf' or '.vtf' (required by VMD)
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output vcf file //{{{
  char *output_vcf = calloc(LINE, sizeof(char));
  char stuff[LINE];
  strcpy(output_vcf, argv[++count]);

  // test if outpuf_vcf has '.vcf' extension - required by vmd //{{{
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}
  //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // read box size //{{{
  char line[LINE], *box[3];
  fgets(line, sizeof(line), fr);
  box[0] = strtok(line, " \t");
  box[1] = strtok(NULL, " \t");
  box[2] = strtok(NULL, " \t");
  Vector BoxLength;
  BoxLength.x = atof(box[0]);
  BoxLength.y = atof(box[1]);
  BoxLength.z = atof(box[2]); //}}}

  // assumes number bead density equal to 3 for number of beads (if '-n' isn't used) //{{{
  if (number_of_beads > 0) {
    Counts.BeadsInVsf = number_of_beads;
  } else {
    Counts.BeadsInVsf = 3 * BoxLength.x * BoxLength.y * BoxLength.z;
  }
  Counts.Beads = Counts.BeadsInVsf; //}}}

  // number of molecules in x and y directions //{{{
  int mols[2];
  if (number_of_mols == 0) {
    mols[0] = round(BoxLength.x / spacing[0]);
    mols[1] = round(BoxLength.y / spacing[1]);
    Counts.Molecules = mols[0] * mols[1] * 2;
  } else {
    double mols_sqrt = sqrt(number_of_mols);
    mols[0] = (int)(mols_sqrt);
    mols[1] = mols[0] + (number_of_mols - SQR(mols[0])) / mols[0];
    if ((mols[0]*mols[1]-number_of_mols) != 0) {
      mols[0]++;
    }
    spacing[0] = BoxLength.x / mols[0];
    spacing[1] = BoxLength.y / mols[1];
    Counts.Molecules = 2 * number_of_mols;
  } //}}}

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
  // number of beads is ignored, because GenLayers only generates the brush;
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

  MoleculeType = calloc(Counts.TypesOfMolecules, sizeof(struct MoleculeType));
  Vector **prototype = malloc(Counts.TypesOfMolecules*sizeof(struct Vector *));

  // read info about molecule types (and calculate numbers of beads and stuff) //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // name //{{{
    fgets(line, sizeof(line), fr);
    strcpy(line, TrimLine(line)); // trim excess whitespace
    strcpy(MoleculeType[i].Name, strtok(line, " \t")); //}}}
    // number of molecules - irrelevant, decided by spacing //{{{
    fgets(line, sizeof(line), fr);
    strtok(line, " \t "); //}}}
    // number of beads //{{{
    fgets(line, sizeof(line), fr);
    strtok(line, " \t ");
    MoleculeType[i].nBeads = atoi(strtok(NULL, " \t"));
    MoleculeType[i].Bead = calloc(MoleculeType[i].nBeads, sizeof(int));
    prototype[i] = calloc(MoleculeType[i].nBeads, sizeof(Vector)); //}}}
    // number & ids of beadtypes & mass //{{{
    MoleculeType[i].nBTypes = 0;
    MoleculeType[i].BType = calloc(1,sizeof(int));
    MoleculeType[i].Mass = 0;
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      bool test = false;
      int type = FindBeadType(strtok(line, " \t "), Counts, BeadType);
      MoleculeType[i].Bead[j] = type;
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
      // read coordinates
      prototype[i][j].x = atof(strtok(NULL, " \t"));
      prototype[i][j].y = atof(strtok(NULL, " \t"));
      prototype[i][j].z = atof(strtok(NULL, " \t"));
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
  fclose(fr);

  // TODO: generalize for more moltypes
  // fill arrays - based on a single molecule type for now //{{{
  Counts.Bonded = 0;
  int i = 0; // at some point, there will be more molecule types
  MoleculeType[i].Number = Counts.Molecules;
  Counts.Bonded += MoleculeType[i].Number * MoleculeType[i].nBeads; //}}}

  // increase number of beads if too small (-n option or default number) //{{{
  if (Counts.Bonded > Counts.Beads) {
    Counts.Beads = Counts.BeadsInVsf = Counts.Bonded;
    fprintf(stderr, "\nWARNING: too few beads to fit the required number of molecules!\n");
    fprintf(stderr, "   Number of beads increased to %d\n", Counts.Beads);
  } //}}}

  // allocate Bead array //{{{
  Bead = malloc(Counts.Beads*sizeof(struct Bead));
  for (int i = 0; i < Counts.Beads; i++) {
    Bead[i].Aggregate = calloc(1,sizeof(double));
  } //}}}

  // fill Molecule array and Bead[].Molecule //{{{
  Counts.Unbonded = Counts.Beads - Counts.Bonded;
  count = Counts.Unbonded;
  for (int j = 0; j < Counts.Molecules; j++) {
    Molecule[j].Bead = calloc(MoleculeType[i].nBeads, sizeof(int));
    Molecule[j].Type = i;
    for (int k = 0; k < MoleculeType[i].nBeads; k++) {
      Molecule[j].Bead[k] = count;
      Bead[count].Molecule = j;
      count++;
    }
  } //}}}

  // set all molecules & beads to write //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = true;
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = true;
  } //}}}

  // put unbonded beads in the middle of the box - they have to be somewhere //{{{
  for (int i = 0; i < Counts.Unbonded; i++) {
    Bead[i].Position.x = BoxLength.x / 2;
    Bead[i].Position.y = BoxLength.y / 2;
    Bead[i].Position.z = BoxLength.z / 2;
    Bead[i].Type = 0;
    Bead[i].Molecule = -1;
    Bead[i].Index = i;
  } //}}}

  // generate brush on the grid & fill remaining Bead[] array //{{{
  count = Counts.Unbonded;
  int count_mols = 0;
  // z=0
  for (int i = 0; i < mols[0]; i++) {
    for (int j = 0; j < mols[1]; j++) {
      Bead[count].Position.x = spacing[0] / 2 + spacing[0] * i;
      Bead[count].Position.y = spacing[1] / 2 + spacing[1] * j;
      Bead[count].Position.z = 0;
      Bead[count].Index = count;
      Bead[count].Type = MoleculeType[0].Bead[0];
      for (int k = 1; k < MoleculeType[0].nBeads; k++) {
        Bead[count+k].Position.x = Bead[count].Position.x + (prototype[0][k].x - prototype[0][0].x);
        Bead[count+k].Position.y = Bead[count].Position.y + (prototype[0][k].y - prototype[0][0].y);
        Bead[count+k].Position.z = Bead[count].Position.z + (prototype[0][k].z - prototype[0][0].z);
        Bead[count+k].Index = count + k;
        Bead[count+k].Type = MoleculeType[0].Bead[k];
      }
      count += MoleculeType[0].nBeads;
      if (++count_mols == number_of_mols) {
        break;
      }
    }
    if (count_mols == number_of_mols) {
      break;
    }
  }
  // z=BoxLength
  count_mols = 0;
  for (int i = 0; i < mols[0]; i++) {
    for (int j = 0; j < mols[1]; j++) {
      Bead[count].Position.x = spacing[0] / 2 + spacing[0] * i;
      Bead[count].Position.y = spacing[1] / 2 + spacing[1] * j;
      Bead[count].Position.z = BoxLength.z;
      Bead[count].Index = count;
      Bead[count].Type = MoleculeType[0].Bead[0];
      for (int k = 1; k < MoleculeType[0].nBeads; k++) {
        Bead[count+k].Position.x = Bead[count].Position.x + (prototype[0][k].x - prototype[0][0].x);
        Bead[count+k].Position.y = Bead[count].Position.y + (prototype[0][k].y - prototype[0][0].y);
        Bead[count+k].Position.z = Bead[count].Position.z - (prototype[0][k].z - prototype[0][0].z);
        Bead[count+k].Index = count + k;
        Bead[count+k].Type = MoleculeType[0].Bead[k];
      }
      count += MoleculeType[0].nBeads;
      if (++count_mols == number_of_mols) {
        break;
      }
    }
    if (count_mols == number_of_mols) {
      break;
    }
  } //}}}

  // open output .vcf file for writing //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  } //}}}

  // string with command to be printed to output vsf //{{{
  strcpy(stuff, "# Generated by: GenLayers ");
  for (int i = 1; i < argc; i++) {
    strcat(stuff, argv[i]);
    strcat(stuff, " ");
  }
  strcat(stuff, "\n"); //}}}

  // TODO: redo so that the None type is added only if needed
  // remove None bead type if there are no extra beads //{{{
  if (Counts.Bonded == Counts.Beads) {
    for (int i = 1; i < Counts.TypesOfBeads; i++) {
      strcpy(BeadType[i-1].Name, BeadType[i].Name);
      BeadType[i-1].Number = BeadType[i].Number;
      BeadType[i-1].Use = BeadType[i].Use;
      BeadType[i-1].Write = BeadType[i].Write;
      BeadType[i-1].Charge = BeadType[i].Charge;
      BeadType[i-1].Mass = BeadType[i].Mass;
    }
    for (int i = 0; i < Counts.Beads; i++) {
      Bead[i].Type--;
    }
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      for (int j = 0; j < MoleculeType[i].nBeads; j++) {
        MoleculeType[i].Bead[j]--;
      }
    }
    Counts.TypesOfBeads--;
  } //}}}

  // write output vsf file
  WriteVsf(output, Counts, BeadType, Bead, MoleculeType, Molecule);

  // adjust for gap value (-g option) //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    int mtype = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
      int id = Molecule[i].Bead[j];
      if (i < (Counts.Molecules/2)) {
        Bead[id].Position.z += gap;
      } else {
        Bead[id].Position.z -= gap;
      }
    }
  } //}}}

  // write output vcf file //{{{
  // pbc
  fprintf(out, "pbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);
  // coordinates
  WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
  fclose(out); //}}}

  // print information - verbose option //{{{
  if (verbose) {
    fprintf(stdout, "\nGrid of %d x %d ", mols[0], mols[1]);
    if (number_of_mols != 0) {
      fprintf(stdout, "(%d molecules) ", number_of_mols);
    }
    fprintf(stdout, "molecules on each wall");
    char null[1] = {'\0'};
    putchar('\n');
    putchar('\n');
    VerboseOutput(null, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    // molecule prototypes
    putchar('\n');
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(stdout, "Prototype of %s\n", MoleculeType[i].Name);
      for (int j = 0; j < MoleculeType[i].nBeads; j++) {
        int btype = MoleculeType[i].Bead[j];
        fprintf(stdout, "%8s %lf %lf %lf\n", BeadType[btype].Name,
                                             prototype[i][j].x,
                                             prototype[i][j].y,
                                             prototype[i][j].z);
      }
      putchar('\n');
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    free(prototype[i]);
  }
  free(prototype);
  free(input);
  free(output_vcf); //}}}

  return 0;
}
