#include "../AnalysisTools.h"

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
  fprintf(ptr, "      --version           print version number and exit\n");
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
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-nm") != 0 &&
        strcmp(argv[i], "-g") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

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
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: option '-s' requires two numeric arguments\n\n");
    fprintf(stderr, "\033[0m");
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
  char *input = calloc(LINE, sizeof(char));
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

  // structures for info from field //{{{
  BEADTYPE *BeadType_field;
  MOLECULETYPE *MoleculeType_field;
  MOLECULE *Molecule_field;
  BEAD *Bead_field;
  COUNTS Counts_field = InitCounts;
  VECTOR BoxLength;
  int *Index_field;
  PARAMS *bond_type;
  PARAMS *angle_type; //}}}

  ReadField(input, &BoxLength, &Counts_field, &BeadType_field, &Bead_field, &Index_field, &MoleculeType_field, &Molecule_field, &bond_type, &angle_type);

  // structures for info for final system //{{{
  BEADTYPE *BeadType;
  MOLECULETYPE *MoleculeType;
  MOLECULE *Molecule;
  BEAD *Bead;
  COUNTS Counts = InitCounts; //}}}

  // TODO: for now, only one molecule type is used
  Counts.TypesOfMolecules = 1;
  MoleculeType = calloc(Counts.TypesOfMolecules, sizeof(struct MoleculeType));

  // set all bead types to not use; to be decided based on what bead types are in molecules
  for (int i = 0; i < Counts_field.TypesOfBeads; i++) {
    BeadType_field[i].Use = false;
  }

  // copy the molecule type to a new struct - only the one for now //{{{
  // for now, only one molecule type is used
  Counts.TypesOfBeads = 0; // count only types in a molecule
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    strcpy(MoleculeType[i].Name, MoleculeType_field[i].Name);
    MoleculeType[i].nBeads = MoleculeType_field[i].nBeads;
    MoleculeType[i].Bead = calloc(MoleculeType[i].nBeads, sizeof(int));
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      int btype = MoleculeType_field[i].Bead[j];
      MoleculeType[i].Bead[j] = btype;
      // increment number of bead types if btype wasn't yet used
      if (!BeadType_field[btype].Use) {
        Counts.TypesOfBeads++;
      }
      BeadType_field[btype].Use = true; // will be copied later to a new struct
    }
    MoleculeType[i].nBonds = MoleculeType_field[i].nBonds;
    MoleculeType[i].Bond = calloc(MoleculeType[i].nBonds, sizeof(int *));
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      MoleculeType[i].Bond[j] = calloc(2, sizeof(int));
      MoleculeType[i].Bond[j][0] = MoleculeType_field[i].Bond[j][0];
      MoleculeType[i].Bond[j][1] = MoleculeType_field[i].Bond[j][1];
    }
    MoleculeType[i].nBTypes = MoleculeType_field[i].nBTypes;
    MoleculeType[i].BType = calloc(MoleculeType[i].nBTypes, sizeof(int));
    for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
      MoleculeType[i].BType[j] = MoleculeType_field[i].BType[j];
    }
    MoleculeType[i].Mass = MoleculeType_field[i].Mass;
  } //}}}

  // copy bead types to be used //{{{
  BeadType = calloc(Counts.TypesOfBeads, sizeof(struct BeadType));
  count = 0;
  for (int i = 0; i < Counts_field.TypesOfBeads; i++) {
    if (BeadType_field[i].Use) {
      strcpy(BeadType[count].Name, BeadType_field[i].Name);
      BeadType[count].Charge = BeadType_field[i].Charge;
      BeadType[count].Mass = BeadType_field[i].Mass;
      count++;
    }
  } //}}}

  // number of molecules in x and y directions //{{{
  int mols[2]; // maximum number of molecules in x and y directions (per one wall)
  if (number_of_mols == 0) { // if not specified via -nm option
    mols[0] = round(BoxLength.x / spacing[0]);
    mols[1] = round(BoxLength.y / spacing[1]);
    Counts.Molecules = mols[0] * mols[1] * 2;
  } else { // if -nm option is used
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

  // assumes number bead density equal to 3 for number of beads (if '-n' isn't used) //{{{
  if (number_of_beads > 0) {
    Counts.BeadsInVsf = number_of_beads;
  } else {
    Counts.BeadsInVsf = 3 * BoxLength.x * BoxLength.y * BoxLength.z;
  }
  Counts.Beads = Counts.BeadsInVsf; //}}}

  // number of molecules of the one possible molecule type (for now)
  MoleculeType[0].Number = Counts.Molecules;
  // count bonded beads - for now, just one molecule type
  Counts.Bonded = MoleculeType[0].Number * MoleculeType[0].nBeads;

  // increase number of beads if too small (-n option or default number) //{{{
  if (Counts.Bonded > Counts.Beads) {
    Counts.Beads = Counts.BeadsInVsf = Counts.Bonded;
    fprintf(stderr, "\033[1;33m");
    fprintf(stderr, "\nWARNING: too few beads to fit the required number of molecules!\n");
    fprintf(stderr, "         Number of beads increased to %d\n", Counts.Beads);
    fprintf(stderr, "\033[0m");
  } //}}}

  // count unbonded beads
  Counts.Unbonded = Counts.Beads - Counts.Bonded;

  // if there are unbonded beads, make them 'None' type //{{{
  if (Counts.Unbonded > 0) {
    int btype = Counts.TypesOfBeads;
    Counts.TypesOfBeads++;
    BeadType = realloc(BeadType, Counts.TypesOfBeads*sizeof(struct BeadType));
    strcpy(BeadType[btype].Name, "None");
    BeadType[btype].Number = Counts.Unbonded;
    BeadType[btype].Charge = 0;
    BeadType[btype].Mass = 1;
  } //}}}

  // allocate arrays
  Molecule = calloc(Counts.Molecules, sizeof(struct Molecule));
  Bead = calloc(Counts.Beads, sizeof(struct Bead));

  // put the possible unbonded beads in the middle of the box - they have to be somewhere //{{{
  for (int i = 0; i < Counts.Unbonded; i++) {
    Bead[i].Position.x = BoxLength.x / 2;
    Bead[i].Position.y = BoxLength.y / 2;
    Bead[i].Position.z = BoxLength.z / 2;
    Bead[i].Type = Counts.TypesOfBeads - 1;
    Bead[i].Molecule = -1;
    Bead[i].Index = i;
  } //}}}

  // fill molecule struct //{{{
  count = Counts.Unbonded;
  for (int i = 0; i < Counts.Molecules; i++) {
    Molecule[i].Type = 0; // for now, there's only one molecule type
    Molecule[i].Bead = calloc(MoleculeType[0].nBeads, sizeof(int));
    for (int j = 0; j < MoleculeType[0].nBeads; j++) {
      Molecule[i].Bead[j] = count;
      Bead[count].Type = MoleculeType_field[0].Bead[j];
      Bead[count].Molecule = i;
      Bead[count].Index = count;
      Bead[count].Aggregate = calloc(1, sizeof(int)); // just to free later
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

  // move the used molecule so that first bead is [0,0,0] //{{{
  for (int i = 1; i < MoleculeType[0].nBeads; i++) {
    int id_field = Molecule_field[0].Bead[i];
    Bead_field[id_field].Position.x -= Bead_field[id_field-i].Position.x;
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
        int id_field = Molecule_field[0].Bead[k];
        Bead[count+k].Position.x = Bead[count].Position.x + Bead_field[id_field].Position.x;
        Bead[count+k].Position.y = Bead[count].Position.y + Bead_field[id_field].Position.y;
        Bead[count+k].Position.z = Bead[count].Position.z + Bead_field[id_field].Position.z;
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
        int id_field = Molecule_field[0].Bead[k];
        Bead[count+k].Position.x = Bead[count].Position.x + Bead_field[id_field].Position.x;
        Bead[count+k].Position.y = Bead[count].Position.y + Bead_field[id_field].Position.y;
        Bead[count+k].Position.z = Bead[count].Position.z - Bead_field[id_field].Position.z;
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

  // count number of beads of each bead type //{{{
  BeadType[0].Number = Counts.Unbonded; // None beads
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      int btype = MoleculeType[i].Bead[j];
      BeadType[btype].Number += MoleculeType[i].Number;
    }
  } //}}}

  // write output vsf file
  WriteVsf(output, Counts, BeadType, Bead, MoleculeType, Molecule, false);

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
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);

  free(BeadType_field);
  FreeMoleculeType(Counts_field, &MoleculeType_field);
  FreeMolecule(Counts_field, &Molecule_field);
  FreeBead(Counts_field, &Bead_field);
  free(Index_field);
  free(bond_type);
  free(angle_type);
  free(output_vcf);
  free(input); //}}}

  return 0;
}
