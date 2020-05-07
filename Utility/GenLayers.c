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
  char line[LINE], split[30][100];
  fgets(line, sizeof(line), fr);
  int words = SplitLine(split, line);
  // first line of the FIELD must be: <double> <double> <double> ...whatever
  // Error if:
  // 1) too few strings
  // 2) enough strings, but not numbers
  if (words < 3 || // 1
      !IsPosDouble(split[0]) || // 2) first <double>
      !IsPosDouble(split[1]) || // 2) first <double>
      !IsPosDouble(split[2])) { // 2) first <double>
    fprintf(stderr, "\nError: %s - first line must start with box size (i.e., three positive numbers)\n", input);
    fprintf(stderr, "       Wrong line:");
    for (int i = 0; i < words; i++) {
      fprintf(stderr, " %s", split[i]);
    }
    fprintf(stderr, "\n\n");
    exit(1);
  }
  Vector BoxLength;
  BoxLength.x = atof(split[0]);
  BoxLength.y = atof(split[1]);
  BoxLength.z = atof(split[2]); //}}}

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
  bool missing = true; // is 'species' keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    words = SplitLine(split, line);
    if (strcasecmp(split[0], "species") == 0) {
      missing = false;
      // check if the next string is a number
      if (words < 2 ||            // missing next string
          !IsInteger(split[1])) { // next string isn't a number
        fprintf(stderr, "\nError: %s - missing number of species\n", input);
        fprintf(stderr, "       Wrong line:");
        for (int j = 0; j < words; j++) {
          fprintf(stderr, " %s", split[j]);
        }
        fprintf(stderr, "\n\n");
        exit(1);
      }
      Counts.TypesOfBeads = atoi(split[1]) + 1;
      break;
    }
  }
  if (missing) {
    fprintf(stderr, "\nError: %s - missing 'species' line\n\n", input);
    exit(1);
  } //}}}

  BeadType = calloc(Counts.TypesOfBeads,sizeof(struct BeadType));

  // read info about bead types //{{{
  // number of beads is ignored, because GenLayers only generates the brush;
  // other beads/molecules can be added using AddToSystem
  strcpy(BeadType[0].Name, "None"); // BeadType[0] is a placeholder for AddToSystem
  BeadType[0].Mass = 1;
  BeadType[0].Charge = 0;
  BeadType[0].Number = 0;
  for (int i = 1; i < Counts.TypesOfBeads; i++) {
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line);
    // Error:
    // 1) empty line
    // 2) less then four strings
    // 3) second string isn't a positive double (mass)
    // 4) third string isn't a double (charge)
    // 5) fifth string isn't an integer (unbonded beads)
    if (words == 1 && split[0][0] == '\0') { // 1)
      fprintf(stderr, "\nError: %s - missing bead type line\n\n", input);
    } else if (words < 4 ||                  // 2)
               !IsPosDouble(split[1]) ||     // 3)
               !IsDouble(split[2]) ||        // 4)
               !IsInteger(split[3])) {       // 5)
      fprintf(stderr, "\nError: %s - wrong bead type line:", input);
      for (int j = 0; j < words; j++) {
        fprintf(stderr, " %s", split[j]);
      }
      fprintf(stderr, "\n\n");
      exit(1);
    }
    strcpy(BeadType[i].Name, split[0]);
    BeadType[i].Mass = atof(split[1]);
    BeadType[i].Charge = atof(split[2]);
    BeadType[i].Number = 0;
  } //}}}

  // TODO: generalize - for now, only the first molecule type is used
  // read number of molecule types //{{{
  missing = true; // is molecule keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    words = SplitLine(split, line);
    if (strncasecmp(split[0], "molecule", 8) == 0) {
      missing = false;
      // error - next string isn't a number
      if (!IsInteger(split[1])) {
        fprintf(stderr, "\nError: %s - missing number of molecule types\n", input);
        fprintf(stderr, "       Wrong line:");
        for (int i = 0; i < words; i++) {
          fprintf(stderr, " %s", split[i]);
        }
        fprintf(stderr, "\n\n");
        exit(1);
      }
      Counts.TypesOfMolecules = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    fprintf(stderr, "\nError: %s - missing 'molecule' line\n\n", input);
    exit(1);
  } //}}}

  // TODO: only for now, as only one molecule type is used for now
  Counts.TypesOfMolecules = 1;

  MoleculeType = calloc(Counts.TypesOfMolecules, sizeof(struct MoleculeType));
  Vector **prototype = malloc(Counts.TypesOfMolecules*sizeof(struct Vector *));

  // read info about molecule types (and calculate numbers of beads and stuff) //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // name //{{{
    fgets(line, sizeof(line), fr);
    SplitLine(split, line);
    if (split[0][0] == '\0') { // empty line
      fprintf(stderr, "\nError: %s - blank line instead of molecule name\n", input);
      fprintf(stderr, "       Beware that a missing 'finish' keyword in multi-molecule systems can make the error somewhat ambiguous.\n");
      exit(1);
    }
    strcpy(MoleculeType[i].Name, split[0]); //}}}
    // number of molecules - irrelevant, decided by spacing
    fgets(line, sizeof(line), fr);
    // number of beads //{{{
    fgets(line, sizeof(line), fr);
    // error if not 'beads <int>'
    if ((words = SplitLine(split, line)) < 2 ||
        strncasecmp(split[0], "beads", 4) != 0 ||
        !IsInteger(split[1])) {
      fprintf(stderr, "\nError: %s - wrong or missing 'beads' line\n", input);
      fprintf(stderr, "       Wrong line:");
      for (int i = 0; i < words; i++) {
        fprintf(stderr, " %s", split[i]);
      }
      fprintf(stderr, "\n\n");
      exit(1);
    }
    MoleculeType[i].nBeads = atoi(split[1]);
    MoleculeType[i].Bead = calloc(MoleculeType[i].nBeads, sizeof(int)); //}}}
    // read bead lines //{{{
    prototype[i] = calloc(MoleculeType[i].nBeads, sizeof(Vector));
    // number & ids of beadtypes & mass
    MoleculeType[i].nBTypes = 0;
    MoleculeType[i].BType = calloc(1,sizeof(int));
    MoleculeType[i].Mass = 0;
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line);
      // error - not enough columns or not three coordinate //{{{
      if (words < 4 || !IsDouble(split[1]) || !IsDouble(split[2]) || !IsDouble(split[3])) {
        fprintf(stderr, "\nError: %s - cannot read coordinates\n", input);
        fprintf(stderr, "       Wrong line:");
        for (int i = 0; i < words; i++) {
          fprintf(stderr, " %s", split[i]);
        }
        fprintf(stderr, "\n\n");
        exit(1);
      } //}}}
      bool test = false;
      int type = FindBeadType(split[0], Counts, BeadType);
      // error - wrong bead name //{{{
      if (type == -1) {
        ErrorBeadType(input, split[0], Counts, BeadType);
        fprintf(stderr, "       Note that 'None' bead type is auto-generated, not read from a file.\n\n");
        exit(1);
      } //}}}
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
      prototype[i][j].x = atof(split[1]);
      prototype[i][j].y = atof(split[2]);
      prototype[i][j].z = atof(split[3]);
    } //}}}
    // number of bonds //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line);
    // error if not 'bonds <int>'
    if (words < 2 ||
        strncasecmp(split[0], "bonds", 4) != 0 ||
        !IsInteger(split[1])) {
      fprintf(stderr, "\nError: %s - wrong or missing 'bonds' line\n", input);
      fprintf(stderr, "       Wrong line:");
      for (int i = 0; i < words; i++) {
        fprintf(stderr, " %s", split[i]);
      }
      fprintf(stderr, "\n\n");
      exit(1);
    }
    MoleculeType[i].nBonds = atoi(split[1]); //}}}
    // connectivity //{{{
    MoleculeType[i].Bond = malloc(MoleculeType[i].nBonds*sizeof(int *));
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      MoleculeType[i].Bond[j] = calloc(2, sizeof(int));
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line);
      // error if not '<string> <int> <int>'
      if (words < 3 ||
          !IsInteger(split[1]) || !IsInteger(split[2])) {
        fprintf(stderr, "\nError: %s - wrong or missing 'bonds' line\n", input);
        fprintf(stderr, "       Wrong line:");
        for (int i = 0; i < words; i++) {
          fprintf(stderr, " %s", split[i]);
        }
        fprintf(stderr, "\n\n");
        exit(1);
      }
      MoleculeType[i].Bond[j][0] = atoi(split[1]) - 1;
      MoleculeType[i].Bond[j][1] = atoi(split[2]) - 1;
    } //}}}
    // skip till 'finish' //{{{
    missing = true;
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line);
      if (strcasecmp(split[0], "finish") == 0) {
        missing = false;
        break;
      }
    }
    if (missing) {
      fprintf(stderr, "\nError: %s - missing 'finish' line\n\n", input);
      exit(1);
    } //}}}
  } //}}}
  fclose(fr);
  free(input);

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

  // count number of beads of each bead type //{{{
  BeadType[0].Number = Counts.Unbonded; // None beads
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      int btype = MoleculeType[i].Bead[j];
      BeadType[btype].Number += MoleculeType[i].Number;
    }
  } //}}}

  // remove bead types with 0 beads //{{{
  count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Number != 0) {
      strcpy(BeadType[count].Name, BeadType[i].Name);
      BeadType[count].Number = BeadType[i].Number;
      BeadType[count].Use = BeadType[i].Use;
      BeadType[count].Write = BeadType[i].Write;
      BeadType[count].Charge = BeadType[i].Charge;
      BeadType[count].Mass = BeadType[i].Mass;
      count++;
    }
  }
  Counts.TypesOfBeads = count; //}}}

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
  free(output_vcf); //}}}

  return 0;
}
