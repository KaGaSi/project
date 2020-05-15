#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Warning: This utility was not extensively tested and is in fact not a \
very good generator of initial configuration.\n\n\
GenSystem reads information from a FIELD-like file to create vsf \
structure file and generate coordinates for all beads (used, e.g., as \
initial configuration for a simulation). This utility only creates linear \
molecules no matter the connectivity in the provided FIELD-like file. \
GenSystem also lacks checking for errors in the FIELD-like file so if an \
incorrect file is provided, the utility will exhibit undefined behaviour \
(it will freeze, crash, or produce wrong results).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <out.vsf> <out.vcf> <options>\n\n", cmd);
  fprintf(ptr, "   <out.vsf>     output structure file (vsf format)\n");
  fprintf(ptr, "   <out.vcf>     output coordinate file (vcf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -f <name>  FIELD-like file (default: FIELD)\n");
  fprintf(ptr, "      -v         verbose output\n");
  fprintf(ptr, "      -h         print this help and exit\n");
  fprintf(ptr, "      --version  print version number and exit\n");
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
  // output verbosity //{{{
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  // }}}

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

  // test if <output.vsf> filename ends with '.vsf' (required by VMD)
  int ext = 1;
  char extension[1][5];
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

  // read system information from FIELD-like //{{{
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // read box size //{{{
  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t");
  fgets(line, sizeof(line), fr);
  int words = SplitLine(split, line, delim);
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

  // read number of bead types //{{{
  bool missing = true; // is 'species' keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    words = SplitLine(split, line, delim);
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
      Counts.TypesOfBeads = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    fprintf(stderr, "\nError: %s - missing 'species' line\n\n", input);
    exit(1);
  } //}}}

  BeadType = calloc(Counts.TypesOfBeads,sizeof(struct BeadType));
  Bead = malloc(1*sizeof(struct Bead));

  // read info about bead types //{{{
  Counts.Unbonded = 0;
  Counts.Beads = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);

    // Error on the line: //{{{
    // 1) empty line
    // 2) less then four strings
    // 3) second string isn't a positive double (mass)
    // 4) third string isn't a double (charge)
    // 5) fifth string isn't an integer (unbonded beads)
    if (words == 1 && split[0][0] == '\0') { // 1)
      fprintf(stderr, "\nError: %s - missing bead type line\n\n", input);
      exit(1);
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
    } //}}}

    strcpy(BeadType[i].Name, split[0]);
    BeadType[i].Mass = atof(split[1]);
    BeadType[i].Charge = atof(split[2]);
    BeadType[i].Number = atoi(split[3]);
    Counts.Unbonded += BeadType[i].Number;

    // realloc & fill Bead array
    Bead = realloc(Bead, Counts.Unbonded*sizeof(struct Bead));
    for (int j = Counts.Beads; j < Counts.Unbonded; j++) {
      Bead[j].Type = i;
      Bead[j].Molecule = -1;
      Bead[j].Index = j;
    }

    Counts.Beads = Counts.Unbonded;
  } //}}}

  // read number of molecule types //{{{
  missing = true; // is molecule keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    words = SplitLine(split, line, delim);
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

  MoleculeType = calloc(Counts.TypesOfMolecules,sizeof(struct MoleculeType));
  Molecule = calloc(1,sizeof(struct Molecule));

  // read info about molecule types (and calculate numbers of beads and stuff) //{{{
  Counts.Bonded = 0;
  Counts.Molecules = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // name //{{{
    fgets(line, sizeof(line), fr);
    SplitLine(split, line, delim);
    if (split[0][0] == '\0') { // empty line
      fprintf(stderr, "\nError: %s - blank line instead of molecule name.\n", input);
      fprintf(stderr, "       Beware that a missing 'finish' keyword in multi-molecule systems can make the error somewhat ambiguous.\n");
      exit(1);
    }
    strcpy(MoleculeType[i].Name, split[0]); //}}}
    // number of molecules 'i' //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // error if not 'nummols <int>'
    if (words < 2 ||
        strcasecmp(split[0], "nummols") != 0 ||
        !IsInteger(split[1])) {
      fprintf(stderr, "\nError: %s - wrong or missing 'nummols' line\n", input);
      fprintf(stderr, "       Wrong line:");
      for (int i = 0; i < words; i++) {
        fprintf(stderr, " %s", split[i]);
      }
      fprintf(stderr, "\n\n");
      exit(1);
    }
    MoleculeType[i].Number = atoi(split[1]); //}}}
    // number of beads //{{{
    fgets(line, sizeof(line), fr);
    // error if not 'beads <int>'
    if ((words = SplitLine(split, line, delim)) < 2 ||
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
    // number of bonded beads
    Counts.Bonded += MoleculeType[i].Number * MoleculeType[i].nBeads;
    // realloc Bead and Molecule arrays //{{{
    Bead = realloc(Bead, (Counts.Unbonded+Counts.Bonded)*sizeof(struct Bead));
    Molecule = realloc(Molecule, (Counts.Molecules+MoleculeType[i].Number)*sizeof(struct Molecule));
    for (int j = 0; j < MoleculeType[i].Number; j++) {
      Molecule[Counts.Molecules+j].Bead = calloc(MoleculeType[i].nBeads, (sizeof(int)));
      Molecule[Counts.Molecules+j].Type = i;
    } //}}}
    // number & ids of beadtypes & mass //{{{
    MoleculeType[i].nBTypes = 0;
    MoleculeType[i].BType = calloc(1,sizeof(int));
    MoleculeType[i].Mass = 0;
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line, delim);
      // error - not enough columns or not three coordinate //{{{
      if (words < 4 || !IsDouble(split[1]) || !IsDouble(split[2]) || !IsDouble(split[3])) {
        fprintf(stderr, "\nError: %s - cannot read coordinates\n", input);
        fprintf(stderr, "       Wrong line:");
        for (int k = 0; k < words; k++) {
          fprintf(stderr, " %s", split[k]);
        }
        fprintf(stderr, "\n\n");
        exit(1);
      } //}}}
      bool test = false;
      int type = FindBeadType(split[0], Counts, BeadType);
      // error - wrong bead name //{{{
      if (type == -1) {
        ErrorBeadType(input, split[0], Counts, BeadType);
        exit(1);
      } //}}}
      for (int k = 0; k < MoleculeType[i].nBTypes; k++) {
        if (MoleculeType[i].BType[k] == type) {
          test = true;
          break;
        }
      }
      MoleculeType[i].Bead[j] = type;
      if (!test) {
        MoleculeType[i].nBTypes++;
        MoleculeType[i].BType = realloc(MoleculeType[i].BType,MoleculeType[i].nBTypes*sizeof(int));
        MoleculeType[i].BType[MoleculeType[i].nBTypes-1] = type;
      }
      MoleculeType[i].Mass += BeadType[type].Mass;
      BeadType[type].Number += MoleculeType[i].Number;
      // fill Bead & Molecule arrays
      for (int k = 0; k < MoleculeType[i].Number; k++) {
        int id = Counts.Beads + MoleculeType[i].nBeads * k + j;
        Bead[id].Type = type;
        Bead[id].Molecule = Counts.Molecules + k;
        Bead[id].Index = id;
        Molecule[Counts.Molecules+k].Bead[j] = id;
      }
    } //}}}
    // number of bonds //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
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
      words = SplitLine(split, line, delim);
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
    // total number of beads
    Counts.Beads = Counts.Unbonded + Counts.Bonded;
    // total number of molecules
    Counts.Molecules = MoleculeType[i].Number;
    // skip till 'finish' //{{{
    missing = true;
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
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

  // allocate Bead Aggregate array - needed only to free() //{{{
  for (int i = 0; i < Counts.Beads; i++) {
    Bead[i].Aggregate = calloc(10,sizeof(double));
  } //}}}

  // calculate total number of beads //{{{
  Counts.Beads = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    Counts.Beads += BeadType[i].Number;
  }
  Counts.Bonded = Counts.Beads - Counts.Unbonded;
  Counts.BeadsInVsf = Counts.Beads; //}}}

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

  // Generate coordinates //{{{
  // open FIELD-like file //{{{
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  // skip to molecule section //{{{
  while(fgets(line, sizeof(line), fr)) {
    SplitLine(split, line, delim);
    if (strncmp(split[0], "molecules", 8) == 0 ||
        strncmp(split[0], "Molecules", 8) == 0 ||
        strncmp(split[0], "MOLECULES", 8) == 0 ) {
      Counts.TypesOfMolecules = atoi(split[1]);
      break;
    }
  } //}}}

  // read bond lenghts for all molecule types and create prototype molecules //{{{
  printf("\nbox = (%lf, %lf, %lf)\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  double *prototype_z[Counts.TypesOfMolecules]; // prototype molecule z coordinate
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // allocate array for bond lengths
    prototype_z[i] = calloc(MoleculeType[i].nBeads,sizeof(double));
    prototype_z[i][0] = 0;

    // skip to bonds section //{{{
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strncmp(split[0], "bonds", 4) == 0 ||
          strncmp(split[0], "Bonds", 4) == 0 ||
          strncmp(split[0], "BONDS", 4) == 0 ) {
        break;
      }
    } //}}}

//  for (int j = 0; j < MoleculeType[i].nBonds; j++) {
    // if we assume linear molecule, there must be nBeads-1 bonds
    for (int j = 0; j < (MoleculeType[i].nBeads-1); j++) {
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      double length = atof(split[4]); // bond length is in fifth columne
      if (fabs(length) < 0.01) {
        length = 0.7;
      }
      prototype_z[i][j+1] = prototype_z[i][j] + length;
    }

    printf("Prototype %s molecule (z coordinate):\n", MoleculeType[i].Name);
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      printf("%2d: %lf\n", j+1, prototype_z[i][j]);
    }
  } //}}}
  fclose(fr);

  // longest molecule type //{{{
  int longest = -1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].nBeads > longest) {
      longest = i;
    }
  }
  double length = prototype_z[longest][MoleculeType[longest].nBeads-1];
  printf("Length of the molecule: %lf\n", length); //}}}

  double dist = 0.7; // distance between layers (and beads)
  int x_n = BoxLength.x / dist, y_n = BoxLength.y / dist, z_n = BoxLength.z / dist; // number of positions per axes
  int x = 0, y = 0; // coordinates that are incrementally increased
  double z = 0.1; // coordinates that are incrementally increased
  int count_free = 0; // count number of unbonded beads that are placed

  int layer[Counts.TypesOfMolecules];
  int total_layers = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    layer[i] = x_n * (int)(BoxLength.z / (prototype_z[i][MoleculeType[i].nBeads-1] + dist));
    printf("Maximum number of %s molecules per layer: %d\n", MoleculeType[i].Name, layer[i]);
    layer[i] = MoleculeType[i].Number / layer[i] + 1;
    total_layers += layer[i];
  }
  int sol_layers = Counts.Unbonded / (x_n * z_n);
  int free_layer_per_mol_layer = sol_layers / total_layers;
  printf("%d molecule layers and %d solvent layers (%d solvent layers per one molecule layer).\n", total_layers, sol_layers, free_layer_per_mol_layer);

  // molecules layers with solvent layers in between //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int id = Molecule[i].Bead[j];
      Bead[id].Position.x = (x % x_n) * dist + 0.1;
      Bead[id].Position.y = (y % y_n) * dist + 0.1;
      Bead[id].Position.z = z + prototype_z[type][j];
    }

    z += prototype_z[type][MoleculeType[type].nBeads-1] + dist;
    // add solvent when molecule do not fit till the end of z direction
    if ((z+prototype_z[type][MoleculeType[type].nBeads-1]) >= BoxLength.z) {
      for (; count_free < Counts.Unbonded; count_free++) {
        Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
        Bead[count_free].Position.y = (y % y_n) * dist + 0.1;
        Bead[count_free].Position.z = z;
        z += dist;
        if (z >= BoxLength.z) {
          break;
        }
      }
      z = 0.1;
    }
    // if z is filled, change also x and y
    if (z == 0.1) {
      x++;
      int y_old = y;
      y = x / x_n;
      // when y changes, whole xz layer of molecules was added -- add xz layers of unbonded beads
      if ((y-y_old) > 0) {
        for (int j = 0; j < free_layer_per_mol_layer; j++) {
          for (; count_free < Counts.Unbonded; count_free++) {
            Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
            Bead[count_free].Position.y = (y % y_n) * dist + 0.1;
            Bead[count_free].Position.z = z;
            z += dist;
            if (z >= BoxLength.z) {
              z = 0.1;
            }
            if (z == 0.1) {
              x++;
              y_old = y;
              y = x / x_n;
              // when y changes, whole xz layer of unbonded beads was added
              if ((y-y_old) > 0) {
                break;
              }
            }
          }
        }
      }
    }
  } //}}}

  // remaining unbonded beads //{{{
  for (; count_free < Counts.Unbonded; count_free++) {
    Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
    Bead[count_free].Position.y = (y % y_n) * dist + 0.1;
    Bead[count_free].Position.z = z;
    z += dist;
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z == 0.1) {
      x++;
      y = x / x_n;
    }
  } //}}}

  // tweak y_dist and do it again //{{{
  double y_dist = ((BoxLength.y - 0.5) / (y * dist)) * dist;
  x = 0; // coordinates that are incrementally increased
  y = 0; // coordinates that are incrementally increased
  z = 0.1; // coordinates that are incrementally increased
  double y_coor = 0.1;
  count_free = 0; // count number of unbonded beads that are placed

  // molecules layers with solvent layers in between //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int id = Molecule[i].Bead[j];
      Bead[id].Position.x = (x % x_n) * dist + 0.1;
      Bead[id].Position.y = y_coor;
      Bead[id].Position.z = z + prototype_z[type][j];
    }

    z += prototype_z[type][MoleculeType[type].nBeads-1] + dist;
    // add solvent when molecule do not fit till the end of z direction
    if ((z+prototype_z[type][MoleculeType[type].nBeads-1]) >= BoxLength.z) {
      for (; count_free < Counts.Unbonded; count_free++) {
        Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
        Bead[count_free].Position.y = y_coor;
        Bead[count_free].Position.z = z;
        z += dist;
        if (z >= BoxLength.z) {
          break;
        }
      }
      z = 0.1;
    }
    // if z is filled, change also x and y
    if (z == 0.1) {
      x++;
      int y_old = y;
      y = x / x_n;
      // when y changes, whole xz layer of molecules was added -- add xz layers of unbonded beads
      if ((y-y_old) > 0) {
        y_coor += y_dist;
        for (int j = 0; j < free_layer_per_mol_layer; j++) {
          for (; count_free < Counts.Unbonded; count_free++) {
            Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
            Bead[count_free].Position.y = y_coor;
            Bead[count_free].Position.z = z;
            z += dist;
            if (z >= BoxLength.z) {
              z = 0.1;
            }
            if (z == 0.1) {
              x++;
              y_old = y;
              y = x / x_n;
              // when y changes, whole xz layer of unbonded beads was added
              if ((y-y_old) > 0) {
                y_coor += y_dist;
                break;
              }
            }
          }
        }
      }
    }
  } //}}}

  // remaining unbonded beads //{{{
  for (; count_free < Counts.Unbonded; count_free++) {
    Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
    Bead[count_free].Position.y = y_coor;
    Bead[count_free].Position.z = z;
    z += dist;
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z == 0.1) {
      x++;
      int y_old = y;
      y = x / x_n;
      if ((y-y_old) > 0) {
        y_coor += y_dist;
      }
    }
  } //}}}
  //}}}
  //}}}

  // open output .vcf file for writing //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  } //}}}

  strcpy(stuff, "# by AddToSystem\n");

  // write pbc
  fprintf(out, "pbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);
  // write coordinates
  WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
  fclose(out);

  // print information - verbose option //{{{
  if (verbose) {
    char null[1] = {'\0'};
    fprintf(stdout, "\n\n");
    VerboseOutput(null, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    free(prototype_z[i]);
  }
  free(input);
  free(output_vcf); //}}}

  return 0;
}
