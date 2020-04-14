#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;

    fprintf(ptr, "\
AddToSystem takes existing coordinates and adds unbonded beads and/or \
molecules to the system. The data for added components is read either \
from a FIELD-like file that has to contain species and molecule sections \
(the same as for DL_MESO simulation program) or from a vsf/vcf files \
(-vtf option). In the first case, new molecules and beads are added \
randomly, while in the second case, the provided coordinates are used. \
The new beads replace either neutral unbonded beads with the lowest \
indices (as written in the vsf file) or option-specified unbonded beads. \
Constraints can be placed on x, y, or z coordinates, as well as on \
distance of the added beads from other specified beads. The constraints \
are ignored for -vtf option. Options -ld and/or -hd must be used in \
accompanied by -bt option.\n\n");
  }
  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input.vcf> ", cmd);
  fprintf(ptr, "<out.vcf> <out.vsf> <options>\n\n");

  fprintf(ptr, "   <input.vcf>         input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <out.vsf>           output structure file (vsf format)\n");
  fprintf(ptr, "   <out.vcf>           output coordinate file (vcf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -f <name>        FIELD-like file with molecules to add (default: FIELD)\n");
  fprintf(ptr, "      -vtf <vsf> <vcf> use vtf file format instead of FIELD (divided to vsf and vcf files)\n");
  fprintf(ptr, "      -st <int>        timestep to add new beads to\n");
  fprintf(ptr, "      -xyz <name>      output xyz file\n");
  fprintf(ptr, "      -ld <float>      specify lowest distance from chosen bead types (default: none)\n");
  fprintf(ptr, "      -hd <float>      specify highest distance from chosen bead types (default: none)\n");
  fprintf(ptr, "      -bt <name(s)>    specify bead types new beads should be far from/near to (default: none)\n");
  fprintf(ptr, "      -cx <num> <num2> constrain x coordinate of randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -cy <num> <num2> constrain y coordinate of randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -cz <num> <num2> constrain z coordinate of randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -gc              use molecule's geometric centre for the distance check instead of its first bead\n");
  fprintf(ptr, "      -xb <name(s)>    specify bead types to exchange new beads for (must be unbonded beads)\n");
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
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-vtf") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-xyz") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-ld") != 0 &&
        strcmp(argv[i], "-hd") != 0 &&
        strcmp(argv[i], "-cx") != 0 &&
        strcmp(argv[i], "-cy") != 0 &&
        strcmp(argv[i], "-cz") != 0 &&
        strcmp(argv[i], "-gc") != 0 &&
        strcmp(argv[i], "-xb") != 0) {
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

  // -f <add> - FIELD-like file with molecules to add //{{{
  char *input_add = calloc(LINE, sizeof(char *));
  if (FileOption(argc, argv, "-f", &input_add)) {
    exit(1);
  }
  if (input_add[0] == '\0') {
    strcpy(input_add, "FIELD");
  } //}}}

  // -vtf <vsf> <vcf> - FIELD-like file with molecules to add //{{{
  char *add_vsf = calloc(LINE, sizeof(char *)),
       *add_vcf = calloc(LINE, sizeof(char *));
  // 1) vsf file
  if (FileOption(argc, argv, "-vtf", &add_vsf)) {
    exit(1);
  }
  // 2) if vsf file exists, look for vcf
  int ext;
  char extension[2][5];
  if (add_vsf[0] != '\0') {
    ext = 2;
    strcpy(extension[0], ".vsf");
    strcpy(extension[1], ".vtf");
    if (ErrorExtension(add_vsf, ext, extension)) {
      Help(argv[0], true);
      exit(1);
    }
    if (add_vsf[strlen(add_vsf)-2] == 't') { // if *.vtf file, use it as a vcf as well
      strcpy(add_vcf, add_vsf);
    } else {
      for (int i = 1; i < argc; i++) { // if *.vsf file, read next file
        if (strcmp(argv[i], "-vtf") == 0) {
          if (argc >= (i+3)) { // is there another cli argument behind '.vsf' file?
            char temp[LINE];
            strcpy(temp, argv[i+1]); // save vsf filename
            strcpy(argv[i+1], argv[i+2]); // copy vcf filename to vsf - required by FileOption()
            ext = 1;
            strcpy(extension[0], ".vcf");
            if (FileOption(argc, argv, "-vtf", &add_vcf)) {
              exit(1);
            }
            if (ErrorExtension(add_vcf, ext, extension)) {
              fprintf(stderr, "       Wrong second filename!\n");
              Help(argv[0], true);
              exit(1);
            }
            strcpy(argv[i+1], temp); // restore vsf filename
            break;
          } else { // there's no more cli arguments
            fprintf(stderr, "\nError: option -vtf is missing second filename!\n");
            Help(argv[0], true);
            exit(1);
          }
        }
      }
    }
  } //}}}

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // save into xyz file? //{{{
  char *output_xyz = calloc(LINE, sizeof(char *));
  if (FileOption(argc, argv, "-xyz", &output_xyz)) {
    exit(1);
  } //}}}

  // specify lowest and/or highest distance from beads of type specified by '-bt' option //{{{
  double lowest_dist = -1, highest_dist = -1;
  if (DoubleOption(argc, argv, "-ld", &lowest_dist)) {
    exit(1);
  }
  if (DoubleOption(argc, argv, "-hd", &highest_dist)) {
    exit(1);
  }
  // if '-ld' and/or '-hd' are present, '-bt' must be too
  if (highest_dist != -1 || lowest_dist != -1) {
    bool bt = false;
    for (int i = 0; i < argc; i++) {
      if (strcmp(argv[i], "-bt") == 0) {
        bt = true;
      }
    }
    if (!bt) {
      fprintf(stderr, "\nError: if '-ld' and/or '-hd' is used, '-bt' must be specified as well\n\n");
      exit(1);
    }
  } //}}}

  // coordinate constraints //{{{
  // x direcion //{{{
  int test = 2;
  double range[2] = {0, 0};
  if (MultiDoubleOption(argc, argv, "-cx", &test, range)) {
    exit(1);
  }
  if (test != 2) {
    fprintf(stderr, "\nError: option '-cx' requires two numeric arguments\n\n");
    exit(1);
  }
  // make sure first number is smaller
  if (range[0] > range[1]) {
    SwapDouble(&range[0], &range[1]);
  }
  Vector constraint[2];
  constraint[0].x = range[0];
  constraint[1].x = range[1]; //}}}
  // y direcion //{{{
  test = 2;
  range[0] = range[1] = 0;
  if (MultiDoubleOption(argc, argv, "-cy", &test, range)) {
    exit(1);
  }
  if (test != 2) {
    fprintf(stderr, "\nError: option '-cy' requires two numeric arguments\n\n");
    exit(1);
  }
  // make sure first number is smaller
  if (range[0] > range[1]) {
    SwapDouble(&range[0], &range[1]);
  }
  constraint[0].y = range[0];
  constraint[1].y = range[1]; //}}}
  // z direcion //{{{
  test = 2;
  range[0] = range[1] = 0;
  if (MultiDoubleOption(argc, argv, "-cz", &test, range)) {
    exit(1);
  }
  if (test != 2) {
    fprintf(stderr, "\nError: option '-cz' requires two numeric arguments\n\n");
    exit(1);
  }
  // make sure first number is smaller
  if (range[0] > range[1]) {
    SwapDouble(&range[0], &range[1]);
  }
  constraint[0].z = range[0];
  constraint[1].z = range[1]; //}}}
  //}}}

  // use centre of mass instead of the first bead for distance check of new molecules //{{{
  bool com = BoolOption(argc, argv, "-gc"); //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input.vcf> - input coordinate file //{{{
  char input_coor[LINE];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  ext = 2;
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vsf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vsf[LINE];
  strcpy(output_vsf, argv[++count]);

  // test if <out.vsf> filename ends with '.vsf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output_vsf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[LINE];
  strcpy(output_vcf, argv[++count]);

  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // -xb <name(s)> - specify what bead types to exchange //{{{
  if (BeadTypeOption(argc, argv, "-xb", true, Counts, &BeadType)) {
    exit(0);
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = BeadType[i].Use; // use Write flag to decide which bead types to use
    BeadType[i].Use = false;
  }//}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", false, Counts, &BeadType)) {
    exit(0);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[LINE];
  // skip till 'pbc' keyword //{{{
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0); //}}}

  // read pbc //{{{
  Vector BoxLength;
  char line[LINE];
  fgets(line, sizeof(line), vcf);
  // split the line into array
  char *split[30];
  split[0] = strtok(line, " \t");
  int i = 0;
  while (split[i] != NULL && i < 29) {
    split[++i] = strtok(NULL, " \t");
  }
  BoxLength.x = atof(split[0]);
  BoxLength.y = atof(split[1]);
  BoxLength.z = atof(split[2]); //}}}
  //}}}

  // print original system //{{{
  if (verbose) {
    fprintf(stdout, "\nORIGINAL SYSTEM\n");
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    if (start > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", start);
    }
  } //}}}

  // define box size using constraints (-c{x,y,z} options)//{{{
  Vector new_box;
  if (constraint[1].x != 0) {
    new_box.x = constraint[1].x - constraint[0].x;
  } else {
    new_box.x = BoxLength.x;
  }
  if (constraint[1].y != 0) {
    new_box.y = constraint[1].y - constraint[0].y;
  } else {
    new_box.y = BoxLength.y;
  }
  if (constraint[1].z != 0) {
    new_box.z = constraint[1].z - constraint[0].z;
  } else {
    new_box.z = BoxLength.z;
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff = calloc(LINE, sizeof(char)); //}}}

  // skip first start-1 steps //{{{
  count = 0;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\nError: premature end of %s file\n\n", input_coor);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent && start > 1) {
    if (script) {
      fprintf(stdout, "Starting step: %d\n", start);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                             ");
      fprintf(stdout, "\rStarting step: %d   \n", start);
    }
  } //}}}
  // is the vcf file continuing?
  if (ErrorDiscard(start, count, input_coor, vcf)) {
    exit(1);
  }
  //}}}

  // read the coordinate timestep to add stuff to //{{{
  if ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      fprintf(stdout, "Using step %6d\n", count);
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
      exit(1);
    } //}}}
  } else {
    fprintf(stderr, "\nWarning: using last step in %s (%d)\n", input_coor, count);
  }
  fclose(vcf); //}}}

  // create structures for new stuff //{{{
  struct Counts Counts_add;
  ZeroCounts(&Counts_add);
  struct Molecule *Molecule_add;
  struct MoleculeType *MoleculeType_add;
  struct BeadType *BeadType_add;
  struct Bead *Bead_add;
  struct Vector **prototype;
  int *Index_add;
  //}}}

  if (add_vsf[0] == '\0') { // read stuff to be added FIELD //{{{
    // open the FIELD-like file for reading //{{{
    FILE *in_add;
    if ((in_add = fopen(input_add, "r")) == NULL) {
      ErrorFileOpen(input_add, 'r');
      exit(1);
    } //}}}

    // read number of bead types //{{{
    // 1) skip till 'species' keyword
    do {
      // get whole line - max 1000 chars
      fgets(line, sizeof(line), in_add);

      // first string of the line
      split[0] = strtok(line, " \t");

    } while (strcmp(split[0], "species") != 0 &&
             strcmp(split[0], "SPECIES") != 0 &&
             strcmp(split[0], "Species") != 0);
    // 2) after 'species' is number of bead types
    split[1] = strtok(NULL, " \t");
    Counts_add.TypesOfBeads = atoi(split[1]); //}}}

    // added bead types
    BeadType_add = malloc(Counts_add.TypesOfBeads*sizeof(struct BeadType));
    // added beads (to be realloc'd later)
    Bead_add = malloc(1*sizeof(struct Bead));

    // read bead type info //{{{
    for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
      fgets(line, sizeof(line), in_add);
      // bead name
      split[0] = strtok(line, " \t");
      strcpy(BeadType_add[i].Name, split[0]);
      // bead mass
      split[1] = strtok(NULL, " \t");
      BeadType_add[i].Mass = atof(split[1]);
      // bead charge
      split[2] = strtok(NULL, " \t");
      BeadType_add[i].Charge = atof(split[2]);
      // number of unbonded beads
      split[3] = strtok(NULL, " \t");
      BeadType_add[i].Number = atoi(split[3]);
      BeadType_add[i].Use = true;
      BeadType_add[i].Write = true;

      // realloc & fill Bead_add
      int total_beads = Counts_add.Unbonded + BeadType_add[i].Number;
      if (total_beads > 0) {
        Bead_add = realloc(Bead_add, total_beads*sizeof(struct Bead));
        for (int j = Counts_add.Unbonded; j < total_beads; j++) {
          Bead_add[j].Type = i;
          Bead_add[j].Molecule = -1;
          Bead_add[j].nAggregates = 0;
          Bead_add[j].Aggregate = calloc(1,sizeof(int));
          Bead_add[j].Index = j; // probably useless here
        }
      }

      Counts_add.Unbonded += BeadType_add[i].Number;
      Counts_add.Beads += BeadType_add[i].Number;
    } //}}}

    // read number of molecul types //{{{
    // 1) skip till 'molecule' keyword
    do {
      // get whole line - max 1000 chars
      fgets(line, sizeof(line), in_add);

      // first string of the line
      split[0] = strtok(line, " \t");

    } while (strncmp(split[0], "molecules", 8) != 0 &&
             strncmp(split[0], "MOLECULES", 8) != 0 &&
             strncmp(split[0], "Molecules", 8) != 0);
    // 2) after 'molecule' is number of molecule types
    split[1] = strtok(NULL, " \t");
    Counts_add.TypesOfMolecules = atoi(split[1]); //}}}

    // allocate structures //{{{
    // added molecule types
    MoleculeType_add = malloc(Counts_add.TypesOfMolecules*sizeof(struct MoleculeType));
    // coordinates of a prototype molecule
    prototype = calloc((Counts_add.TypesOfMolecules+Counts.TypesOfMolecules), sizeof(Vector *));
    // added molecules (to be realloc'd later)
    Molecule_add = malloc(1*sizeof(struct Molecule)); //}}}

    // read molecule info //{{{
    for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
      // molecule name //{{{
      fgets(line, strlen(line), in_add);
      strcpy(line, TrimLine(line)); // trim excess whitespace
      split[0] = strtok(line, " \t");
      strcpy(MoleculeType_add[i].Name, split[0]); //}}}
      // number of molecules of given type //{{{
      fgets(line, sizeof(line), in_add);
      split[0] = strtok(line, " \t");
      split[1] = strtok(NULL, " \t");
      MoleculeType_add[i].Number = atoi(split[1]); //}}}
      // number of beads in molecules of given type //{{{
      fgets(line, sizeof(line), in_add);
      split[0] = strtok(line, " \t");
      split[1] = strtok(NULL, " \t");
      MoleculeType_add[i].nBeads = atoi(split[1]);
      MoleculeType_add[i].Bead = calloc(MoleculeType_add[i].nBeads, sizeof(int)); //}}}
      MoleculeType_add[i].InVcf = true;
      MoleculeType_add[i].Use = true;
      MoleculeType_add[i].Write = true;

      // total number of beads in the molecules of type 'i'
      int beads = MoleculeType_add[i].Number * MoleculeType_add[i].nBeads;
      // total number of molecules and beads so far
      int total_mols = Counts_add.Molecules + MoleculeType_add[i].Number;
      int total_beads = Counts_add.Beads+beads;

      // realloc _add structures //{{{
      Bead_add = realloc(Bead_add, total_beads*sizeof(struct Bead));
      Molecule_add = realloc(Molecule_add, total_mols*sizeof(struct Molecule));
      for (int j = Counts_add.Molecules; j < total_mols; j++) {
        Molecule_add[j].Bead = malloc(MoleculeType_add[i].nBeads*sizeof(int));
      }
      MoleculeType_add[i].BType = malloc(1*sizeof(int));
      MoleculeType_add[i].nBTypes = 0;
      MoleculeType_add[i].Mass = 0; //}}}

      // allocate array for coordinates of prototype molecule of type 'i'
      prototype[i] = malloc(MoleculeType_add[i].nBeads*sizeof(struct Vector));

      // read bead types and coordinates
      Vector zero_first;
      zero_first.x = 0;
      zero_first.y = 0;
      zero_first.z = 0;
      for (int j = 0; j < MoleculeType_add[i].nBeads; j++) {
        fgets(line, sizeof(line), in_add);
        // bead name //{{{
        split[0] = strtok(line, " \t");
        // is the bead type registered in the molecule already?
        int btype = FindBeadType(split[0], Counts_add, BeadType_add);
        bool exists = false;
        for (int k = 0; k < MoleculeType_add[i].nBTypes; k++) {
          if (btype == MoleculeType_add[i].BType[k]) {
            exists = true;
            break;
          }
        }
        // add new bead type to molecule if not yet registered
        if (!exists) {
          MoleculeType_add[i].nBTypes++;
          MoleculeType_add[i].BType = realloc(MoleculeType_add[i].BType, MoleculeType_add[i].nBTypes*sizeof(int));
          MoleculeType_add[i].BType[MoleculeType_add[i].nBTypes-1] = btype;
        } //}}}
        MoleculeType_add[i].Bead[j] = btype;

        // add bead's mass to molecule mass
        MoleculeType_add[i].Mass += BeadType_add[btype].Mass;

        // bead coordinate //{{{
        split[1] = strtok(NULL, " \t");
        prototype[i][j].x = atof(split[1]);
        split[2] = strtok(NULL, " \t");
        prototype[i][j].y = atof(split[2]);
        split[3] = strtok(NULL, " \t");
        prototype[i][j].z = atof(split[3]); //}}}

        // first bead should have coordinates [0,0,0] //{{{
        if (j == 0) {
          zero_first.x = prototype[i][j].x;
          zero_first.y = prototype[i][j].y;
          zero_first.z = prototype[i][j].z;
        } else {
          prototype[i][j].x -= zero_first.x;
          prototype[i][j].y -= zero_first.y;
          prototype[i][j].z -= zero_first.z;
        } //}}}
      }

      // fill _add structures //{{{
      for (int j = 0; j < MoleculeType_add[i].nBeads; j++) {
        int btype = MoleculeType_add[i].Bead[j];
        for (int l = 0; l < MoleculeType_add[i].Number; l++) {
          int index = Counts_add.Beads + l * MoleculeType_add[i].nBeads + j;

          Bead_add[index].Type = btype;
          Bead_add[index].Molecule = Counts_add.Molecules + l;
          Bead_add[index].Index = index;
          Bead_add[index].nAggregates = 0; // useless here
          Bead_add[index].Aggregate = malloc(1*sizeof(int)); // just to free later
          Bead_add[index].Position.x = prototype[i][j].x;
          Bead_add[index].Position.y = prototype[i][j].y;
          Bead_add[index].Position.z = prototype[i][j].z;
          Molecule_add[Counts_add.Molecules+l].Bead[j] = index;
          Molecule_add[Counts_add.Molecules+l].Type = i;
          BeadType_add[btype].Number++;
        }
      } //}}}

      Counts_add.Bonded += beads;
      Counts_add.Beads += beads;
      Counts_add.BeadsInVsf += beads;
      Counts_add.Molecules += MoleculeType_add[i].Number;

      // number of bonds in molecules of given type //{{{
      fgets(line, sizeof(line), in_add);
      split[0] = strtok(line, " \t");
      split[1] = strtok(NULL, " \t");
      MoleculeType_add[i].nBonds = atoi(split[1]); //}}}

      // allocate Bond array //{{{
      MoleculeType_add[i].Bond = malloc(MoleculeType_add[i].nBonds*sizeof(int *));
      for (int j = 0; j < MoleculeType_add[i].nBonds; j++) {
        MoleculeType_add[i].Bond[j] = calloc(2, sizeof(int));
      } //}}}

      // read bond info //{{{
      for (int j = 0; j < MoleculeType_add[i].nBonds; j++) {
        fgets(line, sizeof(line), in_add);
        split[0] = strtok(line, " \t");
        split[1] = strtok(NULL, " \t");
        split[2] = strtok(NULL, " \t");
        MoleculeType_add[i].Bond[j][0] = atoi(split[1]) - 1; // bead ids in FIELD start from 1
        MoleculeType_add[i].Bond[j][1] = atoi(split[2]) - 1;
      } //}}}

      // skip till 'finish' keyword //{{{
      do {
        // get whole line - max 1000 chars
        fgets(line, sizeof(line), in_add);
        strcpy(line, TrimLine(line)); // trim excess whitespace

        // first string of the line
        split[0] = strtok(line, " \t");
      } while (strcmp(split[0], "finish") != 0 &&
               strcmp(split[0], "Finish") != 0 &&
               strcmp(split[0], "FINISH") != 0); //}}}
    } //}}}

    Counts_add.BeadsInVsf = Counts_add.Beads;

    // fill Index_add - it's just 0..Counts_add.Beads
    Index_add = calloc(Counts_add.Beads, sizeof(int)); // not used here; just to later free
    for (int i = 0; i < Counts_add.Beads; i++) {
      Index_add[i] = i;
    }

    fclose(in_add); //}}}
  } else { // read stuff to add from vtf file(s) ('-vtf' option) //{{{
    bool indexed_add = ReadStructure(add_vsf, add_vcf, &Counts_add, &BeadType_add, &Bead_add, &Index_add, &MoleculeType_add, &Molecule_add);

    // open input coordinate file
    if ((vcf = fopen(add_vcf, "r")) == NULL) {
      ErrorFileOpen(add_vcf, 'r');
      exit(1);
    }

    // read coordinates
    if ((test = ReadCoordinates(indexed_add, vcf, Counts_add, Index_add, &Bead_add, &stuff)) != 0) {
      ErrorCoorRead(add_vcf, test, 0, stuff, add_vsf);
      exit(1);
    }
    fclose(vcf);

    int new_moltypes = Counts.TypesOfMolecules + Counts_add.TypesOfMolecules;
    prototype = malloc(new_moltypes*sizeof(Vector *));
    for (int i = 0; i < new_moltypes; i++) {
      if (i < Counts_add.TypesOfMolecules) {
        prototype[i] = calloc(MoleculeType_add[i].nBeads, sizeof(Vector));
      } else {
        prototype[i] = calloc(1, sizeof(Vector));
      }
      prototype[i][0].x = 100000; // just some high number
    }
    for (int i = 0; i < Counts_add.Molecules; i++) {
      int mtype = Molecule_add[i].Type;
      if (prototype[mtype][0].x == 100000) {
        for (int j = 0; j < MoleculeType_add[mtype].nBeads; j++) {
          int id = Molecule_add[i].Bead[j];
          prototype[mtype][j].x = Bead_add[id].Position.x;
          prototype[mtype][j].y = Bead_add[id].Position.y;
          prototype[mtype][j].z = Bead_add[id].Position.z;
        }
      }
    }

    // write all bead types and molecules
    for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
      BeadType_add[i].Use = true;
      BeadType_add[i].Write = true;
    }
    for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
      MoleculeType_add[i].InVcf = true;
      MoleculeType_add[i].Use = true;
      MoleculeType_add[i].Write = true;
    }
  } //}}}

  // if '-gc' is used, change molecular prototypes to have geometric centre (0,0,0) //{{{
  if (com) {
    for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
      int mt_i_nB = MoleculeType_add[i].nBeads;
      // geometric centre
      Vector centre;
      centre.x = 0;
      centre.y = 0;
      centre.z = 0;

      for (int j = 0; j < MoleculeType_add[i].nBeads; j++) {
        centre.x += prototype[i][j].x;
        centre.y += prototype[i][j].y;
        centre.z += prototype[i][j].z;
      }
      centre.x /= mt_i_nB;
      centre.y /= mt_i_nB;
      centre.z /= mt_i_nB;

      for (int j = 0; j < MoleculeType_add[i].nBeads; j++) {
        prototype[i][j].x -= centre.x;
        prototype[i][j].y -= centre.y;
        prototype[i][j].z -= centre.z;
      }
    }
  } //}}}

  // print what is to be added //{{{
  if (verbose) {
    fprintf(stdout, "\nBEADS AND MOLECULES TO ADD\n");
    VerboseOutput(input_coor, Counts_add, BoxLength, BeadType_add, Bead_add, MoleculeType_add, Molecule_add);
  } //}}}

  // join original and added systems //{{{
  // create structures for new stuff
  struct Counts Counts_new;
  ZeroCounts(&Counts_new);
  struct BeadType *BeadType_new;
  struct MoleculeType *MoleculeType_new;
  struct Bead *Bead_new;
  struct Molecule *Molecule_new;
  int *Index_new;

  Counts_new.Beads = Counts.Beads;
  Counts_new.BeadsInVsf = Counts.BeadsInVsf;
  Counts_new.Bonded = Counts.Bonded + Counts_add.Bonded;
  Counts_new.Unbonded = Counts.Beads - Counts_new.Bonded;
  Counts_new.Molecules = Counts.Molecules + Counts_add.Molecules;
  Counts_new.TypesOfBeads = Counts.TypesOfBeads;
  BeadType_new = calloc(Counts_new.TypesOfBeads, sizeof(struct BeadType));
  // copy original BeadType to _new //{{{
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    strcpy(BeadType_new[i].Name, BeadType[i].Name);
    BeadType_new[i].Number = BeadType[i].Number;
    BeadType_new[i].Charge = BeadType[i].Charge;
    BeadType_new[i].Mass = BeadType[i].Mass;
    BeadType_new[i].Use = BeadType[i].Use;
    BeadType_new[i].Write = BeadType[i].Write;
  } //}}}
  // add new bead types - check if their the same based only on Name //{{{
  for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
    bool new = true;
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      if (strcmp(BeadType_add[i].Name, BeadType[j].Name) == 0) {
        new = false;
        BeadType_new[j].Number += BeadType_add[i].Number;
        break;
      }
    }
    if (new) {
      int type = Counts_new.TypesOfBeads;
      BeadType_new = realloc(BeadType_new, (type+1)*sizeof(struct BeadType));
      strcpy(BeadType_new[type].Name, BeadType_add[i].Name);
      BeadType_new[type].Number = BeadType_add[i].Number;
      BeadType_new[type].Use = BeadType_add[i].Use;
      BeadType_new[type].Write = BeadType_add[i].Write;
      BeadType_new[type].Charge = BeadType_add[i].Charge;
      BeadType_new[type].Mass = BeadType_add[i].Mass;
      Counts_new.TypesOfBeads++;
    }
  } //}}}
  Counts_new.TypesOfMolecules = Counts.TypesOfMolecules;
  MoleculeType_new = calloc(Counts_new.TypesOfMolecules, sizeof(struct MoleculeType));
  // copy original MoleculeType to _new //{{{
  for (int i = 0; i < Counts_new.TypesOfMolecules; i++) {
    strcpy(MoleculeType_new[i].Name, MoleculeType[i].Name);
    MoleculeType_new[i].Number = MoleculeType[i].Number;
    MoleculeType_new[i].nBeads = MoleculeType[i].nBeads;
    MoleculeType_new[i].Bead = calloc(MoleculeType_new[i].nBeads, sizeof(int));
    for (int j = 0; j < MoleculeType_new[i].nBeads; j++) {
      MoleculeType_new[i].Bead[j] = MoleculeType[i].Bead[j];
    }
    MoleculeType_new[i].nBonds = MoleculeType[i].nBonds;
    MoleculeType_new[i].Bond = calloc(MoleculeType_new[i].nBonds, sizeof(int *));
    for (int j = 0; j < MoleculeType_new[i].nBonds; j++) {
      MoleculeType_new[i].Bond[j] = calloc(2, sizeof(int));
      MoleculeType_new[i].Bond[j][0] = MoleculeType[i].Bond[j][0];
      MoleculeType_new[i].Bond[j][1] = MoleculeType[i].Bond[j][1];
    }
    MoleculeType_new[i].nBTypes = MoleculeType[i].nBTypes;
    MoleculeType_new[i].BType = calloc(MoleculeType_new[i].nBTypes, sizeof(int));
    for (int j = 0; j < MoleculeType_new[i].nBTypes; j++) {
      MoleculeType_new[i].BType[j] = MoleculeType[i].BType[j];
    }
    MoleculeType_new[i].Mass = MoleculeType[i].Mass;
    MoleculeType_new[i].InVcf = MoleculeType[i].InVcf;
    MoleculeType_new[i].Use = MoleculeType[i].Use;
    MoleculeType_new[i].Write = MoleculeType[i].Write;
  } //}}}
  // add new molecule types - check if their the same based only on Name //{{{
  for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
    bool new = true;
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (strcmp(MoleculeType_add[i].Name, MoleculeType[j].Name) == 0) {
        new = false;
        MoleculeType_new[j].Number += MoleculeType_add[i].Number;
        break;
      }
    }
    if (new) {
      int type = Counts_new.TypesOfMolecules;
      MoleculeType_new = realloc(MoleculeType_new, (type+1)*sizeof(struct MoleculeType));
      strcpy(MoleculeType_new[type].Name, MoleculeType_add[i].Name);
      MoleculeType_new[type].Number = MoleculeType_add[i].Number;
      MoleculeType_new[type].nBeads = MoleculeType_add[i].nBeads;
      MoleculeType_new[type].Bead = calloc(MoleculeType_new[type].nBeads, sizeof(int));
      for (int j = 0; j < MoleculeType_new[type].nBeads; j++) {
        int old_type = FindBeadType(BeadType_add[MoleculeType_add[i].Bead[j]].Name, Counts_new, BeadType_new);
        MoleculeType_new[type].Bead[j] = old_type;
      }
      MoleculeType_new[type].nBonds = MoleculeType_add[i].nBonds;
      MoleculeType_new[type].Bond = calloc(MoleculeType_new[type].nBonds, sizeof(int *));
      for (int j = 0; j < MoleculeType_new[type].nBonds; j++) {
        MoleculeType_new[type].Bond[j] = calloc(2, sizeof(int));
        MoleculeType_new[type].Bond[j][0] = MoleculeType_add[i].Bond[j][0];
        MoleculeType_new[type].Bond[j][1] = MoleculeType_add[i].Bond[j][1];
      }
      MoleculeType_new[type].nBTypes = MoleculeType_add[i].nBTypes;
      MoleculeType_new[type].BType = calloc(MoleculeType_new[type].nBTypes, sizeof(int));
      for (int j = 0; j < MoleculeType_new[type].nBTypes; j++) {
        int old_type = FindBeadType(BeadType_add[MoleculeType_add[i].BType[j]].Name, Counts_new, BeadType_new);
        MoleculeType_new[type].BType[j] = old_type;
      }
      MoleculeType_new[type].Mass = MoleculeType_add[i].Mass;
      MoleculeType_new[type].InVcf = MoleculeType_add[i].InVcf;
      MoleculeType_new[type].Use = MoleculeType_add[i].Use;
      MoleculeType_new[type].Write = MoleculeType_add[i].Write;
      Counts_new.TypesOfMolecules++;
    }
  } //}}}
  Bead_new = calloc(Counts_new.Beads, sizeof(struct Bead));
  Index_new = calloc(Counts_new.Beads, sizeof(int));
  // copy beads that aren't to be exchanged to the start of Bead_new //{{{
  count = 0; // counts copied beads
  for (int i = 0; i < Counts.Beads; i++) {
    int type = Bead[i].Type;
    if (!BeadType[type].Write && // for '-xb' option
        Bead[i].Molecule == -1) {
      int new_type = FindBeadType(BeadType[type].Name, Counts_new, BeadType_new);
      Bead_new[count].Type = new_type;
      Bead_new[count].Molecule = -1;
      Bead_new[count].Index = count;
      Bead_new[count].Position.x = Bead[i].Position.x;
      Bead_new[count].Position.y = Bead[i].Position.y;
      Bead_new[count].Position.z = Bead[i].Position.z;
      Bead_new[count].Flag = false; // do not rewrite
      Bead_new[count].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[count] = count;
      count++;
    }
  } //}}}
  // copy uncharged unbonded beads beyond the charged Bead_new //{{{
  // here, count == number of beads in the original system not to rewrite
  for (int i = 0; i < Counts.Beads; i++) {
    int type = Bead[i].Type;
    if (BeadType[type].Write && // for '-xb' option
        Bead[i].Molecule == -1) {
      int new_type = FindBeadType(BeadType[type].Name, Counts_new, BeadType_new);
      Bead_new[count].Type = new_type;
      Bead_new[count].Molecule = -1;
      Bead_new[count].Index = count;
      Bead_new[count].Position.x = Bead[i].Position.x;
      Bead_new[count].Position.y = Bead[i].Position.y;
      Bead_new[count].Position.z = Bead[i].Position.z;
      Bead_new[count].Flag = false; // do not rewrite
      Bead_new[count].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[count] = count;
      count++;
    }
    if (count == Counts_new.Unbonded) {
      break;
    }
  } //}}}
  // find beads to exchange //{{{
  count = 0; // counts added beads
  for (int i = 0; i < Counts_new.Unbonded; i++) {
    int type = Bead_new[i].Type;
    if (count != Counts_add.Unbonded &&
        BeadType_new[type].Write && // for '-xb' option
        Bead_new[i].Molecule == -1) {
      BeadType_new[type].Number--;
      int old_type = Bead_add[count].Type;
      int new_type = FindBeadType(BeadType_add[old_type].Name, Counts_new, BeadType_new);
      Bead_new[i].Type = new_type;
      Bead_new[i].Flag = true; // exchange this bead
      count++;
    }
    if (count == Counts_add.Unbonded) {
      break;
    }
  } //}}}
  // copy bonded beads from Bead to Bead_new //{{{
  count = Counts_new.Unbonded;
  for (int i = 0; i < Counts.Beads; i++) {
    int old_type = Bead[i].Type;
    if (Bead[i].Molecule != -1) {
      int new_type = FindBeadType(BeadType[old_type].Name, Counts_new, BeadType_new);
      Bead_new[count].Type = new_type;
      Bead_new[count].Molecule = Bead[i].Molecule;
      Bead_new[count].Index = count;
      Bead_new[count].Position.x = Bead[i].Position.x;
      Bead_new[count].Position.y = Bead[i].Position.y;
      Bead_new[count].Position.z = Bead[i].Position.z;
      Bead_new[count].Flag = false; // do not rewrite
      Bead_new[count].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[count] = count;
      count++;
    }
  } //}}}
  // find original unbonded beads to be exchanged for added bonded beads //{{{
  count = Counts_add.Beads;
  for (int i = (Counts.Beads-1); i >= 0; i--) {
    int type = Bead[i].Type;
    if (BeadType[type].Write && // for '-xb' option
        Bead[i].Molecule == -1 && BeadType[type].Charge == 0) {
      count--;
      int new_count = Counts_new.Beads - (Counts_add.Beads - count);
      int new_type = FindBeadType(BeadType[type].Name, Counts_new, BeadType_new);
      BeadType_new[new_type].Number--;
      int old_type = Bead_add[count].Type;
      new_type = FindBeadType(BeadType_add[old_type].Name, Counts_new, BeadType_new);
      Bead_new[new_count].Type = new_type;
      Bead_new[new_count].Molecule = Counts.Molecules + Bead_add[count].Molecule;
      Bead_new[new_count].Index = new_count;
      Bead_new[new_count].Flag = true; // exchange
      Index_new[new_count] = new_count;
    }
    if (count == Counts_add.Unbonded) {
      break;
    }
  } //}}}
  Molecule_new = calloc(Counts_new.Molecules, sizeof(struct Molecule));
  // copy original molecules to _new struct //{{{
  count = -1;
  for (int i = 0; i < Counts.Molecules; i++) {
    int old_type = Molecule[i].Type;
    int new_type = FindMoleculeType(MoleculeType[old_type].Name, Counts_new, MoleculeType_new);
    Molecule_new[i].Type = new_type;
    Molecule_new[i].Bead = calloc(MoleculeType_new[new_type].nBeads, sizeof(int));
    for (int j = 0; j < MoleculeType_new[new_type].nBeads; j++) {
      while (Bead_new[++count].Molecule != i)
        ;
      Molecule_new[i].Bead[j] = count;
    }
  } //}}}
  // put _add molecules into _new struct //{{{
  count = Counts_new.Unbonded - 1;
  for (int i = 0; i < Counts_add.Molecules; i++) {
    int old_type = Molecule_add[i].Type;
    int new_type = FindMoleculeType(MoleculeType_add[old_type].Name, Counts_new, MoleculeType_new);
    int new_i = i + Counts.Molecules;
    Molecule_new[new_i].Type = new_type;
    Molecule_new[new_i].Bead = calloc(MoleculeType_new[new_type].nBeads, sizeof(int));
    for (int j = 0; j < MoleculeType_new[new_type].nBeads; j++) {
      while (count < Counts_new.Beads && Bead_new[++count].Molecule != new_i)
        ;
      Molecule_new[new_i].Bead[j] = count;
    }
  } //}}}
  //}}}

  // print new system //{{{
  if (verbose) {
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(input_coor, Counts_new, BoxLength, BeadType_new, Bead_new, MoleculeType_new, Molecule_new);
  } //}}}

  // create & fill output vsf file
  WriteVsf(output_vsf, Counts_new, BeadType_new, Bead_new, MoleculeType_new, Molecule_new);

  // square lowest/highest distance, if '-ld' and/or '-hd' options used //{{{
  if (lowest_dist != -1) {
    lowest_dist = SQR(lowest_dist);
  }
  if (highest_dist != -1) {
    highest_dist = SQR(highest_dist);
  } //}}}

  // seed random number generator
  srand(time(0));

  // count unbonded neutral beads //{{{
  int can_be_exchanged = 0;
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    int btype = Bead[i].Type;
    if (Bead[i].Molecule == -1 && BeadType[btype].Charge == 0) {
      can_be_exchanged++;
    }
  }
  // count beads to be added
  int to_be_added = 0;
  for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
    to_be_added += BeadType[i].Number;
  }
  if (Counts_add.Beads > can_be_exchanged) {
    fprintf(stderr, "\nError: there are more beads to be added than can be changed in the original system:\n");
    fprintf(stderr, "     Number of unbonded neutral beads in the original system: %d\n", can_be_exchanged);
    fprintf(stderr, "     Number of beads to be added: %d\n\n", Counts_add.Beads);
    exit(1);
  } //}}}

  // add monomeric beads //{{{
  count = 0;
  for (int i = 0; i < Counts_new.Unbonded; i++) {
    if (Bead_new[i].Flag) { // is this an added bead?
      Vector random;

      if (add_vsf[0] == '\0') { // randomly place monomers if FIELD-like file is used
        double min_dist;
        if (lowest_dist != -1 || highest_dist != -1) { // ||
          do {
            random.x = (double)rand() / ((double)RAND_MAX + 1) * new_box.x + constraint[0].x;
            random.y = (double)rand() / ((double)RAND_MAX + 1) * new_box.y + constraint[0].y;
            random.z = (double)rand() / ((double)RAND_MAX + 1) * new_box.z + constraint[0].z;

            min_dist = SQR(BoxLength.x * 100);
            for (int j = 0; j < Counts.Beads; j++) {
              int btype = Bead[j].Type;
              // TODO: is that thing below true? ...anyway, use Counts_new now
              // j can be added monomeric bead, so it's type can be higher than the number of types
              if (btype < Counts.TypesOfBeads && BeadType[btype].Use) {
                Vector dist;
                dist = Distance(Bead[j].Position, random, BoxLength);
                dist.x = SQR(dist.x) + SQR(dist.y) + SQR(dist.z);
                if (dist.x < min_dist) {
                  min_dist = dist.x;
                }
              }
              if ((lowest_dist != -1 && lowest_dist >= min_dist) ||
                  (highest_dist != -1 && highest_dist <= min_dist)) {
                break;
              }
            }
          } while ((lowest_dist != -1 && lowest_dist >= min_dist) ||
                   (highest_dist != -1 && highest_dist <= min_dist));
        } else {
          random.x = (double)rand() / ((double)RAND_MAX + 1) * new_box.x + constraint[0].x;
          random.y = (double)rand() / ((double)RAND_MAX + 1) * new_box.y + constraint[0].y;
          random.z = (double)rand() / ((double)RAND_MAX + 1) * new_box.z + constraint[0].z;
        }
      } else { // place them according to add_vcf file if '-vtf' option is used
        random.x = Bead_add[count].Position.x;
        random.y = Bead_add[count].Position.y;
        random.z = Bead_add[count].Position.z;
      }

      Bead_new[i].Position.x = random.x;
      Bead_new[i].Position.y = random.y;
      Bead_new[i].Position.z = random.z;
      count++;

      // print number of placed beads? //{{{
      if (!silent && !script) {
        fflush(stdout);
        fprintf(stdout, "\rMonomers placed: %d", count);
      } //}}}
    }
    // stop when all unbonded beads are added
    if (count == Counts_add.Unbonded) {
      break;
    }
  }

  // print total number of placed beads? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Monomer placed: %d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                           ");
      fprintf(stdout, "\rMonomer placed: %d\n", count);
    }
  } //}}}
  //}}}

  // add the molecules //{{{
  count = Counts_new.Beads - Counts_add.Bonded; // index of the first added bonded bead
  for (int i = Counts.Molecules; i < Counts_new.Molecules; i++) {
    int mtype = Molecule_new[i].Type;
    int i_nbeads = MoleculeType_new[mtype].nBeads;

    Vector random;
    // if read from FIELD, find rotated
    Vector rotated[i_nbeads];
    if (add_vsf[0] == '\0') {
      // rotate the prototype molecule randomly (only used when FIELD is used) //{{{
      // random rotation axis
      random.x = (double)rand() / ((double)RAND_MAX) * 2 - 1; // random number <-1,1>
      random.y = (double)rand() / ((double)RAND_MAX) * 2 - 1;
      random.z = (double)rand() / ((double)RAND_MAX) * 2 - 1;
      double dist = sqrt(SQR(random.x)+SQR(random.y)+SQR(random.z));
      random.x /= dist;
      random.y /= dist;
      random.z /= dist;
      // random rotation angle
      double angle = (double)rand() / ((double)RAND_MAX) * PI;
      // create rotation matrix
      struct Tensor {
        Vector x, y, z;
      } rot;
      rot.x.x = cos(angle) + SQR(random.x) * (1 - cos(angle));
      rot.x.y = random.x * random.y * (1 - cos(angle)) - random.z * sin(angle);
      rot.x.z = random.x * random.z * (1 - cos(angle)) + random.y * sin(angle);

      rot.y.x = random.x * random.y * (1 - cos(angle)) + random.z * sin(angle);
      rot.y.y = cos(angle) + SQR(random.y) * (1 - cos(angle));
      rot.y.z = random.y * random.z * (1 - cos(angle)) - random.x * sin(angle);

      rot.z.x = random.x * random.z * (1 - cos(angle)) - random.y * sin(angle);
      rot.z.y = random.y * random.z * (1 - cos(angle)) + random.x * sin(angle);
      rot.z.z = cos(angle) + SQR(random.z) * (1 - cos(angle));
      // transform the prototype molecule (rotation matrix x coordinates)
      int number;
      for (int j = Counts_add.Unbonded; j < Counts_add.Beads; j++) {
        if (Molecule_add[Bead_add[j].Molecule].Type == Molecule_add[i-Counts.Molecules].Type) {
          number = j;
          break;
        }
      }
      for (int j = 0; j < i_nbeads; j++) {
        rotated[j].x = rot.x.x * Bead_add[number].Position.x
                     + rot.x.y * Bead_add[number].Position.y
                     + rot.x.z * Bead_add[number].Position.z;
        rotated[j].y = rot.y.x * Bead_add[number].Position.x
                     + rot.y.y * Bead_add[number].Position.y
                     + rot.y.z * Bead_add[number].Position.z;
        rotated[j].z = rot.z.x * Bead_add[number].Position.x
                     + rot.z.y * Bead_add[number].Position.y
                     + rot.z.z * Bead_add[number].Position.z;
        number++;
      } //}}}

      // first bead's distance from specified bead typtes is checked (only used when FIELD is used) //{{{
      // first bead can have coordinates [0,0,0] or such that the molecule's geometric centre is [0,0,0] (if -gc is used)
      double min_dist;
      if (lowest_dist != -1 || highest_dist != -1) {
        do {
          random.x = (double)rand() / ((double)RAND_MAX + 1) * new_box.x + constraint[0].x;
          random.y = (double)rand() / ((double)RAND_MAX + 1) * new_box.y + constraint[0].y;
          random.z = (double)rand() / ((double)RAND_MAX + 1) * new_box.z + constraint[0].z;

          min_dist = SQR(BoxLength.x * 100);
          for (int j = 0; j < Counts.Beads; j++) {
            int btype_j = Bead[j].Type;
            // j can be added monomeric bead, so it's type can be higher than the number of types
            if (btype_j < Counts.TypesOfBeads && BeadType[btype_j].Use) {
              Vector dist;
              dist = Distance(Bead[j].Position, random, BoxLength);
              dist.x = SQR(dist.x) + SQR(dist.y) + SQR(dist.z);
              if (dist.x < min_dist) {
                min_dist = dist.x;
              }
            }
          }
        } while ((lowest_dist != -1 && lowest_dist >= min_dist) ||
                 (highest_dist != -1 && highest_dist <= min_dist));
      } else { // no '-ld' or '-hd' options
        random.x = (double)rand() / ((double)RAND_MAX + 1) * new_box.x + constraint[0].x;
        random.y = (double)rand() / ((double)RAND_MAX + 1) * new_box.y + constraint[0].y;
        random.z = (double)rand() / ((double)RAND_MAX + 1) * new_box.z + constraint[0].z;
      } //}}}
    }

    for (int j = 0; j < MoleculeType_new[mtype].nBeads; j++) {
      int id = Molecule_add[i-Counts.Molecules].Bead[j];

      Vector r;
      if (add_vsf[0] == '\0') { // random placement if FIELD is used
        r.x = random.x + rotated[j].x;
        r.y = random.y + rotated[j].y;
        r.z = random.z + rotated[j].z;
      } else { // place them according to add_vcf file if '-vtf' option is used
        r.x = Bead_add[id].Position.x;
        r.y = Bead_add[id].Position.y;
        r.z = Bead_add[id].Position.z;
      }

      Bead_new[count].Position.x = r.x;
      Bead_new[count].Position.y = r.y;
      Bead_new[count].Position.z = r.z;
      count++;
    }

    // print number of placed molecules? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rMolecules placed: %d", i-Counts.Molecules+1);
    } //}}}
  }

  // print total number of placed molecules? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Molecules placed: %3d\n", Counts_add.Molecules);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                                             ");
      fprintf(stdout, "\rMolecules placed: %3d\n", Counts_add.Molecules);
    }
  } //}}}
  //}}}

  // open output .vcf file //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  } //}}}

  // open xyz file and write coordinates (if -xyz option is present) //{{{
  FILE *xyz = NULL; // just make sure it's initialized - for compiler error purposes
  if (output_xyz[0] != '\0') {
    // open output .xyz file for reading
    if ((xyz = fopen(output_xyz, "w")) == NULL) {
      ErrorFileOpen(output_xyz, 'w');
      exit(1);
    }

    WriteCoorXYZ(xyz, Counts, BeadType, Bead);
  } //}}}

  // print command, bead type names & box size to output .vcf file //{{{
  // command
  fprintf(out, "# Generated by: AddToSystem ");
  for (int i = 1; i < argc; i++) {
    fprintf(out, " %s", argv[i]);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z); //}}}

  // print coordinates to output .vcf file //{{{
  // write all beads (Write flag was used with '-xb' option)
  for (int i = 0; i < Counts_new.TypesOfBeads; i++) {
    BeadType_new[i].Write = true;
  }
  // write all molecules (basically just to make sure)
  for (int i = 0; i < Counts_new.TypesOfMolecules; i++) {
    MoleculeType_new[i].Write = true;
  }
  WriteCoorIndexed(out, Counts_new, BeadType_new, Bead_new, MoleculeType_new, Molecule_new, stuff); //}}}

  // close output files //{{{
  fclose(out);
  if (output_xyz[0] != '\0') {
    fclose(xyz);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(BeadType_add);
  free(BeadType_new);
  free(Index);
  free(Index_add);
  free(Index_new);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMoleculeType(Counts_new, &MoleculeType_new);
  FreeMoleculeType(Counts_add, &MoleculeType_add);
  FreeMolecule(Counts, &Molecule);
  FreeMolecule(Counts_new, &Molecule_new);
  FreeMolecule(Counts_add, &Molecule_add);
  FreeBead(Counts, &Bead);
  FreeBead(Counts_new, &Bead_new);
  FreeBead(Counts_add, &Bead_add);
  free(stuff);
  free(output_xyz);
  free(input_add);
  free(add_vsf);
  free(add_vcf);
  for (int i = 0; i < Counts_new.TypesOfMolecules; i++) {
    free(prototype[i]);
  }
  free(prototype);
  //}}}

  return 0;
}
