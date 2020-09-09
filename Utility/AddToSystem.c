#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;

    fprintf(ptr, "\
AddToSystem takes an existing system and adds unbonded beads and/or \
molecules to that system. The data for added components is read either \
from a FIELD-like file that has to contain the species and molecule sections \
(the same as for the DL_MESO simulation program) or from a vsf/vcf files \
(-vtf option). In the first case, new molecules and beads are added \
randomly (with several possible constraints), while in the second case, \
the provided coordinates are used. The new beads are either added to the \
original system, increasing the total number of beads, or they replace \
unbonded beads from the original system (--switch option).\n\n");
  }
  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input.vcf> ", cmd);
  fprintf(ptr, "<out.vsf> <out.vcf> <options>\n\n");

  fprintf(ptr, "   <input>             input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <out.vsf>           output structure file (vsf format)\n");
  fprintf(ptr, "   <out.vcf>           output coordinate file (vcf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>        timestep to add new beads to (default: 1)\n");
  fprintf(ptr, "      -xyz <name>      also save coordinates to an xyz file\n");
  fprintf(ptr, "      --switch         replace original beads instead of increasing the total number of beads\n");
  fprintf(ptr, "      -b <x> <y> <z>   side lengths of the new simulation box\n");
  fprintf(ptr, "      --centre         place the original simulation box in the middle of the new one\n");
  fprintf(ptr, "   <options for random placement>\n");
  fprintf(ptr, "      -f <name>        FIELD-like file with molecules to add (default: FIELD)\n");
  fprintf(ptr, "      -ld <float>      specify lowest distance from chosen bead types (default: none)\n");
  fprintf(ptr, "      -hd <float>      specify highest distance from chosen bead types (default: none)\n");
  fprintf(ptr, "      -bt <name(s)>    specify bead types new beads should be far from/near to (default: none)\n");
  fprintf(ptr, "      -cx <num> <num2> constrain x coordinate of randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -cy <num> <num2> constrain y coordinate of randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -cz <num> <num2> constrain z coordinate of randomly added beads to interaval (int,int2)\n");
  fprintf(ptr, "      -gc              use molecule's geometric centre for the distance check instead of its first bead\n");
  fprintf(ptr, "   <options for provided coordinates>\n");
  fprintf(ptr, "      -vtf <vsf> <vcf> use vtf file format instead of FIELD (divided to vsf and vcf files)\n");
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
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "--centre") != 0 &&
        strcmp(argv[i], "--switch") != 0 &&
        strcmp(argv[i], "-xb") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input.vcf> - input coordinate file //{{{
  char input_coor[LINE];
  char *input_vsf = calloc(LINE,sizeof(char));
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  // if vtf, copy to input_vsf
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
  } //}}}

  // <out.vsf> - filename of output vsf file (must end with .vsf) //{{{
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
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // -f <add> - FIELD-like file with molecules to add //{{{
  char *input_add = calloc(LINE, sizeof(char));
  if (FileOption(argc, argv, "-f", &input_add)) {
    exit(1);
  }
  if (input_add[0] == '\0') {
    strcpy(input_add, "FIELD");
  } //}}}

  // -vtf <vsf> <vcf> - FIELD-like file with molecules to add //{{{
  char *add_vsf = calloc(LINE, sizeof(char)),
       *input_coor_add = calloc(LINE, sizeof(char));
  // 1) vsf file
  if (FileOption(argc, argv, "-vtf", &add_vsf)) {
    exit(1);
  }
  // 2) if vsf file exists, look for vcf
  if (add_vsf[0] != '\0') {
    ext = 2;
    strcpy(extension[0], ".vsf");
    strcpy(extension[1], ".vtf");
    if (ErrorExtension(add_vsf, ext, extension)) {
      Help(argv[0], true);
      exit(1);
    }
    if (add_vsf[strlen(add_vsf)-2] == 't') { // if *.vtf file, use it as a vcf as well
      strcpy(input_coor_add, add_vsf);
    } else {
      for (int i = 1; i < argc; i++) { // if *.vsf file, read next file
        if (strcmp(argv[i], "-vtf") == 0) {
          if (argc >= (i+3)) { // is there another cli argument behind '.vsf' file?
            char temp[LINE];
            strcpy(temp, argv[i+1]); // save vsf filename
            strcpy(argv[i+1], argv[i+2]); // copy vcf filename to vsf - required by FileOption()
            ext = 1;
            strcpy(extension[0], ".vcf");
            if (FileOption(argc, argv, "-vtf", &input_coor_add)) {
              exit(1);
            }
            if (ErrorExtension(input_coor_add, ext, extension)) {
              fprintf(stderr, "\033[1;31m");
              fprintf(stderr, "       Wrong second filename!\n");
              fprintf(stderr, "\033[0m");
              Help(argv[0], true);
              exit(1);
            }
            strcpy(argv[i+1], temp); // restore vsf filename
            break;
          } else { // there's no more cli arguments
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: option -vtf is missing second filename!\n");
            fprintf(stderr, "\033[0m");
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
  char *output_xyz = calloc(LINE, sizeof(char));
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
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: if '-ld' and/or '-hd' is used, '-bt' must be specified as well\n\n");
      fprintf(stderr, "\033[0m");
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
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: option \033[1;33m-cx\033[1;31m requires two numeric arguments\n\n");
    fprintf(stderr, "\033[0m");
    Help(argv[0], true);
    exit(1);
  }
  // make sure first number is smaller
  if (range[0] > range[1]) {
    SwapDouble(&range[0], &range[1]);
  }
  VECTOR constraint[2];
  constraint[0].x = range[0];
  constraint[1].x = range[1]; //}}}
  // y direcion //{{{
  test = 2;
  range[0] = range[1] = 0;
  if (MultiDoubleOption(argc, argv, "-cy", &test, range)) {
    exit(1);
  }
  if (test != 2) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: option \033[1;33m-cy\033[1;31m requires two numeric arguments\n\n");
    fprintf(stderr, "\033[0m");
    Help(argv[0], true);
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
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: option \033[1;33m-cz\033[1;31m requires two numeric arguments\n\n");
    fprintf(stderr, "\033[0m");
    Help(argv[0], true);
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

  // switch beads instead of adding beads?
  bool sw = BoolOption(argc, argv, "--switch");

  // define new box size //{{{
  double box_option[100] = {-1};
  if (MultiDoubleOption(argc, argv, "-b", &count, box_option)) {
    exit(1);
  }
  if (count != 3) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m-b\033[1;31m option requires three dimensions of the cuboid simulation box.\n");
    fprintf(stderr, "\033[0m");
    Help(argv[0], true);
    exit(1);
  } //}}}

  // centre old box in the new one? (only relevant if -b or -vtf option used)
  bool centre = BoolOption(argc, argv, "--centre");
  //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // -xb <name(s)> - specify what bead types to exchange //{{{
  if (BeadTypeOption(argc, argv, "-xb", false, Counts, &BeadType)) {
    exit(0);
  }
  // use Write flag to decide which bead types to use
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = BeadType[i].Use;
    BeadType[i].Use = false;
  }
  // if all are 'false', then set 'true' to all uncharged bead types
  bool all_false = true;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Write == true) {
      all_false = false;
      break;
    }
  }
  if (all_false) {
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      if (BeadType[i].Charge == 0) {
        BeadType[i].Write = true;
      }
    }
  } //}}}

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

  // read original box length and assign new one if -b is used //{{{
  VECTOR BoxLength = GetPBC(vcf, input_coor);
  VECTOR BoxLength_new = BoxLength;
  if (box_option[0] != -1) {
    BoxLength_new.x = box_option[0];
    BoxLength_new.y = box_option[1];
    BoxLength_new.z = box_option[2];
  } //}}}

  // print original system //{{{
  if (verbose) {
    fprintf(stdout, "\nORIGINAL SYSTEM\n");
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    if (start > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", start);
    }
  } //}}}

  // define box size using constraints (-c{x,y,z} options)//{{{
  VECTOR new_box;
  if (constraint[1].x != 0) {
    new_box.x = constraint[1].x - constraint[0].x;
  } else {
    new_box.x = BoxLength_new.x;
  }
  if (constraint[1].y != 0) {
    new_box.y = constraint[1].y - constraint[0].y;
  } else {
    new_box.y = BoxLength_new.y;
  }
  if (constraint[1].z != 0) {
    new_box.z = constraint[1].z - constraint[0].z;
  } else {
    new_box.z = BoxLength_new.z;
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
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
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
    ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);
  } else {
    fprintf(stderr, "\033[1;33m");
    fprintf(stderr, "\nWarning: using last step in \033[1;36m%s\033[1;33m (\033[1;36m%d\033[1;33m)\n", input_coor, count);
    fprintf(stderr, "\033[0m");
  }
  fclose(vcf); //}}}

  // create structures for added stuff //{{{
  COUNTS Counts_add = InitCounts;
  MOLECULE *Molecule_add;
  MOLECULETYPE *MoleculeType_add;
  BEADTYPE *BeadType_add;
  BEAD *Bead_add;
//struct VECTOR **prototype;
  int *Index_add;
  PARAMS *bond_type; // information about bond types
  PARAMS *angle_type; // information about angle types
  //}}}

  if (add_vsf[0] == '\0') { // read stuff to be added from FIELD
    ReadField(input_add, '\0', &Counts_add, &BeadType_add, &Bead_add, &Index_add, &MoleculeType_add, &Molecule_add, &bond_type, &angle_type);
  } else { // read stuff to add from vtf file(s) ('-vtf' option) //{{{
    bool indexed_add = ReadStructure(add_vsf, input_coor_add, &Counts_add, &BeadType_add, &Bead_add, &Index_add, &MoleculeType_add, &Molecule_add);
    // open input coordinate file
    if ((vcf = fopen(input_coor_add, "r")) == NULL) {
      ErrorFileOpen(input_coor_add, 'r');
      exit(1);
    }
    VECTOR BoxLength_add = GetPBC(vcf, input_coor_add);
    // determine new box size if -b isn't used as the larger of the dimensions
    // from original and to be added systems
    if (box_option[0] != -1) {
      if (BoxLength_add.x > BoxLength_new.x) {
        BoxLength_new.x = BoxLength_add.x;
      }
      if (BoxLength_add.y > BoxLength_new.y) {
        BoxLength_new.y = BoxLength_add.y;
      }
      if (BoxLength_add.z > BoxLength_new.z) {
        BoxLength_new.z = BoxLength_add.z;
      }
    }
    ReadCoordinates(indexed_add, input_coor_add, vcf, Counts_add, Index_add, &Bead_add, &stuff);
    fclose(vcf);
  } //}}}

  // move the original system to the centre of the new box if --centre is used //{{{
  if (centre) {
    VECTOR move;
    move.x = (BoxLength_new.x - BoxLength.x) / 2;
    move.y = (BoxLength_new.y - BoxLength.y) / 2;
    move.z = (BoxLength_new.z - BoxLength.z) / 2;
    for (int i = 0; i < Counts.Beads; i++) {
      Bead[i].Position.x += move.x;
      Bead[i].Position.y += move.y;
      Bead[i].Position.z += move.z;
    }
  } //}}}

  // check number of exchangeable beads //{{{
  int can_be_exchanged = 0;
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    int btype = Bead[i].Type;
    if (Bead[i].Molecule == -1 && BeadType[btype].Write) {
      can_be_exchanged++;
    }
  }
  // count beads to be added
  if (sw && Counts_add.Beads > can_be_exchanged) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: insufficient beads to exchange\n");
    fprintf(stderr, "     Number of exchangeable beads in the original system: \033[1;33m%d\033[1;31m\n", can_be_exchanged);
    fprintf(stderr, "     Number of beads to be added: \033[1;33m%d\033[1;31m\n\n", Counts_add.Beads);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // if '-gc' is used, change molecular prototypes to have geometric centre (0,0,0) //{{{
  if (com) {
    for (int i = 0; i < Counts_add.Molecules; i++) {
      int mtype = Molecule_add[i].Type;
      VECTOR geom_centre;
      geom_centre.x = 0;
      geom_centre.y = 0;
      geom_centre.z = 0;

      for (int j = 0; j < MoleculeType_add[mtype].nBeads; j++) {
        int id = Molecule_add[i].Bead[j];
        geom_centre.x += Bead_add[id].Position.x;
        geom_centre.y += Bead_add[id].Position.y;
        geom_centre.z += Bead_add[id].Position.z;
      }
      geom_centre.x /= MoleculeType_add[mtype].nBeads;
      geom_centre.y /= MoleculeType_add[mtype].nBeads;
      geom_centre.z /= MoleculeType_add[mtype].nBeads;

      for (int j = 0; j < MoleculeType_add[mtype].nBeads; j++) {
        int id = Molecule_add[i].Bead[j];
        Bead_add[id].Position.x -= geom_centre.x;
        Bead_add[id].Position.y -= geom_centre.y;
        Bead_add[id].Position.z -= geom_centre.z;
      }
    }
  } //}}}

  // print what is to be added //{{{
  if (verbose) {
    fprintf(stdout, "\nBEADS AND MOLECULES TO ADD\n");
    VerboseOutput(input_coor, Counts_add, BoxLength, BeadType_add, Bead_add, MoleculeType_add, Molecule_add);
  } //}}}

  // decide which beads to exchange (i.e., give them Bead[].Flag = true); has effect only if --switch is used //{{{
  // zeroize Bead[].Flag //{{{
  for (int i = 0; i < Counts.Beads; i++) {
    Bead[i].Flag = false;
  } //}}}
  count = 0; // counts bead in the original Bead[] struct
  for (int i = 0; i < Counts_add.Beads; i++) {
    for (; count < Counts.Beads; count++) {
      int type = Bead[count].Type;
      if (BeadType[type].Write && Bead[count].Molecule == -1) {
        Bead[count].Flag = true; // exchange bead 'count'
        break;
      }
    }
    count++; // loop didn't update count because of the break
  } //}}}

  // create structures for new system //{{{
  COUNTS Counts_new = InitCounts;
  BEADTYPE *BeadType_new;
  MOLECULETYPE *MoleculeType_new;
  BEAD *Bead_new;
  MOLECULE *Molecule_new;
  int *Index_new; //}}}

  // join original and added systems (depending on '--switch' mode)
  if (sw) { // switch old beads for new ones? //{{{
    Counts_new.Beads = Counts.Beads;
    Counts_new.BeadsInVsf = Counts.BeadsInVsf;
    Counts_new.Bonded = Counts.Bonded + Counts_add.Bonded;
    Counts_new.Unbonded = Counts.Beads - Counts_new.Bonded;
    Counts_new.TypesOfBonds = Counts_add.TypesOfBonds;
    Counts_new.TypesOfAngles = Counts_add.TypesOfAngles;
    Counts_new.Molecules = Counts.Molecules + Counts_add.Molecules;
    // fill BeadType struct for the new system
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
    // add new bead types - check if they're the same based only on Name //{{{
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
    // fill MoleculeType struct for the new system
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
        MoleculeType_new[i].Bond[j] = calloc(3, sizeof(int));
        MoleculeType_new[i].Bond[j][0] = MoleculeType[i].Bond[j][0];
        MoleculeType_new[i].Bond[j][1] = MoleculeType[i].Bond[j][1];
        MoleculeType_new[i].Bond[j][2] = MoleculeType[i].Bond[j][2];
      }
      MoleculeType_new[i].nAngles = 0;
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
          int old_type = MoleculeType_add[i].Bead[j];
          int btype = FindBeadType(BeadType_add[old_type].Name, Counts_new, BeadType_new);
          MoleculeType_new[type].Bead[j] = btype;
        }
        MoleculeType_new[type].nBonds = MoleculeType_add[i].nBonds;
        MoleculeType_new[type].Bond = calloc(MoleculeType_new[type].nBonds, sizeof(int *));
        for (int j = 0; j < MoleculeType_new[type].nBonds; j++) {
          MoleculeType_new[type].Bond[j] = calloc(3, sizeof(int));
          MoleculeType_new[type].Bond[j][0] = MoleculeType_add[i].Bond[j][0];
          MoleculeType_new[type].Bond[j][1] = MoleculeType_add[i].Bond[j][1];
          MoleculeType_new[type].Bond[j][2] = MoleculeType_add[i].Bond[j][2];
        }
        MoleculeType_new[type].nAngles = MoleculeType_add[i].nAngles;
        MoleculeType_new[type].Angle = calloc(MoleculeType_new[type].nAngles, sizeof(int *));
        for (int j = 0; j < MoleculeType_new[type].nAngles; j++) {
          MoleculeType_new[type].Angle[j] = calloc(4, sizeof(int));
          MoleculeType_new[type].Angle[j][0] = MoleculeType_add[i].Angle[j][0];
          MoleculeType_new[type].Angle[j][1] = MoleculeType_add[i].Angle[j][1];
          MoleculeType_new[type].Angle[j][2] = MoleculeType_add[i].Angle[j][2];
          MoleculeType_new[type].Angle[j][3] = MoleculeType_add[i].Angle[j][3];
        }
        MoleculeType_new[type].nBTypes = MoleculeType_add[i].nBTypes;
        MoleculeType_new[type].BType = calloc(MoleculeType_new[type].nBTypes, sizeof(int));
        for (int j = 0; j < MoleculeType_new[type].nBTypes; j++) {
          int old_type = MoleculeType_add[i].BType[j];
          int btype = FindBeadType(BeadType_add[old_type].Name, Counts_new, BeadType_new);
          MoleculeType_new[type].BType[j] = btype;
        }
        MoleculeType_new[type].Mass = MoleculeType_add[i].Mass;
        MoleculeType_new[type].InVcf = MoleculeType_add[i].InVcf;
        MoleculeType_new[type].Use = MoleculeType_add[i].Use;
        MoleculeType_new[type].Write = MoleculeType_add[i].Write;
        Counts_new.TypesOfMolecules++;
      }
    } //}}}
    // fill Bead struct for the new system
    Bead_new = calloc(Counts_new.Beads, sizeof(struct Bead));
    Index_new = calloc(Counts_new.Beads, sizeof(int));
    // copy unbonded beads that aren't to be exchanged to the start of Bead_new //{{{
    // TODO: assumes unbonded beads are before bonded beads
    count = 0; // counts copied beads
    for (int i = 0; i < Counts.Unbonded; i++) {
      // first, copy only beads of the type that's not to be exchange
      if (!Bead[i].Flag) {
        Bead_new[count].Type = Bead[i].Type;
        Bead_new[count].Molecule = -1;
        Bead_new[count].Index = count;
        Bead_new[count].Position.x = Bead[i].Position.x;
        Bead_new[count].Position.y = Bead[i].Position.y;
        Bead_new[count].Position.z = Bead[i].Position.z;
        Bead_new[count].Flag = false; // do not rewrite, obviously
        Bead_new[count].Aggregate = calloc(1, sizeof(int)); // just to free later
        Index_new[count] = count;
        count++;
      }
    }
    // count ended at <number of unbonded original beads> - <added beads> //}}}
    // put unbonded beads to be added beyond the unchanged unbonded beads //{{{
    count = Counts.Unbonded - Counts_add.Beads; // just to be sure
    for (int i = 0; i < Counts_add.Unbonded; i++) {
      int type = Bead_add[i].Type;
      int new_type = FindBeadType(BeadType_add[type].Name, Counts_new, BeadType_new);
      Bead_new[count].Type = new_type;
      Bead_new[count].Molecule = -1;
      Bead_new[count].Index = count;
      Bead_new[count].Position.x = Bead_add[i].Position.x;
      Bead_new[count].Position.y = Bead_add[i].Position.y;
      Bead_new[count].Position.z = Bead_add[i].Position.z;
      Bead_new[count].Flag = true; // coordinates to be rewritten
      Bead_new[count].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[count] = count;
      count++; // use count to make it consistent & easy to read
    } //}}}
    // copy the original bonded beads //{{{
    count = Counts_new.Unbonded;
    for (int i = Counts.Unbonded; i < Counts.Beads; i++) {
      Bead_new[count].Type = Bead[i].Type;
      Bead_new[count].Molecule = Bead[i].Molecule;
      Bead_new[count].Index = count;
      Bead_new[count].Position.x = Bead[i].Position.x;
      Bead_new[count].Position.y = Bead[i].Position.y;
      Bead_new[count].Position.z = Bead[i].Position.z;
      Bead_new[count].Flag = false; // coordinates to be rewritten
      Bead_new[count].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[count] = count;
      count++;
    } //}}}
    // put bonded beads to be added at the very end //{{{
    count = Counts_new.Unbonded + Counts.Bonded;
    for (int i = Counts_add.Unbonded; i < Counts_add.Beads; i++) {
      int type = Bead_add[i].Type;
      int new_type = FindBeadType(BeadType_add[type].Name, Counts_new, BeadType_new);
      Bead_new[count].Type = new_type;
      Bead_new[count].Molecule = Bead_add[i].Molecule + Counts.Molecules;
      Bead_new[count].Index = count;
      Bead_new[count].Position.x = Bead_add[i].Position.x;
      Bead_new[count].Position.y = Bead_add[i].Position.y;
      Bead_new[count].Position.z = Bead_add[i].Position.z;
      Bead_new[count].Flag = true; // coordinates to be rewritten
      Bead_new[count].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[count] = count;
      count++; // use count to make it consistent & easy to read
    } //}}}
    // alocate new molecule struct
    Molecule_new = calloc(Counts_new.Molecules, sizeof(struct Molecule));
    // copy original molecules to _new struct //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int type = Molecule[i].Type;
      Molecule_new[i].Type = type;
      Molecule_new[i].Bead = calloc(MoleculeType_new[type].nBeads, sizeof(int));
      for (int j = 0; j < MoleculeType_new[type].nBeads; j++) {
        Molecule_new[i].Bead[j] = Molecule[i].Bead[j] - Counts_add.Bonded;
      }
    } //}}}
    // put _add molecules into _new struct //{{{
    count = Counts_new.Beads - Counts_add.Bonded;
    for (int i = 0; i < Counts_add.Molecules; i++) {
      int add_type = Molecule_add[i].Type;
      int new_type = FindMoleculeType(MoleculeType_add[add_type].Name, Counts_new, MoleculeType_new);
      int new_i = Counts.Molecules + i;
      Molecule_new[new_i].Type = new_type;
      Molecule_new[new_i].Bead = calloc(MoleculeType_new[new_type].nBeads, sizeof(int));
      for (int j = 0; j < MoleculeType_new[new_type].nBeads; j++) {
        Molecule_new[new_i].Bead[j] = count;
        count++;
      }
    } //}}}
    //}}}
  } else { // or add beads to the system? //{{{
    Counts_new.Beads = Counts.Beads + Counts_add.Beads;
    Counts_new.BeadsInVsf = Counts.BeadsInVsf + Counts_add.BeadsInVsf;
    Counts_new.Bonded = Counts.Bonded + Counts_add.Bonded;
    Counts_new.Unbonded = Counts_new.Beads - Counts_new.Bonded;
    Counts_new.TypesOfBonds = Counts_add.TypesOfBonds;
    Counts_new.TypesOfAngles = Counts_add.TypesOfAngles;
    Counts_new.Molecules = Counts.Molecules + Counts_add.Molecules;
    // fill BeadType struct for the new system
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
    // add new bead types - check if they're the same based only on Name //{{{
    for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
      bool new = true;
      for (int j = 0; j < Counts.TypesOfBeads; j++) {
        if (strcmp(BeadType_add[i].Name, BeadType[j].Name) == 0) {
          new = false;
          BeadType_new[j].Number += BeadType_add[i].Number; // increase old type's number of beads
          break;
        }
      }
      if (new) { // create new type
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
    // fill MoleculeType struct for the new system
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
        MoleculeType_new[i].Bond[j] = calloc(3, sizeof(int));
        MoleculeType_new[i].Bond[j][0] = MoleculeType[i].Bond[j][0];
        MoleculeType_new[i].Bond[j][1] = MoleculeType[i].Bond[j][1];
        MoleculeType_new[i].Bond[j][2] = MoleculeType[i].Bond[j][2];
      }
      MoleculeType_new[i].nAngles = 0;
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
          int old_type = MoleculeType_add[i].Bead[j];
          int btype = FindBeadType(BeadType_add[old_type].Name, Counts_new, BeadType_new);
          MoleculeType_new[type].Bead[j] = btype;
        }
        MoleculeType_new[type].nBonds = MoleculeType_add[i].nBonds;
        MoleculeType_new[type].Bond = calloc(MoleculeType_new[type].nBonds, sizeof(int *));
        for (int j = 0; j < MoleculeType_new[type].nBonds; j++) {
          MoleculeType_new[type].Bond[j] = calloc(3, sizeof(int));
          MoleculeType_new[type].Bond[j][0] = MoleculeType_add[i].Bond[j][0];
          MoleculeType_new[type].Bond[j][1] = MoleculeType_add[i].Bond[j][1];
          MoleculeType_new[type].Bond[j][2] = MoleculeType_add[i].Bond[j][2];
        }
        MoleculeType_new[type].nAngles = MoleculeType_add[i].nAngles;
        MoleculeType_new[type].Angle = calloc(MoleculeType_new[type].nAngles, sizeof(int *));
        for (int j = 0; j < MoleculeType_new[type].nAngles; j++) {
          MoleculeType_new[type].Angle[j] = calloc(4, sizeof(int));
          MoleculeType_new[type].Angle[j][0] = MoleculeType_add[i].Angle[j][0];
          MoleculeType_new[type].Angle[j][1] = MoleculeType_add[i].Angle[j][1];
          MoleculeType_new[type].Angle[j][2] = MoleculeType_add[i].Angle[j][2];
          MoleculeType_new[type].Angle[j][3] = MoleculeType_add[i].Angle[j][3];
        }
        MoleculeType_new[type].nBTypes = MoleculeType_add[i].nBTypes;
        MoleculeType_new[type].BType = calloc(MoleculeType_new[type].nBTypes, sizeof(int));
        for (int j = 0; j < MoleculeType_new[type].nBTypes; j++) {
          int old_type = MoleculeType_add[i].BType[j];
          int btype = FindBeadType(BeadType_add[old_type].Name, Counts_new, BeadType_new);
          MoleculeType_new[type].BType[j] = btype;
        }
        MoleculeType_new[type].Mass = MoleculeType_add[i].Mass;
        MoleculeType_new[type].InVcf = MoleculeType_add[i].InVcf;
        MoleculeType_new[type].Use = MoleculeType_add[i].Use;
        MoleculeType_new[type].Write = MoleculeType_add[i].Write;
        Counts_new.TypesOfMolecules++;
      }
    } //}}}
    // fill Bead struct for the new system
    Bead_new = calloc(Counts_new.Beads, sizeof(struct Bead));
    Index_new = calloc(Counts_new.Beads, sizeof(int));
    // copy original unbonded beads to the start of Bead_new //{{{
    // TODO: assumes unbonded beads are before bonded beads
    for (int i = 0; i < Counts.Unbonded; i++) {
      Bead_new[i].Type = Bead[i].Type;
      Bead_new[i].Molecule = -1;
      Bead_new[i].Index = i;
      Bead_new[i].Position.x = Bead[i].Position.x;
      Bead_new[i].Position.y = Bead[i].Position.y;
      Bead_new[i].Position.z = Bead[i].Position.z;
      Bead_new[i].Flag = false; // do not rewrite, obviously
      Bead_new[i].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[i] = i;
    } //}}}
    // put unbonded beads to be added beyond the original unbonded beads //{{{
    for (int i = Counts.Unbonded; i < Counts_new.Unbonded; i++) {
      int id_add = i - Counts.Unbonded;
      int type = Bead_add[id_add].Type;
      int new_type = FindBeadType(BeadType_add[type].Name, Counts_new, BeadType_new);
      Bead_new[i].Type = new_type;
      Bead_new[i].Molecule = -1;
      Bead_new[i].Index = i;
      Bead_new[i].Position.x = Bead_add[id_add].Position.x;
      Bead_new[i].Position.y = Bead_add[id_add].Position.y;
      Bead_new[i].Position.z = Bead_add[id_add].Position.z;
      Bead_new[i].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[i] = i;
    } //}}}
    // copy the original bonded beads //{{{
    for (int i = Counts.Unbonded; i < Counts.Beads; i++) {
      int id = Counts_new.Unbonded + i - Counts.Unbonded; // goes from Counts_new.Unbonded to (Counts_new.Unbonded+Counts.Bonded)
      Bead_new[id].Type = Bead[i].Type;
      Bead_new[id].Molecule = Bead[i].Molecule;
      Bead_new[id].Index = id;
      Bead_new[id].Position.x = Bead[i].Position.x;
      Bead_new[id].Position.y = Bead[i].Position.y;
      Bead_new[id].Position.z = Bead[i].Position.z;
      Bead_new[id].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[id] = id;
    } //}}}
    // put bonded beads to be added at the very end //{{{
    for (int i = Counts_add.Unbonded; i < Counts_add.Beads; i++) {
      int type = Bead_add[i].Type;
      int new_type = FindBeadType(BeadType_add[type].Name, Counts_new, BeadType_new);
      int id = Counts_new.Beads - Counts_add.Beads + i;
      Bead_new[id].Type = new_type;
      Bead_new[id].Molecule = Bead_add[i].Molecule + Counts.Molecules;
      Bead_new[id].Index = id;
      Bead_new[id].Position.x = Bead_add[i].Position.x;
      Bead_new[id].Position.y = Bead_add[i].Position.y;
      Bead_new[id].Position.z = Bead_add[i].Position.z;
      Bead_new[id].Aggregate = calloc(1, sizeof(int)); // just to free later
      Index_new[id] = id;
    } //}}}
    // alocate new molecule struct
    Molecule_new = calloc(Counts_new.Molecules, sizeof(struct Molecule));
    // copy original molecules to _new struct //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int type = Molecule[i].Type;
      Molecule_new[i].Type = type;
      Molecule_new[i].Bead = calloc(MoleculeType_new[type].nBeads, sizeof(int));
      for (int j = 0; j < MoleculeType_new[type].nBeads; j++) {
        Molecule_new[i].Bead[j] = Molecule[i].Bead[j] + Counts_add.Unbonded;
      }
    } //}}}
    // put _add molecules into _new struct //{{{
    count = Counts_new.Beads - Counts_add.Bonded;
    for (int i = 0; i < Counts_add.Molecules; i++) {
      int add_type = Molecule_add[i].Type;
      int new_type = FindMoleculeType(MoleculeType_add[add_type].Name, Counts_new, MoleculeType_new);
      int new_i = Counts.Molecules + i;
      Molecule_new[new_i].Type = new_type;
      Molecule_new[new_i].Bead = calloc(MoleculeType_new[new_type].nBeads, sizeof(int));
      for (int j = 0; j < MoleculeType_new[new_type].nBeads; j++) {
        Molecule_new[new_i].Bead[j] = count;
        count++;
      }
    } //}}}
  } //}}}

  // print new system //{{{
  if (verbose) {
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(input_coor, Counts_new, BoxLength_new, BeadType_new, Bead_new, MoleculeType_new, Molecule_new);
  } //}}}
  // create & fill output vsf file
  WriteVsf(output_vsf, Counts_new, BeadType_new, Bead_new, MoleculeType_new, Molecule_new, false);

  // square lowest/highest distance, if '-ld' and/or '-hd' options used //{{{
  if (lowest_dist != -1) {
    lowest_dist = SQR(lowest_dist);
  }
  if (highest_dist != -1) {
    highest_dist = SQR(highest_dist);
  } //}}}

  // seed random number generator
  srand(time(0));

  // add beads randomly if FIELD-like file is used //{{{
  if (add_vsf[0] == '\0') {
    count = 0;
    // add monomeric beads //{{{
    for (int i = 0; i < Counts_add.Unbonded; i++) {
      VECTOR random;
      if (lowest_dist != -1 || highest_dist != -1) {
        double min_dist;
        do {
          random.x = (double)rand() / ((double)RAND_MAX + 1) * new_box.x + constraint[0].x;
          random.y = (double)rand() / ((double)RAND_MAX + 1) * new_box.y + constraint[0].y;
          random.z = (double)rand() / ((double)RAND_MAX + 1) * new_box.z + constraint[0].z;

          min_dist = SQR(BoxLength.x * 100);
          for (int j = 0; j < Counts.Beads; j++) {
            int btype = Bead[j].Type;
            // j can be added monomeric bead, so it's type can be higher than the number of types
            if (btype < Counts.TypesOfBeads && BeadType[btype].Use) {
              VECTOR dist;
              dist = Distance(Bead[j].Position, random, BoxLength);
              dist.x = SQR(dist.x) + SQR(dist.y) + SQR(dist.z);
              if (dist.x < min_dist) {
                min_dist = dist.x;
              }
            }
          }
        } while ((lowest_dist != -1 && lowest_dist >= min_dist) ||
                 (highest_dist != -1 && highest_dist <= min_dist));
      } else {
        random.x = (double)rand() / ((double)RAND_MAX + 1) * new_box.x + constraint[0].x;
        random.y = (double)rand() / ((double)RAND_MAX + 1) * new_box.y + constraint[0].y;
        random.z = (double)rand() / ((double)RAND_MAX + 1) * new_box.z + constraint[0].z;
      }

      // determine index of the added bead
      int id;
      if (!sw) { // added beads (no --switch option)
        id = Counts.Unbonded + i;
      } else { // switched beds (--switch option)
        for (int j = count; j < Counts_new.Unbonded; j++) {
          if (Bead_new[j].Flag) { // is this an original bead to be exchanged?
            id = j;
            Bead_new[j].Flag = false; // just exchanged (only pro forma)
            count = j + 1;
            break;
          }
        }
      }

      // add the new coordinate
      Bead_new[id].Position.x = random.x;
      Bead_new[id].Position.y = random.y;
      Bead_new[id].Position.z = random.z;

      // print number of placed beads? //{{{
      if (!silent && !script) {
        fflush(stdout);
        fprintf(stdout, "\rMonomers placed: %d", i+1);
      } //}}}
    } //}}}
    // print total number of placed beads? //{{{
    if (!silent) {
      if (script) {
        fprintf(stdout, "Monomer placed: %d\n", Counts_add.Unbonded);
      } else {
        fflush(stdout);
        fprintf(stdout, "\r                           ");
        fprintf(stdout, "\rMonomer placed: %d\n", Counts_add.Unbonded);
      }
    } //}}}
    // add molecules //{{{
    // doesn't depend on --switch option as it's determined by the Molecule_new
    // array established earlier
    count = 0;
    for (int i = Counts.Molecules; i < Counts_new.Molecules; i++) {
      int mtype = Molecule_new[i].Type;

      VECTOR rotated[MoleculeType_new[mtype].nBeads];
      VECTOR random = {0};

      // rotate the molecule randomly //{{{
      // random rotation axis
      random.x = (double)rand() / ((double)RAND_MAX) * 2 - 1; // random number <-1,1>
      random.y = (double)rand() / ((double)RAND_MAX) * 2 - 1;
      random.z = (double)rand() / ((double)RAND_MAX) * 2 - 1;
      double dist = Length(random);
      random.x /= dist;
      random.y /= dist;
      random.z /= dist;
      // random rotation angle
      double angle = (double)rand() / ((double)RAND_MAX) * PI;
      // create rotation matrix
      struct Tensor {
        VECTOR x, y, z;
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
      // transform the prototype molecule (rotation matrix * coordinates)
      for (int j = 0; j < MoleculeType_new[mtype].nBeads; j++) {
        int id = Molecule_new[i].Bead[j];
        rotated[j].x = rot.x.x * Bead_new[id].Position.x
                     + rot.x.y * Bead_new[id].Position.y
                     + rot.x.z * Bead_new[id].Position.z;
        rotated[j].y = rot.y.x * Bead_new[id].Position.x
                     + rot.y.y * Bead_new[id].Position.y
                     + rot.y.z * Bead_new[id].Position.z;
        rotated[j].z = rot.z.x * Bead_new[id].Position.x
                     + rot.z.y * Bead_new[id].Position.y
                     + rot.z.z * Bead_new[id].Position.z;
      } //}}}

      // first bead's distance from specified bead typtes is checked //{{{
      // first bead can have coordinates [0,0,0] or such that the molecule's geometric centre is [0,0,0] (if -gc is used)
      if (lowest_dist != -1 || highest_dist != -1) {
        double min_dist;
        do {
          random.x = (double)rand() / ((double)RAND_MAX + 1) * new_box.x + constraint[0].x;
          random.y = (double)rand() / ((double)RAND_MAX + 1) * new_box.y + constraint[0].y;
          random.z = (double)rand() / ((double)RAND_MAX + 1) * new_box.z + constraint[0].z;

          min_dist = SQR(BoxLength.x * 100);
          for (int j = 0; j < Counts.Beads; j++) {
            int btype_j = Bead[j].Type;
            // j can be added monomeric bead, so it's type can be higher than the number of types
            if (btype_j < Counts.TypesOfBeads && BeadType[btype_j].Use) {
              dist = Length(Distance(Bead[j].Position, random, BoxLength));
              if (dist < min_dist) {
                min_dist = dist;
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

      // place the molecule //{{{
      for (int j = 0; j < MoleculeType_new[mtype].nBeads; j++) {
        int id = Molecule_new[i].Bead[j];
        Bead_new[id].Position.x = random.x + rotated[j].x;
        Bead_new[id].Position.y = random.y + rotated[j].y;
        Bead_new[id].Position.z = random.z + rotated[j].z;
      } //}}}

      // print number of placed molecules? //{{{
      if (!silent && !script) {
        fflush(stdout);
        fprintf(stdout, "\rMolecules placed: %d", i-Counts.Molecules+1);
      } //}}}
    } //}}}
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
  } //}}}

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
  fprintf(out, "# Generated by:");
  PrintCommand(out, argc, argv);
  fprintf(out, "# AnalysisTools version %s; https://github.com/KaGaSi/AnalysisTools/releases\n", VERSION);

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength_new.x, BoxLength_new.y, BoxLength_new.z); //}}}

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
  free(input_coor_add);
  free(angle_type);
  free(bond_type);
  //}}}

  return 0;
}
