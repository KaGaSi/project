#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

//TODO: add option to use centre of mass to decide on distance

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input add> ", cmd);
  fprintf(stderr, "<out.vcf> <out.vsf> <options>\n\n");

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <input add>       FIELD-like file with molecules to add\n");
  fprintf(stderr, "   <out.vcf>         output coordinate file (vcf format)\n");
  fprintf(stderr, "   <out.vsf>         output structure file (vsf format)\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -w             wrap coordinates (i.e., apply pbc)\n");
  fprintf(stderr, "      -st <start>    number of timestep to start from\n");
  fprintf(stderr, "      -xyz <name>    output xyz file\n");
  fprintf(stderr, "      -ld <float>    specify lowest distance from chosen bead types (default: none)\n");
  fprintf(stderr, "      -hd <float>    specify highest distance from chosen bead types (default: none)\n");
  fprintf(stderr, "      -bt <name(s)>  specify bead types new beads should be far from/near to (default: none)\n");
  CommonHelp(1);
} //}}}

// WriteVsf() //{{{
/*
 * Function creating `.vsf` structure file for use in conjunction with
 * `.vcf` coordinate file for better visualisation via VMD program.
 */
void WriteVsf(char *input_vsf, Counts Counts, BeadType *BeadType, Bead *Bead,
              MoleculeType *MoleculeType, Molecule *Molecule) {

  // opten structure file //{{{
  FILE *fw;
  if ((fw = fopen(input_vsf, "w")) == NULL) {
    ErrorFileOpen(input_vsf, 'w');
    exit(1);
  } //}}}

  // find most common type of bead and make it default //{{{
  int type_def = 0, count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Number >= count) {
      count = BeadType[i].Number;
      type_def = i;
    }
  } //}}}

  // print default bead type //{{{
  fprintf(fw, "atom default name %8s ", BeadType[type_def].Name);
  fprintf(fw, "mass %4.2f ", BeadType[type_def].Mass);
  fprintf(fw, "charge %5.2f\n", BeadType[type_def].Charge); //}}}

  // print beads //{{{
  for (int i = 0; i < Counts.BeadsInVsf; i++) {

    // don't print beads with type 'type_def'
    if (Bead[i].Type != type_def) {
      fprintf(fw, "atom %7d ", i);
      fprintf(fw, "name %8s ", BeadType[Bead[i].Type].Name);
      fprintf(fw, "mass %lf ", BeadType[Bead[i].Type].Mass);
      fprintf(fw, "charge %lf", BeadType[Bead[i].Type].Charge);
      if (Bead[i].Molecule != -1) {
        fprintf(fw, " resname %10s ", MoleculeType[Molecule[Bead[i].Molecule].Type].Name);
        fprintf(fw, "resid %5d", Bead[i].Molecule+1);
      }
      putc('\n', fw);
    }
  } //}}}

  // print bonds //{{{
  putc('\n', fw);
  for (int i = 0; i < Counts.Molecules; i++) {
    fprintf(fw, "# resid %d\n", i+1); // in VMD resid start with 1
    int mol_type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[mol_type].nBonds; j++) {
      fprintf(fw, "bond %6d: %6d\n", Molecule[i].Bead[MoleculeType[mol_type].Bond[j][0]],
                                     Molecule[i].Bead[MoleculeType[mol_type].Bond[j][1]]);
    }
  } //}}}

  // close structure file
  fclose(fw);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
SelectedVcf creates new <output.vcf> file from <input.vcf> containing only \
selected bead types. Also <start> timesteps can be omitted and every <skip> \
timestep can be left out.\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input.vcf> <input add> ", argv[0]);
      fprintf(stdout, "<output.vcf> <options>\n\n");

      fprintf(stdout, "   <input.vcf>       input filename (vcf format)\n");
      fprintf(stdout, "   <input add>       FIELD-like file with molecules to add\n");
      fprintf(stdout, "   <out.vcf>         output coordinate file (vcf format)\n");
      fprintf(stdout, "   <out.vsf>         output structure file (vsf format)\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -w             wrap coordinates (i.e., apply pbc)\n");
      fprintf(stdout, "      -st <start>    number of timestep to start from\n");
      fprintf(stdout, "      -xyz <name>    output xyz file\n");
      fprintf(stdout, "      -ld <float>    specify lowest distance from chosen bead types\n");
      fprintf(stdout, "      -hd <float>    specify highest distance from chosen bead types\n");
      fprintf(stdout, "      -bt <type(s)>  specify bead types new beads should be far from/near to (default: none)\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int req_args = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
//      strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-w") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-xyz") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-ld") != 0 &&
        strcmp(argv[i], "-hd") != 0) {
      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than traject.vsf? //{{{
  char *input_vsf = calloc(1024,sizeof(char *));
  if (FileOption(argc, argv, "-i", &input_vsf)) {
    exit(1);
  }
  if (input_vsf[0] == '\0') {
    strcpy(input_vsf, "traject.vsf");
  }

  // test if structure file ends with '.vsf'
  int ext = 2;
  char **extension;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_vsf, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(1024,sizeof(char *));
  if (FileOption(argc, argv, "-b", &bonds_file)) {
    exit(1);
  } //}}}

  // output verbosity //{{{
  bool verbose2, silent;
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  VerboseLongOption(argc, argv, &verbose, &verbose2); // more verbose output
  SilentOption(argc, argv, &verbose, &verbose2, &silent); // no output
  bool script = BoolOption(argc, argv, "--script"); // do not use \r & co.
  // }}}

  // should output coordinates be wrapped?
  bool wrap = BoolOption(argc, argv, "-w");

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // save into xyz file? //{{{
  char *output_xyz = calloc(1024,sizeof(char *));
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
      fprintf(stderr, "Error - if '-ld' and/or '-hd' is used, '-bt' must be specified as well\n\n");
      fprintf(stderr, "  Use: -bt <name(s)>\n");
      exit(1);
    }
  }
  //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input.vcf> - input coordinate file //{{{
  char input_coor[1024];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  ext = 2;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_coor, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <input add> - FIELD-like file with molecules to add //{{{
  char input_add[1024];
  strcpy(input_add, argv[++count]); //}}}

  // <out.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[1024];
  strcpy(output_vcf, argv[++count]);

  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  extension = malloc(ext*sizeof(char *));
  extension[0] = malloc(5*sizeof(char));
  strcpy(extension[0], ".vcf");
  if (!ErrorExtension(output_vcf, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <out.vsf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vsf[1024];
  strcpy(output_vsf, argv[++count]);

  // test if <out.vsf> filename ends with '.vsf' (required by VMD)
  ext = 1;
  extension = malloc(ext*sizeof(char *));
  extension[0] = malloc(5*sizeof(char));
  strcpy(extension[0], ".vsf");
  if (!ErrorExtension(output_vsf, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, false, Counts, &BeadType)) {
    exit(0);
  } //}}}

  // all beads & molecules are to be written //{{{
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = true;
//  BeadType[i].Use = true;
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = true;
//  MoleculeType[i].Use = true;
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_coor, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "\n   Starting from %d. timestep\n", start);
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[1024];
  // skip till 'pbc' keyword //{{{
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0); //}}}

  // read pbc //{{{
  Vector BoxLength;
  char line[1024];
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

  // print pbc if verbose output //{{{
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = calloc(1024,sizeof(int)); //}}}

  // skip first start-1 steps //{{{
  count = 0;
  int test;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent) {
      if (script) {
        fprintf(stdout, "Discarding step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rDiscarding step: %6d", count);
      }
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\nError: premature end of %s file\n\n", input_coor);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Discarded steps: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rDiscarded steps: %6d\n", count);
    }
  } //}}}

  // is the vcf file continuing? //{{{
  if ((test = getc(vcf)) == EOF) {
    fprintf(stderr, "\nError: %s - number of discard steps is lower (or equal) to the total number of steps\n\n", input_coor);
    exit(1);
  } else {
    ungetc(test,vcf);
  } //}}}
  //}}}

  // read the step to add stuff to //{{{
  if ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      fprintf(stdout, "Using step %6d\n", count);
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      fprintf(stdout, "\n%s", stuff);
  } else {
    fprintf(stderr, "\nWarning - using last step in %s (%d)\n", input_coor, count);
  }
  fclose(vcf); //}}}

  // create and zeroize new Counts structure //{{{
  struct Counts Counts_add;
  Counts_add.Beads = 0;
  Counts_add.Bonded = 0;
  Counts_add.Unbonded = 0;
  Counts_add.BeadsInVsf = 0;
  Counts_add.Molecules = 0;
  Counts_add.TypesOfBeads = 0;
  Counts_add.TypesOfMolecules = 0; //}}}

  // read stuff to be added //{{{
  // open the FIELD-like file for reading //{{{
  FILE *in_add;
  if ((in_add = fopen(input_add, "r")) == NULL) {
    ErrorFileOpen(input_add, 'r');
    exit(1);
  } //}}}

  // skip till 'species' keyword //{{{
  do {
    // get whole line - max 1000 chars
    fgets(line, 1024, in_add);

    // first string of the line
    split[0] = strtok(line, " \t");

  } while (strcmp(split[0], "species") != 0 &&
           strcmp(split[0], "SPECIES") != 0 &&
           strcmp(split[0], "Species") != 0); //}}}

  // after 'species' is number of bead types //{{{
  split[1] = strtok(NULL, " \t");
  Counts_add.TypesOfBeads = atoi(split[1]); //}}}

  // allocate structure for added bead types //{{{
  struct BeadType *BeadType_add;
  BeadType_add = malloc(Counts_add.TypesOfBeads*sizeof(struct BeadType)); //}}}

  // allocate structure for added beads (to be realloc'd later) //{{{
  struct Bead *Bead_add;
  Bead_add = malloc(1*sizeof(struct Bead)); //}}}

  // read bead type info //{{{
  for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
    fgets(line, 1024, in_add);
    // bead name //{{{
    split[0] = strtok(line, " \t");
    strcpy(BeadType_add[i].Name, split[0]); //}}}
    // bead mass //{{{
    split[1] = strtok(NULL, " \t");
    BeadType_add[i].Mass = atof(split[1]); //}}}
    // bead charge //{{{
    split[2] = strtok(NULL, " \t");
    BeadType_add[i].Charge = atof(split[2]); //}}}
    // number of unbonded beads //{{{
    split[3] = strtok(NULL, " \t");
    BeadType_add[i].Number = atoi(split[3]); //}}}

    // realloc & fill Bead_add //{{{
    int total_beads = Counts_add.Unbonded + BeadType_add[i].Number;
    if (total_beads > 0) {
      Bead_add = realloc(Bead_add, total_beads*sizeof(struct Bead));
      for (int j = Counts_add.Unbonded; j < total_beads; j++) {
        Bead_add[j].Type = i;
        Bead_add[j].Molecule = -1;
        Bead_add[j].nAggregates = 0;
        Bead_add[j].Aggregate = calloc(1,sizeof(int));
        Bead_add[j].Index = j; // probably useless here
      } //}}}
    }

    Counts_add.Unbonded += BeadType_add[i].Number;
    Counts_add.Beads += BeadType_add[i].Number;
  } //}}}

  // skip till 'molecule' keyword //{{{
  do {
    // get whole line - max 1000 chars
    fgets(line, 1024, in_add);

    // first string of the line
    split[0] = strtok(line, " \t");

  } while (strcmp(split[0], "molecules") != 0 &&
           strcmp(split[0], "MOLECULES") != 0 &&
           strcmp(split[0], "Molecules") != 0); //}}}

  // after 'molecule' is number of molecule types //{{{
  split[1] = strtok(NULL, " \t");
  Counts_add.TypesOfMolecules = atoi(split[1]); //}}}

  // allocate structure for added molecule types //{{{
  struct MoleculeType *MoleculeType_add;
  MoleculeType_add = malloc(Counts_add.TypesOfMolecules*sizeof(struct MoleculeType));
  struct Vector *prototype[Counts_add.TypesOfMolecules]; //}}}

  // allocate structure for added molecules (to be realloc'd later) //{{{
  struct Molecule *Molecule_add;
  Molecule_add = malloc(1*sizeof(struct Molecule)); //}}}

  // read molecule info //{{{
  for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
    // molecule name //{{{
    fgets(line, 1024, in_add);
    // trim trailing whitespace in line //{{{
    int length = strlen(line);
    // last string character needs to be '\0'
    while (length > 1 &&
           (line[length-1] == ' ' ||
            line[length-1] == '\n' ||
            line[length-1] == '\t')) {
      line[length-1] = '\0';
      length--;
    } //}}}
    split[0] = strtok(line, " \t");
    strcpy(MoleculeType_add[i].Name, split[0]); //}}}
    // number of molecules of given type //{{{
    fgets(line, 1024, in_add);
    split[0] = strtok(line, " \t");
    split[1] = strtok(NULL, " \t");
    MoleculeType_add[i].Number = atoi(split[1]); //}}}
    // number of beads in molecules of given type //{{{
    fgets(line, 1024, in_add);
    split[0] = strtok(line, " \t");
    split[1] = strtok(NULL, " \t");
    MoleculeType_add[i].nBeads = atoi(split[1]); //}}}

    // total number of beads in the molecules of type 'i'
    int beads = MoleculeType_add[i].Number*MoleculeType_add[i].nBeads;

    // realloc _add structures //{{{
    printf("xxx %d %d\n", Counts_add.Beads, beads);
    Bead_add = realloc(Bead_add, (Counts_add.Beads+beads)*sizeof(struct Bead));
    printf("xxx %d %d\n", Counts_add.Molecules, MoleculeType_add[i].Number);
    Molecule_add = realloc(Molecule_add, (Counts_add.Molecules+MoleculeType_add[i].Number)*sizeof(struct Molecule));
    for (int j = Counts_add.Molecules; j < (Counts_add.Molecules+MoleculeType_add[i].Number); j++) {
      Molecule_add[j].Bead = malloc(MoleculeType_add[i].nBeads*sizeof(int));
    } //}}}

    // allocate array for coordinates of prototype molecule of type 'i'
    prototype[i] = malloc(MoleculeType_add[i].nBeads*sizeof(struct Vector));

    MoleculeType_add[i].BType = malloc(1*sizeof(int));
    MoleculeType_add[i].nBTypes = 0;
    MoleculeType_add[i].Mass = 0;

    // read bead types and coordinates
    Vector zero_first;
    zero_first.x = 0;
    zero_first.y = 0;
    zero_first.z = 0;
    for (int j = 0; j < MoleculeType_add[i].nBeads; j++) {
      fgets(line, 1024, in_add);
      // bead name //{{{
      split[0] = strtok(line, " \t");
      // is the bead type registered in the molecule already?
      int type = FindBeadType(split[0], Counts_add, BeadType_add);
      bool exists = false;
      for (int k = 0; k < MoleculeType_add[i].nBTypes; k++) {
        if (type == MoleculeType_add[i].BType[k]) {
          exists = true;
          break;
        }
      }
      // add new bead type to molecule if not yet registered
      if (!exists) {
        MoleculeType_add[i].nBTypes++;
        MoleculeType_add[i].BType = realloc(MoleculeType_add[i].BType, MoleculeType_add[i].nBTypes*sizeof(int));
        MoleculeType_add[i].BType[MoleculeType_add[i].nBTypes-1] = type;
      } //}}}

      // add bead's mass to molecule mass
      MoleculeType_add[i].Mass += BeadType_add[type].Mass;

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

      // fill _add structures //{{{
      for (int l = 0; l < MoleculeType_add[i].Number; l++) {
        int index = Counts_add.Beads + l * MoleculeType_add[i].nBeads + j;

        Bead_add[index].Type = type;
        Bead_add[index].Molecule = Counts_add.Molecules + l;
        Bead_add[index].Index = index;
        Bead_add[index].nAggregates = 0; // useless here
        Bead_add[index].Aggregate = malloc(1*sizeof(int)); // just to free later
        Molecule_add[Counts_add.Molecules+l].Bead[j] = index;
        Molecule_add[Counts_add.Molecules+l].Type = i;
        BeadType_add[type].Number++;
      } //}}}
    }

    Counts_add.Bonded += beads;
    Counts_add.Beads += beads;
    Counts_add.BeadsInVsf += beads;
    Counts_add.Molecules += MoleculeType_add[i].Number;

    // number of bonds in molecules of given type //{{{
    fgets(line, 1024, in_add);
    split[0] = strtok(line, " \t");
    split[1] = strtok(NULL, " \t");
    MoleculeType_add[i].nBonds = atoi(split[1]); //}}}

    // allocate Bond array //{{{
    MoleculeType_add[i].Bond = malloc(MoleculeType_add[i].nBonds*sizeof(int *));
    for (int j = 0; j < MoleculeType_add[i].nBonds; j++) {
      MoleculeType_add[i].Bond[j] = calloc(2,sizeof(int));
    } //}}}

    // read bond info //{{{
    for (int j = 0; j < MoleculeType_add[i].nBonds; j++) {
      fgets(line, 1024, in_add);
      split[0] = strtok(line, " \t");
      split[1] = strtok(NULL, " \t");
      split[2] = strtok(NULL, " \t");
      MoleculeType_add[i].Bond[j][0] = atoi(split[1]) - 1; // bead ids in FIELD start from 1
      MoleculeType_add[i].Bond[j][1] = atoi(split[2]) - 1;
    } //}}}

    // skip till 'finish' keyword //{{{
    do {
      // get whole line - max 1000 chars
      fgets(line, 1024, in_add);
      // trim trailing whitespace in line //{{{
      int length = strlen(line);
      // last string character needs to be '\0'
      while (length > 1 &&
             (line[length-1] == ' ' ||
              line[length-1] == '\n' ||
              line[length-1] == '\t')) {
        line[length-1] = '\0';
        length--;
      } //}}}
      // first string of the line
      split[0] = strtok(line, " \t");

    } while (strcmp(split[0], "finish") != 0 &&
             strcmp(split[0], "Finish") != 0 &&
             strcmp(split[0], "FINISH") != 0); //}}}
  } //}}}

  Counts_add.BeadsInVsf = Counts_add.Beads;

  fclose(in_add); //}}}

  // print what is to be added //{{{
  fprintf(stdout, "To add:\n");

  fprintf(stdout, "   Counts.{");
  fprintf(stdout, "TypesOfBeads =%3d, ", Counts_add.TypesOfBeads);
  fprintf(stdout, "Bonded =%7d, ", Counts_add.Bonded);
  fprintf(stdout, "Unbonded =%7d, ", Counts_add.Unbonded);
  fprintf(stdout, "Beads = %7d, ", Counts_add.Beads);
  fprintf(stdout, "TypesOfMolecules =%3d, ", Counts_add.TypesOfMolecules);
  fprintf(stdout, "Molecules =%4d}\n", Counts_add.Molecules);
  putchar('\n');

  for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
    fprintf(stdout, "   BeadType[%2d].{", i);
    fprintf(stdout, "Name =%10s, ", BeadType_add[i].Name);
    fprintf(stdout, "Number =%7d, ", BeadType_add[i].Number);
    fprintf(stdout, "Charge =%6.2f, ", BeadType_add[i].Charge);
    fprintf(stdout, "Mass =%5.2f}\n", BeadType_add[i].Mass);
  }
  putchar('\n');

  for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
    fprintf(stdout, "   MoleculeType[%d].{", i);
    fprintf(stdout, "Name =%10s", MoleculeType_add[i].Name);
    fprintf(stdout, ", Number =%4d", MoleculeType_add[i].Number);
    fprintf(stdout, ", nBeads =%3d", MoleculeType_add[i].nBeads);
//  fprintf(stdout, ", nBonds =%3d", MoleculeType_add[i].nBonds);

    fprintf(stdout, ", nBTypes =%2d, BType{", MoleculeType_add[i].nBTypes);
    fprintf(stdout, "%8s", BeadType_add[MoleculeType_add[i].BType[0]].Name);
    for (int j = 1; j < MoleculeType_add[i].nBTypes; j++) {
      fprintf(stdout, ",%8s", BeadType_add[MoleculeType_add[i].BType[j]].Name);
    }
    putchar('}');
    fprintf(stdout, ", Mass =%7.2f", MoleculeType_add[i].Mass);
    fprintf(stdout, "}\n");
  }
  putchar('\n');

  fprintf(stdout, "   Molecule prototypes \n");
  for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
    fprintf(stdout, "     %10s: %10.5f %10.5f %10.5f\n", MoleculeType_add[i].Name, prototype[i][0].x, prototype[i][0].y, prototype[i][0].z);
    for (int j = 1; j < MoleculeType_add[i].nBeads; j++) {
      fprintf(stdout, "                 %10.5f %10.5f %10.5f\n", prototype[i][j].x, prototype[i][j].y, prototype[i][j].z);
    }
    putchar('\n');
  }

  if (lowest_dist != -1) {
    fprintf(stdout, "Lowest distance from specied beads: %lf\n", lowest_dist);
  }
  if (highest_dist != -1) {
    fprintf(stdout, "Highest distance from specied beads: %lf\n", highest_dist);
  } putchar('\n'); //}}}

  // square lowest/highest distance, if '-ld' and/or '-hd' options used //{{{
  if (lowest_dist != -1) {
    lowest_dist = SQR(lowest_dist);
  }
  if (highest_dist != -1) {
    highest_dist = SQR(highest_dist);
  } //}}}

  // add monomeric beads //{{{
  count = 0;
  for (int i = 0; i < Counts_add.Beads; i++) {
    if (Bead_add[i].Molecule == -1) {
      Vector random;
      random.x = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.x; // random number <0,BoxLength)
      random.y = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.y;
      random.z = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.z;
      for (int j = count; j < Counts.Beads; j++) {
        if (Bead[j].Molecule == -1 && BeadType[Bead[j].Type].Charge == 0) {
          Bead[j].Position.x = random.x;
          Bead[j].Position.y = random.y;
          Bead[j].Position.z = random.z;
          Bead[j].Type = Counts.TypesOfBeads + Bead_add[i].Type;
          Bead_add[i].Index = Bead[j].Index;
          break;
        }
      }
    }
  } //}}}

  // add the molecules //{{{
  for (int i = 0; i < Counts_add.Molecules; i++) {
    int mol_type = Molecule_add[i].Type;

    Vector random;
//  // rotate the prototype molecule randomly //{{{
//  // random rotation axis
//  random.x = (double)rand() / ((double)RAND_MAX); // random number <0,1>
//  random.y = (double)rand() / ((double)RAND_MAX);
//  random.z = (double)rand() / ((double)RAND_MAX);
//  // random rotation angle
//  double angle = (double)rand() / ((double)RAND_MAX) * PI;
//  // create rotation matrix
//  struct Tensor {
//    Vector x, y, z;
//  } rot;
//  rot.x.x = cos(angle) + SQR(random.x) * (1 - cos(angle));
//  rot.x.y = random.x * random.y * (1 - cos(angle)) - random.z * sin(angle);
//  rot.x.z = random.x * random.z * (1 - cos(angle)) + random.y * sin(angle);
//  rot.y.x = random.x * random.y * (1 - cos(angle)) + random.z * sin(angle);
//  rot.y.y = cos(angle) + SQR(random.y) * (1 - cos(angle));
//  rot.y.z = random.y * random.z * (1 - cos(angle)) - random.x * sin(angle);
//  rot.z.x = random.x * random.z * (1 - cos(angle)) - random.y * sin(angle);
//  rot.z.y = random.y * random.z * (1 - cos(angle)) + random.x * sin(angle);
//  rot.z.z = cos(angle) + SQR(random.z) * (1 - cos(angle));
//  // transform the prototype molecule (rotation matrix x coordinates)
//  Vector rotated[MoleculeType_add[mol_type].nBeads];
//  for (int j = 0; j < MoleculeType_add[mol_type].nBeads; j++) {
//    rotated[j].x = rot.x.x * prototype[mol_type][j].x + rot.x.y * prototype[mol_type][j].y + rot.x.z * prototype[mol_type][j].z;
//    rotated[j].y = rot.y.x * prototype[mol_type][j].x + rot.y.y * prototype[mol_type][j].y + rot.y.z * prototype[mol_type][j].z;
//    rotated[j].z = rot.z.x * prototype[mol_type][j].x + rot.z.y * prototype[mol_type][j].y + rot.z.z * prototype[mol_type][j].z;
//  } //}}}

    // for now: only first bead's distance from stuff is checked //{{{
    double min_dist;
    if (lowest_dist != -1 || highest_dist != -1) {
      do {
        random.x = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.x; // random number <0,BoxLength)
        random.y = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.y;
        random.z = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.z;

        min_dist = SQR(BoxLength.x * 100);
        for (int j = 0; j < Counts.Beads; j++) {
          int btype = Bead[j].Type;
          if (BeadType[btype].Use) {
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
      random.x = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.x; // random number <0,BoxLength)
      random.y = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.y;
      random.z = (double)rand() / ((double)RAND_MAX + 1) * BoxLength.z;
    } //}}}

    fflush(stdout);
    fprintf(stdout, "\rMolecules placed: %d", i);
    for (int j = 0; j < MoleculeType_add[mol_type].nBeads; j++) {
      int id = Molecule_add[i].Bead[j];
      for (int k = count; k < Counts.Beads; k++) {
        count++; // so that the for l doesn't go through all beads every time
        if (Bead[k].Molecule == -1 && BeadType[Bead[k].Type].Charge == 0) {
          Vector r;
//        r.x = random.x + rotated[j].x;
//        r.y = random.y + rotated[j].y;
//        r.z = random.z + rotated[j].z;
          r.x = random.x + prototype[mol_type][j].x;
          r.y = random.y + prototype[mol_type][j].y;
          r.z = random.z + prototype[mol_type][j].z;

          Bead[k].Position.x = r.x;
          Bead[k].Position.y = r.y;
          Bead[k].Position.z = r.z;
          Bead[k].Molecule = Counts.Molecules + Bead_add[id].Molecule;
          Bead[k].Type = Counts.TypesOfBeads + Bead_add[id].Type;
          Bead_add[id].Index = Bead[k].Index;
          break;
        }
      }
    }
  } //}}}

// test print //{{{
//for (int i = 0; i < Counts_add.Molecules; i++) {
//  printf("%s (%d):", MoleculeType_add[Molecule_add[i].Type].Name, MoleculeType_add[Molecule_add[i].Type].nBonds);
//  for (int j = 0; j < MoleculeType_add[Molecule_add[i].Type].nBonds; j++) {
//    printf(" %d-%d", Bead[Molecule_add[i].Bead[MoleculeType_add[Molecule_add[i].Type].Bond[j][0]]].Index,
//                     Bead[Molecule_add[i].Bead[MoleculeType_add[Molecule_add[i].Type].Bond[j][1]]].Index);
//  }
//  putchar('\n');
//  for (int j = 0; j < MoleculeType_add[Molecule_add[i].Type].nBeads; j++) {
//    printf(" %d", Molecule_add[i].Bead[j]);
//  }
//  putchar('\n');
//} //}}}

  // join _add structs to original structs //{{{
  Counts.TypesOfBeads += Counts_add.TypesOfBeads;
  Counts.TypesOfMolecules += Counts_add.TypesOfMolecules;
  Counts.Molecules += Counts_add.Molecules;
  // realloc structs
  BeadType = realloc(BeadType, Counts.TypesOfBeads*sizeof(struct BeadType));
  MoleculeType = realloc(MoleculeType, Counts.TypesOfMolecules*sizeof(struct MoleculeType));
  Bead = realloc(Bead, Counts.Beads*sizeof(struct Bead));
  Molecule = realloc(Molecule, Counts.Molecules*sizeof(struct Molecule));
  // bead types //{{{
  for (int i = 0; i < Counts_add.TypesOfBeads; i++) {
    int new = Counts.TypesOfBeads - Counts_add.TypesOfBeads + i;
    strcpy(BeadType[new].Name, BeadType_add[i].Name);
    BeadType[new].Number = BeadType_add[i].Number;
    BeadType[new].Charge = BeadType_add[i].Charge;
    BeadType[new].Mass = BeadType_add[i].Mass;
    BeadType[new].Write = true;
//  BeadType[new].Use = true;
  } //}}}
  // molecule types //{{{
  for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
    int new = Counts.TypesOfMolecules - Counts_add.TypesOfMolecules + i;
    strcpy(MoleculeType[new].Name, MoleculeType_add[i].Name);
    MoleculeType[new].Number = MoleculeType_add[i].Number;
    MoleculeType[new].nBeads = MoleculeType_add[i].nBeads;
    MoleculeType[new].nBonds = MoleculeType_add[i].nBonds;
    MoleculeType[new].Bond = malloc(MoleculeType[new].nBonds*sizeof(int *));
    for (int j = 0; j < MoleculeType[new].nBonds; j++) {
      MoleculeType[new].Bond[j] = calloc(2,sizeof(int));
      MoleculeType[new].Bond[j][0] = MoleculeType_add[i].Bond[j][0];
      MoleculeType[new].Bond[j][1] = MoleculeType_add[i].Bond[j][1];
    }
    MoleculeType[new].nBTypes = MoleculeType_add[i].nBTypes;
    MoleculeType[new].BType = malloc(MoleculeType[new].nBTypes*sizeof(int));
    for (int j = 0; j < MoleculeType[new].nBTypes; j++) {
      MoleculeType[new].BType[j] = Counts.TypesOfBeads - Counts_add.TypesOfBeads + MoleculeType_add[i].BType[j];
    }
    MoleculeType[new].Mass = MoleculeType_add[i].Mass;
    MoleculeType[new].InVcf = true; // some of those flags are useless (and not just here)
    MoleculeType[new].Write = true;
//  MoleculeType[new].Use = true;
  } //}}}
  // molecules //{{{
  for (int i = 0; i < Counts_add.Molecules; i++) {
    int new = Counts.Molecules - Counts_add.Molecules + i;
    Molecule[new].Type = Counts.TypesOfMolecules - Counts_add.TypesOfMolecules + Molecule_add[i].Type;
    int mol_type = Molecule[new].Type;
    Molecule[new].Bead = malloc(MoleculeType[mol_type].nBeads*sizeof(int));
    for (int j = 0; j < MoleculeType[mol_type].nBeads; j++) {
      Molecule[new].Bead[j] = Bead_add[Molecule_add[i].Bead[j]].Index;
    }
    Molecule[new].Aggregate = 0; // probably useless here
  } //}}}
  //}}}

  // recount bonded & unbonded beads - not necessary, just for accounting //{{{
  Counts.Unbonded = 0;
  Counts.Bonded = 0;
  for (int i = 0; i < Counts.Beads; i++) {
    if (Bead[i].Molecule == -1) {
      Counts.Unbonded++;
    } else {
      Counts.Bonded++;
    }
  } //}}}

  // print overall system //{{{
  fprintf(stdout, "\nOld + new (if added molecules/beads are of already known type, they appear twice):\n");
//VerboseOutput(verbose2, input_coor, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // open xyz file and write coordinates (if -xyz option is present) //{{{
  FILE *xyz;
  if (output_xyz[0] != '\0') {
    // open output .xyz file for reading //{{{
    if ((xyz = fopen(output_xyz, "w")) == NULL) {
      ErrorFileOpen(output_xyz, 'w');
      exit(1);
    } //}}}

    WriteCoorXYZ(xyz, Counts, BeadType, Bead);
  } //}}}

  // open output .vcf file //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  } //}}}

  // print command, bead type names & box size to output .vcf file //{{{
  // command
  putc('#', out);
  for (int i = 0; i < argc; i++) {
    fprintf(out, " %s", argv[i]);
  }
  fprintf(out, "\n\n");

  // bead type names
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(out, "# %s\n", BeadType[i].Name);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z); //}}}

  // print coordinates to output .vcf file
  WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

  // write wrapped coordinates as well? //{{{
  if (wrap) {
    RestorePBC(Counts, BoxLength, &Bead);
    WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

    // write coordinates to xyz (if -xyz option is present)
    if (output_xyz[0] != '\0') {
      WriteCoorXYZ(xyz, Counts, BeadType, Bead);
    }
  } //}}}

  // close output files //{{{
  fclose(out);
  if (output_xyz[0] != '\0') {
    fclose(xyz);
  } //}}}

  // create & fill output vsf file
  WriteVsf(output_vsf, Counts, BeadType, Bead, MoleculeType, Molecule);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(BeadType_add);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMoleculeType(Counts_add, &MoleculeType_add);
  FreeMolecule(Counts, &Molecule);
  FreeMolecule(Counts_add, &Molecule_add);
  FreeBead(Counts, &Bead);
  FreeBead(Counts_add, &Bead_add);
  free(stuff);
  free(output_xyz);
  for (int i = 0; i < Counts_add.TypesOfMolecules; i++) {
    free(prototype[i]);
  }
  //}}}

  return 0;
}
