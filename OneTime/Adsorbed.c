#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
A WORK IN PROGRESS\n\
Adsorbed calculates the number of specified beads near a wall per number of \
molecules these beads are in.\n\n");
  }

  fprintf(ptr, "Usage:\n");
//fprintf(ptr, "   %s <input> <distance> <output> <axis> <bead name(s)> <options>\n\n", cmd);
  fprintf(ptr, "   %s <input> <distance> <axis> <bead name(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename (either vcf or vtf format)\n");
  fprintf(ptr, "   <distance>        distance from the box edges\n");
//fprintf(ptr, "   <output>          output density file (automatic ending '<axis>.rho' added)\n");
  fprintf(ptr, "   <axis>            calculate along x, y, or z axis\n");
  fprintf(ptr, "   <bead name(s)>    names of bead types for closeness calculation\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       number of timestep to end with\n");
//fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
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
//int req_args = 4;
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
//      strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-st") != 0) {
//      strcmp(argv[i], "-x") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
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
  char *input_vsf = calloc(LINE,sizeof(char));
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
  } //}}}

  // <distance> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<distance>");
    Help(argv[0], true);
    exit(1);
  }
  double distance = atof(argv[count]); //}}}

//// <output.rho> - filename with bead densities //{{{
//char output_rho[LINE];
//strcpy(output_rho, argv[++count]); //}}}

  // <axis> - x, y, or z //{{{
  char axis = argv[++count][0];

  if (axis != 'x' && axis != 'y' && axis != 'z') {
    fprintf(stderr, "\nError: <axis> must be 'x', 'y', or 'z'\n\n");
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // ending timestep //{{{
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  } //}}}

  // error if ending step is lower than starging step //{{{
  if (end != -1 && start > end) {
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n", start, end);
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

//// '-x' option //{{{
//if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
//  exit(1);
//} //}}}

  // <type names> - names of bead types to use for closeness calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\nError: %s - non-existent bead name '%s'\n", input_coor, argv[count]);
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    if (BeadType[type].Use) {
      fprintf(stderr, "\nError: bead type %s specified more than once\n\n", argv[count]);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}

//// write initial stuff to output file //{{{
//FILE *out;
//char str[1050];

//sprintf(str, "%s%c.rho", output_rho, axis);
//strcpy(output_rho, str);
//if ((out = fopen(output_rho, "w")) == NULL) {
//  ErrorFileOpen(output_rho, 'w');
//  exit(1);
//}

//// print command to output file //{{{
//putc('#', out);
//for (int i = 0; i < argc; i++)
//  fprintf(out, " %s", argv[i]);
//putc('\n', out); //}}}

//// print bead type names to output file //{{{
//fprintf(out, "# columns: (1) distance;");
//for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  fprintf(out, " (%d) %s", i+2, BeadType[i].Name);
//  if (i != (Counts.TypesOfBeads-1)) {
//    putc(';', out);
//  }
//}
//putc('\n', out);
//for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  fprintf(out, " %d: %s", 4*i+2, BeadType[i].Name);
//}
//fprintf(out, "\n# for each molecule type: rdp | stderr | rnp | stderr\n"); //}}}

//fclose(out); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "Chosen bead types:");
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      if (BeadType[i].Use) {
        fprintf(stdout, " %s", BeadType[i].Name);
      }
    }
    putchar('\n');
  } //}}}

  // skip first start-1 steps //{{{
  count = 0;
  int test;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "Error: premature end of %s file\n\n", input_coor);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step: %d\n", start);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rStarting step: %d\n", start);
    }
  } //}}}
  // is the vcf file continuing?
  if (ErrorDiscard(start, count, input_coor, vcf)) {
    exit(1);
  }
  //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int count_vcf = start - 1;
  long int count_near = 0; // number of specified beads near box edges
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // write step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count_vcf, stuff, input_vsf);
      exit(1);
    } //}}}

    RestorePBC(Counts, BoxLength, &Bead);

    // calculate beads near box edges //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      int btype = Bead[i].Type;
      if (Bead[i].Molecule != -1 && BeadType[btype].Use) {
        if (axis == 'x') {
          if (Bead[i].Position.x <= distance || Bead[i].Position.x >= (BoxLength.x-distance)) {
            count_near++;
          }
        } else if (axis == 'y') {
          if (Bead[i].Position.y <= distance || Bead[i].Position.y >= (BoxLength.y-distance)) {
            count_near++;
          }
        } else {
          if (Bead[i].Position.z <= distance || Bead[i].Position.z >= (BoxLength.z-distance)) {
            count_near++;
          }
        }
      }
    } //}}}

    if (end == count_vcf)
      break;
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rLast Step: %d\n", count_vcf);
    }
  } //}}}

  // count molecules containing used bead types //{{{
  int molecules = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    bool use = false;
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      int btype = Bead[MoleculeType[i].Bead[j]].Type;
      if (BeadType[btype].Use) {
        use = true;
        break;
      }
    }
    if (use) {
      molecules += MoleculeType[i].Number;
    }
  } //}}}

  // number of beads of specified type near box edges per molecule containing those beads
  printf("%lf\n", count_near/(double)(molecules)/count);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff); //}}}

  return 0;
}
