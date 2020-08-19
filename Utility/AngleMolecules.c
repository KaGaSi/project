#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
AngleMolecules utility calculates angles for specified beads (three beads per \
angle) for molecules of specified type(s). The specified beads do not have \
to be connected by bonds.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol name(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <width>           width of a single bin in degrees\n");
  fprintf(ptr, "   <output>          output file with distribution of angles\n");
  fprintf(ptr, "   <mol name(s)>     molecule name(s) to calculate angles for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined coordinates\n");
  fprintf(ptr, "      -n <ints>      bead indices (multiple of 3 <ints>) for angle calculation (default: 1 2 3)\n");
  fprintf(ptr, "      -a <name>      write angles of all molecules in all times to <name>\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
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
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-a") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - filename of input vcf file (must end with .vcf) //{{{
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

  // <width> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (!IsPosDouble(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]);

  // number of bins between 0 and 180 deg
  int bins = ceil(180 / width); //}}}

  // <output> - file name with dihedral angle distribution //{{{
  char output_distr[LINE];
  strcpy(output_distr, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

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
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n", start, end);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
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

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {

    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (mol_type == -1) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m non-existent molecule \033[1;33m%s\033[1;31m\n", input_coor, argv[count]);
      fprintf(stderr, "\033[0m");
      ErrorMoleculeType(Counts, MoleculeType);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }
  } //}}}

  // '-n' option - specify bead ids //{{{
  int bead[100] = {0}, number_of_beads = 3, beads_per_angle = 3, test = 0;
  bead[0] = 1; // default ids for angle
  bead[1] = 2;
  bead[2] = 3;
  if (MultiIntegerOption(argc, argv, "-n", &test, bead)) {
    exit(1);
  }
  if (test != 0) { // -n is present
    number_of_beads = test;
  }

  // Error: wrong number of integers //{{{
  if ((number_of_beads%beads_per_angle) != 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: \033[1;33m-n\033[1;31m option - number of bead ids must be dividable by three\n\n");
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero

    // Error - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && bead[i] >= MoleculeType[j].nBeads) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m-n\033[1;31m option - \033[1;33m%d\033[1;31m is larger", bead[i]);
        fprintf(stderr, " than the number of beads in molecule \033[1;33m%s\033[1;31m\n\n", MoleculeType[j].Name);
        fprintf(stderr, "\033[0m");
        Help(argv[0], true);
        exit(1);
      }
    } //}}}
  } //}}}

  // '-a' option - write angles for all molecules //{{{
  char *output = calloc(LINE,sizeof(char));
  if (FileOption(argc, argv, "-a", &output)) {
    exit(1);
  }

  // write initial stuff to output if '-a' is used
  if (output[0] != '\0') {
    // open output //{{{
    FILE *out;
    if ((out = fopen(output, "w")) == NULL) {
      ErrorFileOpen(output, 'w');
      exit(1);
    } //}}}

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++)
      fprintf(out, " %s", argv[i]);
    putc('\n', out); //}}}

    // print molecule names & bead ids //{{{
    fprintf(out, "# angles between beads:");
    for (int j = 0; j < number_of_beads; j += beads_per_angle) {
      fprintf(out, " (%d) %d-%d-%d;", j/beads_per_angle+1, bead[j]+1, bead[j+1]+1, bead[j+2]+1);
    }
    putc('\n', out);
    fprintf(out, "# columns: (1) step;");
    int j = 2;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        if ((number_of_beads/beads_per_angle*MoleculeType[i].Number) == 1) {
          fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
        } else {
          fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_beads/beads_per_angle*MoleculeType[i].Number-1, MoleculeType[i].Name);
        }
        j += number_of_beads / beads_per_angle * MoleculeType[i].Number;
      }
    }
    putc('\n', out); //}}}

    fclose(out);
  } //}}}

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
    fprintf(stdout, "Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(stdout, " %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
  } //}}}

  // allocate array for distribution of angles //{{{
  double avg_angle[Counts.TypesOfMolecules][number_of_beads/beads_per_angle];
  double *distr[Counts.TypesOfMolecules][number_of_beads/beads_per_angle];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
      avg_angle[i][j] = 0;
      distr[i][j] = calloc(bins, sizeof(double));
    }
  } //}}}

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
  // print starting step? //{{{
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
  count = 0; // count calculated timesteps
  int count_vcf = start - 1;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count_vcf, stuff);
      exit(1);
    } //}}}

    // join molecules if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // calculate angles //{{{
    double angle[Counts.Molecules][number_of_beads/beads_per_angle];
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type_i = Molecule[i].Type;
      if (MoleculeType[mol_type_i].Use) {
        for (int j = 0; j < number_of_beads; j += beads_per_angle) {
          VECTOR u, v;
          // first vector
          u.x = Bead[Molecule[i].Bead[bead[j+0]]].Position.x - Bead[Molecule[i].Bead[bead[j+1]]].Position.x;
          u.y = Bead[Molecule[i].Bead[bead[j+0]]].Position.y - Bead[Molecule[i].Bead[bead[j+1]]].Position.y;
          u.z = Bead[Molecule[i].Bead[bead[j+0]]].Position.z - Bead[Molecule[i].Bead[bead[j+1]]].Position.z;
          // second vector
          v.x = Bead[Molecule[i].Bead[bead[j+2]]].Position.x - Bead[Molecule[i].Bead[bead[j+1]]].Position.x;
          v.y = Bead[Molecule[i].Bead[bead[j+2]]].Position.y - Bead[Molecule[i].Bead[bead[j+1]]].Position.y;
          v.z = Bead[Molecule[i].Bead[bead[j+2]]].Position.z - Bead[Molecule[i].Bead[bead[j+1]]].Position.z;
          // calculate angle between the two vectors
          double size[2];
          size[0] = Length(u);
          size[1] = Length(v);
          double scalar = u.x * v.x + u.y * v.y + u.z * v.z;
          angle[i][j/beads_per_angle] = acos(scalar / (size[0] * size[1])); // in rad
          angle[i][j/beads_per_angle] *= 180 / PI; // in degrees

          // add to average
          avg_angle[mol_type_i][j/beads_per_angle] += angle[i][j/beads_per_angle];

          int k = angle[i][j/beads_per_angle] / width;
          if (k < bins) {
            distr[mol_type_i][j/beads_per_angle][k]++;
          } else {
            fprintf(stderr, "\033[1;33m");
            fprintf(stdout, "\nWarning: weird angle: %lf degrees\n", angle[i][j/beads_per_angle]);
            fprintf(stderr, "\033[0m");
          }
        }
      }
    } //}}}

    // write all angles to output if '-a' is used //{{{
    if (output[0] != '\0') {
      FILE *out;
      if ((out = fopen(output, "a")) == NULL) {
        ErrorFileOpen(output, 'a');
        exit(1);
      }

      fprintf(out, "%6d", count_vcf);
      for (int i = 0; i < Counts.Molecules; i++) {
        int mol_type_i = Molecule[i].Type;
        if (MoleculeType[mol_type_i].Use) {
          for (int j = 0; j < number_of_beads; j += beads_per_angle){
            fprintf(out, " %10.6f", angle[i][j/beads_per_angle]);
          }
        }
      }
      putc('\n', out);

      fclose(out);
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

  // write distribution of angles //{{{
  // open output file for appending //{{{
  FILE *out;
  if ((out = fopen(output_distr, "w")) == NULL) {
    ErrorFileOpen(output_distr, 'w');
    exit(1);
  } //}}}

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print molecule names & bead ids //{{{
  fprintf(out, "# angles between beads:");
  for (int j = 0; j < number_of_beads; j += beads_per_angle) {
    fprintf(out, " (%d) %d-%d-%d;", j/beads_per_angle+1, bead[j]+1, bead[j+1]+1, bead[j+2]+1);
  }
  putc('\n', out);
  fprintf(out, "# columns: (1) angle [deg];");
  count = 2;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if ((number_of_beads/beads_per_angle) == 1) {
        fprintf(out, " (%d) %s molecules;", count, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", count, count+number_of_beads/beads_per_angle-1, MoleculeType[i].Name);
      }
      count += number_of_beads / beads_per_angle;
    }
  }
  putc('\n', out); //}}}

  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%5.1f", width*(2*i+1)/2);

    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {

        for (int k = 0; k < (number_of_beads/beads_per_angle); k++) {
          fprintf(out, "%10f", (double)(distr[j][k][i])/(count*MoleculeType[j].Number));
        }
      }
    }
    putc('\n', out);
  } //}}}

  // write to output average angles //{{{
  fprintf(out, "# simple averages:");
  count = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if ((number_of_beads/beads_per_angle) == 1) {
        fprintf(out, " (%d) %s molecules;", count, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", count, count+number_of_beads/beads_per_angle-1, MoleculeType[i].Name);
      }
      count += number_of_beads / beads_per_angle;
    }
  }
  fprintf(out, "\n#");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
        fprintf(out, " %7.3f", avg_angle[i][j]/(count*MoleculeType[i].Number));
      }
    }
  }
  putc('\n', out); //}}}

  fclose(out);
  //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(output);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
      free(distr[i][j]);
    }
  } //}}}

  return 0;
}
