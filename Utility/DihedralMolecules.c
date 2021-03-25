#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
DihedralMolecules utility calculates angle between planes specified by beads \
in each molecule of specified molecule type(s). Each plane can be specified \
by an arbitrary trio of beads in the molecules.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol name(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>         input filename (either vcf or vtf format)\n");
  fprintf(ptr, "   <width>         width of a single bin in degrees\n");
  fprintf(ptr, "   <output>        output file with distribution of dihedral angles\n");
  fprintf(ptr, "   <mol name(s)>   molecule names to calculate angles for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined     specify that <input> contains joined coordinates\n");
  fprintf(ptr, "      -n <ints>    bead indices (multiple of 6 <ints>) for dihedral calculation (default: 1 2 3 2 3 4)\n");
  fprintf(ptr, "      -a <name>    write angle of all molecules in all times to <name>\n");
  fprintf(ptr, "      -st <int>    starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>     ending timestep for the calculation\n");
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

  int req_args = 4; //}}}

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
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-a") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
  char *input_vsf = calloc(LINE,sizeof(char));
  strcpy(input_coor, argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
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

  // starting & ending timesteps //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  }
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  }
  ErrorStartEnd(start, end); //}}}

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  VECTOR BoxLength; // couboid box dimensions
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &BoxLength, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule);
  free(input_vsf); //}}}

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);
    if (mol_type == -1) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - non-existent molecule");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", argv[count]);
      RedText(STDERR_FILENO);
      fprintf(stderr, "\n");
      ResetColour(STDERR_FILENO);
      ErrorMoleculeType(Counts, MoleculeType);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }
  } //}}}

  // '-n' option - specify bead ids //{{{
  int dihedral[100] = {0}, number_of_beads = 6, beads_per_angle = 6, test = 0;
  dihedral[0] = 1; // default planes: 1-2-3 & 2 3 4
  dihedral[1] = 2;
  dihedral[2] = 3;
  dihedral[3] = 2;
  dihedral[4] = 3;
  dihedral[5] = 4;
  if (MultiIntegerOption(argc, argv, "-n", &test, dihedral)) {
    exit(1);
  }
  if (test != 0) { // -n is present
    number_of_beads = test;
  }

  // Error: wrong number of integers //{{{
  if ((number_of_beads%beads_per_angle) != 0) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "-n");
    RedText(STDERR_FILENO);
    fprintf(stderr, " - number of bead ids must be divisible by six.\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    dihedral[i]--; // ids should start with zero

    // Error - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && dihedral[i] >= MoleculeType[j].nBeads) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "-n");
        RedText(STDERR_FILENO);
        fprintf(stderr, " - index");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%d", dihedral[i]+1);
        RedText(STDERR_FILENO);
        fprintf(stderr, " is large than the number of beads in molecule ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", MoleculeType[j].Name);
        RedText(STDERR_FILENO);
        fprintf(stderr, "\n\n");
        ResetColour(STDERR_FILENO);
        Help(argv[0], true);
        exit(1);
      }
    } //}}}
  }
  // Error - same bead id used in specifying a plane
  for (int i = 0; i < number_of_beads; i += (beads_per_angle/2)) {
    if (dihedral[i] == dihedral[i+1] ||
        dihedral[i] == dihedral[i+2] ||
        dihedral[i+1] == dihedral[i+2]) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "-n");
      RedText(STDERR_FILENO);
      fprintf(stderr, " - a plane must be specified by three different beads (wrong trio: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d %d %d", dihedral[i], dihedral[i+1], dihedral[i+2]);
      RedText(STDERR_FILENO);
      fprintf(stderr, "\n\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
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

    // print command to output file
    putc('#', out);
    PrintCommand(out, argc, argv);

    // print molecule names & bead ids //{{{
    fprintf(out, "# dihedral angles between planes specifief by:");
    for (int j = 0; j < number_of_beads; j += beads_per_angle) {
      fprintf(out, " (%d) %d-%d-%d & %d-%d-%d;", j/beads_per_angle+1, dihedral[j]+1, dihedral[j+1]+1, dihedral[j+2]+1, dihedral[j+3]+1, dihedral[j+4]+1, dihedral[j+5]+1);
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

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    fprintf(stdout, "\nChosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(stdout, " %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
    fprintf(stdout, "\nPlanes for which to calculate angles:\n");
    for (int i = 0; i < number_of_beads; i += beads_per_angle) {
      fprintf(stdout, "  %d-%d-%d &", dihedral[i], dihedral[i+1], dihedral[i+2]);
      fprintf(stdout, " %d-%d-%d\n", dihedral[i+3], dihedral[i+4], dihedral[i+5]);
    }
    putchar('\n');
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vtf, vcf, struct_lines); //}}}

  count = SkipCoorSteps(vcf, input_coor, Counts, start, silent);

  // allocate array for distribution of angles //{{{
  double avg_angle[Counts.TypesOfMolecules][number_of_beads/beads_per_angle];
  double *distr[Counts.TypesOfMolecules][number_of_beads/beads_per_angle];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
      avg_angle[i][j] = 0;
      distr[i][j] = calloc(bins, sizeof(double));
    }
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int count_step = 0; // count calculated timesteps
  int count_vcf = start - 1; // count timesteps from the beginning
  while (true) {
    count_step++;
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    ReadVcfCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

    // join molecules if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // calculate dihedral angles //{{{
    double angle[Counts.Molecules][number_of_beads/beads_per_angle];
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type = Molecule[i].Type;
      if (MoleculeType[mol_type].Use) {

        // calculate normal vectors to specified planes
        // first plane is given by 0 1 2 and second plane by 1 2 3
        for (int j = 0; j < number_of_beads; j += beads_per_angle) {
          // plane normals //{{{
          VECTOR u[2], v[2], n[2];
          // vectors in first plane (points 0 1 2): u=1-2; v=1-2
          u[0].x = Bead[Molecule[i].Bead[dihedral[j+1]]].Position.x - Bead[Molecule[i].Bead[dihedral[j+2]]].Position.x;
          u[0].y = Bead[Molecule[i].Bead[dihedral[j+1]]].Position.y - Bead[Molecule[i].Bead[dihedral[j+2]]].Position.y;
          u[0].z = Bead[Molecule[i].Bead[dihedral[j+1]]].Position.z - Bead[Molecule[i].Bead[dihedral[j+2]]].Position.z;
          v[0].x = Bead[Molecule[i].Bead[dihedral[j+0]]].Position.x - Bead[Molecule[i].Bead[dihedral[j+2]]].Position.x;
          v[0].y = Bead[Molecule[i].Bead[dihedral[j+0]]].Position.y - Bead[Molecule[i].Bead[dihedral[j+2]]].Position.y;
          v[0].z = Bead[Molecule[i].Bead[dihedral[j+0]]].Position.z - Bead[Molecule[i].Bead[dihedral[j+2]]].Position.z;
          // normal in first plane
          n[0].x = u[0].y * v[0].z - u[0].z * v[0].y;
          n[0].y = u[0].z * v[0].x - u[0].x * v[0].z;
          n[0].z = u[0].x * v[0].y - u[0].y * v[0].x;
          // vectors in second plane (points 3 4 5): u=4-3; v=5-3
          u[1].x = Bead[Molecule[i].Bead[dihedral[j+4]]].Position.x - Bead[Molecule[i].Bead[dihedral[j+3]]].Position.x;
          u[1].y = Bead[Molecule[i].Bead[dihedral[j+4]]].Position.y - Bead[Molecule[i].Bead[dihedral[j+3]]].Position.y;
          u[1].z = Bead[Molecule[i].Bead[dihedral[j+4]]].Position.z - Bead[Molecule[i].Bead[dihedral[j+3]]].Position.z;
          v[1].x = Bead[Molecule[i].Bead[dihedral[j+5]]].Position.x - Bead[Molecule[i].Bead[dihedral[j+3]]].Position.x;
          v[1].y = Bead[Molecule[i].Bead[dihedral[j+5]]].Position.y - Bead[Molecule[i].Bead[dihedral[j+3]]].Position.y;
          v[1].z = Bead[Molecule[i].Bead[dihedral[j+5]]].Position.z - Bead[Molecule[i].Bead[dihedral[j+3]]].Position.z;
          // normal in second plane
          n[1].x = u[1].y * v[1].z - u[1].z * v[1].y;
          n[1].y = u[1].z * v[1].x - u[1].x * v[1].z;
          n[1].z = u[1].x * v[1].y - u[1].y * v[1].x; //}}}

          // calculate angle between the two normals //{{{
          double size[2];
          size[0] = Length(n[0]);
          size[1] = Length(n[1]);
          double scalar = n[0].x * n[1].x + n[0].y * n[1].y + n[0].z * n[1].z;
          angle[i][j/beads_per_angle] = acos(scalar / (size[0] * size[1])); // in rad

          // too close to 180 or 0 degrees to call //{{{
          if (scalar < 0 && fabs((scalar/(size[0]*size[1])+1)) <= 0.00000001) {
            angle[i][j/beads_per_angle] = PI - 0.00000001;
          } else if (scalar > 0 && fabs((scalar/(size[0]*size[1])-1)) <= 0.00000001) {
            angle[i][j/beads_per_angle] = 0;
          } //}}}

          angle[i][j/beads_per_angle] *= 180 / PI; // in degrees //}}}

          // add to average
          avg_angle[mol_type][j/beads_per_angle] += angle[i][j/beads_per_angle];

          int k = angle[i][j/beads_per_angle] / width;
          if (k < bins) {
            distr[mol_type][j/beads_per_angle][k]++;
          } else {
            YellowText(STDERR_FILENO);
            fprintf(stdout, "\nWarning - weird angle: ");
            CyanText(STDERR_FILENO);
            fprintf(stdout, "%lf", angle[i][j/beads_per_angle]);
            YellowText(STDERR_FILENO);
            fprintf(stdout, " degrees\n");
            ResetColour(STDERR_FILENO);
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
        int mol_type = Molecule[i].Type;
        if (MoleculeType[mol_type].Use) {
          for (int j = 0; j < number_of_beads; j += beads_per_angle){
            fprintf(out, " %10.6f", 180-angle[i][j/beads_per_angle]); // write angle between planes, not normals
          }
        }
      }
      putc('\n', out);

      // write stuff
      fclose(out);
    } //}}}

    if (end == count_vcf) {
      break;
    }
    // if there's no additional timestep, exit the while loop
    bool rubbish; // not used
    if (ReadTimestepPreamble(&rubbish, input_coor, vcf, &stuff, false) == -1) {
      break;
    }
  }
  fclose(vcf);

  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // write distribution of angles //{{{
  // open output file for appending //{{{
  FILE *out;
  if ((out = fopen(output_distr, "w")) == NULL) {
    ErrorFileOpen(output_distr, 'w');
    exit(1);
  } //}}}

  // print command to output file
  putc('#', out);
  PrintCommand(out, argc, argv);

  // print molecule names & bead ids //{{{
  fprintf(out, "# dihedral angles between planes specifief by:");
  for (int j = 0; j < number_of_beads; j += beads_per_angle) {
    fprintf(out, " (%d) %d-%d-%d & %d-%d-%d;", j/beads_per_angle+1, dihedral[j]+1, dihedral[j+1]+1, dihedral[j+2]+1, dihedral[j+3]+1, dihedral[j+4]+1, dihedral[j+5]+1);
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
  fprintf(out, "# dihedral angles between planes specifief by:");
  for (int j = 0; j < number_of_beads; j += beads_per_angle) {
    fprintf(out, " (%d) %d-%d-%d & %d-%d-%d;", j/beads_per_angle+1, dihedral[j]+1, dihedral[j+1]+1, dihedral[j+2]+1, dihedral[j+3]+1, dihedral[j+4]+1, dihedral[j+5]+1);
  }
  putc('\n', out);
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
  putc('\n', out);
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
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  free(output);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
      free(distr[i][j]);
    }
  } //}}}

  return 0;
}
