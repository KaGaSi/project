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
        strcmp(argv[i], "--script") != 0 &&
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
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
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

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {

    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (mol_type == -1) {
      fprintf(stderr, "Error: molecule '%s' does not exist in FIELD\n\n", argv[count]);
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
    fprintf(stderr, "\nError: '-n' option - number of bead ids must be divisible by six.\n");
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    dihedral[i]--; // ids should start with zero

    // Error - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && dihedral[i] >= MoleculeType[j].nBeads) {
        fprintf(stderr, "\nError: '-n' option - %d is larger than the number of beads in molecule %s (%d beads)\n\n", dihedral[i]+1, MoleculeType[j].Name, MoleculeType[j].nBeads);
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
      fprintf(stderr, "\nError: '-n' option - a plane must be specified by three different beads (wrong trio: %d %d %d)\n", dihedral[i], dihedral[i+1], dihedral[i+2]);
      exit(1);
    }
  } //}}}

  // '-a' option - write angles for all molecules //{{{
  char *output = calloc(LINE,sizeof(char *));
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
      ErrorCoorRead(input_coor, test, count_vcf, stuff, input_vsf);
      exit(1);
    } //}}}

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
            fprintf(stdout, "\nWarning - weird angle: %lf degrees\n", angle[i][j/beads_per_angle]);
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

    if (end == count_vcf)
      break;
  }
  fclose(vcf);

  // print last sep? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rLast Step: %d\n", count_vcf);
    }
  } //}}}
  //}}}

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
  fprintf(out, "# dihedral angles between planes specifief by:");
  for (int j = 0; j < number_of_beads; j += beads_per_angle) {
    fprintf(out, " (%d) %d-%d-%d & %d-%d-%d;", j/beads_per_angle+1, dihedral[j]+1, dihedral[j+1]+1, dihedral[j+2]+1, dihedral[j+3]+1, dihedral[j+4]+1, dihedral[j+5]+1);
  }
  putc('\n', out);
  fprintf(out, "# columns: (1) angle [deg];");
  int j = 2;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if ((number_of_beads/beads_per_angle) == 1) {
        fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_beads/beads_per_angle-1, MoleculeType[i].Name);
      }
      j += number_of_beads / beads_per_angle;
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
  j = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if ((number_of_beads/beads_per_angle) == 1) {
        fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_beads/beads_per_angle-1, MoleculeType[i].Name);
      }
      j += number_of_beads / beads_per_angle;
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
