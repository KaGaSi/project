#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
OrientOrder utility calculates orientational order parameter for specified \
molecules and specified bead pairs in the molecule. Besides average value, \
it calculates distribution of the parameter. For now, the angle is \
calculated relative to xy plane.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol name(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <width>           width of a single distribution bin\n");
  fprintf(ptr, "   <output>          output file with distribution of the orientation parameter\n");
  fprintf(ptr, "   <mole name(s)>    molecule name(s) to calculate orientational parameter for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -a             axis to use as the normal vector (default: z)\n");
  fprintf(ptr, "      -n <ints>      bead indices (multiple of 2 <ints>; default: 1 2)\n");
  CommonHelp(error);
} //}}}

// function to test if provided indices are too high for a given molecule type //{{{
int test_ids(int id1, int id2, int mbeads) {
  if (id1 >= mbeads &&
      id2 >= mbeads) { // both higher than number beads - skip
    return 0;
  } else if (id1 >= mbeads &&
             id2 < mbeads) { // 'id1' too high, 'id2' okay
    return 1;
  } else if (id1 < mbeads &&
             id2 >= mbeads) { // 'id1' okay, 'id2' too high
    return 2;
  } else { // both are okay
    return 3;
  }
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

  if (argc < req_args) {
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
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-a") != 0 &&
        strcmp(argv[i], "-n") != 0) {

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

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - file name with orientational parameter distribution //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

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
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n\n", start, end);
    exit(1);
  } //}}}

  // axis to use as a normal //{{{
  int test = 'z';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-a") == 0) {
      // Error - missing argument //{{{
      if ((i+1) >= argc) {
        fprintf(stderr, "Missing numeric argument for '-a' option!\n");
        Help(argv[0], true);
        exit(1);
      } //}}}
      test = argv[i+1][0];
    }
  }
  Vector normal;
  switch(test) {
    case 120: // x
      normal.x = 1;
      normal.y = 0;
      normal.z = 0;
      break;
    case 121: // y
      normal.x = 0;
      normal.y = 1;
      normal.z = 0;
      break;
    case 122: // z
      normal.x = 0;
      normal.y = 0;
      normal.z = 1;
      break;
    default:
      fprintf(stderr, "\nError: argument to -a option must be 'x', 'y', or 'z'\n\n");
      Help(argv[0], true);
      exit(1);
  }
//printf("normal=(%lf,%lf,%lf)\n", normal.x, normal.y, normal.z);
//// use later, when arbitrary vector can be used
//double normal_length = sqrt(SQR(normal.x) + SQR(normal.y) + SQR(normal.z));
//normal.x /= normal_length;
//normal.y /= normal_length;
//normal.z /= normal_length; //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <molecule names> - names of molecule types to use //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (type == -1) {
      fprintf(stderr, "\nError: molecule type '%s' is not in %s file\n\n", argv[count], input_coor);
      exit(1);
    }

    MoleculeType[type].Use = true;
  } //}}}

  // '-n' option - specify bead ids //{{{
  int bead[100] = {0}, // specified bead indices
      number_of_angles, // total number of angles to calculate
      number_of_beads = 2, // number of parameters to -n
      beads_per_angle = 2; // the numbers must come in pairs
  test = 0;
  bead[0] = 1; // default ids for angle
  bead[1] = 2;
  if (MultiIntegerOption(argc, argv, "-n", &test, bead)) {
    exit(1);
  }
  if (test != 0) { // -n is present
    number_of_beads = test;
  }
  number_of_angles = number_of_beads / beads_per_angle;

  // Error: wrong number of integers //{{{
  if ((number_of_beads%beads_per_angle) != 0) {
    fprintf(stderr, "\nError: '-n' option - number of bead ids must be divisible by two\n\n");
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    // Error - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && bead[i] > MoleculeType[j].nBeads) {
        if ((i%2) == 1 || bead[i+1] <= MoleculeType[j].nBeads) { // warning if one is higher
          fprintf(stderr, "\nWarning: '-n' option - %d is larger than the number of beads in %s; using %d instead\n", bead[i], MoleculeType[j].Name, MoleculeType[j].nBeads);
        } else if ((i%2) == 0 && bead[i+1] > MoleculeType[j].nBeads) {
          fprintf(stderr, "\nWarning: '-n' option - both %d and %d are larger than the number of beads in %s; this pair is ignored\n", bead[i], bead[i+1], MoleculeType[j].Name);
        }
      }
    } //}}}
    // Error - two same beads //{{{
    if ((i%2) == 0 && bead[i] == bead[i+1]) {
      fprintf(stderr, "\nError: '-n' option - each pair must contain two different bead ids (invalid: %d %d)\n", bead[i], bead[i+1]-1);
      Help(argv[0], true);
      exit(1);
    } //}}}
    bead[i]--; // ids should start with zero
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  Vector BoxLength = GetPBC(vcf, input_coor);

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // number of bins
  int bins = 1.5 / width;

  // array for distributions //{{{
  double *orientation[Counts.TypesOfMolecules][number_of_angles],
         avg_orientation[Counts.TypesOfMolecules][number_of_angles];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_angles; j++) {
      orientation[i][j] = calloc(bins, sizeof(double));
      avg_orientation[i][j] = 0;
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

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
  count = 0;
  int count_vcf = start - 1;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // print step?  //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      ErrorCoorRead(input_coor, test, count_vcf, stuff, input_vsf);
      exit(1);
    } //}}}

    // join all molecules
    RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

    // calculate orientation parameter //{{{
    // go through all molecules
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      if (MoleculeType[mtype].Use) { // use only specified molecule types
        for (int j = 0; j < number_of_beads; j += 2) {
          int mbeads = MoleculeType[mtype].nBeads;
          int id1 = -2, id2 = -1;
          // check provided ids for given molecule type
          switch(test_ids(bead[j], bead[j+1], mbeads)) {
            case 0: // bead[k] & bead[k+1] too high
              continue;
            case 1: // bead[k] too high & bead[k+1] okay
              id1 = Molecule[i].Bead[mbeads-1];
              id2 = Molecule[i].Bead[bead[j+1]];
              break;
            case 2: // bead[k] okay & bead[k+1] too high
              id1 = Molecule[i].Bead[bead[j]];
              id2 = Molecule[i].Bead[mbeads-1];
              break;
            case 3: // bead[k] & bead[k+1] okay
              id1 = Molecule[i].Bead[bead[j]];
              id2 = Molecule[i].Bead[bead[j+1]];
              break;
          }
          // if 'j' or 'j+1' is too high, it can lead to both being the highest id in molecule
          if (id1 == id2) {
            continue;
          }
          double length = SQR(Bead[id1].Position.x - Bead[id2].Position.x) +
                          SQR(Bead[id1].Position.y - Bead[id2].Position.y) +
                          SQR(Bead[id1].Position.z - Bead[id2].Position.z);
          length = sqrt(length);
          Vector unit;
          unit.x = (Bead[id1].Position.x - Bead[id2].Position.x) / length;
          unit.y = (Bead[id1].Position.y - Bead[id2].Position.y) / length;
          unit.z = (Bead[id1].Position.z - Bead[id2].Position.z) / length;

          double value = unit.x * normal.x + unit.y * normal.y + unit.z * normal.z; // cos(angle)
          value = 0.5 * (3 * SQR(value) - 1);

          // add to overall averages
          avg_orientation[mtype][j/2] += value;

          // place into appropriate distribution bin
          int k = (value + 0.5) / width; // +0.5, because its range is <-0.5,1>
          if (k < bins) {
            orientation[mtype][j/2][k]++;
          } else if (k == bins) { // in case of value=1 exactly
            orientation[mtype][j/2][bins-1]++;
          }
        }
      }
    } //}}}

    if (end == count_vcf)
      break;
  }

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rLast Step: %d\n", count_vcf);
    }
  }

  fclose(vcf); //}}}

  // count number of calculated steps //{{{
  int steps;
  if (end != -1) {
    steps = end - start + 1;
  } else {
    steps = count_vcf - start + 1;
  } //}}}

  // write distribution of orientation parameters //{{{
  // open output file for writing //{{{
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

  // print first line of output file - molecule names and beadtype pairs //{{{
  fprintf(out, "# columns: (1) Orientation order parameter");
  count = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, "\n#          %s:", MoleculeType[i].Name);
      for (int j = 0; j < number_of_beads; j += 2) {
        // are ids right?
        switch(test_ids(bead[j], bead[j+1], MoleculeType[i].nBeads)) {
          case 0: // bead[k] & bead[k+1] too high
            continue;
          case 1: // bead[k] too high & bead[k+1] okay
            if ((bead[j+1]+1) != MoleculeType[i].nBeads) {
              fprintf(out, " (%d) %d-%d,", ++count, MoleculeType[i].nBeads, bead[j+1]+1);
            }
            break;
          case 2: // bead[k] okay & bead[k+1] too high
            if ((bead[j]+1) != MoleculeType[i].nBeads) {
              fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, MoleculeType[i].nBeads);
            }
            break;
          case 3: // bead[k] & bead[k+1] okay
            fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, bead[j+1]+1);
            break;
        }
      }
    }
  }
  putc('\n', out); //}}}

  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%7.4f", width*(2*i+1)/2-0.5); // orientation order parameter: <-0.5,1>
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {
        for (int k = 0; k < number_of_beads; k += 2) {
          switch(test_ids(bead[k], bead[k+1], MoleculeType[j].nBeads)) {
            case 0: // bead[k] & bead[k+1] too high
              continue;
            case 1: // bead[k] too high & bead[k+1] okay
              if ((bead[k+1]+1) == MoleculeType[j].nBeads) {
                continue;
              }
            case 2: // bead[k] okay & bead[k+1] too high
              if ((bead[k]+1) == MoleculeType[j].nBeads) {
                continue;
              }
            case 3: // bead[k] & bead[k+1] okay
              fprintf(out, "%10f", orientation[j][k/2][i]/(steps*MoleculeType[j].Number));
              break;
          }
        }
      }
    }
    putc('\n', out);
  } //}}}

  // write averages //{{{
  // legend line
  fprintf(out, "# Average orientation order parameter");
  count = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, "\n#    %s:", MoleculeType[i].Name);
      for (int j = 0; j < number_of_beads; j += 2) {
        if (bead[j] < MoleculeType[i].nBeads && bead[j+1] < MoleculeType[i].nBeads) {
            fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, bead[j+1]+1);
        } else if (bead[j] >= MoleculeType[i].nBeads && bead[j+1] < MoleculeType[i].nBeads) {
            fprintf(out, " (%d) %d-%d,", ++count, MoleculeType[i].nBeads, bead[j+1]+1);
        } else if (bead[j+1] >= MoleculeType[i].nBeads && bead[j] < MoleculeType[i].nBeads) {
            fprintf(out, " (%d) %d-%d,", ++count, bead[j]+1, MoleculeType[i].nBeads);
        }
      }
    }
  }
  putc('\n', out);
  // data line
  putc('#', out);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      for (int j = 0; j < number_of_beads; j += 2) {
        int mbeads = MoleculeType[i].nBeads;
        // check provided ids for given molecule type
        if (bead[j] >= mbeads &&
            bead[j+1] >= mbeads) { // both higher than number beads - skip
          continue;
        }
        fprintf(out, "%10f", avg_orientation[i][j/2]/(steps*MoleculeType[i].Number));
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_angles; j++) {
      free(orientation[i][j]);
    }
  } //}}}

  return 0;
}
