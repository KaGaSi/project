#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
BondLength utility calculates distribution of bond lengths for all bead type \
pairs in specified molecule(s). It can also calculate distribution of distances \
between any two beads in a molecule.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> <mol name(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>                 input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <width>                 width of a single distribution bin\n");
  fprintf(ptr, "   <output>                output file with distribution of bond lengths\n");
  fprintf(ptr, "   <mole name(s)>          molecule name(s) to calculate bond lengths for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>            starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>             ending timestep for calculation\n");
  fprintf(ptr, "      -d <out> <ints>      write distribution of distances between specified beads in the molecule to file <out>\n");
  fprintf(ptr, "      -w <double>          warn if bond length exceeds <double> (default: half a box length)\n");
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
        strcmp(argv[i], "-d") != 0 &&
        strcmp(argv[i], "-w") != 0) {

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
  if (!IsPosDouble(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - file name with bond length distribution //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

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
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n\n", start, end);
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

  // <molecule names> - names of molecule types to use //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindMoleculeType(argv[count], Counts, MoleculeType);
    if (type == -1) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", input_coor);
      fprintf(stderr, " - non-existent molecule type \033[1;33m%s\033[1;31m\n\n", argv[count]);
      fprintf(stderr, "\033[0m");
      ErrorMoleculeType(Counts, MoleculeType);
      exit(1);
    }
    MoleculeType[type].Use = true;
  } //}}}

  // '-d' option - specify bead ids to calculate distance between //{{{
  int bead[100] = {0}, number_of_beads = 0;
  char output_d[LINE];
  output_d[0] = '\0';
  if (FileIntsOption(argc, argv, "-d", bead, &number_of_beads, output_d)) {
    exit(1);
  }

  // if '-d' is present, but without numbers - use first and last for each molecule
  if (output_d[0] != '\0' && number_of_beads == 0) {
    number_of_beads = 2;
    bead[0] = 1;
    bead[1] = 1000000;
  }

  // Error: wrong number of integers //{{{
  if (output_d[0] != '\0' && (number_of_beads%2) != 0) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: '-d' option - number of bead ids must be even\n\n");
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  // Error: same bead ids //{{{
  for (int i = 0; i < number_of_beads; i += 2) {
    if (bead[i] == bead[i+1] || bead[i] == 0 || bead[i+1] == 0) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: '-d' option - the bead indices must be non-zero and different\n");
      fprintf(stderr, "\033[0m");
      exit(1);
    }
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero (or -1 if none specified)
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // '-w' option - bond length warning //{{{
  double warn = Min3(BoxLength.x, BoxLength.y, BoxLength.z);
  if (DoubleOption(argc, argv, "-w", &warn)) {
    exit(1);
  } //}}}

  // number of bins
  int bins = Max3(BoxLength.x, BoxLength.y, BoxLength.z) / (2 * width);

  // arrays for distributions //{{{
  double *length[Counts.TypesOfMolecules][Counts.TypesOfBeads][Counts.TypesOfBeads];
  double min_max[Counts.TypesOfMolecules][Counts.TypesOfBeads][Counts.TypesOfBeads][2];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        length[i][j][k] = calloc(bins,sizeof(double));
        min_max[i][j][k][0] = 10 * Max3(BoxLength.x, BoxLength.y, BoxLength.z);
        min_max[i][j][k][1] = 0;
      }
    }
  }
  // extra arrays for -d option
  double *distance[Counts.TypesOfMolecules][number_of_beads/2];
  double min_max_d_option[Counts.TypesOfMolecules][number_of_beads/2][2];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/2); j++) {
      distance[i][j] = calloc(bins, sizeof(double));
      min_max_d_option[i][j][0] = 10 * Max3(BoxLength.x, BoxLength.y, BoxLength.z);
      min_max_d_option[i][j][1] = 0;
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

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
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
      exit(1);
    }
  }
  // print starting step? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step: %6d\n", start);
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

    ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

    // join all molecules
    RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

    // calculate bond length //{{{
    // go through all molecules
    for (int i = 0; i < Counts.Molecules; i++) {
      int mtype = Molecule[i].Type;
      if (MoleculeType[mtype].Use) { // use only specified molecule types
        for (int j = 0; j < MoleculeType[mtype].nBonds; j++) {

          // bead ids in the bond //{{{
          int id1 = Molecule[i].Bead[MoleculeType[mtype].Bond[j][0]];
          int id2 = Molecule[i].Bead[MoleculeType[mtype].Bond[j][1]]; //}}}
          int btype1 = Bead[id1].Type;
          int btype2 = Bead[id2].Type;

          // bond length //{{{
          VECTOR bond;
          bond.x = Bead[id1].Position.x - Bead[id2].Position.x;
          bond.y = Bead[id1].Position.y - Bead[id2].Position.y;
          bond.z = Bead[id1].Position.z - Bead[id2].Position.z;

          bond.x = Length(bond); //}}}

          // warn if bond is too long //{{{
          if (bond.x > warn) {
            fprintf(stderr, "\033[1;33m");
            fprintf(stderr, "\nWarning: bond longer than %lf\n", warn);
            fprintf(stderr, " Step: %d\n", count_vcf);
            fprintf(stderr, " Beads: %6d (%s): %lf %lf %lf\n", Bead[id1].Index,
                                                               BeadType[btype1].Name,
                                                               Bead[id1].Position.x,
                                                               Bead[id1].Position.y,
                                                               Bead[id1].Position.z);
            fprintf(stderr, "        %6d (%s): %lf %lf %lf\n", Bead[id2].Index,
                                                               BeadType[btype2].Name,
                                                               Bead[id2].Position.x,
                                                               Bead[id2].Position.y,
                                                               Bead[id2].Position.z);
            fprintf(stderr, " Bond length: %lf\n", bond.x);
            fprintf(stderr, "\033[0m");
          } //}}}

          // btype1 must be lower then btype2
          if (btype1 > btype2) {
            int swap = btype1;
            btype1 = btype2;
            btype2 = swap;
          }

          // mins & maxes //{{{
          if (bond.x < min_max[mtype][btype1][btype2][0]) {
            min_max[mtype][btype1][btype2][0] = bond.x;
          } else if (bond.x > min_max[mtype][btype1][btype2][1]) {
            min_max[mtype][btype1][btype2][1] = bond.x;
          } //}}}

          int k = bond.x / width;
          if (k < bins) {
            length[mtype][btype1][btype2][k]++;
          }
        }
      }
    } //}}}

    // calculate distance (-d option) //{{{
    if (output_d[0] != '\0') {
      for (int i = 0; i < Counts.Molecules; i++) {
        int mtype = Molecule[i].Type;
        if (MoleculeType[mtype].Use) { // use only specified molecule types
          for (int j = 0; j < number_of_beads; j += 2) {

            // bead ids the distance //{{{
            int id1, id2;
            // use first molecule bead if bead index too high or -1
            if (bead[j] == -1 || bead[j] >= MoleculeType[mtype].nBeads) {
              id1 = Molecule[i].Bead[0];
            } else{ // use specified index otherwise
              id1 = Molecule[i].Bead[bead[j]];
            }
            // use last molecule bead if bead index too high or -1
            if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[mtype].nBeads) {
              id2 = Molecule[i].Bead[MoleculeType[mtype].nBeads-1];
            } else{ // use specified index otherwise
              id2 = Molecule[i].Bead[bead[j+1]];
            } //}}}

            // distance //{{{
            VECTOR dist;
            dist.x = Bead[id1].Position.x - Bead[id2].Position.x;
            dist.y = Bead[id1].Position.y - Bead[id2].Position.y;
            dist.z = Bead[id1].Position.z - Bead[id2].Position.z;

            dist.x = Length(dist); //}}}

            // mins & maxes //{{{
            if (dist.x < min_max_d_option[mtype][j/2][0]) {
              min_max_d_option[mtype][j/2][0] = dist.x;
            } else if (dist.x > min_max_d_option[mtype][j/2][1]) {
              min_max_d_option[mtype][j/2][1] = dist.x;
            } //}}}

            int k = dist.x / width;
            if (k < bins) {
              distance[mtype][j/2][k]++;
            }
          }
        }
      }
    } //}}}

    if (end == count_vcf)
      break;
  }

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count_vcf);
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

  // count total number of bonds in molecules //{{{
  int bonds[Counts.TypesOfMolecules][Counts.TypesOfBeads][Counts.TypesOfBeads]; //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] = 0;
        }
      }
    }
  } //}}}

  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] += length[i][j][k][l];
        }
      }
    }
  } //}}}

  // write distribution of bond lengths //{{{
  // open output file for writing //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    ErrorFileOpen(output, 'w');
    exit(1);
  } //}}}

  // print command to output file
  putc('#', out);
  PrintCommand(out, argc, argv);

  // print first line of output file - molecule names and beadtype pairs //{{{
  fprintf(out, "# (1) distance;");
  count = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, " %s molecule:", MoleculeType[i].Name);
      for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
        for (int k = j; k < MoleculeType[i].nBTypes; k++) {
          if (bonds[i][MoleculeType[i].BType[j]][MoleculeType[i].BType[k]] > 0) {
            fprintf(out, " (%d) %s-%s", ++count, BeadType[MoleculeType[i].BType[j]].Name, BeadType[MoleculeType[i].BType[k]].Name);
            // add semicolon if this is the last pair for this molecule, add comma otherwise
            if (k == (MoleculeType[i].nBTypes-1)) {
              putc(';', out);
            } else {
              putc(',', out);
            }
          }
        }
      }
    }
  }
  putc('\n', out); //}}}

  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%7.4f", width*(2*i+1)/2);

    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {

        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < MoleculeType[j].nBTypes; k++) {
          for (int l = k; l < MoleculeType[j].nBTypes; l++) {

            int btype1 = MoleculeType[j].BType[k];
            int btype2 = MoleculeType[j].BType[l];

            // btype1 must be lower then btype2
            if (btype1 > btype2) {
              int swap = btype1;
              btype1 = btype2;
              btype2 = swap;
            }

            if (bonds[j][btype1][btype2] > 0) {
              fprintf(out, "%10f", length[j][btype1][btype2][i]/bonds[j][btype1][btype2]);
            }
          }
        }
      }
    }
    putc('\n', out);
  } //}}}

  // write mins and maxes //{{{
  // legend line
  fprintf(out, "# mins(odd columns)/maxes(even columns) -");
  count = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, " %s molecule:", MoleculeType[i].Name);

      for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
        for (int k = j; k < MoleculeType[i].nBTypes; k++) {
          if (bonds[i][MoleculeType[i].BType[j]][MoleculeType[i].BType[k]] > 0) {
            fprintf(out, " (%d) %s-%s", count, BeadType[MoleculeType[i].BType[j]].Name, BeadType[MoleculeType[i].BType[k]].Name);
            count += 2;
            if (k == (MoleculeType[i].nBTypes-1)) {
              if (i != (Counts.TypesOfMolecules-1)) {
                putc(';', out);
              }
            } else {
              putc(',', out);
            }
          }
        }
      }
    }
  }
  putc('\n', out);
  // data line
  putc('#', out);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        if (min_max[i][j][k][1] > 0) {
          fprintf(out, " %lf %lf", min_max[i][j][k][0], min_max[i][j][k][1]);
        }
      }
    }
  }
  putc('\n', out); //}}}
  fclose(out); //}}}

  // write distribution of distances from '-d' option //{{{
  if (output_d[0] != '\0') {
    // open output file for appending //{{{
    FILE *out;
    if ((out = fopen(output_d, "w")) == NULL) {
      ErrorFileOpen(output_d, 'w');
      exit(1);
    } //}}}

    // print command to output file
    putc('#', out);
    PrintCommand(out, argc, argv);

    // print the first line of output file - molecule names with bead order //{{{
    fprintf(out, "# bead order in molecule(s) -");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %s:", MoleculeType[i].Name);
        for (int j = 0; j < MoleculeType[i].nBeads; j++) {
          fprintf(out, " %s", BeadType[MoleculeType[i].Bead[j]].Name);
        }
        putc(';', out);
      }
    }
    putc('\n', out); //}}}

    // print the second line of output file - molecule names and indices with column numbers //{{{
    fprintf(out, "# (1) distance;");
    count = 1;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %s:", MoleculeType[i].Name);

        for (int j = 0; j < number_of_beads; j += 2) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (bead[j] == -1 || bead[j] >= MoleculeType[i].nBeads) {
            id1 = 1;
          } else{ // use specified index otherwise
            id1 = bead[j]+1;
          }
          // use last molecule bead if bead index too high or -1
          if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[i].nBeads) {
            id2 = MoleculeType[i].nBeads;
          } else{ // use specified index otherwise
            id2 = bead[j+1]+1;
          } //}}}

          fprintf(out, " (%d) %d-%d", ++count, id1, id2);
          // add semicolon if this is the last pair for this molecule, add comma otherwise
          if (j == (number_of_beads-2)) {
            putc(';', out);
          } else {
            putc(',', out);
          }
        }
      }
    }
    putc('\n', out); //}}}

    // write the distribution to output file //{{{
    for (int i = 0; i < bins; i++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2);

      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        if (MoleculeType[j].Use) {
          for (int k = 0; k < number_of_beads; k += 2) {
            fprintf(out, " %10f", (double)(distance[j][k/2][i])/(steps*MoleculeType[j].Number));
          }
        }
      }
      putc('\n', out);
    } //}}}

    // write mins and maxes //{{{
    // legend line
    fprintf(out, "# mins(odd columns)/maxes(even columns) -");
    count = 1;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %s:", MoleculeType[i].Name);

        for (int j = 0; j < number_of_beads; j += 2) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (bead[j] == -1 || bead[j] >= MoleculeType[i].nBeads) {
            id1 = 1;
          } else{ // use specified index otherwise
            id1 = bead[j]+1;
          }
          // use last molecule bead if bead index too high or -1
          if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[i].nBeads) {
            id2 = MoleculeType[i].nBeads;
          } else{ // use specified index otherwise
            id2 = bead[j+1]+1;
          } //}}}

          fprintf(out, " (%d) %d-%d", count, id1, id2);
          count += 2;
          // add semicolon if this is the last pair for this molecule, add comma otherwise
          if (j == (number_of_beads-2)) {
            putc(';', out);
          } else {
            putc(',', out);
          }
        }
      }
    }
    putc('\n', out);

    // data line
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        for (int j = 0; j < number_of_beads; j+=2) {
          fprintf(out, " %lf %lf", min_max_d_option[i][j/2][0], min_max_d_option[i][j/2][1]);
        }
      }
    }
    putc('\n', out); //}}}

    fclose(out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        free(length[i][j][k]);
      }
    }
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/2); j++) {
      free(distance[i][j]);
    }
  } //}}}

  return 0;
}
