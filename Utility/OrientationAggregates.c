#include "../AnalysisTools.h"

// TODO: orientation for aggregates? Is there a point? Well, there is - e.g., under shear...

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
        REWRITE!!!\n\n\
DensityAggregates utility calculates radial bead density for aggregates \
of given size(s) from their centre of mass. For beads in molecules, \
it takes into account only beads from the current aggregate, not from \
any other aggregate. The utility does not check what molecule type a \
given bead belongs to; therefore, if the same bead type appears in \
more molecule types, the resulting density for that bead type will be \
averaged without regard for the various types of molecule it comes from \
(i.e., if 'mol1' and 'mol2' both both contain bead 'A', there will be \
only one column for 'A' bead type).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <input.agg> <width> <output> <agg size(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <input.agg>       input agg file\n");
  fprintf(ptr, "   <width>           width of a single bin of the distribution\n");
  fprintf(ptr, "   <output>          output file with orientation (automatic '#.txt' ending)\n");
  fprintf(ptr, "   <agg size(s)>     aggregate size(s) to calculate orientation for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -n <ints>      bead indices (multiple of 2 <ints>; default: 1 2)\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> molecule types in an aggregate\n");
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
  int req_args = 5; //}}}

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
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-m") != 0) {

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

  // <input.agg> - input agg file //{{{
  char input_agg[LINE];
  strcpy(input_agg, argv[++count]);

  // test if <input.agg> ends with '.agg'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".agg");
  if (ErrorExtension(input_agg, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (!IsPosDouble(argv[++count])) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.txt> - filename with bead densities //{{{
  char output_rho[LINE];
  strcpy(output_rho, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // number of bins to average //{{{
  int avg = 1;
  if (IntegerOption(argc, argv, "-n", &avg)) {
    exit(1);
  } //}}}

  int start, end;
  StartEndTime(argc, argv, &start, &end); //}}}

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

  // '-m' option //{{{
  int *specific_moltype_for_size = malloc(Counts.TypesOfMolecules*sizeof(int));
  // all are to be used without '-m' option
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", &specific_moltype_for_size, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // '-n' option - specify bead ids //{{{
  int bead[100] = {0}, // specified bead indices
      number_of_angles, // total number of angles to calculate TODO: why is it called this?
      number_of_beads = 2, // number of parameters to -n
      beads_per_angle = 2, // the numbers must come in pairs
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
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "-n");
    RedText(STDERR_FILENO);
    fprintf(stderr, " - number of bead ids must be divisible by two\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    // Warning - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && bead[i] >= MoleculeType[j].nBeads) {
//      if ((i%2) == 1 || bead[i+1] <= MoleculeType[j].nBeads) { // warning if one is higher
//        fprintf(stderr, "\033[1;33m");
//        fprintf(stderr, "\nWarning: '-n' option - %d is larger than the number of beads in %s; using %d instead\n", bead[i], MoleculeType[j].Name, MoleculeType[j].nBeads);
//        fprintf(stderr, "\033[0m");
//      } else if ((i%2) == 0 && bead[i+1] > MoleculeType[j].nBeads) {
//        fprintf(stderr, "\033[1;33m");
//        fprintf(stderr, "\nWarning: '-n' option - both %d and %d are larger than the number of beads in %s; this pair is ignored\n", bead[i], bead[i+1], MoleculeType[j].Name);
//        fprintf(stderr, "\033[0m");
//      }
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "-n");
        RedText(STDERR_FILENO);
        fprintf(stderr, " - index ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%d", bead[i]);
        RedText(STDERR_FILENO);
        fprintf(stderr, " is larger than the number of beads in molecule ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", MoleculeType[j].Name);
        fprintf(stderr, "\n\n");
        ResetColour(STDERR_FILENO);
      }
    } //}}}
    // Error - two same beads //{{{
    if ((i%2) == 0 && bead[i] == bead[i+1]) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m-n\033[1;31m option - each pair must contain two different bead ids (invalid: %d %d)\n", bead[i], bead[i+1]);
      fprintf(stderr, "\033[0m");
      Help(argv[0], true);
      exit(1);
    } //}}}
    bead[i]--; // ids should start with zero
  } //}}}

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  int **agg_mols = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
    agg_mols[i] = calloc(Counts.TypesOfMolecules, sizeof(int));
  }

  int aggs = 0; // number of aggregate sizes

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (!IsInteger(argv[count]) || atoi(argv[count]) == 0) {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}

    agg_sizes[aggs][0] = atoi(argv[count]);

    // write initial stuff to output density file //{{{
    FILE *out;
    char str[LINE];
    strcpy(str, output_rho);

    char str2[1030];
    sprintf(str2, "%s%d.txt", str, agg_sizes[aggs][0]);
    strcpy(str, str2);
    if ((out = fopen(str, "w")) == NULL) {
      ErrorFileOpen(str, 'w');
      exit(1);
    }

    // print command to output file
    putc('#', out);
    PrintCommand(out, argc, argv);

    // print first lines of output file - molecule names and beadtype pairs //{{{
    fprintf(out, "# columns: (1) Orientation order parameter");
    int count_print = 1;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, "\n#          %s:", MoleculeType[i].Name);
      for (int j = 0; j < number_of_beads; j += 2) {
        // are ids right?
        switch(test_ids(bead[j], bead[j+1], MoleculeType[i].nBeads)) {
          case 0: // bead[k] & bead[k+1] too high
            continue;
          case 1: // bead[k] too high & bead[k+1] okay
            if ((bead[j+1]+1) != MoleculeType[i].nBeads) {
              fprintf(out, " (%d) %d-%d,", ++count_print, MoleculeType[i].nBeads, bead[j+1]+1);
            }
            break;
          case 2: // bead[k] okay & bead[k+1] too high
            if ((bead[j]+1) != MoleculeType[i].nBeads) {
              fprintf(out, " (%d) %d-%d,", ++count_print, bead[j]+1, MoleculeType[i].nBeads);
            }
            break;
          case 3: // bead[k] & bead[k+1] okay
            fprintf(out, " (%d) %d-%d,", ++count_print, bead[j]+1, bead[j+1]+1);
            break;
        }
      }
    }
    putc('\n', out); //}}}

    fclose(out); //}}}

    aggs++; // number of aggregate sizes
  } //}}}

  double distance; // <distance> parameter from Aggregate command
  int contacts; // <contacts> parameter from Aggregate command - not used here
  ReadAggCommand(BeadType, Counts, input_coor, input_agg, &distance, &contacts);

  // open input aggregate file and skip the first lines (Aggregate command & blank line) //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  char line[LINE];
  fgets(line, sizeof(line), agg);
  fgets(line, sizeof(line), agg); //}}}

  // number of bins
  int bins = 1.5 / width;

  // array for distributions //{{{
  double *orientation[Counts.TypesOfMolecules][number_of_angles][aggs],
         avg_orientation[Counts.TypesOfMolecules][number_of_angles][aggs];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < number_of_angles; j++) {
      for (int k = 0; k < aggs; k++) {
        orientation[i][j][k] = calloc(bins, sizeof(double));
        avg_orientation[i][j][k] = 0;
      }
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
    fprintf(stdout, "Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      fprintf(stdout, " %d", agg_sizes[i][0]);
    }
    putchar('\n');
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  count = SkipCoorAggSteps(vcf, input_coor, agg, input_agg, Counts, start, silent);

  // use z axis as a normal TODO: option for axis //{{{
  VECTOR normal;
  normal.x = 0;
  normal.y = 0;
  normal.z = 1; //}}}

  // main loop //{{{
  count = 0; // count timesteps in the main loop
  int count_vcf = start - 1; // count timesteps from the beginning
  while (true) {
    count++;
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    ReadAggregates(agg, input_agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index);
    ReadVcfCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);
    RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

    // calculate orientation //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {
      // test if aggregate 'i' should be used //{{{
      int size = 0;
      // agg size = number of molecules of type 'specific_moltype_for_size'
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type]) {
          size++;
        }
      }
      // is 'size' in provided list?
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j][0] == size) {
          correct_size = j;
        }
      } //}}}
      if (correct_size != -1) {
        agg_sizes[correct_size][1]++;

        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol = Aggregate[i].Molecule[j];
          int mtype = Molecule[mol].Type;

          for (int k = 0; k < number_of_beads; k += 2) {
            int mbeads = MoleculeType[mtype].nBeads;
            int id1 = -2, id2 = -1;
            // check provided ids for given molecule type //{{{
            switch(test_ids(bead[k], bead[k+1], mbeads)) {
              case 0: // bead[k] & bead[k+1] too high
                continue;
              case 1: // bead[k] too high & bead[k+1] okay
                id1 = Molecule[mol].Bead[mbeads-1];
                id2 = Molecule[mol].Bead[bead[k+1]];
                break;
              case 2: // bead[k] okay & bead[k+1] too high
                id1 = Molecule[mol].Bead[bead[k]];
                id2 = Molecule[mol].Bead[mbeads-1];
                break;
              case 3: // bead[k] & bead[k+1] okay
                id1 = Molecule[mol].Bead[bead[k]];
                id2 = Molecule[mol].Bead[bead[k+1]];
                break;
            } //}}}
            // if 'k' or 'k+1' is too high, it can lead to both being the highest id in molecule
            if (id1 == id2) {
              continue;
            }
            double length = SQR(Bead[id1].Position.x - Bead[id2].Position.x) +
                            SQR(Bead[id1].Position.y - Bead[id2].Position.y) +
                            SQR(Bead[id1].Position.z - Bead[id2].Position.z);
            length = sqrt(length);
            VECTOR unit;
            unit.x = (Bead[id1].Position.x - Bead[id2].Position.x) / length;
            unit.y = (Bead[id1].Position.y - Bead[id2].Position.y) / length;
            unit.z = (Bead[id1].Position.z - Bead[id2].Position.z) / length;

            double value = unit.x * normal.x + unit.y * normal.y + unit.z * normal.z; // cos(angle)
            value = 0.5 * (3 * SQR(value) - 1);

            // add to overall averages
            avg_orientation[mtype][k/2][correct_size] += value;

            // place into appropriate distribution bin
            int bin = (value + 0.5) / width; // +0.5, because its range is <-0.5,1>
            if (k < bins) {
              orientation[mtype][k/2][correct_size][bin]++;
            } else if (k == bins) { // in case of value=1 exactly
              orientation[mtype][k/2][correct_size][bins-1]++;
            }
          }
        }

        // count mol types in the aggregate
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          agg_mols[correct_size][Molecule[Aggregate[i].Molecule[j]].Type]++;
        }
      }
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
  fclose(agg);

  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // write orientation to output file(s) //{{{
  for (int temp = 0; temp < aggs; temp++) {
    FILE *out;
    char str[LINE];
    // assemble correct name
    sprintf(str, "%s", output_rho);
    char str2[1030];
    sprintf(str2, "%s%d.txt", str, agg_sizes[temp][0]);
    strcpy(str, str2);
    if ((out = fopen(str, "a")) == NULL) {
      ErrorFileOpen(str, 'a');
      exit(1);
    }

    // calculate normalisation factor
    long int normalisation[Counts.TypesOfMolecules][number_of_beads/2];
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      for (int j = 0; j < number_of_beads; j+=2) {
        normalisation[i][j/2] = 0;
        for (int k = 0; k < bins; k++) {
          normalisation[i][j/2] += orientation[i][j/2][temp][k];
        }
      }
    }

    // write distribution to output file //{{{
    for (int i = 0; i < bins; i++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2-0.5); // orientation order parameter: <-0.5,1>
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
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
              fprintf(out, "%10f", orientation[j][k/2][temp][i]/normalisation[j][k/2]);
              break;
          }
        }
      }
      putc('\n', out);
    } //}}}

    // write overall averages //{{{
    // legend line
    fprintf(out, "# Average orientation order parameter");
    count = 0;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
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
    putc('\n', out);
    // data line
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      for (int j = 0; j < number_of_beads; j += 2) {
        int mbeads = MoleculeType[i].nBeads;
        // check provided ids for given molecule type
        if (bead[j] >= mbeads &&
            bead[j+1] >= mbeads) { // both higher than number beads - skip
          continue;
        }
        fprintf(out, "%10f", avg_orientation[i][j/2][temp]/agg_sizes[temp][1]);
      }
    }
    putc('\n', out); //}}}

//  // calculate rdf //{{{
//  for (int j = 0; j < (bins-avg); j++) {

//    // calculate volume of every shell that will be averaged
//    double shell[avg];
//    for (int k = 0; k < avg; k++) {
//      shell[k] = 4 * PI * CUBE(width) * (CUBE(j+k+1) - CUBE(j+k)) / 3;
//    }

//    fprintf(out, "%.2f", width*(j+0.5*avg));

//    for (int k = 0; k < Counts.TypesOfBeads; k++) {
//      double temp_rdp = 0, temp_number = 0,
//             temp_rdp_err = 0, temp_number_err = 0;

//      // sum up rdfs from all shells to be averaged
//      for (int l = 0; l < avg; l++) {
//        temp_rdp += rho[k][i][j+l] / (shell[l] * agg_sizes[i][1]);
//        temp_rdp_err += rho_2[k][i][j+l] / (shell[l] * agg_sizes[i][1]);
//        temp_number += rho[k][i][j+l] / agg_sizes[i][1];
//        temp_number_err += rho_2[k][i][j+l] / agg_sizes[i][1];
//      }

//      temp_rdp_err = sqrt(temp_rdp_err - temp_rdp);
//      temp_number_err = sqrt(temp_number_err - temp_number);

//      // print average value to output file
//      fprintf(out, " %10f %10f", temp_rdp/avg, temp_rdp_err/avg);
//      fprintf(out, " %10f %10f", temp_number/avg, temp_number_err/avg);
//    }
//    putc('\n',out);
//  } //}}}

//  fprintf(out, "# (1) Number of aggregates; Average numbers of molecules:"); //{{{
//  count = 2;
//  for (int j = 0; j < Counts.TypesOfMolecules; j++) {
//    fprintf(out," (%d) %s", count++, MoleculeType[j].Name);
//    if (i < (Counts.TypesOfMolecules-1)) {
//      putc(';', out);
//    }
//  }
//  putc('\n', out);
//  putc('#', out);
//  fprintf(out," %d", agg_sizes[i][1]);
//  for (int j = 0; j < Counts.TypesOfMolecules; j++) {
//    fprintf(out," %lf", (double)(agg_mols[i][j])/agg_sizes[i][1]);
//  }
//  putc('\n', out); //}}}

    fclose(out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  FreeAggregate(Counts, &Aggregate);
  free(stuff);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
    free(agg_mols[i]);
  }
  free(agg_sizes);
  free(agg_mols);
  free(specific_moltype_for_size); //}}}

  return 0;
}
