#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
DistrAgg calculates weight and number average aggregation numbers during the \
simulation run as well as overall weight and number distributions and volume \
fractions of aggregates.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input.agg> <output distr file> <output avg file> <options>\n\n", cmd);

  fprintf(ptr, "   <input.agg>            input agg file\n");
  fprintf(ptr, "   <distr file>           output file with distribution of aggregate sizes\n");
  fprintf(ptr, "   <avg file>             output file with per-timestep averages\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -st <int>           starting timestep for calculation (default: 1)\n");
  fprintf(ptr, "      -n <int> <int>      use aggregate sizes in a given range\n");
  fprintf(ptr, "      -m <name(s)>        use number of specified molecule(s) as aggrete size\n");
  fprintf(ptr, "      -x <name(s)>        exclude aggregates containing only specified molecule(s)\n");
  fprintf(ptr, "      --only <name(s)>    use only aggregates composed of a specified molecule\n");
  fprintf(ptr, "      -c <name> <int(s)>  composition distribution of specific aggregate size(s) to <name>\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
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
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "--only") != 0 &&
        strcmp(argv[i], "-c") != 0 ) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
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
    Help(argv[0], true);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(1024,sizeof(char *));
  if (FileOption(argc, argv, "-b", &bonds_file)) {
    exit(0);
  } //}}}

  // output verbosity //{{{
  bool verbose2, silent;
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  VerboseLongOption(argc, argv, &verbose, &verbose2); // more verbose output
  SilentOption(argc, argv, &verbose, &verbose2, &silent); // no output
  bool script = BoolOption(argc, argv, "--script"); // do not use \r & co.
  // }}}

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input.agg> - input agg file //{{{
  char input_agg[1024];
  strcpy(input_agg, argv[++count]);

  // test if <output.agg> filename ends with '.agg' (required by VMD)
  ext = 1;
  extension = malloc(ext*sizeof(char *));
  extension[0] = malloc(5*sizeof(char));
  strcpy(extension[0], ".agg");
  if (!ErrorExtension(input_agg, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // open input file and skip the first two lines //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }

  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // <output distr file> - filename with weight and number distributions //{{{
  char output_distr[1024];
  strcpy(output_distr, argv[++count]); //}}}

  // <output avg file> - filename with weight and number average aggregation numbers //{{{
  char output_avg[1024];
  strcpy(output_avg, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // revers of Bead[].Index
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  char input_coor[1] = {'\0'};
  ReadStructure(input_vsf, input_coor, bonds_file, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // '-n' option - range of aggregation numbers //{{{
  int range_As[2], test = 2;
  range_As[0] = 1;
  range_As[1] = Counts.Molecules;
  if (MultiIntegerOption(argc, argv, "-n", &test, range_As)) {
    exit(1);
  }
  if (test != 2) {
    fprintf(stderr, "\nError: option '-n' requires two numeric arguments\n\n");
    exit(1);
  }

  // make sure first number is smaller
  if (range_As[0] > range_As[1]) {
    int tmp = range_As[0];
    range_As[0] = range_As[1];
    range_As[1] = tmp;
  } //}}}

  // '-m' option //{{{
  int *specific_moltype_for_size = calloc(Counts.TypesOfMolecules,sizeof(int *));
  // all are to be used without '-m' option
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", &specific_moltype_for_size, Counts, &MoleculeType)) {
    exit(1);
  }

  // is '-m' used?
  bool m_option = false;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (specific_moltype_for_size[i] == 0) {
      m_option = true;
      break;
    }
  }
  // mass is only one; -m decides just what do we consider as A_S
  m_option = false; //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }

  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  }

  // count total number of chains in excluded aggs
  long int exclude_count_chains = 0;
  // count total number of excluded aggs
  long int exclude_count_agg = 0; //}}}

  // '-c' option - aggregate sizes //{{{
  int multiply = 1000;
  int composition[100] = {0}, comp_number_of_sizes = 0,
      types[2][2] = {{-1},{-1}}; // [x][0]: mol type; [x][1]: number of mols
  char output_comp[1024];
  if (FileIntsOption(argc, argv, "-c", composition, &comp_number_of_sizes, output_comp)) {
    exit(1);
  }
  if (comp_number_of_sizes > 0) {
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (specific_moltype_for_size[i] && types[0][0] == -1) {
        types[0][0] = i;
      } else if (specific_moltype_for_size[i] && types[0][0] != -1 && types[1][0] == -1) {
        types[1][0] = i;
      } else if (specific_moltype_for_size[i]) {
        fprintf(stderr, "\nError: '-c' option - more than two molecule types for composition distribution\n\n");
        exit(1);
      }
    }
    if (types[0][0] == -1 || types[1][0] == -1) {
      fprintf(stderr, "Error: '-c' option - less than two molecule types for composition distribution\n");
      exit(1);
    }
    for (int i = 0; i < comp_number_of_sizes; i++) {
      composition[i]+=0;
    }
  } //}}}

  // '--only' option //{{{
  int *only_specific_moltype_aggregates = calloc(Counts.TypesOfMolecules,sizeof(int));
  if (MoleculeTypeOption2(argc, argv, "--only", &only_specific_moltype_aggregates, Counts, &MoleculeType)) {
    exit(1);
  }

  // is '--only' used?
  bool only = false;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (only_specific_moltype_aggregates[i] == 1) {
      only = true;
      break;
    }
  }

  // count total number chains in specified aggs
  long int only_count_chains = 0;
  // count total number of specified aggs
  long int only_count_agg = 0; //}}}

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));

  for (int i = 0; i < Counts.Molecules; i++) {
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    fprintf(stdout, "Since no coordinates are used, no structure information is available and therefore the data is for the whole simulated system!\n\n");
    char null[1] = {'\0'};
    VerboseOutput(verbose2, null, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // arrays for distribution //{{{
  // number distribution
  long double ndistr[Counts.Molecules];
  // weight and z distributions - [][0] = mass of mols according to options; [][1] = mass of whole agg
  long double wdistr[Counts.Molecules][2];
  long double zdistr[Counts.Molecules][2];
  // volume distribution - probably works, but not really needed, so not sure
  // [][0] = volume of mols according to options; [][1] = volume of whole agg
  long double voldistr[Counts.Molecules][2];
  // number distribution of agg composition
  long double *ndistr_comp[Counts.Molecules];
  for (int i = 0; i < Counts.Molecules; i++) {
    ndistr_comp[i] = malloc((multiply*Counts.Molecules)*sizeof(long int));
  }
  // number of aggregates throughout simulation
  int count_agg[Counts.Molecules];
  // molecule typs in aggregates: [agg size][mol type][number or SQR(number)]
  int **molecules_sum = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    molecules_sum[i] = calloc(Counts.TypesOfMolecules,sizeof(int));
  }

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    ndistr[i] = 0;
    wdistr[i][0] = 0;
    wdistr[i][1] = 0;
    zdistr[i][0] = 0;
    zdistr[i][1] = 0;
    voldistr[i][0] = 0;
    voldistr[i][1] = 0;
    count_agg[i] = 0;
  } //}}}

  // print the first two lines to output file with per-step averages //{{{
  FILE *out;
  if ((out = fopen(output_avg, "w")) == NULL) {
    ErrorFileOpen(output_avg, 'w');
    exit(1);
  }

  // print command
  putc('#', out);
  for (int i = 0; i < argc; i++){
    fprintf(out, " %s", argv[i]);
  }
  putc('\n', out);

  count = 1;
  fprintf(out, "# column: (%d) step, ", count++);
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out, "(%d) <M>_n, ", count++);
    count++;
    fprintf(out, "(%d) <M>_w, ", count++);
    count++;
    fprintf(out, "(%d) <M>_z, ", count++);
    count++;
    fprintf(out, "(%d) <As>_n, ", count++);
    fprintf(out, "(%d) <As>_w, ", count++);
    count++;
    fprintf(out, "(%d) <As>_z, ", count++);
    count++;
  } else {
    fprintf(out, "(%d) <M>_n, ", count++);
    fprintf(out, "(%d) <M>_w, ", count++);
    fprintf(out, "(%d) <M>_z, ", count++);
    fprintf(out, "(%d) <As>_n, ", count++);
    fprintf(out, "(%d) <As>_w, ", count++);
    fprintf(out, "(%d) <As>_z, ", count++);
  }
  fprintf(out, "(%d) number of aggregates", count++);
  putc('\n', out);
  fclose(out); //}}}

  // main loop //{{{
  int count_step = 0;
  // [0][] = simple sum, [1][] = sum of squares, [2][] = sume of cubes
  // [][0] = mass of mols in agg from options, [][1] = mass of the whole aggregate
  double mass_sum[3][2] = {{0}}, As_sum[3][2] = {{0}};
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    // print (or not) step number //{{{
    count_step++;
    if (count_step < start) {
      if (!silent) {
        if (script) {
          fprintf(stdout, "Discarding Step: %6d\n", count_step);
        } else {
          fflush(stdout);
          fprintf(stdout, "\rDiscarding Step: %6d", count_step);
        }
      }
    } else if (count_step == start) {
      if (!silent) {
        if (script) {
          fprintf(stdout, "Discarded Steps: %6d\n", count_step-1);
        } else {
          fflush(stdout);
          fprintf(stdout, "\rDiscarded Steps: %6d\n", count_step-1);
        }
      }
    } else {
      if (!silent) {
        if (script) {
          fprintf(stdout, "Step: %6d\n", count_step);
        } else {
          fflush(stdout);
          fprintf(stdout, "\rStep: %6d", count_step);
        }
      }
    } //}}}

    ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index);

    // print info about aggregates if '-V' is used //{{{
    if (verbose2) {
      for (int i = 0; i < Counts.Aggregates; i++) {
        fprintf(stdout, "\nAggregate[%3d].{Mass = %6.2f,\nnMolecules = %3d:", i+1, Aggregate[i].Mass, Aggregate[i].nMolecules);
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          fprintf(stdout, " %d", Aggregate[i].Molecule[j]+1);
        }
        fprintf(stdout, ",\n nBeads = %4d:", Aggregate[i].nBeads);
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          fprintf(stdout, " %d", Aggregate[i].Bead[j]);
        }
        fprintf(stdout, ",\n nMonomers = %4d:", Aggregate[i].nMonomers);
        for (int j = 0; j < Aggregate[i].nMonomers; j++) {
          fprintf(stdout, " %d", Aggregate[i].Monomer[j]);
        }
        fprintf(stdout, "}\n");
      }
      putchar('\n');
    } //}}}

    // go through all aggregates
    int aggs_step = 0; // number of eligible aggregates per step
    double avg_mass_n_step[2] = {0}, avg_mass_w_step[2] = {0}, avg_mass_z_step[2] = {0}, // per-step mass averages
           avg_As_n_step[2] = {0}, avg_As_w_step[2] = {0},  avg_As_z_step[2] = {0}; // per-step As averages
    for (int i = 0; i < Counts.Aggregates; i++) {

      /* determine aggregate size and mass: //{{{
       * if `-m` is used, the size is the number of 'specific_moltype_for_size'
       * mols and the mass is the mass only of those molecules */
      int size = 0;
      double agg_mass = 0, agg_vol = 0; // mass and volume according to options
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type]) {
          size++;
          agg_mass += MoleculeType[mol_type].Mass;
          agg_vol += MoleculeType[mol_type].nBeads; // all DPD beads have same volume
        }
      }
      // make calculations only if agg size is well defined and within given range
      if (size == 0 || size < range_As[0] || size > range_As[1]) {
        continue;
      } //}}}

      // if '--only' is used, use only aggregates composed the specified molecule(s) //{{{
      bool test = true;
      if (only) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (only_specific_moltype_aggregates[Molecule[id].Type] == 0) {
            test = false; // a molecule is not of the required type
            break;
          }
        }
        if (!test) { // should the rest of the for loop agg i be skipped?
          continue;
        }
      } //}}}

      // if '-x' option is used, discount aggregates with only specified molecule type(s) //{{{
      test = false;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
        if (MoleculeType[moltype].Write) {
          test = true; // a molecule that shouldn't be in agg 'i' is there
          break;
        }
      }
      if (!test) { // should the rest of the for loop agg i be skipped?
        exclude_count_chains += Aggregate[i].nMolecules;
        exclude_count_agg++;
        continue;
      } //}}}

      // if '-c' is used, calculate composition distribution //{{{
      if (comp_number_of_sizes > 0) {
        // use the size?
        bool use_size_comp = false;
        for (int j = 0; j < comp_number_of_sizes; j++) {
          if (size == composition[j]) {
            use_size_comp = true;
            break;
          }
        }

        if (use_size_comp) {
          // count numbers of the two mol types in aggregate 'i'
          types[0][1] = 0;
          types[1][1] = 0;
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            // molecule type
            int id = Molecule[Aggregate[i].Molecule[j]].Type;
            if (types[0][0] == id) {
              types[0][1]++;
            } else if (types[1][0] == id) {
              types[1][1]++;
            }
          }
          double ratio;
          if (types[1][1] == 0) {
            ratio = 1;
          } else {
            ratio = (double)(types[0][1]) / types[1][1];
          }
          ndistr_comp[size-1][(int)(ratio*multiply)]++;
        }
      } //}}}

      // for average aggregate mass during the step
      avg_mass_n_step[0] += agg_mass;
      avg_mass_w_step[0] += SQR(agg_mass);
      avg_mass_z_step[0] += CUBE(agg_mass);
      avg_mass_n_step[1] += Aggregate[i].Mass;
      avg_mass_w_step[1] += SQR(Aggregate[i].Mass);
      avg_mass_z_step[1] += CUBE(Aggregate[i].Mass);

      // for average aggregation number during the step
      avg_As_n_step[0] += size;
      avg_As_w_step[0] += size * agg_mass;
      avg_As_z_step[0] += size * SQR(agg_mass);
      avg_As_n_step[1] += Aggregate[i].nMolecules;
      avg_As_w_step[1] += Aggregate[i].nMolecules * Aggregate[i].Mass;
      avg_As_z_step[1] += Aggregate[i].nMolecules * SQR(Aggregate[i].Mass);

      aggs_step++;

      // start calculation of averages from specified 'start' timestep
      if (count_step >= start) {

        // distribution //{{{
        ndistr[size-1]++;
        wdistr[size-1][0] += agg_mass;
        wdistr[size-1][1] += Aggregate[i].Mass;
        zdistr[size-1][0] += SQR(agg_mass);
        zdistr[size-1][1] += SQR(Aggregate[i].Mass);
        voldistr[size-1][0] += agg_vol;
        voldistr[size-1][1] += Aggregate[i].nBeads; //}}}

        // number of various species in the aggregate //{{{
        count_agg[size-1]++;

        As_sum[0][0] += size;
        As_sum[1][0] += size * agg_mass;
        As_sum[2][0] += size * SQR(agg_mass);
        As_sum[0][1] += Aggregate[i].nMolecules;
        As_sum[1][1] += Aggregate[i].nMolecules * Aggregate[i].Mass;
        As_sum[2][1] += Aggregate[i].nMolecules * SQR(Aggregate[i].Mass);

        mass_sum[0][0] += agg_mass;
        mass_sum[1][0] += SQR(agg_mass);
        mass_sum[2][0] += CUBE(agg_mass);
        mass_sum[0][1] += Aggregate[i].Mass;
        mass_sum[1][1] += SQR(Aggregate[i].Mass);
        mass_sum[2][1] += CUBE(Aggregate[i].Mass);

        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
          molecules_sum[size-1][mol_type]++;
        } //}}}

        // count chains if `--only` is used //{{{
        if (only) {
          only_count_chains += Aggregate[i].nMolecules; // superfluous
          only_count_agg++;
        } //}}}
      }
    }

    // print averages to output file (only if non-zero eligible aggregates) //{{{
    if (aggs_step > 0) {
      if ((out = fopen(output_avg, "a")) == NULL) {
        ErrorFileOpen(output_avg, 'a');
        exit(1);
      }

      fprintf(out, "%5d", count_step); // step
      fprintf(out, " %10.5f", avg_mass_n_step[0]/aggs_step); // <mass>_n
      if (m_option) {
        fprintf(out, " %10.5f", avg_mass_n_step[1]/aggs_step); // <mass>_n (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_mass_w_step[0]/avg_mass_n_step[0]); // <mass>_w
      if (m_option) {
        fprintf(out, " %10.5f", avg_mass_w_step[1]/avg_mass_n_step[1]); // <mass>_w (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_mass_z_step[0]/avg_mass_w_step[0]); // <mass>_z
      if (m_option) {
        fprintf(out, " %10.5f", avg_mass_z_step[1]/avg_mass_w_step[1]); // <mass>_z (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_As_n_step[0]/aggs_step); // <As>_n
      if (m_option) {
        fprintf(out, " %10.5f", avg_As_n_step[1]/aggs_step); // <As>_n (whole agg)
      }
      fprintf(out, " %10.5f", avg_As_w_step[0]/avg_mass_n_step[0]); // <As>_w
      if (m_option) {
        fprintf(out, " %10.5f", avg_As_w_step[1]/avg_mass_n_step[1]); // <As>_w (whole agg mass)
      }
      fprintf(out, " %10.5f", avg_As_z_step[0]/avg_mass_w_step[0]); // <As>_z
      if (m_option) {
        fprintf(out, " %10.5f", avg_As_z_step[1]/avg_mass_w_step[1]); // <As>_z (whole agg mass)
      }
      fprintf(out, " %5d", aggs_step); // number of aggregates in the step

      putc('\n', out);
      fclose(out); //}}}
    }
  }
  fclose(agg);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count_step);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count_step);
    }
  } //}}}

  // print the first two lines to output file with distributions //{{{
  if ((out = fopen(output_distr, "w")) == NULL) {
    ErrorFileOpen(output_distr, 'w');
    exit(1);
  }

  // print command
  putc('#', out);
  for (int i = 0; i < argc; i++){
    fprintf(out, " %s", argv[i]);
  }
  putc('\n', out);

  count = 1;
  fprintf(out, "# column: ");
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out, "(%d) As, ", count++);
    fprintf(out, "(%d) F_n(As), ", count++);
    fprintf(out, "(%d) F_w(As), ", count++);
    count++;
    fprintf(out, "(%d) F_z(As), ", count++);
    count++;
    fprintf(out, "(%d) <volume distribution>, ", count++);
    count++;
  } else {
    fprintf(out, "(%d) As, ", count++);
    fprintf(out, "(%d) F_n(As), ", count++);
    fprintf(out, "(%d) F_w(As), ", count++);
    fprintf(out, "(%d) F_z(As), ", count++);
    fprintf(out, "(%d) <volume distribution>, ", count++);
  }
  fprintf(out, "(%d) number of aggs,", count++);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " (%d) <%s>_n", i+count, MoleculeType[i].Name);
    if (i != (Counts.TypesOfMolecules-1)) {
      putc(',', out);
    }
  }
  putc('\n', out);
  fclose(out); //}}}

  // print distributions to output file //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    ErrorFileOpen(output_distr, 'a');
    exit(1);
  }

  count_step -= start - 1;

  // normalization factors
  long int ndistr_norm = 0, wdistr_norm[2] = {0}, zdistr_norm[2] = {0}, voldistr_norm[2] = {0};
  for (int i = 0; i < Counts.Molecules; i++) {
    ndistr_norm += ndistr[i];

    wdistr_norm[0] += wdistr[i][0];
    wdistr_norm[1] += wdistr[i][1];

    zdistr_norm[0] += zdistr[i][0];
    zdistr_norm[1] += zdistr[i][1];

    voldistr_norm[0] += voldistr[i][0];
    voldistr_norm[1] += voldistr[i][1];
  }

  for (int i = 0; i < Counts.Molecules; i++) {
    if (count_agg[i] > 0) {
      fprintf(out, "%4d", i+1); // As
      fprintf(out, " %10.5f", (double)(ndistr[i])/ndistr_norm); // number As distr
      fprintf(out, " %10.5f", (double)(wdistr[i][0])/wdistr_norm[0]); // weight distr
      if (m_option) {
        fprintf(out, " %10.5f", (double)(wdistr[i][1])/wdistr_norm[1]); // weight distr (whole agg mass)
      }
      fprintf(out, " %10.5f", (double)(zdistr[i][0])/zdistr_norm[0]); // z distr
      if (m_option) {
        fprintf(out, " %10.5f", (double)(zdistr[i][1])/zdistr_norm[1]); // z distr (whole agg mass)
      }
      fprintf(out, " %10.5f", (double)(voldistr[i][0])/voldistr_norm[0]); // volume distribution
      if (m_option) {
        fprintf(out, " %10.5f", (double)(voldistr[i][1])/voldistr_norm[1]); // volume distribution (whole agg mass)
      }
      fprintf(out, " %6d", count_agg[i]); // number of aggregates
      // print average number of molecule types in aggregates
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        if (count_agg[i] == 0) {
          fprintf(out, "%8s", "?");
        } else {
          fprintf(out, " %10.5f", (double)(molecules_sum[i][j])/count_agg[i]);
        }
      }
      putc('\n', out);
    }
  }
  fclose(out); //}}}

  // print overall averages to avg and distr output files //{{{
  for (int i = 1; i < Counts.Molecules; i++) {
    count_agg[0] += count_agg[i]; // total number of aggregates
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules_sum[0][j] += molecules_sum[i][j]; // total number molecules of each type
    }
  }

  // open output files //{{{
  // distr file
  if ((out = fopen(output_distr, "a")) == NULL) {
    ErrorFileOpen(output_distr, 'a');
    exit(1);
  }
  // avg file
  FILE *out2;
  if ((out2 = fopen(output_avg, "a")) == NULL) {
    ErrorFileOpen(output_avg, 'a');
    exit(1);
  } //}}}

  // print legend (with column numbers) //{{{
  putc('#', out);
  putc('#', out2);
  // distr file
  count = 1;
  fprintf(out, " (%d) <As>_n,", count++);
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out, " (%d) <As>_w,", count++);
    count++;
    fprintf(out, " (%d) <As>_z,", count++);
    count++;
    fprintf(out, " (%d) <M>_n,", count++);
    count++;
    fprintf(out, " (%d) <M>_w,", count++);
    count++;
    fprintf(out, " (%d) <M>_z,", count++);
    count++;
  } else {
    fprintf(out, " (%d) <As>_w,", count++);
    fprintf(out, " (%d) <As>_z,", count++);
    fprintf(out, " (%d) <M>_n,", count++);
    fprintf(out, " (%d) <M>_w,", count++);
    fprintf(out, " (%d) <M>_z,", count++);
  }
//count++;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " (%d) <%s>_n", i+count, MoleculeType[i].Name);
    if (i != (Counts.TypesOfMolecules-1) && !only) {
      putc(',', out);
    }
  }
  // avg file
  count = 1;
  fprintf(out2, " (%d) <As>_n,", count++);
  if (m_option) { // -m is used, two agg masses are considered
    fprintf(out2, " (%d) <As>_w,", count++);
    count++;
    fprintf(out2, " (%d) <As>_z,", count++);
    count++;
    fprintf(out2, " (%d) <M>_n,", count++);
    count++;
    fprintf(out2, " (%d) <M>_w,", count++);
    count++;
    fprintf(out2, " (%d) <M>_z,", count++);
    count++;
  } else {
    fprintf(out2, " (%d) <As>_w,", count++);
    fprintf(out2, " (%d) <As>_z,", count++);
    fprintf(out2, " (%d) <M>_n,", count++);
    fprintf(out2, " (%d) <M>_w,", count++);
    fprintf(out2, " (%d) <M>_z,", count++);
  }
//count++;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out2, " (%d) <%s>_n", i+count, MoleculeType[i].Name);
    if (i != (Counts.TypesOfMolecules-1) && !only) {
      putc(',', out2);
    }
  }
//if (only) {
//  // distr file
//  fprintf(out, " (%d) <fraction of '--only' mols not touching other mols>", Counts.TypesOfMolecules+count);
//  // avg file
//  fprintf(out2, " (%d) <fraction of '--only' mols not touching other mols>", Counts.TypesOfMolecules+count);
//}
  putc('\n', out);
  putc('\n', out2); //}}}

  // print the averages //{{{
  // distr file
  fprintf(out, "#");
  fprintf(out2, "#");

  // distr file
  fprintf(out, " %lf", As_sum[0][0]/count_agg[0]); // <As>_n
  fprintf(out, " %lf", As_sum[1][0]/mass_sum[0][0]); // <As>_w
  if (m_option) {
    fprintf(out, " %lf", As_sum[1][1]/mass_sum[0][1]); // <As>_w (whole agg)
  }
  fprintf(out, " %lf", As_sum[2][0]/mass_sum[1][0]); // <As>_z
  if (m_option) {
    fprintf(out, " %lf", As_sum[2][1]/mass_sum[1][1]); // <As>_z (whole agg)
  }

  fprintf(out, " %lf", mass_sum[0][0]/count_agg[0]); // <M>_n
  if (m_option) {
    fprintf(out, " %lf", mass_sum[0][1]/count_agg[0]); // <M>_n (whole agg)
  }
  fprintf(out, " %lf", mass_sum[1][0]/mass_sum[0][0]); // <M>_w
  if (m_option) {
    fprintf(out, " %lf", mass_sum[1][1]/mass_sum[0][1]); // <M>_w (whole agg)
  }
  fprintf(out, " %lf", mass_sum[2][0]/mass_sum[1][0]); // <M>_z
  if (m_option) {
    fprintf(out, " %lf", mass_sum[2][1]/mass_sum[1][1]); // <M>_z (whole agg)
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %lf", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
  }
  // avg file
  fprintf(out2, " %lf", As_sum[0][0]/count_agg[0]); // <As>_n
  fprintf(out2, " %lf", As_sum[1][0]/mass_sum[0][0]); // <As>_w
  if (m_option) {
    fprintf(out2, " %lf", As_sum[1][1]/mass_sum[0][1]); // <As>_w (whole agg)
  }
  fprintf(out2, " %lf", As_sum[2][0]/mass_sum[1][0]); // <As>_z
  if (m_option) {
    fprintf(out2, " %lf", As_sum[2][1]/mass_sum[1][1]); // <As>_z (whole agg)
  }

  fprintf(out2, " %lf", mass_sum[0][0]/count_agg[0]); // <M>_n
  if (m_option) {
    fprintf(out2, " %lf", mass_sum[0][1]/count_agg[0]); // <M>_n (whole agg)
  }
  fprintf(out2, " %lf", mass_sum[1][0]/mass_sum[0][0]); // <M>_w
  if (m_option) {
    fprintf(out2, " %lf", mass_sum[1][1]/mass_sum[0][1]); // <M>_w (whole agg)
  }
  fprintf(out2, " %lf", mass_sum[2][0]/mass_sum[1][0]); // <M>_z
  if (m_option) {
    fprintf(out2, " %lf", mass_sum[2][1]/mass_sum[1][1]); // <M>_z (whole agg)
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out2, " %lf", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
  }

//if (only) {
//  int only_molecules_sum = 0;
//  int only_moltype_number = 0;
//  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//    if (only_specific_moltype_aggregates[i]) {
//      only_molecules_sum += molecules_sum[0][i];
//      only_moltype_number += MoleculeType[i].Number;
//    }
//  }
//  // distr file
//  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//    if (only_specific_moltype_aggregates[i]) {
//      fprintf(out, "%10.3f", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
//    }
//  }
//  fprintf(out, " %10.3f", (double)(only_count_chains)/(count_step*only_moltype_number));
//  // avg file
//  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//    if (only_specific_moltype_aggregates[i]) {
//      fprintf(out2, "%10.3f", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
//    }
//  }
//  fprintf(out2, " %10.3f", (double)(only_count_chains)/(count_step*only_moltype_number));
//}

  putc('\n', out);
  putc('\n', out2); //}}}

  fclose(out);
  fclose(out2);
  //}}}

  // print composition distribution //{{{
  if (comp_number_of_sizes > 0) {
    // open file with composition distribution //{{{
    if ((out = fopen(output_comp, "w")) == NULL) {
      ErrorFileOpen(output_comp, 'w');
      exit(1);
    } //}}}

    // print command //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++){
      fprintf(out, " %s", argv[i]);
    }
    putc('\n', out); //}}}

    // print header line //{{{
    fprintf(out, "# column: (1) ");
    count = 0;
    // molecule names
    fprintf(out, "%s/%s; agg sizes - ", MoleculeType[types[0][0]].Name, MoleculeType[types[1][0]].Name);
    // aggregate sizes
    for (int i = 0; i < comp_number_of_sizes; i++) {
      fprintf(out, " (%d) %d;", i+2, composition[i]);
    }
    putc('\n', out); //}}}

    // normalization factor (simple sum)
    long int *sum = calloc(Counts.Molecules,sizeof(long int));
    int low = multiply * Counts.Molecules;
    bool *include = calloc((multiply*Counts.Molecules),sizeof(bool));
    for (int i = 0; i < (multiply*Counts.Molecules); i++) {
      for (int j = 0; j < comp_number_of_sizes; j++) {
        if (ndistr_comp[composition[j]-1][i] > 0 && low == (multiply*Counts.Molecules)) {
          low = i;
        }
        if (ndistr_comp[composition[j]-1][i] > 0) {
          include[i] = true;
        }
        sum[composition[j]-1] += ndistr_comp[composition[j]-1][i];
      }
    }

    // print data
    for (int i = low; i < (multiply*Counts.Molecules); i++) {
      if (include[i]) {
        fprintf(out, "%10.5f", (double)(i)/multiply);
        for (int j = 0; j < comp_number_of_sizes; j++) {
          if (ndistr_comp[composition[j]-1][i] != 0) {
            fprintf(out, " %7.5f", (double)(ndistr_comp[composition[j]-1][i])/sum[composition[j]-1]);
          } else {
            fprintf(out, "       ?");
          }
        }
        putc('\n', out);
      }
    }
    fclose(out);

    free(sum);
    free(include);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(specific_moltype_for_size);
  free(only_specific_moltype_aggregates);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(molecules_sum[i]);
    free(ndistr_comp[i]);
  }
  free(molecules_sum); //}}}

  return 0;
}
