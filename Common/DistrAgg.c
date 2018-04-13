#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input> <output distr file> <output avg file> <options>\n\n", cmd);

  fprintf(stderr, "   <input>           input filename (agg format)\n");
  fprintf(stderr, "   <distr file>      filename with weight and number distributions\n");
  fprintf(stderr, "   <avg file>        filename with weight and number averages throughout simulation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -st <int>      start distribution calculation with <int>-th step\n");
  fprintf(stderr, "      -n <int> <int> calculate for aggregate sizes in given range\n");
  fprintf(stderr, "      -m <name>      agg size means number of <name> molecule types in an aggregate\n");
  fprintf(stderr, "      -x <name(s)>   exclude aggregates containing only specified molecule(s)\n");
  fprintf(stderr, "      --only <name>  use just aggregates composed of a specified molecule\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
DistrAgg calculates weight and number average aggregation numbers during the \
simulation run as well as overall weight and number distributions and volume \
fractions of aggregates.\n\n");

      fprintf(stdout, "\
The utility uses dl_meso.vsf (or other input structure file) and FIELD \
(along with optional bond file) files to determine all information about \
the system.\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input> <output distr file> <output avg file> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input>           input filename (agg format)\n");
      fprintf(stdout, "   <distr file>      filename with weight and number distributions\n");
      fprintf(stdout, "   <avg file>        filename with weight and number averages throughout simulation\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -st <int>      start distribution calculation with <int>-th step\n");
      fprintf(stdout, "      -n <int> <int> calculate for aggregate sizes in given range\n");
      fprintf(stdout, "      -m <name>      agg size means number of <name> molecule types in an aggregate\n");
      fprintf(stdout, "      -x <name(s)>   exclude aggregates containing only specified molecule(s)\n");
      fprintf(stdout, "      --only <name>  use just aggregates composed of a specified molecule\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 3; //}}}

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
        strcmp(argv[i], "--only") != 0 ) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < options) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than dl_meso.vsf? //{{{
  char *vsf_file = calloc(32,sizeof(char *));
  if (VsfFileOption(argc, argv, &vsf_file)) {
    exit(1);
  } //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(32,sizeof(char *));
  if (BondsFileOption(argc, argv, &bonds_file)) {
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

  // <input> - filename of input agg file //{{{
  char input_agg[32];
  strcpy(input_agg, argv[++count]); //}}}

  // open input file and skip the first two lines //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // <output distr file> - filename with weight and number distributions //{{{
  char output_distr[123];
  strcpy(output_distr, argv[++count]); //}}}

  // <output avg file> - filename with weight and number average aggregation numbers //{{{
  char output_avg[123];
  strcpy(output_avg, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information //{{{
  char vcf[1];
  vcf[0] = '\0';
  ReadStructure(vsf_file, vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file); //}}}

  // '-n' option - range of aggregation numbers //{{{
  int range_As[2];
  range_As[0] = 1;
  range_As[1] = Counts.Molecules;
  if (TwoIntegerOption(argc, argv, "-n", range_As)) {
    exit(1);
  }
  range_As[0]--;

  // make sure first number is larger
  if (range_As[0] > range_As[1]) {
    int tmp = range_As[0];
    range_As[0] = range_As[1];
    range_As[1] = tmp;
  } //}}}

  // '-m' option //{{{
  int *specific_moltype_for_size;
  specific_moltype_for_size = malloc(Counts.TypesOfMolecules*sizeof(int *));
  // all are to be used without '-m' option
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", &specific_moltype_for_size, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }

  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  }
  //}}}

  // '--only' option //{{{
  int only_specific_moltype_aggregates = -1;
  if (MoleculeTypeOption(argc, argv, "--only", &only_specific_moltype_aggregates, Counts, &MoleculeType)) {
    exit(1);
  }

  // to calculate fraction of chains not in contact with other moltypes
  long int only_count_chains = 0; //}}}

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
  long int ndistr[Counts.Molecules];
  // weight and z distributions - [][0] = mass of mols according to options; [][1] = mass of whole agg
  long int wdistr[Counts.Molecules][2];
  long int zdistr[Counts.Molecules][2];
  // volume distribution - probably works, but not really needed, so not sure
  // [][0] = volume of mols according to options; [][1] = volume of whole agg
  long int voldistr[Counts.Molecules][2];
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
    fprintf(stderr, "Cannot open file %s!\n", output_avg);
    exit(1);
  }

  // print command
  putc('#', out);
  for (int i = 0; i < argc; i++){
    fprintf(out, " %s", argv[i]);
  }
  putc('\n', out);

  fprintf(out, "# 1:step ");
  fprintf(out, "2:<M>_n (options' mass) ");
  fprintf(out, "3:<M>_n (whole agg mass) ");
  fprintf(out, "4:<M>_w (options' mass) ");
  fprintf(out, "5:<M>_w (whole agg mass)\n");
  fclose(out); //}}}

  // main loop //{{{
  int test; // sum_mass = 0;
  count = 0;
  // [0][] = simple sum, [1][] = sum of squares, [2][] = sume of cubes
  // [][0] = mass of mols in agg from options, [][1] = mass of the whole aggregate
  double mass_sum[3][2] = {0};
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    // print (or not) step number //{{{
    count++;
    if (count < start) {
      if (!silent) {
        if (script) {
          fprintf(stdout, "Discarding Step: %6d\n", count);
        } else {
          fflush(stdout);
          fprintf(stdout, "\rDiscarding Step: %6d", count);
        }
      }
    } else if (count == start) {
      if (!silent) {
        if (script) {
          fprintf(stdout, "Discarded Steps: %6d\n", count-1);
        } else {
          fflush(stdout);
          fprintf(stdout, "\rDiscarded Steps: %6d\n", count-1);
        }
      }
    } else {
      if (!silent) {
        if (script) {
          fprintf(stdout, "Step: %6d\n", count);
        } else {
          fflush(stdout);
          fprintf(stdout, "\rStep: %6d", count);
        }
      }
    } //}}}

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

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
    int aggs_step = 0, // number of aggregates (w/o unimers if --no-unimers)
        avg_mass_n_step[2] = {0}, avg_mass_w_step[2] = {0}, // per-step mass averages
        avg_As_n_step[2] = {0}, avg_As_w_step[2] = {0}; // per-step As averages
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
      if (size == 0 || size <= range_As[0] || size >= range_As[1]) {
        continue;
      } //}}}

      // for average aggregate mass during the step //{{{
      avg_mass_n_step[0] += agg_mass;
      avg_mass_w_step[0] += SQR(agg_mass);
      avg_mass_n_step[1] += Aggregate[i].Mass;
      avg_mass_w_step[1] += SQR(Aggregate[i].Mass); //}}}

      // for average aggregation number during the step //{{{
      avg_As_n_step[0] += size;
      avg_As_w_step[0] += SQR(size);
      avg_As_n_step[1] += Aggregate[i].nMolecules;
      avg_As_w_step[1] += SQR(Aggregate[i].nMolecules); //}}}

      // if '--only' is used, use only aggregates composed the specified molecule //{{{
      bool test = true;
      if (only_specific_moltype_aggregates != -1) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (only_specific_moltype_aggregates != Molecule[id].Type) {
            test = false; // a molecule is not of the required type
            break;
          }
        }
        if (!test) { // should the rest of the for loop agg i be skipped?
          continue;
        }
      } //}}}

      // if '-x' option is used, discount aggregates with only specified molecule types //{{{
      test = false;
      for (int j = 0; j < size; j++) {
        int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
        if (MoleculeType[moltype].Write) {
          test = true; // a molecule that shouldn't be in agg 'i' is there
          break;
        }
      }
      if (!test) { // should the rest of the for loop agg i be skipped?
        continue;
      } //}}}

      // number of eligible aggregate for a step
      aggs_step++;

      // start calculation of averages from specified 'start' timestep
      if (count >= start) {

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

        mass_sum[0][0] += agg_mass;
        mass_sum[0][1] += Aggregate[i].Mass;
        mass_sum[1][0] += SQR(agg_mass);
        mass_sum[1][1] += SQR(Aggregate[i].Mass);
        mass_sum[2][0] += CUBE(agg_mass);
        mass_sum[2][1] += CUBE(Aggregate[i].Mass);

        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
          molecules_sum[size-1][mol_type]++;
        } //}}}

        // count chains if `--only` is used //{{{
        if (only_specific_moltype_aggregates != -1) {
          only_count_chains += Aggregate[i].nMolecules;
        } //}}}
      }
    }

    // print averages to output file //{{{
    // REDO CORRECTLY
    if ((out = fopen(output_avg, "a")) == NULL) {
      // print newline to stdout if Step... doesn't end with one
      if (!script && !silent) {
        putchar('\n');
      }
      fprintf(stderr, "Cannot open file %s!\n", output_avg);
      exit(1);
    }

    fprintf(out, "%5d", count); // step
    fprintf(out, " %8.3f", (double)(avg_mass_n_step[0])/aggs_step); // <M>_n (options' mass)
    fprintf(out, " %8.3f", (double)(avg_mass_n_step[1])/aggs_step); // <M>_n (whole agg mass)
    fprintf(out, " %8.3f", (double)(avg_mass_w_step[0])/avg_mass_n_step[0]); // <M>_w (options' mass)
    fprintf(out, " %8.3f", (double)(avg_mass_w_step[1])/avg_mass_n_step[1]); // <M>_w (whole agg mass)
    fprintf(out, " %8.3f", (double)(avg_As_n_step[0])/aggs_step); // <A_s>_n (options' mass)
    fprintf(out, " %8.3f", (double)(avg_As_n_step[1])/aggs_step); // <A_s>_n (whole agg mass)
    fprintf(out, " %8.3f", (double)(avg_As_w_step[0])/avg_As_n_step[0]); // <A_s>_w (options' mass)
    fprintf(out, " %8.3f", (double)(avg_As_w_step[1])/avg_As_n_step[1]); // <A_s>_w (whole agg mass)

    putc('\n', out);
    fclose(out); //}}}
  }
  fclose(agg);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  } //}}}

  // print the first two lines to output file with distributions //{{{
  if ((out = fopen(output_distr, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  // print command
  putc('#', out);
  for (int i = 0; i < argc; i++){
    fprintf(out, " %s", argv[i]);
  }
  putc('\n', out);

  fprintf(out, "# ");
  fprintf(out, "1:A_s ");
  fprintf(out, "2:F_n(As) ");
  fprintf(out, "3:F_w(As) (options' mass) ");
  fprintf(out, "4:F_w(As) (whole agg mass) ");
  fprintf(out, "5:F_z(As) (options' mass) ");
  fprintf(out, "6:F_z(As) (whole agg mass) ");
  fprintf(out, "7:<volume distribution> (options' volume) ");
  fprintf(out, "8:<volume distribution> (whole agg volume) ");
  fprintf(out, "9:number of aggs ");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %d:<%s>_n", i+10, MoleculeType[i].Name);
  }
  putc('\n', out);
  fclose(out); //}}}

  // print distributions to output file //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  count -= start -1;

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
      fprintf(out, "%4d %lf %lf %lf %lf %lf %lf %lf %6d", i+1, // A_s
                  (double)(ndistr[i])/ndistr_norm, // number As distr
                  (double)(wdistr[i][0])/wdistr_norm[0], // weight distr (options mass)
                  (double)(wdistr[i][1])/wdistr_norm[1], // weight distr (whole mass)
                  (double)(zdistr[i][0])/zdistr_norm[0], // z distr (options mass)
                  (double)(zdistr[i][1])/zdistr_norm[1], // z distr (whole agg mass)
                  (double)(voldistr[i][0])/voldistr_norm[0], // volume distribution (options agg)
                  (double)(voldistr[i][1])/voldistr_norm[1], // volume distribution (whole agg)
                  count_agg[i]);
      // print average number of molecule types in aggregates
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        if (count_agg[i] == 0) {
          fprintf(out, "%8s", "?");
        } else {
          fprintf(out, " %7.3f", (double)(molecules_sum[i][j])/count_agg[i]);
        }
      }
      putc('\n', out);
    }
  }
  fclose(out); //}}}

  // print overall averages //{{{
  for (int i = 1; i < Counts.Molecules; i++) {
    count_agg[0] += count_agg[i]; // total number of aggregates (possibly without unimers)
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules_sum[0][j] += molecules_sum[i][j]; // total number of each molecular species
    }
  }

  // open file with distribution //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_avg);
    exit(1);
  } //}}}

  // print legend (with column numbers) //{{{
  fprintf(out, "# ");
  if (only_specific_moltype_aggregates != -1) {
    fprintf(out, "1:<M>_n ");
    fprintf(out, "2:<M>_w ");
    fprintf(out, "3:<M>_z ");
    fprintf(out, "4:<%s>_n ", MoleculeType[only_specific_moltype_aggregates].Name);
    fprintf(out, "5:<fraction of %s chains not in contact with other moltypes>", MoleculeType[only_specific_moltype_aggregates].Name);
  } else {
    fprintf(out, "1:<M>_n (options' mass) ");
    fprintf(out, "2:<M>_n (whole agg mass) ");
    fprintf(out, "3:<M>_w (options' mass) ");
    fprintf(out, "4:<M>_w (whole agg mass) ");
    fprintf(out, "5:<M>_z (options' mass) ");
    fprintf(out, "6:<M>_z (whole agg mass)");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, " %d:<%s>_n", i+7, MoleculeType[i].Name);
    }
  }
  putc('\n', out); //}}}

  // print the averages //{{{
  fprintf(out, "# ");
  if (only_specific_moltype_aggregates != -1) {
    // --only option doesn't distinguish between options's stuff and whole agg stuff
    fprintf(out, "%10.3f ", mass_sum[0][0]/count_agg[0]); // <M>_n
    fprintf(out, "%10.3f ", mass_sum[1][0]/mass_sum[0][0]); // <M>_w
    fprintf(out, "%10.3f ", mass_sum[2][0]/mass_sum[1][0]); // <M>_z
    fprintf(out, "%10.3f", (double)(molecules_sum[0][only_specific_moltype_aggregates])/count_agg[0]); // <species>_n
    fprintf(out, " %10.3f", (double)(only_count_chains)/(count*MoleculeType[only_specific_moltype_aggregates].Number));
  } else {
    fprintf(out, "%10.3f ", mass_sum[0][0]/count_agg[0]); // <M>_n (options' mass)
    fprintf(out, "%10.3f ", mass_sum[0][1]/count_agg[0]); // <M>_n (whole agg mass)
    fprintf(out, "%10.3f ", mass_sum[1][0]/mass_sum[0][0]); // <M>_w (options' mass)
    fprintf(out, "%10.3f ", mass_sum[1][1]/mass_sum[0][1]); // <M>_w (whole agg mass)
    fprintf(out, "%10.3f ", mass_sum[2][0]/mass_sum[1][0]); // <M>_z (options' mass)
    fprintf(out, "%10.3f ", mass_sum[2][1]/mass_sum[1][1]); // <M>_z (whole agg mass)
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, "%10.3f", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
    }
  }
  putc('\n', out);
  fclose(out); //}}}
  //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(molecules_sum[i]);
  }
  free(molecules_sum); //}}}

  return 0;
}
