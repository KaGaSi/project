#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
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
  char output_distr[32];
  strcpy(output_distr, argv[++count]); //}}}

  // <output avg file> - filename with weight and number average aggregation numbers //{{{
  char output_avg[32];
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
  }

//// testing output
//for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//  printf("%d %s\n", specific_moltype_for_size[i], MoleculeType[i].Name);
//} //}}}

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
  } //}}}

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
  // distributions of agg mass - [][0] = mass of mols according to options; [][1] = mass of whole agg
  long int ndistr_mass[Counts.Molecules][2];
  long int wdistr_mass[Counts.Molecules][2];
  // number distribution of As
  long int ndistr_As[Counts.Molecules];
  // weight distribution of As - [][0] = As as number of mols according to options; [][1] = all mols in aggregate
  long int wdistr_As[Counts.Molecules][2];
  // volume distribution - probably works, but not really needed, so not sure
  long int voldistr[Counts.Molecules];
  // number of aggregates throughout simulation
  int count_agg[Counts.Molecules];
  // molecule typs in aggregates: [agg size][mol type][number or SQR(number)]
  int **molecules_sum = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    molecules_sum[i] = calloc(Counts.TypesOfMolecules,sizeof(int));
  }

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    ndistr_mass[i][0] = 0;
    wdistr_mass[i][0] = 0;
    ndistr_mass[i][1] = 0;
    wdistr_mass[i][1] = 0;

    ndistr_As[i] = 0;

    wdistr_As[i][0] = 0;
    wdistr_As[i][1] = 0;

    voldistr[i] = 0;
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
  // add date
  char date[26];
  Date(date);
  fprintf(out, " on %s\n", date);

  fprintf(out, "# 1:step ");
  fprintf(out, "2:<M>_n (options' mass) ");
  fprintf(out, "3:<M>_w (options' mass) ");
  fprintf(out, "4:<M>_n (whole agg mass) ");
  fprintf(out, "5:<M>_w (whole agg mass)\n");
  fprintf(out, "6:<As>_n (options' As) ");
  fprintf(out, "7:<As>_w (options' As) ");
  fprintf(out, "8:<As>_w (whole agg As)\n");
  fclose(out); //}}}

  // main loop //{{{
  int test; // sum_mass = 0;
  count = 0;
  // [0][] = simple sum, [1][] = sum of squares
  // [][0] = mass of mols in agg from options, [][1] = mass of the whole aggregate
  double mass_sum[2][2] = {0};
  // [0][] = simple sum, [1][] = sum of squares
  // [][0] = As (mols in agg from options), [][1] = total As
  double As_sum[2][2] = {0};
  // total mass disregarding all options -- regardless of options used, it should always be the same
  double volume_sum = 0;
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
      double agg_mass = 0;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type] == 1) {
          size++;
          agg_mass += MoleculeType[mol_type].Mass;
        }
      }
      // make calculations only if agg size is well defined and within given range
      if (size == 0 || size < range_As[0] || size > range_As[1]) {
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
        ndistr_mass[size-1][0] += agg_mass;
        wdistr_mass[size-1][0] += SQR(agg_mass);
        ndistr_mass[size-1][1] += Aggregate[i].Mass;
        wdistr_mass[size-1][1] += SQR(Aggregate[i].Mass);

        // ndistr_As is the number of aggregates with As=size as a function
        // of As (hopefully)
        ndistr_As[size-1]++;

        // really? this is weight average?
        wdistr_As[size-1][0] += size;
        wdistr_As[size-1][1] += Aggregate[i].nMolecules;

        voldistr[size-1] += Aggregate[i].nBeads; //}}}

        // number of various species in the aggregate //{{{
        count_agg[size-1]++;

        mass_sum[0][0] += agg_mass;
        mass_sum[1][0] += SQR(agg_mass);
        mass_sum[0][1] += Aggregate[i].Mass;
        mass_sum[1][1] += SQR(Aggregate[i].Mass);

        As_sum[0][0] += size;
        As_sum[1][0] += SQR(size);
        As_sum[0][1] += Aggregate[i].nMolecules;
        As_sum[1][1] += SQR(Aggregate[i].nMolecules);

        volume_sum += Aggregate[i].nBeads;

        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
          molecules_sum[size-1][mol_type]++;
        } //}}}
      }
    }

    // print averages to output file //{{{
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
    fprintf(out, " %8.3f", (double)(avg_mass_w_step[0])/avg_mass_n_step[0]); // <M>_w (options' mass)
    fprintf(out, " %8.3f", (double)(avg_mass_n_step[1])/aggs_step); // <M>_n (whole agg mass)
    fprintf(out, " %8.3f", (double)(avg_mass_w_step[1])/avg_mass_n_step[1]); // <M>_w (whole agg mass)

    fprintf(out, " %8.3f", (double)(avg_As_n_step[0])/aggs_step); // <As>_n (options' As)
    fprintf(out, " %8.3f", (double)(avg_As_w_step[0])/avg_As_n_step[0]); // <As>_w (options' As)
    fprintf(out, " %8.3f", (double)(avg_As_n_step[1])/aggs_step); // <As>_n (whole agg As)
    fprintf(out, " %8.3f", (double)(avg_As_w_step[1])/avg_As_n_step[1]); // <As>_w (whole agg As)

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
  // add date
  fprintf(out, " on %s\n", date);

  fprintf(out, "# ");
  fprintf(out, "1:A_s ");
  fprintf(out, "2:M_n(As) (options' mass) ");
  fprintf(out, "3:M_w(As) (options' mass) ");
  fprintf(out, "4:M_n(As) (whole agg mass) ");
  fprintf(out, "5:M_w(As) (whole agg mass) ");
  fprintf(out, "6:F_n(As) ");
  fprintf(out, "7:F_w(As) (options' As) ");
  fprintf(out, "8:F_w(As) (whole agg As) ");
  fprintf(out, "9:<volume distribution>_n ");
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

  // normalization factors for As distributions are calculated simply as a
  // sum of all the ndistr_As or wdistr_As values
  int norm_n = 0, norm_w[2] = {0};
  for (int i = 0; i < Counts.Molecules; i++) {
    norm_n += ndistr_As[i];

    norm_w[0] += wdistr_As[i][0];
    norm_w[1] += wdistr_As[i][1];
  }

  for (int i = 0; i < Counts.Molecules; i++) {
    /* number distribution of aggregate masses, M_n(As) =
     * ndistr_mass[i][]/agg_sum, is normalised by number average of mass, <M>_n
     * = mass_sum[0]/agg_sum, so M_n(As) = ndistr_mass[i][]/mass_sum[0][] to
     * give a sum of 1 over the whole of M_n(As) */
    /* weighed distribution of aggregate masses, M_w(As) =
     * wdistr_mass[i][]/mass_sum[0][], is normalised by weighed average of
     * mass, <M>_w = mass_sum[1][]/mass_sum[0][], so M_w(As) =
     * wdistr_mass[i][]/mass_sum[1][] to give a sum of 1 over the whole of
     * M_w(As) */
    fprintf(out, "%4d %lf %lf %lf %lf %lf %lf %lf %lf", i+1, // A_s
                (double)(ndistr_mass[i][0])/mass_sum[0][0], // number mass distr (options mass)
                (double)(wdistr_mass[i][0])/mass_sum[1][0], // weight mass distr (whole agg mass)
                (double)(ndistr_mass[i][1])/mass_sum[0][1], // number mass distr (options mass)
                (double)(wdistr_mass[i][1])/mass_sum[1][1], // weight mass distr (whole agg mass)
                (double)(ndistr_As[i])/norm_n, // number As distr
                (double)(wdistr_As[i][0])/norm_w[0], // number As distr
                (double)(wdistr_As[i][1])/norm_w[1], // weight As distr
                voldistr[i]/(volume_sum)); // volume distribution
    // print average number of molecule types in aggregates
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (count_agg[i] == 0) {
        fprintf(out, "%8s", "?");
      } else {
        fprintf(out, " %7.3f", (double)(molecules_sum[i][j])/count_agg[i]);
      }
    }
    // print total number of aggregates with given size
    fprintf(out, " %6d", count_agg[i]);
    putc('\n', out);
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
  fprintf(out, "1:<M>_n (options' mass) ");
  fprintf(out, "2:<M>_w (options' mass) ");
  fprintf(out, "3:<M>_n (whole agg mass) ");
  fprintf(out, "4:<M>_w (whole agg mass)");
  fprintf(out, "5:<As>_n (options' As) ");
  fprintf(out, "6:<As>_w (options' As) ");
  fprintf(out, "7:<As>_n (whole agg As) ");
  fprintf(out, "8:<As>_w (whole agg As)");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %d:<%s>_n", i+9, MoleculeType[i].Name);
  }
  putc('\n', out); //}}}

  // print the averages //{{{
  fprintf(out, "# ");
  fprintf(out, "%10.3f ", mass_sum[0][0]/count_agg[0]); // <M>_n (options' mass)
  fprintf(out, "%10.3f ", mass_sum[1][0]/mass_sum[0][0]); // <M>_w (whole agg mass)
  fprintf(out, "%10.3f ", mass_sum[0][1]/count_agg[0]); // <M>_w (options' mass)
  fprintf(out, "%10.3f ", mass_sum[1][1]/mass_sum[0][1]); // <M>_w (whole agg mass)
  fprintf(out, "%10.3f ", As_sum[0][0]/count_agg[0]); // <M>_n (options' mass)
  fprintf(out, "%10.3f ", As_sum[1][0]/As_sum[0][0]); // <M>_w (whole agg mass)
  fprintf(out, "%10.3f ", As_sum[0][1]/count_agg[0]); // <M>_w (options' mass)
  fprintf(out, "%10.3f ", As_sum[1][1]/As_sum[0][1]); // <M>_w (whole agg mass)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, "%7.3f", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
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
