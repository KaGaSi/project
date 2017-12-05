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

  fprintf(stderr, "   <input>              input filename (agg format)\n");
  fprintf(stderr, "   <output distr file>  filename with weight and number distributions\n");
  fprintf(stderr, "   <output avg file>    filename with weight and number averages throughout simulation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -st <int>         start distribution calculation with <int>-th step\n");
  fprintf(stderr, "      --no-unimers      do not count unimers into averages\n");
  fprintf(stderr, "      -m <name>         agg size means number of <name> molecule types in an aggregate\n");
  fprintf(stderr, "      -x <name(s)>      exclude aggregates containing only specified molecule(s)\n");
  fprintf(stderr, "      --only <name>     use just aggregates composed of a specified molecule\n");
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

      fprintf(stdout, "   <input>              input filename (agg format)\n");
      fprintf(stdout, "   <output distr file>  filename with weight and number distributions\n");
      fprintf(stdout, "   <output avg file>    filename with weight and number averages throughout simulation\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -st <int>         start distribution calculation with <int>-th step\n");
      fprintf(stdout, "      --no-unimers      do not count unimers into averages\n");
      fprintf(stdout, "      -m <name>         agg size means number of <name> molecule types in an aggregate\n");
      fprintf(stdout, "      -x <name(s)>      exclude aggregates containing only specified molecule(s)\n");
      fprintf(stdout, "      --only <name>     use just aggregates composed of a specified molecule\n");
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
        strcmp(argv[i], "--no-unimers") != 0 &&
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

  // do not count unimers towards averages //{{{
  bool no_uni = BoolOption(argc, argv, "--no-unimers"); //}}}
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

  // '-m' option //{{{
  int specific_moltype_for_size;
  if (MoleculeTypeOption(argc, argv, "-m", &specific_moltype_for_size, Counts, &MoleculeType)) {
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

  // disused, really - but it's still in the main loop and may some day come in useful //{{{
  int max_mass = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    max_mass += MoleculeType[i].Number * MoleculeType[i].Mass;
  }
  long int ndistr_mass[max_mass];
  long int wdistr_mass[max_mass];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    ndistr_mass[i] = 0;
    wdistr_mass[i] = 0;
  } //}}}

  // arrays for distribution //{{{
  long int ndistr_As[Counts.Molecules];
  long int wdistr_As[Counts.Molecules];
  long int voldistr_As[Counts.Molecules];
  // number of aggregates throughout simulation
  int count_agg[Counts.Molecules];
  // molecule typs in aggregates: [agg size][mol type][number or SQR(number)]
  int **molecules_sum = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    molecules_sum[i] = calloc(Counts.TypesOfMolecules,sizeof(int));
  }

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    wdistr_As[i] = 0;
    ndistr_As[i] = 0;
    voldistr_As[i] = 0;
    count_agg[i] = 0;
  }
  //}}}

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

  fprintf(out, "# 1:step 2:<M>_n 3:<M>_w\n");
  fclose(out); //}}}

  // main loop //{{{
  int test, agg_sum = 0;// sum_mass = 0;
  count = 0;
  double mass_sum[2] = {0}; // [0] = simple sum, [1] = sum of squares
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
        avg_n_step = 0, avg_w_step = 0; // per-step averages
    for (int i = 0; i < Counts.Aggregates; i++) {

      // determine aggregate size -- if `-m` is used size is the number of 'specific_moltype_for_size' mols //{{{
      int size = 0;
      double agg_mass = 0;
      if (specific_moltype_for_size != -1) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
          if (specific_moltype_for_size == mol_type) {
            size++;
            agg_mass += MoleculeType[mol_type].Mass;
          }
        }
      } else { // agg size = total number of all molecules
        size = Aggregate[i].nMolecules;
        agg_mass = Aggregate[i].Mass;
      } //}}}

      // if '--only' is used, use only aggregates composed the specified molecule
      bool use = true;
      if (only_specific_moltype_aggregates != -1) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (only_specific_moltype_aggregates != Molecule[id].Type) {
            use = false;
            break;
          }
        }
      }

      // make calculations only if agg size is well defined
      if (size > 0 && use == true) {
        // start calculation of averages from specified 'start' timestep
        if (count >= start) {
//        printf("%d %lf\n", size, agg_mass);

          // if '-x' option is used, discount aggregates with only specified molecule types //{{{
          bool test = true; // aggregate with only excluded molecules?
          for (int j = 0; j < size; j++) {
            if (MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Use) {
              test = false;
              break;
            }
          }
          // should the rest of the for loop be skipped?
          if (test) {
            continue;
          } //}}}

          // distribution //{{{
          ndistr_As[size-1] += agg_mass;
          wdistr_As[size-1] += SQR(agg_mass);
          voldistr_As[size-1] += Aggregate[i].nBeads;
          int mass = agg_mass;
          ndistr_mass[mass] += agg_mass;
          wdistr_mass[mass] += SQR(agg_mass); //}}}

          // number of various species in the aggregate //{{{
          if (!no_uni || size != 1) {
            count_agg[size-1]++;
            mass_sum[0] += agg_mass;
            mass_sum[1] += SQR(agg_mass);

            volume_sum += Aggregate[i].nBeads;

            int *mol_count = calloc(Counts.TypesOfMolecules,sizeof(int));

            for (int j = 0; j < Aggregate[i].nMolecules; j++) {
              mol_count[Molecule[Aggregate[i].Molecule[j]].Type]++;
            }

            for (int j = 0; j < Counts.TypesOfMolecules; j++) {
              molecules_sum[size-1][j] += mol_count[j];
            }

            free(mol_count);
          } //}}}
        }

        // average aggregation number during the step //{{{
        if (!no_uni || size != 1) {
          aggs_step++;

          avg_n_step += agg_mass;
          avg_w_step += SQR(agg_mass);
        } //}}}
      }
//    else {
//      fprintf(stdout, "ok\n");
//    }
    }

    // add to total number of aggregates
    agg_sum += aggs_step;

    // print averages to output file //{{{
    if ((out = fopen(output_avg, "a")) == NULL) {
      // print newline to stdout if Step... doesn't end with one
      if (!script && !silent) {
        putchar('\n');
      }
      fprintf(stderr, "Cannot open file %s!\n", output_avg);
      exit(1);
    }

    fprintf(out, "%5d %8.3f %8.3f\n", count, (double)(avg_n_step)/aggs_step, (double)(avg_w_step)/avg_n_step);
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

  fprintf(out, "# 1:A_s 2:<M>_n 3:<M>_w 3:<volume distribution>_n ");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %d:<%s>_n", i+4, MoleculeType[i].Name);
  }
  putc('\n', out);
  fclose(out); //}}}

  // print distributions to output file //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  count -= start -1;
  for (int i = 0; i < Counts.Molecules; i++) {
    /* number distribution of aggregate masses, F(M_i)_n =
     * ndistr_As[i]/agg_sum, is normalised by number average of mass, <M>_n =
     * mass_sum[0]/agg_sum, so F(M_i)_n = ndistr_As[i]/mass_sum[0] to give a
     * sum of 1 overall the whole of F(M_i)_w */
    /* weighed distribution of aggregate masses, F(M_i)_w =
     * wdistr_As[i]/mass_sum[0], is normalised by weighed average of mass, <M>_w =
     * mass_sum[1]/mass_sum[0], so F(M_i)_w = wdistr_As[i]/mass_sum[1] to give a
     * sum of 1 overall the whole of F(M_i)_w */
    fprintf(out, "%4d %lf %lf %lf", i+1, // A_s
                                    (double)(ndistr_As[i])/mass_sum[0], // number distribution
                                    (double)(wdistr_As[i])/mass_sum[1], // weight distribution
                                    voldistr_As[i]/(volume_sum)); // volume distribution
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

  // disused - still something in main loop; calculates F as function of mass //{{{
  //out = fopen("distr-mass.txt", "w");
  //for (int i = 0; i < max_mass; i++) {
  //  if (ndistr_mass[i] > 0) {
  //    fprintf(out, "%5d %lf %lf\n", i, // mass
  //                                  (double)(ndistr_mass[i])/mass_sum[0], // number distribution
  //                                  (double)(wdistr_mass[i])/mass_sum[1]); // weight distribution
  //  }
  //}
  //fclose(out); //}}}

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
  fprintf(out, "# 1:<M>_n 2:<M>_w");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %d:<%s>_n", i+3, MoleculeType[i].Name);
  }
  putc('\n', out); //}}}

  // print the averages //{{{
  fprintf(out, "# %10.3f", (double)(mass_sum[0])/count_agg[0]); // <A_s>_n
  fprintf(out, "%10.3f", mass_sum[1]/mass_sum[0]); // <Mass>_w
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %7.3f", (double)(molecules_sum[0][i])/count_agg[0]); // <species>_n
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
