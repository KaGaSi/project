#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input.agg> <output> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <input.agg>         input filename with information about aggregates (agg format)\n");
  fprintf(stderr, "   <output>            output file with data during simulation run\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined         specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -bt              specify bead types to be used for calculation (default is all)\n");
  fprintf(stderr, "      -m <name>        agg size means number of <name> molecule types in an aggregate\n");
  fprintf(stderr, "      --no-unimers      do not count unimers into averages\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
GyrationAggregates calculates radii of gyration, acylindricities, \
asphericities and relative shape anisotropies during the simulation for all \
aggregates izes. The shape descriptors are calculated from eigenvalues of \
gyration tensor. It also prints averages to the stdout. Instead of aggregate \
size, a number of specified molecular species in an aggregate can be used and \
only specified bead types can be used for all calculations.\n\n");

      fprintf(stdout, "\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input.vcf> <input.agg> <output> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input.vcf>         input filename (vcf format)\n");
      fprintf(stdout, "   <input.agg>         input filename with information about aggregates (agg format)\n");
      fprintf(stdout, "   <output>            output file with data during simulation run\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      --joined         specify that aggregates with joined coordinates are used\n");
      fprintf(stdout, "      -bt              specify bead types to be used for calculation (default is all)\n");
      fprintf(stdout, "      -m <name>        agg size means number of <name> molecule types in an aggregate\n");
      fprintf(stdout, "      --no-unimers      do not count unimers into averages\n");
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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "--no-unimers") != 0) {

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
    fprintf(stderr, "Too little mandatory arguments (%d instead of at least %d)!\n\n", count, options);
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

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

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

  // <input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf[32];
  strcpy(input_vcf, argv[++count]);

  // test if <input.vcf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(input_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<input.vcf> '%s' does not have .vcf ending!\n", input_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <input.agg> - filename of input file with aggregate information //{{{
  char input_agg[32];
  strcpy(input_agg, argv[++count]); //}}}

  // <output> - filename with data during simulation run //{{{
  char output[32];
  strcpy(output, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file);

  // -m <name> option - specify MoleculeType that is used for determining agg sizes //{{{
  int specific_molecule = -1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-m") == 0) {

      // Error - missing or wrong argument //{{{
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument for '-m' option (or molecule name beginning with a dash)!\n");
        ErrorHelp(argv[0]);
        exit(1);
      } //}}}

      specific_molecule = FindMoleculeType(argv[i+1], Counts, MoleculeType);
      if (specific_molecule == -1) {
        fprintf(stderr, "Molecule '%s' does not exist in FIELD ('-m' option)!\n", argv[i+1]);
        exit(1);
      }
    }
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, Counts, &BeadType)) {
    exit(0);
  } //}}}

  // write initial stuff to output file //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output);
    exit(1);
  }

  // print command to output file
  putc('#', out);
  for (int i = 0; i < argc; i++) {
    fprintf(out, " %s", argv[i]);
  }
  // add date
  char date[26];
  Date(date);
  fprintf(out, " on %s\n", date);

  // print legend line to output file
  fprintf(out, "# 1:dt 2:<Rg>_n 3:<Rg>_w 4:<Rg>_z ");
  fprintf(out, "5:<Rg^2>_n 6:<Rg^2>_w 7:<Rg^2>_z ");
  fprintf(out, "8:<Anis>_n 9:<Acyl>_n 10:<Aspher>_n\n");

  fclose(out); //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  // read minimum distance for closeness check (<distance> argument in Aggregates utility)
  double distance;
  fscanf(agg, "%*s %*s %lf", &distance);

  // skip <contacts> and <output.agg> in Aggregates command
  fscanf(agg, "%*s %*s");

  // read <type names> from Aggregates command //{{{
  int test;
  // reading ends if next argument (beginning with '-') or the following empty line is read
  while ((test = getc(agg)) != '-' && test != '\n') {
    ungetc(test, agg);

    char name[10];
    fscanf(agg, "%s", name);
    int type = FindBeadType(name, Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", name, input_vcf);
      exit(1);
    }

    while ((test = getc(agg)) == ' ')
      ;
    ungetc(test, agg);
  } //}}}
  fclose(agg);

  // open again for production run - to ensure the pointer position in file is correct (at first 'Step')
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf);
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read string from '%s' file!\n", input_vcf);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_vcf);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;
  // skip blank line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
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
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // allocate memory for sum of various things //{{{

  // numbers of aggregates of all possibe sizes
  int *agg_counts_sum = calloc(Counts.Molecules,sizeof(int *));
  // total radius of gyration: [size][0] normal sum, [size][1] sum of Rg*mass, [size][2] Rg*mass^2
  double **Rg_sum = malloc(Counts.Molecules*sizeof(double *));
  // total square of radius of gyration: [size][0] normal sum, [size][1] sum of Rg^2*mass, [size][2] Rg^2*mass^2
  double **sqrRg_sum = malloc(Counts.Molecules*sizeof(double *));
  // relative shape anisotropy: only normal sum
  double *Anis_sum = calloc(Counts.Molecules,sizeof(double));
  // acylindricity: only normal sum
  double *Acyl_sum = calloc(Counts.Molecules,sizeof(double));
  // asphericity: only normal sum
  double *Aspher_sum = calloc(Counts.Molecules,sizeof(double));
  // total mass of aggregates: [size][0] normal sum, [size][1] sum of squares
  long int **mass_sum = malloc(Counts.Molecules*sizeof(int *));
  // number of molecule types in aggregates: [size][mol type] only normal sum
  int **molecules_sum = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    Rg_sum[i] = calloc(3,sizeof(double));
    sqrRg_sum[i] = calloc(3,sizeof(double));
    mass_sum[i] = calloc(2,sizeof(long int));
    molecules_sum[i] = calloc(Counts.TypesOfMolecules,sizeof(int));
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
      }
    }

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf, Counts, &Bead, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      }
    } //}}}

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate arrays for the timestep //{{{
    int *agg_counts_step = calloc(Counts.Molecules,sizeof(int));
    double **Rg_step = malloc(Counts.Molecules*sizeof(double *));
    double **sqrRg_step = malloc(Counts.Molecules*sizeof(double *));
    double *Anis_step = calloc(Counts.Molecules,sizeof(double));
    double *Acyl_step = calloc(Counts.Molecules,sizeof(double));
    double *Aspher_step = calloc(Counts.Molecules,sizeof(double));
    for (int i = 0; i < Counts.Molecules; i++) {
      Rg_step[i] = calloc(3,sizeof(double));
      sqrRg_step[i] = calloc(3,sizeof(double));
    } //}}}

    // calculate shape descriptors //{{{
    double mass_step[2] = {0}; // total mass of aggregates in a step: [0] normal, [1] sum of squares
    for (int i = 0; i < Counts.Aggregates; i++) {
      if (!no_uni || Aggregate[i].nMolecules != 1) {

        // test if aggregate 'i' should be used //{{{
        int mols = 0; // agg size
        if (specific_molecule != -1) { // agg size = number of molecules of type 'specific_molecule'
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int id = Aggregate[i].Molecule[j];
            if (specific_molecule == Molecule[id].Type) {
              mols++;
            }
          }
        } else { // agg size = total number of all molecules
          mols = Aggregate[i].nMolecules;
        }
        // is 'mols' agg size in provided list?
        int correct_size = -1;
        for (int j = 0; j < Counts.Molecules; j++) {
          if ((j+1) == mols) {
            correct_size = j;
          }
        } //}}}

        if (correct_size != -1) {
          agg_counts_step[correct_size]++;
          agg_counts_sum[correct_size]++;

          // copy bead ids to a separate array //{{{
          int *list = malloc(Aggregate[i].nBeads*sizeof(int));
          int n = 0;
          for (int j = 0; j < Aggregate[i].nBeads; j++) {
            int id = Aggregate[i].Bead[j];
            if (BeadType[Bead[id].Type].Use) {
              list[n] = id;
              n++;
            }
          } //}}}

          Vector eigen = Gyration(n, list, Counts, BoxLength, BeadType, &Bead);

  //      // calcule Rg the 'usual way' -- for testing purposes //{{{
  //      double Rg2 = 0;
  //      Vector com = CentreOfMass(n, list, Bead, BeadType);
  //      for (int j = 0; j < n; j++) {
  //        Vector rij = Distance(Bead[list[j]].Position, com, BoxLength);
  //        Rg2 += SQR(rij.x) + SQR(rij.y) + SQR(rij.z);
  //      }
  //      Rg2 /= n; //}}}

          free(list); // free array of bead ids for gyration calculation

          double Rgi = sqrt(eigen.x + eigen.y + eigen.z);

          // agg masses
          mass_step[0] += Aggregate[i].Mass; // for this timestep
          mass_step[1] += SQR(Aggregate[i].Mass); // for this timestep
          // radius of gyration
          Rg_step[correct_size][0] += Rgi; // for number avg
          Rg_step[correct_size][1] += Rgi * Aggregate[i].Mass; // for weight average
          Rg_step[correct_size][2] += Rgi * SQR(Aggregate[i].Mass); // for z-average
          // squared radius of gyration
          sqrRg_step[correct_size][0] += SQR(Rgi); // for number avg
          sqrRg_step[correct_size][1] += SQR(Rgi) * Aggregate[i].Mass; // for weight average
          sqrRg_step[correct_size][2] += SQR(Rgi) * SQR(Aggregate[i].Mass); // for z-average
          // relative shape anisotropy
          Anis_step[correct_size] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) / SQR(eigen.x + eigen.y + eigen.z) - 0.5;
          // acylindricity
          Acyl_step[correct_size] += eigen.y - eigen.x;
          // asphericity
          Aspher_step[correct_size] += eigen.z - 0.5 * (eigen.x + eigen.y);
          // aggregate size
          mass_sum[correct_size][0] += Aggregate[i].Mass;
          mass_sum[correct_size][1] += SQR(Aggregate[i].Mass);

          // count number of Diblocks and Surfacts in aggrate
          // TODO: generalise
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
            molecules_sum[correct_size][mol_type]++;
          }
        }
      }
    } //}}}

    // add values to sums //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      Rg_sum[i][0] += Rg_step[i][0];
      Rg_sum[i][1] += Rg_step[i][1];
      Rg_sum[i][2] += Rg_step[i][2];
      sqrRg_sum[i][0] += sqrRg_step[i][0];
      sqrRg_sum[i][1] += sqrRg_step[i][1];
      sqrRg_sum[i][2] += sqrRg_step[i][2];
      Anis_sum[i] += Anis_step[i];
      Acyl_sum[i] += Acyl_step[i];
      Aspher_sum[i] += Aspher_step[i];
    } //}}}

    // print data to output file //{{{
    FILE *out;
    if ((out = fopen(output, "a")) == NULL) { // out file opened fine? //{{{
      // print newline to stdout if Step... doesn't end with one
      if (!script && !silent) {
        putchar('\n');
      }
      fprintf(stderr, "Cannot open file %s!\n", output);
      exit(1);
    } //}}}

    fprintf(out, "%5d", count); // timestep
    // sum up contributions from all aggregate sizes
    for (int i = 1; i < Counts.Molecules; i++) {
      Rg_step[0][0] += Rg_step[i][0];
      Rg_step[0][1] += Rg_step[i][1];
      Rg_step[0][2] += Rg_step[i][2];
      sqrRg_step[0][0] += sqrRg_step[i][0];
      sqrRg_step[0][1] += sqrRg_step[i][1];
      sqrRg_step[0][2] += sqrRg_step[i][2];
      Anis_step[0] += Anis_step[i];
      Acyl_step[0] += Acyl_step[i];
      Aspher_step[0] += Aspher_step[i];

      agg_counts_step[0] += agg_counts_step[i];
    }
    // <R_G>
    fprintf(out, " %8.5f %8.5f %8.5f", Rg_step[0][0]/agg_counts_step[0], Rg_step[0][1]/mass_step[0], Rg_step[0][2]/mass_step[1]);
    // <R_G^2>
    fprintf(out, " %8.5f %8.5f %8.5f", sqrRg_step[0][0]/agg_counts_step[0], sqrRg_step[0][1]/mass_step[0], sqrRg_step[0][2]/mass_step[1]);
    // relative shape anisotropy
    fprintf(out, " %8.5f", Anis_step[0]/agg_counts_step[0]);
    // acylindricity
    fprintf(out, " %8.5f", Acyl_step[0]/agg_counts_step[0]);
    // asphericity
    fprintf(out, " %8.5f", Aspher_step[0]/agg_counts_step[0]);
    putc('\n', out);

    fclose(out); //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      fprintf(stdout, "\n%s", stuff);
    } //}}}

    // free memory //{{{
    free(agg_counts_step);
    for (int i = 0; i < Counts.Molecules; i++) {
      free(Rg_step[i]);
      free(sqrRg_step[i]);
    }
    free(Rg_step);
    free(sqrRg_step);
    free(Anis_step);
    free(Acyl_step);
    free(Aspher_step); //}}}
  }
  fclose(vcf);
  fclose(agg);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  } //}}}

  // calculate per-size averages //{{{
  fprintf(stdout, "# 1:As");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(stdout, " %d:<%s>", i+2, MoleculeType[i].Name);
  }
  fprintf(stdout, " %d:<Rg>", Counts.TypesOfMolecules+2);
  fprintf(stdout, " %d:<Rg^2>", Counts.TypesOfMolecules+3);
  fprintf(stdout, " %d:<Acyl>", Counts.TypesOfMolecules+4);
  fprintf(stdout, " %d:<Ashper>", Counts.TypesOfMolecules+5);
  fprintf(stdout, " %d:<Anis>", Counts.TypesOfMolecules+6);
  fprintf(stdout, " %d:number of aggs", Counts.TypesOfMolecules+7);
  putchar('\n');
  for (int i = 0; i < Counts.Molecules; i++) {
    if (agg_counts_sum[i] > 0) {
      fprintf(stdout, "%4d", i+1);
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        fprintf(stdout, " %7.3f", (double)(molecules_sum[i][j])/agg_counts_sum[i]);
      }
      fprintf(stdout, " %7.3f", Rg_sum[i][0]/agg_counts_sum[i]);
      fprintf(stdout, " %7.3f", sqrRg_sum[i][0]/agg_counts_sum[i]);
      fprintf(stdout, " %7.3f", Anis_sum[i]/agg_counts_sum[i]);
      fprintf(stdout, " %7.3f", Acyl_sum[i]/agg_counts_sum[i]);
      fprintf(stdout, " %7.3f", Aspher_sum[i]/agg_counts_sum[i]);
      fprintf(stdout, " %d", agg_counts_sum[i]);
      putchar('\n');
    }
  } //}}}

  // total averages //{{{
  for (int i = 1; i < Counts.Molecules; i++) {
    Rg_sum[0][0] += Rg_sum[i][0];
    Rg_sum[0][1] += Rg_sum[i][1];
    Rg_sum[0][2] += Rg_sum[i][2];
    sqrRg_sum[0][0] += sqrRg_sum[i][0];
    sqrRg_sum[0][1] += sqrRg_sum[i][1];
    sqrRg_sum[0][2] += sqrRg_sum[i][2];
    Anis_sum[0] += Anis_sum[i];
    Acyl_sum[0] += Acyl_sum[i];
    Aspher_sum[0] += Aspher_sum[i];

    agg_counts_sum[0] += agg_counts_sum[i];

    mass_sum[0][0] += mass_sum[i][0];
    mass_sum[0][1] += mass_sum[i][1];

    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules_sum[0][j] += molecules_sum[i][j];
    }
  }

  // print to output file
  if ((out = fopen(output, "a")) == NULL) { //{{{
    fprintf(stderr, "Cannot open file %s!\n", output);
    exit(1);
  } //}}}

  fprintf(out, "# 1:<M>_n 2:_w ");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, "%d:<%s> ", i+3, MoleculeType[i].Name);
  }
  fprintf(out, "%d:<Rg>_n ", Counts.TypesOfMolecules+3);
  fprintf(out, "%d:_w ", Counts.TypesOfMolecules+4);
  fprintf(out, "%d:_z ", Counts.TypesOfMolecules+5);
  fprintf(out, "%d:<Rg^2>_n ", Counts.TypesOfMolecules+6);
  fprintf(out, "%d:_w ", Counts.TypesOfMolecules+7);
  fprintf(out, "%d:_z ", Counts.TypesOfMolecules+8);
  fprintf(out, "%d:<Anisotropy> ", Counts.TypesOfMolecules+9);
  fprintf(out, "%d:<Acylindricity> ", Counts.TypesOfMolecules+10);
  fprintf(out, "%d:<Asphericity>\n", Counts.TypesOfMolecules+11);
  fprintf(out, "# %8.3f", (double)(mass_sum[0][0])/agg_counts_sum[0]); //<M_As>_n
  fprintf(out, " %8.3f", (double)(mass_sum[0][1])/mass_sum[0][0]); //<M_As>_w
  // molecule types
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %8.3f", (double)(molecules_sum)[0][i]/agg_counts_sum[0]);
  }
  fprintf(out, " %8.3f", Rg_sum[0][0]/agg_counts_sum[0]); // <Rg>_n
  fprintf(out, " %8.3f", Rg_sum[0][1]/mass_sum[0][0]); // <Rg>_w
  fprintf(out, " %8.3f", Rg_sum[0][2]/mass_sum[0][1]); // <Rg>_z
  fprintf(out, " %8.3f", sqrRg_sum[0][0]/agg_counts_sum[0]); // <Rg^2>_n
  fprintf(out, " %8.3f", sqrRg_sum[0][1]/mass_sum[0][0]); // <Rg^2>_w
  fprintf(out, " %8.3f", sqrRg_sum[0][2]/mass_sum[0][1]); // <Rg^2>_z
  fprintf(out, " %8.3f", Anis_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %8.3f", Acyl_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %8.3f\n", Aspher_sum[0]/agg_counts_sum[0]);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(Rg_sum[i]);
    free(sqrRg_sum[i]);
    free(mass_sum[i]);
    free(molecules_sum[i]);
  }
  free(molecules_sum);
  free(mass_sum);
  free(agg_counts_sum);
  free(Rg_sum);
  free(sqrRg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(stuff); //}}}

  return 0;
}
