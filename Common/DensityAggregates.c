#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <input.agg> <width> <output.rho> <agg sizes> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>     input filename (vcf format)\n");
  fprintf(stderr, "   <input.agg>     input filename with information about aggregates (agg format)\n");
  fprintf(stderr, "   <width>         width of a single bin\n");
  fprintf(stderr, "   <output.rho>    output density file (automatic ending '#.rho' added)\n");
  fprintf(stderr, "   <agg sizes>     aggregate sizes to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined     specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -n <int>     number of bins to average\n");
  fprintf(stderr, "      -st <int>    starting timestep for calculation\n");
  fprintf(stderr, "      -m <name>    agg size means number of <name> molecule types in an aggregate\n");
  fprintf(stderr, "      -x <name(s)> exclude specified molecule(s)\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
AggDensity utility calculates bead density for aggregates of given size(s) from \
their centre of mass. Beside unbonded beads it takes into account only beads \
from the current aggregate, not from any other aggregate.\n\n");

      fprintf(stdout, "\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input.vcf> <input.agg> <width> <output.rho> <agg sizes> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input.vcf>     input filename (vcf format)\n");
      fprintf(stdout, "   <input.agg>     input filename with information about aggregates (agg format)\n");
      fprintf(stdout, "   <width>         width of a single bin\n");
      fprintf(stdout, "   <output.rho>    output density file (automatic ending '#.rho' added)\n");
      fprintf(stdout, "   <agg sizes>     aggregate sizes to calculate density for\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      --joined     specify that aggregates with joined coordinates are used\n");
      fprintf(stdout, "      -n <int>     number of bins to average\n");
      fprintf(stdout, "      -st <int>    starting timestep for calculation\n");
      fprintf(stdout, "      -m <name>    agg size means number of <name> molecule types in an aggregate\n");
      fprintf(stdout, "      -x <name(s)> exclude specified molecule(s)\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 5; //}}}

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
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 ) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
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

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // number of bins to average //{{{
  int avg = 1;
  if (IntegerOption(argc, argv, "-n", &avg)) {
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

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <width>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.rho> - filename with bead densities //{{{
  char output_rho[32];
  strcpy(output_rho, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information{{{
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file); //}}}

  // '-m' option //{{{
  int specific_moltype_for_size;
  if (AggSizeSpecificMolType(argc, argv, &specific_moltype_for_size, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  int *other_mols = malloc(Counts.Molecules*sizeof(int));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
    other_mols[i] = 0;
  }

  int aggs = 0;

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (argv[count][0] < '1' || argv[count][0] > '9') {
      fprintf(stderr, "Non-numeric option in <agg sizes>!\n");
      exit(1);
    } //}}}

    agg_sizes[aggs][0] = atoi(argv[count]);

    // write initial stuff to output density file //{{{
    FILE *out;
    char str[32];

    if (specific_moltype_for_size == -1) {
      sprintf(str, "%s%d.rho", output_rho, agg_sizes[aggs][0]);
    } else {
      sprintf(str, "%s%s%d.rho", output_rho, MoleculeType[specific_moltype_for_size].Name, agg_sizes[aggs][0]);
    }
    if ((out = fopen(str, "w")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", str);
      exit(1);
    }

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++)
      fprintf(out, " %s", argv[i]);
    putc('\n', out); //}}}

    // print bead type names to output file //{{{
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      fprintf(out, " %s", BeadType[i].Name);
    }
    putc('\n', out); //}}}

    fclose(out); //}}}

    aggs++; // number of aggregate sizes
  } //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  // open for the first time to read the distance and type names
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

    BeadType[type].Use = true;

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
  // skip line with command to produce the agg file
  while (getc(agg) != '\n')
    ;
  // skip empty line
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

  // number of bins //{{{
  double max_dist = 0.5 * Min3(BoxLength.x, BoxLength.y, BoxLength.z);
  int bins = ceil(max_dist / width); //}}}

  // allocate memory for density arrays //{{{
  double ***rho = malloc(Counts.TypesOfBeads*sizeof(double **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    rho[i] = malloc(aggs*sizeof(double *));
    for (int j = 0; j < aggs; j++) {
      rho[i][j] = calloc(bins,sizeof(double));
    }
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

    fprintf(stdout, "Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      fprintf(stdout, " %d", agg_sizes[i][0]);
    }
    putchar('\n');
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // skip first start-1 steps //{{{
  count = 0;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent) {
      if (script) {
        fprintf(stdout, "Discarding step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rDiscarding step: %6d", count);
      }
    } //}}}

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);
    if (SkipCoor(vcf, Counts, &stuff) == 1) {
      fprintf(stderr, "Premature end of %s file!\n", input_vcf);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Discarded steps: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rDiscarded steps: %6d\n", count);
    }
  } //}}}
  //}}}

  // main loop //{{{
  count = 0; // count timesteps
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
      }
    } //}}}

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

    // calculate densities //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {

      // test if aggregate 'i' should be used //{{{
      int size = 0;
      if (specific_moltype_for_size != -1) { // agg size = number of molecules of type 'specific_moltype_for_size'
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int id = Aggregate[i].Molecule[j];
          if (specific_moltype_for_size == Molecule[id].Type) {
            size++;
          }
        }
      } else { // agg size = total number of all molecules
        size = Aggregate[i].nMolecules;
      }
      // is 'size' in provided list?
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j][0] == size) {
          correct_size = j;
        }
      } //}}}

      if (correct_size != -1) {
        other_mols[correct_size] += Aggregate[i].nMolecules - size;

        Vector com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, Bead, BeadType);

        // aggregate beads //{{{
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int bead = Aggregate[i].Bead[j];
          int mol = Bead[bead].Molecule;
          int moltype = Molecule[mol].Type;

          if (MoleculeType[moltype].Use) {
//          fprintf(stdout, "%s: %d\n", MoleculeType[moltype].Name, MoleculeType[moltype].Use);

            Vector dist = Distance(Bead[Aggregate[i].Bead[j]].Position, com, BoxLength);
            dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

            if (dist.x < max_dist) {
              int k = dist.x / width;

              rho[Bead[Aggregate[i].Bead[j]].Type][correct_size][k]++;
            }
          }
        } //}}}

        // monomeric beads //{{{
        for (int j = 0; j < Counts.Unbonded; j++) {
          Vector dist = Distance(Bead[j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            rho[Bead[j].Type][correct_size][k]++;
          }
        } //}}}

        agg_sizes[correct_size][1]++;
      }
    } //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      fprintf(stdout, "\n%s", stuff);
    } //}}}
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

  // write densities to output file(s) //{{{
  for (int i = 0; i < aggs; i++) {
    FILE *out;
    char str[32];

    if (specific_moltype_for_size == -1) {
      sprintf(str, "%s%d.rho", output_rho, agg_sizes[i][0]);
    } else {
      sprintf(str, "%s%s%d.rho", output_rho, MoleculeType[specific_moltype_for_size].Name, agg_sizes[i][0]);
    }
    if ((out = fopen(str, "a")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", str);
      exit(1);
    }

    // calculate rdf
    for (int j = 0; j < bins; j++) {

      // calculate volume of every shell that will be averaged
      double shell[avg];
      for (int k = 0; k < avg; k++) {
        shell[k] = 4 * PI * CUBE(width) *(CUBE(j+k+1) - CUBE(j+k)) / 3;
      }

      fprintf(out, "%.2f", width*(j+0.5*avg));

      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        double temp = 0;

        // sum up rdfs from all shells to be averaged
        for (int l = 0; l < avg; l++) {
          temp += rho[k][i][j+l] / (shell[l] * agg_sizes[i][1]);
        }

        // print average value to output file
        fprintf(out, " %10f", temp/avg);
      }
      putc('\n',out);
    }

    if (specific_moltype_for_size != -1) {
      fprintf(out, "# %d %s molecules in aggregate and %.2f other molecules (%d aggregates)\n",
          agg_sizes[i][0], MoleculeType[specific_moltype_for_size].Name, (double)(other_mols[i])/agg_sizes[i][1], agg_sizes[i][1]);
    } else {
      fprintf(out, "# %d molecules in aggregate (%d aggregates)\n", agg_sizes[i][0], agg_sizes[i][1]);
    }
    fclose(out);
  } //}}}

  // print to stdout number of unspecified molecules if '-m' option was used //{{{
  if (specific_moltype_for_size != -1 && !silent) {
    for (int i = 0; i < aggs; i++) {
      fprintf(stdout, "%d %s molecules in aggregate and %.2f other molecules (%d aggregates)\n",
          agg_sizes[i][0], MoleculeType[specific_moltype_for_size].Name, (double)(other_mols[i])/agg_sizes[i][1], agg_sizes[i][1]);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(other_mols);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
  }
  free(agg_sizes);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < aggs; j++) {
      free(rho[i][j]);
    }
    free(rho[i]);
  }
  free(rho); //}}}

  return 0;
}
