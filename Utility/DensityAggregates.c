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
    fprintf(ptr, "\
DensityAggregates utility calculates radial bead density for aggregates \
of given size(s) from their centre of mass. For beads in molecules, \
it takes into account only beads from the current aggregate, not from \
any other aggregate. The utility does not check what molecule type a \
given bead belongs to, therefore if the same bead type appears in \
more molecule types, the resulting density for that bead type will be \
averaged without regard for the various types of molecule it comes from \
(i.e., if 'mol1' and 'mol2' both both contain bead 'A', there will not be \
a separate column for 'A' from 'mol1' and from 'mol2').\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <input.agg> <width> <output.rho> <agg size(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <input.agg>       input agg file\n");
  fprintf(ptr, "   <width>           width of a single bin of the distribution\n");
  fprintf(ptr, "   <output>          output density file with automatic '#.rho' ending (# is aggregate size)\n");
  fprintf(ptr, "   <agg size(s)>     aggregate size(s) to calculate density for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined       specify that aggregates with joined coordinates are used\n");
  fprintf(ptr, "      -n <int>       number of bins to average\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> molecule types in an aggregate\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
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
//      strcmp(argv[i], "-b") != 0 &&
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
  char extension[2][5];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_vsf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
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

  // <input> - input coordinate file //{{{
  char input_coor[1024];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  ext = 2;
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <input.agg> - input agg file with aggregate information //{{{
  char input_agg[1024];
  strcpy(input_agg, argv[++count]);

  // test if <input.gg> ends with '.agg'
  ext = 1;
  strcpy(extension[0], ".agg");
  if (!ErrorExtension(input_agg, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.rho> - filename with bead densities //{{{
  char output_rho[1024];
  strcpy(output_rho, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

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
  } //}}}

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  int **agg_mols = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
    agg_mols[i] = calloc(Counts.TypesOfMolecules, sizeof(int));
  }

  int aggs = 0;

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (argv[count][0] < '1' || argv[count][0] > '9') {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}

    agg_sizes[aggs][0] = atoi(argv[count]);

    // write initial stuff to output density file //{{{
    FILE *out;
    char str[1024];
    strcpy(str, output_rho);

    char str2[1030];
    sprintf(str2, "%s%d.rho", str, agg_sizes[aggs][0]);
    strcpy(str, str2);
    if ((out = fopen(str, "w")) == NULL) {
      ErrorFileOpen(str, 'w');
      exit(1);
    }

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++)
      fprintf(out, " %s", argv[i]);
    putc('\n', out); //}}}

    // print bead type names to output file //{{{
    fprintf(out, "# for each bead type: (1) rdp; (2) stderr; (3) rnp; (4) stderr\n");
    fprintf(out, "# columns: (1) distance;");
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      fprintf(out, " (%d) %s", 4*i+2, BeadType[i].Name);
      if (i != (Counts.TypesOfBeads-1)) {
        putc(';', out);
      }
    }
    putc('\n', out); //}}}

    fclose(out); //}}}

    aggs++; // number of aggregate sizes
  } //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  // open for the first time to read the distance and type names
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
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

    char name[1024];
    fscanf(agg, "%s", name);
    int type = FindBeadType(name, Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", name, input_coor);
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
    ErrorFileOpen(input_agg, 'r');
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
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[1024];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "\nError: cannot read pbc from %s\n\n", input_coor);
    exit(1);
  }

  // skip remainder of pbc line
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
  double ***rho_2 = malloc(Counts.TypesOfBeads*sizeof(double **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    rho[i] = malloc(aggs*sizeof(double *));
    rho_2[i] = malloc(aggs*sizeof(double *));
    for (int j = 0; j < aggs; j++) {
      rho[i][j] = calloc(bins,sizeof(double));
      rho_2[i][j] = calloc(bins,sizeof(double));
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(1024*sizeof(int));

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
    VerboseOutput(verbose2, input_coor, Counts, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      fprintf(stdout, " %d", agg_sizes[i][0]);
    }
    putchar('\n');
  } //}}}

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

    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      fprintf(stderr, "Error: cannot read coordinates from %s (%d. step - '%s'; %d. bead)\n\n", input_coor, count, stuff, test);
      test = '\0';
      exit(1);
    }
    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "Error: cannot read coordinates from %s (%d. step - '%s'; %d. bead)\n\n", input_coor, count, stuff, test);
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

    // read aggregates //{{{
    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      test = start - 1 + count; // total number of processed steps in agg file
      fprintf(stderr, "\nError: premature end of %s file (after %d. step - '%s')\n\n", input_agg, test, stuff);
      break;
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate memory for temporary density arrays //{{{
    double ***temp_rho = malloc(Counts.TypesOfBeads*sizeof(double **));
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      temp_rho[i] = malloc(aggs*sizeof(double *));
      for (int j = 0; j < aggs; j++) {
        temp_rho[i][j] = calloc(bins,sizeof(double));
      }
    } //}}}

    // calculate densities //{{{
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
        Vector com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, Bead, BeadType);

        // free temporary density array //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < aggs; k++) {
            for (int l = 0; l < bins; l++) {
              temp_rho[j][k][l] = 0;
            }
          }
        } //}}}

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

              temp_rho[Bead[Aggregate[i].Bead[j]].Type][correct_size][k]++;
            }
          }
        } //}}}

        // monomeric beads //{{{
        for (int j = 0; j < Counts.Unbonded; j++) {
          Vector dist = Distance(Bead[j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            temp_rho[Bead[j].Type][correct_size][k]++;
          }
        } //}}}

        agg_sizes[correct_size][1]++;

        // add from temporary density array to global density arrays //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < bins; k++) {
            rho[j][correct_size][k] += temp_rho[j][correct_size][k];
            rho_2[j][correct_size][k] += SQR(temp_rho[j][correct_size][k]);
          }
        } //}}}

        // count mol types in the aggregate
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          agg_mols[correct_size][Molecule[Aggregate[i].Molecule[j]].Type]++;
        }
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      for (int j = 0; j < aggs; j++) {
        free(temp_rho[i][j]);
      }
      free(temp_rho[i]);
    }
    free(temp_rho); //}}}

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
    // assemble correct name
    sprintf(str, "%s", output_rho);
    char str2[1030];
    sprintf(str2, "%s%d.rho", str, agg_sizes[i][0]);
    strcpy(str, str2);
    if ((out = fopen(str, "a")) == NULL) {
      ErrorFileOpen(str, 'a');
      exit(1);
    }

    // calculate rdf
    for (int j = 0; j < (bins-avg); j++) {

      // calculate volume of every shell that will be averaged
      double shell[avg];
      for (int k = 0; k < avg; k++) {
        shell[k] = 4 * PI * CUBE(width) * (CUBE(j+k+1) - CUBE(j+k)) / 3;
      }

      fprintf(out, "%.2f", width*(j+0.5*avg));

      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        double temp_rdp = 0, temp_number = 0,
               temp_rdp_err = 0, temp_number_err = 0;

        // sum up rdfs from all shells to be averaged
        for (int l = 0; l < avg; l++) {
          temp_rdp += rho[k][i][j+l] / (shell[l] * agg_sizes[i][1]);
          temp_rdp_err += rho_2[k][i][j+l] / (shell[l] * agg_sizes[i][1]);
          temp_number += rho[k][i][j+l] / agg_sizes[i][1];
          temp_number_err += rho_2[k][i][j+l] / agg_sizes[i][1];
        }

        temp_rdp_err = sqrt(temp_rdp_err - temp_rdp);
        temp_number_err = sqrt(temp_number_err - temp_number);

        // print average value to output file
//      fprintf(out, " %10f", temp_rdp/avg);
        fprintf(out, " %10f %10f", temp_rdp/avg, temp_rdp_err/avg);
        fprintf(out, " %10f %10f", temp_number/avg, temp_number_err/avg);
      }
      putc('\n',out);
    }

    fprintf(out, "# (1) Number of aggregates; Average numbers of molecules:");
    count = 2;
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      fprintf(out," (%d) %s", count++, MoleculeType[j].Name);
      if (i < (Counts.TypesOfMolecules-1)) {
        putc(';', out);
      }
    }
    putc('\n', out);
    putc('#', out);
    fprintf(out," %d", agg_sizes[i][1]);
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      fprintf(out," %lf", (double)(agg_mols[i][j])/agg_sizes[i][1]);
    }
    putc('\n', out);

    fclose(out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
    free(agg_mols[i]);
  }
  free(agg_sizes);
  free(agg_mols);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < aggs; j++) {
      free(rho[i][j]);
      free(rho_2[i][j]);
    }
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho);
  free(rho_2);
  free(specific_moltype_for_size); //}}}

  return 0;
}
