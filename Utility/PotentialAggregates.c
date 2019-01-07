#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input> <input.agg> <width> <output.txt> <agg sizes> <options>\n\n", cmd);

  fprintf(stderr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(stderr, "   <input.agg>       input agg file\n");
  fprintf(stderr, "   <width>           width of a single bin\n");
  fprintf(stderr, "   <output.txt>      output file (automatic ending '#.txt' added)\n");
  fprintf(stderr, "   <agg size(s)>     aggregate size(s) to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined       specify that <input> contains joined coordinates\n");
  fprintf(stderr, "      -st <int>      starting timestep for calculation\n");
  fprintf(stderr, "      -m <name(s)>   agg size means number of <name(s)> molecule types in an aggregate\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
/*      fprintf(stdout, "\
DensityAggregates utility calculates bead density for aggregates of given \
size(s) from \
their centre of mass. Beside unbonded beads it takes into account only beads \
from the current aggregate, not from any other aggregate. \
Care must be taken with beadtype names in molecule types, because if \
one beadtype appears in more molecule types, the resulting density for that \
beadtype will be averaged without regard for the various types of molecule it \
comes from (in that case, use -x option with SelectedVcf utility).\n\n");*/

/*      fprintf(stdout, "\
The utility uses traject.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n"); */

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input> <input.agg> <width> <output.txt> <agg size(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input.vcf>       input coordinate file (either vcf or vtf format)\n");
      fprintf(stdout, "   <input.agg>       input agg file\n");
      fprintf(stdout, "   <width>           width of a single bin\n");
      fprintf(stdout, "   <output.txt>      output file (automatic ending '#.txt' added)\n");
      fprintf(stdout, "   <agg size(s)>     aggregate size(s) to calculate density for\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      --joined       specify that <input> contains joined coordinates\n");
      fprintf(stdout, "      -st <int>      starting timestep for calculation\n");
      fprintf(stdout, "      -m <name(s)>   agg size means number of <name(s)> molecule types in an aggregate\n");
      CommonHelp(0);
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
    ErrorHelp(argv[0]);
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
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-m") != 0) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // elstat parameters
  double bjerrum = 1.1,
         lambda = 0.2,
         r_c = 3.0,
         beta = (5 * r_c) / (8 * lambda),
         images = 5;
  int point_density = 2;
  int max_points = 100;

  // options before reading system data //{{{
  // use .vsf file other than traject.vsf? //{{{
  char *input_vsf = calloc(32,sizeof(char *));
  if (FileOption(argc, argv, "-i", &input_vsf)) {
    exit(1);
  }
  if (input_vsf[0] == '\0') {
    strcpy(input_vsf, "traject.vsf");
  }

  // test if structure file ends with '.vsf' or '.vtf'
  int ext = 2;
  char **extension;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_vsf, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(32,sizeof(char *));
  if (FileOption(argc, argv, "-b", &bonds_file)) {
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
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_vcf[32];
  strcpy(input_vcf, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  ext = 2;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_vcf, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <input.agg> - filename of input file with aggregate information //{{{
  char input_agg[32];
  strcpy(input_agg, argv[++count]);

  // test if <input.agg> ends with '.agg'
  ext = 1;
  extension = malloc(ext*sizeof(char *));
  extension[0] = malloc(5*sizeof(char));
  strcpy(extension[0], ".agg");
  if (!ErrorExtension(input_agg, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <width> - width of single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.txt> - filename with bead densities //{{{
  char output_elstat[32];
  strcpy(output_elstat, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information{{{
  bool indexed = ReadStructure(input_vsf, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf); //}}}

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

  // <agg sizes> - aggregate sizes for calculation //{{{
  int **agg_sizes = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    agg_sizes[i] = calloc(2,sizeof(int));
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
    char str[128];
    strcpy(str, output_elstat);

    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (specific_moltype_for_size[i]) {
        char str2[256];
        sprintf(str2, "%s%s", str, MoleculeType[i].Name);
        strcpy(str, str2);
      }
    }
    char str2[256];
    sprintf(str2, "%s_%d.txt", str, agg_sizes[aggs][0]);
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

    char name[10];
    fscanf(agg, "%s", name);
    int type = FindBeadType(name, Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\nError: bead type '%s' is not in %s file\n\n", name, input_vcf);
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
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    ErrorFileOpen(input_vcf, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[128];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_vcf);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "\nError: cannot read pbc from %s\n\n", input_vcf);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // number of bins
  int bins = ceil(Min3(BoxLength.x, BoxLength.y, BoxLength.z) / (3 * width));

  double max_dist = (0.5 + images) * Min3(BoxLength.x, BoxLength.y, BoxLength.z);

  // allocate memory for density arrays //{{{
//double **elstat = malloc(aggs*sizeof(double *));
  double **elstat_potential = malloc(aggs*sizeof(double *));
  double **elstat_potential_sqr = malloc(aggs*sizeof(double *));
  for (int i = 0; i < aggs; i++) {
//  elstat[i] = calloc(bins,sizeof(double));
    elstat_potential[i] = calloc(bins,sizeof(double));
    elstat_potential_sqr[i] = calloc(bins,sizeof(double));
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

    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      fprintf(stderr, "\nError: premature end of %s file (after %d. step - '%s')\n\n", input_agg, count, stuff);
      test = '\0';
      break;
    }
    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\nError: premature end of %s file (%d. step - '%s')\n\n", input_vcf, --count, stuff);
      exit(1);
    }
  }
  if (test == '\0') {
    exit(1);
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

  // create array of bead indices
  int *index = calloc(Counts.BeadsInVsf,sizeof(int));
  for (int i = 0; i < Counts.Beads; i++) {
    index[Bead[i].Index] = i;
  }

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
    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      test = start - 1 + count; // total number of processed steps in agg file
      fprintf(stderr, "\nError: premature end of %s file (after %d. step - '%s')\n", input_agg, test, stuff);
      break;
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // calculate potential //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {

      // allocate memory for temporary elstat array
      double *temp_elstat = calloc(bins,sizeof(double));

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
//      printf("\nAggregate[%d]\n", i);

        // zeroize temporary array //{{{
        for (int j = 0; j < bins; j++) {
          temp_elstat[j] = 0;
        } //}}}

        Vector com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, Bead, BeadType);

        // move beads so that com is in the box's centre
        for (int j = 0; j < Counts.Beads; j++) {
          Bead[j].Position.x -= com.x;
          Bead[j].Position.y -= com.y;
          Bead[j].Position.z -= com.z;
        }

        for (int j = 1; j < bins; j++) {
          double r = width * j;

          // start calculations from distance 0.1
          if (r > 0.1) {
            int N = point_density * 4 * PI * SQR(r);
            if (N > max_points) {
              N = max_points;
            }
//          N = 1000;

            int N_count = 0; // number of points

            double a = 4 * PI / N;
            double d = sqrt(a);
            int M1 = round(PI / d);
            double d1 = PI / M1;
            double d2 = a / d1;

            for (int mm = 0; mm < M1; mm++) {
              double angle1 = PI * (mm + 0.5) / M1;
              int M2 = round(2 * PI * sin(angle1) / d2);
              for (int nn = 0; nn < M2; nn++) {
                double angle2 = 2 * PI * nn / M2;
                Vector point;
                point.x = r * sin(angle1) * cos(angle2);
                point.y = r * sin(angle1) * sin(angle2);
                point.z = r * cos(angle1);
                N_count++;

                for (int l = 0; l < Counts.Beads; l++) {
                  if (BeadType[Bead[l].Type].Charge != 0) {
                    Vector dist;
                    dist = Distance(Bead[l].Position, point, BoxLength);
                    double rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                    double coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                    if (rij < r_c) { // short ranged part
                      temp_elstat[j] += coulomb * (1 - (1 + beta) * exp(-2 * beta * rij));
                    } else { // long ranged part
                      temp_elstat[j] += coulomb;
                    }
                    if (max_dist < 0);

                    // periodic images //{{{
                    Vector dist_orig;
                    dist_orig.x = dist.x;
                    dist_orig.y = dist.y;
                    dist_orig.z = dist.z;
                    for (int m = 1; m <= images; m++) {
                      // +mz
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +0y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -1y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -1y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -1y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z + m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0z
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -1y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -1y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -1y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mz
                      // -mx +my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y + m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx +0y //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x +0y //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx +0y //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // -mx -my //{{{
                      dist.x = dist_orig.x - m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +0x -my //{{{
                      dist.x = dist_orig.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                      // +mx -my //{{{
                      dist.x = dist_orig.x + m * BoxLength.x;
                      dist.y = dist_orig.y - m * BoxLength.y;
                      dist.z = dist_orig.z - m * BoxLength.z;
                      rij = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
                      coulomb = bjerrum * BeadType[Bead[l].Type].Charge / rij;
                      if (rij < max_dist) {
                        temp_elstat[j] += coulomb;
                      } //}}}
                    } //}}}
                  }
                }
              }
            }

            // add from temporary elstat array to global elstat array
            elstat_potential[correct_size][j] += temp_elstat[j] / N_count;
            elstat_potential_sqr[correct_size][j] += SQR(temp_elstat[j]) / N_count; // for error calculation
          }
        }

        agg_sizes[correct_size][1]++;
      }

      // free temporary elstat array
      free(temp_elstat);
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

  // write elstat to output file(s) //{{{
  for (int i = 0; i < aggs; i++) {
    FILE *out;
    // assemble correct name
    sprintf(str, "%s", output_elstat);
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (specific_moltype_for_size[j]) {
        char str2[256];
        sprintf(str2, "%s%s", str, MoleculeType[j].Name);
        strcpy(str, str2);
      }
    }
    char str2[256];
    sprintf(str2, "%s_%d.txt", str, agg_sizes[i][0]);
    strcpy(str, str2);
    if ((out = fopen(str, "a")) == NULL) {
      ErrorFileOpen(str, 'a');
      exit(1);
    }

    // average elstat
    for (int j = 0; j < bins; j++) {
      if ((width*j) > 0.1) {
        // print distance from aggregate com
        fprintf(out, "%.2f", width*(j+0.5));

        // print average value and standard deviation to output file
        fprintf(out, " %20.15f", elstat_potential[i][j]/agg_sizes[i][1]);
        putc('\n',out);
      }
    }

    fprintf(out, "# %d molecules in aggregate (sum of", agg_sizes[i][0]);
    int types = 0;
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (specific_moltype_for_size[j]) {
        if (++types == 1) { // first mol type
          fprintf(out, " %s", MoleculeType[j].Name);
        } else { // second and higher mol type
          fprintf(out, ", %s", MoleculeType[j].Name);
        }
      }
    }
    fprintf(out, " - %d aggregates)\n", agg_sizes[i][1]);

    fclose(out);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(agg_sizes[i]);
  }
  free(agg_sizes);
  for (int i = 0; i < aggs; i++) {
//  free(elstat[i]);
    free(elstat_potential[i]);
    free(elstat_potential_sqr[i]);
  }
//free(elstat);
  free(elstat_potential_sqr);
  free(specific_moltype_for_size);
  free(index); //}}}

  return 0;
}
