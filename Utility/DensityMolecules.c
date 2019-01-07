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
  fprintf(stderr, "   %s <input> <width> <output.rho> <mol name(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input>           input filename (either vcf or vtf format)\n");
  fprintf(stderr, "   <width>           width of a single bin\n");
  fprintf(stderr, "   <output.rho>      output density file (automatic ending 'molecule_name.rho' added)\n");
  fprintf(stderr, "   <mol name(s)>     molecule names to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined       specify that <input> contains joined coordinates\n");
  fprintf(stderr, "      -n <int>       number of bins to average\n");
  fprintf(stderr, "      -st <int>      starting timestep for calculation\n");
  fprintf(stderr, "      -c 'x's/<ints> use <int>-th molecule bead instead of centre of mass\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
DensityMolecules utility calculates number \
density for all bead types from the \
centre of mass (or specified bead number in a molecule) of specified molecules. \
\n\n");

/*      fprintf(stdout, "\
The utility uses traject.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");
*/

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input> <width> <output.rho> <mol name(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input>           input filename (either vcf or vtf format)\n");
      fprintf(stdout, "   <width>           width of a single bin\n");
      fprintf(stdout, "   <output.rho>      output density file (automatic ending 'molecule_name.rho' added)\n");
      fprintf(stdout, "   <mol name(s)>     molecule names to calculate density for\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      --joined       specify that <input> contains joined coordinates\n");
      fprintf(stdout, "      -n <int>       number of bins to average\n");
      fprintf(stdout, "      -st <int>      starting timestep for calculation\n");
      fprintf(stdout, "      -c 'x's/<ints> use <int>-th molecule bead instead of centre of mass\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int req_args = 4; //}}}

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
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-c") != 0) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than traject.vsf? //{{{
  char *input_vsf = calloc(32,sizeof(char *));
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

  // <input> - input coordinate file //{{{
  char input_coor[32];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  ext = 2;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(8*sizeof(char));
  }
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_coor, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <width> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.rho> - filename with bead densities //{{{
  char output_rho[16];
  strcpy(output_rho, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {

    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (mol_type == -1) {
      fprintf(stderr, "\nError: molecule '%s' does not exist in FIELD\n\n", argv[count]);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }

    // write initial stuff to output density file //{{{
    FILE *out;
    char str[128];

    sprintf(str, "%s%s.rho", output_rho, argv[count]);
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
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      fprintf(out, " %d: %s", 4*i+2, BeadType[i].Name);
    }
    fprintf(out, "\n# for each molecule type: rdp | stderr | rnp | stderr\n"); //}}}

    fclose(out); //}}}
  } //}}}

  // -c 'x's/<int(s)> option - specify which bead to use as a molecule centre //{{{

  // array for considering whether to use COM or specified bead number //{{{
  int centre[Counts.TypesOfMolecules];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    centre[i] = -1; // helper value
  } //}}}

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-c") == 0) {
      int j = 1; // count extra arguments
      while ((i+j) < argc && argv[i+j][0] != '-') {

        int mol_type = FindMoleculeType(argv[i+j], Counts, MoleculeType);

        if (mol_type == -1) {
          fprintf(stderr, "\nError: molecule '%s' does not exist in FIELD ('-c' option)\n\n", argv[i+j]);
          exit(1);
        } else {
          // Error - non-numeric argument //{{{
          if (argv[i+j+1][0] < '0' || argv[i+j+1][0] > '9') {
            ErrorNaN("-c");
            ErrorHelp(argv[0]);
            exit(1);
          } //}}}

          MoleculeType[mol_type].Use = true;

          centre[mol_type] = atoi(argv[i+j+1]) - 1; // bead indices start from 0

          // Error - too high bead number //{{{
          if (centre[mol_type] > MoleculeType[mol_type].nBeads) {
            fprintf(stderr, "Error: incorrect number in '-c' option (%dth bead in molecule "
              "%s containing only %d beads)\n", centre[mol_type]+1, MoleculeType[mol_type].Name,
              MoleculeType[mol_type].nBeads);
            exit(1);
          } //}}}
        }

        j += 2; // +2 because there's "<mol name> number"
      }
    }
  } //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[128];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_coor);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_coor);
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
    rho[i] = malloc(Counts.TypesOfMolecules*sizeof(double *));
    rho_2[i] = malloc(Counts.TypesOfMolecules*sizeof(double *));
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      rho[i][j] = calloc(bins,sizeof(double));
      rho_2[i][j] = calloc(bins,sizeof(double));
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
    VerboseOutput(verbose2, input_coor, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(stdout, " %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // skip first start-1 steps //{{{
  count = 0;
  int test;
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

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "Error: premature end of %s file\n\n", input_coor);
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
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
      }
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // join molecules if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate memory for temporary density arrays //{{{
    double ***temp_rho = malloc(Counts.TypesOfBeads*sizeof(double **));
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      temp_rho[i] = malloc(Counts.TypesOfMolecules*sizeof(double *));
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        temp_rho[i][j] = calloc(bins,sizeof(double));
      }
    } //}}}

    // calculate densities //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type_i = Molecule[i].Type;
      if (MoleculeType[mol_type_i].Use) {

        // determine centre to calculate densities from //{{{
        Vector com;
        if (centre[mol_type_i] == -1 ) { // use molecule's centre of mass
          com = CentreOfMass(MoleculeType[mol_type_i].nBeads, Molecule[i].Bead, Bead, BeadType);
        } else { // use centre[mol_type_i]-th molecule's bead as com
          com.x = Bead[Molecule[i].Bead[centre[mol_type_i]]].Position.x;
          com.y = Bead[Molecule[i].Bead[centre[mol_type_i]]].Position.y;
          com.z = Bead[Molecule[i].Bead[centre[mol_type_i]]].Position.z;
        } //}}}

        // free temporary density array //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < Counts.TypesOfMolecules; k++) {
            for (int l = 0; l < bins; l++) {
              temp_rho[j][k][l] = 0;
            }
          }
        } //}}}

        // molecule beads //{{{
        for (int j = 0; j < MoleculeType[mol_type_i].nBeads; j++) {
          int bead_j = Molecule[i].Bead[j];

          Vector dist = Distance(Bead[bead_j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            temp_rho[Bead[bead_j].Type][mol_type_i][k]++;
          }
        } //}}}

        // monomeric beads //{{{
        for (int j = 0; j < Counts.Unbonded; j++) {
          Vector dist = Distance(Bead[j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            temp_rho[Bead[j].Type][mol_type_i][k]++;
          }
        } //}}}

        // add from temporary density array to global density arrays //{{{
        for (int j = 0; j < Counts.TypesOfBeads; j++) {
          for (int k = 0; k < bins; k++) {
            rho[j][mol_type_i][k] += temp_rho[j][mol_type_i][k];
            rho_2[j][mol_type_i][k] += SQR(temp_rho[j][mol_type_i][k]);
          }
        } //}}}
      }
    } //}}}

    // free temporary density array //{{{
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
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

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  } //}}}

  // write densities to output file(s) //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      FILE *out;

      sprintf(str, "%s%s.rho", output_rho, MoleculeType[i].Name);
      if ((out = fopen(str, "a")) == NULL) {
        ErrorFileOpen(str, 'a');
        exit(1);
      }

      // calculate rdf
      for (int j = 1; j < (bins-avg); j++) {

        // calculate volume of every shell that will be averaged
        double shell[avg];
        for (int k = 0; k < avg; k++) {
          shell[k] = 4 * PI * CUBE(width) *(CUBE(j+k+1) - CUBE(j+k)) / 3;
        }

        fprintf(out, "%.2f", width*(j+0.5*avg));

        for (int k = 0; k < Counts.TypesOfBeads; k++) {
          double temp_rdp = 0, temp_number = 0,
                 temp_rdp_err = 0, temp_number_err = 0;

          // sum rdfs from all shells to be averaged
          for (int l = 0; l < avg; l++) {
            temp_rdp += rho[k][i][j+l] / (shell[l] * MoleculeType[i].Number * count);
            temp_rdp_err += rho_2[k][i][j+l] / (shell[l] * MoleculeType[i].Number * count);
            temp_number += rho[k][i][j+l] / MoleculeType[i].Number;
            temp_number_err += rho_2[k][i][j+l] / MoleculeType[i].Number;
          }

          temp_rdp_err = sqrt(temp_rdp_err - temp_rdp);
          temp_number_err = sqrt(temp_number_err - temp_number);

          // print average value to output file
          fprintf(out, " %10f %10f", temp_rdp/avg, temp_rdp_err/avg);
          fprintf(out, " %10f %10f", temp_number/avg, temp_number_err/avg);
        }
        putc('\n',out);
      }

      fclose(out);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      free(rho[i][j]);
      free(rho_2[i][j]);
    }
    free(rho[i]);
    free(rho_2[i]);
  }
  free(rho_2);
  free(rho); //}}}

  return 0;
}
