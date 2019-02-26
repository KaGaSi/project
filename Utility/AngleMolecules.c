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
  fprintf(stderr, "   %s <input> <width> <output> <mol name(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input>         input coordinate file (either vcf or vtf format)\n");
  fprintf(stderr, "   <width>         width of a single bin in degrees\n");
  fprintf(stderr, "   <output>        output file with distribution of dihedral angles\n");
  fprintf(stderr, "   <mol name(s)>   molecule name(s) to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined     specify that <input> contains joined coordinates\n");
  fprintf(stderr, "      -n <ints>    bead indices (multiple of 3 <ints>) for dihedral calculation (default: 1 2 3)\n");
  fprintf(stderr, "      -a <name>    write angle of all molecules in all times to <name>\n");
  fprintf(stderr, "      -st <int>    starting timestep for calculation\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
AngleMolecules utility calculates angles \
for specified beads (three beads per angel) in each molecule of specified \
molecule type(s). The specified beads do not to be connected by bonds. \
\n\n");

/*      fprintf(stdout, "\
The utility uses traject.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");
*/

      fprintf(stdout, "Usage:\n");

      fprintf(stdout, "   %s <input> <width> <output> <mol name(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input>         input coordinate file (either vcf or vtf format)\n");
      fprintf(stdout, "   <width>         width of a single bin in degrees\n");
      fprintf(stdout, "   <output>        output file with distribution of dihedral angles\n");
      fprintf(stdout, "   <mol name(s)>   molecule name(s) to calculate density for\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -n <ints>    bead indices (multiple of 3 <ints>) for dihedral calculation (default: 1 2 3)\n");
      fprintf(stdout, "      -a <name>    write angle of all molecules in all times to <name>\n");
      fprintf(stdout, "      -st <int>    starting timestep for calculation\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int req_args = 2; //}}}

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
        strcmp(argv[i], "-a") != 0 &&
        strcmp(argv[i], "-st") != 0) {

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
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - filename of input vcf file (must end with .vcf) //{{{
  char input_coor[32];
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD) //{{{
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
  //}}}

  // <width> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]);

  // number of bins between 0 and 180 deg
  int bins = ceil(180 / width); //}}}

  // <output> - file name with dihedral angle distribution //{{{
  char output_distr[32];
  strcpy(output_distr, argv[++count]); //}}}

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
      fprintf(stderr, "Error: molecule '%s' does not exist in FIELD\n\n", argv[count]);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }
  } //}}}

  // '-n' option - specify bead ids //{{{
  int bead[100] = {0}, number_of_beads = 3, beads_per_angle = 3, test;
  bead[0] = 1; // default ids for angle
  bead[1] = 2;
  bead[2] = 3;
  if (MultiIntegerOption(argc, argv, "-n", &test, bead)) {
    exit(1);
  }
  if (test != 0) { // -a is present
    number_of_beads = test;
  }

  // Error: wrong number of integers //{{{
  if ((number_of_beads%beads_per_angle) != 0) {
    fprintf(stderr, "\nError: '-a' option - number of bead ids must be dividable by three.\n");
    exit(1);
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero

    // Error - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && bead[i] >= MoleculeType[j].nBeads) {
        fprintf(stderr, "\nError: '-a' option - %d is larger than the number of beads in molecule %s\n\n", bead[i], MoleculeType[j].Name);
        ErrorHelp(argv[0]);
        exit(1);
      }
    } //}}}
  } //}}}

  // '-a' option - write angles for all molecules //{{{
  char *output = calloc(32,sizeof(char *));
  if (FileOption(argc, argv, "-a", &output)) {
    exit(1);
  }

  // write initial stuff to output if '-a' is used
  if (output[0] != '\0') {
    // open output //{{{
    FILE *out;
    if ((out = fopen(output, "w")) == NULL) {
      ErrorFileOpen(output, 'w');
      exit(1);
    } //}}}

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++)
      fprintf(out, " %s", argv[i]);
    putc('\n', out); //}}}

    // print molecule names & bead ids //{{{
    fprintf(out, "# angles between beads:");
    for (int j = 0; j < number_of_beads; j += beads_per_angle) {
      fprintf(out, " (%d) %d-%d-%d;", j/beads_per_angle+1, bead[j]+1, bead[j+1]+1, bead[j+2]+1);
    }
    putc('\n', out);
    fprintf(out, "# columns: (1) step;");
    int j = 2;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        if ((number_of_beads/beads_per_angle*MoleculeType[i].Number) == 1) {
          fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
        } else {
          fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_beads/beads_per_angle*MoleculeType[i].Number-1, MoleculeType[i].Name);
        }
        j += number_of_beads / beads_per_angle * MoleculeType[i].Number;
      }
    }
    putc('\n', out); //}}}

    fclose(out);
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

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
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

  // allocate array for distribution of angles //{{{
  double avg_angle[Counts.TypesOfMolecules][number_of_beads/beads_per_angle];
  double *distr[Counts.TypesOfMolecules][number_of_beads/beads_per_angle];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
      avg_angle[i][j] = 0;
      distr[i][j] = calloc(bins, sizeof(double));
    }
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

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\nError: premature end of %s file\n\n", input_coor);
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

    // calculate angles //{{{
    double angle[Counts.Molecules][number_of_beads/beads_per_angle];
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type_i = Molecule[i].Type;
      if (MoleculeType[mol_type_i].Use) {
        for (int j = 0; j < number_of_beads; j += beads_per_angle) {
          Vector u, v;
          // first vector
          u.x = Bead[Molecule[i].Bead[bead[j+0]]].Position.x - Bead[Molecule[i].Bead[bead[j+1]]].Position.x;
          u.y = Bead[Molecule[i].Bead[bead[j+0]]].Position.y - Bead[Molecule[i].Bead[bead[j+1]]].Position.y;
          u.z = Bead[Molecule[i].Bead[bead[j+0]]].Position.z - Bead[Molecule[i].Bead[bead[j+1]]].Position.z;
          // second vector
          v.x = Bead[Molecule[i].Bead[bead[j+2]]].Position.x - Bead[Molecule[i].Bead[bead[j+1]]].Position.x;
          v.y = Bead[Molecule[i].Bead[bead[j+2]]].Position.y - Bead[Molecule[i].Bead[bead[j+1]]].Position.y;
          v.z = Bead[Molecule[i].Bead[bead[j+2]]].Position.z - Bead[Molecule[i].Bead[bead[j+1]]].Position.z;
          // calculate angle between the two vectors
          double size[2];
          size[0] = sqrt(SQR(u.x) + SQR(u.y) + SQR(u.z));
          size[1] = sqrt(SQR(v.x) + SQR(v.y) + SQR(v.z));
          double scalar = u.x * v.x + u.y * v.y + u.z * v.z;
          angle[i][j/beads_per_angle] = acos(scalar / (size[0] * size[1])); // in rad
          angle[i][j/beads_per_angle] *= 180 / PI; // in degrees

          // add to average
          avg_angle[mol_type_i][j/beads_per_angle] += angle[i][j/beads_per_angle];

          int k = angle[i][j/beads_per_angle] / width;
          if (k < bins) {
            distr[mol_type_i][j/beads_per_angle][k]++;
          } else {
            fprintf(stdout, "\nWarning - weird angle: %lf degrees\n", angle[i][j/beads_per_angle]);
          }
        }
      }
    } //}}}

    // write all angles to output if '-a' is used //{{{
    if (output[0] != '\0') {
      FILE *out;
      if ((out = fopen(output, "a")) == NULL) {
        ErrorFileOpen(output, 'a');
        exit(1);
      }

      fprintf(out, "%6d", count);
      for (int i = 0; i < Counts.Molecules; i++) {
        int mol_type_i = Molecule[i].Type;
        if (MoleculeType[mol_type_i].Use) {
          for (int j = 0; j < number_of_beads; j += beads_per_angle){
            fprintf(out, " %10.6f", angle[i][j/beads_per_angle]);
          }
        }
      }
      putc('\n', out);

      fclose(out);
    } //}}}

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

  // write distribution of angles //{{{
  // open output file for appending //{{{
  FILE *out;
  if ((out = fopen(output_distr, "w")) == NULL) {
    ErrorFileOpen(output_distr, 'w');
    exit(1);
  } //}}}

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print molecule names & bead ids //{{{
  fprintf(out, "# angles between beads:");
  for (int j = 0; j < number_of_beads; j += beads_per_angle) {
    fprintf(out, " (%d) %d-%d-%d;", j/beads_per_angle+1, bead[j], bead[j+1], bead[j+2]);
  }
  putc('\n', out);
  fprintf(out, "# columns: (1) angle [deg];");
  int j = 2;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if ((number_of_beads/beads_per_angle) == 1) {
        fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_beads/beads_per_angle-1, MoleculeType[i].Name);
      }
      j += number_of_beads / beads_per_angle;
    }
  }
  putc('\n', out); //}}}

  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%5.1f", width*(2*i+1)/2);

    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {

        for (int k = 0; k < (number_of_beads/beads_per_angle); k++) {
          fprintf(out, "%10f", (double)(distr[j][k][i])/(count*MoleculeType[j].Number));
        }
      }
    }
    putc('\n', out);
  } //}}}

  // write to output average angles //{{{
  fprintf(out, "# simple averages:");
  j = 1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      if ((number_of_beads/beads_per_angle) == 1) {
        fprintf(out, " (%d) %s molecules;", j, MoleculeType[i].Name);
      } else {
        fprintf(out, " (%d) to (%d) %s molecules;", j, j+number_of_beads/beads_per_angle-1, MoleculeType[i].Name);
      }
      j += number_of_beads / beads_per_angle;
    }
  }
  fprintf(out, "\n#");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
        fprintf(out, " %7.3f", avg_angle[i][j]/(count*MoleculeType[i].Number));
      }
    }
  }
  putc('\n', out); //}}}

  fclose(out);
  //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(output);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/beads_per_angle); j++) {
      free(distr[i][j]);
    }
  } //}}}

  return 0;
}