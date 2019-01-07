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

  fprintf(stderr, "   <input>                 input filename (vcf format)\n");
  fprintf(stderr, "   <width>                 width of a single bin\n");
  fprintf(stderr, "   <output>                output file with distribution of bond lengths\n");
  fprintf(stderr, "   <mole name(s)>          names of molecule type(s) to use for calculation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -st <int>            starting timestep for calculation\n");
  fprintf(stderr, "      -d <name> <int(s)>   calculate distribution of distances between specified bead indices\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
BondLength utility calculates distribution of bond lengths for all bead type \
pairs in specified molecule(s).\n\n");

/*      fprintf(stdout, "\
The utility uses traject.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");
*/

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input> <width> <output> <mol name(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input>                 input filename (vcf format)\n");
      fprintf(stdout, "   <width>                 width of a single bin\n");
      fprintf(stdout, "   <output>                output file with distribution of bond length\n");
      fprintf(stdout, "   <mol name(s)>           names of molecule type(s) to use for calculation\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -st <int>            starting timestep for calculation\n");
      fprintf(stdout, "      -d <name> <int(s)>   calculate distribution of distances between specified bead indices\n");
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

  if (argc < req_args) {
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
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-d") != 0) {

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

  // <input.vcf> - filename of input vcf file (must end with .vcf) //{{{
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

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<width>");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output> - file name with bond length distribution //{{{
  char output[32];
  strcpy(output, argv[++count]); //}}}

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

  // <molecule names> - names of molecule types to use //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (type == -1) {
      fprintf(stderr, "Error: molecule type '%s' is not in %s file\n\n", argv[count], input_coor);
      exit(1);
    }

    MoleculeType[type].Use = true;
  } //}}}

  // '-d' option - specify for bead ids to calculate distance between //{{{
  int bead[100] = {0}, number_of_beads = 0;
  char output_d[32];
  output_d[0] = '\0';
  if (FileIntsOption(argc, argv, "-d", bead, &number_of_beads, output_d)) {
    exit(1);
  }

  // if '-d' is present, but without numbers - use first and last for each molecule
  if (output_d[0] != '\0' && number_of_beads == 0) {
    number_of_beads = 2;
  }

  // Error: wrong number of integers //{{{
  if (output_d[0] != '\0' && (number_of_beads%2) != 0) {
    fprintf(stderr, "\nError: '-d' option - number of bead ids must be even\n");
    exit(1);
  } //}}}

  // Error: same bead ids //{{{
  for (int i = 0; i < number_of_beads; i += 2) {
    if (bead[i] == bead[i+1] || bead[i] == 0 || bead[i+1] == 0) {
      fprintf(stderr, "\nError: '-d' option - the bead indices must be non-zero and different\n");
      exit(1);
    }
  } //}}}

  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero (or -1 if none specified)
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_coor, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "Error: cannot read a string from '%s' file\n\n", input_coor);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Error: cannot read pbc from %s\n\n", input_coor);
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
  int bins = BoxLength.x / (2 * width);

  // arrays for distributions //{{{
  double *distance[Counts.TypesOfMolecules][number_of_beads/2];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/2); j++) {
      distance[i][j] = calloc(bins, sizeof(double));
    }
  }
  double *length[Counts.TypesOfMolecules][Counts.TypesOfBeads][Counts.TypesOfBeads];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        length[i][j][k] = calloc(bins,sizeof(double));
      }
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = calloc(128,sizeof(int)); //}}}

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
  count = 0;
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

    // join all molecules
    RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

    // calculate bond length //{{{

    // go through all molecules
    for (int i = 0; i < Counts.Molecules; i++) {
      int type = Molecule[i].Type;
      if (MoleculeType[type].Use) { // use only specified molecule types
        for (int j = 0; j < MoleculeType[type].nBonds; j++) {

          // bead ids in the bond //{{{
          int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
          int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]]; //}}}

          // bond length //{{{
          Vector bond;
          bond.x = Bead[id1].Position.x - Bead[id2].Position.x;
          bond.y = Bead[id1].Position.y - Bead[id2].Position.y;
          bond.z = Bead[id1].Position.z - Bead[id2].Position.z;

          bond.x = sqrt(SQR(bond.x) + SQR(bond.y) + SQR(bond.z)); //}}}

          int k = bond.x / width;
          if (k < bins) {
            length[type][Bead[id1].Type][Bead[id2].Type][k]++;
          }
        }
      }
    } //}}}

    // calculate distance (-d option) //{{{
    if (output_d[0] != '\0') {
      for (int i = 0; i < Counts.Molecules; i++) {
        int type = Molecule[i].Type;
        if (MoleculeType[type].Use) { // use only specified molecule types
          for (int j = 0; j < number_of_beads; j += 2) {

            // bead ids the distance //{{{
            int id1, id2;
            // use first molecule bead if bead index too high or -1
            if (bead[j] == -1 || bead[j] >= MoleculeType[type].nBeads) {
              id1 = Molecule[i].Bead[0];
            } else{ // use specified index otherwise
              id1 = Molecule[i].Bead[bead[j]];
            }
            // use last molecule bead if bead index too high or -1
            if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[type].nBeads) {
              id2 = Molecule[i].Bead[MoleculeType[type].nBeads-1];
            } else{ // use specified index otherwise
              id2 = Molecule[i].Bead[bead[j+1]];
            } //}}}

            // distance //{{{
            Vector dist;
            dist.x = Bead[id1].Position.x - Bead[id2].Position.x;
            dist.y = Bead[id1].Position.y - Bead[id2].Position.y;
            dist.z = Bead[id1].Position.z - Bead[id2].Position.z;

            dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z)); //}}}

            int k = dist.x / width;
            if (k < bins) {
              distance[type][j/2][k]++;
            }
          }
        }
      }
    } //}}}

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      fprintf(stdout, "\n%s", stuff);
  }

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  }

  fclose(vcf); //}}}

  // count total number of bonds in molecules //{{{
  int bonds[Counts.TypesOfMolecules][Counts.TypesOfBeads][Counts.TypesOfBeads]; //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] = 0;
        }
      }
    }
  } //}}}

  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = j; k < Counts.TypesOfBeads; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] += length[i][j][k][l] + length[i][k][j][l];
        }
      }
    }
  } //}}}

  // write distribution of bond lengths //{{{
  // open output file for appending //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    ErrorFileOpen(output, 'w');
    exit(1);
  } //}}}

  // print first line of output file - molecule names and beadtype pairs //{{{
  fprintf(out, "#distance ");

  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, "%s: ", MoleculeType[i].Name);

      for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
        for (int k = j; k < MoleculeType[i].nBTypes; k++) {
          fprintf(out, " %s-%s", BeadType[MoleculeType[i].BType[j]].Name, BeadType[MoleculeType[i].BType[k]].Name);
        }
      }
      putc(';', out);
    }
  }
  putc('\n', out); //}}}

  // write to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(out, "%7.4f", width*(2*i+1)/2);

    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use) {

        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < MoleculeType[j].nBTypes; k++) {
          for (int l = k; l < MoleculeType[j].nBTypes; l++) {

            int btype1 = MoleculeType[j].BType[k];
            int btype2 = MoleculeType[j].BType[l];

            // btype1 must be lower then btype2 - TODO: check why
            if (btype1 > btype2) {
              int swap = btype1;
              btype1 = btype2;
              btype2 = swap;
            }

            fprintf(out, "%10f", (double)(length[j][btype1][btype2][i]+length[j][btype2][btype1][i])/bonds[j][btype1][btype2]);
          }
        }
      }
    }
    putc('\n', out);
  }
  fclose(out); //}}}
  //}}}

  // write distribution of distances from '-d' option //{{{
  if (output_d[0] != '\0') {
    // open output file for appending //{{{
    FILE *out;
    if ((out = fopen(output_d, "w")) == NULL) {
      ErrorFileOpen(output_d, 'w');
      exit(1);
    } //}}}

    // print first line of output file - molecule names and indices //{{{
    fprintf(out, "# distance ");

    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(out, " %s: ", MoleculeType[i].Name);

        for (int j = 0; j < number_of_beads; j += 2) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (bead[j] == -1 || bead[j] >= MoleculeType[i].nBeads) {
            id1 = 1;
          } else{ // use specified index otherwise
            id1 = bead[j]+1;
          }
          // use last molecule bead if bead index too high or -1
          if (bead[j+1] == -1 || bead[j+1] >= MoleculeType[i].nBeads) {
            id2 = MoleculeType[i].nBeads;
          } else{ // use specified index otherwise
            id2 = bead[j+1]+1;
          } //}}}

          fprintf(out, " %d-%d", id1, id2);
        }
        putc(';', out);
      }
    }
    putc('\n', out); //}}}

    // write to output file //{{{
    for (int i = 0; i < bins; i++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2);

      for (int j = 0; j < Counts.TypesOfMolecules; j++) {
        if (MoleculeType[j].Use) {
          for (int k = 0; k < number_of_beads; k += 2) {
            fprintf(out, " %10f", (double)(distance[j][k/2][i])/(count*MoleculeType[j].Number));
          }
        }
      }
      putc('\n', out);
    }
    fclose(out); //}}}
  }
  //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < Counts.TypesOfBeads; j++) {
      for (int k = 0; k < Counts.TypesOfBeads; k++) {
        free(length[i][j][k]);
      }
    }
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (number_of_beads/2); j++) {
      free(distance[i][j]);
    }
  }
  //}}}

  return 0;
}
