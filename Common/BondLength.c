#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <width> <output file> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <width>           width of a single bin\n");
  fprintf(stderr, "   <output file>     name of output file with end-to-end distances\n");
  fprintf(stderr, "   <molecule names>  names of molecule type(s) to use for calculation\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("BondLength utility calculates distribution of bond lengths for all bead type\n");
      printf("pairs in specified molecule(s).                                             \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <width> <output file> <molecule names> <options>\n\n", argv[0]);

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <width>           width of a single bin\n");
      printf("   <output file>     name of output file with end-to-end distances\n");
      printf("   <molecule names>  names of molecule type(s) to use for calculation\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 4; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (argc < options) {
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
        strcmp(argv[i], "--script") != 0) {

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
  bool verbose, verbose2, silent;
  SilentOption(argc, argv, &verbose, &verbose2, &silent); // no output
  VerboseShortOption(argc, argv, &verbose); // verbose output
  VerboseLongOption(argc, argv, &verbose, &verbose2); // more verbose output
  bool script = BoolOption(argc, argv, "--script"); // do not use \r & co.
  // }}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      printf(" %s", argv[i]);
    printf("\n\n");
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

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <width>!\n");
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
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file);

  // <molecule names> - names of molecule types to use //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (type == -1) {
      fprintf(stderr, "Molecule type '%s' is not in %s coordinate file!\n", argv[count], input_vcf);
      exit(1);
    }

    MoleculeType[type].Use = true;
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

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
    printf("   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // number of bins
  int bins = BoxLength.x / (2 * width);

  // array for average bond length //{{{
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

  // main loop //{{{
  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      if (script) {
        printf("Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rStep: %6d", count);
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

    // join all molecules
    RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

    // calculate bond distance //{{{

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

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      printf("\n%s", stuff);
  }

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
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

  // open output file for appending //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output);
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
  //}}}

  return 0;
}
