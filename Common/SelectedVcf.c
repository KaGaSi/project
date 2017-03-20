#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> ", cmd);
  fprintf(stderr, "<output.vcf> <type names> <options>\n\n");

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <output.vcf>      output filename (vcf format)\n");
  fprintf(stderr, "   <type names>      names of bead types to save\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -j             join molecules (remove pbc)\n");
  fprintf(stderr, "      -st <start>    number of timestep to start from\n");
  fprintf(stderr, "      -sk <skip>     leave out every 'skip' steps\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("SelectedVcf creates new <output.vcf> file from <input.vcf> containing only  \n");
      printf("selected bead types. Also <start> timesteps can be omitted and every <skip> \n");
      printf("timestep can be left out.                                                   \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> ", argv[0]);
      printf("<output.vcf> <type names> <options>\n\n");

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <output.vcf>      output filename (vcf format)\n");
      printf("   <type names>      names of bead types to save\n");
      printf("   <options>\n");
      printf("      -j             join molecules (remove pbc)\n");
      printf("      -st <start>    number of timestep to start from\n");
      printf("      -sk <skip>     leave out every 'skip' steps\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 3; //}}}

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
        strcmp(argv[i], "-j") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-sk") != 0) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // standard options //{{{
  char *vsf_file = calloc(32,sizeof(char *));
  char *bonds_file = calloc(32,sizeof(char *));
  bool verbose, verbose2, silent, script;
  bool error = CommonOptions(argc, argv, &vsf_file, &bonds_file, &verbose, &verbose2, &silent, &script);

  // was there error during CommonOptions()?
  if (error) {
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -j option - join molecules //{{{
  bool join = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      join = true;

      break;
    }
  } //}}}

  // -st <start> - number of starting timestep //{{{
  int start = 1;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-st") == 0) {

      // Error - non-numeric or missing argument //{{{
      if ((i+1) >= argc) {
        fprintf(stderr, "Missing numeric argument for '-st' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }
      if ((i+1) >= argc || argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '-st' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      } //}}}

      start = atoi(argv[i+1]);
    }
  } //}}}

  // -sk <skip> - number of steps to skip per one used //{{{
  int skip = 0;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-sk") == 0) {

      // Error - non-numeric or missing argument //{{{
      if ((i+1) >= argc) {
        fprintf(stderr, "Missing numeric argument for '-sk' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }
      if (argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '-sk' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      } //}}}

      skip = atoi(argv[i+1]);
    }
  } //}}}

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

  // <output.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[32];
  strcpy(output_vcf, argv[++count]);

  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  dot = strrchr(output_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<output.vcf> '%s' does not have .vcf ending!\n", output_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

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

  // <type names> - names of bead types to save //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);

    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", argv[count], input_vcf);
      exit(1);
    }

    BeadType[type].Write = true;
  }

//// Error - does not make sense to use all bead types
//if ((count-5) == Counts.TypesOfBeads) {
//  fprintf(stderr, "All beadtypes are selected which would only copy %s,\n", input_vcf);
//  fprintf(stderr, "therefore this is not allowed!\n\n");
//  ErrorHelp(argv[0]);
//  exit(1);
//} //}}}

  // print selected bead type names to output .vcf file //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_vcf);
    exit(1);
  }

  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Write) {
      fprintf(out, "# %s\n", BeadType[i].Name);
    }
  }

  fclose(out); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    printf("\n   Starting from %d. timestep\n", start);
    printf("   Every %d. timestep used\n", skip+1);
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

  // print pbc to output .vcf file //{{{
  if ((out = fopen(output_vcf, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_vcf);
    exit(1);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  fclose(out); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{ //{{{
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
        printf("Discarding step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rDiscarding step: %6d", count);
      }
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff) == 1) {
      fprintf(stderr, "Premature end of %s file!\n", input_vcf);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      printf("Discarded steps: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rDiscarded steps: %6d\n", count);
    }
  } //}}}
  //}}} //}}}

  // main loop //{{{
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

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

    count++;
    if (!silent) {
      if (script) {
        printf("Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rStep: %6d", count);
      }
    }

    // join molecules? //{{{
    if (join) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } else { // if rounding leads to BoxLength, move it bead to other side of box
      for (int i = 0; i < (Counts.Bonded+Counts.Unbonded); i++) {
        char check[8];
        char box[8];
        // x direction
        sprintf(check, "%.3f", Bead[i].Position.x);
        sprintf(box, "%.3f", BoxLength.x);
        if (strcmp(check, box) == 0) {
          Bead[i].Position.x = 0;
        }
        // y direction
        sprintf(check, "%.3f", Bead[i].Position.y);
        sprintf(box, "%.3f", BoxLength.y);
        if (strcmp(check, box) == 0) {
          Bead[i].Position.y = 0;
        }
        // z direction
        sprintf(check, "%.3f", Bead[i].Position.z);
        sprintf(box, "%.3f", BoxLength.z);
        if (strcmp(check, box) == 0) {
          Bead[i].Position.z = 0;
        }
      }
    } //}}}

    // open output .vcf file for appending //{{{
    if ((out = fopen(output_vcf, "a")) == NULL) {
      // print newline to stdout if Step... doesn't end with one
      if (!script && !silent) {
        putchar('\n');
      }
      fprintf(stderr, "Cannot open file %s!\n", output_vcf);
      exit(1);
    } //}}}

    WriteCoorIndexed(out, Counts, BeadType, Bead, stuff);

    fclose(out);

    // skip every 'skip' steps //{{{
    for (int i = 0; i < skip; i++) {
      // test whether at vcf's eof //{{{
      if ((test = getc(vcf)) == EOF) {
        break;
      }
      ungetc(test, vcf); //}}}

      if (!silent) {
        fflush(stdout);
        if (script) {
          printf("Step: %6d\n", count);
        } else {
          printf("\rStep: %6d", ++count);
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

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  //}}}

  return 0;
}
