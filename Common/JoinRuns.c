#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <1st input.vcf> <2nd input.vcf> <2nd input.vsf> ", cmd);
  fprintf(stderr, "<output.vcf> <type names> <options>\n\n");

  fprintf(stderr, "   <1st input.vcf>   input filename of 1st run (vcf format)\n");
  fprintf(stderr, "   <2nd input.vcf>   input filename of 2nd run (vcf format)\n");
  fprintf(stderr, "   <2nd input.vsf>   input filename of 2nd run (vsf format)\n");
  fprintf(stderr, "   <output.vcf>      output filename (vcf format)\n");
  fprintf(stderr, "   <type names>      names of bead types to save\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --join         join molecules (remove pbc)\n");
  fprintf(stderr, "      -st1 <int>     starting timestep from 1st run\n");
  fprintf(stderr, "      -st2 <int>     starting timestep from 2nd run\n");
  fprintf(stderr, "      -sk1 <int>     skip every <int> steps from 1st run\n");
  fprintf(stderr, "      -sk2 <int>     skip every <int> steps from 2nd run\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help nd exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
JoinRuns joins two simulation runs with different .vsf files. The first .vsf is \
assumed to be dl_meso.vsf (if not, use '-i' option) and the FIELD file has to \
be the same for both simulation runs. Bead types in both .vcf files must be the \
same, but only selected bead types are saved to output.vcf file.\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <1st input.vcf> <2nd input.vcf> <2nd input.vsf> ", argv[0]);
      fprintf(stdout, "<output.vcf> <type names> <options>\n\n");

      fprintf(stdout, "   <1st input.vcf>   input filename of 1st run (vcf format)\n");
      fprintf(stdout, "   <2nd input.vcf>   input filename of 2nd run (vcf format)\n");
      fprintf(stdout, "   <2nd input.vsf>   input filename of 2nd run (vsf format)\n");
      fprintf(stdout, "   <output.vcf>      output filename (vcf format)\n");
      fprintf(stdout, "   <type names>      names of bead types to save\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      --join         join molecules (remove pbc)\n");
      fprintf(stdout, "      -st1 <int>     starting timestep from 1st run\n");
      fprintf(stdout, "      -st2 <int>     starting timestep from 2nd run\n");
      fprintf(stdout, "      -sk1 <int>     skip every <int> steps from 1st run\n");
      fprintf(stdout, "      -sk2 <int>     skip every <int> steps from 2nd run\n");
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
        strcmp(argv[i], "--join") != 0 &&
        strcmp(argv[i], "-st1") != 0 &&
        strcmp(argv[i], "-st2") != 0 &&
        strcmp(argv[i], "-sk1") != 0 &&
        strcmp(argv[i], "-sk2") != 0) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // save coordinates of joined aggregates //{{{
  char joined_vcf[32];
  bool error = JoinCoorOption(argc, argv, joined_vcf);
  if (error) {
    exit(1);
  } //}}}

  // use .vsf file other than dl_meso.vsf? //{{{
  char *vsf_file_1 = calloc(32,sizeof(char *));
  if (VsfFileOption(argc, argv, &vsf_file_1)) {
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

  // should output coordinates be joined? //{{{
  bool join = BoolOption(argc, argv, "--join"); //}}}

  // starting timesteps for both simulations //{{{
  int start_1 = 1, start_2 = 1;
  if (IntegerOption(argc, argv, "-st1", &start_1)) {
    exit(1);
  }
  if (IntegerOption(argc, argv, "-st2", &start_2)) {
    exit(1);
  } //}}}

  // skipped timesteps per used step //{{{
  int skip_1 = 1, skip_2 = 1;
  if (IntegerOption(argc, argv, "-sk1", &skip_1)) {
    exit(1);
  }
  if (IntegerOption(argc, argv, "-sk2", &skip_2)) {
    exit(1);
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  count = 0; // count mandatory arguments

  // <1st input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf_1[32];
  strcpy(input_vcf_1, argv[++count]);

  // test if <1st input.vcf> filename ends with '.vcf' (required by VMD)
  char *dot = strrchr(input_vcf_1, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<input.vcf> '%s' does not have .vcf ending!\n", input_vcf_1);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <2nd input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf_2[32];
  strcpy(input_vcf_2, argv[++count]);

  // test if <1st input.vcf> filename ends with '.vcf' (required by VMD)
  dot = strrchr(input_vcf_2, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<input.vcf> '%s' does not have .vcf ending!\n", input_vcf_2);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <2nd input.vsf> - second structure file (must end with .vsf) //{{{
  char *vsf_file_2 = calloc(32,sizeof(char *));
  strcpy(vsf_file_2, argv[++count]);

  // test if <1st input.vcf> filename ends with '.vcf' (required by VMD)
  dot = strrchr(vsf_file_2, '.');
  if (!dot || strcmp(dot, ".vsf")) {
    fprintf(stderr, "<2nd input.vsf> '%s' does not have .vsf ending!\n", vsf_file_2);
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
  // data from 1st run
  BeadType *BeadType1; // structure with info about all bead types
  MoleculeType *MoleculeType1; // structure with info about all molecule types
  Bead *Bead1; // structure with info about every bead
  Molecule *Molecule1; // structure with info about every molecule
  // data from 2nd run
  BeadType *BeadType2; // structure with info about all bead types
  MoleculeType *MoleculeType2; // structure with info about all molecule types
  Bead *Bead2; // structure with info about every bead
  Molecule *Molecule2; // structure with info about every molecule
  // Counts is the same for both runs
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(vsf_file_1, input_vcf_1, bonds_file, &Counts, &BeadType1, &Bead1, &MoleculeType1, &Molecule1);
  ReadStructure(vsf_file_2, input_vcf_2, bonds_file, &Counts, &BeadType2, &Bead2, &MoleculeType2, &Molecule2);

  // vsf files are not needed anymore
  free(vsf_file_1);
  free(vsf_file_2);

  // <type names> - names of bead types to save //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType1);

    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s or %s coordinate file!\n", argv[count], input_vcf_1, input_vcf_2);
      exit(1);
    }

    BeadType1[type].Write = true;
  } //}}}

  // print selected bead type names to output .vcf file //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_vcf);
    exit(1);
  }

  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType1[i].Write) {
      fprintf(out, "# %s\n", BeadType1[i].Name);
    }
  }

  fclose(out); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf_1, bonds_file, Counts, BeadType1, Bead1, MoleculeType1, Molecule1);

    fprintf(stdout, "\n   Starting from %d. (%d.) timestep\n", start_1, start_2);
    fprintf(stdout, "   Every %d. (%d.) timestep used\n", skip_1+1, skip_2+1);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // open input coordinate files //{{{
  FILE *vcf_1, *vcf_2;
  if ((vcf_1 = fopen(input_vcf_1, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf_1);
    exit(1);
  }
  if ((vcf_2 = fopen(input_vcf_2, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf_2);
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // 1st vcf file - skip till 'pbc' keyword //{{{
  do {
    if (fscanf(vcf_1, "%s", str) != 1) {
      fprintf(stderr, "Cannot read string from '%s' file!\n", input_vcf_1);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf_1, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_vcf_1);
    exit(1);
  }
  // skip remainder of pbc line
  while (getc(vcf_1) != '\n')
    ;
  // skip blank line
  while (getc(vcf_1) != '\n')
    ; //}}}

  // 2nd vcf file - skip till 'pbc' keword //{{{
  do {
    if (fscanf(vcf_2, "%s", str) != 1) {
      fprintf(stderr, "Cannot read string from '%s' file!\n", input_vcf_2);
    }
  } while (strcmp(str, "pbc") != 0);
  // skip remainder of pbc line
  while (getc(vcf_2) != '\n')
    ;
  // skip blank line
  while (getc(vcf_2) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}
  //}}}

  // print pbc to output .vcf file //{{{
  if ((out = fopen(output_vcf, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_vcf);
    exit(1);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  fclose(out); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = calloc(128,sizeof(int)); //}}}

  // start first run with start-th step //{{{
  int test;
  // first run
  count = 0;
  for (int i = 1; i < start_1; i++) {
    count++;

    if (!silent) {
      if (script) {
        fprintf(stdout, "Discarded from 1st coordinate file: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rDiscarded from 1st coordinate file: %6d", count);
      }
    }

    SkipCoor(vcf_1, Counts, &stuff);
  }

  putchar('\n'); //}}}

  // main loop - 1st run //{{{
  while ((test = getc(vcf_1)) != EOF) {
    ungetc(test, vcf_1);

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf_1, Counts, &Bead1, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_1, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf_1, Counts, &Bead1, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_1, count, test);
        exit(1);
      }
    } //}}}

    count++;
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep from 1st run: %6d", count);
      }
    }

    // join molecules? //{{{
    if (join) {
      RemovePBCMolecules(Counts, BoxLength, BeadType1, &Bead1, MoleculeType1, Molecule1);
    } else { // if rounding leads to BoxLength, move it bead to other side of box
      for (int i = 0; i < (Counts.Bonded+Counts.Unbonded); i++) {
        char check[8];
        char box[8];
        // x direction
        sprintf(check, "%.3f", Bead1[i].Position.x);
        sprintf(box, "%.3f", BoxLength.x);
        if (strcmp(check, box) == 0) {
          Bead1[i].Position.x = 0;
        }
        // y direction
        sprintf(check, "%.3f", Bead1[i].Position.y);
        sprintf(box, "%.3f", BoxLength.y);
        if (strcmp(check, box) == 0) {
          Bead1[i].Position.y = 0;
        }
        // z direction
        sprintf(check, "%.3f", Bead1[i].Position.z);
        sprintf(box, "%.3f", BoxLength.z);
        if (strcmp(check, box) == 0) {
          Bead1[i].Position.z = 0;
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

    WriteCoorIndexed(out, Counts, BeadType1, Bead1, MoleculeType1, Molecule1, stuff);

    fclose(out);

    // skip every 'skip' steps //{{{
    for (int i = 0; i < skip_1; i++) {
      // test whether at vcf's eof //{{{
      if ((test = getc(vcf_1)) == EOF) {
        break;
      }
      ungetc(test, vcf_1); //}}}

      fflush(stdout);
      fprintf(stdout, "\rStep from 1st run: %6d", ++count);

      // read indexed timestep from input .vcf file //{{{
      if (indexed) {
        if ((test = ReadCoorIndexed(vcf_1, Counts, &Bead1, &stuff)) != 0) {
          // print newline to stdout if Step... doesn't end with one
          if (!script && !silent) {
            putchar('\n');
          }
          fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_1, count, test);
          exit(1);
        } //}}}
      // or read ordered timestep from input .vcf file //{{{
      } else {
        if ((test = ReadCoorOrdered(vcf_1, Counts, &Bead1, &stuff)) != 0) {
          // print newline to stdout if Step... doesn't end with one
          if (!script && !silent) {
            putchar('\n');
          }
          fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_1, count, test);
          exit(1);
        }
      } //}}}
    } //}}}

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      fprintf(stdout, "\n%s", stuff);
  }

  if (!silent) {
    fflush(stdout);
    fprintf(stdout, "\rLast Step from first run: %6d\n", count);
  }

  fclose(vcf_1); //}}}

  // start second run with start-th step //{{{
  count = 0;
  if (!silent)
    fprintf(stdout, "\rDiscarded from 2nd coordinate file: %6d", count);
  for (int i = 1; i < start_2; i++) {
    count++;

    if (!silent) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarded from 2nd coordinate file: %6d", count);
    }

    SkipCoor(vcf_2, Counts, &stuff);
  }

  putchar('\n'); //}}}

  // main loop - 2nd run //{{{
  while ((test = getc(vcf_2)) != EOF) {
    ungetc(test, vcf_2);

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf_2, Counts, &Bead2, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_2, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf_2, Counts, &Bead2, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_2, count, test);
        exit(1);
      }
    } //}}}

    count++;
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step from 2nd run: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep from 2nd run: %6d", count);
      }
    }

    // join molecules? //{{{
    if (join) {
      RemovePBCMolecules(Counts, BoxLength, BeadType2, &Bead2, MoleculeType2, Molecule2);
    } else { // if rounding leads to BoxLength, move it bead to other side of box
      for (int i = 0; i < (Counts.Bonded+Counts.Unbonded); i++) {
        char check[8];
        char box[8];
        // x direction
        sprintf(check, "%.3f", Bead2[i].Position.x);
        sprintf(box, "%.3f", BoxLength.x);
        if (strcmp(check, box) == 0) {
          Bead2[i].Position.x = 0;
        }
        // y direction
        sprintf(check, "%.3f", Bead2[i].Position.y);
        sprintf(box, "%.3f", BoxLength.y);
        if (strcmp(check, box) == 0) {
          Bead2[i].Position.y = 0;
        }
        // z direction
        sprintf(check, "%.3f", Bead2[i].Position.z);
        sprintf(box, "%.3f", BoxLength.z);
        if (strcmp(check, box) == 0) {
          Bead2[i].Position.z = 0;
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

    // copy Bead2 coordinates to Bead1 //{{{
    bool used[Counts.Unbonded];
    for (int i = 0; i < Counts.Unbonded; i++) {
      used[i] = false;
    }

    for (int i = 0; i < Counts.Unbonded; i++) {
      for (int j = 0; j < Counts.Unbonded; j++) {
        if (strcmp(BeadType2[Bead2[j].Type].Name, BeadType1[Bead1[i].Type].Name) == 0
            && !used[j]) {
          Bead1[i].Position.x = Bead2[j].Position.x;
          Bead1[i].Position.y = Bead2[j].Position.y;
          Bead1[i].Position.z = Bead2[j].Position.z;

          used[j] = true;

          break;
        }
      }
    }

    for (int i = 0; i < Counts.Unbonded; i++) {
      if (!used[i]) {
        fprintf(stdout, "ERROR - used[%d] = false\n", i);
      }
    }

    for (int i = Counts.Unbonded; i < (Counts.Unbonded+Counts.Bonded); i++) {
      Bead1[i].Position.x = Bead2[i].Position.x;
      Bead1[i].Position.y = Bead2[i].Position.y;
      Bead1[i].Position.z = Bead2[i].Position.z;
    } //}}}

    // TODO: there probably shouldn't be Molecule(Type)1, but Molecule(Type)2
    // or something different -- this program isn't used anyway
    // TODO: implement -x or maybe the write should be hardcoded in?
    WriteCoorIndexed(out, Counts, BeadType1, Bead1, MoleculeType1, Molecule1, stuff);

    fclose(out);

    // skip every 'skip' steps //{{{
    for (int i = 0; i < skip_2; i++) {
      // test whether at vcf's eof //{{{
      if ((test = getc(vcf_2)) == EOF) {
        break;
      }
      ungetc(test, vcf_2); //}}}

      fflush(stdout);
      fprintf(stdout, "\rStep from 2st run: %6d", ++count);

      // read indexed timestep from input .vcf file //{{{
      if (indexed) {
        if ((test = ReadCoorIndexed(vcf_2, Counts, &Bead2, &stuff)) != 0) {
          // print newline to stdout if Step... doesn't end with one
          if (!script && !silent) {
            putchar('\n');
          }
          fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_2, count, test);
          exit(1);
        } //}}}
      // or read ordered timestep from input .vcf file //{{{
      } else {
        if ((test = ReadCoorOrdered(vcf_2, Counts, &Bead2, &stuff)) != 0) {
          // print newline to stdout if Step... doesn't end with one
          if (!script && !silent) {
            putchar('\n');
          }
          fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf_2, count, test);
          exit(1);
        }
      } //}}}
    } //}}}

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      fprintf(stdout, "\n%s", stuff);
  }

  if (!silent) {
    fflush(stdout);
    fprintf(stdout, "\rLast Step from 2nd run: %6d\n", count);
  }

  fclose(vcf_2); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType1);
  free(BeadType2);
  FreeMoleculeType(Counts, &MoleculeType1);
  FreeMoleculeType(Counts, &MoleculeType2);
  FreeMolecule(Counts, &Molecule1);
  FreeMolecule(Counts, &Molecule2);
  FreeBead(Counts, &Bead1);
  FreeBead(Counts, &Bead2);
  free(stuff);
  //}}}

  return 0;
}
