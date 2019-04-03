#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
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
SelectedVcf creates new <output> file from <input> containing only \
selected bead types. Also <start> timesteps can be omitted and every <skip> \
timestep can be left out.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> ", cmd);
  fprintf(ptr, "<output> <type names> <options>\n\n");

  fprintf(ptr, "   <input>           input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <output>          output filename (vcf format)\n");
  fprintf(ptr, "   <type names>      names of bead types to save (optional if '-r' used)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -r             reverse <type name(s)>, i.e., exclude the specified bead types\n");
  fprintf(ptr, "      --join         join molecules (remove pbc)\n");
  fprintf(ptr, "      -w             wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "      -st <start>    number of timestep to start from\n");
  fprintf(ptr, "      -e <end>       number of timestep to end with\n");
  fprintf(ptr, "      -sk <skip>     leave out every 'skip' steps\n");
  fprintf(ptr, "      -n <int(s)>    save only specified timesteps\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  fprintf(ptr, "      -xyz <name>    output xyz file\n");
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
  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  // reverse bead type selection? ...do now to check correct number of arguments //{{{
  bool reverse = BoolOption(argc, argv, "-r");
  // }}}

  // possible to exclude <type name(s)> if '-r' is used
  if (count < (req_args-1) || (count == (req_args-1) && !reverse)) {
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
        strcmp(argv[i], "--join") != 0 &&
        strcmp(argv[i], "-w") != 0 &&
        strcmp(argv[i], "-r") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-xyz") != 0 &&
        strcmp(argv[i], "-x") != 0) {

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
  char **extension;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_vsf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(1024,sizeof(char *));
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

  // should output coordinates be joined?
  bool join = BoolOption(argc, argv, "--join");

  // should output coordinates be wrapped?
  bool wrap = BoolOption(argc, argv, "-w");

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // ending timestep //{{{
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  } //}}}

  // error if ending step is lower than starging step //{{{
  if (start > end) {
    fprintf(stderr, "\nErorr: Starting step (%d) is higher than ending step (%d)\n", start, end);
    exit(1);
  } //}}}

  // number of steps to skip per one used //{{{
  int skip = 0;
  if (IntegerOption(argc, argv, "-sk", &skip)) {
    exit(1);
  } //}}}

  // save into xyz file? //{{{
  char *output_xyz = calloc(1024,sizeof(char *));
  if (FileOption(argc, argv, "-xyz", &output_xyz)) {
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
  char input_vcf[1024];
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
    Help(argv[0], true);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // <output.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[1024];
  strcpy(output_vcf, argv[++count]);

  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  extension = malloc(ext*sizeof(char *));
  extension[0] = malloc(5*sizeof(char));
  strcpy(extension[0], ".vcf");
  if (!ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <type names> - names of bead types to save //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);

    if (type == -1) {
      fprintf(stderr, "\nError: bead type '%s' is not in %s file\n\nPresent bead types:\n", argv[count], input_vcf);
      for (int i = 0; i < Counts.TypesOfBeads; i++) {
        fprintf(stderr, "   %s\n", BeadType[i].Name);
      }
      exit(1);
    }

    BeadType[type].Write = true;
  } //}}}

  // if '-r' is used, switch Write bools for all bead types //{{{
  if (reverse) {
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      if (BeadType[i].Write) {
        BeadType[i].Write = false;
      } else {
        BeadType[i].Write = true;
      }
    }
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }

  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  }
  //}}}

  // '-n' option - specify bead ids //{{{
  int save_step[100] = {0}, number_of_steps = 0, test = 0;
  if (MultiIntegerOption(argc, argv, "-n", &test, save_step)) {
    exit(1);
  }
  if (test != 0) { // -n is present
    number_of_steps = test;

    // sort from lowest to highest
    for (int i = 0; i < (number_of_steps-1); i++) {
      for (int j = (i+1); j < number_of_steps; j++) {
        if (save_step[i] > save_step[j]) {
          int swap = save_step[i];
          save_step[i] = save_step[j];
          save_step[j] = swap;
        }
      }
    }
  }
  //}}}

  // print selected bead type names to output .vcf file //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
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

    fprintf(stdout, "\n   Starting from %d. timestep\n", start);
    fprintf(stdout, "   Every %d. timestep used\n", skip+1);
    if (end != -1) {
      fprintf(stdout, "   Ending with %d. timestep\n", end);
    }

    if (number_of_steps > 0) {
      fprintf(stdout, "   Save %d timesteps:", number_of_steps);
      for (int i = 0; i < number_of_steps; i++) {
        fprintf(stdout, " %d", save_step[i]);
        if (i != (number_of_steps-1)) {
          putchar(',');
        } else {
          putchar('\n');
        }
      }
    }
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    ErrorFileOpen(input_vcf, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[1024];
  // skip till 'pbc' keyword //{{{
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_vcf);
      exit(1);
    }
  } while (strcmp(str, "pbc") != 0); //}}}

  // read pbc //{{{
  Vector BoxLength;
  char line[1024];
  fgets(line, sizeof(line), vcf);
  // split the line into array
  char *split[30];
  split[0] = strtok(line, " \t");
  int i = 0;
  while (split[i] != NULL && i < 29) {
    split[++i] = strtok(NULL, " \t");
  }
  BoxLength.x = atof(split[0]);
  BoxLength.y = atof(split[1]);
  BoxLength.z = atof(split[2]); //}}}
  //}}}

  // print pbc if verbose output //{{{
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // print pbc to output .vcf file //{{{
  if ((out = fopen(output_vcf, "a")) == NULL) {
    ErrorFileOpen(output_vcf, 'a');
    exit(1);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  fclose(out); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = calloc(1024,sizeof(int)); //}}}

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
      fprintf(stderr, "\nError: premature end of %s file\n\n", input_vcf);
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

  // is the vcf file continuing? //{{{
  if ((test = getc(vcf)) == EOF) {
    fprintf(stderr, "\nError: %s - number of discard steps is lower (or equal) to the total number of steps\n\n", input_vcf);
    exit(1);
  } else {
    ungetc(test,vcf);
  } //}}}
  //}}}

  // main loop //{{{
  int count_n_opt = 0;
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

    if (number_of_steps != 0) {
      if (count_n_opt < number_of_steps) {
        if (save_step[count_n_opt] == count) {
          // read coordinates //{{{
          if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
            // print newline to stdout if Step... doesn't end with one
            ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
            exit(1);
          } //}}}

          // wrap coordinates? //{{{
          if (wrap) {
            RestorePBC(Counts, BoxLength, &Bead);
          } // }}}

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
            ErrorFileOpen(output_vcf, 'a');
            exit(1);
          } //}}}

          WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

          fclose(out);

          // save to xyz file? //{{{
          if (output_xyz[0] != '\0') {

            // open output .xyz file for appending //{{{
            if ((out = fopen(output_xyz, "a")) == NULL) {
              ErrorFileOpen(output_xyz, 'a');
              exit(1);
            } //}}}

            WriteCoorXYZ(out, Counts, BeadType, Bead);

            fclose(out);
          } //}}}

          count_n_opt++;
        } else {
          if (SkipCoor(vcf, Counts, &stuff)) {
            fprintf(stderr, "\nError: premature end of %s file\n\n", input_vcf);
            exit(1);
          }
        }
      } else { // all required steps saved
        break;
      }
    } else {
      // read coordinates //{{{
      if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
        exit(1);
      } //}}}

      // wrap coordinates? //{{{
      if (wrap) {
        RestorePBC(Counts, BoxLength, &Bead);
      } // }}}

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
        ErrorFileOpen(output_vcf, 'a');
        exit(1);
      } //}}}

      WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

      fclose(out);

      // save to xyz file? //{{{
      if (output_xyz[0] != '\0') {

        // open output .xyz file for appending //{{{
        if ((out = fopen(output_xyz, "a")) == NULL) {
          ErrorFileOpen(output_xyz, 'a');
          exit(1);
        } //}}}

        WriteCoorXYZ(out, Counts, BeadType, Bead);

        fclose(out);
      } //}}}

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
            fprintf(stdout, "Step: %6d\n", count);
          } else {
            fprintf(stdout, "\rStep: %6d", ++count);
          }
        }

  //    // read coordinates //{{{
  //    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
  //      // print newline to stdout if Step... doesn't end with one
  //      ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
  //      exit(1);
  //    } //}}}
        if (SkipCoor(vcf, Counts, &stuff)) {
          fprintf(stderr, "\nError: premature end of %s file\n\n", input_vcf);
          exit(1);
        }
      } //}}}
    }

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      fprintf(stdout, "\n%s", stuff);

    if (end == count)
      break;
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

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(output_xyz);
  //}}}

  return 0;
}
