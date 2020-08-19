#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
SelectedVcf creates new <output.vcf> file (and possibly xyz file) from <input> \
containing only selected bead types. Periodic boundary conditions can be either \
stripped away or applied (which happens first if both '--join' and '-w' options \
are used). Also, specified molecules can be excluded. However, AnalysisTools \
utilities can only read coordinate files containing all beads of any given \
type, so the usefulness is very limited (for, e.g., visualization using \
vmd).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> ", cmd);
  fprintf(ptr, "<output> <type names> <options>\n\n");

  fprintf(ptr, "   <input>           input filename (vcf or vtf format)\n");
  fprintf(ptr, "   <output.vcf>      output filename (vcf format)\n");
  fprintf(ptr, "   <type names>      names of bead types to save (optional if '-r' used)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -r             reverse <type name(s)>, i.e., exclude the specified bead types\n");
  fprintf(ptr, "      --join         join molecules (remove pbc)\n");
  fprintf(ptr, "      -w             wrap coordinates (i.e., apply pbc)\n");
  fprintf(ptr, "      -st <start>    starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  fprintf(ptr, "      -sk <skip>     leave out every 'skip' steps\n");
  fprintf(ptr, "      -n <int(s)>    save only specified timesteps\n");
  fprintf(ptr, "      -x <name(s)>   exclude specified molecule(s)\n");
  fprintf(ptr, "      -xyz <name>    output xyz file\n");
  fprintf(ptr, "      --last         use only the last step\n");
  CommonHelp(error);
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
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

  // reverse bead type selection? ...do now to check correct number of arguments
  bool reverse = BoolOption(argc, argv, "-r");

  // possible to omit <type name(s)> if '-r' is used
  if (count < (req_args-1) || (count == (req_args-1) && !reverse)) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-r") != 0 &&
        strcmp(argv[i], "--join") != 0 &&
        strcmp(argv[i], "-w") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-xyz") != 0 &&
        strcmp(argv[i], "--last") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE];
  char *input_vsf = calloc(LINE,sizeof(char));
  strcpy(input_coor, argv[++count]);

  // test if <input> filename ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  // if vtf, copy to input_vsf
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
  } //}}}

  // <output.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[LINE];
  strcpy(output_vcf, argv[++count]);

  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

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
  if (end != -1 && start > end) {
    fprintf(stderr, "\nErorr: starting step (%d) is higher than ending step (%d)\n", start, end);
    exit(1);
  } //}}}

  // number of steps to skip per one used //{{{
  int skip = 0;
  if (IntegerOption(argc, argv, "-sk", &skip)) {
    exit(1);
  } //}}}

  // save into xyz file? //{{{
  char *output_xyz = calloc(LINE,sizeof(char));
  if (FileOption(argc, argv, "-xyz", &output_xyz)) {
    exit(1);
  } //}}}

  // use only the last step?
  bool last = BoolOption(argc, argv, "--last");
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);
  free(input_vsf);

  // <type names> - names of bead types to save //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);
    if (type == -1) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - non-existent bead name \033[1;33m%s\033[1;31m\n", input_coor, argv[count]);
      fprintf(stderr, "\033[0m");
      ErrorBeadType(Counts, BeadType);
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

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);

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
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // print pbc to output .vcf file //{{{
  if ((out = fopen(output_vcf, "a")) == NULL) {
    ErrorFileOpen(output_vcf, 'a');
    exit(1);
  }

  fprintf(out, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

  fclose(out); //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // skip first start-1 steps //{{{
  count = 0;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
      exit(1);
    }
  }
  // print starting step? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Starting step: %d\n", start);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rStarting step: %d\n", start);
    }
  } //}}}
  // is the vcf file continuing?
  if (ErrorDiscard(start, count, input_coor, vcf)) {
    exit(1);
  } //}}}

  // main loop //{{{
  int count_n_opt = 0; // count saved steps if -n option is used
  count = 0; // count timesteps in the main loop
  int count_vcf = start - 1; // count timesteps from the beginning
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // if --last is used, just read coordinates and move on //{{{
    if (last) {
      // read coordinates
      ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);
      continue;
    } //}}}

    if (number_of_steps != 0) {
      if (count_n_opt < number_of_steps) {
        if (save_step[count_n_opt] == count_vcf) {
          ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

          // wrap coordinates? //{{{
          if (wrap) {
            RestorePBC(Counts, BoxLength, &Bead);
          } // }}}

          // join molecules? //{{{
          if (join) {
            RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
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
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
            fprintf(stderr, "\033[0m");
            exit(1);
          }
        }
      } else { // all required steps saved
        break;
      }
    } else {
      ReadCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

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
        // test whether at vcf's eof or reached end due to -e option //{{{
        if ((test = getc(vcf)) == EOF ||
            end == count_vcf) {
          break;
        }
        ungetc(test, vcf); //}}}

        count_vcf++;
        count++;
        if (!silent && !script) {
          fflush(stdout);
          fprintf(stdout, "\rStep: %d", count_vcf);
        }

        if (SkipCoor(vcf, Counts, &stuff)) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
          fprintf(stderr, "\033[0m");
          exit(1);
        }
      } //}}}
    }

    if (end == count_vcf)
      break;
  }

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rLast Step: %d\n", count_vcf);
    }
  }

  fclose(vcf); //}}}

  // save last step if --last is used //{{{
  if (last) {
    // open .vcf file for appending
    if ((out = fopen(output_vcf, "a")) == NULL) {
      ErrorFileOpen(output_vcf, 'a');
      exit(1);
    }

    // wrap coordinates? //{{{
    if (wrap) {
      RestorePBC(Counts, BoxLength, &Bead);
    } // }}}

    // join molecules? //{{{
    if (join) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
    fclose(out);

    // save to xyz file?
    if (output_xyz[0] != '\0') {
      // open output .xyz file for appending
      if ((out = fopen(output_xyz, "a")) == NULL) {
        ErrorFileOpen(output_xyz, 'a');
        exit(1);
      }
      WriteCoorXYZ(out, Counts, BeadType, Bead);
      fclose(out);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(output_xyz);
  //}}}

  return 0;
}
