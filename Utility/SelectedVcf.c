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
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[LINE];
  strcpy(output_vcf, argv[++count]);
  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // should output coordinates be joined?
  bool join = BoolOption(argc, argv, "--join");

  // should output coordinates be wrapped?
  bool wrap = BoolOption(argc, argv, "-w");

  // starting & ending timesteps //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  }
  int end = -1;
  if (IntegerOption(argc, argv, "-e", &end)) {
    exit(1);
  }
  ErrorStartEnd(start, end); //}}}

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

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  VECTOR BoxLength; // couboid box dimensions
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &BoxLength, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule);
  free(input_vsf); //}}}

  // <type names> - names of bead types to save //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);
    if (type == -1) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - non-existent bead name ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s\n", argv[count]);
      ResetColour(STDERR_FILENO);
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
  } //}}}

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
  }
  SkipVtfStructure(vtf, vcf, struct_lines); //}}}

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

  // create array for the timestep preamble
  char *stuff = calloc(LINE, sizeof(char));

  count = SkipCoorSteps(vcf, input_coor, Counts, start, silent);

  // main loop //{{{
  int count_n_opt = 0; // count saved steps if -n option is used
  count = 0; // count timesteps in the main loop
  int count_vcf = start - 1; // count timesteps from the beginning
  while (true) {
    count++;
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // read coordinates & wrap/join molecules //{{{
    ReadVcfCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);
    // wrap and/or join molecules?
    if (wrap) {
      RestorePBC2(Counts.Beads, BoxLength, &Bead);
    }
    if (join) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    if (!last) { // save coordinates only if --last option is not used
      if (count_n_opt < number_of_steps) { // if -n option is used
        if (save_step[count_n_opt] == count_vcf) {
          // write to output .vcf file //{{{
          if ((out = fopen(output_vcf, "a")) == NULL) {
            ErrorFileOpen(output_vcf, 'a');
            exit(1);
          }
          WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
          fclose(out); //}}}
          // write to xyz file? //{{{
          if (output_xyz[0] != '\0') {
            if ((out = fopen(output_xyz, "a")) == NULL) {
              ErrorFileOpen(output_xyz, 'a');
              exit(1);
            }
            WriteCoorXYZ(out, Counts, BeadType, Bead);
            fclose(out);
          } //}}}
          count_n_opt++;
        }
        if (count_n_opt == number_of_steps) {
          break;
        }
      } else { // if -n option is not used
        // write to output .vcf file //{{{
        if ((out = fopen(output_vcf, "a")) == NULL) {
          ErrorFileOpen(output_vcf, 'a');
          exit(1);
        }
        WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
        fclose(out); //}}}
        // write to xyz file? //{{{
        if (output_xyz[0] != '\0') {
          if ((out = fopen(output_xyz, "a")) == NULL) {
            ErrorFileOpen(output_xyz, 'a');
            exit(1);
          }
          WriteCoorXYZ(out, Counts, BeadType, Bead);
          fclose(out);
        } //}}}
        // skip every 'skip' steps //{{{
        for (int i = 0; i < skip; i++) {
          bool rubbish; // not used
          if (ReadTimestepPreamble(&rubbish, input_coor, vcf, &stuff, false) == -1) {
            break;
          }
          count_vcf++;
          count++;
          // -e option
          if (end == count_vcf) {
            break;
          }
          if (!silent && isatty(STDOUT_FILENO)) {
            fflush(stdout);
            fprintf(stdout, "\rStep: %d", count_vcf);
          }
          SkipVcfCoor(vcf, input_coor, Counts, &stuff);
        } //}}}
      }
    }

    // -e option - exit main loop if last step is done
    if (end == count_vcf) {
      break;
    }
    // if there's no additional timestep, exit the while loop
    bool rubbish; // not used
    if (ReadTimestepPreamble(&rubbish, input_coor, vcf, &stuff, false) == -1) {
      break;
    }
  }
  fclose(vcf);

  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // save last step if --last is used //{{{
  if (last) {
    // write to output .vcf file //{{{
    if ((out = fopen(output_vcf, "a")) == NULL) {
      ErrorFileOpen(output_vcf, 'a');
      exit(1);
    }
    WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
    fclose(out); //}}}
    // write to xyz file? //{{{
    if (output_xyz[0] != '\0') {
      // open output .xyz file for appending
      if ((out = fopen(output_xyz, "a")) == NULL) {
        ErrorFileOpen(output_xyz, 'a');
        exit(1);
      }
      WriteCoorXYZ(out, Counts, BeadType, Bead);
      fclose(out);
    } //}}}
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  free(output_xyz);
  //}}}

  return 0;
}
