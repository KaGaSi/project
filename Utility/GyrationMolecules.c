#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
GyrationMolecules calculates the gyration tensor for molecules and \
determines their shape descriptors like the radius of gyration, \
acylindricity, asphericity, or relative shape anisotropy. It writes \
per-timestep averages to the output file and appends overall averages to \
that file.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <output> <mol name(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <output>          output file with shape descriptors (automatic ending '-<name>.txt')\n");
  fprintf(ptr, "   <mol name(s)>     molecule types to calculate shape descriptors for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined coordinates\n");
  fprintf(ptr, "      -bt            specify bead types to be used for calculation (default is all)\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
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

  if (count < req_args) {
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
        strcmp(argv[i], "--joined") &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0) {

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

  // <output> - filename with shape descriptors //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

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

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}
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

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);
    if (mol_type == -1) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - non-existent molecule");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", argv[count]);
      RedText(STDERR_FILENO);
      fprintf(stderr, "\n");
      ResetColour(STDERR_FILENO);
      ErrorMoleculeType(Counts, MoleculeType);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }
  } //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", true, Counts, &BeadType)) {
    exit(0);
  } //}}}

  // write initial stuff to output file //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      char str2[1050];
      sprintf(str2, "%s-%s.txt", output, MoleculeType[i].Name);

      FILE *out;
      if ((out = fopen(str2, "w")) == NULL) {
        ErrorFileOpen(str2, 'w');
        exit(1);
      }

      // print command to output file
      putc('#', out);
      PrintCommand(out, argc, argv);

      fprintf(out, "# %s\n", MoleculeType[i].Name);
      fprintf(out, "# (1) dt, (2) <Rg>, (3) <Rg^2>, ");
      fprintf(out, "(4) <Anis>, (5) <Acyl>, (6) <Aspher>");
      fprintf(out, "(7) <eigen.x>, (8) <eigen.y>, (9) <eigen.z>");
      putc('\n', out);

      fclose(out);
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // allocate memory for sums of shape descriptors //{{{
  double *Rg_sum = calloc(Counts.TypesOfMolecules, sizeof(double));
  double *sqrRg_sum = calloc(Counts.TypesOfMolecules, sizeof(double));
  double *Anis_sum = calloc(Counts.TypesOfMolecules, sizeof(double));
  double *Acyl_sum = calloc(Counts.TypesOfMolecules, sizeof(double));
  double *Aspher_sum = calloc(Counts.TypesOfMolecules, sizeof(double));
  VECTOR *eigen_sum = calloc(Counts.TypesOfMolecules, sizeof(VECTOR)); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vtf, vcf, struct_lines); //}}}

  // main loop //{{{
  count = 0; // count timesteps
  while (true) {
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}

    ReadVcfCoordinates(indexed, input_coor, vcf, Counts, Index, &Bead, &stuff);

    // join molecules if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate arrays for the timestep //{{{
    double *Rg_step = calloc(Counts.TypesOfMolecules, sizeof(double));
    double *sqrRg_step = calloc(Counts.TypesOfMolecules, sizeof(double));
    double *Anis_step = calloc(Counts.TypesOfMolecules,sizeof(double));
    double *Acyl_step = calloc(Counts.TypesOfMolecules,sizeof(double));
    double *Aspher_step = calloc(Counts.TypesOfMolecules,sizeof(double));
    VECTOR *eigen_step = calloc(Counts.TypesOfMolecules,sizeof(VECTOR)); //}}}

    // calculate shape descriptors //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type = Molecule[i].Type;

      if (MoleculeType[mol_type].Use) {
        // copy bead ids to a separate array //{{{
        int *list = malloc(MoleculeType[mol_type].nBeads*sizeof(int));
        int n = 0;
        for (int j = 0; j < MoleculeType[mol_type].nBeads; j++) {
          int bead_id = Molecule[i].Bead[j];
          if (BeadType[Bead[bead_id].Type].Use) {
            list[n] = Molecule[i].Bead[j];
            n++;
          }
        } //}}}

        VECTOR eigen = Gyration(n, list, Counts, BoxLength, BeadType, &Bead);

        free(list); // free array of bead ids for gyration calculation

        double Rgi = sqrt(eigen.x + eigen.y + eigen.z);

        if (eigen.x < 0 || eigen.y < 0 || eigen.z < 0) {
          YellowText(STDERR_FILENO);
          fprintf(stderr, "\nWarning: negative eigenvalues (%lf, %lf, %lf)\n", eigen.x, eigen.y, eigen.z);
          ResetColour(STDERR_FILENO);
        }
        // radius of gyration
        Rg_step[mol_type] += Rgi; // for number avg
        // squared radius of gyration
        sqrRg_step[mol_type] += SQR(Rgi); // for number avg
        // relative shape anisotropy
        Anis_step[mol_type] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) / SQR(eigen.x + eigen.y + eigen.z) - 0.5;
        // acylindricity
        Acyl_step[mol_type] += eigen.y - eigen.x;
        // asphericity
        Aspher_step[mol_type] += eigen.z - 0.5 * (eigen.x + eigen.y);
        // eigenvalues
        eigen_step[mol_type].x += eigen.x;
        eigen_step[mol_type].y += eigen.y;
        eigen_step[mol_type].z += eigen.z;
      }
    } //}}}

    // add values to sums //{{{
    if (count >= start && (end == -1 || count <= end)) {
      for (int i = 0; i < Counts.TypesOfMolecules; i++) {
        Rg_sum[i] += Rg_step[i];
        sqrRg_sum[i] += sqrRg_step[i];
        Anis_sum[i] += Anis_step[i];
        Acyl_sum[i] += Acyl_step[i];
        Aspher_sum[i] += Aspher_step[i];
        eigen_sum[i].x += eigen_step[i].x;
        eigen_sum[i].y += eigen_step[i].y;
        eigen_sum[i].z += eigen_step[i].z;
      }
    } //}}}

    // print shape descriptors to output file(s) //{{{
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        char str2[1050];
        sprintf(str2, "%s-%s.txt", output, MoleculeType[i].Name);

        FILE *out;
        if ((out = fopen(str2, "a")) == NULL) {
          ErrorFileOpen(str2, 'a');
          exit(1);
        }

        fprintf(out, "%5d", count);
        fprintf(out, " %8.5f", Rg_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", sqrRg_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Anis_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Acyl_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", Aspher_step[i]/MoleculeType[i].Number);
        fprintf(out, " %8.5f", eigen_step[i].x/MoleculeType[i].Number);
        fprintf(out, " %8.5f", eigen_step[i].y/MoleculeType[i].Number);
        fprintf(out, " %8.5f", eigen_step[i].z/MoleculeType[i].Number);
        putc('\n', out);

        fclose(out);
      }
    } //}}}

    // free memory //{{{
    free(Rg_step);
    free(sqrRg_step);
    free(Anis_step);
    free(Acyl_step);
    free(Aspher_step);
    free(eigen_step); //}}}

    // -e option - exit main loop if last step is done
    if (end == count) {
      break;
    }
    // if there's no additional timestep, exit the while loop
    bool rubbish; // not used
    if (ReadTimestepPreamble(&rubbish, input_coor, vcf, &stuff, false) == -1) {
      break;
    }
  }
  fclose(vcf);

  // print last step? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count);
  } //}}}
  //}}}

  // write simple averages to <output> //{{{
  if (end != -1) {
    count = end - start + 1;
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      // open file //{{{
      char str2[1050];
      sprintf(str2, "%s-%s.txt", output, MoleculeType[i].Name);
      FILE *out;
      if ((out = fopen(str2, "a")) == NULL) {
        ErrorFileOpen(str2, 'a');
        exit(1);
      } //}}}

      // write legend line to output file //{{{
      fprintf(out, "# %s\n", MoleculeType[i].Name);
      fprintf(out, "# simple averages (from steps %d to %d): ", start, start+count);
      int counter = 1;
      fprintf(out, "(%d) <Rg>, ", counter++);
      fprintf(out, "(%d) <Anis>, ", counter++);
      fprintf(out, "(%d) <Acyl>, ", counter++);
      fprintf(out, "(%d) <Aspher>", counter++);
//    fprintf(out, "(%d) <eigen.x>, ", counter++);
//    fprintf(out, "(%d) <eigen.y>, ", counter++);
//    fprintf(out, "(%d) <eigen.z>, ", counter++);
      putc('\n', out); //}}}

      // write averages to output file //{{{
      putc('#', out);
      fprintf(out, " %lf", Rg_sum[i]/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", Anis_sum[i]/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", Acyl_sum[i]/(count*MoleculeType[i].Number));
      fprintf(out, " %lf", Aspher_sum[i]/(count*MoleculeType[i].Number));
      putc('\n', out); //}}}

      fclose(out);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(Rg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(stuff); //}}}

  return 0;
}
