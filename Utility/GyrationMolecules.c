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
        strcmp(argv[i], "--script") != 0 &&
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
  char *input_vsf = calloc(LINE,sizeof(char));
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
  } //}}}

  // <output> - filename with shape descriptors //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

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
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n\n", start, end);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {

    bool test = false;
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (strcmp(argv[count], MoleculeType[i].Name) == 0) {
        MoleculeType[i].Use = true;

        test = true;

        break;
      }
    }

    // wrong molecule name //{{{
    if (!test) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - non-existent molecule name %s\n\n", input_vsf, argv[count]);
      fprintf(stderr, "\033[0m");
      ErrorMoleculeType(Counts, MoleculeType);
      exit(1);
    } //}}}
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

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

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

  // main loop //{{{
  count = 0; // count timesteps
  int test;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      ErrorCoorRead(input_coor, test, count, stuff);
      exit(1);
    } //}}}

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
          fprintf(stderr, "\033[1;33m");
          fprintf(stderr, "\nWarning: negative eigenvalues (%lf, %lf, %lf)\n", eigen.x, eigen.y, eigen.z);
          fprintf(stderr, "\033[0m");
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
  }
  fclose(vcf);

  // print last step? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                      ");
      fprintf(stdout, "\rLast Step: %d\n", count);
    }
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
  free(input_vsf);
  free(BeadType);
  free(Index);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(Rg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(stuff); //}}}

  return 0;
}
