#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
GyrationAggregates calculates gyration tensor for aggregates and determines \
shape descriptors like radius of gyration, acylindricity, asphericity, or \
relative shape anisotropy. By default, it calculates per-timestep averages, \
but per-size averages can also be determined. Overall averages are appended \
to the output file. The definition of aggregate size is quite flexible and \
the calculation can also be made only for aggregate sizes in a given \
range.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <input.agg> <output> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <input.agg>       input agg file\n");
  fprintf(ptr, "   <output>          output file with data during simulation run\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined coordinates\n");
  fprintf(ptr, "      -bt            specify bead types to be used for calculation (default is all)\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> molecule types in an aggregate\n");
  fprintf(ptr, "      -x <name(s)>   exclude aggregates containing only specified molecule(s)\n");
  fprintf(ptr, "      -ps <file>     save per-size averages to a file\n");
  fprintf(ptr, "      -n <int> <int> calculate for aggregate sizes in given range\n");
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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-bt") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-ps") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
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

  // <input.agg> - filename of input agg file //{{{
  char input_agg[LINE];
  strcpy(input_agg, argv[++count]);

  // test if <input.agg> filename ends with '.agg'
  ext = 1;
  strcpy(extension[0], ".agg");
  if (ErrorExtension(input_agg, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output> - filename with data during simulation run //{{{
  char output[LINE];
  strcpy(output, argv[++count]); //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

  // write per-agg averages to a file? //{{{
  char *per_size_file = calloc(LINE,sizeof(char));
  if (FileOption(argc, argv, "-ps", &per_size_file)) {
    exit(0);
  } //}}}

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
    fprintf(stderr, "\nError: Starting step (%d) is higher than ending step (%d)\n\n", start, end);
    exit(1);
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // '-n' option - range of aggregation numbers //{{{
  int range_As[2], test = 2;
  range_As[0] = 1;
  range_As[1] = Counts.Molecules;
  if (MultiIntegerOption(argc, argv, "-n", &test, range_As)) {
    exit(1);
  }
  if (test != 2) {
    fprintf(stderr, "\nError: option '-n' needs two numeric arguments\n\n");
    exit(1);
  }

  // make sure first number is larger
  if (range_As[0] > range_As[1]) {
    int tmp = range_As[0];
    range_As[0] = range_As[1];
    range_As[1] = tmp;
  } //}}}

  // '-m' option //{{{
  int *specific_moltype_for_size;
  specific_moltype_for_size = malloc(Counts.TypesOfMolecules*sizeof(int *));
  // all are to be used without '-m' option
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", &specific_moltype_for_size, Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }

  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  }

  // count total number of chains in excluded aggs
  long int exclude_count_chains = 0;
  // count total number of excluded aggs
  long int exclude_count_agg = 0; //}}}

  // -bt <name(s)> - specify what bead types to use //{{{
  if (BeadTypeOption(argc, argv, "-bt", true, Counts, &BeadType)) {
    exit(0);
  } //}}}

  bool types_for_gyr[Counts.TypesOfBeads];
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    types_for_gyr[i] = BeadType[i].Use;
    BeadType[i].Use = false;
  }

  // write initial stuff to output file //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    ErrorFileOpen(output, 'w');
    exit(1);
  }

  // print command to output file
  putc('#', out);
  for (int i = 0; i < argc; i++) {
    fprintf(out, " %s", argv[i]);
  }
  putc('\n', out);

  // print legend line to output file
  count = 1;
  fprintf(out, "# column: ");
  fprintf(out, "(%d) timestep", count++);
  fprintf(out, ", (%d) <Rg>_n", count++);
  fprintf(out, ", (%d) <Rg>_w", count++);
  fprintf(out, ", (%d) <Rg>_z", count++);
  fprintf(out, ", (%d) <Rg^2>_n", count++);
  fprintf(out, ", (%d) <Rg^2>_w", count++);
  fprintf(out, ", (%d) <Rg^2>_z", count++);
  fprintf(out, ", (%d) <Anis>_n", count++);
  fprintf(out, ", (%d) <Acyl>_n", count++);
  fprintf(out, ", (%d) <Aspher>_n", count++);
  fprintf(out, ", (%d) <eigen.x>_n", count++);
  fprintf(out, ", (%d) <eigen.y>_n", count++);
  fprintf(out, ", (%d) <eigen.z>_n", count++);
  putc('\n', out);

  fclose(out); //}}}

  // open input aggregate file and read info from first line (Aggregates command) //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }

  // read minimum distance for closeness check (<distance> argument in Aggregates utility)
  double distance;
  fscanf(agg, "%*s %*s %lf", &distance);

  // skip <contacts> and <output.agg> in Aggregates command
  fscanf(agg, "%*s %*s");

  // read <type names> from Aggregates command //{{{
  // reading ends if next argument (beginning with '-') or the following empty line is read
  while ((test = getc(agg)) != '-' && test != '\n') {
    ungetc(test, agg);

    char name[LINE];
    fscanf(agg, "%s", name);
    int type = FindBeadType(name, Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\nError: %s - non-existent bead name '%s'\n", input_coor, name);
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }

    BeadType[type].Use = true;

    while ((test = getc(agg)) == ' ')
      ;
    ungetc(test, agg);
  } //}}}
  fclose(agg);

  // open again for production run - to ensure the pointer position in file is correct (at first 'Step')
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }

  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  Vector BoxLength = GetPBC(vcf, input_coor);

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // allocate memory for sum of various things //{{{
  // numbers of aggregates of all possibe sizes
  int *agg_counts_sum = calloc(Counts.Molecules,sizeof(int *));
  // total radius of gyration: [size][0] normal sum, [size][1] sum of Rg*mass, [size][2] Rg*mass^2
  double **Rg_sum = malloc(Counts.Molecules*sizeof(double *));
  // total square of radius of gyration: [size][0] normal sum, [size][1] sum of Rg^2*mass, [size][2] Rg^2*mass^2
  double **sqrRg_sum = malloc(Counts.Molecules*sizeof(double *));
  // relative shape anisotropy: only normal sum
  double *Anis_sum = calloc(Counts.Molecules,sizeof(double));
  // acylindricity: only normal sum
  double *Acyl_sum = calloc(Counts.Molecules,sizeof(double));
  // asphericity: only normal sum
  double *Aspher_sum = calloc(Counts.Molecules,sizeof(double));
  // gyration tensor eigenvalues
  struct Vector *eigen_sum = calloc(Counts.Molecules,sizeof(struct Vector));
  // total mass of aggregates: [size][0] normal sum, [size][1] sum of squares
  long int **mass_sum = malloc(Counts.Molecules*sizeof(int *));
  // number of molecule types in aggregates: [size][mol type] only normal sum
  int **molecules_sum = malloc(Counts.Molecules*sizeof(int *));
  for (int i = 0; i < Counts.Molecules; i++) {
    Rg_sum[i] = calloc(3,sizeof(double));
    sqrRg_sum[i] = calloc(3,sizeof(double));
    mass_sum[i] = calloc(2,sizeof(long int));
    molecules_sum[i] = calloc(Counts.TypesOfMolecules,sizeof(int));
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;

    // print steps? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}

    // read aggregates //{{{
    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      fprintf(stderr, "\nError: premature end of %s file (after %d. step - '%s')\n\n", input_agg, count, stuff);
      break;
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // allocate arrays for the timestep //{{{
    int *agg_counts_step = calloc(Counts.Molecules,sizeof(int));
    double **Rg_step = malloc(Counts.Molecules*sizeof(double *));
    double **sqrRg_step = malloc(Counts.Molecules*sizeof(double *));
    double *Anis_step = calloc(Counts.Molecules,sizeof(double));
    double *Acyl_step = calloc(Counts.Molecules,sizeof(double));
    double *Aspher_step = calloc(Counts.Molecules,sizeof(double));
    struct Vector *eigen_step = calloc(Counts.Molecules,sizeof(struct Vector));
    for (int i = 0; i < Counts.Molecules; i++) {
      Rg_step[i] = calloc(3,sizeof(double));
      sqrRg_step[i] = calloc(3,sizeof(double));
    } //}}}

    // calculate shape descriptors //{{{
    double mass_step[2] = {0}; // total mass of aggregates in a step: [0] normal, [1] sum of squares
    for (int i = 0; i < Counts.Aggregates; i++) {

      // test if aggregate 'i' should be used //{{{
      int size = 0;
      // agg size = number of molecules of type 'specific_moltype_for_size'
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type]) {
          size++;
        }
      }
      // is 'size' agg size in provided list?
      int correct_size = -1;
      for (int j = 0; j < Counts.Molecules; j++) {
        if ((j+1) == size) {
          correct_size = j;
        }
      } //}}}

      // if '-x' option is used, discount aggregates with only specified molecule type(s) //{{{
      test = false;
      for (int j = 0; j < size; j++) {
        int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
        if (MoleculeType[moltype].Write) {
          test = true; // a molecule that shouldn't be in agg 'i' is there
          break;
        }
      }
      if (!test) { // should the rest of the for loop agg i be skipped?
        exclude_count_chains += Aggregate[i].nMolecules;
        exclude_count_agg++;
        continue;
      } //}}}

      if (size != 0 && size >= range_As[0] && size <= range_As[1]) {
        if (correct_size != -1) {
          // copy bead ids to a separate array //{{{
          int *list = malloc(Aggregate[i].nBeads*sizeof(int));
          int n = 0;
          double agg_mass = 0;
          for (int j = 0; j < Aggregate[i].nBeads; j++) {
            int id = Aggregate[i].Bead[j];
            if (types_for_gyr[Bead[id].Type]) {
              list[n] = id;
              n++;
              agg_mass += BeadType[Bead[id].Type].Mass;
            }
          } //}}}

          Vector eigen = Gyration(n, list, Counts, BoxLength, BeadType, &Bead);

//        // calcule Rg the 'usual way' -- for testing purposes //{{{
//        double Rg2 = 0;
//        Vector test_com;
//        test_com.x = 0;
//        test_com.y = 0;
//        test_com.z = 0;
//        for (int j = 0; j < n; j++) {
//          Vector rij = Distance(Bead[list[j]].Position, com, BoxLength);
//          Rg2 += SQR(rij.x) + SQR(rij.y) + SQR(rij.z);
//          Rg2 += SQR(Bead[list[j]].Position.x) + SQR(Bead[list[j]].Position.y) + SQR(Bead[list[j]].Position.z);
//          test_com.x += Bead[list[j]].Position.x;
//          test_com.y += Bead[list[j]].Position.y;
//          test_com.z += Bead[list[j]].Position.z;
//        }
//        Rg2 /= n; //}}}

          free(list); // free array of bead ids for gyration calculation

          double Rgi = sqrt(eigen.x + eigen.y + eigen.z);

          if (eigen.x < 0 || eigen.y < 0 || eigen.z < 0) {
            fprintf(stderr, "Error: negative eigenvalues (%lf, %lf, %lf)\n\n", eigen.x, eigen.y, eigen.z);
          }
          // agg masses
          mass_step[0] += agg_mass; // for this timestep
          mass_step[1] += SQR(agg_mass); // for this timestep
          // radius of gyration
          Rg_step[correct_size][0] += Rgi; // for number avg
          Rg_step[correct_size][1] += Rgi * agg_mass; // for weight average
          Rg_step[correct_size][2] += Rgi * SQR(agg_mass); // for z-average
          // squared radius of gyration
          sqrRg_step[correct_size][0] += SQR(Rgi); // for number avg
          sqrRg_step[correct_size][1] += SQR(Rgi) * agg_mass; // for weight average
          sqrRg_step[correct_size][2] += SQR(Rgi) * SQR(agg_mass); // for z-average
          // relative shape anisotropy
          Anis_step[correct_size] += 1.5 * (SQR(eigen.x) + SQR(eigen.y) + SQR(eigen.z)) / SQR(eigen.x + eigen.y + eigen.z) - 0.5;
          // acylindricity
          Acyl_step[correct_size] += eigen.y - eigen.x;
          // asphericity
          Aspher_step[correct_size] += eigen.z - 0.5 * (eigen.x + eigen.y);
          // gyration vector eigenvalues
          eigen_step[correct_size].x += eigen.x;
          eigen_step[correct_size].y += eigen.y;
          eigen_step[correct_size].z += eigen.z;
          // aggregate count
          agg_counts_step[correct_size]++;

          // sum molecules and aggregates (if step between start and end)
          if (count >= start && (end == -1 || count <= end)) {
            agg_counts_sum[correct_size]++;
            // aggregate mass
            mass_sum[correct_size][0] += agg_mass;
            mass_sum[correct_size][1] += SQR(agg_mass);
            for (int j = 0; j < Aggregate[i].nMolecules; j++) {
              int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
              molecules_sum[correct_size][mol_type]++;
            }
          }
        }
      }
    } //}}}

    // add values to sums //{{{
    if (count >= start && (end == -1 || count <= end)) {
      for (int i = 0; i < Counts.Molecules; i++) {
        Rg_sum[i][0] += Rg_step[i][0];
        Rg_sum[i][1] += Rg_step[i][1];
        Rg_sum[i][2] += Rg_step[i][2];
        sqrRg_sum[i][0] += sqrRg_step[i][0];
        sqrRg_sum[i][1] += sqrRg_step[i][1];
        sqrRg_sum[i][2] += sqrRg_step[i][2];
        Anis_sum[i] += Anis_step[i];
        Acyl_sum[i] += Acyl_step[i];
        Aspher_sum[i] += Aspher_step[i];
        eigen_sum[i].x += eigen_step[i].x;
        eigen_sum[i].y += eigen_step[i].y;
        eigen_sum[i].z += eigen_step[i].z;
      }
    } //}}}

    // print data to output file //{{{
    FILE *out;
    if ((out = fopen(output, "a")) == NULL) { // out file opened fine? //{{{
      ErrorFileOpen(output, 'a');
      exit(1);
    } //}}}

    fprintf(out, "%5d", count); // timestep
    // sum up contributions from all aggregate sizes
    for (int i = 1; i < Counts.Molecules; i++) {
      Rg_step[0][0] += Rg_step[i][0];
      Rg_step[0][1] += Rg_step[i][1];
      Rg_step[0][2] += Rg_step[i][2];
      sqrRg_step[0][0] += sqrRg_step[i][0];
      sqrRg_step[0][1] += sqrRg_step[i][1];
      sqrRg_step[0][2] += sqrRg_step[i][2];
      Anis_step[0] += Anis_step[i];
      Acyl_step[0] += Acyl_step[i];
      Aspher_step[0] += Aspher_step[i];
      eigen_step[0].x += eigen_step[i].x;
      eigen_step[0].y += eigen_step[i].y;
      eigen_step[0].z += eigen_step[i].z;

      agg_counts_step[0] += agg_counts_step[i];
    }
    // <R_G>
    fprintf(out, " %8.5f %8.5f %8.5f", Rg_step[0][0]/agg_counts_step[0], Rg_step[0][1]/mass_step[0], Rg_step[0][2]/mass_step[1]);
    // <R_G^2>
    fprintf(out, " %8.5f %8.5f %8.5f", sqrRg_step[0][0]/agg_counts_step[0], sqrRg_step[0][1]/mass_step[0], sqrRg_step[0][2]/mass_step[1]);
    // relative shape anisotropy
    fprintf(out, " %8.5f", Anis_step[0]/agg_counts_step[0]);
    // acylindricity
    fprintf(out, " %8.5f", Acyl_step[0]/agg_counts_step[0]);
    // asphericity
    fprintf(out, " %8.5f", Aspher_step[0]/agg_counts_step[0]);
    // eigenvalues
    fprintf(out, " %8.5f %8.5f %8.5f", eigen_step[0].x/agg_counts_step[0], eigen_step[0].y/agg_counts_step[0], eigen_step[0].z/agg_counts_step[0]);
    putc('\n', out);

    fclose(out); //}}}

    // free memory //{{{
    free(agg_counts_step);
    for (int i = 0; i < Counts.Molecules; i++) {
      free(Rg_step[i]);
      free(sqrRg_step[i]);
    }
    free(Rg_step);
    free(sqrRg_step);
    free(Anis_step);
    free(Acyl_step);
    free(Aspher_step);
    free(eigen_step); //}}}
  }
  fclose(vcf);
  fclose(agg);

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

  // calculate per-size averages? //{{{
  if (per_size_file[0] != '\0') {

    // open file //{{{
    FILE *out;
    if ((out = fopen(per_size_file, "w")) == NULL) {
      ErrorFileOpen(per_size_file, 'w');
      exit(1);
    } //}}}

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++) {
      fprintf(out, " %s", argv[i]);
    }
    putc('\n', out); //}}}

    fprintf(out, "# column: (1) agg size");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      fprintf(out, " (%d) <%s>", i+2, MoleculeType[i].Name);
    }
    fprintf(out, " (%d) <Rg>, ", Counts.TypesOfMolecules+2);
    fprintf(out, " (%d) <Rg^2>, ", Counts.TypesOfMolecules+3);
    fprintf(out, " (%d) <Anis>, ", Counts.TypesOfMolecules+4);
    fprintf(out, " (%d) <Acyl>, ", Counts.TypesOfMolecules+5);
    fprintf(out, " (%d) <Aspher>, ", Counts.TypesOfMolecules+6);
    fprintf(out, " (%d) <eigen.x>, ", Counts.TypesOfMolecules+7);
    fprintf(out, " (%d) <eigen.y>, ", Counts.TypesOfMolecules+8);
    fprintf(out, " (%d) <eigen.z>, ", Counts.TypesOfMolecules+9);
    fprintf(out, " (%d) number of aggs", Counts.TypesOfMolecules+10);
    putc('\n', out);
    for (int i = 0; i < Counts.Molecules; i++) {
      if (agg_counts_sum[i] > 0) {
        fprintf(out, "%4d", i+1);
        for (int j = 0; j < Counts.TypesOfMolecules; j++) {
          fprintf(out, " %7.3f", (double)(molecules_sum[i][j])/agg_counts_sum[i]);
        }
        fprintf(out, " %7.3f", Rg_sum[i][0]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", sqrRg_sum[i][0]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", Anis_sum[i]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", Acyl_sum[i]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", Aspher_sum[i]/agg_counts_sum[i]);
        fprintf(out, " %7.3f", eigen_sum[i].x/agg_counts_sum[i]);
        fprintf(out, " %7.3f", eigen_sum[i].y/agg_counts_sum[i]);
        fprintf(out, " %7.3f", eigen_sum[i].z/agg_counts_sum[i]);
        fprintf(out, " %d", agg_counts_sum[i]);
        putc('\n', out);
      }
    }

    fclose(out);
  } //}}}

  // total averages //{{{
  for (int i = 1; i < Counts.Molecules; i++) {
    Rg_sum[0][0] += Rg_sum[i][0];
    Rg_sum[0][1] += Rg_sum[i][1];
    Rg_sum[0][2] += Rg_sum[i][2];
    sqrRg_sum[0][0] += sqrRg_sum[i][0];
    sqrRg_sum[0][1] += sqrRg_sum[i][1];
    sqrRg_sum[0][2] += sqrRg_sum[i][2];
    Anis_sum[0] += Anis_sum[i];
    Acyl_sum[0] += Acyl_sum[i];
    Aspher_sum[0] += Aspher_sum[i];
    eigen_sum[0].x += eigen_sum[i].x;
    eigen_sum[0].y += eigen_sum[i].y;
    eigen_sum[0].z += eigen_sum[i].z;

    agg_counts_sum[0] += agg_counts_sum[i];

    mass_sum[0][0] += mass_sum[i][0];
    mass_sum[0][1] += mass_sum[i][1];

    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      molecules_sum[0][j] += molecules_sum[i][j];
    }
  }

  // print to output file
  if ((out = fopen(output, "a")) == NULL) { //{{{
    ErrorFileOpen(output, 'a');
    exit(1);
  } //}}}

  fprintf(out, "# (1) <M>_n, (2) <M>_w ");
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, "(%d) <%s>, ", i+3, MoleculeType[i].Name);
  }
  fprintf(out, "(%d) <Rg>_n, ", Counts.TypesOfMolecules+3);
  fprintf(out, "(%d) <Rg>_w, ", Counts.TypesOfMolecules+4);
  fprintf(out, "(%d) <Rg>_z, ", Counts.TypesOfMolecules+5);
  fprintf(out, "(%d) <Rg^2>_n, ", Counts.TypesOfMolecules+6);
  fprintf(out, "(%d) <Rg^2>_w, ", Counts.TypesOfMolecules+7);
  fprintf(out, "(%d) <Rg^2>_z, ", Counts.TypesOfMolecules+8);
  fprintf(out, "(%d) <Anis>, ", Counts.TypesOfMolecules+9);
  fprintf(out, "(%d) <Acyl>, ", Counts.TypesOfMolecules+10);
  fprintf(out, "(%d) <Aspher>, ", Counts.TypesOfMolecules+11);
  fprintf(out, "(%d) <eigen.x>, ", Counts.TypesOfMolecules+12);
  fprintf(out, "(%d) <eigen.y>, ", Counts.TypesOfMolecules+13);
  fprintf(out, "(%d) <eigen.z>, ", Counts.TypesOfMolecules+14);
  putc('\n', out);
  fprintf(out, "# %lf", (double)(mass_sum[0][0])/agg_counts_sum[0]); //<M_As>_n
  fprintf(out, " %lf", (double)(mass_sum[0][1])/mass_sum[0][0]); //<M_As>_w
  // molecule types
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(out, " %lf", (double)(molecules_sum)[0][i]/agg_counts_sum[0]);
  }
  fprintf(out, " %lf", Rg_sum[0][0]/agg_counts_sum[0]); // <Rg>_n
  fprintf(out, " %lf", Rg_sum[0][1]/mass_sum[0][0]); // <Rg>_w
  fprintf(out, " %lf", Rg_sum[0][2]/mass_sum[0][1]); // <Rg>_z
  fprintf(out, " %lf", sqrRg_sum[0][0]/agg_counts_sum[0]); // <Rg^2>_n
  fprintf(out, " %lf", sqrRg_sum[0][1]/mass_sum[0][0]); // <Rg^2>_w
  fprintf(out, " %lf", sqrRg_sum[0][2]/mass_sum[0][1]); // <Rg^2>_z
  fprintf(out, " %lf", Anis_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %lf", Acyl_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %lf", Aspher_sum[0]/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0].x/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0].y/agg_counts_sum[0]);
  fprintf(out, " %lf", eigen_sum[0].z/agg_counts_sum[0]);
  putc('\n', out);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(Rg_sum[i]);
    free(sqrRg_sum[i]);
    free(mass_sum[i]);
    free(molecules_sum[i]);
  }
  free(molecules_sum);
  free(mass_sum);
  free(agg_counts_sum);
  free(Rg_sum);
  free(sqrRg_sum);
  free(Anis_sum);
  free(Acyl_sum);
  free(Aspher_sum);
  free(eigen_sum);
  free(stuff);
  free(specific_moltype_for_size); //}}}

  return 0;
}
