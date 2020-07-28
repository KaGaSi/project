#include "../AnalysisTools.h"

/*
 * TODO: move monomeric beads towards the aggregate
  // put monomeric beads in contact with their aggregates //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    // go through monomeric beads in the aggregate
    for (int j = 0; j < Aggregate[i].nMonomers; j++) {
      int id1 = Aggregate[i].Monomer[j]; // id_move_to = -1;
      // find smallest distance between the monomeric bead and bonded beads //{{{
      double min_dist = 1000000;
      for (int k = 0; k < Aggregate[i].nBeads; k++) {
        int id2 = Aggregate[i].Bead[k];
        if (BeadType[(*Bead)[id2].Type].Use) {
          VECTOR dist;
          dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
          dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
          dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
          double d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
          if (d < min_dist) {
            min_dist = d;
          }
          // stop if the monomeric bead is confirmed close to the aggregate
          if (min_dist <= distance) {
            break;
          }
        }
      } //}}}

      // move monomer if it's too far from aggregate //{{{
      while (min_dist > distance) {
        double d, min_dist_2;
        if (id1 == 10) {
          printf("before -x %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by -BoxLength.x //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.x -= BoxLength.x;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.x += BoxLength.x;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before +x %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by +BoxLength.x //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.x += BoxLength.x;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.x -= BoxLength.x;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before -y %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by -BoxLength.y //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.y -= BoxLength.y;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.y += BoxLength.y;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before +y %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by +BoxLength.y //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.y += BoxLength.y;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.y -= BoxLength.y;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before -z %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by -BoxLength.z //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.z -= BoxLength.z;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.z += BoxLength.z;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before +z %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by +BoxLength.z //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.z += BoxLength.z;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.z -= BoxLength.z;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("after +z %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }

        // find smallest distance between the monomeric bead and bonded beads
        min_dist = 1000000;
        for (int k = 0; k < Aggregate[i].nBeads; k++) {
          int id2 = Aggregate[i].Bead[k];
          VECTOR dist;
          dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
          dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
          dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
          double d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
          if (d < min_dist) {
            min_dist = d;
          }
        }
        if (id1 == 10) {
          printf("%3d %4d %lf (%lf, %lf, %lf)\n", i, id1, min_dist,
                                                  (*Bead)[id1].Position.x,
                                                  (*Bead)[id1].Position.y,
                                                  (*Bead)[id1].Position.z);
        }
        if (id1 == 10) {
          printf("%3d %4d %lf (%lf, %lf, %lf)\n", i, id1, min_dist,
                                                  (*Bead)[id1].Position.x,
                                                  (*Bead)[id1].Position.y,
                                                  (*Bead)[id1].Position.z);
        }
      } //}}}
    }
  } //}}}
printf("(%lf, %lf, %lf)\n", (*Bead)[10].Position.x,
                            (*Bead)[10].Position.y,
                            (*Bead)[10].Position.z);

  // test monomeric beads //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    for (int j = 0; j < Aggregate[i].nMonomers; j++) {
      int id1 = Aggregate[i].Monomer[j];
      int near;
      double min_dist = 1000000;
      VECTOR near_d;
      for (int k = 0; k < Aggregate[i].nBeads; k++) {
        int id2 = Aggregate[i].Bead[k];
        VECTOR dist;
        dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
        dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
        dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
        double d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
        if (d < min_dist) {
          min_dist = d;
          near = id2;
          near_d.x = dist.x;
          near_d.y = dist.y;
          near_d.z = dist.z;
        }
      }
      if (min_dist > distance) {
        printf("%3d %4d %4d %lf (%lf, %lf, %lf); (%lf, %lf, %lf)\n", i, id1, near, min_dist, near_d.x, near_d.y, near_d.z,
                                                                     (*Bead)[id1].Position.x,
                                                                     (*Bead)[id1].Position.y,
                                                                     (*Bead)[id1].Position.z);
      }
    }
  } //}}}
*/

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
VisualizeAgg writes specified aggregates into a new vcf file. This file does \
not have to contain all beads of any type, so it cannot be used for further \
analysis using AnalysisTools utilities.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <input.agg> <output> <agg size(s)> <options>\n\n", cmd);

  fprintf(ptr, "   <input>           input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <input.agg>       input agg file\n");
  fprintf(ptr, "   <output>          output vcf file(s) - one per aggregate size (automatic ending <size>.vcf)\n");
  fprintf(ptr, "   <agg size(s)>     aggregate size(s) to calculate density for\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined coordinates\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> molecule types in an aggregate\n");
  fprintf(ptr, "      -x <name(s)>   exclude aggregates containing only specified molecule(s)\n");
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
  int req_args = 4; //}}}

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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
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

  // <input.agg> - input agg file //{{{
  char input_agg[LINE];
  strcpy(input_agg, argv[++count]);

  // test if <input.agg> ends with '.agg'
  ext = 1;
  strcpy(extension[0], ".agg");
  if (ErrorExtension(input_agg, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output> - vcf file(s) for aggregates //{{{
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
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
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

  // vsf file is not needed anymore
  free(input_vsf);

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
  } //}}}

  // <agg sizes> - aggregate sizes to write //{{{
  int *agg_sizes = calloc(Counts.Molecules, sizeof(int));

  int aggs = 0;

  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (!IsInteger(argv[count])) {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}

    agg_sizes[aggs] = atoi(argv[count]);

    aggs++; // number of aggregate sizes
  } //}}}

  double distance; // <distance> parameter from Aggregate command
  int contacts; // <contacts> parameter from Aggregate command - not used here
  ReadAggCommand(BeadType, Counts, input_coor, input_agg, &distance, &contacts);

  // open input aggregate file and skip the first lines (Aggregate command & blank line) //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  char line[LINE];
  fgets(line, sizeof(line), agg);
  fgets(line, sizeof(line), agg); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  VECTOR BoxLength = GetPBC(vcf, input_coor);

  // write initial stuff to output density file //{{{
  for (int i = 0; i < aggs; i++) {
    FILE *out;
    char str[LINE];
    strcpy(str, output);

    char str2[LINE+6];
    sprintf(str2, "%s%d.vcf", str, agg_sizes[i]);
    strcpy(str, str2);
    if ((out = fopen(str, "w")) == NULL) {
      ErrorFileOpen(str, 'w');
      exit(1);
    }

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++)
      fprintf(out, " %s", argv[i]);
    putc('\n', out); //}}}

    // print agg size and periodic boundary conditions
    fprintf(out, "# aggregate size: %d\n", agg_sizes[i]);
    fprintf(out, "pbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

    fclose(out);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

  // initialize the array
  for (int i = 0; i < LINE; i++) {
    stuff[i] = '\0';
  } //}}}

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
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

    fprintf(stdout, "Chosen aggregate sizes:");
    for (int i = 0; i < aggs; i++) {
      fprintf(stdout, " %d", agg_sizes[i]);
    }
    putchar('\n');
  } //}}}

  // skip first start-1 steps //{{{
  count = 0;
  int test;
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}

    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_agg);
      fprintf(stderr, "\033[0m");
      test = '\0';
      exit(1);
    }

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_coor);
      fprintf(stderr, "\033[0m");
      exit(1);
    }
  }
  // print number of starting step? //{{{
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
  }
  //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int count_vcf = start - 1;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    count_vcf++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // read aggregates //{{{
    if (ReadAggregates(agg, &Counts, &Aggregate, BeadType, &Bead, MoleculeType, &Molecule, Index)) {
      if (!silent && !script) { // end of line if \r is used for printing step number
        putchar('\n');
      }
      count--; // because last step isn't processed
      count_vcf--;
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: premature end of \033[1;33m%s\033[1;31m file\n\n", input_agg);
      fprintf(stderr, "\033[0m");
      break;
    } //}}}

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count_vcf, stuff);
      exit(1);
    } //}}}

    // join agggregates if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // find correct aggregates and save them //{{{
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
      // is 'size' in provided aggregate sizes?
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j] == size) {
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
        continue;
      } //}}}

      if (correct_size != -1) {
        FILE *out;
        char str[LINE];
        strcpy(str, output);
        char str2[LINE+6];
        sprintf(str2, "%s%d.vcf", str, agg_sizes[correct_size]);
        strcpy(str, str2);
        if ((out = fopen(str, "a")) == NULL) {
          ErrorFileOpen(str, 'a');
          exit(1);
        }
        fprintf(out, "%s\nindexed\n", stuff);
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int id = Aggregate[i].Bead[j];
          fprintf(out, "%6d %7.3f %7.3f %7.3f\n", Bead[id].Index,
                                                  Bead[id].Position.x,
                                                  Bead[id].Position.y,
                                                  Bead[id].Position.z);
        }
        fclose(out);
      }
    } //}}}

    if (end == count_vcf)
      break;
  }
  fclose(vcf);
  fclose(agg);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count_vcf);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                          ");
      fprintf(stdout, "\rLast Step: %d\n", count_vcf);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  free(Index);
  free(agg_sizes);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  free(specific_moltype_for_size); //}}}

  return 0;
}
