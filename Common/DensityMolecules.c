#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <width> <output.rho> <molecule name(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <width>           width of a single bin\n");
  fprintf(stderr, "   <output.rho>      output density file (automatic ending 'molecule_name.rho' added)\n");
  fprintf(stderr, "   <molecule names>  molecule names to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined       specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -n <int>       number of bins to average\n");
  fprintf(stderr, "      -st <int>      starting timestep for calculation\n");
  fprintf(stderr, "      -c 'x's/<ints> use <int>-th molecule bead instead of centre of mass\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("\
MolDensity utility calculates number beads density for all bead types from the \
centre of mass (or specified bead number in a molecule) of specified molecules. \
Care must be taken with beadtype names in various molecules types, because if \
one beadtype appears in more molecule types, the resulting density for that \
beadtype will be averaged without regard for the various types of molecule it \
appears in.\n\n");

      printf("\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <width> <output.rho> <molecule name(s)> <options>\n\n", argv[0]);

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <width>           width of a single bin\n");
      printf("   <output.rho>      output density file (automatic ending 'molecule_name.rho' added)\n");
      printf("   <molecule names>  molecule names to calculate density for\n");
      printf("   <options>\n");
      printf("      --joined       specify that aggregates with joined coordinates are used\n");
      printf("      -n <int>       number of bins to average\n");
      printf("      -st <int>      starting timestep for calculation\n");
      printf("      -c 'x's/<ints> use <int>-th molecule bead instead of centre of mass\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 4; //}}}

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
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-n") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-c") != 0) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < 4) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of at least %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
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
  bool verbose2, silent;
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  VerboseLongOption(argc, argv, &verbose, &verbose2); // more verbose output
  SilentOption(argc, argv, &verbose, &verbose2, &silent); // no output
  bool script = BoolOption(argc, argv, "--script"); // do not use \r & co.
  // }}}

  // are provided coordinates joined? //{{{
  bool joined = BoolOption(argc, argv, "--joined"); //}}}

  // starting timestep //{{{
  int start = 1;
  if (IntegerOption(argc, argv, "-st", &start)) {
    exit(1);
  } //}}}

  // number of bins to average //{{{
  int avg = 1;
  if (IntegerOption(argc, argv, "-n", &avg)) {
    exit(1);
  } //}}}
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

  // <width> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <width>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double width = atof(argv[count]); //}}}

  // <output.rho> - filename with bead densities //{{{
  char output_rho[16];
  strcpy(output_rho, argv[++count]); //}}}

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

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {

    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);
    printf("%s\n", argv[count]);

    if (mol_type == -1) {
      fprintf(stderr, "Molecule '%s' does not exist in FIELD!\n", argv[count]);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }

    // write initial stuff to output density file //{{{
    FILE *out;
    char str[32];

    sprintf(str, "%s%s.rho", output_rho, argv[count]);
    if ((out = fopen(str, "w")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", str);
      exit(1);
    }

    // print command to output file //{{{
    putc('#', out);
    for (int i = 0; i < argc; i++)
      fprintf(out, " %s", argv[i]);
    putc('\n', out); //}}}

    // print bead type names to output file //{{{
    putc('#', out);
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      fprintf(out, " %s", BeadType[i].Name);
    }
    putc('\n', out); //}}}

    fclose(out); //}}}
  } //}}}

  // -c 'x's/<int(s)> option - specify which bead to use as a molecule centre //{{{

  // array for considering whether to use COM or specified bead number //{{{
  int centre[Counts.TypesOfMolecules];
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    centre[i] = -1; // use COM for this molecule type
  } //}}}

  int used_mols = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-c") == 0) {
      int j = 0;
      while ((++j+i) < argc && argv[i+j][0] != '-') {
        printf("%s\n", argv[i+j]);

        // Error - non-numeric or non-"x" argument //{{{
        if ((argv[i+j][0] < '0' || argv[i+j][0] > '9') && argv[i+j][0] != 'x') {
          fprintf(stderr, "Wrong argement for '-c' option!\n");
          ErrorHelp(argv[0]);
          exit(1);
        } //}}}

        if (argv[i+j][0] == 'x') {
          centre[j-1] = -1;
        } else {

          int k;
          for (k = 0; k < Counts.TypesOfMolecules && used_mols < j; k++) {
            if (MoleculeType[k].Use) {
              used_mols++;
            }
          }
          printf("k = %d; used_mols = %d\n", k, used_mols);

          centre[used_mols-1] = atoi(argv[i+j]) - 1;

          // Error - too high bead number //{{{
          if (centre[used_mols-1] > MoleculeType[j-1].nBeads) {
            fprintf(stderr, "Incorrect number in '-c' option (%dth bead in molecule"
              "%s containing only %d beads)\n", centre[j-1]+1, MoleculeType[j-1].Name,
              MoleculeType[j-1].nBeads);
            exit(1);
          } //}}}
        }
      }
    }
  } //}}}

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

  // number of bins //{{{
  double max_dist = 0.5 * Min3(BoxLength.x, BoxLength.y, BoxLength.z);
  int bins = ceil(max_dist / width); //}}}

  // allocate memory for density arrays //{{{
  double ***rho = malloc(Counts.TypesOfBeads*sizeof(double **));
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    rho[i] = malloc(Counts.TypesOfMolecules*sizeof(double *));
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      rho[i][j] = calloc(bins,sizeof(double));
    }
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

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
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    printf("Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        printf(" %s", MoleculeType[i].Name);
      }
    }
    putchar('\n');
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

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
  //}}}

  // main loop //{{{
  count = 0; // count timesteps
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

    // join molecules if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // calculate densities //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type_i = Molecule[i].Type;
      if (MoleculeType[mol_type_i].Use) {

        // determine centre to calculate densities from //{{{
        Vector com;
        if (centre[mol_type_i] == -1 ) { // use molecule's centre of mass
          com = CentreOfMass(MoleculeType[mol_type_i].nBeads, Molecule[i].Bead, Bead, BeadType);
        } else { // use centre[mol_type_i]-th molecule's bead as com
          com.x = Bead[Molecule[i].Bead[centre[mol_type_i]]].Position.x;
          com.y = Bead[Molecule[i].Bead[centre[mol_type_i]]].Position.y;
          com.z = Bead[Molecule[i].Bead[centre[mol_type_i]]].Position.z;
        } //}}}

        // molecule beads //{{{
        for (int j = 0; j < MoleculeType[mol_type_i].nBeads; j++) {
          int bead_j = Molecule[i].Bead[j];

          Vector dist = Distance(Bead[bead_j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            rho[Bead[bead_j].Type][mol_type_i][k]++;
          }
        } //}}}

        // beads from other molecules //{{{
        for (int j = 0; j < Counts.Molecules; j++) {
          int mol_type_j = Molecule[j].Type;

          if (strcmp(MoleculeType[mol_type_i].Name, MoleculeType[mol_type_j].Name) != 0) {
            for (int k = 0; k < MoleculeType[mol_type_j].nBeads; k++) {
              int bead_k = Molecule[j].Bead[k];

              Vector dist = Distance(Bead[bead_k].Position, com, BoxLength);
              dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

              if (dist.x < max_dist) {
                int l = dist.x / width;

                rho[Bead[bead_k].Type][mol_type_i][l]++;
              }
            }
          }
        } //}}}

        // monomeric beads //{{{
        for (int j = 0; j < Counts.Unbonded; j++) {
          Vector dist = Distance(Bead[j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            rho[Bead[j].Type][mol_type_i][k]++;
          }
        } //}}}
      }
    } //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
    }
  } //}}}

  // write densities to output file(s) //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      FILE *out;
      char str[32];

      sprintf(str, "%s%s.rho", output_rho, MoleculeType[i].Name);
      if ((out = fopen(str, "a")) == NULL) {
        fprintf(stderr, "Cannot open file %s!\n", str);
        exit(1);
      }

      // calculate rdf
      for (int j = 0; j < bins; j++) {

        // calculate volume of every shell that will be averaged
        double shell[avg];
        for (int k = 0; k < avg; k++) {
          shell[k] = 4 * PI * CUBE(width) *(CUBE(j+k+1) - CUBE(j+k)) / 3;
        }

        fprintf(out, "%.2f", width*(j+0.5*avg));

        for (int k = 0; k < Counts.TypesOfBeads; k++) {
          double temp = 0;

          // sum rdfs from all shells to be averaged
          for (int l = 0; l < avg; l++) {
            temp += rho[k][i][j+l] / (shell[l] * MoleculeType[i].Number * count);
          }

          // print average value to output file
          fprintf(out, " %10f", temp/avg);
        }
        putc('\n',out);
      }

      fclose(out);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      free(rho[i][j]);
    }
    free(rho[i]);
  }
  free(rho); //}}}

  return 0;
}
