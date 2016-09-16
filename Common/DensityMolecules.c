#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <width> <output.rho> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <width>           width of a single bin\n");
  fprintf(stderr, "   <output.rho>      output density file (automatic ending 'agg#.rho' added)\n");
  fprintf(stderr, "   <molecule names>  molecule names to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -j             specify that aggregates with joined coordinates are used\n");
  fprintf(stderr, "      -n <average>   number of bins to average\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("MolDensity utility calculates number beads density for all bead types from  \n");
      printf("the center of mass of specified molecules. Care must be taken with beadtype \n");
      printf("names in various molecules types, because if one beadtype appears in more   \n");
      printf("molecule types, the resulting density for that beadtype will be averaged    \n");
      printf("without regard for the various types of molecule it appears in.             \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <width> <output.rho> <molecule names> <options>\n\n", argv[0]);

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <width>           width of a single bin\n");
      printf("   <output.rho>      output density file (automatic ending 'agg#.rho' added)\n");
      printf("   <molecule names>  molecule names to calculate density for\n");
      printf("   <options>\n");
      printf("      -j             specify that aggregates with joined coordinates are used\n");
      printf("      -n <average>   number of bins to average\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 4; //}}}

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

  // -j option - coordinates are joined //{{{
  bool joined = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {
      joined = true;
    }
  } //}}}

  // -n option - number of bins to average //{{{
  int avg = 1;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0) {

      // Error - non-numeric argument
      if (argv[i+1][0] < '0' || argv[i+1][0] > '9') {
        fprintf(stderr, "Non-numeric argement for '-n' option!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      avg = atoi(argv[i+1]);
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
      fprintf(stderr, "Non-existent molecule name:%s!\n", argv[count]);
      exit(1);
    } //}}}

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

  // main loop //{{{
  count = 0; // count timesteps
  int test;
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
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
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
      if (MoleculeType[Molecule[i].Type].Use) {
        Vector com = CenterOfMass(MoleculeType[Molecule[i].Type].nBeads, Molecule[i].Bead, Bead, BeadType);

        // molecule beads
        for (int j = 0; j < MoleculeType[Molecule[i].Type].nBeads; j++) {
          Vector dist = Distance(Bead[Molecule[i].Bead[j]].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            rho[Bead[Molecule[i].Bead[j]].Type][Molecule[i].Type][k]++;
          }
        }

        // beads from other molecules
        for (int j = 0; j < Counts.Molecules; j++) {
          if (strcmp(MoleculeType[Molecule[i].Type].Name, MoleculeType[Molecule[j].Type].Name) != 0) {
            for (int k = 0; k < MoleculeType[Molecule[j].Type].nBeads; k++) {
              Vector dist = Distance(Bead[Molecule[j].Bead[k]].Position, com, BoxLength);
              dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

              if (dist.x < max_dist) {
                int l = dist.x / width;

                rho[Bead[Molecule[j].Bead[k]].Type][Molecule[i].Type][l]++;
              }
            }
          }
        }

        // monomeric beads
        for (int j = 0; j < Counts.Unbonded; j++) {
          Vector dist = Distance(Bead[j].Position, com, BoxLength);
          dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

          if (dist.x < max_dist) {
            int k = dist.x / width;

            rho[Bead[j].Type][Molecule[i].Type][k]++;
          }
        }
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
