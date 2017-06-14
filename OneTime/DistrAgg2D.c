#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input> <output distr file> <options>\n\n", cmd);

  fprintf(stderr, "   <input>              input filename (agg format)\n");
  fprintf(stderr, "   <output distr file>  filename with weight and number distributions\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -i <name>         use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>         file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v                verbose output\n");
  fprintf(stderr, "      -V                verbose output with more information\n");
  fprintf(stderr, "      -h                print this help and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("                                                                            \n");
      printf("                                                                            \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input> <output distr file> <options>\n\n", argv[0]);

      printf("   <input>              input filename (agg format)\n");
      printf("   <output distr file>  filename with weight and number distributions\n");
      printf("   <options>\n");
      printf("      -i <name>         use input .vsf file different from dl_meso.vsf\n");
      printf("      -b <name>         file containing bond alternatives to FIELD\n");
      printf("      -v                verbose output\n");
      printf("      -V                verbose output with more information\n");
      printf("      -h                print this help and exit\n");
      exit(0);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  printf("\n\n"); //}}}

  // check if correct number of arguments //{{{
  if (argc < 4) {
    fprintf(stderr, "Too little arguments!\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -i <name> option - filename of input structure file //{{{
  char vsf_file[32];
  vsf_file[0] = '\0'; // check if -i option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");
        exit(1);
      }

      // check if .vsf ending is present
      char *vsf = strrchr(argv[i+1], '.');
      if (!vsf || strcmp(vsf, ".vsf")) {
        fprintf(stderr, "'-i' arguments does not have .vsf ending!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(vsf_file, argv[i+1]);
    }
  }

  // -i option is not used
  if (vsf_file[0] == '\0') {
    strcpy(vsf_file, "dl_meso.vsf");
  } //}}}

  // -b <name> option - filename of input bond file //{{{
  char bonds_file[32];
  bonds_file[0] = '\0'; // check if -b option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-b") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-b' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(bonds_file, argv[i+1]);
    }
  } //}}}

  // -v option - verbose output //{{{
  bool verbose = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      verbose = true;

      break;
    }
  } //}}}

  // -V option - verbose output with comments from input .vcf file //{{{
  bool verbose2 = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-V") == 0) {
      verbose = true;
      verbose2 = true;

      break;
    }
  } //}}}

  int count = 0; // count mandatory arguments

  // <input> - filename of input agg file //{{{
  char input_agg[32];
  strcpy(input_agg, argv[++count]); //}}}

  // open input file and skip the first two lines //{{{
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_agg);
    exit(1);
  }

  while (getc(agg) != '\n')
    ;
  while (getc(agg) != '\n')
    ; //}}}

  // <output distr file> - filename with weight and number distributions //{{{
  char output_distr[32];
  strcpy(output_distr, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  char vcf[1];
  vcf[0] = '\0';
  ReadStructure(vsf_file, vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));

  for (int i = 0; i < Counts.Molecules; i++) {
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
  } //}}}

  // open output files and print first line //{{{
  FILE *out;
  if ((out = fopen(output_distr, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  fprintf(out, "# N_Homo N_Diblock wdistr\n");
  fclose(out); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    char null[1] = {'\0'};
    VerboseOutput(verbose2, null, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // arrays for distribution //{{{
  double wdistr[Counts.Molecules][Counts.Molecules];
  int HomosForDiblocks[Counts.Molecules][4];

  // zeroize arrays
  for (int i = 0; i < Counts.Molecules; i++) {
    HomosForDiblocks[i][0] = 0; // number of aggregates
    HomosForDiblocks[i][1] = 0; // number of homopolymers
    HomosForDiblocks[i][2] = 0; // CI+
    HomosForDiblocks[i][3] = 0; // CI-
    for (int j = 0; j < Counts.Molecules; j++) {
      wdistr[i][j] = 0;
    }
  } //}}}

  // main loop //{{{
  int test;
  count = 0;
  while ((test = getc(agg)) != 'L') { // cycle ends with 'Last Step' line in agg file
    ungetc(test, agg);

    fflush(stdout);
    printf("\rStep: %6d", ++count);

    ReadAggregates(agg, &Counts, &Aggregate, MoleculeType, Molecule);

    // print info about aggregates if '-V' is used //{{{
    if (verbose2) {
      for (int i = 0; i < Counts.Aggregates; i++) {
        printf("\nAggregate[%3d].{Mass = %6.2f,\nnMolecules = %3d:", i+1, Aggregate[i].Mass, Aggregate[i].nMolecules);
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          printf(" %d", Aggregate[i].Molecule[j]+1);
        }
        printf(",\n nBeads = %4d:", Aggregate[i].nBeads);
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          printf(" %d", Aggregate[i].Bead[j]);
        }
        printf(",\n nMonomers = %4d:", Aggregate[i].nMonomers);
        for (int j = 0; j < Aggregate[i].nMonomers; j++) {
          printf(" %d", Aggregate[i].Monomer[j]);
        }
        printf("}\n");
      }
      putchar('\n');
    } //}}}

    // go through all aggregates
    for (int i = 0; i < Counts.Aggregates; i++) {
      int Diblock = 0, Homo = 0;

      // go through all molecules in an aggregate
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name,"Diblock") == 0) {
          Diblock++;
        } else if (strcmp(MoleculeType[Molecule[Aggregate[i].Molecule[j]].Type].Name,"Homo") == 0) {
          Homo++;
        }
      }

      if ((Diblock+Homo) != Aggregate[i].nMolecules) { //{{{
        fprintf(stderr, "Error: Aggregate[%d].nMolecules = %d\n", i, Aggregate[i].nMolecules);
        fprintf(stderr, "Error: Diblock + Homo = %d + %d = %d\n", Diblock, Homo, Diblock+Homo);

        exit(1);
      } //}}}

      // distribution
      wdistr[Diblock][Homo] += Aggregate[i].nMolecules;

      HomosForDiblocks[Diblock][0]++;
      HomosForDiblocks[Diblock][1] += Homo;

      // find number of CI's near aggregate core
      for (int j = 0; j < Aggregate[i].nMonomers; j++) {
        if (BeadType[Bead[Aggregate[i].Monomer[j]].Type].Charge != 0) { // use only charged monomers
          bool near = false;

          for (int k = 0; k < Aggregate[i].nBeads && !near; k++) {
            if (BeadType[Bead[Aggregate[i].Bead[k]].Type].Charge != 0) { // use only charged aggregate beads

              Vector box;
              box.x = 30;
              box.y = 30;
              box.z = 30;

              // distance between two beads
              Vector rij = Distance(Bead[Aggregate[i].Monomer[j]].Position, Bead[Aggregate[i].Bead[k]].Position, box);

              rij.x = sqrt(SQR(rij.x) + SQR(rij.y) + SQR(rij.z));

              if (rij.x < 1) {
                near = true;
                if (BeadType[Bead[Aggregate[i].Monomer[j]].Type].Charge > 0) {
                  HomosForDiblocks[Diblock][2]++;
                } else {
                  HomosForDiblocks[Diblock][3]++;
                }
              }
            }
          }
        }
      }
    }
  }
  fclose(agg);

  fflush(stdout);
  printf("\rLast Step: %6d\n", count); //}}}

  // print distributions to output file //{{{
  if ((out = fopen(output_distr, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  int data = 80;
  for (int i = 0; i < data; i++) {
    for (int j = 0; j <= i; j++) {
      fprintf(out, "%4d %4d %lf\n", j, j+data-i, (double)(wdistr[j][j+data-i])/(Counts.Molecules*count));
    }
    putc('\n', out);
  }
  for (int i = 0; i <= data; i++) {
    fprintf(out, "%4d %4d %lf\n", i, i, (double)(wdistr[i][i])/(Counts.Molecules*count));
  }
  putc('\n', out);
  for (int i = 0; i < data; i++) {
    for (int j = 0; j <= i; j++) {
      fprintf(out, "%4d %4d %lf\n", j+data-i, j, (double)(wdistr[j+data-i][j])/(Counts.Molecules*count));
    }
    putc('\n', out);
  }
  fclose(out); //}}}

  if ((out = fopen("another_distr.txt", "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_distr);
    exit(1);
  }

  fprintf(out, "# N_Diblock <N_Homo> <N_CI+> <N_CI->\n");

  for (int i = 1; i <= data; i++) {
    if (HomosForDiblocks[i][0] != 0) {
      fprintf(out, "%3d %lf %lf %lf\n", i, (double)(HomosForDiblocks[i][1])/HomosForDiblocks[i][0], (double)(HomosForDiblocks[i][2])/HomosForDiblocks[i][0], (double)(HomosForDiblocks[i][3])/HomosForDiblocks[i][0]);
    }
  }

  fclose(out);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead); //}}}

  return 0;
}
