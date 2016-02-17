#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <output file> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <output file>     name of output file with persistence length data\n");
  fprintf(stderr, "   <molecule names>  names of molecule type(s) to use for calculation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -i <name>      use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>      file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v             verbose output\n");
  fprintf(stderr, "      -V             verbose output with comments from input .vcf file\n");
  fprintf(stderr, "      -h             print this help and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("PersistenceLength utility calculates persistence length of linear chains (no    \n");
      printf("check whether the molecules are linear is performed). It calculates distance    \n");
      printf("between first and last bead in a molecule The utility uses dl_meso.vsf (or      \n");
      printf("other input structure file) and FIELD (along with optional bond file) files     \n");
      printf("to determine all information about the system.                                \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <output file> <molecule names> <options>\n\n", argv[0]);

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <output file>     name of output file with persistence length data\n");
      printf("   <molecule names>  names of molecule type(s) to use for calculation\n");
      printf("   <options>\n");
      printf("      -i <name>      use input .vsf file different from dl_meso.vsf\n");
      printf("      -b <name>      file containing bond alternatives to FIELD\n");
      printf("      -v             verbose output\n");
      printf("      -V             verbose output with comments from input .vcf file\n");
      printf("      -h             print this help and exit\n");
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

  // <output> - file name with end-to-end distances //{{{
  char output[32];
  strcpy(output, argv[++count]); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // <molecule names> - names of molecule types to use //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindMoleculeType(argv[count], Counts, MoleculeType);

    if (type == -1) {
      fprintf(stderr, "Molecule type '%s' is not in %s coordinate file!\n", argv[count], input_vcf);
      exit(1);
    }

    MoleculeType[type].Use = true;
  } //}}}

  // open output file and print molecule names //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output);
    exit(1);
  }

  fprintf(out, "#timestep ");

  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, "%10s", MoleculeType[i].Name);
    }
  }
  putc('\n', out);

  fclose(out); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
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

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = calloc(128,sizeof(int)); //}}}

  // main loop //{{{
  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    fflush(stdout);
    printf("\rStep: %6d", ++count);

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

    // join all molecules
    RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

    // calculate & write persistence length //{{{
    double bond_dist[Counts.TypesOfMolecules];
    double *cos_angle[Counts.TypesOfMolecules];
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      bond_dist[i] = 0;

      cos_angle[i] = calloc(MoleculeType[i].nBonds, sizeof(double));
    }

    // go through all molecules //{{{
    for (int i = 0; i < Counts.Molecules; i++) {
      int type = Molecule[i].Type;

      // use only specified molecules
      if (MoleculeType[type].Use) {

        // go through bond pairs in molecule 'i'
        for (int j = 0; j < MoleculeType[type].nBonds; j++) {
          // beads in first bond //{{{
          int id_j0 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
          int id_j1 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]]; //}}}

          // first bond vector //{{{
          Vector bond_j;
          bond_j.x = Bead[id_j0].Position.x - Bead[id_j1].Position.x;
          bond_j.y = Bead[id_j0].Position.y - Bead[id_j1].Position.y;
          bond_j.z = Bead[id_j0].Position.z - Bead[id_j1].Position.z; //}}}

          // bond length
          double dist_j = sqrt(SQR(bond_j.x) + SQR(bond_j.y) + SQR(bond_j.z));

          // add bond length for averaging
          bond_dist[type] += dist_j;

          for (int k = (j+1); k < MoleculeType[type].nBonds; k++) {
            // beads in second bond //{{{
            int id_k0 = Molecule[i].Bead[MoleculeType[type].Bond[k][0]];
            int id_k1 = Molecule[i].Bead[MoleculeType[type].Bond[k][1]]; //}}}

            // second bond vector //{{{
            Vector bond_k;
            bond_k.x = Bead[id_k0].Position.x - Bead[id_k1].Position.x;
            bond_k.y = Bead[id_k0].Position.y - Bead[id_k1].Position.y;
            bond_k.z = Bead[id_k0].Position.z - Bead[id_k1].Position.z; //}}}

            // bond length
            double dist_k = sqrt(SQR(bond_k.x) + SQR(bond_k.y) + SQR(bond_k.z));

            // dot product of two vectors
            double dot = bond_j.x * bond_k.x + bond_j.y * bond_k.y + bond_j.z * bond_k.z;

            // cos(angle) between two vectors
            cos_angle[type][k-j-1] += dot / (dist_j * dist_k);
          }
        }
      }
    } //}}}

    // open output file for appending //{{{
    if ((out = fopen(output, "a")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", output);
      exit(1);
    } //}}}

    // write persistence length to output file //{{{
    fprintf(out, "%6d", count);
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {

      // use only specified molecules
      if (MoleculeType[i].Use) {

        double persistence = 0;
        for (int j = 0; j < MoleculeType[i].nBonds; j++) {
          persistence += cos_angle[i][j] / (MoleculeType[i].nBonds - j);
        }
        persistence *= bond_dist[i] / (MoleculeType[i].nBonds * MoleculeType[i].Number);

        fprintf(out, "%10.5f", persistence/MoleculeType[i].Number);
      }
    }
    putc('\n', out); //}}}

    fclose(out); //}}}

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      printf("\n%s", stuff);

    // free calloc'd memory
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      free(cos_angle[i]);
    }
  }

  fflush(stdout);
  printf("\rLast Step: %6d\n", count);

  fclose(vcf); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      free(MoleculeType[i].Bond[j]);
    }
    free(MoleculeType[i].Bond);
  }
  free(MoleculeType);
  free(Bead);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(Molecule[i].Bead);
  }
  free(Molecule);
  free(stuff);
  //}}}

  return 0;
}