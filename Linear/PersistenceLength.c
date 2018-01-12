#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <output> <molecule names> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <output>          name of output file with persistence length data\n");
  fprintf(stderr, "   <molecule names>  names of molecule type(s) to use for calculation\n");
  fprintf(stderr, "   <options>\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
WARNING! Probably does not work as it should.\n\n");

      fprintf(stdout, "\
PersistenceLength utility calculates persistence length of linear chains (no \
check whether the molecules are linear is performed). It calculates distance \
between first and last bead in a molecule.\n\n");

      fprintf(stdout, "\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input.vcf> <output> <molecule names> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input.vcf>       input filename (vcf format)\n");
      fprintf(stdout, "   <output>          name of output file with persistence length data\n");
      fprintf(stdout, "   <molecule names>  names of molecule type(s) to use for calculation\n");
      fprintf(stdout, "   <options>\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 3; //}}}
  fprintf(stdout, "WARNING! Probably does not work as it should.\n\n");

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count][0] != '-'; i++) {
    count++;
  }

  if (count < options) {
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
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
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

  // vsf file is not needed anymore
  free(vsf_file);

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
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = calloc(128,sizeof(int)); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // allocate memory for sums of persistence length
  double *sum_persistence = calloc(Counts.TypesOfMolecules,sizeof(double));

  // main loop //{{{
  int test;
  count = 0;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      if (script) {
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
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
        sum_persistence[i] += persistence; // overall sum

        fprintf(out, "%10.5f", persistence/MoleculeType[i].Number);
      }
    }
    putc('\n', out); //}}}

    fclose(out); //}}}

    // if -V option used, print comment at the beginning of a timestep
    if (verbose2)
      fprintf(stdout, "\n%s", stuff);

    // free calloc'd memory
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      free(cos_angle[i]);
    }
  }

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  }

  fclose(vcf); //}}}

  // write overall averages to output file //{{{
  // open output file for appending //{{{
  if ((out = fopen(output, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output);
    exit(1);
  } //}}}

  putc('#', out);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      fprintf(out, " %s", MoleculeType[i].Name);
    }
  }
  putc('\n', out);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].Use) {
      sum_persistence[i] /= count * MoleculeType[i].Number;

      fprintf(out, " %lf", sum_persistence[i]);
    }
  }
  putc('\n', out);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  //}}}

  return 0;
}
