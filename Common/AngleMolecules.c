#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <molecule(s)> <options>\n\n", cmd);

  fprintf(stderr, "   <input.vcf>            input filename (vcf format)\n");
  fprintf(stderr, "   <molecule(s)>          molecule names to calculate density for\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      --joined            specify that molecules with joined coordinates are used\n");
  fprintf(stderr, "      -a <name> <ints>    filename and bead ids for angle calculation ");
  fprintf(stderr, "(default: 1 2 3 for each molecule type and filename 'angle.txt')\n");
  fprintf(stderr, "      -st <int>           starting timestep for calculation\n");
  CommonHelp(1);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
AngleMolecules utility calculates angles \
for specified beads (three beads per angel) in each molecule of specified \
molecule type(s). The specified beads do not to be connected by bonds. \
\n\n");

/*      fprintf(stdout, "\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");
*/

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <input.vcf> <molecule(s)> <options>\n\n", argv[0]);

      fprintf(stdout, "   <input.vcf>            input filename (vcf format)\n");
      fprintf(stdout, "   <molecule(s)>          molecule names to calculate density for\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      --joined            specify that molecules with joined coordinates are used\n");
      fprintf(stdout, "      -a <name> <ints>    filename and bead ids for angle calculation ");
      fprintf(stdout, "(default: 1 2 3 for each molecule type and filename 'angle.txt')\n");
      fprintf(stdout, "      -st <int>           starting timestep for calculation\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int req_args = 2; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
//      strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-a") != 0 &&
        strcmp(argv[i], "-st") != 0) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // options before reading system data //{{{
  // use .vsf file other than dl_meso.vsf? //{{{
  char *input_vsf = calloc(32,sizeof(char *));
  if (VsfFileOption(argc, argv, &input_vsf)) {
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
    ErrorExtension(input_vcf, ".vcf");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <molecule names> - types of molecules for calculation //{{{
  while (++count < argc && argv[count][0] != '-') {

    int mol_type = FindMoleculeType(argv[count], Counts, MoleculeType);
    fprintf(stdout, "%s\n", argv[count]);

    if (mol_type == -1) {
      fprintf(stderr, "Error: molecule '%s' does not exist in FIELD\n\n", argv[count]);
      exit(1);
    } else {
      MoleculeType[mol_type].Use = true;
    }
  } //}}}

  // '-a' option - specify for bead ids //{{{
  int bead[100] = {0}, number_of_beads = 0;
  char output[32];
  strcpy(output, "angle.txt"); // default output filename (if '-a' is missing)
  if (HundredIntegerOption(argc, argv, "-a", bead, &number_of_beads, output)) {
    exit(1);
  }
  // Error: wrong number of integers //{{{
  if ((number_of_beads%3) != 0) {
    fprintf(stderr, "\nError: '-a' option - number of bead ids must be divisible by three.\n");
    exit(1);
  } //}}}
  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero
    // Warning - too high id for specific molecule //{{{
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      if (MoleculeType[j].Use && bead[i] >= MoleculeType[j].nBeads) {
        fprintf(stderr, "\nWarning: '-a' option - %d is larger than the number of beads in molecule %s.\n", bead[i], MoleculeType[j].Name);
        fprintf(stderr, "Series with this bead id will not be used to calculate angle in %s.\n\n", MoleculeType[j].Name);
      }
    } //}}}
  } //}}}

  // open output file and write initial stuff //{{{
  FILE *out;
  if ((out = fopen(output, "w")) == NULL) {
    ErrorFileOpen(output, 'w');
    exit(1);
  }

  // print command to output file //{{{
  putc('#', out);
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out); //}}}

  // print molecule names & bead ids //{{{
  count = 1;
  fprintf(out, "# (1) step;");
  for (int i = 0; i < Counts.Molecules; i++) {
    int mol_type_i = Molecule[i].Type;
    if (MoleculeType[mol_type_i].Use) {
      for (int j = 0; j < number_of_beads; j += 3) {
        if (bead[j+0] < MoleculeType[mol_type_i].nBeads &&
            bead[j+1] < MoleculeType[mol_type_i].nBeads &&
            bead[j+2] < MoleculeType[mol_type_i].nBeads) {
          fprintf(out, " (%d) %d %s:%d-%d-%d;", ++count, i+1, MoleculeType[mol_type_i].Name, bead[j]+1, bead[j+1]+1, bead[j+2]+1);
        }
      }
    }
  }
  putc('\n', out); //}}}

  fclose(out); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    ErrorFileOpen(input_vcf, 'r');
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[128];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "\nError: cannot read a string from '%s' file\n\n", input_vcf);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "\nError: cannot read pbc from %s\n\n", input_vcf);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    fprintf(stdout, "   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
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

    fprintf(stdout, "Chosen molecule types:");
    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      if (MoleculeType[i].Use) {
        fprintf(stdout, " %s", MoleculeType[i].Name);
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
        fprintf(stdout, "Discarding step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rDiscarding step: %6d", count);
      }
    } //}}}

    if (SkipCoor(vcf, Counts, &stuff)) {
      fprintf(stderr, "\nError: premature end of %s file\n\n", input_vcf);
      exit(1);
    }
  }
  // print number of discarded steps? //{{{
  if (!silent) {
    if (script) {
      fprintf(stdout, "Discarded steps: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rDiscarded steps: %6d\n", count);
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
        fprintf(stdout, "Step: %6d\n", count);
      } else {
        fflush(stdout);
        fprintf(stdout, "\rStep: %6d", count);
      }
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_vcf, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    // join molecules if un-joined coordinates provided //{{{
    if (!joined) {
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
    } //}}}

    // calculate angles //{{{
    double angle[Counts.Molecules][number_of_beads/3];
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type_i = Molecule[i].Type;
      if (MoleculeType[mol_type_i].Use) {

        for (int j = 0; j < number_of_beads; j += 3) {
          if (bead[j+0] < MoleculeType[mol_type_i].nBeads &&
              bead[j+1] < MoleculeType[mol_type_i].nBeads &&
              bead[j+2] < MoleculeType[mol_type_i].nBeads) {
            Vector u, v;
            u.x = Bead[Molecule[i].Bead[bead[j+0]]].Position.x - Bead[Molecule[i].Bead[bead[j+1]]].Position.x;
            u.y = Bead[Molecule[i].Bead[bead[j+0]]].Position.y - Bead[Molecule[i].Bead[bead[j+1]]].Position.y;
            u.z = Bead[Molecule[i].Bead[bead[j+0]]].Position.z - Bead[Molecule[i].Bead[bead[j+1]]].Position.z;
            v.x = Bead[Molecule[i].Bead[bead[j+2]]].Position.x - Bead[Molecule[i].Bead[bead[j+1]]].Position.x;
            v.y = Bead[Molecule[i].Bead[bead[j+2]]].Position.y - Bead[Molecule[i].Bead[bead[j+1]]].Position.y;
            v.z = Bead[Molecule[i].Bead[bead[j+2]]].Position.z - Bead[Molecule[i].Bead[bead[j+1]]].Position.z;
            // calculate angle between the two vectors
            double size[2];
            size[0] = sqrt(SQR(u.x) + SQR(u.y) + SQR(u.z));
            size[1] = sqrt(SQR(v.x) + SQR(v.y) + SQR(v.z));
            double scalar = u.x * v.x + u.y * v.y + u.z * v.z;
//          printf("n[0]=( %lf %lf %lf )\n", n[0].x, n[0].y, n[0].z);
//          printf("n[1]=( %lf %lf %lf )\n", n[1].x, n[1].y, n[1].z);
            angle[i][j/3] = acos(scalar / (size[0] * size[1])); // in rad
            angle[i][j/3] *= 180 / PI; // in degrees
          }
        }
      }
    } //}}}

    // write angles to output file //{{{
    if ((out = fopen(output, "a")) == NULL) {
      ErrorFileOpen(output, 'a');
      exit(1);
    }

    fprintf(out, "%6d", count);
    for (int i = 0; i < Counts.Molecules; i++) {
      int mol_type_i = Molecule[i].Type;
      if (MoleculeType[mol_type_i].Use) {
        for (int j = 0; j < number_of_beads; j += 3){
          if (bead[j+0] < MoleculeType[mol_type_i].nBeads &&
              bead[j+1] < MoleculeType[mol_type_i].nBeads &&
              bead[j+2] < MoleculeType[mol_type_i].nBeads) {
            fprintf(out, " %10.6f", angle[i][j/3]);
          }
        }
      }
    }
    putc('\n', out);

    // write stuff
    fclose(out); //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      fprintf(stdout, "\n%s", stuff);
    } //}}}
  }
  fclose(vcf);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\rLast Step: %6d\n", count);
    }
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff); //}}}

  return 0;
}
