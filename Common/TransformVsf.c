#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../AnalysisTools.h"
#include "../Options.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <output.vsf> <options>\n\n", cmd);
  fprintf(stderr, "   <output.vsf>    output structure file (*.vsf)\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -i <name>  use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>  file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v         verbose output\n");
  fprintf(stderr, "      -h         print this help and exit\n");
} //}}}

// WriteVsf() //{{{
/*
 * Function creating `.vsf` structure file for use in conjunction with
 * `.vcf` coordinate file for better visualisation via VMD program.
 */
void WriteVsf(char *vsf_file, Counts Counts, BeadType *BeadType, Bead *Bead,
              MoleculeType *MoleculeType, Molecule *Molecule) {

  // opten structure file //{{{
  FILE *fw;
  if ((fw = fopen(vsf_file, "w")) == NULL) {
    fprintf(stderr, "Cannot open file input.vsf for writing!\n");
    exit(1);
  } //}}}

  // find most common type of bead and make it default //{{{
  int type_def = 0, count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Number >= count) {
      count = BeadType[i].Number;
      type_def = i;
    }
  } //}}}

  // print default bead type //{{{
  fprintf(fw, "atom default name %8s ", BeadType[type_def].Name);
  fprintf(fw, "mass %4.2f ", BeadType[type_def].Mass);
  fprintf(fw, "charge %5.2f\n", BeadType[type_def].Charge); //}}}

  // print unbonded beads //{{{
  for (int i = 0; i < Counts.Unbonded; i++) {

    // don't print beads with type 'type_def'
    if (Bead[i].Type != type_def) {
      fprintf(fw, "atom %7d ", i);
      fprintf(fw, "name %8s ", BeadType[Bead[i].Type].Name);
      fprintf(fw, "mass %4.2f ", BeadType[Bead[i].Type].Mass);
      fprintf(fw, "charge %5.2f\n", BeadType[Bead[i].Type].Charge);
    }
  } //}}}

  // print bonded beads //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;

    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int bead = Molecule[i].Bead[j];
      fprintf(fw, "atom %7d ", bead);
      fprintf(fw, "name %8s ", BeadType[Bead[bead].Type].Name);
      fprintf(fw, "mass %4.2f ", BeadType[Bead[bead].Type].Mass);
      fprintf(fw, "charge %5.2f ", BeadType[Bead[bead].Type].Charge);
      fprintf(fw, "segid %10s ", MoleculeType[type].Name);
      fprintf(fw, "resid %5d\n", i+1);
    }
  } //}}}

  // print bonds //{{{
  putc('\n', fw);
  count = 0;
  // go through all molecule types
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // go through all molecules of type 'i'
    for (int j = 0; j < MoleculeType[i].Number; j++) {
      // go through all bonds of 'j'-th molecule of type 'i'
      fprintf(fw, "# resid %d\n", count+1); // in VMD resid start with 1
      for (int k = 0; k < MoleculeType[i].nBonds; k++) {
        fprintf(fw, "bond %6d: %6d\n", Molecule[count].Bead[MoleculeType[i].Bond[k][0]],
                                       Molecule[count].Bead[MoleculeType[i].Bond[k][1]]);
      }

      count++;
    }
  } //}}}

  // close structure file
  fclose(fw);
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      fprintf(stdout, "\
TransformVsf reads information from FIELD and dl_meso.vsf files and creates \
.vsf structure file used for visualisation of trajectory (.vcf files) via VMD \
visualisation tool.\n\n");

/*      fprintf(stdout, "\
The utility uses dl_meso.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");
*/

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <output.vsf>\n\n", argv[0]);
      fprintf(stdout, "   <output.vsf>    output structure file (*.vsf)\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -i   use input .vsf file different from dl_meso.vsf\n");
      fprintf(stdout, "      -b   file containing bond alternatives to FIELD\n");
      fprintf(stdout, "      -v   verbose output\n");
      fprintf(stdout, "      -h   print this help and exit\n");
      exit(0);
    }
  }

  int options = 1; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count][0] != '-'; i++) {
    count++;
  }

  if (count < options) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    fprintf(stdout, " %s", argv[i]);
  putchar('\n'); //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
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
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  // }}}
  //}}}

  count = 0; // count arguments

  // <output.vsf> - file name of output structure file (must end with .vsf) //{{{
  char output[32];
  strcpy(output, argv[++count]);

  // test if <output.vsf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(output, '.');
  if (!dot || strcmp(dot, ".vsf")) {
    fprintf(stderr, "<output.vsf> %s does not have .vsf ending!\n", output);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information{{{
  char vcf[1];
  vcf[0] = '\0';
  ReadStructure(vsf_file, vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file); //}}}

  // print information - verbose option //{{{
  if (verbose) {
    char null[1] = {'\0'};
    VerboseOutput(false, null, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);
  }

  free(bonds_file); //}}}

  // create & fill output vsf file
  WriteVsf(output, Counts, BeadType, Bead, MoleculeType, Molecule);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead); //}}}

  return 0;
}
