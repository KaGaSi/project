#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "../Errors.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <output.vsf> <options>\n\n", cmd);
  fprintf(stderr, "   <output.vsf>    output structure file (*.vsf)\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -i <name>  use input .vsf file different from traject.vsf\n");
//fprintf(stderr, "      -b <name>  file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v         verbose output\n");
  fprintf(stderr, "      -h         print this help and exit\n");
} //}}}

// WriteVsf() //{{{
/*
 * Function creating `.vsf` structure file for use in conjunction with
 * `.vcf` coordinate file for better visualisation via VMD program.
 */
void WriteVsf(char *input_vsf, Counts Counts, BeadType *BeadType, Bead *Bead,
              MoleculeType *MoleculeType, Molecule *Molecule) {

  // opten structure file //{{{
  FILE *fw;
  if ((fw = fopen(input_vsf, "w")) == NULL) {
    ErrorFileOpen(input_vsf, 'w');
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
  fprintf(fw, "mass %9.5f ", BeadType[type_def].Mass);
  fprintf(fw, "charge %9.5f\n", BeadType[type_def].Charge); //}}}

  // print beads //{{{
  for (int i = 0; i < Counts.BeadsInVsf; i++) {

    // don't print beads with type 'type_def'
    if (Bead[i].Type != type_def) {
      fprintf(fw, "atom %7d ", i);
      fprintf(fw, "name %8s ", BeadType[Bead[i].Type].Name);
      fprintf(fw, "mass %9.5f ", BeadType[Bead[i].Type].Mass);
      fprintf(fw, "charge %9.5f", BeadType[Bead[i].Type].Charge);
      if (Bead[i].Molecule != -1) {
        fprintf(fw, " resname %10s ", MoleculeType[Molecule[Bead[i].Molecule].Type].Name);
        fprintf(fw, "resid %5d", Bead[i].Molecule+1);
      }
      putc('\n', fw);
    // print highest bead id even if it's default type
    } else if (i == (Counts.BeadsInVsf-1)) {
      fprintf(fw, "atom %7d ", i);
      fprintf(fw, "name %8s ", BeadType[Bead[i].Type].Name);
      fprintf(fw, "mass %lf ", BeadType[Bead[i].Type].Mass);
      fprintf(fw, "charge %lf", BeadType[Bead[i].Type].Charge);
      putc('\n', fw);
    }
  } //}}}

  // print bonds //{{{
  putc('\n', fw);
  for (int i = 0; i < Counts.Molecules; i++) {
    fprintf(fw, "# resid %d\n", i+1); // in VMD resid start with 1
    int mol_type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[mol_type].nBonds; j++) {
      fprintf(fw, "bond %6d: %6d\n", Molecule[i].Bead[MoleculeType[mol_type].Bond[j][0]],
                                     Molecule[i].Bead[MoleculeType[mol_type].Bond[j][1]]);
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
TransformVsf reads information from FIELD and traject.vsf files and creates \
.vsf structure file used for visualisation of trajectory (.vcf files) via VMD \
visualisation tool.\n\n");

/*      fprintf(stdout, "\
The utility uses traject.vsf (or other input structure file) and FIELD (along \
with optional bond file) files to determine all information about the \
system.\n\n");
*/

      fprintf(stdout, "Usage:\n");
      fprintf(stdout, "   %s <output.vsf>\n\n", argv[0]);
      fprintf(stdout, "   <output.vsf>    output structure file (*.vsf)\n");
      fprintf(stdout, "   <options>\n");
      fprintf(stdout, "      -i   use input .vsf file different from traject.vsf\n");
//    fprintf(stdout, "      -b   file containing bond alternatives to FIELD\n");
      fprintf(stdout, "      -v   verbose output\n");
      fprintf(stdout, "      -h   print this help and exit\n");
      exit(0);
    }
  }

  int req_args = 1; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count][0] != '-'; i++) {
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
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0) {

      ErrorOption(argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    fprintf(stdout, " %s", argv[i]);
  putchar('\n'); //}}}

  // options before reading system data //{{{
  // use .vsf file other than traject.vsf? //{{{
  char *input_vsf = calloc(1024,sizeof(char *));
  if (FileOption(argc, argv, "-i", &input_vsf)) {
    exit(1);
  }
  if (input_vsf[0] == '\0') {
    strcpy(input_vsf, "traject.vsf");
  }

  // test if structure file ends with '.vsf'
  int ext = 2;
  char **extension;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(input_vsf, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // use bonds file? //{{{
  char *bonds_file = calloc(1024,sizeof(char *));
  if (FileOption(argc, argv, "-b", &bonds_file)) {
    exit(0);
  } //}}}

  // output verbosity //{{{
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  // }}}
  //}}}

  count = 0; // count arguments

  // <output.vsf> - output structure file (must end with .vsf) //{{{
  char output[1024];
  strcpy(output, argv[++count]);

  // test if <output.vsf> filename ends with '.vsf' or '.vtf' (required by VMD)
  ext = 2;
  extension = malloc(ext*sizeof(char *));
  for (int i = 0; i < ext; i++) {
    extension[i] = malloc(5*sizeof(char));
  }
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (!ErrorExtension(output, ext, extension)) {
    ErrorHelp(argv[0]);
    exit(1);
  }
  for (int i = 0; i < ext; i++) {
    free(extension[i]);
  }
  free(extension); //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information{{{
  char vcf[1];
  vcf[0] = '\0';
  ReadStructure(input_vsf, vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf); //}}}
//for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//  printf("%d\n", MoleculeType[i].nBTypes);
//}

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
