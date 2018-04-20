#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "../AnalysisTools.h"
#include "../Options.h"

int main(int argc, char *argv[]) {

  char *vsf_file = "dl_meso.vsf";
  char *bonds_file = "\0";

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  char vcf[32];
  strcpy(vcf, argv[1]);
  printf("vcf file: %s\n", vcf);
  ReadStructure(vsf_file, vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  printf("Bonded=%d Unbonded=%d\n", Counts.Bonded, Counts.Unbonded);
  printf("Beads in vsf=%d Molecules in vsf=%d\n", Counts.BeadsInVsf, Counts.MoleculesInVsf);
  printf("Beads=%d Molecules=%d\n", Counts.Beads, Counts.Molecules);
  printf("Bead types=%d Molecule types=%d\n", Counts.TypesOfBeads, Counts.TypesOfMolecules);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    printf("%s %d %d\n", BeadType[i].Name, i, BeadType[i].Number);
  }
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    printf("%s %d %d %d\n", MoleculeType[i].Name, i, MoleculeType[i].Number, MoleculeType[i].nBeads);
  }

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead); //}}}

  return 0;
}
