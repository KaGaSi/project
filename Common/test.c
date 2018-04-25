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

//// prints for debugging //{{{
//// Counts struct
//printf("Counts.{TypesOfMolecules = %d,", Counts.TypesOfMolecules);
//printf(" Molecules = %d,", Counts.Molecules);
//printf(" TypesOfBeads = %d,", Counts.TypesOfBeads);
//printf(" BeadsInVsf = %d,", Counts.BeadsInVsf);
//printf(" Beads = %d}\n", Counts.Beads);
//// MoleculeType struct
//for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//  printf("MoleculeType[%d].{Name = %s, ", i, MoleculeType[i].Name);
//  printf("Number = %d, ", MoleculeType[i].Number);
//  printf("nBeads = %d, ", MoleculeType[i].nBeads);
//  printf("nBTypes = %d, BTypes = {", MoleculeType[i].nBTypes);
//  for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
//    if (j == 0 && MoleculeType[i].nBTypes == 1) {
//      printf("%d} ", MoleculeType[i].BType[j]);
//    } else if (j == 0) {
//      printf("%d,", MoleculeType[i].BType[j]);
//    } else if (j < (MoleculeType[i].nBTypes-1)) {
//      printf(" %d,", MoleculeType[i].BType[j]);
//    } else {
//      printf(" %d}, ", MoleculeType[i].BType[j]);
//    }
//  }
//  printf("nBonds = %d, Bond = {", MoleculeType[i].nBonds);
//  for (int j = 0; j < MoleculeType[i].nBonds; j++) {
//    if (j == 0 && MoleculeType[i].nBonds == 1) {
//      printf("%d-%d}", MoleculeType[i].Bond[j][0], MoleculeType[i].Bond[j][1]);
//    } else if (j == 0) {
//      printf("%d-%d,", MoleculeType[i].Bond[j][0], MoleculeType[i].Bond[j][1]);
//    } else if (j < (MoleculeType[i].nBonds-1)) {
//      printf(" %d-%d,", MoleculeType[i].Bond[j][0], MoleculeType[i].Bond[j][1]);
//    } else {
//      printf(" %d-%d}", MoleculeType[i].Bond[j][0], MoleculeType[i].Bond[j][1]);
//    }
//  }
//  printf("}\n");
//}
//// BeadType struct
//for (int i = 0; i < Counts.TypesOfBeads; i++) {
//  printf("BeadType[%d].{Name = %s, Number = %d, Use = %s}\n", i,
//                                                             BeadType[i].Name,
//                                                             BeadType[i].Number,
//                                                             BeadType[i].Use? "Yes":"No");
//}
//// Molecule struct
//for (int i = 0; i < Counts.Molecules; i++) {
//  int type = Molecule[i].Type;
//  printf("Molecule[%d].{Type = %d", i, Molecule[i].Type);
//  printf(" (%s),", MoleculeType[type].Name);
//  printf(" Bead (%d) = {", MoleculeType[type].nBeads);
//  for (int j = 0; j < MoleculeType[type].nBeads; j++) {
//    if (j == 0 && j == MoleculeType[type].nBeads) {
//      printf("%d}", Molecule[i].Bead[j]);
//    } else if (j == 0) {
//      printf("%d,", Molecule[i].Bead[j]);
//    } else if (j < (MoleculeType[type].nBeads-1)) {
//      printf(" %d,", Molecule[i].Bead[j]);
//    } else {
//      printf(" %d}", Molecule[i].Bead[j]);
//    }
//  }
//  printf("}\n");
//}
//// Bead struct
//for (int i = 0; i < Counts.Beads; i++) {
//  printf("Bead[%d].{Type = %d", i, Bead[i].Type);
//  printf(" (%s),", BeadType[Bead[i].Type].Name);
//  printf(" Molecule = %d,", Bead[i].Molecule);
//  printf(" Index = %d}\n", Bead[i].Index);
//} //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead); //}}}

  return 0;
}
