#include "Read.h"

// WriteCoorIndexed() //{{{
/**
 * Function writing coordinates to a `.vcf` file. According to the Write flag
 * in BeadType and MoleculeType structures only certain bead types will be
 * saved into the indexed timestep in .vcf file.
 */
void WriteCoorIndexed(FILE *vcf_file, Counts Counts,
                      BeadType *BeadType, Bead *Bead,
                      MoleculeType *MoleculeType, Molecule *Molecule,
                      char *stuff) {

  // print blank line
  putc('\n', vcf_file);
  // print comment at the beginning of a timestep if present in initial vcf file
  fprintf(vcf_file, "%s", stuff);
  // print 'indexed' on the next
  fprintf(vcf_file, "indexed\n");

  for (int i = 0; i < Counts.Beads; i++) {
    int btype = Bead[i].Type;
    if (BeadType[btype].Write) {
      if (Bead[i].Molecule != -1) { // bead in a molecule
        int mtype = Molecule[Bead[i].Molecule].Type;
        if (MoleculeType[mtype].Write) {
          fprintf(vcf_file, "%8d %7.3f %7.3f %7.3f\n", Bead[i].Index,
                                                       Bead[i].Position.x,
                                                       Bead[i].Position.y,
                                                       Bead[i].Position.z);
        }
      } else { // monomer bead
        fprintf(vcf_file, "%8d %7.3f %7.3f %7.3f\n", Bead[i].Index,
                                                     Bead[i].Position.x,
                                                     Bead[i].Position.y,
                                                     Bead[i].Position.z);
      }
    }
  }
} //}}}

// WriteCoorXYZ() //{{{
void WriteCoorXYZ(FILE *xyz_file, Counts Counts,
                  BeadType *BeadType, Bead *Bead) {

  // count beads to write
  int count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Write) {
      count += BeadType[i].Number;
    }
  }

  // print number of beads to file
  fprintf(xyz_file, "%d\n\n", count);

  // print coordinates
  for (int i = 0; i < Counts.Beads; i++) {
    int type = Bead[i].Type;
    if (BeadType[type].Write) {
      fprintf(xyz_file, "%8s %7.3f %7.3f %7.3f\n", BeadType[type].Name, Bead[i].Position.x, Bead[i].Position.y, Bead[i].Position.z);
    }
  }
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
  int type_def = -1, count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    bool use = true;
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      for (int k = 0; k < MoleculeType[j].nBTypes; k++) {
        if (MoleculeType[j].BType[k] == i) {
          use = false;
        }
      }
    }
    if (use && BeadType[i].Number >= count) {
      count = BeadType[i].Number;
      type_def = i;
    }
  } //}}}

// type_def = -1; // uncomment if no default is to be in vsf (for debugging)
  // print default bead type //{{{
  if (type_def != -1) {
    fprintf(fw, "atom default name %8s ", BeadType[type_def].Name);
    fprintf(fw, "mass %9.5f ", BeadType[type_def].Mass);
    fprintf(fw, "charge %9.5f\n", BeadType[type_def].Charge);
  } //}}}

  // print beads //{{{
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    // don't print beads with type 'type_def'
    if (Bead[i].Type != type_def || Bead[i].Molecule != -1) {
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
