/**
 * \file
 * \brief Functions reading files.
 */

#ifndef _READ_H_
#define _READ_H_

#include "AnalysisTools.h"

// GetPBC() //{{{
/*
 * \brief Function to get box dimensions.
 *
 * \param [in] vcf         opened coordinate file
 * \param [in] input_coor  name of the coordinate file
 * \return vector with box dimensions
 */
VECTOR GetPBC(FILE *vcf, char *input_coor); //}}}

// ReadAggCommand() //{{{
/**
 * \brief Function reading Aggregate command from agg file.
 *
 * \param [in]  BeadType      information about bead types
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [in]  input_coor    coordinate file
 * \param [in]  input_agg     aggregate file
 * \param [out] distance      <distance> parameter from Aggregate command
 * \param [out] contacts      <contacts> parameter from Aggregate command
 */
void ReadAggCommand(BEADTYPE *BeadType, COUNTS Counts,
                    char *input_coor, char *input_agg,
                    double *distance, int *contacts); //}}}

// ReadStructure() //{{{
/**
 * \brief Function reading information from dl_meso FIELD and vsf
 * structure files.
 *
 * \param [in]  vsf_file      .vsf structure file
 * \param [in]  vcf_file      .vcf coordinate file
 * \param [out] Counts        numbers of beads, molecules, etc.
 * \param [out] BeadType      information about bead types
 * \param [out] Bead          informationn about individual beads
 * \param [out] Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * \return 'true' or 'false' for .vcf file with indexed or ordered
 * timesteps, respectively
 * */
bool ReadStructure(char *vsf_file, char *vcf_file, COUNTS *Counts,
                   BEADTYPE **BeadType, BEAD **Bead, int **Index,
                   MOLECULETYPE **MoleculeType, MOLECULE **Molecule); //}}}

int ReadTimestepPreamble(bool indexed, char *input_coor, FILE *vcf_file, char **stuff);

// ReadCoordinates() //{{{
/**
 * \brief Function reading ordered coordinates from .vcf coordinate file.
 *
 * \param [in]  indexed    is the vcf indexed?
 * \param [in]  input_coor name of input coordinate file
 * \param [in]  vcf_file   pointer to the open coordinate file
 * \param [in]  Counts     numbers of beads, molecules, etc.
 * \param [in]  Index      bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] Bead       coordinates of individual beads
 * \param [out] stuff      first line of a timestep
 */
void ReadCoordinates(bool indexed, char *input_coor, FILE *vcf_file, COUNTS Counts, int *Index, BEAD **Bead, char **stuff); //}}}

// SkipCoor() //{{{
/**
 * \brief Function to skip one timestep in coordinates file.
 *
 * \param [in]  vcf_file   file with vcf coordinates
 * \param [in]  Counts     number of beads in vcf file
 * \param [out] stuff      first line of a timestep
 * \return 1 if premature end of file or 0 for no error
 */
bool SkipCoor(FILE *vcf_file, COUNTS Counts, char **stuff); //}}}

// ReadAggregates() //{{{
/**
 * \brief Function reading information about aggregates from `.agg` file
 *
 * \param [in]  fr            pointer to open aggregate file
 * \param [in]  agg_file      name of aggregate file
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [in]  BeadType      information about bead types
 * \param [out] Bead          information about individual beads
 * \param [out] Aggregate     information about aggregates
 * \param [in]  MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * \return 1 if 'Last Step' detected or 0 for no error
 */
bool ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts, AGGREGATE **Aggregate,
                    BEADTYPE *BeadType, BEAD **Bead,
                    MOLECULETYPE *MoleculeType, MOLECULE **Molecule, int *Index); //}}}

// ReadField() //{{{
/**
 * \brief Function reading structure information from FIELD-like file
 *
 * \param [in]  field         input FIELD-like file
 * \param [out] BoxLength     simulation box size
 * \param [out] Counts        numbers of beads, molecules, etc.
 * \param [out] BeadType      information about bead types
 * \param [out] Bead          informationn about individual beads
 * \param [out] Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * \param [out] bond_types    information abount bond types
 * \param [out] angle_types   information abount angle types
 * */
void ReadField(char *field, VECTOR *BoxLength, COUNTS *Counts,
               BEADTYPE **BeadType, BEAD **Bead, int **Index,
               MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
               PARAMS **bond_types, PARAMS **angle_types); //}}}

// ReadLmpData() //{{{
/**
 * \brief Function reading all information from lammps data file
 *
 * \param [in]  data_field    input data file file
 * \param [out] bonds         number of bonds
 * \param [out] bond_type     information about bond types
 * \param [out] angles        number of angles
 * \param [out] angle_type    information about angle types
 * \param [out] BoxLength     simulation box size
 * \param [out] box_lo        minimum box coordinates
 * \param [out] Counts        numbers of beads, molecules, etc.
 * \param [out] BeadType      information about bead types
 * \param [out] Bead          informationn about individual beads
 * \param [out] Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * */
void ReadLmpData(char *data_file, int *bonds, PARAMS **bond_type,
                 int *angles, PARAMS **angle_type,
                 VECTOR *BoxLength, VECTOR *box_lo, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule); //}}}
#endif
