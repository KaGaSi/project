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
Vector GetPBC(FILE *vcf, char *input_coor); //}}}

// ReadAggCommand() //{{{
/**
 * \brief Function reading Aggregate command from agg file.
 *
 * \param [in]  BeadType      information about bead types
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [in]  input_coor    coordinate file
 * \param [in]  input_agg     aggregate file
 * \param [in]  agg           opened aggregate file
 * \param [out] distance      <distance> parameter from Aggregate command
 * \param [out] contacts      <contacts> parameter from Aggregate command
 */
void ReadAggCommand(BeadType *BeadType, Counts Counts,
                    char *input_coor, char *input_agg,
                    FILE *agg, double *distance, int *contacts); //}}}

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
bool ReadStructure(char *vsf_file, char *vcf_file, Counts *Counts,
                   BeadType **BeadType, Bead **Bead, int **Index,
                   MoleculeType **MoleculeType, Molecule **Molecule); //}}}

// ReadCoordinates() //{{{
/**
 * \brief Function reading ordered coordinates from .vcf coordinate file.
 *
 * \param [in]  indexed    is the vcf indexed?
 * \param [in]  vcf_file   name of input .vcf coordinate file
 * \param [in]  Counts     numbers of beads, molecules, etc.
 * \param [in]  Index         bead indices between program and vsf (i.e., opposite of Bead[].Index)
 * \param [out] Bead       coordinates of individual beads
 * \param [out] stuff      first line of a timestep
 * \return 0 for no errors or index number of bead (starting from 1) for which coordinates cannot be read
 */
int ReadCoordinates(bool indexed, FILE *vcf_file, Counts Counts, int *Index, Bead **Bead, char **stuff); //}}}

// SkipCoor() //{{{
/**
 * \brief Function to skip one timestep in coordinates file.
 *
 * \param [in]  vcf_file   file with vcf coordinates
 * \param [in]  Counts     number of beads in vcf file
 * \param [out] stuff      first line of a timestep
 * \return 1 if premature end of file or 0 for no error
 */
bool SkipCoor(FILE *vcf_file, Counts Counts, char **stuff); //}}}

// ReadAggregates() //{{{
/**
 * \brief Function reading information about aggregates from `.agg` file
 *
 * \param [in]  agg_file      name of input aggregate file
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [in]  BeadType      information about bead types
 * \param [out] Bead          information about individual beads
 * \param [out] Aggregate     information about aggregates
 * \param [in]  MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * \return 1 if 'Last Step' detected or 0 for no error
 */
bool ReadAggregates(FILE *agg_file, Counts *Counts, Aggregate **Aggregate,
                    BeadType *BeadType, Bead **Bead,
                    MoleculeType *MoleculeType, Molecule **Molecule, int *Index); //}}}
#endif
