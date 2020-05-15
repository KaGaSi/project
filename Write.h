/**
 * \file
 * \brief Functions writing to files.
 */

#ifndef _WRITE_H_
#define _WRITE_H_

#include "AnalysisTools.h"

// WriteCoorIndexed //{{{
/**
 * \brief Function writing indexed coordinates to a .vcf file.
 *
 * \param [in] vcf_file      name of output .vcf coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          coordinates of individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      coordinates of individual molecules
 * \param [in] stuff         array of chars containing comment line to place at the beginning
 */
void WriteCoorIndexed(FILE *vcf_file, Counts Counts, BeadType *BeadType, Bead *Bead, MoleculeType *MoleculeType, Molecule *Molecule, char *stuff); //}}}

// WriteCoorXYZ() //{{{
/**
 * \brief Function for writing xyz coordinates
 *
 * \param [in] xyz_file      output .xyz coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          coordinates of individual beads
 */
void WriteCoorXYZ(FILE *xyz_file, Counts Counts,
                  BeadType *BeadType, Bead *Bead); //}}}

// WriteVsf() //{{{
/**
 * \brief Function writing vsf file
 *
 * \param [in] vsf_file      name of output .vsf structure file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          coordinates of individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      coordinates of individual molecules
 */
void WriteVsf(char *input_vsf, Counts Counts, BeadType *BeadType, Bead *Bead,
              MoleculeType *MoleculeType, Molecule *Molecule); //}}}
#endif
