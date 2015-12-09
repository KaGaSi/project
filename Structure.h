/**
 * \file
 * \brief functions to manipulate DL_MESO FIELD and vsf structure files.
 */

#ifndef _TRANSFORMVSF_H_
#define _TRANSFORMVSF_H_

/** \brief Function reading information from dl_meso FIELD file.
 *
 * \param [out] Counts       numbers of beads, molecules, etc.
 * \param [out] BeadType     information about bead types
 * \param [out] MoleculeType    information about molecule types
 * */
void ReadFIELD(Counts *Counts, BeadType **BeadType, MoleculeType **MoleculeType);

/** \brief Function reading bead numbers from provided vsf file.
 *
 * \param [in]  file      name of vsf input file
 * \param [in]  Counts    numbers of beads, molecules, etc.
 * \param [in]  BeadType  informationn about bead types
 * \param [out] Bead      informationn about individual beads
 * */
void ReadVsf(char *file, Counts Counts, BeadType *BeadType, Bead **Bead);

// WriteVsf()
/**
 * \brief Function to create structure file input.vsf.
 *
 * \param [in] file        name of file to be created
 * \param [in] Counts      numbers of beads, molecules, etc.
 * \param [in] BeadType    informationn about bead types
 * \param [in] Bead        informationn about individual beads
 * \param [in] MoleculeType   information about molecule types
 * \param [in] Molecule     information about individual molecules
 */
void WriteVsf(char *file, Counts Counts, BeadType *BeadType, Bead *Bead,
              MoleculeType *MoleculeType, Molecule *Molecule);

#endif
