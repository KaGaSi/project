#ifndef _TRANSFORMVSF_H_
#define _TRANSFORMVSF_H_

// WriteVsf() //{{{
/*
 * \brief Function to create structure file input.vsf.
 *
 * \param [in] vsf_file      name of file to be created
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      informationn about bead types
 * \param [in] Bead          informationn about individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      information about individual molecules
 */
void WriteVsf(char *vsf_file, Counts Counts, BeadType *BeadType, Bead *Bead,
              MoleculeType *MoleculeType, Molecule *Molecule); //}}}
#endif
