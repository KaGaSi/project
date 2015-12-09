/**
 * \file
 * \brief Auxiliary functions that are often used.
 */

#ifndef _AUX_H_
#define _AUX_H_

/** \brief Function to identify type of bead from its name
 *
 * \param [in]  name      bead name
 * \param [in]  Counts    numbers of beads, residues, etc.
 * \param [in]  BeadType  informationn about bead types
 * \return bead type id corresponding to index in BeadType struct
 */
int FindType(char *name, Counts Counts, BeadType *BeadType);
#endif
