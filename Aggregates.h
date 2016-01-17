#ifndef _AGGREGATES_H_
#define _AGGREGATES_H_

// CalculateAggregates() //{{{
/* Function to determine which molecules belong to which aggregate
 * according to minimum number of close pairs criterium.
 *
 * \param [out] Aggregate     information about aggregate
 * \param [out] Counts        number of aggregates (and molecules, beads, etc.)
 * \param [in]  sqdist        squared distance for closeness check
 * \param [in]  contacs       minimum number of contact pairs for aggregate check
 * \param [in]  BoxLength     vector of simulation box length
 * \param [in]  BeadType      information about bead types
 * \param [in]  Bead          informationn about individual beads
 * \param [in]  MoleculeType  information about molecule types
 * \param [in]  Molecule      information about individual molecules
 */
void CalculateAggregates(Aggregate **Aggregate, Counts *Counts, int sqdist, int contacts,
                         Vector BoxLength, BeadType *BeadType, Bead *Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule); //}}}

#endif
