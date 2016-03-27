#ifndef _AGGREGATES_H_
#define _AGGREGATES_H_

// CalculateAggregates() //{{{
/*
 * Function determine which molecules are in what aggregate.
 *
 * \param [out] Aggregate     information about aggregates
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [in]  sqdist        squared distance for closeness check
 * \param [in]  contacts      minimal number of contacts aggregate check
 * \param [in]  BoxLength     dimensions of simulation box
 * \param [in]  BeadType      informationn about bead types
 * \param [in]  Bead          informationn about individual beads
 * \param [in]  MoleculeType  information about molecule types
 * \param [in]  Molecule      information about individual molecules
 */
void CalculateAggregates(Aggregate **Aggregate, Counts *Counts, int sqdist, int contacts,
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule); //}}}
#endif
