/**
 * \file
 * \brief Functions common to all analysis utilities.
 */

#ifndef _ANALYSISTOOLS_H_
#define _ANALYSISTOOLS_H_

#include "Structs.h"

// CommonHelp() //{{{
/**
 * \brief Function printing help for common options.
 *
 * \param [in] error   `true` or `false` whether to use stderr or stdout
 */
void CommonHelp(bool error); //}}}

// VerboseOutput() //{{{
/**
 * \brief Function printing basic information about system if `-v` or `-V`
 * option is provided
 *
 * \param [in] input_vcf     .vcf coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          informationn about individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      information about individual molecules
 */
void VerboseOutput(char *input_vcf, Counts Counts,
                   BeadType *BeadType, Bead *Bead,
                   MoleculeType *MoleculeType, Molecule *Molecule); //}}}

// PrintCounts()  //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts   numbers of beads, molecules, etc.
 */
void PrintCounts(Counts Counts);
//}}}

// PrintBeadType() //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts     numbers of beads, molecules, etc.
 * \param [in] BeadType   information about bead types
 */
void PrintBeadType(Counts Counts, BeadType *BeadType); //}}}

// PrintMoleculeTypeType()  //{{{
/**
 * \brief Function printing MoleculeType structure.
 *
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] MoleculeType  information about molecule types
 */
void PrintMoleculeType(Counts Counts, BeadType *BeadType, MoleculeType *MoleculeType); //}}}

// PrintBead() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead(Counts Counts, int *Index, BeadType *BeadType, Bead *Bead); //}}}

// PrintMolecule() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule(Counts Counts, int *Index, MoleculeType *MoleculeType, Molecule *Molecule, Bead *Bead, BeadType *BeadType); //}}}

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

// FindBeadType() //{{{
/** \brief Function to identify type of bead from its name
 *
 * \param [in]  name      bead name
 * \param [in]  Counts    numbers of beads, residues, etc.
 * \param [in]  BeadType  information about bead types
 * \return bead type id corresponding to index in BeadType struct (or -1 if non-existent bead name)
 */
int FindBeadType(char *name, Counts Counts, BeadType *BeadType); //}}}

// FindMoleculeType() //{{{
/** \brief Function to identify type of molecule from its name
 *
 * \param [in]  name          molecule name
 * \param [in]  Counts        numbers of beads, residues, etc.
 * \param [in]  MoleculeType  information about bead types
 * \return molecule type      id corresponding to index in BeadType struct (or -1 for non-existent molecule)
 */
int FindMoleculeType(char *name, Counts Counts, MoleculeType *MoleculeType); //}}}

// Distance between two beads //{{{
/**
 * \brief Function to calculate distance vector between two beads.
 *
 * \param [in] id1         first coordinate vector
 * \param [in] id2         second coordinate vector
 * \param [in] BoxLength   dimensions of simulation box
 * \return distance vector between the two provided beads (without pbc)
 */
Vector Distance(Vector id1, Vector id2, Vector BoxLength); //}}}

// RemovePBCMolecules() //{{{
/**
 * \brief Function to join all molecules.
 *
 * \param [in]  Counts         numbers of beads, molecules, etc.
 * \param [in]  BoxLength      dimension of the simulation box
 * \param [in]  BeadType       information about bead types
 * \param [out] Bead           information about individual beads (coordinates)
 * \param [in]  MoleculeType   information about molecule types
 * \param [in]  Molecule       information about individual molecules
 */
void RemovePBCMolecules(Counts Counts, Vector BoxLength,
                        BeadType *BeadType, Bead **Bead,
                        MoleculeType *MoleculeType, Molecule *Molecule); //}}}

// RemovePBCAggregates() //{{{
/**
 * \brief Funcion to join all aggregates.
 *
 * \param [in]  distance       distance for closeness check (taken from agg file)
 * \param [in]  Aggregate      information about aggregates
 * \param [in]  Counts         number of beads, molecu.es, etc.
 * \param [in]  BoxLength      dimensions of the simulation box
 * \param [in]  BeadType       information about bead types
 * \param [out] Bead           information about individual beads (coordinates)
 * \param [in]  MoleculeType   information about molecule types
 * \param [in]  Molecule       information about individual molecules
 */
void RemovePBCAggregates(double distance, Aggregate *Aggregate, Counts Counts,
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule); //}}}

// RestorePBC() //{{{
/**
 * \brief Function to restore pbc.
 *
 * \param [in]  Counts         numbers of beads, molecules, etc.
 * \param [in]  BoxLength      dimension of the simulation box
 * \param [out] Bead           information about individual beads (coordinates)
 */
void RestorePBC(Counts Counts, Vector BoxLength, Bead **Bead); //}}}

// CentreOfMass() //{{{
/**
 * \brief Function to calculate centre of mass for a collection of beads.
 *
 * \param [in] n              number of beads
 * \param [in] list           list of bead ids (corresponding to indices in Bead struct)
 * \param [in] Bead           information about individual beads (coordinates)
 * \param [in] BeadType       information about beadtypes (masses)
 * \return coordinates of centre of mass of a given aggregate
 */
Vector CentreOfMass(int n, int *list, Bead *Bead, BeadType *BeadType); //}}}

// Gyration() //{{{
/*
 * \brief Function calculating principal moments of the gyration tensor.
 *
 * \param [in] n             number of beads
 * \param [in] list          array of bead ids
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BoxLength     dimensions of simulation box
 * \param [in] BeadType      informationn about bead types
 * \param [in] Bead          informationn about individual beads
 * \return vector with principal moments of gyration tensor (sorted as x<y<z)
 */
Vector Gyration(int n, int *list, Counts Counts, Vector BoxLength,
                BeadType *BeadType, Bead **Bead); //}}}

// Min3() //{{{
/**
 * \brief Function returning the lowest number from three floats.
 *
 * \param [in] x   first double precision number
 * \param [in] y   second double precision number
 * \param [in] z   third double precision number
 * \return lowest of the supplied numbers
 */
double Min3(double x, double y, double z); //}}}

// Max3() //{{{
/**
 * \brief Function returning the highest number from three floats.
 *
 * \param [in] x   first double precision number
 * \param [in] y   second double precision number
 * \param [in] z   third double precision number
 * \return highest of the supplied numbers
 */
double Max3(double x, double y, double z); //}}}

// Sort3() //{{{
/**
 * \brief Function returning sorted numbers x < y < z.
 *
 * \param [in] in   first double precision number
 * \return sorted vector
 */
Vector Sort3(Vector in); //}}}

// Swap() //{{{
/**
 * \brief Function to swap two integers.
 *
 * \param [in] a   first integer to swap
 * \param [in] b   second integer to swap
 */
void Swap(int *a, int *b);
// }}}

// SwapDouble() //{{{
/**
 * \brief Function to swap two doubles.
 *
 * \param [in] a   first integer to swap
 * \param [in] b   second integer to swap
 */
void SwapDouble(double *a, double *b);
// }}}

// SwapBool() //{{{
/**
 * \brief Function to swap two booleans.
 *
 * \param [in] a   first integer to swap
 * \param [in] b   second integer to swap
 */
void SwapBool(bool *a, bool *b);
// }}}

// SortArray() //{{{
/**
 * \brief Function to sort an integer array.
 *
 * \param [out] array   integer array to sort
 * \param [in]  length  array length
 * \param [in]  mode    0 for ascending order, 1 for descending order
 */
void SortArray(int **array, int length, int mode); //}}}

// SortAggStruct() //{{{
/**
 * \brief Function to sort Aggregate struct.
 *
 * \param [out] Aggregate  Aggregate struct to sort
 * \param [in]  Counts     numbers of beads, molecules, etc.
 */
void SortAggStruct(Aggregate **Aggregate, Counts Counts); //}}}

// SortBonds() //{{{
/**
 * \brief Function to sort 2D array of bonds.
 *
 * \param [out] bond    2D array of bonds
 * \param [in]  length  number of bonds
 */
void SortBonds(int **bond, int length); //}}}

// ZeroCounts() //{{{
/**
 * \brief Zeroize Counts structure.
 *
 * \param[in] Counts   Counts structure to zeroize
 */
void ZeroCounts(Counts *Counts);

// FreeBead() //{{{
/**
 * \brief Free memory allocated for Bead struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Bead        information about individual beads
 */
void FreeBead(Counts Counts, Bead **Bead); //}}}

// FreeMolecule() //{{{
/**
 * \brief Free memory allocated for Molecule struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Molecule    information about individual molecules
 */
void FreeMolecule(Counts Counts, Molecule **Molecule); //}}}

// FreeMoleculeType() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  Counts         number of beads, molecu.es, etc.
 * \param [out] MoleculeType   information about individual molecules
 */
void FreeMoleculeType(Counts Counts, MoleculeType **MoleculeType); //}}}

// FreeAggregate() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Aggregate   information about individual molecules
 */
void FreeAggregate(Counts Counts, Aggregate **Aggregate); //}}}
#endif
