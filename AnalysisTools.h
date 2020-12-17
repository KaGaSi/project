/**
 * \file
 * \brief Functions common to all analysis utilities.
 */

#ifndef _ANALYSISTOOLS_H_
#define _ANALYSISTOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include "General.h"
#include "Errors.h"
#include "Structs.h"
#include "Options.h"
#include "Read.h"
#include "Write.h"

// InputCoor() //{{{
/**
 * \brief Function test input coordinate file is correct
 *
 * \param [out] vtf          is the coordinate file vtf or vcf?
 * \param [in]  file_coor    name of coordinate file
 * \param [in]  file_struct  name of structure file
 * \return false if file_coor has wrong extension, true otherwise
 */
bool InputCoor(bool *vtf, char *file_coor, char *file_struct); //}}}

// VerboseOutput() //{{{
/**
 * \brief Function printing basic information about system if `-v` or `-V`
 * option is provided
 *
 * \param [in] input_vcf     .vcf coordinate file
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BoxLength     dimension of the simulation box
 * \param [in] BeadType      information about bead types
 * \param [in] Bead          informationn about individual beads
 * \param [in] MoleculeType  information about molecule types
 * \param [in] Molecule      information about individual molecules
 */
void VerboseOutput(char *input_vcf, COUNTS Counts, VECTOR BoxLength,
                   BEADTYPE *BeadType, BEAD *Bead,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}

// PrintCounts()  //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts   numbers of beads, molecules, etc.
 */
void PrintCounts(COUNTS Counts);
//}}}

// PrintBeadType() //{{{
/**
 * \brief Function printing Counts structure.
 *
 * \param [in] Counts     numbers of beads, molecules, etc.
 * \param [in] BeadType   information about bead types
 */
void PrintBeadType(COUNTS Counts, BEADTYPE *BeadType); //}}}
void PrintBeadType2(int number, BEADTYPE *BeadType);

// PrintMoleculeTypeType()  //{{{
/**
 * \brief Function printing MoleculeType structure.
 *
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] BeadType      information about bead types
 * \param [in] MoleculeType  information about molecule types
 */
void PrintMoleculeType(COUNTS Counts, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType); //}}}
void PrintMoleculeType2(int number_of_types, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType);

// PrintBead() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead(COUNTS Counts, int *Index, BEADTYPE *BeadType, BEAD *Bead); //}}}
void PrintBead2(int number_of_beads, int *Index, BEADTYPE *BeadType, BEAD *Bead);

// PrintMolecule() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule(int number_of_molecules,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                   BEADTYPE *BeadType, BEAD *Bead); //}}}

// PrintAggregate() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintAggregate(COUNTS Counts, int *Index,
                    MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                    BEAD *Bead, BEADTYPE *BeadType, AGGREGATE *Aggregate); //}}}

// PrintBondTypes() //{{{
void PrintBondTypes(COUNTS Counts, PARAMS *bond_type);  //}}}

// PrintAngleTypes() //{{{
void PrintAngleTypes(COUNTS Counts, PARAMS *angle_type);  //}}}

// FindBeadType() //{{{
/** \brief Function to identify type of bead from its name
 *
 * \param [in]  name      bead name
 * \param [in]  Counts    numbers of beads, residues, etc.
 * \param [in]  BeadType  information about bead types
 * \return bead type id corresponding to index in BeadType struct (or -1 if non-existent bead name)
 */
int FindBeadType(char *name, COUNTS Counts, BEADTYPE *BeadType); //}}}
int FindBeadType2(char *name, int types_of_beads, BEADTYPE *BeadType);

// FindMoleculeType() //{{{
/** \brief Function to identify type of molecule from its name
 *
 * \param [in]  name          molecule name
 * \param [in]  Counts        numbers of beads, residues, etc.
 * \param [in]  MoleculeType  information about bead types
 * \return molecule type      id corresponding to index in BeadType struct (or -1 for non-existent molecule)
 */
int FindMoleculeType(char *name, COUNTS Counts, MOLECULETYPE *MoleculeType); //}}}
int FindMoleculeType2(char *name, int number_of_types, MOLECULETYPE *MoleculeType);

// Distancet() //{{{
/**
 * \brief Function to calculate distance vector between two beads.
 *
 * \param [in] id1         first coordinate vector
 * \param [in] id2         second coordinate vector
 * \param [in] BoxLength   dimensions of simulation box
 * \return distance vector between the two provided beads (without pbc)
 */
VECTOR Distance(VECTOR id1, VECTOR id2, VECTOR BoxLength); //}}}

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
void RemovePBCMolecules(COUNTS Counts, VECTOR BoxLength,
                        BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}

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
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNTS Counts,
                         VECTOR BoxLength, BEADTYPE *BeadType, BEAD **Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE *Molecule); //}}}

// RestorePBC() //{{{
/**
 * \brief Function to restore pbc.
 *
 * \param [in]  Counts         numbers of beads, molecules, etc.
 * \param [in]  BoxLength      dimension of the simulation box
 * \param [out] Bead           information about individual beads (coordinates)
 */
void RestorePBC(COUNTS Counts, VECTOR BoxLength, BEAD **Bead); //}}}
// RestorePBC2() //{{{
/**
 * \brief Function to restore pbc.
 *
 * \param [in]  number_of_beads   number of beads in the system
 * \param [in]  BoxLength         dimension of the simulation box
 * \param [out] Bead              information about individual beads (coordinates)
 */
void RestorePBC2(int number_of_beads, VECTOR BoxLength, BEAD **Bead); //}}}

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
VECTOR CentreOfMass(int n, int *list, BEAD *Bead, BEADTYPE *BeadType); //}}}

// GeomCentre() //{{{
/**
 * \brief Function to calculate geometric centre for a collection of beads.
 *
 * \param [in] n              number of beads
 * \param [in] list           list of bead ids (corresponding to indices in Bead struct)
 * \param [in] Bead           information about individual beads (coordinates)
 * \return coordinates of geometric centre of a given aggregate
 */
VECTOR GeomCentre(int n, int *list, BEAD *Bead); //}}}

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
VECTOR Gyration(int n, int *list, COUNTS Counts, VECTOR BoxLength,
                BEADTYPE *BeadType, BEAD **Bead); //}}}

// EvaluateContacts() //{{{
/**
 * \brief Function evaluating contacts for aggregate detection
 *
 * \param [in]  Counts        numbers of beads, molecules, etc.
 * \param [out] Aggregate     information about aggregates
 * \param [in]  Molecule      information about individual molecules
 * \param [in]  contacts      number of contacts for aggregate check
 * \param [in]  contact       2D array containing number of contacts between molecules
 */
void EvaluateContacts(COUNTS *Counts, AGGREGATE **Aggregate,
                      MOLECULE **Molecule,
                      int contacts, int **contact);
//}}}

// SortAggStruct() //{{{
/**
 * \brief Function to sort Aggregate struct.
 *
 * \param [out] Aggregate  Aggregate struct to sort
 * \param [in]  Counts     numbers of beads, molecules, etc.
 */
void SortAggStruct(AGGREGATE **Aggregate, COUNTS Counts,
                   MOLECULE *Molecule, MOLECULETYPE *MoleculeType,
                   BEAD **Bead, BEADTYPE *BeadType); //}}}

// LinkedList() //{{{
void LinkedList(VECTOR BoxLength, COUNTS Counts, BEAD *Bead,
                int **Head, int **Link, double cell_size, INTVECTOR *n_cells,
                int *Dcx, int *Dcy, int *Dcz); //}}}

// SortBonds() //{{{
/**
 * \brief Function to sort 2D array of bonds.
 *
 * \param [out] bond    2D array of bonds
 * \param [in]  length  number of bonds
 */
void SortBonds(int **bond, int length); //}}}

// SortAngles() //{{{
/**
 * \brief Function to sort 2D array of bonds.
 *
 * \param [out] angle   2D array of angles
 * \param [in]  length  number of angles
 */
void SortAngles(int **angle, int length); //}}}

// FreeBead() //{{{
/**
 * \brief Free memory allocated for Bead struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Bead        information about individual beads
 */
void FreeBead(COUNTS Counts, BEAD **Bead); //}}}
void FreeBead2(int number_of_beads, BEAD **Bead);

// FreeMolecule() //{{{
/**
 * \brief Free memory allocated for Molecule struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Molecule    information about individual molecules
 */
void FreeMolecule(COUNTS Counts, MOLECULE **Molecule); //}}}
void FreeMolecule2(int number_of_molecules, MOLECULE **Molecule);

// FreeMoleculeType() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  Counts         number of beads, molecu.es, etc.
 * \param [out] MoleculeType   information about individual molecules
 */
void FreeMoleculeType(COUNTS Counts, MOLECULETYPE **MoleculeType); //}}}
void FreeMoleculeType2(int number_of_types, MOLECULETYPE **MoleculeType);

// FreeAggregate() //{{{
/**
 * \brief Free memory allocated for MoleculeType struct array.
 *
 * \param [in]  Counts      number of beads, molecu.es, etc.
 * \param [out] Aggregate   information about individual molecules
 */
void FreeAggregate(COUNTS Counts, AGGREGATE **Aggregate); //}}}
#endif
