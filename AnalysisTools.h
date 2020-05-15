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
#include "General.h"
#include "Errors.h"
#include "Structs.h"
#include "Options.h"
#include "Read.h"

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
void VerboseOutput(char *input_vcf, Counts Counts, Vector BoxLength,
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

// PrintAggregate() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintAggregate(Counts Counts, int *Index, MoleculeType *MoleculeType, Molecule *Molecule, Bead *Bead, BeadType *BeadType, Aggregate *Aggregate); //}}}

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
void EvaluateContacts(Counts *Counts, Aggregate **Aggregate,
                      Molecule **Molecule,
                      int contacts, int **contact);
//}}}

// SortAggStruct() //{{{
/**
 * \brief Function to sort Aggregate struct.
 *
 * \param [out] Aggregate  Aggregate struct to sort
 * \param [in]  Counts     numbers of beads, molecules, etc.
 */
void SortAggStruct(Aggregate **Aggregate, Counts Counts); //}}}

// LinkedList() //{{{
void LinkedList(Vector BoxLength, Counts Counts, Bead *Bead,
                int **Head, int **Link, double cell_size, IntVector *n_cells,
                int *Dcx, int *Dcy, int *Dcz); //}}}

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
void ZeroCounts(Counts *Counts); //}}}

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
