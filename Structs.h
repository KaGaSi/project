/**
 * \file
 * \brief Structures for utilities
 */

#ifndef _STRUCTS_H_
#define _STRUCTS_H_

#define PI 3.141593 ///< value of pi

#define SQR(x) ((x)*(x)) ///< macro for algebraic square
#define CUBE(x) ((x)*(x)*(x)) ///< macro for algebraic cube

// struct Vector //{{{
/**
 * \brief 3D vector of floats.
 */
typedef struct Vector {
  double x, y, z;
} Vector; //}}}

// struct LongVector //{{{
/**
 * \brief 3D vector of floats.
 */
typedef struct LongVector {
  long double x, y, z;
} LongVector; //}}}

// struct IntVector //{{{
/**
 * \brief 3D vector of integers.
 */
typedef struct IntVector {
  int x, y, z;
} IntVector; //}}}

// struct Counts //{{{
/**
 * \brief Total numbers of various things.
 */
typedef struct Counts {
  int TypesOfBeads, ///< number of bead types
      TypesOfMolecules, ///< number of molecule types
      Beads, ///< total number of beads in all molecules
      Bonded, ///< total number of beads in all molecules (TO BE REMOVED)
      Unbonded, ///< total number of monomeric beads (TO BE REMOVED)
      BeadsInVsf, ///< total number of all beads in .vsf file (not necessarily in .vcf)
      Molecules, ///< total number of molecules
      MoleculesInVsf, ///< total number of all molecules in .vsf file (not necessarily in .vcf)
      Aggregates; ///< total number of aggregates
} Counts; //}}}

// struct BeadType //{{{
/**
 * \brief Information about bead types.
 */
typedef struct BeadType {
  char Name[16]; ///< name of given bead type

  int Number; ///< number of beads of given type

  bool Use, ///< should bead type in .vcf file be used for calculation?
       Write; ///< should bead type in .vcf file be written to output .vcf?

  double Charge, ///< charge of every bead of given type
         Mass; ///< mass of every bead of given type
} BeadType; //}}}

// struct MoleculeType //{{{
/**
 * \brief Information about molecule types.
 */
typedef struct MoleculeType {
  char Name[16]; ///< name of given molecule type

  int Number, ///< number of molecules of given type
      nBeads, ///< number of beads in every molecule of given type
      nBonds, ///< number of bonds in every molecule of given type
      **Bond, ///< pair of ids for every bond (with relative bead numbers from 0 to nBeads)
               // has to be sorted; size: [MoleculeType[i].Bonds][2]
      nBTypes, ///< number of bead types in every molecule of given type
      *BType; ///< ids of bead types in every molecule of given type (corresponds to indices in BeadType struct)

  double Mass; ///< total mass of every molecule of given type

  bool InVcf, ///< is molecule type in vcf file?
       Use, ///< should molecule type be used for calculation?
       Write; ///< should molecule type be used for calculation?
} MoleculeType; //}}}

// struct Bead //{{{
/**
 * \brief Information about every bead.
 */
typedef struct Bead {
  int Type, ///< type of bead corresponding to index in BeadType struct
      Molecule, ///< index number of molecule corresponding to Molecule struct (-1 for monomeric bead)
      nAggregates, ///< number of aggregates the bead is in (only monomeric beads can be in more aggregates - allocated memory for 10)
      *Aggregate, ///< index numbers of aggregates corresponding to Aggregate struct (-1 for bead in no aggregate)
      Index; ///< index of the bead according to .vsf file (needed for indexed timesteps)

  Vector Position; ///< cartesian coordinates of the bead
} Bead; //}}}

// struct Molecule //{{{
/**
 * \brief Information about every molecule.
 */
typedef struct Molecule {
  int Type, ///< type of molecule corresponding to index in MoleculeType struct
      *Bead, ///< ids of beads in the molecule
      Aggregate; ///< id of aggregate molecule is in (corresponding to index in Aggregate struct)
} Molecule; //}}}

// struct Aggregate //{{{
/**
 * \brief Information about every aggregate.
 */
typedef struct Aggregate {
  int nMolecules, ///< number of molecules in aggregate
      *Molecule, ///< ids of molecules in aggregate
      nBeads, ///< number of bonded beads in aggregate
      *Bead, ///< ids of bonded beads in aggregate
      nMonomers, ///< number of monomeric beads in aggregate
      *Monomer; ///< ids of monomeric beads in aggregate

  double Mass; ///< total mass of the aggregate

  bool Use; ///< should aggregate be used for calculation?
} Aggregate; //}}}
#endif
