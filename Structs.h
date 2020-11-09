/**
 * \file
 * \brief Structures for utilities
 */

#ifndef _STRUCTS_H_
#define _STRUCTS_H_

#define CHARGE 10000 // 'impossible' charge to define a given bead type has charge specified in an input file
#define MASS 0 // 'impossible' mass to define a given bead type has charge specified in an input file
#define RADIUS 0 // 'impossible' radius to define a given bead type has charge specified in an input file

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
      Aggregates, ///< total number of aggregates
      TypesOfBonds, ///< number of bond types; -1 if not read from anywhere
      TypesOfAngles; ///< number of bond types; -1 if not read from anywhere
} COUNTS;

// Initialize Counts
static const COUNTS InitCounts = {
  .TypesOfBeads = 0,
  .TypesOfMolecules = 0,
  .Beads = 0,
  .Bonded = 0,
  .Unbonded = 0,
  .BeadsInVsf = 0,
  .Molecules = 0,
  .Aggregates = 0,
  .TypesOfBonds = -1,
  .TypesOfAngles = -1,
}; //}}}

// struct Params //{{{
/**
 * \brief Parameters for various things (bonds and angles for now)
 */
typedef struct Params {
  double a, b;
} PARAMS; //}}}

// struct BeadType //{{{
/**
 * \brief Information about bead types.
 */
typedef struct BeadType {
  char Name[17]; ///< name of given bead type

  int Number; ///< number of beads of given type

  bool Use, ///< should bead type in .vcf file be used for calculation?
       Write; ///< should bead type in .vcf file be written to output .vcf?

  double Charge, ///< charge of every bead of given type
         Mass, ///< mass of every bead of given type
         Radius; ///< radius of every bead of the given type
} BEADTYPE; //}}}

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

  VECTOR Position; ///< cartesian coordinates of the bead

  bool Flag; ///< some flag for, e.g., use/not use
} BEAD; //}}}

// struct MoleculeType //{{{
/**
 * \brief Information about molecule types.
 */
typedef struct MoleculeType {
  char Name[17]; ///< name of given molecule type

  int Number, ///< number of molecules of given type
      nBeads, ///< number of beads in every molecule of given type
      *Bead, ///< ids of bead types of every molecule bead
      nBonds, ///< number of bonds in every molecule of given type
      **Bond, ///< pair of ids for every bond (with relative bead numbers from 0 to nBeads)
               // has to be sorted; size: [MoleculeType[i].Bond[2]
      nAngles, ///< number of Angles in every molecule of given type
      **Angle, ///< trio of ids for every angle (with relative bead numbers from 0 to nBeads)
               // has to be sorted; size: [MoleculeType[i].Angle[3]
      nBTypes, ///< number of bead types in every molecule of given type
      *BType; ///< ids of bead types in every molecule of given type (corresponds to indices in BeadType struct)

  double Mass, ///< total mass of every molecule of given type
         Charge; ///< total charge of every molecule of given type

  bool InVcf, ///< is molecule type in vcf file?
       Use, ///< should molecule type be used for calculation?
       Write; ///< should molecule type be used for calculation?
} MOLECULETYPE; //}}}

// struct Molecule //{{{
/**
 * \brief Information about every molecule.
 */
typedef struct Molecule {
  int Type, ///< type of molecule corresponding to index in MoleculeType struct
      *Bead, ///< ids of beads in the molecule
      Aggregate; ///< id of aggregate molecule is in (corresponding to index in Aggregate struct)
} MOLECULE; //}}}

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
} AGGREGATE; //}}}
#endif
