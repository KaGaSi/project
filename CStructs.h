#ifndef _CSTRUCTS_H_
#define _CSTRUCTS_H_
/**
 * \file
 * \brief Structures used by all analysis programs
 */

// Vector //{{{
/**
 * \brief 3D vector
 */
typedef struct Vector {
  double x, y, z;
} Vector; //}}}

// Counts //{{{
/**
 * \brief Total numbers of various things */
typedef struct Counts {
  int TypesOfBeads, ///< number of bead types
      TypesOfMolecules, ///< number of molecule types
      Bonded, ///< total number of beads in all molecules
      Unbonded, ///< total number of monomeric beads
      Molecules; ///< total number of molecules
} Counts; //}}}

// BeadType //{{{
/**
 * \brief Information about bead types
 */
typedef struct BeadType {
  char Name[16]; ///< name of given bead type

  int Number; ///< number of beads of given type

  double Charge, ///< charge of every bead of given type
         Mass; ///< mass of every bead of given type
} BeadType; //}}}

// MoleculeType //{{{
/**
 * \brief Information about molecule types
 */
typedef struct MoleculeType {
  char Name[16]; ///< name of given molecule type

  int Number, ///< number of molecules of given type
      nBeads, ///< number of beads in every molecule of given type
      nBonds, ///< number of bonds in every molecule of given type
      **Bond; ///< pair of ids for every bond (with relative bead numbers from 0 to nBeads)
               // has to be sorted; size: [MoleculeType[i].Bonds][2]
} MoleculeType; //}}}

// Bead //{{{
/**
 * \brief Information about every bead
 */
typedef struct Bead {
  int Type; ///< type of bead corresponding to index in BeadType struct

  Vector Position; ///< cartesian coordinates of the bead
} Bead; //}}}

// Molecule //{{{
/**
 * \brief Information about every molecule
 */
typedef struct Molecule {
  int Type, ///< type of molecule corresponding to index in MoleculeType struct
      *Bead; ///< ids of beads in the molecule
} Molecule; //}}}
#endif
