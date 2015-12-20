#ifndef _ANALYSISTOOLS_H_
#define _ANALYSISTOOLS_H_

#include <stdbool.h>

/**
 * \file
 * \brief Structures used by all analysis programs
 */

// struct Vector //{{{
/**
 * \brief 3D vector
 */
typedef struct Vector {
  double x, y, z;
} Vector; //}}}

// struct Counts //{{{
/**
 * \brief Total numbers of various things */
typedef struct Counts {
  int TypesOfBeads, ///< number of bead types
      TypesOfMolecules, ///< number of molecule types
      Bonded, ///< total number of beads in all molecules
      Unbonded, ///< total number of monomeric beads
      Molecules; ///< total number of molecules
} Counts; //}}}

// struct BeadType //{{{
/**
 * \brief Information about bead types
 */
typedef struct BeadType {
  char Name[16]; ///< name of given bead type

  int Number; ///< number of beads of given type

  double Charge, ///< charge of every bead of given type
         Mass; ///< mass of every bead of given type

  bool Use; ///< indication if beads of given type are in coordinate .vcf file
} BeadType; //}}}

// struct MoleculeType //{{{
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

// struct Bead //{{{
/**
 * \brief Information about every bead
 */
typedef struct Bead {
  int Type; ///< type of bead corresponding to index in BeadType struct

  Vector Position; ///< cartesian coordinates of the bead
} Bead; //}}}

// struct Molecule //{{{
/**
 * \brief Information about every molecule
 */
typedef struct Molecule {
  int Type, ///< type of molecule corresponding to index in MoleculeType struct
      *Bead; ///< ids of beads in the molecule
} Molecule; //}}}

// ReadStructure() //{{{
/**
 * \brief Function reading information from dl_meso FIELD and vsf
 * structure files.
 *
 * \param [in]  vsf_file      .vsf structure file
 * \param [in]  bonds_file    filename with bonds
 * \param [out] Counts        numbers of beads, molecules, etc.
 * \param [out] BeadType      information about bead types
 * \param [out] Bead          informationn about individual beads
 * \param [out] MoleculeType  information about molecule types
 * \param [out] Molecule      information about individual molecules
 * */
void ReadStructure(char *vsf_file, char *bonds_file, Counts *Counts,
                   BeadType **BeadType, Bead **Bead,
                   MoleculeType **MoleculeType, Molecule **Molecule); //}}}

// WriteVsf() //{{{
/**
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

// ReadCoorOrdered() //{{{
/**
 * \brief Function reading ordered coordinates from .vcf coordinate file.
 *
 * \param [in]  vcf_file   name of input .vcf coordinate file
 * \param [in]  Counts     numbers of beads, molecules, etc.
 * \param [out] Bead       coordinates of individual beads
 * \param [out] stuff      first two line of a timestep
 * \return 0 for no errors or index number of bead (starting from 1) for which coordinates cannot be read
 */
int ReadCoorOrdered(FILE *vcf_file, Counts Counts, Bead **Bead, char **stuff); //}}}

// WriteCoorIndexed //{{{
/**
 * \brief Function writing indexed coordinates to a .vcf file.
 *
 * \param [in] vcf_file   name of output .vcf coordinate file
 * \param [in] Counts     numbers of beads, molecules, etc.
 * \param [in] BeadType   information about bead types
 * \param [in] Bead       coordinates of individual beads
 * \param [in] stuff      array of chars containing comment line to place at the beginning
 */
void WriteCoorIndexed(FILE *vcf_file, Counts Counts, BeadType *BeadType, Bead *Bead, char *stuff); //}}}

// FindType() //{{{
/** \brief Function to identify type of bead from its name
 *
 * \param [in]  name      bead name
 * \param [in]  Counts    numbers of beads, residues, etc.
 * \param [in]  BeadType  informationn about bead types
 * \return bead type id corresponding to index in BeadType struct
 */
int FindType(char *name, Counts Counts, BeadType *BeadType); //}}}
#endif
