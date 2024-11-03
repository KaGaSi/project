#ifndef READWRITE_H
#define READWRITE_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

// reading functions
void CopyMoleculeTypeBeadsToMoleculeBeads(SYSTEM *System);
// fill MoleculeType[].Bond arrays with data from a separate bond array
void FillMoleculeTypeBonds(SYSTEM *System, const int (*bond)[3],
                           const int n);
// remove molecule/bead types with 0 molecules/beads
void RemoveExtraTypes(SYSTEM *System);
// SYSTEM ReadStructure(int struct_type, char struct_file[], int coor_type,
//                      char coor_file[], bool detailed);
SYSTEM ReadStructure(SYS_FILES f, bool detailed);
bool ReadTimestep(SYS_FILES f, FILE *fr, SYSTEM *System, int *line_count);
bool SkipTimestep(SYS_FILES f, FILE *fr, int *line_count);
int ReadAggregates(FILE *fr, const char *file, SYSTEM *System,
                   AGGREGATE *Aggregate, int *line_count);
// fill MoleculeType[].Angle arrays with data from a separate angle array
void FillMoleculeTypeAngles(SYSTEM *System, const int (*angle)[4],
                            const int n);
// fill MoleculeType[].Dihedral arrays with data from a separate dihedral array
void FillMoleculeTypeDihedral(SYSTEM *System, const int (*dihedral)[5],
                              const int n);
// fill MoleculeType[].Improper arrays with data from a separate improper array
void FillMoleculeTypeImproper(SYSTEM *System, const int (*improper)[5],
                              const int n);
// writing functions
void InitCoorFile(FILE_TYPE fout, SYSTEM System, int argc, char **argv);
void WriteOutput(SYSTEM System, const bool *write, FILE_TYPE fw,
                 bool lmp_mass, int vsf_def, int argc, char **argv);
void WriteOutputAll(SYSTEM System, FILE_TYPE fw, bool lmp_mass,
                    int vsf_def, int argc, char **argv);
void WriteTimestep(FILE_TYPE f, SYSTEM System, int count_step,
                   const bool *write, int argc, char **argv);
void WriteTimestepAll(FILE_TYPE f, SYSTEM System, int count_step,
                      int argc, char **argv);
void WriteStructure(FILE_TYPE f, SYSTEM System, int vsf_def_type,
                    bool lmp_mass, int argc, char **argv);
void WriteAggregates(int step_count, char *agg_file, SYSTEM System,
                     AGGREGATE *Aggregate);

// verbose output (print various structures and some such)
void VerboseOutput(SYSTEM System);
void PrintCount(COUNT Count);
void PrintBeadType(SYSTEM System);
void PrintOneMolType(SYSTEM System, int n);
void PrintAllMolTypes(SYSTEM System);
void Print1Molecule(SYSTEM System, int n);
void PrintMolecules(SYSTEM System);
void PrintBead(SYSTEM System);
void PrintBondType(SYSTEM System);
void PrintAngleType(SYSTEM System);
void PrintDihedralType(SYSTEM System);
void PrintImproperType(SYSTEM System);
// TODO: use SYSTEM
void PrintBondTypes(COUNT Counts, PARAMS *bond_type);
// TODO: use SYSTEM
void PrintAngleTypes(COUNT Counts, PARAMS *angle_type);
void PrintBox(BOX Box);
void PrintByline(char *file, int argc, char **argv);
void PrintStep(int *count_coor, int start, bool silent);
void PrintAggregate(SYSTEM System, AGGREGATE *Aggregate);
#endif
