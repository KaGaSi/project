/**
 * \file
 * \brief Options usable in utilities
 */

#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include "Structs.h"

// ExcludeOption() //{{{
bool ExcludeOption(int argc, char **argv, Counts Counts,
                   MoleculeType **MoleculeType); //}}}

bool JoinCoorOption(int argc, char **argv, char *joined_vcf);
#endif
