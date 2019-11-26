/**
 * \file
 * \brief Options usable in utilities
 */

#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include "Structs.h"

// Help() //{{{
/**
 * \brief Function to print help.
 *
 * \param [in]  cmd    utility name
 * \param [out] error  is it help due to error or not?
 */
void Help(char cmd[50], bool error); //}}}

// CommonOptions() //{{{
/**
 * \brief Function for options common to most of the utilities.
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [out] input.vsf    .vsf structure file
 * \param [out] verbose      verbose output?
 * \param [out] silent       no output?
 * \apram [out] script       output to file?
 */
void CommonOptions(int argc, char **argv, char **vsf_file,
                   bool *verbose, bool *silent, bool *script); //}}}

// VerboseLongOption() //{{{
/**
 * \brief Option whether to use long verbose output (overrides
 * VerboseShortOutput) (`-V`).
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [out] verbose      bool for `-v` option (verbose output)
 * \return `true` or `false` for error on common options
 */
void VerboseLongOption(int argc, char **argv, bool *verbose); //}}}

// SilentOption() //{{{
/**
 * \brief Option whether not to print to stdout (overrides Verbose options)
 * (`-s`).
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [out] verbose      bool for `-v` option (verbose output)
 * \param [out] silent       bool for this option
 * \return `true` or `false` for error on common options
 */
void SilentOption(int argc, char **argv, bool *verbose, bool *silent); //}}}

// ExcludeOption() //{{{
/**
 * \brief Option whether to exclude molecule types (`-x <name(s)>`).
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [in]  Counts       numbers of beads, molecules, etc.
 * \param [out] MoleculeType information about molecule types
 * \return `true` or `false` error or not error
 */
bool ExcludeOption(int argc, char **argv, Counts Counts,
                   MoleculeType **MoleculeType); //}}}

// JoinCoorOption() //{{{
/**
 * \brief Option whether to write joined aggregate coordinates to file (`-j
 * <joined.vcf>`).
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [in]  Counts       numbers of beads, molecules, etc.
 * \param [out] MoleculeType information about molecule types
 * \return `true` or `false` error or not error
 */
bool JoinCoorOption(int argc, char **argv, char *joined_vcf); //}}}

// BeadTypeOption() //{{{
/**
 * \brief Option to choose which bead types to use in calculations (`-bt
 * <name(s)>`).
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [in]  opt          option switich (e.g., string '-bt')
 * \param [in]  use          if the option is not present, set all BeadType[].Use flags to 'use'
 * \param [in]  Counts       numbers of beads, molecules, etc.
 * \param [out] BeadType     information about bead types
 * \return `true` or `false` error or not error
 */
bool BeadTypeOption(int argc, char **argv, char *opt, bool use,
                    Counts Counts, BeadType **BeadType); //}}}

// BoolOption() //{{{
/**
 * \brief Option whether not to print rewrite stdout line (`--script`).
 *
 * \param [in] argc  number of program's arguments
 * \param [in] argv  program's arguments
 * \param [in] opt   option switch (e.g. array containing `-n`)
 * \return `true` if `opt` present, `false` otherwise
 */
bool BoolOption(int argc, char **argv, char *opt); //}}}

// IntegerOption() //{{{
/**
 * \brief Function for any option with integer argument.
 *
 * \param [in]  argc  number of program's arguments
 * \param [in]  argv  program's arguments
 * \param [in]  opt   option switch (e.g. array containing `-n`)
 * \param [out] value integer value of given option
 * \return `true` or `false` for error
 */
bool IntegerOption(int argc, char **argv, char *opt, int *value);
// }}}

// DoubleOption() //{{{
/**
 * \brief Function for any option with integer argument.
 *
 * \param [in]  argc  number of program's arguments
 * \param [in]  argv  program's arguments
 * \param [in]  opt   option switch (e.g. array containing `-n`)
 * \param [out] value double value of given option
 * \return `true` or `false` for error
 */
bool DoubleOption(int argc, char **argv, char *opt, double *value);
// }}}

// MultiIntegerOption() //{{{
/**
 * \brief Function for any option with two integer arguments.
 *
 * \param [in]  argc   number of program's arguments
 * \param [in]  argv   program's arguments
 * \param [in]  opt    option switch (e.g. array containing `-n`)
 * \param [out] count  number of numeric arguments
 * \param [out] values array of integer values of given option
 * \return `true` or `false` for error
 */
bool MultiIntegerOption(int argc, char **argv, char *opt, int *count, int *values);
// }}}

// MultiDoubleOption() //{{{
/**
 * \brief Function for any option with two integer arguments.
 *
 * \param [in]  argc   number of program's arguments
 * \param [in]  argv   program's arguments
 * \param [in]  opt    option switch (e.g. array containing `-n`)
 * \param [out] count  number of numeric arguments
 * \param [out] values array of double values of given option
 * \return `true` or `false` for error
 */
bool MultiDoubleOption(int argc, char **argv, char *opt, int *count, double *values);
// }}}

// FileIntsOption() //{{{
/**
 * \brief Function for any option with filename and up to a hundred
 * integer arguments.
 *
 * \param [in]  argc   number of program's arguments
 * \param [in]  argv   program's arguments
 * \param [in]  opt    option switch (e.g. array containing `-c`)
 * \param [out] values array of two integer values of given option
 * \param [out] count  number of numeric arguments
 * \param [out] file   file name (first argument of option)
 * \return `true` or `false` for error
 */
bool FileIntsOption(int argc, char **argv, char *opt, int *values, int *count, char *file);
 //}}}

// FileOption() //{{{
/**
 * \brief Function for any option with filename.
 *
 * \param [in]  argc  number of program's arguments
 * \param [in]  argv  program's arguments
 * \param [in]  opt   option switch (e.g. array containing `-n`)
 * \param [out] name  array containing the filename
 * \return `true` or `false` for error
 */
bool FileOption(int argc, char **argv, char *opt, char **name); //}}}

// MoleculeTypeOption() //{{{
/**
 * \brief Function for any option with molecule type name.
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [in]  opt          option switch (e.g. array containing `-n`)
 * \param [in]  Counts       numbers of beads, molecules, etc.
 * \param [out] moltype      id of the molecule type
 * \param [in]  MoleculeType information about molecule types
 * \return `true` or `false` for error
 */
bool MoleculeTypeOption(int argc, char **argv, char *opt, int *moltype, Counts Counts,
                        MoleculeType **MoleculeType); //}}}

// MoleculeTypeOption2() //{{{
/**
 * \brief Function for any option with molecule type name(s).
 *
 * \param [in]  argc         number of program's arguments
 * \param [in]  argv         program's arguments
 * \param [in]  opt          option switch (e.g. array containing `-n`)
 * \param [out] moltype      array for the molecule type
 * \param [in]  Counts       numbers of beads, molecules, etc.
 * \param [in]  MoleculeType information about molecule types
 * \return `true` or `false` for error
 */
bool MoleculeTypeOption2(int argc, char **argv, char *opt, int **moltype, Counts Counts,
                         MoleculeType **MoleculeType); //}}}
#endif
