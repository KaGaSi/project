#ifndef OPTIONS_H
#define OPTIONS_H

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdbool.h>
#include "AnalysisTools.h"

// version/help printing and initial check of provided options
int OptionCheck(const int argc, char **argv, const int req, const int common,
                const int all, const bool check_extra,
                char opt[all][OPT_LENGTH], ...);
// print help for common options
void CommonHelp(const bool error, const int n,
                const char option[n][OPT_LENGTH]);
// detect options common for most utilities
COMMON_OPT CommonOptions(const int argc, char **argv, const SYS_FILES f);
// exclude specified molecule names (-x <mol name(s)>)
bool ExcludeOption(const int argc, char **argv, SYSTEM *System);
// tag which bead types to use (if not present, set to specified value)
bool BeadTypeOption(const int argc, char **argv, const char *opt,
                    const bool use, bool *flag, const SYSTEM System);
// tag which molecule types to use (if not present, set to specified value)
bool MoleculeTypeOption(const int argc, char **argv, const char *opt,
                        const bool use, bool *flag, const SYSTEM System);
// general boolean option
bool BoolOption(const int argc, char **argv, const char *opt);
// general option with multiple integer arguments (up to 'max')
bool IntegerOption(const int argc, char **argv, const int max,
                   const char *opt, int *count, int *values);
bool IntegerOption1(const int argc, char **argv, const char *opt, int *value);
bool IntegerOption2(const int argc, char **argv,
                    const char *opt, int value[2]);
// general option with multiple double arguments (up to 'max')
bool DoubleOption(const int argc, char **argv, const int max,
                  const char *opt, int *count, double values[max]);
bool DoubleOption1(const int argc, char **argv,
                   const char *opt, double *value);
bool DoubleOption2(const int argc, char **argv,
                   const char *opt, double value[2]);
bool DoubleOption3(const int argc, char **argv,
                   const char *opt, double value[3]);
// general option with filename and integer(s) arguments
bool FileIntegerOption(const int argc, char **argv, const int min,
                       const int max, const char *opt, int *values,
                       int *count, char *file);
bool FileOption(const int argc, char **argv, const char *opt, char *file);
// general option with filename and double(s) arguments
bool FileDoubleOption(const int argc, char **argv, const int max,
                      const char *opt, double *values, int *count, char *file);
// print help - function body in each utility
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]);

#if 0 //{{{
// TODO redo
bool MoleculeTypeOption(int argc, char *argv[], char opt[], int *moltype,
                        COUNTS counts, MOLECULETYPE **MoleculeType);
bool MoleculeTypeOption2(int argc, char *argv[], char opt[], int *moltype,
                         COUNTS Counts, MOLECULETYPE **MoleculeType);
bool MoleculeTypeIntOption(int argc, int i, char *argv[], char opt[],
                           int *moltype, int *value, COUNTS Counts,
                           MOLECULETYPE *MoleculeType);
#endif //}}}
#endif
