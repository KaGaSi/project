#ifndef READWRITEFIELD_H
#define READWRITEFIELD_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

SYSTEM FieldRead(const char *file);
void WriteField(SYSTEM System, char *file_field, int argc, char **argv);

#endif

