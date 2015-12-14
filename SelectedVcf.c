#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CStructs.h"
#include "Structure.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   TransformVsf <output.vsf> <options>\n\n");
  fprintf(stderr, "   <start>         number of timestep to start from\n");
  fprintf(stderr, "   <skip>          leave-out every 'skip' steps\n");
  fprintf(stderr, "   <output.vcf>    output filename (vcf format)\n");
  fprintf(stderr, "   <type names> name of bead types to save\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -v   verbose output\n");
  fprintf(stderr, "      -h   print this help and exit\n");
  fprintf(stderr, "      -i   use input .vsf file different from dl_meso.vsf\n\n");
} //}}}
