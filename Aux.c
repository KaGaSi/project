#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CStructs.h"

int FindType(char *name, Counts Counts, BeadType *BeadType) { //{{{
  int type;

  // compare give 'name' with all known bead types & return bead type id
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (strcmp(name, BeadType[i].Name) == 0) {
      type = i;
      return (type);
    }
  }

  // name isn't in BeadType struct
  fprintf(stderr, "Bead type %s doesn't exist!\n", name);
  exit(1);
} //}}}
