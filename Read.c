#include "Read.h"

void SkipStructVtf(FILE *vtf, char *name_vtf) { //{{{
  char line[LINE];
  while (fgets(line, sizeof(line), vtf)) {
    char split[30][100], delim[8];
    strcpy(delim, " \t");
    int words = SplitLine(split, line, delim);
    // only certain keywords besides pbc can be present before the first coordinate block
    // 1) t(imestep) ... starting the coordinate block
    // 2) t(imestep) i(ndexed)/o(ordered) ... starting the coordinate block
    // 3) i(ndexed)/o(ordered) ... starting the coordinate block
    // 4) a(tom) ... in case of vtf file
    // 5) b(ond) ... in case of vtf file
    // 6) empty line
    // 7) comment
    // 8) pbc line
    if (!(split[0][0] == 't' && words == 1) && // 1)
        !(split[0][0] == 't' && words > 1 && (split[1][0] == 'o' || split[1][0] =='i')) && // 2)
        !(split[0][0] == 'o' || split[0][0] == 'i') && // 3)
        split[0][0] != 'a' && // 4)
        split[0][0] != 'b' && // 5)
        split[0][0] != '\0' && // 6)
        split[0][0] != '#' && // 7)
        strcasecmp(split[0], "pbc") != 0) { // 8)
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError:");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", name_vtf);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - unrecognised line in structure section (or first timestep's preamble)\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }

    // save file pointer at the last bond line
    fpos_t position;
    if (split[0][0] == 'b') { // 5)
      // save pointer (i.e., possible last bond line)
      fgetpos(vtf, &position);
    }
    // restore file pointer to after the last bond line if timestep line encountered
    if ((split[0][0] == 't' && words == 1) || // 1)
        (split[0][0] == 't' && words > 1 && (split[1][0] == 'o' || split[1][0] =='i')) || // 2)
        (split[0][0] == 'o' || split[0][0] == 'i')) { // 3)
      fsetpos(vtf, &position); // restore pointer position
      break;
    }
  };
} //}}}

// GetPBC() //{{{
/*
 * Function to get box dimensions from the provided coordinate file.
 */
VECTOR GetPBC(FILE *vcf, char *vcf_file) {

  VECTOR BoxLength;

  char line[LINE];
  while (fgets(line, sizeof(line), vcf)) {
    char split[30][100], delim[8];
    strcpy(delim, " \t");
    int words = SplitLine(split, line, delim);

    if (strcasecmp(split[0], "pbc") == 0) {
      BoxLength.x = atof(split[1]);
      BoxLength.y = atof(split[2]);
      BoxLength.z = atof(split[3]);
      break;
    // only certain keywords besides pbc can be present before the first coordinate block
    // 1) t(imestep) ... starting the coordinate block
    // 2) t(imestep) i(ndexed)/o(ordered) ... starting the coordinate block
    // 3) i(ndexed)/o(ordered) ... starting the coordinate block
    // 4) a(tom) ... in case of vtf file
    // 5) b(ond) ... in case of vtf file
    // 6) empty line
    // 7) comment
    } else if (!(split[0][0] == 't' && words == 1) && // 1)
               !(split[0][0] == 't' && words > 1 && (split[1][0] == 'o' || split[1][0] =='i')) && // 2)
               !(split[0][0] == 'o' || split[0][0] == 'i') && // 3)
               split[0][0] != 'a' && // 4)
               split[0][0] != 'b' && // 5)
               split[0][0] != '\0' && // 6)
               split[0][0] != '#') { // 7)
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", vcf_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - unrecognised line\n");
      ErrorPrintLine(split, words);
      RedText(STDERR_FILENO);
      putc('\n', stderr);
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  };
  return BoxLength;
} //}}}

// GetPBC2() //{{{
/*
 * Function to get box dimensions from the provided coordinate file.
 */
VECTOR GetPBC2(char *vcf_file) {

  VECTOR BoxLength;

  FILE *vcf;
  if ((vcf = fopen(vcf_file, "r")) == NULL) {
    ErrorFileOpen(vcf_file, 'r');
    exit(1);
  }
  char line[LINE];
  while (fgets(line, sizeof(line), vcf)) {
    char split[30][100], delim[8];
    strcpy(delim, " \t");
    int words = SplitLine(split, line, delim);

    if (strcasecmp(split[0], "pbc") == 0) {
      BoxLength.x = atof(split[1]);
      BoxLength.y = atof(split[2]);
      BoxLength.z = atof(split[3]);
      break;
    // only certain keywords besides pbc can be present before the first coordinate block
    // 1) t(imestep) ... starting the coordinate block
    // 2) t(imestep) i(ndexed)/o(ordered) ... starting the coordinate block
    // 3) i(ndexed)/o(ordered) ... starting the coordinate block
    // 4) a(tom) ... in case of vtf file
    // 5) b(ond) ... in case of vtf file
    // 6) empty line
    // 7) comment
    } else if (!(split[0][0] == 't' && words == 1) && // 1)
               !(split[0][0] == 't' && words > 1 && (split[1][0] == 'o' || split[1][0] =='i')) && // 2)
               !(split[0][0] == 'o' || split[0][0] == 'i') && // 3)
               split[0][0] != 'a' && // 4)
               split[0][0] != 'b' && // 5)
               split[0][0] != '\0' && // 6)
               split[0][0] != '#') { // 7)
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", vcf_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - unrecognised line\n");
      ErrorPrintLine(split, words);
      RedText(STDERR_FILENO);
      putc('\n', stderr);
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  };

  fclose(vcf);

  return BoxLength;
} //}}}

// ReadAggCommand() //{{{
void ReadAggCommand(BEADTYPE *BeadType, COUNTS Counts,
                    char *input_coor, char *input_agg,
                    double *distance, int *contacts) {

  // open input aggregate file
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  // read first line (Aggregate command)
  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t");
  fgets(line, sizeof(line), agg);
  int words = SplitLine(split, line, delim);
  // errof if not enough striings for an Aggregate command //{{{
  if (words < 6) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - first line must contain a valid Aggregates command\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  } //}}}
  // read minimum distance for closeness check (<distance> argument in Aggregates utility) //{{{
  if (!IsPosDouble(split[2])) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - <distance> in Aggregate command must be a non-negative real number\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  }
  *distance = atof(split[2]); //}}}
  // read number of contacts aggregate check (<contacts> argument in Aggregates utility) //{{{
  if (!IsInteger(split[3])) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - <contacts> in Aggregate command must be a non-negative whole number\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  }
  *contacts = atof(split[3]); //}}}
  // warn if a differently named vcf file is used than the one in agg file //{{{
  if (strcmp(split[1], input_coor) != 0) {
    YellowText(STDERR_FILENO);
    fprintf(stderr, "\nWarning: coordinate file ");
    CyanText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    YellowText(STDERR_FILENO);
    fprintf(stderr, " is different to the one in the aggregate file (");
    CyanText(STDERR_FILENO);
    fprintf(stderr, "%s", split[1]);
    YellowText(STDERR_FILENO);
    fprintf(stderr, ")\n");
    fprintf(stderr, "         Mismatch between beads present in both files can lead to undefined behaviour.\n\n");
    ResetColour(STDERR_FILENO);
  } //}}}
  // read <type names> from Aggregates command //{{{
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType(split[i], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m -", input_agg);
      fprintf(stderr, " - non-existent bead name \033[1;33m%s\033[1;31m in Aggregate command\n", split[i]);
      fprintf(stderr, "\033[0m");
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}
  fclose(agg);
} //}}}

// CountVtfStructLines() //{{{
// Assumes vtf; counts number of structure part lines (i.e., the number of the
// last line containing atom/bond).
int CountVtfStructLines(bool vtf, char *input) {
  if (vtf) {
    // open input file
    FILE *vtf;
    if ((vtf = fopen(input, "r")) == NULL) {
      ErrorFileOpen(input, 'r');
      exit(1);
    }
    /*
     * i) check for errors
     * ii) find number of beads
     * iii) save bead and molecule names
     */
    int count_lines = 0, last_line;
    while(true) {
      char split[30][100], line[LINE];
      fgets(line, sizeof line, vtf);
      if (feof(vtf)) { // exit while loop on vsf file finish
        break;
      }
      int words = SplitLine(split, line, "\t ");
      // finish reading at the 'timestep' line
      if (CheckVtfTimestepLine(words, split)) {
        break;
      }
      count_lines++;
      // skip blank, comment, and pbc lines
      if (words == 0 || split[0][0] == '#' || strcasecmp(split[0], "pbc") == 0) {
        continue;
      }
      if (strncmp(split[0], "atom", 1) == 0 ||
          strncmp(split[0], "bond", 1) == 0) {
        last_line = count_lines;
      }
    }
    fclose(vtf);
    return last_line;
  } else {
    return -1; // if not vtf, the number doesn't matter
  }
} //}}}

// SkipVtfStructure() //{{{
void SkipVtfStructure(bool vtf, FILE *vcf, int struct_lines) {
  // skip structure part of a vtf file
  if (vtf) {
    for (int i = 0; i < struct_lines; i++) {
      char line[LINE];
      fgets(line, sizeof line, vcf);
    }
  }
} //}}}

// CheckVtfTimestepLine() //{{{
/*
 * Function to check if the provided line is a timestep line.
 */
bool CheckVtfTimestepLine(int words, char split[30][100]) {
  // there are several possibilities how the timestep line looks
  // 1) 't[imestep]'
  if (words == 1 && strncasecmp(split[0], "timestep", 1) == 0) {
    return true;
  }
  // 2) 't[imestep] i[indexed]/o[rdered] any extra stuff'
  if (words > 1 && strncasecmp(split[0], "timestep", 1) == 0 &&
      (strncasecmp(split[1], "indexed", 1) == 0 || strncasecmp(split[1], "ordered", 1) == 0)) {
    return true;
  }
  // 3) 'i[indexed]/o[rdered] any extra stuff'
  if (words > 0 &&
      (strncasecmp(split[0], "indexed", 1) == 0 || strncasecmp(split[0], "ordered", 1) == 0)) {
    return true;
  }
  return false;
} //}}}

// CheckVtfAtomLine() //{{{
// the line must be: a[tom] <id>/'default' ... n[ame] <char(16)>... [resid <char> ... resname <char>]
bool CheckVtfAtomLine(int words, char split[30][100], char *error) {
  // is there even number of strings?
  if ((words%2) != 0) {
    strcpy(error, "odd number of strings in atom line");
    return false;
  }
  // is the line starting with a[tom] <id>/default?
  if (strcmp(split[1], "default") != 0 && !IsInteger(split[1])) {
    strcpy(error, "'atom' must be followed by <int> or 'default'");
    return false;
  }
  // is the mandatory n[ame] keyword present?
  bool name = false;
  for (int i = 2; i < words; i+=2) {
    if (split[i][0] == 'n') {
      name = true;
    }
  }
  if (!name) {
    strcpy(error, "missing 'name' in atom line");
    return false;
  }
  // TODO: some default naming if resname is missing
  // if resid is present, there must be resname as well
  for (int i = 2; i < words; i+=2) {
    if (strcmp(split[i], "resid") == 0) {
      // check that resid is followed by <int>
      if (!IsInteger(split[i+1])) {
        strcpy(error, "'resid' must be followed by <int>");
        return false;
      }
      // search for res[name]
      name = false;
      for (int j = 2; j < words; j+=2) {
        if (strncmp(split[j], "resname", 3) == 0 && j != i) { // 'resid' also satisfies the strncmp()
          name = true;
          break;
        }
      }
      if (!name) {
        strcpy(error, "if 'resid' is present, 'resname' must be too");
        return false;
      }
      break;
    }
  }
  // if mass, charge, or radius keywoards are present, the must be followed by <float>
  for (int i = 2; i < words; i+=2) {
    if ((strcmp(split[i], "charge") == 0 || split[i][0] == 'q') && !IsDouble(split[i+1])) {
      strcpy(error, "'charge|q' must be followed by <float>");
      return false;
    } else if (split[i][0] == 'r' && strncmp(split[i], "res", 3) != 0 && (!IsPosDouble(split[i+1]) || atof(split[i+1]) <= 0)) {
      strcpy(error, "bead radius must be positive");
      return false;
    } else if (split[i][0] == 'm' && (!IsPosDouble(split[i+1]) || atof(split[i+1]) <= 0)) {
      strcpy(error, "bead mass must positive");
      return false;
    }
  }
  return true;
} //}}}

// CheckVtfBondLine() //{{{
// the line must be: b[ond] <id>:<id> (with possible blanks around ':') 
bool CheckVtfBondLine(int words, char split[30][100], char *error) {
  switch(words) {
    case 1: // there must be at least 1 more string after 'bond'
      strcpy(error, "missing a string after 'bond'");
      return false;
    case 2: // if there's exactly 1 string after 'bond', it must be exactly '<int>:<int>'
      for (int i = 1; i < (strlen(split[1])-2); i++) {
        if ((split[1][i] == ':' && split[1][i+1] == ':') ||
            split[1][i] == ',' || split[1][i+1] == ',') {
          strcpy(error, "for now, only bond lines in the format <int>:<int> are allowed (i.e., no comma separated list or 'two-colon' string of bonds)");
          return false;
        }
      }
      char index[30][100], string[100];
      strcpy(string, split[1]);
      int strings = SplitLine(index, string, ":");
      if (strings < 2 || !IsInteger(index[0]) || !IsInteger(index[1]) ||
          split[1][0] == ':' || split[1][strlen(split[1])-1] == ':') { // check for preceding/ending colon stripped away by SplitLine()
        strcpy(error, "missing two colon separated integers in a bond line");
        return false;
      }
      return true;
    case 3: // if there are exactly 2 strings after 'bond', it must be exactly '<int>: <int>'
      if (split[1][strlen(split[1])-1] == ':' && IsInteger(split[2])) {
        if (split[1][strlen(split[1])-2] == ':') { // two colons aren't allowed
          strcpy(error, "for now, only bond lines in the format <int>:<int> are allowed (i.e., no comma separated list or 'two-colon' string of bonds)");
          return false;
        }
        // test that <int> preceeds the colon
        strcpy(string, split[1]);
        string[strlen(string)-1] = '\0';
        if (!IsInteger(string)) {
          strcpy(error, "missing an integer before colon in a bond line");
          return false;
        }
      } else {
        strcpy(error, "unrecognised bond line (note that whitespace is not allowed before the colon; i.e., '<int> : <int>' is not allowed)");
        return false;
      }
      return true;
    default: // more then 2 strings after 'bond' isn't allowed (at least for now)
      strcpy(error, "unrecognised bond line (note that for now, only the format <int>:<int> is allowed; i.e., no comma separated list or 'two-colon' string of bonds)");
      return false;
  }
} //}}}

// CheckVtfCoordinateLine() //{{{
// the line must be: <x> <y> <z>, preceded by <int> for indexed timestep
bool CheckVtfCoordinateLine(int words, char split[30][100], char *error, bool indexed) {
  // not enough strings
  if ((words < 3 && !indexed) ||
      (words < 4 && indexed)) {
    strcpy(error, "wrong coordinate line");
    return false;
  }
  // incorrect string types
  if (indexed && (!IsInteger(split[0]) || !IsDouble(split[1]) || !IsDouble(split[2]) || !IsDouble(split[3]))) {
    strcpy(error, "wrong coordinate line for an indexed timestep");
    return false;
  }
  if (!indexed && (!IsDouble(split[0]) || !IsDouble(split[1]) || !IsDouble(split[2]))) {
    strcpy(error, "wrong coordinate line for an ordered timestep");
    return false;
  }
  // correct line
  return true;
} //}}}

// CopyBeadType() //{{{
void CopyBeadType(int number_of_types, BEADTYPE **bt_out, BEADTYPE *bt_in) {
  for (int i = 0; i < number_of_types; i++) {
    strcpy((*bt_out)[i].Name, bt_in[i].Name);
    (*bt_out)[i].Number = bt_in[i].Number;
    (*bt_out)[i].Charge = bt_in[i].Charge;
    (*bt_out)[i].Mass = bt_in[i].Mass;
    (*bt_out)[i].Radius = bt_in[i].Radius;
    (*bt_out)[i].Use = bt_in[i].Use;
    (*bt_out)[i].Write = bt_in[i].Write;
  }
} //}}}

// CopyMoleculeType() //{{{
void CopyMoleculeType(int number_of_types, MOLECULETYPE **mt_out, MOLECULETYPE *mt_in) {
  FreeMoleculeType2(number_of_types, mt_out);
  *mt_out = calloc(number_of_types, sizeof(MOLECULETYPE));
  for (int i = 0; i < number_of_types; i++) {
    strcpy((*mt_out)[i].Name, mt_in[i].Name);
    (*mt_out)[i].Number = mt_in[i].Number;
    (*mt_out)[i].Charge = mt_in[i].Charge;
    (*mt_out)[i].Mass = mt_in[i].Mass;
    (*mt_out)[i].Use = mt_in[i].Use;
    (*mt_out)[i].Write = mt_in[i].Write;
    (*mt_out)[i].nBeads = mt_in[i].nBeads;
    (*mt_out)[i].Bead = calloc((*mt_out)[i].nBeads, sizeof(int));
    for (int j = 0; j < (*mt_out)[i].nBeads; j++) {
      (*mt_out)[i].Bead[j] = mt_in[i].Bead[j];
    }
    (*mt_out)[i].nBTypes = mt_in[i].nBTypes;
    (*mt_out)[i].BType = calloc((*mt_out)[i].nBTypes, sizeof(int));
    for (int j = 0; j < (*mt_out)[i].nBTypes; j++) {
      (*mt_out)[i].BType[j] = mt_in[i].BType[j];
    }
    if (mt_in[i].nBonds > 0) {
      (*mt_out)[i].nBonds = mt_in[i].nBonds;
      (*mt_out)[i].Bond = calloc((*mt_out)[i].nBonds, sizeof(int *));
      for (int j = 0; j < (*mt_out)[i].nBonds; j++) {
        (*mt_out)[i].Bond[j] = calloc(3, sizeof(int));
        (*mt_out)[i].Bond[j][0] = mt_in[i].Bond[j][0];
        (*mt_out)[i].Bond[j][1] = mt_in[i].Bond[j][1];
        (*mt_out)[i].Bond[j][2] = mt_in[i].Bond[j][2];
      }
    }
    if (mt_in[i].nAngles > 0) {
      (*mt_out)[i].nAngles = mt_in[i].nAngles;
      (*mt_out)[i].Angle = calloc((*mt_out)[i].nAngles, sizeof(int *));
      for (int j = 0; j < (*mt_out)[i].nBonds; j++) {
        (*mt_out)[i].Angle[j] = calloc(4, sizeof(int));
        (*mt_out)[i].Angle[j][0] = mt_in[i].Angle[j][0];
        (*mt_out)[i].Angle[j][1] = mt_in[i].Angle[j][1];
        (*mt_out)[i].Angle[j][2] = mt_in[i].Angle[j][2];
        (*mt_out)[i].Angle[j][3] = mt_in[i].Angle[j][3];
      }
    }
  }
} //}}}

// CopyMolecule() //{{{
void CopyMolecule(int number_of_molecules, MOLECULETYPE *mt, MOLECULE **m_out, MOLECULE *m_in) {
  for (int i = 0; i < number_of_molecules; i++) {
    (*m_out)[i].Type = m_in[i].Type;
    (*m_out)[i].Aggregate = m_in[i].Aggregate;
    int mtype = m_in[i].Type;
    (*m_out)[i].Bead = realloc((*m_out)[i].Bead, mt[mtype].nBeads*sizeof(int));
    for (int j = 0; j < mt[mtype].nBeads; j++) {
      (*m_out)[i].Bead[j] = m_in[i].Bead[j];
    }
  }
} //}}}

// CopyBead() //{{{
void CopyBead(int number_of_beads, BEAD **b_out, BEAD *b_in) {
  for (int i = 0; i < number_of_beads; i++) {
    (*b_out)[i].Type = b_in[i].Type;
    (*b_out)[i].nAggregates = b_in[i].nAggregates;
    (*b_out)[i].Index = b_in[i].Index;
    (*b_out)[i].Position.x = b_in[i].Position.x;
    (*b_out)[i].Position.y = b_in[i].Position.y;
    (*b_out)[i].Position.z = b_in[i].Position.z;
    (*b_out)[i].Flag = b_in[i].Flag;

    // TODO: this useless stuff will be removed at some point
    if ((*b_out)[i].nAggregates > 0) {
      (*b_out)[i].Aggregate = realloc((*b_out)[i].Aggregate, 10*sizeof(int));
      for (int j = 0; j < 10; j++) {
        (*b_out)[i].Aggregate[j] = b_in[i].Aggregate[j];
      }
    }
  }
} //}}}

// NewBeadType() //{{{
/*
 * Function to create a new bead type in a BEADTYPE struct.
 */
void NewBeadType(char *name, int *number_of_types, BEADTYPE **BeadType, char *vsf_file) {
  int btype = (*number_of_types)++;
  *BeadType = realloc(*BeadType, (*number_of_types)*sizeof(struct BeadType));
  (*BeadType)[btype].Number = 0;
  // TODO: make a safe copy function
  // copy new name to BeadType[].Name
  if (strlen(name) > (sizeof (*BeadType)[btype].Name)) {
    name[sizeof (*BeadType)[btype].Name - 1] = '\0';
  }
  strcpy((*BeadType)[btype].Name, name);
  // initialize charge and mass by 'impossible' value
  (*BeadType)[btype].Charge = 1000;
  (*BeadType)[btype].Mass = 0;
}; //}}}

// NewMolType() //{{{
/*
 * Function to create a new bead type in a BEADTYPE struct.
 */
void NewMolType(char *name, int *number_of_types, MOLECULETYPE **MoleculeType, char *vsf_file) {
  int mtype = (*number_of_types)++;
  *MoleculeType = realloc(*MoleculeType, (*number_of_types)*sizeof(struct MoleculeType));
  // TODO: make a safe copy function
  // copy new name to MoleculeType[].Name
  if (strlen(name) > (sizeof (*MoleculeType)[mtype].Name)) {
    name[sizeof (*MoleculeType)[mtype].Name - 1] = '\0';
  }
  strcpy((*MoleculeType)[mtype].Name, name);
  // initialize struct members
  (*MoleculeType)[mtype].Number = 0;
  (*MoleculeType)[mtype].nBonds = 0;
  (*MoleculeType)[mtype].nBeads = 0;
  (*MoleculeType)[mtype].nBTypes = 0;
}; //}}}

// FillMolBTypes //{{{
/*
 * Function to fill MoleculeType[].BType array based on MoleculeType[].Bead array.
 */
void FillMolBTypes(int number_of_types, MOLECULETYPE **MoleculeType) {
  for (int i = 0; i < number_of_types; i++) {
    (*MoleculeType)[i].nBTypes = 0;
    (*MoleculeType)[i].BType = calloc(1, sizeof(int));
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      bool new = true;
      for (int k = 0; k < (*MoleculeType)[i].nBTypes; k++) {
        if ((*MoleculeType)[i].Bead[j] == (*MoleculeType)[i].BType[k]) {
          new = false;
          break;
        }
      }
      if (new) {
        int type = (*MoleculeType)[i].nBTypes++;
        (*MoleculeType)[i].BType = realloc((*MoleculeType)[i].BType, (*MoleculeType)[i].nBTypes*sizeof(int));
        (*MoleculeType)[i].BType[type] = (*MoleculeType)[i].Bead[j];
      }
    }
  }
} //}}}

// FillMolMass //{{{
/*
 * Function to calculate mass of all molecules. If at least one bead has
 * undefined mass, the mass of the molecule is also undefined
 */
void FillMolMass(int number_of_types, BEADTYPE *BeadType, MOLECULETYPE **MoleculeType) {
  for (int i = 0; i < number_of_types; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int type = (*MoleculeType)[i].Bead[j];
      if (BeadType[type].Mass == MASS) {
        (*MoleculeType)[i].Mass = MASS;
        break;
      } else {
        (*MoleculeType)[i].Mass += BeadType[type].Mass;
      }
    }
    if ((*MoleculeType)[i].Mass == MASS) {
      fprintf(stderr, "\033[1;33m");
      fprintf(stderr, "\nWarning: molecule type \033[1;36m%s\033[1;33m has undefined mass\n", (*MoleculeType)[i].Name);
      fprintf(stderr, "\033[0m");
    }
  }
} //}}}

// FillMolCharge //{{{
/*
 * Function to calculate charge of all molecules. If at least one bead has
 * undefined charge, the charge of the molecule is also undefined
 */
void FillMolCharge(int number_of_types, BEADTYPE *BeadType, MOLECULETYPE **MoleculeType) {
  for (int i = 0; i < number_of_types; i++) {
    (*MoleculeType)[i].Charge = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int type = (*MoleculeType)[i].Bead[j];
      if (BeadType[type].Charge == CHARGE) {
        (*MoleculeType)[i].Charge = CHARGE;
        break;
      } else {
        (*MoleculeType)[i].Charge += BeadType[type].Charge;
      }
    }
    if ((*MoleculeType)[i].Charge == CHARGE) {
      fprintf(stderr, "\033[1;33m");
      fprintf(stderr, "\nWarning: molecule type \033[1;36m%s\033[1;33m has undefined charge\n", (*MoleculeType)[i].Name);
      fprintf(stderr, "\033[0m");
    }
  }
} //}}}

// FillMolType //{{{
/*
 * Function to fill BType array and mass and charge for each molecule type
 */
void FillMolType(int number_of_types, BEADTYPE *BeadType, MOLECULETYPE **MoleculeType) {
  FillMolBTypes(number_of_types, MoleculeType);
  FillMolMass(number_of_types, BeadType, MoleculeType);
  FillMolCharge(number_of_types, BeadType, MoleculeType);
} //}}}

// CheckVtfTimestep() //{{{
bool CheckVtfTimestep(FILE *vcf, char *vcf_file, COUNTS *Counts,
                      BEADTYPE **BeadType, BEAD **Bead, int **Index,
                      MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {

//PrintMoleculeType2((*Counts).TypesOfMolecules, *BeadType, *MoleculeType);
  bool indexed; // is the timestep indexed? ...returned by this function
  // skip timestep preamble, determining if timesteps are ordered/indexed //{{{
  char *stuff = calloc(LINE, sizeof(char));
  int lines = ReadTimestepPreamble(&indexed, vcf_file, vcf, &stuff, true); // number of preamble lines
  for (int i = 0; i < lines; i++) {
    char line[100];
    fgets(line, sizeof line, vcf);
  } //}}}
  // count beads in the timestep & save their ids //{{{
  lines = 0;
  // use Bead[].Flag to determine which beads are in the timestep (true=present)
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    (*Bead)[i].Flag = false;
  }
  while (true) {
    char split[30][100], line[LINE];
    fgets(line, sizeof line, vcf);
    if (feof(vcf)) { // exit while loop on coordinate file finish
      break;
    }
    int words = SplitLine(split, line, "\t ");
    // test if the line is a proper coordinate line
    if (!CheckVtfCoordinateLine(words, split, stuff, indexed)) {
      // exit while loop when the first line of the next timestep preamble is found
      if (words == 0 || split[0][0] == '#' || strcasecmp(split[0], "pbc") == 0) {
        break;
      } else { // error if not proper coordinate line & not the start of a new timestep
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", vcf_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - %s\n", stuff);
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
    }
    int id;
    if (indexed) { // indexed timestep => 1st is bead id
      id = (*Index)[atoi(split[0])];
    } else { // ordered timestep => bead id is based on the line count
      id = lines;
    }
    (*Bead)[id].Flag = true;
    lines++;
  }
  free(stuff); //}}}
  // error - if ordered timestep, all beads must be present //{{{
  if (!indexed && lines != (*Counts).Beads) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", vcf_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - 'ordered' timestep requires all beads present in the coordinate file\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // if all beads are in the coordinate file, nothing more to do //{{{
  if (lines == (*Counts).Beads) {
    return indexed;
  } //}}}
  // otherwise, identify what's present according to the saved Bead[].Flag
  COUNTS c_new = InitCounts; // initialize new COUNTS struct
  c_new.Beads = lines; // number of coordinate lines in the coordinate file
  c_new.BeadsInVsf = (*Counts).BeadsInVsf;
  // copy beads present in the coordinate file to a new BEAD struct & Index array //{{{
  int count = 0;
  BEAD *b_new = calloc(c_new.Beads, sizeof(BEAD));
  int *index_new = calloc((*Counts).BeadsInVsf, sizeof(int));
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if ((*Bead)[i].Flag) {
      b_new[count].Type = (*Bead)[i].Type;
      b_new[count].Molecule = (*Bead)[i].Molecule;
      b_new[count].Index = (*Bead)[i].Index;
      b_new[count].nAggregates = 0;
      index_new[b_new[count].Index] = count;
      // count bonded & unbonded beads
      if (b_new[count].Molecule == -1) {
        c_new.Unbonded++;
      } else {
        c_new.Bonded++;
      }
      count++;
    }
  } //}}}
  // count numbers of beads of different types present in the coordinate file //{{{
  int *numbers = calloc((*Counts).TypesOfBeads, sizeof(int));
  for (int i = 0; i < c_new.Beads; i++) {
    int btype = b_new[i].Type;
    numbers[btype]++;
  }
  // error - not all beads of a given type are in the vcf file
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if (numbers[i] != 0 && numbers[i] != (*BeadType)[i].Number) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", vcf_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - not all ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", (*BeadType)[i].Name);
      RedText(STDERR_FILENO);
      fprintf(stderr, " beads are in the coordinate file\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
    if (numbers[i] > 0) {
      c_new.TypesOfBeads++;
    }
  } //}}}
  // copy bead types present in the coordinate file to a new struct //{{{
  BEADTYPE *bt_new = calloc(c_new.TypesOfBeads, sizeof(BEADTYPE));
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if (numbers[i] != 0) {
      strcpy(bt_new[count].Name, (*BeadType)[i].Name);
      bt_new[count].Number = (*BeadType)[i].Number;
      bt_new[count].Charge = (*BeadType)[i].Charge;
      bt_new[count].Mass = (*BeadType)[i].Mass;
      bt_new[count].Radius = (*BeadType)[i].Radius;
      count++;
    }
  }
  free(numbers);
  c_new.TypesOfBeads = count; //}}}
  // update bead types the new BEAD according to the new BEADTYPE //{{{
  for (int i = 0; i < c_new.Beads; i++) {
    int btype = b_new[i].Type;
    btype = FindBeadType2((*BeadType)[btype].Name, c_new.TypesOfBeads, bt_new);
    b_new[i].Type = btype;
  } //}}}
  // determine molecule types (and number of molecules) present in the vcf //{{{
  /*
   * A molecule type is present in the 'new' system if it contains at least one
   * bead type present in the coordinate file.
   * i) count number of molecule types and of molecules
   * ii) copy the molecule types present in the coordinate file to a new MOLECULETYPE struct
   * iii) determine what bead types are present in the 'new' system (only parts
   *      of a molecule may be in the vcf, so adjust bonds too)
   *      TODO: in future, vcf coordinate file might be combined with a
   *            different structure filetype, so angles should in principle be
   *            adjusted too
   * iv) determine BTypes, charge & mass of the 'new' molecule types (taking
   *     into account that only parts of the original molecules may be present)
   */
  // i) //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int old_type = (*MoleculeType)[i].BType[j];
      if (FindBeadType2((*BeadType)[old_type].Name, c_new.TypesOfBeads, bt_new) != -1) {
        c_new.TypesOfMolecules++;
        c_new.Molecules += (*MoleculeType)[i].Number;
        break;
      }
    }
  } //}}}
  // ii) & iii) //{{{
  MOLECULETYPE *mt_new = calloc(c_new.TypesOfMolecules, sizeof(MOLECULETYPE));
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int old_type = (*MoleculeType)[i].BType[j];
      if (FindBeadType2((*BeadType)[old_type].Name, c_new.TypesOfBeads, bt_new) != -1) {
        strcpy(mt_new[count].Name, (*MoleculeType)[i].Name);
        mt_new[count].Number = (*MoleculeType)[i].Number;
        // copy beads present in the vcf to the new MOLECULETYPE.Bead array
        mt_new[count].nBeads = 0;
        mt_new[count].Bead = calloc(1, sizeof(int));
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) { // go through all original types, testing each one
          int btype = FindBeadType2((*BeadType)[(*MoleculeType)[i].Bead[k]].Name, c_new.TypesOfBeads, bt_new);
          if (btype != -1) {
            int id = mt_new[count].nBeads++;
            mt_new[count].Bead = realloc(mt_new[count].Bead, mt_new[count].nBeads*sizeof(int));
            mt_new[count].Bead[id] = btype;
          }
        }
        // copy bonds so that only those between beads present in the coordinate file remain
        mt_new[count].nBonds = 0;
        mt_new[count].Bond = calloc(1, sizeof(int *));
        for (int k = 0; k < (*MoleculeType)[i].nBonds; k++) { // go through all original types, testing each pair of beads
          int old_type1 = (*MoleculeType)[i].Bead[(*MoleculeType)[i].Bond[k][0]];
          int old_type2 = (*MoleculeType)[i].Bead[(*MoleculeType)[i].Bond[k][1]];
          int btype1 = FindBeadType2((*BeadType)[old_type1].Name, c_new.TypesOfBeads, bt_new);
          int btype2 = FindBeadType2((*BeadType)[old_type2].Name, c_new.TypesOfBeads, bt_new);
          if (btype1 != -1 && btype2 != -1) {
            int id = mt_new[count].nBonds++;
            mt_new[count].Bond = realloc(mt_new[count].Bond, mt_new[count].nBonds*sizeof(int *));
            mt_new[count].Bond[id] = calloc(3,sizeof(int));
            // these indices correspond to the old MOLECULETYPE that contains all beads - will be adjusted later
            mt_new[count].Bond[id][0] = (*MoleculeType)[i].Bond[k][0];
            mt_new[count].Bond[id][1] = (*MoleculeType)[i].Bond[k][1];
            mt_new[count].Bond[id][2] = (*MoleculeType)[i].Bond[k][2];
          }
        }
        if (mt_new[count].nBonds == 0) { // if there are no bonds, free the array
          free(mt_new[count].Bond);
        }
        // adjust the bonded beads' indices based on the numbering in the new MOLECULETYPE //{{{
        /*
         * i) count how much to subtract from old indices in the new
         *    MOLECULETYPE's Bond array by going through the old MOLECULETYPE's
         *    index of bead types and testing if they are in the coordinate file
         * ii) subtract the required amount from each index in every bond
         */
        int subtract[(*MoleculeType)[i].nBeads];
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          subtract[k] = 0;
        }
        // i)
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          int old_type = (*MoleculeType)[i].Bead[k];
          int btype = FindBeadType2((*BeadType)[old_type].Name, c_new.TypesOfBeads, bt_new);
          if (btype == -1) {
            for (int l = k; l < (*MoleculeType)[i].nBeads; l++) {
              subtract[l]++;
            }
          }
        }
        // ii)
        for (int k = 0; k < mt_new[count].nBonds; k++) {
          int id1 = mt_new[count].Bond[k][0];
          int id2 = mt_new[count].Bond[k][1];
          mt_new[count].Bond[k][0] -= subtract[id1];
          mt_new[count].Bond[k][1] -= subtract[id2];
        } //}}}
        count++; // increment number of molecule types copied
        break;
      }
    }
  } //}}}
  // iv)
  FillMolBTypes(c_new.TypesOfMolecules, &mt_new);
  FillMolMass(c_new.TypesOfMolecules, bt_new, &mt_new);
  FillMolCharge(c_new.TypesOfMolecules, bt_new, &mt_new);
  //}}}
  // copy molecules present in the coordinate file to a new MOLECULE struct //{{{
  MOLECULE *m_new = calloc(c_new.Molecules, sizeof(MOLECULE));
  count = 0; // count molecules
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int old_mtype = (*Molecule)[i].Type;
    int new_mtype = FindMoleculeType2((*MoleculeType)[old_mtype].Name, c_new.TypesOfMolecules, mt_new);
    if (new_mtype != -1) {
      m_new[count].Type = new_mtype;
      m_new[count].Bead = calloc(mt_new[new_mtype].nBeads, sizeof(int));
      int count2 = 0; // counts beads in count-th molecule
      for (int j = 0; j < (*MoleculeType)[old_mtype].nBeads; j++) {
        int old_btype = (*MoleculeType)[old_mtype].Bead[j];
        if (FindBeadType2((*BeadType)[old_btype].Name, c_new.TypesOfBeads, bt_new) != -1) {
          int id = (*Molecule)[i].Bead[j];
          int new_id = index_new[(*Bead)[id].Index];
          m_new[count].Bead[count2] = new_id;
          b_new[new_id].Molecule = count;
          count2++;
        }
      }
      count++;
    }
  } //}}}
  // copy data back to original structs & realloc the structs //{{{
  // BEAD struct - realloc and copy; copy Index array
  *Bead = realloc(*Bead, c_new.Beads*sizeof(BEAD));
  for (int i = 0; i < c_new.Beads; i++) {
    (*Bead)[i].Type = b_new[i].Type;
    (*Bead)[i].Molecule = b_new[i].Molecule;
    (*Bead)[i].Index = b_new[i].Index;
    (*Index)[(*Bead)[i].Index] = i;
  }
  // BEADTYPE - realloc and copy
  *BeadType = realloc(*BeadType, c_new.TypesOfBeads*sizeof(BEADTYPE));
  CopyBeadType(c_new.TypesOfBeads, BeadType, bt_new);
  // MOLECULETYPE - free memory struct[].arrays for extra molecule types, realloc MOLECULETYPE, and copy
  if (c_new.TypesOfMolecules < (*Counts).TypesOfMolecules) {
    for (int i = c_new.TypesOfMolecules; i < (*Counts).TypesOfMolecules; i++) {
      for (int j = 0; (*MoleculeType)[i].nBonds > 0 && j < (*MoleculeType)[i].nBonds; j++) {
        free((*MoleculeType)[i].Bond[j]);
      }
      for (int j = 0; (*MoleculeType)[i].nAngles > 0 && j < (*MoleculeType)[i].nAngles; j++) {
        free((*MoleculeType)[i].Angle[j]);
      }
      free((*MoleculeType)[i].Angle);
      free((*MoleculeType)[i].BType);
      free((*MoleculeType)[i].Bead);
      free((*MoleculeType)[i].Bond);
    }
    *MoleculeType = realloc(*MoleculeType, c_new.TypesOfMolecules*sizeof(MOLECULETYPE));
  }
  CopyMoleculeType(c_new.TypesOfMolecules, MoleculeType, mt_new);
  // MOLECULE - free struct[].arrays for extra molecules, realloc MOLECULE, and copy
  if (c_new.Molecules < (*Counts).Molecules) {
    for (int i = c_new.Molecules; i < (*Counts).Molecules; i++) {
      free((*Molecule)[i].Bead);
    }
    *Molecule = realloc(*Molecule, c_new.Molecules*sizeof(MOLECULE));
  }
  CopyMolecule(c_new.Molecules, *MoleculeType, Molecule, m_new);
  // copy the new Counts struct back to the original one
  *Counts = c_new; //}}}
  // free memory //{{{
  free(b_new);
  free(bt_new);
  free(index_new);
  FreeMoleculeType2((*Counts).TypesOfMolecules, &mt_new);
  FreeMolecule2((*Counts).Molecules, &m_new); //}}}
  return indexed;
} //}}}

// ReadVtfStructure() //{{{
void ReadVtfStructure(char *vsf_file, bool detailed, COUNTS *Counts,
                      BEADTYPE **BeadType, BEAD **Bead, int **Index,
                      MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {

  (*Counts) = InitCounts; // zeroize
  // test if vsf_file is vtf or vsf //{{{
  bool vtf = false;
  char *dot = strrchr(vsf_file, '.');
  if (strcmp(dot, ".vtf") == 0) {
    vtf = true;
  } //}}}
  // open vsf_file //{{{
  FILE *vsf;
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}
  fpos_t pos;
  fgetpos(vsf, &pos); // save file pointer
  // 1) read through the structure part to find basic info: //{{{
  /*
   * i) check for errors
   * ii) find number of beads
   * iii) save bead and molecule names
   */
  int count_atom_lines = 0, // number of atom lines
      default_atom_line = -1, // if there is a default atom line, what line is it (the first one, if multiple are present)?
      count_bond_lines = 0, // number of bond lines
      count_comment_lines = 0, // number of comments (#) or blank lines
      atom_names = 0, res_names = 0, // number of different bead and molecule names
      last_struct = 0; // number of the last structure line (i.e., bond or atom); used in case of vtf
  char **atom_name = calloc(1, sizeof(char *)); // names of different atoms
  atom_name[0] = calloc(17, sizeof(char));
  char **res_name = calloc(1, sizeof(char *)); // names of different molecules
  res_name[0] = calloc(9, sizeof(char));
  while(true) {
    char error[LINE] = {'\0'}, // error message; if (strlen(error) == 0), then no error
         split[30][100], line[LINE];
    fgets(line, sizeof line, vsf);
    if (feof(vsf)) { // exit while loop on vsf file finish
      break;
    }
    int words = SplitLine(split, line, "\t ");
    // if the structure part is in a vtf file, finish reading at the 'timestep' line
    if (vtf && CheckVtfTimestepLine(words, split)) {
      break;
    }
    // skip blank, comment, and pbc lines
    if (words == 0 || split[0][0] == '#' || strcasecmp(split[0], "pbc") == 0) {
      count_comment_lines++;
      continue;
    }
    switch(split[0][0]) {
      case 'a': // atom line //{{{
        if (CheckVtfAtomLine(words, split, error)) {
          count_atom_lines++;
          last_struct = count_atom_lines + count_bond_lines + count_comment_lines;
        } else {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError: ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", vsf_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - %s\n", error);
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        }
        if (strcmp(split[1], "default") == 0) {
          if (default_atom_line != -1) {
            fprintf(stderr, "\033[1;33m");
            fprintf(stderr, "\nWarning: \033[1;36m%s\033[1;33m - multiple 'atom default' lines\n", vsf_file);
            fprintf(stderr, "         Using line number \033[1;36m%d\033[1;33\n", default_atom_line+1);
            fprintf(stderr, "\033[0m");
          } else {
            default_atom_line = count_atom_lines + count_comment_lines + count_bond_lines - 1; // so it starts from 0
          }
        } else {
          // check for highest index (i.e., the number of beads in vsf)
          if (atoi(split[1]) >= (*Counts).BeadsInVsf) {
            (*Counts).BeadsInVsf = atoi(split[1]) + 1; // bead indices in vtf start with 0
          }
        }
        // save bead & molecule names (if not saved already)
        for (int i = 0; i < words; i+= 2) {
          if (split[i][0] == 'n') {
            int in = -1;
            for (int j = 0; j  < atom_names; j++) {
              if (strncmp(split[i+1], atom_name[j], 15) == 0) {
                in = j;
                break;
              }
            }
            if (in == -1) {
              if (atom_names > 0) {
                atom_name = realloc(atom_name, (atom_names+1)*sizeof(char *));
                atom_name[atom_names] = calloc(16, sizeof(char));
              }
              strcpy(atom_name[atom_names], split[i+1]);
              atom_names++;
            }
          } else if (strncmp(split[i], "res", 3) == 0 && strcmp(split[i], "resid") != 0) {
            int in = -1;
            for (int j = 0; j  < res_names; j++) {
              if (strncmp(split[i+1], res_name[j], 7) == 0) {
                in = j;
                break;
              }
            }
            if (in == -1) {
              if (res_names > 0) {
                res_name = realloc(res_name, (res_names+1)*sizeof(char *));
                res_name[res_names] = calloc(9, sizeof(char));
              }
              strcpy(res_name[res_names], split[i+1]);
              res_names++;
            }
          }
        }
        break; //}}}
      case 'b': // bond line //{{{
        if (CheckVtfBondLine(words, split, error)) {
          count_bond_lines++;
          last_struct = count_atom_lines + count_bond_lines + count_comment_lines;
        } else {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError: ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", vsf_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - %s\n", error);
          ErrorPrintLine(split, words);
          exit(1);
        }
        break; //}}}
      default: // something unrecognised //{{{
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        RedText(STDERR_FILENO);
        fprintf(stderr, "%s", vsf_file);
        YellowText(STDERR_FILENO);
        fprintf(stderr, " - unrecognised line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1); //}}}
    }
  } //}}}
  // error - missing default line but have too few atom lines //{{{
  if (default_atom_line == -1 && count_atom_lines != (*Counts).BeadsInVsf) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", vsf_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing atom line(s)\n");
    fprintf(stderr, "       (i.e., when 'atom default' is omitted, there must be a line for each atom)\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // total number of lines in the structure section (vtf/vsf file dependent) //{{{
  int total_lines; total_lines = last_struct;
  if (vtf) {
    total_lines = last_struct;
  } else {
    total_lines = count_atom_lines + count_bond_lines + count_comment_lines;
  } //}}}
  // declare & allocate & initilize stuff for holding all vsf info //{{{
    // structure to save info from atom lines
    struct atom {
      int index, name, resid, resname;
      double charge, mass, radius;
    } *atom = calloc(count_atom_lines, sizeof(struct atom));
    for (int i = 0; i < count_atom_lines; i++) {
      atom[i].resid = -1; // i.e., not in a molecule
      atom[i].charge = CHARGE; // i.e., missing charge in the atom line
      atom[i].mass = MASS; // i.e., missing mass in the atom line
      atom[i].radius = RADIUS; // i.e., missing radius in the atom line
    }
    // array to connect atom ids in vsf with struct atom numbering
    int *atom_id = calloc((*Counts).BeadsInVsf, sizeof(int));
    for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
      atom_id[i] = -1; // if it stays -1, it's default atom
    }
    // structure to save info from bond lines
    struct bond {
      int index1, index2;
    } *bond = calloc(count_bond_lines, sizeof(struct bond)); //}}}
  // 2) save all vsf lines into atom & bond structs & and count number of molecules //{{{
  fsetpos(vsf, &pos); // restore file pointer (to the beginning of vsf_file)
  int count_atoms = 0, count_bonds = 0;
  for (int line_no = 0; line_no < total_lines; line_no++) {
    char split[30][100], line[LINE];
    fgets(line, sizeof line, vsf);
    int words = SplitLine(split, line, "\t ");
    // if the structure part is in a vtf file, finish reading at the 'timestep' line
    if (vtf && CheckVtfTimestepLine(words, split)) {
      break;
    }
    // skip blank, comment, and pbc lines, (i.e., count them as comment lines)
    if (words == 0 || split[0][0] == '#' || strcasecmp(split[0], "pbc") == 0) {
      continue;
    }
    switch (split[0][0]) {
      case 'a': // save atom line into atom struct //{{{
        if (strcmp(split[1], "default") != 0 || line_no == default_atom_line) { // use only the designated atom default line, not extra ones
          int id;
          if (line_no == default_atom_line) { // first atom default line in vsf
            id = count_atom_lines - 1;
            atom[id].index = -1;
          } else { // use atom <id> lines, but not extra atom default lines
            id = count_atoms++;
            atom_id[atoi(split[1])] = id;
            atom[id].index = atoi(split[1]);
          }
          for (int i = 0; i < words; i+=2) {
            // if the bead is in a molecule, count it as bonded
            if (strncmp(split[i], "resid", 3) == 0 && strcmp(split[i], "resname") != 0) {
              // find highest id molecule
              if (atoi(split[i+1]) > (*Counts).Molecules) {
                (*Counts).Molecules = atoi(split[i+1]);
              }
              atom[id].resid = atoi(split[i+1]) - 1; // resid in vsf starts with 1
            } else if (split[i][0] == 'n') {
              for (int j = 0; j < atom_names; j++) {
                if (strncmp(split[i+1], atom_name[j], 16) == 0) {
                  atom[id].name = j;
                  break;
                }
              }
            } else if (strcmp(split[i], "resname") == 0) {
              for (int j = 0; j < res_names; j++) {
                if (strncmp(split[i+1], res_name[j], 8) == 0) {
                  atom[id].resname = j;
                  break;
                }
              }
            } else if (strcmp(split[i], "charge") == 0 || strcmp(split[i], "q") == 0) {
              atom[id].charge = atof(split[i+1]);
            } else if (split[i][0] == 'm') {
              atom[id].mass = atof(split[i+1]);
            } else if (split[i][0] == 'r' && strncmp(split[i], "res", 3) != 0) {
              atom[id].radius = atof(split[i+1]);
            }
          }
        }
        break; //}}}
      case 'b': // save bond line into bond struct //{{{
        if (words == 2) {
          char index[30][100];
          SplitLine(index, split[1], ":");
          bond[count_bonds].index1 = atoi(index[0]);
          bond[count_bonds].index2 = atoi(index[1]);
        } else {
          split[1][strlen(split[1])-1] = '\0';
          bond[count_bonds].index1 = atoi(split[1]);
          bond[count_bonds].index2 = atoi(split[2]);
        }
        count_bonds++;
        break; //}}}
    }
  } //}}}
  fclose(vsf); // everything is saved already
  // fill atom_id with default beads //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if (atom_id[i] == -1) {
      atom_id[i] = count_atom_lines - 1;
    }
  } //}}}
  // 3) fill BEADTYPE struct //{{{
  BEADTYPE *bt = calloc(1, sizeof(BEADTYPE));
  bt[0].Charge = CHARGE; //
  bt[0].Mass = MASS;     // initial values for atom default
  bt[0].Radius = RADIUS; //
  // save default type - if there is one //{{{
  if (atom[count_atom_lines-1].index == -1) {
    (*Counts).TypesOfBeads = 1;
    bt = realloc(bt, (*Counts).TypesOfBeads*sizeof(BEADTYPE));
    strcpy(bt[0].Name, atom_name[atom[count_atom_lines-1].name]);
    bt[0].Charge = atom[count_atom_lines-1].charge;
    bt[0].Mass = atom[count_atom_lines-1].mass;
    bt[0].Radius = atom[count_atom_lines-1].radius;
  } //}}}
  if (detailed) { // check other stuff besides name //{{{
    /*
     * i) identify bead types according to all data not just name, e.g., lines
     *      atom 0 n x q 1 m 1
     *      atom 0 n x q 2 m 1
     *    will be of different types
     * ii) however, should keyword be missing in one line but present in
     *     another, that does not create a different type, e.g., lines
     *       atom 0 n x q 1 m 1
     *       atom 0 n x m 1
     *     will be of the same type
     * iii) however, there can be ambiguities, so e.g., lines
     *        atom 0 n x q 1 m 1
     *        atom 0 n x     m 1
     *        atom 0 n x q 0 m 1
     *      shouldn't be joined - see the merge part below
     */
    // go through all atoms (besides atom default line) //{{{
    for (int i = 0; i < count_atoms; i++) {
      int btype = -1;
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[j].Name, atom_name[atom[i].name]) == 0 &&
            bt[j].Charge == atom[i].charge &&
            bt[j].Mass == atom[i].mass &&
            bt[j].Radius == atom[i].radius) {
          btype = j;
        }
      }
      if (btype == -1) { // new bead type
        btype = (*Counts).TypesOfBeads++;
        bt = realloc(bt, (*Counts).TypesOfBeads*sizeof(BEADTYPE));
        strcpy(bt[btype].Name, atom_name[atom[i].name]);
        bt[btype].Number = 0;
        bt[btype].Charge = atom[i].charge;
        bt[btype].Mass = atom[i].mass;
        bt[btype].Radius = atom[i].radius;
      } else { // check possible missing charge/mass/radius in an exsting bead type
        if (bt[btype].Charge == CHARGE) {
          bt[btype].Charge = atom[i].charge;
        }
        if (bt[btype].Mass == MASS) {
          bt[btype].Mass = atom[i].mass;
        }
        if (bt[btype].Radius == RADIUS) {
          bt[btype].Radius = atom[i].radius;
        }
      }
      // count number of beads of given type (except for 'default' beads that are counted later)
      if (atom[i].index != (count_atom_lines-1)) {
        bt[btype].Number++;
      }
    } //}}}
    // count number of beads of default type if atom default present //{{{
    if (atom[count_atom_lines-1].index == -1) {
      bt[0].Number = (*Counts).BeadsInVsf;
      for (int i = 1; i < (*Counts).TypesOfBeads; i++) {
        bt[0].Number -= bt[i].Number;
      }
    } //}}}
    // decide what bead types to merge & and do the merge //{{{
    /*
     * i) lines 'atom <id> name x mass 1'
     *          'atom <id> name x'
     *          'atom <id> name x radius 1'
     *    are to be merged as they only differ in that one has undefined mass
     * ii) lines 'atom <id> name x mass 1'
     *           'atom <id> name x'
     *           'atom <id> name x mass 2 radius 1'
     *     are not to be merged, as it's uncertain if the second atom should have mass 1 or 2
     */
    // 1) assume all bead types are to be merged //{{{
    /* helper variable merge[i][j] matrix
     * i) if i!=j; 0 means do not merge i and j, 1 means merge i and j
     * ii) if i=j; 0 means i was merged with a type x>i; 1 means there's x>i type that i was merged with
     */
    bool merge[(*Counts).TypesOfBeads][(*Counts).TypesOfBeads];
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        merge[i][j] = true;
      }
    } //}}}
    // 2) bead types that definitely aren't to be merged: //{{{
    //    i) a different Name
    //    ii) different well defined values of charge/radius/mass
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[i].Name, bt[j].Name) != 0) { // i)
          merge[i][j] = false;
        } else if ((bt[i].Charge != bt[j].Charge && bt[i].Charge != CHARGE && bt[j].Charge != CHARGE) || //
                   (bt[i].Radius != bt[j].Radius && bt[i].Radius != RADIUS && bt[j].Radius != RADIUS) || // ii)
                   (bt[i].Mass   != bt[j].Mass   && bt[i].Mass   != MASS   && bt[j].Mass   != MASS  )) { //
          merge[i][j] = false;
        }
      }
    } //}}}
/*
// test print merge matrix //{{{
printf("   ");
for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
  printf("%s ", bt[i].Name);
}
putchar('\n');
for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
  printf("%s ", bt[i].Name);
  for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
    if (j >= i) {
      printf(" %d", merge[i][j]);
    } else {
      printf("  ");
    }
  }
  putchar('\n');
} //}}}
*/
    // 3) compare bead type pairs with a third bead type //{{{
    /*
     * i) if the first two are to be merged, check that they're both also to
     *    be merged with the third one (if the third one has the same name);
     *    otherwise, nothing is to be merged
     * ii) if the first two aren't to be merged, nohing of the same name is
     *     to be merged
     */
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        if (merge[i][j]) { // i)
          for (int k = (j+1); k < (*Counts).TypesOfBeads; k++) {
            if (strcmp(bt[i].Name, bt[k].Name) == 0 &&
                (!merge[i][k] || !merge[j][k])) {
              merge[i][j] = false;
              merge[i][k] = false;
              merge[j][k] = false;
            }
          }
        } else if (strcmp(bt[i].Name, bt[j].Name) == 0) { // ii)
          for (int k = (j+1); k < (*Counts).TypesOfBeads; k++) {
            if (strcmp(bt[i].Name, bt[k].Name) == 0) {
              merge[i][j] = false;
              merge[i][k] = false;
              merge[j][k] = false;
            }
          }
        }
      }
    } //}}}
    // 4) fill diagonal merge matrix elements and make it symmetric //{{{
    /*
     * They have 0 if any non-diagonal elements have 1 (i.e., when that type
     * is to be merged) and 1 if all non-diagonal elements have 0
     */
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        if (merge[i][j]) {
          merge[i][i] = false;
          break;
        }
      }
    }
    // symmetricise
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        merge[j][i] = merge[i][j];
      }
    } //}}}
/*
// test print merge matrix //{{{
printf("   ");
for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
  printf("%s ", bt[i].Name);
}
putchar('\n');
for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
  printf("%s ", bt[i].Name);
  for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
    printf(" %d", merge[i][j]);
  }
  putchar('\n');
} //}}}
PrintBeadType2((*Counts).TypesOfBeads, bt);
*/
    // 5) merge the bead types to be merged //{{{
    /*
     * Copy values from a temporary BEADTYPE struct if a diagonal merge
     * element have 1; then, if there are any 0 diagonal elements for a bead
     * type sharing its name with the original type, copy values from the
     * first 0 diagonal element, as that contains the merged values.
     */
    // copy all bead types to a temporary struct //{{{
    BEADTYPE *temp = calloc((*Counts).TypesOfBeads, sizeof(BEADTYPE));
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      strcpy(temp[i].Name, bt[i].Name);
      temp[i].Number = bt[i].Number;
      temp[i].Charge = bt[i].Charge;
      temp[i].Mass = bt[i].Mass;
      temp[i].Radius = bt[i].Radius;
      if (!merge[i][i]) {
        for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
          if (merge[i][j]) {
            temp[i].Number += bt[j].Number;
            if (temp[i].Charge == CHARGE) {
              temp[i].Charge = bt[j].Charge;
            }
            if (temp[i].Mass == MASS) {
              temp[i].Mass = bt[j].Mass;
            }
            if (temp[i].Radius == RADIUS) {
              temp[i].Radius = bt[j].Radius;
            }
          }
        }
      }
    } //}}}
    int count = 0;
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      if (merge[i][i]) {
        strcpy(bt[count].Name, temp[i].Name);
        bt[count].Number = temp[i].Number;
        bt[count].Charge = temp[i].Charge;
        bt[count].Mass = temp[i].Mass;
        bt[count].Radius = temp[i].Radius;
        for (int j = 0; j < i; j++) {
          if (!merge[j][j] && merge[i][j]) {
            bt[count].Number = temp[j].Number;
            bt[count].Charge = temp[j].Charge;
            bt[count].Mass = temp[j].Mass;
            bt[count].Radius = temp[j].Radius;
            break;
          }
        }
        count++;
      }
    }
    (*Counts).TypesOfBeads = count; //}}}
    //}}}
    // reorder the bead types so that those with the same name are next to each other //{{{
    count = 0;
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      // copy all bead types to a temporary struct
      strcpy(temp[i].Name, bt[i].Name);
      temp[i].Number = bt[i].Number;
      temp[i].Charge = bt[i].Charge;
      temp[i].Mass = bt[i].Mass;
      temp[i].Radius = bt[i].Radius;
      temp[i].Use = false; // to ascertain whether this was already copied back to bt
    }
    // copy the bead types back in a proper order
    count = 0;
    for (int i = 0; i < ((*Counts).TypesOfBeads-1); i++) {
      if (!temp[i].Use) {
        strcpy(bt[count].Name, temp[i].Name);
        bt[count].Number = temp[i].Number;
        bt[count].Charge = temp[i].Charge;
        bt[count].Mass = temp[i].Mass;
        bt[count].Radius = temp[i].Radius;
        temp[i].Use = true;
        for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
          if (strcmp(temp[i].Name, temp[j].Name) == 0 && !temp[j].Use) {
            count++;
            strcpy(bt[count].Name, temp[j].Name);
            bt[count].Number = temp[j].Number;
            bt[count].Charge = temp[j].Charge;
            bt[count].Mass = temp[j].Mass;
            bt[count].Radius = temp[j].Radius;
            temp[j].Use = true;
          }
        }
        count++;
      }
    } //}}}
    free(temp);
    //}}}
  } else { // check only bead type name //{{{
    // charge/mass/radius is taken from the first bead with the given name (even if undefined)
    // go through all atoms (besides atom default line)
    for (int i = 0; i < count_atoms; i++) {
      int btype = -1;
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[j].Name, atom_name[atom[i].name]) == 0) {
          btype = j;
        }
      }
      if (btype == -1) { // new bead type
        btype = (*Counts).TypesOfBeads++;
        bt = realloc(bt, (*Counts).TypesOfBeads*sizeof(BEADTYPE));
        strcpy(bt[btype].Name, atom_name[atom[i].name]);
        bt[btype].Number = 0;
        bt[btype].Charge = atom[i].charge;
        bt[btype].Mass = atom[i].mass;
        bt[btype].Radius = atom[i].radius;
      }
      // count number of beads of given type (except for 'default' beads that are counted later)
      if (atom[i].index != (count_atom_lines-1)) {
        bt[btype].Number++;
      }
    }
    // count number of beads of default type if atom default present
    if (atom[count_atom_lines-1].index == -1) {
      bt[0].Number = (*Counts).BeadsInVsf;
      for (int i = 1; i < (*Counts).TypesOfBeads; i++) {
        bt[0].Number -= bt[i].Number;
      }
    }
  } //}}}
  //}}}
  // 4) fill BEAD struct //{{{
  BEAD *bead_all = calloc((*Counts).BeadsInVsf, sizeof(BEAD));
  int *index_all = calloc((*Counts).BeadsInVsf, sizeof(int)); // connect to internal indexing (atom_id)
  // unbonded beads
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int internal_id = atom_id[i];
    if (atom[internal_id].resid == -1) {
      int id = (*Counts).Unbonded++;
      bead_all[id].Molecule = -1;
      index_all[id] = internal_id;
      if (detailed) { // find bead type based on all information //{{{
        /*
         * i) find a bead type that shares a name
         * ii) check whether there are other bead types with the same name
         *     a) if so, go through all bead types and check name as well as
         *        charge, mass, and radius
         *     b) if not, no need to go through the bead types again
         *     ...necessary as when there's only one bead type with given
         *     name, the test bead may have unspecified value of
         *     charge/mass/radius, whereas the bead type may have the value
         *     specified
         */
        for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
          if (strcmp(bt[j].Name, atom_name[atom[internal_id].name]) == 0) {
            bool name2 = false;
            for (int k = 0; k < (*Counts).TypesOfBeads; k++) {
              if (j != k && strcmp(bt[j].Name, bt[k].Name) == 0) {
                name2 = true;
                break;
              }
            }
            if (name2) {
              for (int k = 0; k < (*Counts).TypesOfBeads; k++) {
                if (strcmp(bt[k].Name, atom_name[atom[internal_id].name]) == 0 &&
                    bt[k].Charge == atom[internal_id].charge &&
                    bt[k].Mass == atom[internal_id].mass &&
                    bt[k].Radius == atom[internal_id].radius) {
                  bead_all[id].Type = k;
                }
              }
            } else {
              bead_all[id].Type = j;
            }
          }
        }
        //}}}
      } else { // find bead based only on name //{{{
        bead_all[id].Type = FindBeadType2(atom_name[atom[internal_id].name], (*Counts).TypesOfBeads, bt);
      } //}}}
      if (internal_id == (count_atom_lines-1)) { // default bead
        bead_all[id].Index = -1; // to be filled later
      } else {
        bead_all[id].Index = atom[internal_id].index;
      }
    }
  }
  // bonded beads
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int internal_id = atom_id[i];
    if (atom[internal_id].resid != -1) {
      int id = (*Counts).Unbonded + (*Counts).Bonded++;
      bead_all[id].Molecule = atom[internal_id].resid;
      if (detailed) { // find bead type based on all information //{{{
        /*
         * Same as for unbonded
         * i) find a bead type that shares a name
         * ii) check whether there are other bead types with the same name
         *     a) if so, go through all bead types and check name as well as
         *        charge, mass, and radius
         *     b) if not, no need to go through the bead types again
         *     ...necessary as when there's only one bead type with given
         *     name, the test bead may have unspecified value of
         *     charge/mass/radius, whereas the bead type may have the value
         *     specified
         */
        for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
          if (strcmp(bt[j].Name, atom_name[atom[internal_id].name]) == 0) {
            bool name2 = false;
            for (int k = 0; k < (*Counts).TypesOfBeads; k++) {
              if (j != k && strcmp(bt[j].Name, bt[k].Name) == 0) {
                name2 = true;
                break;
              }
            }
            if (name2) {
              for (int k = 0; k < (*Counts).TypesOfBeads; k++) {
                if (strcmp(bt[k].Name, atom_name[atom[internal_id].name]) == 0 &&
                    bt[k].Charge == atom[internal_id].charge &&
                    bt[k].Mass == atom[internal_id].mass &&
                    bt[k].Radius == atom[internal_id].radius) {
                  bead_all[id].Type = k;
                }
              }
            } else {
              bead_all[id].Type = j;
            }
          }
        }
        //}}}
      } else { // find bead based only on name //{{{
        bead_all[id].Type = FindBeadType2(atom_name[atom[internal_id].name], (*Counts).TypesOfBeads, bt);
      } //}}}
      bead_all[id].Index = atom[internal_id].index;
      index_all[id] = internal_id;
    }
  }
  // fill BEAD[].Index array for default beads //{{{
  // i) identify already used indices (i.e., indices specifically written in the vsf file)
  bool *filled = calloc((*Counts).BeadsInVsf, sizeof(bool));
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if (bead_all[i].Index != -1) {
      filled[bead_all[i].Index] = true;
    }
  }
  // ii) assign unused indices to beads of default type
  int count = 0;
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if (bead_all[i].Index == -1) {
      for (; count < (*Counts).BeadsInVsf; count++) {
        if (!filled[count]) {
          bead_all[i].Index = count;
          count++;
          break;
        }
      }
    }
  }
  free(filled); //}}}
  // construct 'proper' Index array //{{{
  /*
   * until now, index_all was connected to atom struct, but from now on it will
   * be connected to the BEAD bead_all struct, because atom struct is no longer
   * used anywhere
   */
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    index_all[bead_all[i].Index] = i;
  } //}}}
  // rename the bead types with the same name (does nothing if !detailed) //{{{
  // ...must be after 4) as that uses the original, possibly duplicate names
  for (int i = 0; i < ((*Counts).TypesOfBeads-1); i++) {
    count = 0;
    for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
      if (strcmp(bt[i].Name, bt[j].Name) == 0) {
        count++;
        char name[17];
        // maximum strlen is 16, so if 15th or 14th place in Name (i.e.,
        // Name[14] or Name[13]) must contain '\0' (if the name is short, it
        // changes nothing)
        if (count < 10) {
          bt[j].Name[14] = '\0';
        } else if (count < 100) {
          bt[j].Name[13] = '\0';
        }
        strcpy(name, bt[j].Name);
        sprintf(bt[j].Name, "%s_%d", name, count);
      }
    }
  } //}}}
  //}}}
  // 5) identify molecule types base on all data //{{{
  /*
   * Molecules of one type must share:
   * i) molecule name
   * ii) number of beads and bonds
   * iii) order of bead types
   * iv) connectivity
   */
  // count beads in each molecule for ii) //{{{
  int *atoms_per_mol = calloc ((*Counts).Molecules, sizeof(int));
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int internal_id = atom_id[i];
    if (atom[internal_id].resid != -1) {
      atoms_per_mol[atom[internal_id].resid]++;
    }
  } //}}}
// test print number of beads in each molecule //{{{
/*
for (int i = 0; i < (*Counts).Molecules; i++) {
  printf("%d: %d\n", i, atoms_per_mol[i]);
}
putchar('\n');
*/ //}}}
  // allocate MOLECULE struct //{{{
  MOLECULE *mol = calloc((*Counts).Molecules, sizeof(MOLECULE));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    mol[i].Bead = calloc(atoms_per_mol[i], sizeof(int));
//printf("%d\n", atoms_per_mol[i]);
  } //}}}
  // add bead indices to Molecule[].Bead array for iv) //{{{
  int *count_in_mol = calloc((*Counts).Molecules, sizeof(int));
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int bead_id = atom_id[i];
    int mol_id = atom[bead_id].resid;
    if (mol_id != -1) {
      mol[mol_id].Bead[count_in_mol[mol_id]] = index_all[atom[bead_id].index];
      count_in_mol[mol_id]++;
    }
  } //}}}
// test print beads in molecules //{{{
/*
for (int i = 0; i < (*Counts).Molecules; i++) {
  printf("%d:", i);
  for (int j = 0; j < atoms_per_mol[i]; j++) {
    printf(" %d", mol[i].Bead[j]);
  }
  putchar('\n');
}
putchar('\n');
*/ //}}}
  // error - resid numbering in vsf must be continuous //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    if (atoms_per_mol[i] == 0) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", vsf_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - 'resid' numbering is not continuous\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  } //}}}
  // count bonds in all molecules for ii) //{{{
  int *bonds_per_mol = calloc ((*Counts).Molecules, sizeof(int));
  for (int i = 0; i < count_bond_lines; i++) {
    int mol_id1 = atom[atom_id[bond[i].index1]].resid;
    int mol_id2 = atom[atom_id[bond[i].index2]].resid;
    // error - beads from one bond are in different molecules
    if (mol_id1 != mol_id2) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", vsf_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - bonded beads in different molecules (beads ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", bond[i].index1);
      RedText(STDERR_FILENO);
      fprintf(stderr, " and ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", bond[i].index2);
      RedText(STDERR_FILENO);
      fprintf(stderr, ")\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
    bonds_per_mol[mol_id1]++;
  } //}}}
  // save connectivity for each molecule for iv) //{{{
  // allocate 3D connectivity array //{{{
  /*
   * in connectivity[i][j][k]:
   *   i ... molecule id
   *   j ... bond id
   *   k ... 0 and 1 contain connected bead indices; 2 is not used as there are no bond types
   */
  int ***connectivity = calloc((*Counts).Molecules, sizeof(int **));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    count_in_mol[i] = 0;
    connectivity[i] = calloc(bonds_per_mol[i], sizeof(int *));
    for (int j = 0; j < bonds_per_mol[i]; j++) {
      connectivity[i][j] = calloc(3, sizeof(int));
      connectivity[i][j][0] = -1;
      connectivity[i][j][1] = -1;
    }
  } //}}}
  // save ids of bonded beads
  for (int i = 0; i < count_bond_lines; i++) {
    int mol_id = atom[atom_id[bond[i].index1]].resid;
    int id1 = index_all[bond[i].index1];
    int id2 = index_all[bond[i].index2];
    bool done[2] = {false};
    for (int j = 0; j < atoms_per_mol[mol_id]; j++) {
      if (!done[0] && mol[mol_id].Bead[j] == id1) {
        connectivity[mol_id][count_in_mol[mol_id]][0] = j;
      }
      if (!done[1] && mol[mol_id].Bead[j] == id2) {
        connectivity[mol_id][count_in_mol[mol_id]][1] = j;
      }
      if (done[0] && done[1]) {
        break;
      }
    }
    count_in_mol[mol_id]++;
  }
  // sort bonds -- pro forma as bead ids in atom struct are sorted
  for (int i = 0; i < (*Counts).Molecules; i++) {
    SortBonds(connectivity[i], bonds_per_mol[i]);
  }
  //}}}
// test print bonds in molecules //{{{
/*
for (int i = 0; i < (*Counts).Molecules; i++) {
  printf("%d:", i);
  for (int j = 0; j < bonds_per_mol[i]; j++) {
    printf(" %d-%d", connectivity[i][j][0]+1, connectivity[i][j][1]+1);
  }
  putchar('\n');
}
putchar('\n');
*/ //}}}
  // count molecule types based on i) through to iv) //{{{
  // ...or just on i) if !detailed (and exit if a molecule contains too few/too many beads)
  MOLECULETYPE *mt = calloc(1, sizeof(MOLECULETYPE));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    bool add = false;
    int first_bead = bead_all[mol[i].Bead[0]].Index;
    int name = atom[atom_id[first_bead]].resname;
//  printf("%5d %5d %5d %s\n", first_bead, atom_id[first_bead], name, res_name[name]);
    // go through all molecule types to find a match for molecule i //{{{
    for (int j = 0; j < (*Counts).TypesOfMolecules; j++) {
      if (detailed) { // check name, numbers of beads & bonds, and connectivity
        if (strcmp(res_name[name], mt[j].Name) == 0 &&
            atoms_per_mol[i] == mt[j].nBeads &&
            bonds_per_mol[i] == mt[j].nBonds) {
          bool same_beads = true;
          for (int k = 0; k < mt[j].nBeads; k++) {
            if (bead_all[mol[i].Bead[k]].Type != mt[j].Bead[k]) {
              same_beads = false;
              break;
            }
          }
          bool same_bonds = true;
          for (int k = 0; k < mt[j].nBonds; k++) {
            if (connectivity[i][k][0] != mt[j].Bond[k][0] ||
                connectivity[i][k][1] != mt[j].Bond[k][1]) {
              same_bonds = false;
              break;
            }
          }
          if (same_beads && same_bonds) {
            add = true;
          }
        }
      } else { // check name only
        if (strcmp(res_name[name], mt[j].Name) == 0) {
          add = true;
          // error - too few/too many beads in a molecule
          if (atoms_per_mol[i] != mt[j].nBeads) {
            RedText(STDERR_FILENO);
            fprintf(stderr, "\nError: ");
            YellowText(STDERR_FILENO);
            fprintf(stderr, "%s", vsf_file);
            RedText(STDERR_FILENO);
            fprintf(stderr, " - molecule ");
            YellowText(STDERR_FILENO);
            fprintf(stderr, "resid %d", i);
            RedText(STDERR_FILENO);
            fprintf(stderr, " contains too ");
            if (atoms_per_mol[i] > mt[j].nBeads) {
              fprintf(stderr, "many");
            } else {
              fprintf(stderr, "few");
            }
            fprintf(stderr, "beads (%d instead of %d)\n", atoms_per_mol[i], mt[j].nBeads);
            ResetColour(STDERR_FILENO);
            exit(1);
          }
        }
      }
      if (add) { // found matching molecule type
        mt[j].Number++;
        mol[i].Type = j;
        break;
      }
    } //}}}
    // didn't find matching type (or no type yet exists), so create molecule type //{{{
    if (!add || (*Counts).TypesOfMolecules == 0) {
      int type = (*Counts).TypesOfMolecules++;
      mol[i].Type = type;
      mt = realloc(mt, (*Counts).TypesOfMolecules*sizeof(MOLECULETYPE));
      strcpy(mt[type].Name, res_name[name]);
      mt[type].Number = 1;
      mt[type].nBeads = atoms_per_mol[i];
      mt[type].Bead = calloc(mt[type].nBeads, sizeof(int));
      for (int j = 0; j < mt[type].nBeads; j++) {
        int bead = mol[i].Bead[j];
        int id = bead;
        mt[type].Bead[j] = bead_all[id].Type;
      }
      mt[type].nBonds = bonds_per_mol[i];
      if (bonds_per_mol[i] > 0) {
        mt[type].Bond = calloc(bonds_per_mol[i], sizeof(int *));
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          mt[type].Bond[j] = calloc(3, sizeof(int));
          mt[type].Bond[j][0] = connectivity[i][j][0];
          mt[type].Bond[j][1] = connectivity[i][j][1];
          mt[type].Bond[j][2] = -1; // no bond types
        }
      }
      mt[type].nAngles = 0;
    } //}}}
  } //}}}
  // rename the molecule types with the same name (does nothing if !detailed) //{{{
  for (int i = 0; i < ((*Counts).TypesOfMolecules-1); i++) {
    count = 0;
    for (int j = (i+1); j < (*Counts).TypesOfMolecules; j++) {
      if (strcmp(mt[i].Name, mt[j].Name) == 0) {
        count++;
        char name[17];
        // maximum strlen is 16, so if 15th or 14th place in Name (i.e.,
        // Name[14] or Name[13]) must contain '\0' (if the name is short, it
        // changes nothing)
        if (count < 10) {
          mt[j].Name[14] = '\0';
        } else if (count < 100) {
          mt[j].Name[13] = '\0';
        }
        strcpy(name, mt[j].Name);
        sprintf(mt[j].Name, "%s_%d", name, count);
      }
    }
  } //}}}
  // calculate molecules' mass and charge and fill their BType array
  FillMolType((*Counts).TypesOfMolecules, bt, &mt); //}}}
  // 6) end by copying everything back to their 'proper' arrays and structures //{{{
  (*Counts).Beads = (*Counts).BeadsInVsf;
  // indexing
  *Index = calloc((*Counts).Beads, sizeof(int));
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Index)[i] = index_all[i];
  }
  // bead types
  *BeadType = calloc((*Counts).TypesOfBeads, sizeof(BEADTYPE));
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    strcpy((*BeadType)[i].Name, bt[i].Name);
    (*BeadType)[i].Number = bt[i].Number;
    (*BeadType)[i].Charge = bt[i].Charge;
    (*BeadType)[i].Mass = bt[i].Mass;
    (*BeadType)[i].Radius = bt[i].Radius;
  }
  // individual beads
  *Bead = calloc((*Counts).Beads, sizeof(BEAD));
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Type = bead_all[i].Type;
    (*Bead)[i].Molecule = bead_all[i].Molecule;
    (*Bead)[i].nAggregates = 0;
    (*Bead)[i].Index = bead_all[i].Index;
  }
  // molecule types
  *MoleculeType = calloc((*Counts).TypesOfMolecules, sizeof(MOLECULETYPE));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    strcpy((*MoleculeType)[i].Name, mt[i].Name);
    (*MoleculeType)[i].Number = mt[i].Number;
    (*MoleculeType)[i].nBeads = mt[i].nBeads;
    (*MoleculeType)[i].Bead = calloc((*MoleculeType)[i].nBeads, sizeof(int));
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      (*MoleculeType)[i].Bead[j] = mt[i].Bead[j];
    }
    (*MoleculeType)[i].nBonds = mt[i].nBonds;
    if ((*MoleculeType)[i].nBonds > 0) {
      (*MoleculeType)[i].Bond = calloc((*MoleculeType)[i].nBonds, sizeof(int *));
      for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
        (*MoleculeType)[i].Bond[j] = calloc(3, sizeof(int));
        (*MoleculeType)[i].Bond[j][0] = mt[i].Bond[j][0];
        (*MoleculeType)[i].Bond[j][1] = mt[i].Bond[j][1];
        (*MoleculeType)[i].Bond[j][2] = mt[i].Bond[j][2];
      }
    }
    (*MoleculeType)[i].nAngles = 0;
//  (*MoleculeType)[i].Angle = calloc(1, sizeof(int *)); // ...why should this be allocated when there's no angles?
    (*MoleculeType)[i].nBTypes = mt[i].nBTypes;
    (*MoleculeType)[i].BType = calloc((*MoleculeType)[i].nBTypes, sizeof(int));
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      (*MoleculeType)[i].BType[j] = mt[i].BType[j];
    }
    (*MoleculeType)[i].Mass = mt[i].Mass;
    (*MoleculeType)[i].Charge = mt[i].Charge;
  }
  // individual molecules
  *Molecule = calloc((*Counts).Molecules, sizeof(MOLECULE));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Molecule)[i].Type = mol[i].Type;
    (*Molecule)[i].Bead = calloc((*MoleculeType)[(*Molecule)[i].Type].nBeads, sizeof(int));
    for (int j = 0; j < (*MoleculeType)[(*Molecule)[i].Type].nBeads; j++) {
      (*Molecule)[i].Bead[j] = mol[i].Bead[j];
    }
    (*Molecule)[i].Aggregate = -1;
  } //}}}
  // free memory //{{{
  free(bt);
  free(index_all);
  FreeBead2((*Counts).BeadsInVsf, &bead_all);
  FreeMolecule2((*Counts).Molecules, &mol);
  FreeMoleculeType2((*Counts).TypesOfMolecules, &mt);
  for (int i = 0; i < (*Counts).Molecules; i++) {
    for (int j = 0; j < bonds_per_mol[i]; j++) {
      free(connectivity[i][j]);
    }
    free(connectivity[i]);
  }
  free(connectivity);
  for (int i = 0; i < atom_names; i++) {
    free(atom_name[i]);
  }
  free(atom_name);
  free(atom_id);
  for (int i = 0; i < res_names; i++) {
    free(res_name[i]);
  }
  free(res_name);
  free(atom);
  free(bond);
  free(atoms_per_mol);
  free(bonds_per_mol);
  free(count_in_mol); //}}}
//PrintMoleculeType2((*Counts).TypesOfMolecules, *BeadType, *MoleculeType);
//PrintMolecule((*Counts).Molecules, *MoleculeType, *Molecule, *BeadType, *Bead);
//PrintBead2((*Counts).Beads, *Index, *BeadType, *Bead);
} //}}}

void FullVtfRead(char *vsf_file, char *vcf_file, bool detailed, bool vtf, bool *indexed, int *struct_lines,
                 VECTOR *BoxLength, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {
  // get box size
  *BoxLength = GetPBC2(vcf_file);
  // read the whole structure section
  ReadVtfStructure(vsf_file, detailed, Counts, BeadType, Bead, Index, MoleculeType, Molecule);
  // count structure lines if vtf as the coordinate file (-1 otherwise)
  *struct_lines = CountVtfStructLines(vtf, vcf_file);
  // determine timestep type & what coordinate file contains from the first timestep
  FILE *vcf;
  if ((vcf = fopen(vcf_file, "r")) == NULL) {
    ErrorFileOpen(vcf_file, 'r');
    exit(1);
  }
  SkipVtfStructure(vtf, vcf, *struct_lines);
  *indexed = CheckVtfTimestep(vcf, vcf_file, Counts, BeadType, Bead, Index, MoleculeType, Molecule);
  fclose(vcf);

  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, vsf_file);

  // allocate memory for aggregates
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(1, sizeof(int));
  }
}

// ReadStructure() //{{{
/** Function reading information about beads and molecules from DL_MESO
 * `FIELD` file, `.vsf` structure file, and `.vcf` coordinate file. Charge
 * and mass of beads is read from `FIELD` file, but all other information
 * is read from `vsf` structure file. If `vcf` coordinate is passed to
 * `ReadStructure()`, bead types not present in the `vcf` file get ignored.
 *
 * Overly complicated, some things done more than once, needs complete
 * overhaul. Thankfully, it works.
 */
bool ReadStructure(char *vsf_file, char *vcf_file, COUNTS *Counts,
                   BEADTYPE **BeadType, BEAD **Bead, int **Index,
                   MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {

  FILE *vsf;

  (*Counts) = InitCounts; // zeroize

  // initial allocations (realloced later)
  *BeadType = calloc(1, sizeof(struct BeadType));
  *MoleculeType = calloc(1, sizeof(struct MoleculeType));

  // delimiters for SplitLine - whitespace in general, whitespace+colon in bond lines
  char delim[8], bond_delim[8];
  strcpy(delim, " \t");
  strcpy(bond_delim, " \t:");

  // first, read through vsf to find highest bead/mol id and all bead/mol types and identify default bead type (if exist) //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}
  int max_bead = 0, // highest bead id detected
      max_mol = 0, // highest mol id detected
      total_atom_lines = 0, // number of atom lines in vsf (including comments and blanks)
      atom_lines = 0, // number of atom lines in vsf (excluding comments and blanks)
      type_default = -1, // default bead type; -1 remains if no default exists
      missing = 1000; // placeholder for testing whether bead charge and/or mass is present
  char line[LINE];
  while(fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    int words = SplitLine(split, line, delim);
    if (strncmp(split[0], "atom", 1) == 0) { // go through an atom line
      total_atom_lines++;
      atom_lines++;
      // error - odd number of values: a(tom) lines are composed of 'keyword <value>' pairs //{{{
      if ((words % 2) == 1) { // because it's split[0] ... split[2n-1]
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - odd number of strings an atom line\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // line must begin with 'a(tom) default' or 'a(tom) <int>' //{{{
      bool def = false;
      if (strcmp("default", split[1]) == 0) { // if 'default' keyword, assign type_default
        def = true; // default line - do not count this atom to the total number of atoms
        // warn if there has already been a default line
        if (type_default != -1) {
          fprintf(stderr, "\033[1;33m");
          fprintf(stderr, "\nWarning: \033[1;36m%s\033[1;33m - multiple 'atom default' lines\n", vsf_file);
          fprintf(stderr, "\033[0m");
        }
        // scan the line for 'name' to identify the default bead
        bool name = false;
        for (int j = 2; j < words; j += 2) {
          if (strncmp(split[j], "name", 1) == 0) {
            type_default = FindBeadType(split[j+1], *Counts, *BeadType);
            if (type_default == -1) {
              type_default = (*Counts).TypesOfBeads;
            }
            name = true;
          }
        }
        // error: missing name
        if (!name) {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom line must contain 'name'\n", vsf_file);
          fprintf(stderr, "\033[0m");
          ErrorPrintLine(split, words);
          exit(1);
        }
      } else if (!IsInteger(split[1])) { // if not 'default' line, there must be integer after 'atom'
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'a(tom)' must be followed by 'default' or a whole number\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } else if (atoi(split[1]) > max_bead) { // if there's atom id, test if it's the highest id yet
        max_bead = atoi(split[1]);
      } //}}}
      // check that line contains either both 'resid' and 'resname' or neither //{{{
      // TODO: only resid is necessary, if not resname, add some autogenerated name, creating a new MoleculeType
      int test = 0;
      for (int i = 2; i < words && test < 2; i+= 2) {
        if (strcmp("resname", split[i]) == 0) {
          test++;
        } else if (strcmp("resid", split[i]) == 0) {
          test++;
          if (!IsInteger(split[i+1]) || atoi(split[i+1]) == 0) {
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'resid' must be followed by a positive whole number\n", vsf_file);
            fprintf(stderr, "\033[0m");
            ErrorPrintLine(split, words);
            exit(1);
          }
        }
      }
      if (test == 1) { // test == 0: unbonded bead; test == 1: only one keyword present; test == 2: both keywords present
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom line must contain both 'resid' and 'resname'\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // go through the rest of the line by <keyword> <value> pairs
      int new = -1, btype; // test whether this is a new bead type
      double charge = missing, mass = missing; // just a placeholder to know if bead types have charge and/or mass
      bool name = false; // is the necessary name keyword present?
      for (int i = 2; i < words; i += 2) {
        if (strncmp("name", split[i], 1) == 0) { // bead type name //{{{
          name = true;
          // if the bead name doesn't exist, add it to the structures
          btype = FindBeadType(split[i+1], *Counts, *BeadType);
          if (btype == -1) {
            // increment number of bead types
            btype = (*Counts).TypesOfBeads++;
            // realloc BeadType array
            *BeadType = realloc(*BeadType, (*Counts).TypesOfBeads*sizeof(struct BeadType));
            (*BeadType)[btype].Number = 0;
            // TODO: make a safe copy function
            // copy new name to BeadType[].Name
            if (strlen(split[i+1]) > (sizeof (*BeadType)[btype].Name)) {
              split[i+1][sizeof (*BeadType)[btype].Name - 1] = '\0';
            }
            strcpy((*BeadType)[btype].Name, split[i+1]);
            // initialize charge and mass by 'impossible' value
            (*BeadType)[btype].Charge = missing;
            (*BeadType)[btype].Mass = missing;
            new = btype;
          }
          // increment number of beads of given type (except if it's a 'default' line)
          if (!def) {
            (*BeadType)[btype].Number++;
          } //}}}
        } else if (strcmp("resname", split[i]) == 0) { // molecule type name //{{{
          // if the molecule name doesn't exist, add it to the structures
          int mtype = FindMoleculeType(split[i+1], *Counts, *MoleculeType);
          if (mtype == -1) {
            // increment number of molecule types
            mtype = (*Counts).TypesOfMolecules++;
            // realloc MoleculeType array
            *MoleculeType = realloc(*MoleculeType, (*Counts).TypesOfMolecules*sizeof(struct MoleculeType));
            // initialize stuff
            (*MoleculeType)[mtype].Number = 0;
            (*MoleculeType)[mtype].nBonds = 0;
            (*MoleculeType)[mtype].nBeads = 0;
            (*MoleculeType)[mtype].nBTypes = 0;
            (*MoleculeType)[mtype].Use = false;
            // copy new name to MoleculeType[].Name
            strcpy((*MoleculeType)[mtype].Name, split[i+1]);
          } //}}}
        } else if (strcmp("resid", split[i]) == 0) { // molecule id //{{{
          if (atoi(split[i+1]) > max_mol) {
            if (!IsInteger(split[i+1]) || atoi(split[i+1]) == 0) {
              fprintf(stderr, "\033[1;31m");
              fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'resid' must be followed by a positive whole number\n", vsf_file);
              fprintf(stderr, "\033[0m");
              ErrorPrintLine(split, words);
              exit(1);
            }
            max_mol = atoi(split[i+1]);
          }
        //}}}
        } else if (strcmp("charge", split[i]) == 0 || strcmp("q", split[i]) == 0) { //{{{
          if (!IsDouble(split[i+1])) {
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom charge must be a real number\n", vsf_file);
            fprintf(stderr, "\033[0m");
            ErrorPrintLine(split, words);
            exit(1);
          }
          charge = atof(split[i+1]);
        //}}}
        } else if (strncmp("mass", split[i], 1) == 0) { //{{{
          if (!IsPosDouble(split[i+1])) {
            fprintf(stderr, "\033[1;31m");
            fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom mass must be a non-negative real number\n", vsf_file);
            fprintf(stderr, "\033[0m");
            ErrorPrintLine(split, words);
            exit(1);
          }
          mass = atof(split[i+1]);
        } //}}}
      }
      // error - no 'name' //{{{
      if (!name) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - atom line must contain 'name'\n", vsf_file);
        fprintf(stderr, "\033[0m");
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // assign charge/mass only if it's the first bead of that type or if it doesn't have proper charge/mass yet
      if (new != -1 || (*BeadType)[btype].Charge == missing) {
        (*BeadType)[btype].Charge = charge;
      }
      if (new != -1 || (*BeadType)[btype].Mass == missing) {
        (*BeadType)[btype].Mass = mass;
      }
    } else if (split[0][0] == '#' || split[0][0] == '\0' || strcasecmp(split[0], "pbc") == 0) { // count the line if comment, blank, or pbc
      total_atom_lines++;
    } else { // end of vsf structure block
      break;
    }
  }
  fclose(vsf);
  // test whether there's enough atom lines if there's no 'default' line
  if (type_default == -1 && atom_lines != (max_bead+1)) {
    fprintf(stderr, "\033[1;31m");
    fprintf(stderr, "\nError - \033[1;33m%s\033[1;31m: too few atom lines (or 'default' atom missing)\n", vsf_file);
    fprintf(stderr, "\033[0m");
    exit(1);
  } //}}}

  (*Counts).BeadsInVsf = max_bead + 1; // bead ids start from 0 in vsf
  (*Counts).Molecules = max_mol; // mol ids start from 1 in vsf
  // allocate Bead and Molecule structures
  *Bead = calloc((*Counts).BeadsInVsf, sizeof(struct Bead));
  *Molecule = calloc(max_mol, sizeof(struct Molecule));
  // reverse of Bead[].Index
  *Index = malloc((*Counts).BeadsInVsf*sizeof(int));

  // assign type -1 to all beads to later identify 'default' beads
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    (*Bead)[i].Type = -1;
  }
  // assigne type -1 to all molecules to later check for possible errors
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Molecule)[i].Type = -1;
  }

  // second, read through vsf to find stuff about beads and mols //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}
  // helper array to identify which molecule is used to count beads in which molecule type
  int first_mol[(*Counts).TypesOfMolecules];
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    first_mol[i] = -1;
  }
  // go through atom lines - no error checking,
  // because it's been done in the first read-through
  for (int count = 0; count < total_atom_lines; count++) {
    fgets(line, sizeof(line), vsf);
    char split[30][100];
    int words = SplitLine(split, line, delim);
    // go through the line (ignoring the 'default' one)
    if (split[0][0] == 'a' && strcmp("default", split[1]) != 0) {
      int bead_type = -1, mol_id = -1, mol_type = -1;
      int bead_id = atoi(split[1]);
      // by <keyword> <value> pairs
      for (int i = 2; i < words; i += 2) {
        if (strncmp("name", split[i], 1) == 0) { // bead name
          bead_type = FindBeadType(split[i+1], *Counts, *BeadType);
        } else if (strcmp("resname", split[i]) == 0) { // molecule name
          mol_type = FindMoleculeType(split[i+1], *Counts, *MoleculeType);
        } else if (strcmp("resid", split[i]) == 0) { // molecule id
          mol_id = atoi(split[i+1]) - 1; // mol ids start with 1 in vsf
        }
      }
      // assign values to Bead & Molecule
      (*Index)[bead_id] = bead_id; // reverse of Bead[].Index
      (*Bead)[bead_id].Type = bead_type;
      (*Bead)[bead_id].Index = bead_id; // same as Index[] for now as it assumes all beads are present
      (*Bead)[bead_id].Molecule = mol_id;
      if (mol_id > -1) {
        (*Molecule)[mol_id].Type = mol_type;
        // if this molecule type wasn't used yet, use mol_id to calculate its number of beads
        if (first_mol[mol_type] == -1) {
          first_mol[mol_type] = mol_id;
        }
        if (first_mol[mol_type] == mol_id) {
          (*MoleculeType)[mol_type].nBeads++;
        }
      }
    }
  } //}}}

  // assign 'type_default' to default beads //{{{
  if (type_default != -1) {
    for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
      if ((*Bead)[i].Type == -1) {
        (*Index)[i] = i; // reverse of Bead[].Index
        (*Bead)[i].Type = type_default;
        (*Bead)[i].Index = i;
        (*Bead)[i].Molecule = -1; // default beads aren't in molecules
        (*BeadType)[type_default].Number++;
      } else if ((*Bead)[i].Type == type_default) { // default type beads explicitly specified by 'atom' line
        (*BeadType)[type_default].Number++;
      }
    }
  } //}}}

  // calculate number of molecules for each molecule type
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    (*MoleculeType)[mtype].Number++;
  }

  // third, go through the bonds section of vsf to find the number of bonds in molecules //{{{
  int bonds[(*Counts).TypesOfMolecules]; // helper array
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  }
  // error checking as this is the first time bonds are reead
  // TODO: check colon; possibly double colon and assign bonds properly
  while (fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    int words = SplitLine(split, line, bond_delim);
    if (strncmp(split[0], "bond", 1) == 0) {
      if (!IsInteger(split[1]) || !IsInteger(split[2])) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bond' must be followed by two colon-separated numbers\n", vsf_file);
        fprintf(stderr, "\033[0m");
        exit(1);
      }
      int mol_id = (*Bead)[atoi(split[1])].Molecule;
      int mol_type = (*Molecule)[mol_id].Type;
      // is this mol_type already in use in bond numbers calculation?
      if (bonds[mol_type] == -1) {
        bonds[mol_type] = mol_id;
      }
      // increment number of bonds if correct molecule
      if (bonds[mol_type] == mol_id) {
        (*MoleculeType)[mol_type].nBonds++;
      }
    /* Break while loop if at end of bond section;
     * Excepting comments or white lines, behind bonds section can be either coordinates (vtf file):
     * timestep starting line can be
     * 1) t(imestep) ...and nothing behind it
     * 2) t(imestep) o(rdered)/i(ndexed)
     * 3) o(rdered)/i(ndexed)
     * or pbc line:
     * 4) pbc <number> <number> <number>
     */
    } else if ((split[0][0] == 't' && words == 1) || // 1)
               (split[0][0] == 't' && words > 1 && // 2) 1st string: t(imestep)
                (split[1][0] == 'o'|| split[1][0] == 'i')) || // 2) 2nd string: o/i
               (split[0][0] == 'o' || split[0][0] == 'i') || // 3)
               (strcasecmp("pbc", split[0]) == 0 && words > 3 && // 4)
                IsPosDouble(split[1]) && // 4) 1st <number>
                IsPosDouble(split[2]) && // 4) 2nd <number>
                IsPosDouble(split[3]))) { // 4) 3rd <number>
      break;
    // exit with error if not empty line or comment
    } else if (split[0][0] != '\0' && // empty line
               split[0][0] != '#') { // comment
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError - \033[1;33m%s\033[1;31m: unrecognised line\n", vsf_file);
      fprintf(stderr, "\033[0m");
      ErrorPrintLine(split, words);
      exit(1);
    }
  } //}}}

  fclose(vsf); // will be later reopened to read from the beginning again

  // TODO: Is this needed? Why not alloc precise number of bonds? Maybe for later when stuff is kicked out due to vcf?
  // find the highest number of bonds
  int max_bonds = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    if ((*MoleculeType)[i].nBonds > max_bonds) {
      max_bonds = (*MoleculeType)[i].nBonds;
    }
  }
  // allocate MoleculeType[].Bond array
  int moltype_alloced = (*Counts).TypesOfMolecules;
  for (int i = 0; i < moltype_alloced; i++) {
    (*MoleculeType)[i].Bond = calloc(max_bonds, sizeof(int *));
    for (int j = 0; j < max_bonds; j++) {
      (*MoleculeType)[i].Bond[j] = calloc(3, sizeof(int));
      (*MoleculeType)[i].Bond[j][2] = -1; // no bond type assigned
    }
  }

  // fourth, go through vsf and assign bead ids to molecules //{{{
  // open vsf structure file //{{{
  if ((vsf = fopen(vsf_file, "r")) == NULL) {
    ErrorFileOpen(vsf_file, 'r');
    exit(1);
  } //}}}
  // go through atom lines - no error checking,
  // because it's been done in the first read-through
  int *beads = calloc((*Counts).BeadsInVsf, sizeof(int)); // helper array
  for (int count = 0; count < total_atom_lines; count++) {
    fgets(line, sizeof(line), vsf);
    char split[30][100];
    int words = SplitLine(split, line, delim);
    // go through the line
    if (strncmp(split[0], "atom", 1) == 0 && strcmp("default", split[1]) != 0) { // non-default a(tom) line
      int bead_id = atoi(split[1]);
      // by <keyword> <value> pairs
      int mol_id = -1, mtype = -1;
      for (int i = 2; i < words; i += 2) {
        if (strcmp("resid", split[i]) == 0) { // molecule id
          mol_id = atoi(split[i+1]) - 1; // mol ids start with 1 in vsf
        } else if (strcmp("resname", split[i]) == 0) { // molecule name (i.e., type)
          mtype = FindMoleculeType(split[i+1], *Counts, *MoleculeType);
        }
      }
      if (mol_id > -1) {
        // allocate Molecule[].Bead array if bead_id is the first bead of the molecule
        if (beads[mol_id] == 0) {
          (*Molecule)[mol_id].Bead = calloc((*MoleculeType)[mtype].nBeads, sizeof(int));
        }
        (*Molecule)[mol_id].Bead[beads[mol_id]] = bead_id;
        beads[mol_id]++;
      }
    }
  }
  free(beads);

  // sort ids in molecule according to ascending id in vsf //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    SortArray(&(*Molecule)[i].Bead, (*MoleculeType)[mtype].nBeads, 0);
  } //}}}
  //}}}

  // fill MoleculeType[].Bead array according to the first molecule with the name MoleculeType[].Name //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Bead = calloc((*MoleculeType)[i].nBeads, sizeof(int));
    (*MoleculeType)[i].Bead[0] = -1; // just to test whether this molecule type was used already
  }
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    if ((*MoleculeType)[mtype].Bead[0] == -1) {
      for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
        int id = (*Molecule)[i].Bead[j];
        int btype = (*Bead)[id].Type;
        (*MoleculeType)[mtype].Bead[j] = btype;
      }
    }
  } //}}}

  // fifth, go again through the bonds and assign ids of bonded beads //{{{
  // initialize helper arrays //{{{
  int *count_bonds;
  count_bonds = calloc((*Counts).TypesOfMolecules, sizeof(int));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    bonds[i] = -1;
  } //}}}
  // again, no error checking, as it was done earlier already
  while (fgets(line, sizeof(line), vsf)) {
    char split[30][100];
    SplitLine(split, line, bond_delim);
    if (strncmp(split[0], "bond", 1) == 0) {
      int mol_id = (*Bead)[atoi(split[1])].Molecule;
      int mol_type = (*Molecule)[mol_id].Type;
      // is this mol_type already in use?
      if (bonds[mol_type] == -1) {
        bonds[mol_type] = mol_id;
      }
      // increment number of bonds if correct molecule
      if (bonds[mol_type] == mol_id) {
        int id1 = atoi(split[1]);
        int id2 = atoi(split[2]);
        bool test_id1 = false, test_id2 = false;
        for (int i = 0; i < (*MoleculeType)[mol_type].nBeads; i++) {
          if (!test_id1 && id1 == (*Molecule)[mol_id].Bead[i]) {
            (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][0] = i;
            test_id1 = true;
          }
          if (!test_id2 && id2 == (*Molecule)[mol_id].Bead[i]) {
            (*MoleculeType)[mol_type].Bond[count_bonds[mol_type]][1] = i;
            test_id2 = true;
          }
          if (test_id1 && test_id2) {
            break;
          }
        }
        count_bonds[mol_type]++;
      }
    } else if (split[0][0] != '\0' && split[0][0] != '#') {
      break;
    }
  }
  free(count_bonds);
  fclose(vsf); //}}}

  // sort the Bond arrays //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    SortBonds((*MoleculeType)[i].Bond, (*MoleculeType)[i].nBonds);
  } //}}}

  // fill BTypes array //{{{
  // allocate & initialize arrays //{{{
  bool used[(*Counts).TypesOfMolecules];
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].BType = calloc((*Counts).TypesOfBeads, sizeof(int));
    (*MoleculeType)[i].nBTypes = 0;
    used[i] = false;
  } //}}}
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mol_type = (*Molecule)[i].Type;
    if (!used[mol_type]) { // go only through the first molecule of the given type
      for (int j = 0; j < (*MoleculeType)[mol_type].nBeads; j++) {
        int id = (*Molecule)[i].Bead[j];
        // test if bead id's type is present in BType array //{{{
        bool in_mol = false;
        for (int k = 0; k < (*MoleculeType)[mol_type].nBTypes; k++) {
          if ((*MoleculeType)[mol_type].BType[k] == (*Bead)[id].Type) {
            in_mol = true;
            break;
          }
        } //}}}
        // if bead id is of a type not yet present in BType array, add it there
        if (!in_mol) {
          (*MoleculeType)[mol_type].BType[(*MoleculeType)[mol_type].nBTypes] = (*Bead)[id].Type;
          (*MoleculeType)[mol_type].nBTypes++;
        }
      }
    }
  } //}}}

  // sixth, go through vcf to find which bead types are there //{{{
  bool indexed = true;
  if (vcf_file[0] != '\0') {
    // open vcf coordinate file //{{{
    FILE *vcf;
    if ((vcf = fopen(vcf_file, "r")) == NULL) {
      ErrorFileOpen(vcf_file, 'r');
      exit(1);
    } //}}}

    // skip initial stuff //{{{
    // skip lines including 't(imestep)' line //{{{
    char line2[LINE], *split[30], str[30];
    str[0] = '\0';
    do {
      fgets(line, sizeof(line), vcf);
      strcpy(line2, line);
      split[0] = strtok(line, " \t");

      // 't(imestep)' line
      if (split[0][0] == 't' || split[0][0] == 'T') {
        split[1] = strtok(NULL, " \t");
        // if only 't(imestep)' present, assume ordered timestep
        if (split[1] == NULL) {
          str[0] = 'o';
          break;
        // otherwise, error if 'i(ndexed)' or 'o(ordered)' not present
        } else if (split[1][0] != 'o' && split[1][0] != 'O' &&
            split[1][0] != 'i' && split[1][0] != 'I') {
          fprintf(stderr, "\033[1;31m");
          fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m", vsf_file);
          fprintf(stderr, " - unrecognised keywords '\033[1;33m%s %s\033[1;31m'", split[0], split[1]);
          fprintf(stderr, "\033[0m");
          exit(1);
        // if no error, find if 'i(ndexed)' or 'o(rdered)' timestep
        } else {
          str[0] = split[1][0];
          break;
        }

      // 'o(rdered)' or 'i(ndexed)' line
      } else if (split[0][0] == 'o' || split[0][0] == 'O' ||
                 split[0][0] == 'i' || split[0][0] == 'I') {
        str[0] = split[0][0];
        break;
      }
    } while (true); //}}}
    // skip the rest until first coordinates //{{{
    fpos_t pos;
    do {
      fgetpos(vcf, &pos); // save vcf file pointer
      fgets(line, sizeof(line), vcf);
      split[0] = strtok(line, " \t");
    } while (split[0][0] < '0' || split[0][0] > '9');
    fsetpos(vcf, &pos); // restore vcf file pointer to before the first coordinate line //}}}
    //}}}

    if (str[0] == 'o' || str[0] == 'O') { // ordered timesteps //{{{
      indexed = false;

      // set all bead types to 'use'
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        (*BeadType)[i].Use = true;
      }
      (*Counts).Beads = (*Counts).BeadsInVsf;
      //}}}
    } else if (str[0] == 'i' || str[0] == 'I') { // indexed timesteps //{{{
      // TODO: check that if one beads with given name is there, all with the same name are
      indexed = true;
      // set all bead types to 'do not use' //{{{
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        (*BeadType)[i].Use = false;
      } //}}}

      // read data
      // the first coordinate line
      fgets(line, sizeof(line), vcf);
      strcpy(line, TrimLine(line)); // trim excess whitespace
      // for lines containing only whitespace
      if (strlen(line) == 1) {
        fprintf(stderr, "\033[1;31m");
        fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - blank line instead of a first coordinate line\n\n", vcf_file);
        fprintf(stderr, "\033[0m");
        exit(1);
      }

      // split the line into array //{{{
      split[0] = strtok(line, " \t");
      int i = 0;
      while (i < 29 && split[i] != NULL) {
        i++;
        split[i] = strtok(NULL, " \t");
      } //}}}

      // first split is bead index; read data till there is <int> at the beginning of the line //{{{
      while (split[0][0] >= '0' && split[0][0] <= '9') {
        int type = (*Bead)[atoi(split[0])].Type;
        (*BeadType)[type].Use = true;
        char test;
        if ((test = getc(vcf)) == EOF) {
          break;
        }
        ungetc(test, vcf);
        // read line
        fgets(line, sizeof(line), vcf);
        strcpy(line, TrimLine(line)); // trim excess whitespace
        // split the line into array //{{{
        split[0] = strtok(line, " \t");
        i = 0;
        while (i < 29 && split[i] != NULL) {
          i++;
          split[i] = strtok(NULL, " \t");
        } //}}}
      } //}}}

      // count the number of beads in vcf //{{{
      (*Counts).Beads = 0;
      for (i = 0; i < (*Counts).TypesOfBeads; i++) {
        if ((*BeadType)[i].Use) {
          (*Counts).Beads += (*BeadType)[i].Number;
        }
      } //}}}
    //}}}
    } else { // error //{{{
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - missing 'i(ndexed)' or 'o(rdered)' keyword\n\n", vcf_file);
      fprintf(stderr, "\033[0m");
      exit(1);
    } //}}}
    fclose(vcf);

  // use all beads if no vcf file provided
  } else {
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      (*BeadType)[i].Use = true;
    }
    (*Counts).Beads = (*Counts).BeadsInVsf;
  } //}}}

  // determine which molecule types are present in vcf //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Use = false;
  }
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // find what bead types are in vcf
    int test = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int type = (*MoleculeType)[i].BType[j];
      if ((*BeadType)[type].Use) {
        test++;
      }
    }
    if (test != 0) {
      (*MoleculeType)[i].Use = true;
    }
  } //}}}

  // remove unused molecules from Molecule struct  //{{{
  int mols_in_vcf = 0;
  int count = 0;
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int type = (*Molecule)[i].Type;
    if ((*MoleculeType)[type].Use) {
      (*Molecule)[count].Type = (*Molecule)[i].Type;
      for (int j = 0; j < (*MoleculeType)[type].nBeads; j++) {
        (*Molecule)[count].Bead[j] = (*Molecule)[i].Bead[j];
      }
      mols_in_vcf++;
      count++;
    }
  }
  (*Counts).Molecules = mols_in_vcf; //}}}

  // copy molecule type names //{{{
  char **name;
  name = calloc(moltype_alloced, sizeof(char *));
  for (int i = 0; i < moltype_alloced; i++) {
    name[i] = calloc(16, sizeof(char));
  }
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    strcpy(name[i], (*MoleculeType)[i].Name);
  } //}}}

  // remove unused molecule types //{{{
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    if ((*MoleculeType)[i].Use) {
      if (i != count) {
        strcpy((*MoleculeType)[count].Name, (*MoleculeType)[i].Name);
        (*MoleculeType)[count].Number = (*MoleculeType)[i].Number;
        (*MoleculeType)[count].nBeads = (*MoleculeType)[i].nBeads;
        for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
          (*MoleculeType)[count].Bead[j] = (*MoleculeType)[i].Bead[j];
        }
        (*MoleculeType)[count].nBonds = (*MoleculeType)[i].nBonds;
        for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
          (*MoleculeType)[count].Bond[j][0] = (*MoleculeType)[i].Bond[j][0];
          (*MoleculeType)[count].Bond[j][1] = (*MoleculeType)[i].Bond[j][1];
        }
        (*MoleculeType)[count].nBTypes = (*MoleculeType)[i].nBTypes;
        for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
          (*MoleculeType)[count].BType[j] = (*MoleculeType)[i].BType[j];
        }
      }
      count++;
    }
  }
  (*Counts).TypesOfMolecules = count; //}}}

  // correct molecule ids (in case of some molecule type absent from vcf) //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int type = (*Molecule)[i].Type;
    (*Molecule)[i].Type = FindMoleculeType(name[type], *Counts, *MoleculeType);
  } //}}}

  // remove unused beads from bonds in molecules //{{{
  // helper array indicating if a given type was already changed //{{{
  int mol_type[(*Counts).TypesOfMolecules];
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    mol_type[i] = -1;
  } //}}}
  // go through all molecules
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int mtype = (*Molecule)[i].Type;
    // is this molecule type done?
    if (mol_type[mtype] == -1) {
      count = 0;
      // go through all bonds and also find initial minimal id //{{{
      int min1 = 1000000;
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        if ((*MoleculeType)[mtype].Bond[j][0] < min1) {
          min1 = (*MoleculeType)[mtype].Bond[j][0];
        }
        int id1 = (*Molecule)[i].Bead[(*MoleculeType)[mtype].Bond[j][0]];
        int btype1 = (*Bead)[id1].Type;
        int id2 = (*Molecule)[i].Bead[(*MoleculeType)[mtype].Bond[j][1]];
        int btype2 = (*Bead)[id2].Type;
        if ((*BeadType)[btype1].Use && (*BeadType)[btype2].Use) {
          (*MoleculeType)[mtype].Bond[count][0] = (*MoleculeType)[mtype].Bond[j][0];
          (*MoleculeType)[mtype].Bond[count][1] = (*MoleculeType)[mtype].Bond[j][1];
          count++;
        }
      } //}}}
      // change number of bonds
      (*MoleculeType)[mtype].nBonds = count;
      // find minimum id in Bond array //{{{
      int min2 = 1000000;
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        if ((*MoleculeType)[mtype].Bond[j][0] < min2) {
          min2 = (*MoleculeType)[mtype].Bond[j][0];
        }
      } //}}}
      // minimize ids in Bond array //{{{
      for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
        (*MoleculeType)[mtype].Bond[j][0] += min1 - min2;
        (*MoleculeType)[mtype].Bond[j][1] += min1 - min2;
      } //}}}
      // go through all beads to correct MoleculeType[].Bead array //{{{
      count = 0;
      for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
        int btype = (*MoleculeType)[mtype].Bead[j];
        if ((*BeadType)[btype].Use) {
          (*MoleculeType)[mtype].Bead[count] = (*MoleculeType)[mtype].Bead[j];
          count++;
        }
      } //}}}
      // mark this molecule type as done
      mol_type[mtype] = i;
    }
  } //}}}

  // remove unused bead types in molecules (Molecule[].Bead arrays) //{{{
  int beads_in_mols[(*Counts).TypesOfMolecules];
  for (int i = 0; i < (*Counts).Molecules; i++) {
    count = 0;
    int mtype = (*Molecule)[i].Type;
    for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
      int type = (*Bead)[(*Molecule)[i].Bead[j]].Type;
      if ((*BeadType)[type].Use) {
        (*Molecule)[i].Bead[count] = (*Molecule)[i].Bead[j];
        (*Bead)[(*Molecule)[i].Bead[j]].Molecule = i;
        count++;
      }
    }
    beads_in_mols[mtype] = count;
  }
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nBeads = beads_in_mols[i];
  } //}}}

  // remove unused bead types in molecules (MoleculeType[].BType arrays) //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    count = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int type = (*MoleculeType)[i].BType[j];
      if ((*BeadType)[type].Use) {
        (*MoleculeType)[i].BType[count] = type;
        count++;
      }
    }
    (*MoleculeType)[i].nBTypes = count;
  } //}}}

  // remove unused beads from Bead struct //{{{
  count = 0;
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int type = (*Bead)[i].Type;
    if ((*BeadType)[type].Use) {
      (*Bead)[count].Type = (*Bead)[i].Type;
      (*Bead)[count].Molecule = (*Bead)[i].Molecule;
      (*Bead)[count].Index = (*Bead)[i].Index;
      (*Index)[(*Bead)[i].Index] = count;
      int id = (*Bead)[count].Molecule;
      if (id != -1) {
        type = (*Molecule)[id].Type;
        for (int j = 0; j < (*MoleculeType)[type].nBeads; j++) {
          if ((*Molecule)[id].Bead[j] == (*Bead)[count].Index) {
            (*Molecule)[id].Bead[j] = count;
          }
        }
      }
      count++;
    }
  } //}}}

  // copy bead type names //{{{
  // first, free name[] ...to prevent memory leak
  for (int i = 0; i < moltype_alloced; i++) {
    free(name[i]);
  }
  // second, realloc name
  int beadtype_alloced = (*Counts).TypesOfMolecules;
  name = realloc(name, beadtype_alloced*sizeof(char *));
  // third, calloc name[]
  for (int i = 0; i < beadtype_alloced; i++) {
    name[i] = calloc(16, sizeof(char));
  }
  // fourth, copy names
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    strcpy(name[i], (*MoleculeType)[i].Name);
  } //}}}

  // remove unused bead types from BeadType struct and molecule bonds //{{{
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if ((*BeadType)[i].Use) {
      if (count != i) {
        strcpy((*BeadType)[count].Name, (*BeadType)[i].Name);
        (*BeadType)[count].Number = (*BeadType)[i].Number;
        (*BeadType)[count].Charge = (*BeadType)[i].Charge;
        (*BeadType)[count].Mass = (*BeadType)[i].Mass;
        (*BeadType)[count].Use = (*BeadType)[i].Use;
        for (int j = 0; j < (*Counts).TypesOfMolecules; j++) {
          for (int k = 0; k < (*MoleculeType)[j].nBTypes; k++) {
            if ((*MoleculeType)[j].BType[k] == i) {
              (*MoleculeType)[j].BType[k] = count;
            }
          }
          for (int k = 0; k < (*MoleculeType)[j].nBeads; k++) {
            if ((*MoleculeType)[j].Bead[k] == i) {
              (*MoleculeType)[j].Bead[k] = count;
            }
          }
        }
        for (int j = 0; j < (*Counts).Beads; j++) {
          if ((*Bead)[j].Type == i) {
            (*Bead)[j].Type = count;
          }
        }
      }
      count++;
    }
  } //}}}

  (*Counts).TypesOfBeads = count;

  // free name array //{{{
  for (int i = 0; i < beadtype_alloced; i++) {
    free(name[i]);
  }
  free(name); //}}}

  // free extraneous arrays to prevent memory leak //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = (*MoleculeType)[i].nBonds; j < max_bonds; j++) {
      free((*MoleculeType)[i].Bond[j]);
    }
  }
  for (int i = (*Counts).TypesOfMolecules; i < moltype_alloced; i++) {
    for (int j = 0; j < max_bonds; j++) {
      free((*MoleculeType)[i].Bond[j]);
    }
    free((*MoleculeType)[i].Bond);
    free((*MoleculeType)[i].BType);
  }
  for (int i = (*Counts).Molecules; i < max_mol; i++) {
    free((*Molecule)[i].Bead);
  } //}}}

  // assign default charge (0) and mass (1) to bead types if not read from vsf //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if ((*BeadType)[i].Charge == 1000) {
      (*BeadType)[i].Charge = 0;
    }
    if ((*BeadType)[i].Mass == 1000) {
      (*BeadType)[i].Mass = 1;
    }
  } //}}}

  // calculate molecule mass //{{{
  count = 0; // to store id of the first molecule id of the given type
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int id = (*Molecule)[count].Bead[j];
      int type = (*Bead)[id].Type;
      (*MoleculeType)[i].Mass += (*BeadType)[type].Mass;
    }
    count += (*MoleculeType)[i].Number;
  } //}}}

  // calculate number of (un)bonded beads //{{{
  (*Counts).Unbonded = 0;
  (*Counts).Bonded = 0;
  for (int i = 0; i < (*Counts).Beads; i++) {
    if ((*Bead)[i].Molecule == -1) {
      (*Counts).Unbonded++;
    } else {
      (*Counts).Bonded++;
    }
  } //}}}

  // allocate Bead[].Aggregate array //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(1, sizeof(int));
    (*Bead)[i].Aggregate[0] = -1;
  } //}}}

  // fill Molecule[].Aggregate //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Molecule)[i].Aggregate = -1;
  } //}}}

  // set all molecule & bead types to be unused //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Use = false;
    (*MoleculeType)[i].Write = false;
  }
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i].Use = false;
    (*BeadType)[i].Write = false;
  } //}}}

  // allocate memory for MoleculeType[].Angles (just to free later) //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nAngles = 0;
    (*MoleculeType)[i].Angle = calloc(1, sizeof(int *));
  } //}}}

  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, vsf_file);

  return indexed;
} //}}}

// ReadTimestepPreamble() //{{{
// TODO: add proper checking for the first coor line via CheckVtfCoordinateLine()
int ReadTimestepPreamble(bool *indexed, char *input_coor, FILE *vcf_file, char **stuff, bool quit) {
  // save pointer position in vcf_file
  fpos_t position;
  fgetpos(vcf_file, &position);

  (*stuff)[0] = '\0'; // empty the array
  int words, count_lines = 0;
  char split[30][100];
  bool timestep = false;
  do {
    char line[LINE], line2[LINE];
    fgets(line, sizeof(line), vcf_file);
    if (feof(vcf_file)) {
      return -1;
    }
    strcpy(line2, line);
    words = SplitLine(split, line, " \t");
    // comment line - copy to stuff array
    if (split[0][0] == '#') {
      char temp[LINE];
      strcpy(temp, *stuff);
      sprintf(temp, "%s%s", *stuff, line2);
      strcpy(*stuff, temp);
    // error if not timestep or pbc line, blank line, or a double (i.e., start of the coordinate lines)
    } else if (quit && (split[0][0] != 't' && split[0][0] != 'i' && split[0][0] != 'o' &&
               words != 0 && strcasecmp(split[0], "pbc") != 0 && !IsDouble(split[0]))) {
      fprintf(stderr, "\033[1;31m");
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - unrecognised line in the timestep preamble\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    // test for a t(imestep) i(ndexed)/o(rdered) line //{{{
    if ((words > 1 && split[0][0] == 't' && split[1][0] == 'i') ||
        split[0][0] == 'i') {
      timestep = true;
      *indexed = true; // indexed timestep present
    } else if ((words > 1 && split[0][0] == 't' && split[1][0] == 'o') ||
               split[0][0] == 'o' || (words == 1 && split[0][0] == 't')){
      timestep = true;
      *indexed = false; // ordered timestep present
    } //}}}
    count_lines++;
  } while (words == 0 || !IsDouble(split[0]));
  count_lines--; // the last counted line contained the first coordinate line
  // error - missing timestep line //{{{
  if (quit && !timestep) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing [t(imestep)] o(rdered)/i(indexed) line\n");
    if (indexed) {
      fprintf(stderr, "       or more coordinate lines than in the first timestep\n");
    }
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  fsetpos(vcf_file, &position); // restore pointer position
  return count_lines;
} //}}}

// ReadCoordinates() //{{{
/**
 * Function reading coordinates from .vcf file with indexed timesteps (\ref IndexedCoorFile).
 */
void ReadCoordinates(bool indexed, char *input_coor, FILE *vcf_file, COUNTS Counts, int *Index, BEAD **Bead, char **stuff) {
  bool test_indexed; // to test if the present timestep type is the same as detected by ReadStructure()
  int count_lines = ReadTimestepPreamble(&test_indexed, input_coor, vcf_file, stuff, true);
  // error - wrong type of step (indexed vs. ordered) //{{{
  if (test_indexed != indexed) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    if (test_indexed) {
      fprintf(stderr, " - indexed timestep instead of an ordered one\n");
    } else {
      fprintf(stderr, " - ordered timestep instead of an indexed one\n");
    }
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // skip the preamble lines //{{{
  for (int i = 0; i < count_lines; i++) {
    char line[LINE];
    fgets(line, sizeof(line), vcf_file);
  } //}}}
  if (indexed) { // indexed timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[30][100];
      fgets(line, sizeof(line), vcf_file);
      // error - end of file
      if (feof(vcf_file)) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - unexpected end of file\n");
        fprintf(stderr, "       possibly fewer beads in a timestep than in the first timestep\n");
        ResetColour(STDERR_FILENO);
        exit(1);
      }
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <int> <double> <double> <double>
      if (words < 4 || !IsInteger(split[0]) || !IsDouble(split[1]) ||
          !IsDouble(split[2]) || !IsDouble(split[3])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read a coordinate line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      int index = atoi(split[0]);
      // bead coordinates
      (*Bead)[Index[index]].Position.x = atof(split[1]);
      (*Bead)[Index[index]].Position.y = atof(split[2]);
      (*Bead)[Index[index]].Position.z = atof(split[3]);
    } //}}}
  } else { // ordered timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[30][100];
      fgets(line, sizeof(line), vcf_file);
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <double> <double> <double>
      if (words < 3 || !IsDouble(split[0]) ||
          !IsDouble(split[1]) || !IsDouble(split[2])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read coordinates\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      // bead coordinates
      (*Bead)[i].Position.x = atof(split[0]);
      (*Bead)[i].Position.y = atof(split[1]);
      (*Bead)[i].Position.z = atof(split[2]);
    }
  } //}}}
} //}}}

// ReadVcfCoordinates() //{{{
/**
 * Function reading coordinates from .vcf file with indexed timesteps (\ref IndexedCoorFile).
 */
void ReadVcfCoordinates(bool indexed, char *input_coor, FILE *vcf_file, COUNTS Counts, int *Index, BEAD **Bead, char **stuff) {
  bool test_indexed; // to test if the present timestep type is the same as detected by ReadStructure()
  int count_lines = ReadTimestepPreamble(&test_indexed, input_coor, vcf_file, stuff, true);
  // error - wrong type of step (indexed vs. ordered) //{{{
  if (test_indexed != indexed) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    if (test_indexed) {
      fprintf(stderr, " - indexed timestep instead of an ordered one\n");
    } else {
      fprintf(stderr, " - ordered timestep instead of an indexed one\n");
    }
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // skip the preamble lines //{{{
  for (int i = 0; i < count_lines; i++) {
    char line[LINE];
    fgets(line, sizeof(line), vcf_file);
  } //}}}
  if (indexed) { // indexed timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[30][100];
      fgets(line, sizeof(line), vcf_file);
      // error - end of file
      if (feof(vcf_file)) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - unexpected end of file\n");
        fprintf(stderr, "       possibly fewer beads in a timestep than in the first timestep\n");
        ResetColour(STDERR_FILENO);
        exit(1);
      }
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <int> <double> <double> <double>
      if (words < 4 || !IsInteger(split[0]) || !IsDouble(split[1]) ||
          !IsDouble(split[2]) || !IsDouble(split[3])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read a coordinate line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      int index = atoi(split[0]);
      // bead coordinates
      (*Bead)[Index[index]].Position.x = atof(split[1]);
      (*Bead)[Index[index]].Position.y = atof(split[2]);
      (*Bead)[Index[index]].Position.z = atof(split[3]);
    } //}}}
  } else { // ordered timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[30][100];
      fgets(line, sizeof(line), vcf_file);
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <double> <double> <double>
      if (words < 3 || !IsDouble(split[0]) ||
          !IsDouble(split[1]) || !IsDouble(split[2])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read coordinates\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      // bead coordinates
      (*Bead)[i].Position.x = atof(split[0]);
      (*Bead)[i].Position.y = atof(split[1]);
      (*Bead)[i].Position.z = atof(split[2]);
    }
  } //}}}
} //}}}

// SkipCoor() //{{{
/**
 * Function to skip one timestep in coordinates file. It works with both
 * indexed and ordered vcf files.
 */
bool SkipCoor(FILE *vcf_file, COUNTS Counts, char **stuff) {

  bool error = false;
  // initial stuff //{{{
  (*stuff)[0] = '\0'; // no comment line
  char line[LINE];
  fpos_t position;
  do {
    fgetpos(vcf_file, &position); // save pointer position
    fgets(line, sizeof(line), vcf_file);

    while (line[0] == '\t' || line[0] == ' ') {
      int i;
      for (i = 1; line[i] != '\n'; i++) {
        line[i-1] = line[i];
      }
      line[i-1] = '\n';
      line[i] = '\0';
    }

    if (line[0] == '#') {
      strcat(*stuff, line);
    }
  } while (line[0] < '0' || line[0] > '9');

  // return file pointer to before the first coordinate line
  fsetpos(vcf_file, &position); //}}}
  for (int i = 0; i < Counts.Beads; i++) {
    int test;
    while ((test = getc(vcf_file)) != '\n' && test != EOF)
      ;
    // premature end of file
    if (test == EOF) {
      error = true;
      break;
    }
  }
  if (!error) {
    getc(vcf_file);
  }
  return error;
} //}}}

// SkipVcfCoor() //{{{
/**
 * Function to skip one timestep in coordinates file. It works with both
 * indexed and ordered vcf files.
 */
bool SkipVcfCoor(FILE *vcf_file, char *input_coor, COUNTS Counts, char **stuff) {

  bool error = false;
//// initial stuff //{{{
//(*stuff)[0] = '\0'; // no comment line
//char line[LINE];
//fpos_t position;
//do {
//  fgetpos(vcf_file, &position); // save pointer position
//  fgets(line, sizeof(line), vcf_file);

//  while (line[0] == '\t' || line[0] == ' ') {
//    int i;
//    for (i = 1; line[i] != '\n'; i++) {
//      line[i-1] = line[i];
//    }
//    line[i-1] = '\n';
//    line[i] = '\0';
//  }

//  if (line[0] == '#') {
//    strcat(*stuff, line);
//  }
//} while (line[0] < '0' || line[0] > '9');

//// return file pointer to before the first coordinate line
//fsetpos(vcf_file, &position); //}}}
  bool rubbish; // testing timestep type - not used here
  int count_lines = ReadTimestepPreamble(&rubbish, input_coor, vcf_file, stuff, false);
  // skip the preamble lines //{{{
  for (int i = 0; i < count_lines; i++) {
    char line[LINE];
    fgets(line, sizeof(line), vcf_file);
  } //}}}
  for (int i = 0; i < Counts.Beads; i++) {
    int test;
    while ((test = getc(vcf_file)) != '\n' && test != EOF)
      ;
    // premature end of file
    if (test == EOF) {
      error = true;
      break;
    }
  }
  if (!error) {
    getc(vcf_file);
  }
  return error;
} //}}}

// ReadAggregates() //{{{
/**
 * Function reading information about aggregates from `.agg` file (\ref AggregateFile) generated by Aggregates utility.
 */
bool ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts, AGGREGATE **Aggregate,
                    BEADTYPE *BeadType, BEAD **Bead,
                    MOLECULETYPE *MoleculeType, MOLECULE **Molecule, int *Index) {

  bool error = false;
  char line[LINE], split[30][100], delim[8];
  // read (Last) Step line
  fgets(line, sizeof(line), fr);
  strcpy(delim, " \t"); // delimiters for SplitLine()
  int words = SplitLine(split, line, delim);

  // is there a Step? i.e., isn't this the line 'Last Step'?
  if (split[0][0] == 'L') {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: premature end of ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " file");
    ResetColour(STDERR_FILENO);
    exit(1);
  // Step line must be 'Step: <int>'
  } else if (words < 2 || !IsInteger(split[1])) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - wrong 'Step' line\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  } else {
    // initialize array of number of aggregates per bead //{{{
    for (int i = 0; i < (*Counts).Beads; i++) {
      (*Bead)[i].nAggregates = 0;
    } //}}}

    // get number of aggregates //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // number of aggregates must be <int>
    if (words != 0 && !IsInteger(split[0])) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", agg_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - number of aggregates must be a whole number\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*Counts).Aggregates = atoi(split[0]); //}}}

    // skip blank line
    fgets(line, sizeof(line), fr);

    // go through all aggregates
    for (int i = 0; i < (*Counts).Aggregates; i++) {
      // read molecules in Aggregate 'i' //{{{
      fscanf(fr, "%d :", &(*Aggregate)[i].nMolecules);
      for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
        int mol;
        fscanf(fr, "%d", &mol);
        mol--; // in agg file the numbers correspond to vmd

        (*Aggregate)[i].Molecule[j] = mol;
        (*Molecule)[mol].Aggregate = i;
      }

      while (getc(fr) != '\n')
       ; //}}}

      // read monomeric beads in Aggregate 'i' //{{{
      int count;
      fscanf(fr, "%d :", &count);
      (*Aggregate)[i].nMonomers = 0; // their number will be counted according to which beads are in vcf
      for (int j = 0; j < count; j++) {
        int vsf_id;
        fscanf(fr, "%d", &vsf_id); // monomer index in vsf file
        int id = Index[vsf_id]; // monomer index in Bead structure
        if (id > -1) {
          int beads = (*Aggregate)[i].nMonomers++;
          (*Aggregate)[i].Monomer[beads] = id;
          (*Bead)[id].nAggregates++;
          (*Bead)[id].Aggregate = realloc((*Bead)[id].Aggregate, (*Bead)[id].nAggregates*sizeof(int));
          (*Bead)[id].Aggregate[(*Bead)[id].nAggregates-1] = i;
        }
      }

      while (getc(fr) != '\n')
       ; //}}}
    }

    // skip blank line at the end of every entry
    while (getc(fr) != '\n')
      ;

    // fill Aggregate[].Bead, Bead[].Aggregate, and Molecule[].Aggregate arrays //{{{
    for (int i = 0; i < (*Counts).Aggregates; i++) {
      (*Aggregate)[i].nBeads = 0;

      for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
        int mol = (*Aggregate)[i].Molecule[j];
        (*Molecule)[mol].Aggregate = i;
        for (int k = 0; k < MoleculeType[(*Molecule)[mol].Type].nBeads; k++) {
          int bead = (*Molecule)[mol].Bead[k];
          (*Aggregate)[i].Bead[(*Aggregate)[i].nBeads+k] = bead;
          (*Bead)[bead].Aggregate[(*Bead)[bead].nAggregates] = i;
          (*Bead)[bead].nAggregates++;
        }

        // increment number of beads in aggregate 'i'
        (*Aggregate)[i].nBeads += MoleculeType[(*Molecule)[mol].Type].nBeads;
      }
    }

    // fill Bead[].Aggregate for monomeric Beads
//  for (int i = 0, i < (*Counts).Bead)
    //}}}

    // calculate aggregates' masses //{{{
    for (int i = 0; i < (*Counts).Aggregates; i++) {
      (*Aggregate)[i].Mass = 0;

      for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
        int type = (*Molecule)[(*Aggregate)[i].Molecule[j]].Type;

        (*Aggregate)[i].Mass += MoleculeType[type].Mass;
      }
    } //}}}
  }

  return error;
} //}}}

// ReadFieldPbc() //{{{
/**
 * Function reading box size from the first line of the FIELD-like file.
 */
void ReadFieldPbc(char *field, VECTOR *BoxLength) {
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}

  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t"); // delimiters for SplitLine()

  // read box size //{{{
  fgets(line, sizeof(line), fr);
  int words = SplitLine(split, line, delim);
  // first line of the FIELD must be: <double> <double> <double> ...whatever
  // Error if:
  // 1) too few strings
  // 2) enough strings, but not numbers
  if (words < 3 || // 1
      !IsPosDouble(split[0]) || // 2) first <double>
      !IsPosDouble(split[1]) || // 2) first <double>
      !IsPosDouble(split[2])) { // 2) first <double>
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - first line must start with box size (i.e., three positive numbers)\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  }
  (*BoxLength).x = atof(split[0]);
  (*BoxLength).y = atof(split[1]);
  (*BoxLength).z = atof(split[2]); //}}}

  fclose(fr);
} //}}}

// ReadFieldBeadType() //{{{
/**
 * Function reading the species section of the FIELD-like file.
 */
void ReadFieldBeadType(char *field, COUNTS *Counts, BEADTYPE **BeadType, BEAD **Bead) {

  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}

  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t"); // delimiters for SplitLine()

  // read number of bead types //{{{
  bool missing = true; // is 'species' keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    int words = SplitLine(split, line, delim);
    if (strcasecmp(split[0], "species") == 0) {
      missing = false;
      // check if the next string is a number
      if (words < 2 ||            // missing next string
          !IsInteger(split[1])) { // next string isn't a number
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - missing number of species\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'species' keyword\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}

  *BeadType = calloc((*Counts).TypesOfBeads,sizeof(struct BeadType));

  // read info about bead types //{{{
  (*Counts).Unbonded = 0;
  (*Counts).Beads = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    fgets(line, sizeof(line), fr);
    int words = SplitLine(split, line, delim);

    // Error on the line: //{{{
    // 1) empty line
    // 2) less then four strings
    // 3) second string isn't a positive double (mass)
    // 4) third string isn't a double (charge)
    // 5) fifth string isn't an integer (unbonded beads)
    if (words == 1 && split[0][0] == '\0') { // 1)
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - no blank spaces permitted in the species section\n\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    } else if (words < 4 ||                  // 2)
               !IsPosDouble(split[1]) ||     // 3)
               !IsDouble(split[2]) ||        // 4)
               !IsInteger(split[3])) {       // 5)
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - wrong species line\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}

    strcpy((*BeadType)[i].Name, split[0]);
    (*BeadType)[i].Mass = atof(split[1]);
    (*BeadType)[i].Charge = atof(split[2]);
    (*BeadType)[i].Number = atoi(split[3]);
    (*Counts).Unbonded += (*BeadType)[i].Number;
  } //}}}

  (*Counts).Beads = (*Counts).Unbonded;

  // allocate & fill Bead array //{{{
  *Bead = calloc((*Counts).Unbonded, sizeof(struct Bead));
  int count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    for (int j = 0; j < (*BeadType)[i].Number; j++) {
      (*Bead)[count].Type = i;
      (*Bead)[count].Molecule = -1;
      (*Bead)[count].Index = count;
      count++;
    }
  } //}}}

  fclose(fr);
} //}}}

// ReadFieldMolecules() //{{{
/**
 * Function reading the molecules section of the FIELD-like file.
 */
void ReadFieldMolecules(char *field, COUNTS *Counts,
                        BEADTYPE **BeadType, BEAD **Bead,
                        MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
                        PARAMS **bond_type, PARAMS **angle_type) {

  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}

  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t"); // delimiters for SplitLine()

  // read number of molecule types //{{{
  bool missing = true; // is molecule keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    int words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "molecule", 8) == 0) {
      missing = false;
      // error - next string isn't a number
      if (!IsInteger(split[1])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - missing number of molecule types\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfMolecules = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'molecule' line\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}

  // allocate molecule type struct //{{{
  *MoleculeType = calloc((*Counts).TypesOfMolecules,sizeof(struct MoleculeType));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Number = 0;
    (*MoleculeType)[i].nBeads = 0;
    (*MoleculeType)[i].nBonds = 0;
    (*MoleculeType)[i].BType = calloc((*Counts).TypesOfBeads, sizeof(int));
  } //}}}

  // test proper number of 'finish' keywords //{{{
  fpos_t position;
  fgetpos(fr, &position); // save pointer position
  int count = 0;
  while(fgets(line, sizeof(line), fr)) {
    SplitLine(split, line, delim);
    if (strcmp(split[0], "finish") == 0) {
      count++;
    }
  }
  if (count < (*Counts).TypesOfMolecules) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'finish' keyword\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  fsetpos(fr, &position); //}}}

  // read info about molecule types //{{{
  fgetpos(fr, &position); // save pointer position
  (*Counts).TypesOfBonds = 0;
  (*Counts).TypesOfAngles = 0;
  (*bond_type) = calloc(1, sizeof(struct Params));
  (*angle_type) = calloc(1, sizeof(struct Params));
  // stored bond & angle types - temporary
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // read i-th molecule name //{{{
    fgets(line, sizeof(line), fr);
    SplitLine(split, line, delim);
    if (split[0][0] == '\0') {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - missing molecule name\n\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
    strcpy((*MoleculeType)[i].Name, split[0]); //}}}
    // read number of the molecules //{{{
    fgets(line, sizeof(line), fr);
    int words = SplitLine(split, line, delim);
    if (strcasecmp(split[0], "nummols") != 0) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - missing 'nummols' keyword\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - 'nummols' must be followed by a non-negative integer\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].Number = atoi(split[1]); //}}}
    // read number of beads in the molecule //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "beads", 4) != 0) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - missing 'beads' keyword\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - 'beads' must be followed by a non-negative integer\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].nBeads = atoi(split[1]);
    (*MoleculeType)[i].Bead = calloc((*MoleculeType)[i].nBeads, sizeof(int)); //}}}
    // read bead info //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line, delim);
      // bead line must be 'name <double> <double> <double>'
      if (words < 4 ||
          !IsDouble(split[1]) || !IsDouble(split[2]) || !IsDouble(split[3])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - wrong bead line in molecule ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      int type = FindBeadType(split[0], *Counts, *BeadType);
      if (type == -1) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - non-existent bead name ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", split[0]);
        RedText(STDERR_FILENO);
        fprintf(stderr, "in molecule ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
        ResetColour(STDERR_FILENO);
        ErrorBeadType(*Counts, *BeadType);
        exit(1);
      }
      (*MoleculeType)[i].Bead[j] = type;
    } //}}}
    // read number of bonds in the molecule //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    if (strncasecmp(split[0], "bonds", 4) != 0) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - missing 'bonds' keyword\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    if (words < 2 || !IsInteger(split[1])) {
      RedText(STDERR_FILENO);
      YellowText(STDERR_FILENO);
      RedText(STDERR_FILENO);
      ResetColour(STDERR_FILENO);
      fprintf(stderr, "\nError: \033[1;33m%s\033[1;31m - 'bonds' must be followed by a non-negative integer\n", field);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].nBonds = atoi(split[1]);
    (*MoleculeType)[i].Bond = calloc((*MoleculeType)[i].nBonds, sizeof(int *));
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      (*MoleculeType)[i].Bond[j] = calloc(3, sizeof(int));
      (*MoleculeType)[i].Bond[j][2] = -1; // no bond type assigned
    } //}}}
    // read bond info //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      fgets(line, sizeof(line), fr);
      words = SplitLine(split, line, delim);
      // error - bead line must be '<type> <bead id> <bead id> <bond strenthe> <distance>' //{{{
      if (words < 5 ||
          !IsInteger(split[1]) || atoi(split[1]) == 0 || // bead ids must start from 1
          !IsInteger(split[2]) || atoi(split[2]) == 0 ||
          !IsPosDouble(split[3]) || !IsPosDouble(split[4])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - wrong bond line in molecule\n");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      (*MoleculeType)[i].Bond[j][0] = atoi(split[1]) - 1;
      (*MoleculeType)[i].Bond[j][1] = atoi(split[2]) - 1;
      // find if the bond type is new
      bool known = false;
      for (int k = 0; k < (*Counts).TypesOfBonds; k++) {
        if ((*bond_type)[k].a == atof(split[3]) && (*bond_type)[k].b == atof(split[4])) {
          known = true;
          break;
        }
      }
      if (!known) {
        (*Counts).TypesOfBonds++;
        (*bond_type) = realloc((*bond_type), (*Counts).TypesOfBonds*sizeof(struct Params));
        (*bond_type)[(*Counts).TypesOfBonds-1].a = atof(split[3]);
        (*bond_type)[(*Counts).TypesOfBonds-1].b = atof(split[4]);
      }
    } //}}}
    // read angle info if present //{{{
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // are angles present?
    if (strncasecmp(split[0], "angles", 5) == 0) {
      // get number of angles //{{{
      if (words < 2 || !IsInteger(split[1])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'angles' must be followed by a non-negative integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*MoleculeType)[i].nAngles = atoi(split[1]);
      (*MoleculeType)[i].Angle = calloc((*MoleculeType)[i].nAngles, sizeof(int *));
      for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
        (*MoleculeType)[i].Angle[j] = calloc(4, sizeof(int));
        (*MoleculeType)[i].Angle[j][3] = -1; // no angle type assigned
      } //}}}
      // get info about angles
      for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // bead line must be '<type> <bead id> <bead id> <bead id> <angle strength> <angle>'
        if (words < 6 ||
            !IsInteger(split[1]) || atoi(split[1]) == 0 || // bead ids must start from 1
            !IsInteger(split[2]) || atoi(split[2]) == 0 ||
            !IsInteger(split[3]) || atoi(split[3]) == 0 ||
            !IsPosDouble(split[4]) || !IsPosDouble(split[5])) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError: ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", field);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - wrong angle line in molecule ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
          ErrorPrintLine(split, words);
          exit(1);
        }
        (*MoleculeType)[i].Angle[j][0] = atoi(split[1]) - 1;
        (*MoleculeType)[i].Angle[j][1] = atoi(split[2]) - 1;
        (*MoleculeType)[i].Angle[j][2] = atoi(split[3]) - 1;
        bool known = false;
        for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
          if ((*angle_type)[k].a == atof(split[4]) && (*angle_type)[k].b == atof(split[5])) {
            known = true;
            break;
          }
        }
        if (!known) {
          (*Counts).TypesOfAngles++;
          (*angle_type) = realloc((*angle_type), (*Counts).TypesOfAngles*sizeof(struct Params));
          (*angle_type)[(*Counts).TypesOfAngles-1].a = atof(split[4]);
          (*angle_type)[(*Counts).TypesOfAngles-1].b = atof(split[5]);
        }
      }
    // error - isn't there more bonds, by any chance?
    } else if (strncasecmp(split[0], "harm", 4) == 0) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - extra bond line in molecule ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}
    (*Counts).Bonded += (*MoleculeType)[i].Number * (*MoleculeType)[i].nBeads;
    (*Counts).Molecules += (*MoleculeType)[i].Number;
    // skip till 'finish' //{{{
    if (strcasecmp(split[0], "finish") == 0) {
      continue;
    }
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  } //}}}

  // return the file pointer to the beginning of the molecules section
  fsetpos(fr, &position);

  // calculate molecule masses //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      (*MoleculeType)[i].Mass += (*BeadType)[btype].Mass;
    }
  } //}}}

  // update the number beads in the system //{{{
  (*Counts).Beads += (*Counts).Bonded;
  (*Counts).BeadsInVsf = (*Counts).Beads;
  *Bead = realloc(*Bead, (*Counts).Beads*sizeof(struct Bead));
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      (*BeadType)[btype].Number += (*MoleculeType)[i].Number;
    }
  } //}}}

  // read coordinates of bonded beads & assign bond and angle types //{{{
  // get the coordinates from FIELD to each molecule
  count = (*Counts).Unbonded;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // skip name, nummols, and beads lines //{{{
    fgets(line, sizeof(line), fr);
    fgets(line, sizeof(line), fr);
    fgets(line, sizeof(line), fr); //}}}
    // read bead lines //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      for (int k = 0; k < (*MoleculeType)[i].Number; k++) {
        int id = count + k * (*MoleculeType)[i].nBeads;
        (*Bead)[id].Position.x = atof(split[1]);
        (*Bead)[id].Position.y = atof(split[2]);
        (*Bead)[id].Position.z = atof(split[3]);
      }
      count++;
    } //}}}
    // set count to the first bead of the next molecule type
    count += ((*MoleculeType)[i].Number - 1) * (*MoleculeType)[i].nBeads;
    // skip bonds line
    fgets(line, sizeof(line), fr);
    // read bond types //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      for (int k = 0; k < (*Counts).TypesOfBonds; k++) {
        if ((*bond_type)[k].a == atof(split[3]) && (*bond_type)[k].b == atof(split[4])) {
          (*MoleculeType)[i].Bond[j][2] = k;
          break;
        }
      }
    } //}}}
    // read angle types //{{{
    for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
      // skip 'angles' line at the beginning of the angle section //{{{
      if (j == 0) {
        fgets(line, sizeof(line), fr);
      } //}}}
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
        if ((*angle_type)[k].a == atof(split[4]) && (*angle_type)[k].b == atof(split[5])) {
          (*MoleculeType)[i].Angle[j][3] = k;
          break;
        }
      }
    } //}}}
    // skip till 'finish' //{{{
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  } //}}}

  // return the file pointer to the beginning of the molecules section
  fsetpos(fr, &position);

  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    SortBonds((*MoleculeType)[i].Bond, (*MoleculeType)[i].nBonds);
    SortAngles((*MoleculeType)[i].Angle, (*MoleculeType)[i].nAngles);
  }

  fclose(fr);
} //}}}

// ReadField() //{{{
/**
 * Function reading the FIELD-like file; it completely fills provided structs,
 * overriding any possible data in there. If the FIELD-like file is a source of
 * some additional structure information, new structs must be used, and the
 * data copied from there afterwards.
 */
void ReadField(char *field, VECTOR *BoxLength, COUNTS *Counts,
               BEADTYPE **BeadType, BEAD **Bead, int **Index,
               MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
               PARAMS **bond_type, PARAMS **angle_type) {

  // read pbc only if required
  if (BoxLength != NULL) {
    ReadFieldPbc(field, BoxLength);
  }
  ReadFieldBeadType(field, Counts, BeadType, Bead);
  ReadFieldMolecules(field, Counts, BeadType, Bead, MoleculeType, Molecule, bond_type, angle_type);

  // allocate Bead[].Aggregate array - needed only to free() //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(10, sizeof(int));
  } //}}}

  // fill MoleculeType[].BType array //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].nBTypes = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int type = (*MoleculeType)[i].Bead[j];
      bool present = false;
      for (int k = 0; k < (*MoleculeType)[i].nBTypes; k++) {
        if (type == (*MoleculeType)[i].BType[k]) {
          present = true;
          break;
        }
      }
      if (!present) {
        (*MoleculeType)[i].BType[(*MoleculeType)[i].nBTypes] = type;
        (*MoleculeType)[i].nBTypes++;
      }
    }
  } //}}}

  // fill Molecule & Bead structs //{{{
  *Molecule = calloc((*Counts).Molecules,sizeof(struct Molecule));
  int count_mol = 0, count_bead = (*Counts).Unbonded;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].Number; j++) {
      (*Molecule)[count_mol].Type = i;
      (*Molecule)[count_mol].Bead = calloc((*MoleculeType)[i].nBeads, sizeof(int));
      for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
        int btype = (*MoleculeType)[i].Bead[k];
        (*Molecule)[count_mol].Bead[k] = count_bead;
        (*Bead)[count_bead].Molecule = count_mol;
        (*Bead)[count_bead].Type = btype;
        (*Bead)[count_bead].Index = count_bead;
        count_bead++;
      }
      count_mol++;
    }
  } //}}}

  // fill Index - I don't thinks it's used right now //{{{
  *Index = calloc((*Counts).Beads, sizeof(int));
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Index)[i] = i;
  } //}}}

  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, field);
} //}}}

// ReadLmpData() //{{{
void ReadLmpData(char *data_file, int *bonds, PARAMS **bond_type,
                 int *angles, PARAMS **angle_type,
                 VECTOR *BoxLength, VECTOR *box_lo, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {
  // open data file //{{{
  FILE *fr;
  if ((fr = fopen(data_file, "r")) == NULL) {
    ErrorFileOpen(data_file, 'r');
    exit(1);
  } //}}}

  // ignore first line (comment) //{{{
  char line[LINE];
  fgets(line, sizeof(line), fr); //}}}

  // read data file header //{{{
  // data file header lines must start with a number (or '#' for comment),
  // therefore read until something else is encountered
  char split[30][100], delim[8];
  int words;
  strcpy(delim, " \t");
  do {
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
    // number of atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsInteger(split[0])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'atoms' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).BeadsInVsf = atoi(split[0]);
      (*Counts).Beads = atoi(split[0]);
    } //}}}
    // number of bonds //{{{
    if (words > 1 && strcmp(split[1], "bonds") == 0) {
      if (!IsInteger(split[0])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'bonds' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      *bonds = atoi(split[0]);
    } //}}}
    // number of angles //{{{
    if (words > 1 && strcmp(split[1], "angles") == 0) {
      if (!IsInteger(split[0])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'angles' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      *angles = atoi(split[0]);
    } //}}}
    // number of bead types //{{{
    if (words > 2 && strcmp(split[1], "atom") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'atom types' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[0]);
    } //}}}
    // number of bond types //{{{
    if (words > 2 && strcmp(split[1], "bond") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'bond types' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBonds = atoi(split[0]);
    } //}}}
    // number of angle types //{{{
    if (words > 2 && strcmp(split[1], "angle") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'angle types' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfAngles = atoi(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "xlo") == 0 && strcmp(split[3], "xhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'xlo xhi' keyword must be preceded by two floats\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).x = atof(split[1]) - atof(split[0]);
      (*box_lo).x = atof(split[0]);
    } //}}}
    // box length in y //{{{
    if (words > 3 && strcmp(split[2], "ylo") == 0 && strcmp(split[3], "yhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'ylo yhi' keyword must be preceded by two floats\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).y = atof(split[1]) - atof(split[0]);
      (*box_lo).y = atof(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "zlo") == 0 && strcmp(split[3], "zhi") == 0) {
      if (!IsDouble(split[0]) || !IsDouble(split[1])) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'zlo zhi' keyword must be preceded by two floats\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).z = atof(split[1]) - atof(split[0]);
      (*box_lo).z = atof(split[0]);
    } //}}}
  } while (words == 0 ||
           split[0][0] == '#' ||
           IsDouble(split[0]) ||
           IsInteger(split[0])); //}}}

  // some error checking //{{{
  if ((*Counts).TypesOfBeads == 0) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'atom types' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*Counts).BeadsInVsf == 0) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'atoms' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).x == 0) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'xlo xhi' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).y == 0) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'ylo yhi' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).z == 0) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'zlo zhi' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}

  // fill something in BeadType struct //{{{
  *BeadType = calloc((*Counts).TypesOfBeads, sizeof(struct BeadType));
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    sprintf((*BeadType)[i].Name, "bead%d", i+1);
    (*BeadType)[i].Use = true;
    (*BeadType)[i].Write = true;
  } //}}}

  // bead struct memory allocation //{{{
  *Bead = calloc((*Counts).Beads, sizeof(struct Bead));
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(1, sizeof(int));
  } //}}}

  *Index = calloc((*Counts).Beads, sizeof(int)); // link between indices in vsf and in program (i.e., opposite of Bead[].Index)
  *bond_type = calloc((*Counts).TypesOfBonds, sizeof(PARAMS));
  *angle_type = calloc((*Counts).TypesOfAngles, sizeof(PARAMS));

  // read body of data file //{{{
  int test,
      *mols, // number of beads in each molecule
      monomer = 0; // monomer beads designated by mol_ID = 0 in data file
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);
    // atom masses //{{{
    if (words > 0 && strcmp(split[0], "Masses") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 2 || !IsInteger(split[0]) || !IsPosDouble(split[1])) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each line in 'Masses' section must start with '<int> <float>'\n");
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*BeadType)[i].Mass = atof(split[1]);
        // if there's a comment at the end of the line, consider it bead name
        if (words > 2 && split[2][0] == '#') {
          if (strlen(split[2]) == 1 && words > 3) { // comment form '# name'
            strcpy((*BeadType)[i].Name, split[3]);
          } else if (strlen(split[2]) > 1) { // comment form '#name'
            for (int j = 0; j < strlen(split[2]); j++) {
              split[2][j] = split[2][j+1];
            }
            strcpy((*BeadType)[i].Name, split[2]);
          }
        }
      }
    } //}}}
    // bond coefficients //{{{
    if (words > 1 && strcmp(split[0], "Bond") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfBonds; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosDouble(split[1]) || !IsPosDouble(split[2])) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each line in 'Bond Coeffs' section must start with '<int> <float> <float>'\n");
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*bond_type)[i].a = atof(split[1]);
        (*bond_type)[i].b = atof(split[2]);
      }
    } //}}}
    // angle coefficients //{{{
    if (words > 1 && strcmp(split[0], "Angle") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfAngles; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosDouble(split[1]) || !IsPosDouble(split[2])) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each line in 'Angle Coeffs' section must start with '<int> <float> <float>'\n");
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*angle_type)[i].a = atof(split[1]);
        (*angle_type)[i].b = atof(split[2]);
      }
    } //}}}
    // atoms section //{{{
    if (words > 0 && strcmp(split[0], "Atoms") == 0) {
      // array for number of beads in each molecule
      mols = calloc((*Counts).Beads, sizeof(int));
      // array for list of beads in each molecule
      // bead_mols[i] ... molecule's id; bead_mols[][i] ... molecule's beads
      int **bead_mols = calloc((*Counts).Beads, sizeof(int *));
      for (int i = 0; i < (*Counts).Beads; i++) {
        bead_mols[i] = calloc(1, sizeof(int));
      }
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // go through atom section to get basic info //{{{
      fpos_t pos; // set file counter
      fgetpos(fr, &pos); // save file pointer
      for (int i = 0; i < (*Counts).Beads; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <id> <mol_id> <btype> <charge> <x> <y> <z>
        // Error - incorrect format //{{{
        if (words < 7 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) || !IsInteger(split[2]) ||
            !IsDouble(split[3]) || !IsDouble(split[4]) || !IsDouble(split[5]) || !IsDouble(split[6])) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each 'Atoms' line must be <id> <mol_id> <bead type> <charge> <x> <y> <z>\n");
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int id = atoi(split[0]) - 1, // in lammps, these start at 1
            mol_id = atoi(split[1]) - 1, // in lammps, molecules start with 1; unbonded atoms can be 0
            btype = atoi(split[2]) - 1; // in lammps, these start at 1

        if (mol_id == -1) { // corresponds to 0 in the data file
          monomer++;
          (*Bead)[id].Molecule = -1;
        } else { // possibly in a molecule (if more beads share its mol_id)
          mols[mol_id]++;
          bead_mols[mol_id] = realloc(bead_mols[mol_id], mols[mol_id]*sizeof(int));
          bead_mols[mol_id][mols[mol_id]-1] = id;
          (*Bead)[id].Molecule = mol_id;
        }
        (*BeadType)[btype].Charge = atof(split[3]);
        (*BeadType)[btype].Number++;
        (*Bead)[id].Position.x = atof(split[4]) - (*box_lo).x;
        (*Bead)[id].Position.y = atof(split[5]) - (*box_lo).y;
        (*Bead)[id].Position.z = atof(split[6]) - (*box_lo).z;
        (*Bead)[id].Type = btype;
        (*Bead)[id].Index = id;
        (*Index)[id] = id;
      } //}}}
      (*Counts).Unbonded = monomer;
      // go through possible molecules and remove 1-bead molecules //{{{
      int count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < (*Counts).Beads; i++) {
        if (mols[i] == 1) {
          int bead = bead_mols[i][0];
          (*Bead)[bead].Molecule = -1;
          (*Counts).Unbonded++;
        } else if (mols[i] > 1){
          (*Counts).Molecules++;
          (*Counts).Bonded += mols[i];
          for (int j = 0; j < mols[i]; j++) {
            int bead = bead_mols[i][j];
            (*Bead)[bead].Molecule = count;
          }
          count++;
        }
      } //}}}
      // TODO: is the 'remove 1-bead...' necessary? Join with 'allocate Molecule struct...'
      // remove 1-bead molecules from mols array and sort bead_mols arrayy according to ascending bead index //{{{
      count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < (*Counts).Beads; i++) {
        if (mols[i] > 1) {
          mols[count] = mols[i];
          // realloc bead_mols[count] to bead_mols[i] as it realloc'd in the first place
          bead_mols[count] = realloc(bead_mols[count], mols[count]*sizeof(int));
          for (int j = 0; j < mols[i]; j++) {
            bead_mols[count][j] = bead_mols[i][j];
          }
          // sort molecules in bead_mols[count][] according to ascending id
          SortArray(&bead_mols[count], mols[i], 0);
          count++;
        }
      } //}}}
      // zeroize unused part of mols array - just to be on the save side //{{{
      for (int i = (*Counts).Molecules; i < (*Counts).Beads; i++) {
        mols[i] = 0;
      } //}}}
      // allocate Molecule struct and fill Molecule[].Bead array //{{{
      (*Molecule) = calloc((*Counts).Molecules, sizeof(struct Molecule));
      for (int i = 0; i < (*Counts).Molecules; i++) {
        (*Molecule)[i].Bead = calloc(mols[i], sizeof(int));
        for (int j = 0; j < mols[i]; j++) {
          (*Molecule)[i].Bead[j] = bead_mols[i][j];
        }
      } //}}}
      // free helper array  //{{{
      for (int i = 0; i < (*Counts).Beads; i++) {
        free(bead_mols[i]);
      }
      free(bead_mols); //}}}
    } //}}}
    // bonds section //{{{
    if (words > 0 && strcmp(split[0], "Bonds") == 0) {
      // allocate helper arrays to hold bond info //{{{
      // number of bonds in each molecule
      int *bonds_per_mol = calloc((*Counts).Molecules, sizeof(int));
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] ... ids of connected beads in molecule 'i'
      // [i][][2] ... bond type
      int ***mol_bonds = calloc((*Counts).Molecules, sizeof(int **));
      for (int i = 0; i < (*Counts).Molecules; i++) {
         mol_bonds[i] = calloc(1, sizeof(int *));
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // read all bonds //{{{
      for (int i = 0; i < *bonds; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <bond id> <bond type> <bead1> <bead2>
        // Error - incorrect format //{{{
        if (words < 4 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3])) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each 'Bonds' line must be <bond id> <bond type> <bead1d> <bead2>\n");
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        // assign molecule to the bond
        int mol = (*Bead)[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != (*Bead)[bead2].Molecule) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - beads in ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "bond %d", atoi(split[0]));
          RedText(STDERR_FILENO);
          fprintf(stderr, "are in different molecules\n\n");
          ResetColour(STDERR_FILENO);
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        bonds_per_mol[mol]++;
        int bond = bonds_per_mol[mol];
        // add bonded beads to the molecule they belong to
        mol_bonds[mol] = realloc(mol_bonds[mol], bond*sizeof(int *));
        mol_bonds[mol][bond-1] = calloc(3, sizeof(int));
        mol_bonds[mol][bond-1][0] = bead1;
        mol_bonds[mol][bond-1][1] = bead2;
        mol_bonds[mol][bond-1][2] = type;
      } //}}}
      // sort bonds according to the id of the first bead in a bond //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        SortBonds(mol_bonds[i], bonds_per_mol[i]);
      } //}}}
      // minimize mol_bonds based on lowest id in each molecule //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int lowest = (*Counts).Beads; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if ((*Molecule)[i].Bead[j] < lowest) {
            lowest = (*Molecule)[i].Bead[j];
          }
        }
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          mol_bonds[i][j][0] -= lowest;
          mol_bonds[i][j][1] -= lowest;
        }
      } //}}}
      // identify molecule type based on bead order and connectivity //{{{
      *MoleculeType = calloc(1, sizeof(struct MoleculeType));
      // number of molecule types differing in connectivity but not in beads - naming purposes
      int *diff_conn = calloc(1, sizeof(int));
      // number of molecule types differing in bead order - naming purposes
      int mtype_bead = 0;
      for (int i = 0; i < (*Counts).Molecules; i++) {
        // is molecule 'i' of known type?
        bool exists = false;
        // do 'i' and given molecule type share connectivity and bead order?
        bool same_conn = true, same_bead = true;
        // molecule type with which 'i' shares bead order - naming purposes
        int type_bead = -1;
        for (int j = 0; j < (*Counts).TypesOfMolecules; j++) { //{{{
          if ((*MoleculeType)[j].nBeads == mols[i] && // same number of molecules?
              (*MoleculeType)[j].nBonds == bonds_per_mol[i]) { // same number of bonds?
            // same connectivity?
            same_conn = true;
            for (int k = 0; k < (*MoleculeType)[j].nBonds; k++) {
              if ((*MoleculeType)[j].Bond[k][0] != mol_bonds[i][k][0] ||
                  (*MoleculeType)[j].Bond[k][1] != mol_bonds[i][k][1] ||
                  (*MoleculeType)[j].Bond[k][2] != mol_bonds[i][k][2]) {
                same_conn = false;
                break;
              }
            }
            // same bead types?
            same_bead = true;
            for (int k = 0; k < (*MoleculeType)[j].nBeads; k++) {
              int btype = (*Bead)[(*Molecule)[i].Bead[k]].Type;
              if ((*MoleculeType)[j].Bead[k] != btype) {
                same_bead = false;
                break;
              }
            }
            if (same_bead && type_bead == -1) {
              type_bead = j;
            }
            // if the molecule has the same connectivity and bead order, it's of known type
            if (same_conn && same_bead) {
              exists = true;
              (*MoleculeType)[j].Number++;
              (*Molecule)[i].Type = j;
              break;
            }
          }
        } //}}}
        // add new type? //{{{
        if (!exists) {
          int mtype = (*Counts).TypesOfMolecules;
          *MoleculeType = realloc(*MoleculeType, (mtype+1)*sizeof(struct MoleculeType));
          diff_conn = realloc(diff_conn, (mtype+1)*sizeof(int));
          diff_conn[mtype] = 0;
          // molecule name
          if (same_bead && !same_conn) { // same beads - mol<type_bead>-b#
            diff_conn[type_bead]++;
            char name[5];
            strcpy(name, (*MoleculeType)[type_bead].Name);
            sprintf((*MoleculeType)[mtype].Name, "%s-b%d", name, diff_conn[type_bead]);
          } else { // same connectivity or both different - mol#
            mtype_bead++;
            sprintf((*MoleculeType)[mtype].Name, "mol%d", mtype_bead);
          }
          (*Molecule)[i].Type = mtype;
          (*MoleculeType)[mtype].Number = 1;
          // copy bead sequence and determine BTypes stuff
          (*MoleculeType)[mtype].nBeads = mols[i];
          (*MoleculeType)[mtype].Bead = calloc((*MoleculeType)[mtype].nBeads, sizeof(int));
          (*MoleculeType)[mtype].nBTypes = 0;
          (*MoleculeType)[mtype].BType = calloc(1, sizeof(int));
          for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
            int btype = (*Bead)[(*Molecule)[i].Bead[j]].Type;
            (*MoleculeType)[mtype].Bead[j] = btype;
            exists = false; // recycling the bool to check if btype is already in BType[]
            for (int k = 0; k < (*MoleculeType)[mtype].nBTypes; k++) {
              if ((*MoleculeType)[mtype].BType[k] == btype) {
                exists = true;
              }
            }
            if (!exists) { // recycled exists
              int types = (*MoleculeType)[mtype].nBTypes;
              (*MoleculeType)[mtype].nBTypes++;
              (*MoleculeType)[mtype].BType = realloc((*MoleculeType)[mtype].BType, (types+1)*sizeof(int));
              (*MoleculeType)[mtype].BType[types] = btype;
            }
          }
          // copy bonds
          (*MoleculeType)[mtype].nBonds = bonds_per_mol[i];
          (*MoleculeType)[mtype].Bond = calloc((*MoleculeType)[mtype].nBonds, sizeof(int *));
          for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
            (*MoleculeType)[mtype].Bond[j] = calloc(3, sizeof(int));
            (*MoleculeType)[mtype].Bond[j][0] = mol_bonds[i][j][0];
            (*MoleculeType)[mtype].Bond[j][1] = mol_bonds[i][j][1];
            (*MoleculeType)[mtype].Bond[j][2] = mol_bonds[i][j][2];
          }
          (*MoleculeType)[mtype].nAngles = 0;
          (*MoleculeType)[mtype].Angle = calloc(1, sizeof(int *));
          (*MoleculeType)[mtype].Write = true;
          (*Counts).TypesOfMolecules++;
        } //}}}
      } //}}}
      // free helper arrays //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          free(mol_bonds[i][j]);
        }
        free(mol_bonds[i]);
      }
      free(mol_bonds);
      free(bonds_per_mol);
      free(diff_conn); //}}}
    } //}}}
    // angles section //{{{
    if (words > 0 && strcmp(split[0], "Angles") == 0) {
      // allocate helper arrays to hold angle info //{{{
      // number of angles in each molecule
      int *angles_per_mol = calloc((*Counts).Molecules, sizeof(int));
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] & [i][][2] ... ids of connected beads in molecule 'i'
      // [i][][3] ... angle type
      int ***mol_angles = calloc((*Counts).Molecules, sizeof(int **));
      for (int i = 0; i < (*Counts).Molecules; i++) {
         mol_angles[i] = calloc(1, sizeof(int *));
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof(line), fr);
      // read all angles //{{{
      for (int i = 0; i < *angles; i++) {
        fgets(line, sizeof(line), fr);
        words = SplitLine(split, line, delim);
        // format of each line: <angle id> <angle type> <bead1> <bead2> <bead3>
        // Error - incorrect format //{{{
        if (words < 5 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3]) || !IsInteger(split[4])) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each 'Angles' line must be <angle id> <angle type> <bead1d> <bead2> <bead3>\n");
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        int bead3 = atoi(split[4]) - 1;
        // assign molecule to the bond
        int mol = (*Bead)[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != (*Bead)[bead2].Molecule || mol != (*Bead)[bead3].Molecule) {
          RedText(STDERR_FILENO);
          fprintf(stderr, "\nError ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - atoms in ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "angle %d", atoi(split[0]));
          RedText(STDERR_FILENO);
          fprintf(stderr, " are in different molecules\n\n");
          ResetColour(STDERR_FILENO);
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        angles_per_mol[mol]++;
        int angle = angles_per_mol[mol];
        // add angle beads to the molecule they belong to
        mol_angles[mol] = realloc(mol_angles[mol], angle*sizeof(int *));
        mol_angles[mol][angle-1] = calloc(4, sizeof(int));
        mol_angles[mol][angle-1][0] = bead1;
        mol_angles[mol][angle-1][1] = bead2;
        mol_angles[mol][angle-1][2] = bead3;
        mol_angles[mol][angle-1][3] = type;
      } //}}}
      // sort angles according to the id of the first bead in a bond //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        SortAngles(mol_angles[i], angles_per_mol[i]);
      } //}}}
      // minimize mol_angles based on lowest id in each molecule //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int lowest = (*Counts).Beads; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if ((*Molecule)[i].Bead[j] < lowest) {
            lowest = (*Molecule)[i].Bead[j];
          }
        }
        for (int j = 0; j < angles_per_mol[i]; j++) {
          mol_angles[i][j][0] -= lowest;
          mol_angles[i][j][1] -= lowest;
          mol_angles[i][j][2] -= lowest;
        }
      } //}}}
      // add angles to existing molecule types (or create a new one differing only in angles)
      // ...'basic' molecule types already generated in bonds section
      int *extra = calloc((*Counts).TypesOfMolecules, sizeof(int)); // number of molecule types differing in angles
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int mtype = (*Molecule)[i].Type;
        if ((*MoleculeType)[mtype].nAngles == 0) { // add angles if there are no angles in the molecule type //{{{
          (*MoleculeType)[mtype].nAngles = angles_per_mol[i];
          (*MoleculeType)[mtype].Angle = realloc((*MoleculeType)[mtype].Angle, (*MoleculeType)[mtype].nAngles*sizeof(int *));
          for (int j = 0; j < (*MoleculeType)[mtype].nAngles; j++) {
            (*MoleculeType)[mtype].Angle[j] = calloc(4, sizeof(int));
            (*MoleculeType)[mtype].Angle[j][0] = mol_angles[i][j][0];
            (*MoleculeType)[mtype].Angle[j][1] = mol_angles[i][j][1];
            (*MoleculeType)[mtype].Angle[j][2] = mol_angles[i][j][2];
            (*MoleculeType)[mtype].Angle[j][3] = mol_angles[i][j][3];
          } //}}}
        } else { // if angles already present, check whether they're the same //{{{
          bool exists = true;
          // same number of angles? //{{{
          if ((*MoleculeType)[mtype].nAngles != angles_per_mol[i]) {
            exists = false;
          } //}}}
          // same angles as in mtype? //{{{
          for (int j = 0; j < angles_per_mol[i] && j < (*MoleculeType)[mtype].nAngles; j++) {
            if ((*MoleculeType)[mtype].Angle[j][0] != mol_angles[i][j][0] ||
                (*MoleculeType)[mtype].Angle[j][1] != mol_angles[i][j][1] ||
                (*MoleculeType)[mtype].Angle[j][2] != mol_angles[i][j][2] ||
                (*MoleculeType)[mtype].Angle[j][3] != mol_angles[i][j][3] ) {
              exists = false;
            }
          } //}}}
          // check against angle-generated molecule types //{{{
          // If its angles aren't the same as in mtype, check against other
          // types (i.e., against newly generated thanks to different angles)
          if (!exists) {
            for (int j = (mtype+1); j < (*Counts).TypesOfMolecules; j++) {
              if ((*MoleculeType)[j].nAngles == angles_per_mol[i]) {
                int count = 0;
                for (int k = 0; k < (*MoleculeType)[j].nAngles; k++) {
                  if ((*MoleculeType)[j].Angle[k][0] == mol_angles[i][k][0] &&
                      (*MoleculeType)[j].Angle[k][1] == mol_angles[i][k][1] &&
                      (*MoleculeType)[j].Angle[k][2] == mol_angles[i][k][2] &&
                      (*MoleculeType)[j].Angle[k][3] == mol_angles[i][k][3] ) {
                    count++;
                  }
                }
                if (count == angles_per_mol[i]) {
                  exists = true;
                  (*MoleculeType)[j].Number++;
                  (*MoleculeType)[mtype].Number--;
                  (*Molecule)[i].Type = j;
                  break;
                }
              }
            }
          } //}}}
          // create a new molecule type //{{{
          if (!exists) {
            int new = (*Counts).TypesOfMolecules;
            extra[mtype]++;
            (*Molecule)[i].Type = new;
            (*Counts).TypesOfMolecules++;
            *MoleculeType = realloc(*MoleculeType, (new+1)*sizeof(struct MoleculeType));
            char name[20];
            sprintf(name, "%s-a%d", (*MoleculeType)[mtype].Name, extra[mtype]);
            strcpy((*MoleculeType)[new].Name, name);
            (*MoleculeType)[new].Number = 1;
            (*MoleculeType)[mtype].Number--;
            (*MoleculeType)[new].nBeads = (*MoleculeType)[mtype].nBeads;
            (*MoleculeType)[new].Bead = calloc((*MoleculeType)[new].nBeads, sizeof(int));
            for (int j = 0; j < (*MoleculeType)[new].nBeads; j++) {
              (*MoleculeType)[new].Bead[j] = (*MoleculeType)[mtype].Bead[j];
            }
            (*MoleculeType)[new].nBonds = (*MoleculeType)[mtype].nBonds;
            (*MoleculeType)[new].Bond = calloc((*MoleculeType)[new].nBonds, sizeof(int *));
            for (int j = 0; j < (*MoleculeType)[new].nBonds; j++) {
              (*MoleculeType)[new].Bond[j] = calloc(3, sizeof(int));
              (*MoleculeType)[new].Bond[j][0] = (*MoleculeType)[mtype].Bond[j][0];
              (*MoleculeType)[new].Bond[j][1] = (*MoleculeType)[mtype].Bond[j][1];
              (*MoleculeType)[new].Bond[j][2] = (*MoleculeType)[mtype].Bond[j][2];
            }
            (*MoleculeType)[new].nAngles = (*MoleculeType)[mtype].nAngles;
            (*MoleculeType)[new].Angle = calloc((*MoleculeType)[new].nAngles, sizeof(int *));
            for (int j = 0; j < (*MoleculeType)[new].nAngles; j++) {
              (*MoleculeType)[new].Angle[j] = calloc(4, sizeof(int));
              (*MoleculeType)[new].Angle[j][0] = mol_angles[i][j][0];
              (*MoleculeType)[new].Angle[j][1] = mol_angles[i][j][1];
              (*MoleculeType)[new].Angle[j][2] = mol_angles[i][j][2];
              (*MoleculeType)[new].Angle[j][3] = mol_angles[i][j][3];
            }
            (*MoleculeType)[new].nBTypes = (*MoleculeType)[mtype].nBTypes;
            (*MoleculeType)[new].BType = calloc((*MoleculeType)[new].nBTypes, sizeof(int));
            for (int j = 0; j < (*MoleculeType)[new].nBTypes; j++) {
              (*MoleculeType)[new].BType[j] = (*MoleculeType)[mtype].BType[j];
            }
            (*MoleculeType)[new].Mass = (*MoleculeType)[mtype].Mass;
            (*MoleculeType)[new].Write = true;
          } //}}}
        } //}}}
      }
      // free helper arrays //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        for (int j = 0; j < angles_per_mol[i]; j++) {
          free(mol_angles[i][j]);
        }
        free(mol_angles[i]);
      }
      free(mol_angles);
      free(angles_per_mol);
      free(extra); //}}}
    } //}}}
    // read and split next line
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, delim);
  }
  free(mols); //}}}

  // calculate molecular mass //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      (*MoleculeType)[i].Mass += (*BeadType)[(*MoleculeType)[i].Bead[j]].Mass;
    }
  } //}}}

  WarnElNeutrality(*Counts, *BeadType, data_file);

  fclose(fr);
} //}}}
