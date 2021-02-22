#include "Read.h"

void SkipStructVtf(FILE *vtf, char *name_vtf) { //{{{
  char line[LINE];
  while (fgets(line, sizeof(line), vtf)) {
    char split[30][100];
    int words = SplitLine(split, line, " \t");
    // only certain keywords can be present before the first coordinate block
    // 1) t(imestep) ... starting the coordinate block
    // 2) t(imestep) i(ndexed)/o(ordered) ... starting the coordinate block
    // 3) i(ndexed)/o(ordered) ... starting the coordinate block
    // 4) a(tom) ... in case of vtf file
    // 5) b(ond) ... in case of vtf file
    // 6) blank line
    // 7) # ... comment
    // 8) pbc line
    if (!(split[0][0] == 't' && words == 1) && // 1)
        !(split[0][0] == 't' && words > 1 &&
          (split[1][0] == 'o' || split[1][0] =='i')) && // 2)
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
      fprintf(stderr, " - unrecognised line in structure section");
      fprintf(stderr, "or the first timestep's preamble\n");
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
    // restore file pointer to the last bond line when timestep line encountered
    if ((split[0][0] == 't' && words == 1) || // 1)
        (split[0][0] == 't' && words > 1 &&
         (split[1][0] == 'o' || split[1][0] =='i')) || // 2)
        (split[0][0] == 'o' || split[0][0] == 'i')) { // 3)
      fsetpos(vtf, &position);
      break;
    }
  };
} //}}}

// GetPBC() - TO DELETE //{{{
/*
 * Function to get box dimensions from the provided coordinate file.
 */
VECTOR GetPBC(FILE *vcf, char *vcf_file) {

  VECTOR BoxLength;

  char line[LINE];
  while (fgets(line, sizeof(line), vcf)) {
    char split[30][100];
    int words = SplitLine(split, line, " \t");

    if (strcasecmp(split[0], "pbc") == 0) {
      BoxLength.x = atof(split[1]);
      BoxLength.y = atof(split[2]);
      BoxLength.z = atof(split[3]);
      break;
    // only certain keywords can be present before the first coordinate block
    // 1) t(imestep) ... starting the coordinate block
    // 2) t(imestep) i(ndexed)/o(ordered) ... starting the coordinate block
    // 3) i(ndexed)/o(ordered) ... starting the coordinate block
    // 4) a(tom) ... in case of vtf file
    // 5) b(ond) ... in case of vtf file
    // 6) blank line
    // 7) # ... comment
    // 8) pbc line
    } else if (!(split[0][0] == 't' && words == 1) && // 1)
               !(split[0][0] == 't' && words > 1 &&
                 (split[1][0] == 'o' || split[1][0] =='i')) && // 2)
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

// GetPBC2() - RENAME TO GetPBC() //{{{
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
    char split[30][100];
    int words = SplitLine(split, line, " \t");

    if (strcasecmp(split[0], "pbc") == 0) {
      if (!IsPosDouble(split[1]) ||
          !IsPosDouble(split[2]) ||
          !IsPosDouble(split[3]) ) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", vcf_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - incorrect pbc line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      BoxLength.x = atof(split[1]);
      BoxLength.y = atof(split[2]);
      BoxLength.z = atof(split[3]);
      break;
    // possible keywords (besides pbc) before the first coordinate block
    // 1) t(imestep) ... starting the coordinate block
    // 2) t(imestep) i(ndexed)/o(ordered) ... starting the coordinate block
    // 3) i(ndexed)/o(ordered) ... starting the coordinate block
    // 4) a(tom) ... in case of vtf file
    // 5) b(ond) ... in case of vtf file
    // 6) blank line
    // 7) # ... comment
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
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
  };

  fclose(vcf);

  return BoxLength;
} //}}}

// ReadAggCommand() //{{{
/*
 * Function to read the Aggregate command from an agg file.
 */
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
  char line[LINE], split[30][100];
  fgets(line, sizeof(line), agg);
  int words = SplitLine(split, line, " \t");
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
  // read <distance> argument from Aggregates command //{{{
  if (!IsPosDouble(split[2])) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - <distance> in the Aggregate command");
    fprintf(stderr, " must be non-negative real number\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  }
  *distance = atof(split[2]); //}}}
  // read <contacts> argument from Aggregates command //{{{
  if (!IsInteger(split[3])) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - <contacts> in the Aggregate command");
    fprintf(stderr, "must be a non-negative integer\n");
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
    fprintf(stderr, "         Mismatch between beads present in both files");
    fprintf(stderr, "can lead to undefined behaviour.\n\n");
    ResetColour(STDERR_FILENO);
  } //}}}
  // read <type names> from Aggregates command //{{{
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType(split[i], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_agg);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - non-existent bead name (");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", split[i]);
      RedText(STDERR_FILENO);
      fprintf(stderr, ") in Aggregate command\n");
      ResetColour(STDERR_FILENO);
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}
  fclose(agg);
} //}}}

// CountVtfStructLines() //{{{
/**
 * Function to count lines in the structure part of a vtf file (i.e., it
 * returns the number of the last line containing atom/bond keyword).
 */
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
    int count_lines = 0, last_line = -1;
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
    if (last_line == -1) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "CountVtfStructLines()");
      RedText(STDERR_FILENO);
      fprintf(stderr, " - something weird\n\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
    return last_line;
  } else {
    return -1; // if not vtf, the number doesn't matter
  }
} //}}}

// SkipVtfStructure() //{{{
// TODO: move to General.c; change to while (getc(..)!='\n');
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
  // there are several possibilities how the timestep line can look
  // 1) 't[imestep]'
  if (words == 1 && strncasecmp(split[0], "timestep", 1) == 0) {
    return true;
  }
  // 2) 't[imestep] i[indexed]/o[rdered] any extra stuff'
  if (words > 1 && strncasecmp(split[0], "timestep", 1) == 0 &&
      (strncasecmp(split[1], "indexed", 1) == 0 ||
       strncasecmp(split[1], "ordered", 1) == 0)) {
    return true;
  }
  // 3) 'i[indexed]/o[rdered] any extra stuff'
  if (words > 0 &&
      (strncasecmp(split[0], "indexed", 1) == 0 ||
       strncasecmp(split[0], "ordered", 1) == 0)) {
    return true;
  }
  return false; // not a timestep line
} //}}}

// CheckVtfAtomLine() //{{{
/**
 * Function to check if the provided line is a vtf atom line.
 */
bool CheckVtfAtomLine(int words, char split[30][100], char *error) {
  // is there even number of strings?
  if ((words%2) != 0) {
    strcpy(error, "odd number of strings in an atom line");
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
        if (strncmp(split[j], "resname", 3) == 0 && j != i) {
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
  // if mass, charge, or radius are present, they must be followed by <float>
  for (int i = 2; i < words; i+=2) {
    if ((strcmp(split[i], "charge") == 0 || split[i][0] == 'q') &&
        !IsDouble(split[i+1])) {
      strcpy(error, "'charge|q' must be followed by <float>");
      return false;
    } else if (split[i][0] == 'r' && strncmp(split[i], "res", 3) != 0 &&
               (!IsPosDouble(split[i+1]) || atof(split[i+1]) <= 0)) {
      strcpy(error, "'radius' must be followed by positive <float>");
      return false;
    } else if (split[i][0] == 'm' &&
               (!IsPosDouble(split[i+1]) || atof(split[i+1]) <= 0)) {
      strcpy(error, "'mass' must be followed by positive <float>");
      return false;
    }
  }
  return true;
} //}}}

// CheckVtfBondLine() //{{{
/**
 * Function to check if the provided line is a vtf bond line.
 */
bool CheckVtfBondLine(int words, char split[30][100], char *error) {
  // the line must be: b[ond] <id>:<id> (with possible blanks around ':')
  switch(words) {
    case 1: // there must be at least 1 more string after 'bond'
      strcpy(error, "missing a string after 'bond'");
      return false;
    case 2: // if there's only 1 string after 'bond', it must be '<int>:<int>'
      for (int i = 1; i < (strlen(split[1])-2); i++) {
        if ((split[1][i] == ':' && split[1][i+1] == ':') ||
            split[1][i] == ',' || split[1][i+1] == ',') {
          strcpy(error, "for now, only bond lines in the format <int>:<int> \
are allowed (i.e., no comma separated list or 'two-colon' string of bonds)");
          return false;
        }
      }
      char index[30][100], string[100];
      strcpy(string, split[1]);
      // split '<int>:<int>' into two strings and test those
      int strings = SplitLine(index, string, ":");
      if (strings < 2 || !IsInteger(index[0]) || !IsInteger(index[1]) ||
          split[1][0] == ':' || split[1][strlen(split[1])-1] == ':') {
        strcpy(error, "missing two colon separated integers in a bond line");
        return false;
      }
      return true;
    case 3: // if there are 2 strings after 'bond', they must be '<int>: <int>'
      if (split[1][strlen(split[1])-1] == ':' && IsInteger(split[2])) {
        if (split[1][strlen(split[1])-2] == ':') { // two colons aren't allowed
          strcpy(error, "for now, only bond lines in the format <int>:<int> \
are allowed (i.e., no comma separated list or 'two-colon' string of bonds)");
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
        strcpy(error, "unrecognised bond line (note that whitespace \
is not allowed before the colon; i.e., '<int> : <int>' is not allowed)");
        return false;
      }
      return true;
    default: // more then 2 strings after 'bond' isn't allowed (for now)
      strcpy(error, "unrecognised bond line (note that for now, only \
the format <int>:<int> is allowed; i.e., no comma separated list \
or 'two-colon' string of bonds)");
      return false;
  }
} //}}}

// CheckVtfCoordinateLine() //{{{
/**
 * Function to check if the provided line is a vtf coordinate line.
 */
bool CheckVtfCoordinateLine(int words, char split[30][100],
                            char *error, bool indexed) {
// the line must be: <x> <y> <z>, preceded by <int> for indexed timestep
  // not enough strings
  if ((words < 3 && !indexed) ||
      (words < 4 && indexed)) {
    strcpy(error, "wrong coordinate line");
    return false;
  }
  // incorrect string types
  if (indexed && (!IsInteger(split[0]) ||
                  !IsDouble(split[1]) ||
                  !IsDouble(split[2]) ||
                  !IsDouble(split[3]))) {
    strcpy(error, "wrong coordinate line for an indexed timestep");
    return false;
  }
  if (!indexed && (!IsDouble(split[0]) ||
                   !IsDouble(split[1]) ||
                   !IsDouble(split[2]))) {
    strcpy(error, "wrong coordinate line for an ordered timestep");
    return false;
  }
  // correct line
  return true;
} //}}}

// CopyMoleculeType() //{{{
/**
 * Function to copy a MOLECULETYPE struct into a new one. Requires that the
 * output struct is not alloc'd.
 */
void CopyMoleculeType(int number_of_types,
                      MOLECULETYPE **mt_out, MOLECULETYPE *mt_in) {
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
/**
 * Function to copy a MOLECULE struct into a new one. Requires that the output
 * struct is not alloc'd.
 */
void CopyMolecule(int number_of_molecules, MOLECULETYPE *mt,
                  MOLECULE **m_out, MOLECULE *m_in) {
  *m_out = calloc(number_of_molecules, sizeof(MOLECULE));
  for (int i = 0; i < number_of_molecules; i++) {
    (*m_out)[i].Type = m_in[i].Type;
    (*m_out)[i].Aggregate = m_in[i].Aggregate;
    int mtype = m_in[i].Type;
    (*m_out)[i].Bead = calloc(mt[mtype].nBeads, sizeof(int));
    for (int j = 0; j < mt[mtype].nBeads; j++) {
      (*m_out)[i].Bead[j] = m_in[i].Bead[j];
    }
  }
} //}}}

// CopyBead() //{{{
/**
 * Function to copy a BEAD struct into a new one. Requires that the new
 * BEAD is properly alloc'd.
 * TODO: when (n)Aggregate(s) stuff is removed, add bool whether output struct
 * is alloc'd (akin to CopyBeadType).
 */
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
 * Function to add a new bead type to a BEADTYPE struct (and increment the
 * number of bead types).
 */
void NewBeadType(BEADTYPE **BeadType, int *number_of_types, char *name,
                 double charge, double mass, double radius) {
  int btype = *number_of_types;
  (*number_of_types)++;
  *BeadType = realloc(*BeadType, (btype+1)*sizeof(BEADTYPE));
  strcpy((*BeadType)[btype].Name, name);
  (*BeadType)[btype].Number = 0;
  (*BeadType)[btype].Charge = charge;
  (*BeadType)[btype].Mass = mass;
  (*BeadType)[btype].Radius = radius;
}; //}}}

// NewMolType() //{{{
/*
 * Function to create a new molecule type in a MOLECULETYPE struct.
 * TODO: not used anywhere
 */
void NewMolType(char *name, int *number_of_types,
                MOLECULETYPE **MoleculeType, char *vsf_file) {
  int mtype = (*number_of_types)++;
  *MoleculeType = realloc(*MoleculeType,
                          (*number_of_types)*sizeof(MOLECULETYPE));
  // copy new name to MoleculeType[].Name
  strncpy((*MoleculeType)[mtype].Name, name, MOL_NAME);
  // initialize struct members
  (*MoleculeType)[mtype].Number = 0;
  (*MoleculeType)[mtype].nBonds = 0;
  (*MoleculeType)[mtype].nBeads = 0;
  (*MoleculeType)[mtype].nBTypes = 0;
}; //}}}

// FillMolMass //{{{
/*
 * Function to calculate mass of all molecules. If at least one bead has
 * undefined mass, the mass of the molecule is also undefined
 */
void FillMolMass(int number_of_types,
                 BEADTYPE *BeadType, MOLECULETYPE **MoleculeType) {
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
      YellowText(STDERR_FILENO);
      fprintf(stderr, "\nWarning: molecule type ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%s", (*MoleculeType)[i].Name);
      YellowText(STDERR_FILENO);
      fprintf(stderr, " has undefined mass\n");
      ResetColour(STDERR_FILENO);
    }
  }
} //}}}

// FillMolCharge //{{{
/*
 * Function to calculate charge of all molecules. If at least one bead has
 * undefined charge, the charge of the molecule is also undefined
 */
void FillMolCharge(int number_of_types, BEADTYPE *BeadType,
                   MOLECULETYPE **MoleculeType) {
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
      YellowText(STDERR_FILENO);
      fprintf(stderr, "\nWarning: molecule type ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%s", (*MoleculeType)[i].Name);
      YellowText(STDERR_FILENO);
      fprintf(stderr, " has undefined charge\n");
      ResetColour(STDERR_FILENO);
    }
  }
} //}}}

// FillMolType //{{{
/*
 * Function to fill BType array and mass and charge for each molecule type
 */
void FillMolType(int number_of_types, BEADTYPE *BeadType,
                 MOLECULETYPE **MoleculeType) {
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
  int lines = ReadTimestepPreamble(&indexed, vcf_file, vcf, &stuff, true);
  free(stuff);
  for (int i = 0; i < lines; i++) {
    char line[LINE];
    fgets(line, sizeof line, vcf);
  } //}}}
  // use Bead[].Flag to determine which beads are in the timestep //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    (*Bead)[i].Flag = false;
  } //}}}
  // count beads in the timestep & save their ids //{{{
  lines = 0;
  while (true) {
    char split[30][100], line[LINE];
    fgets(line, sizeof line, vcf);
    // break loop on the end of the coordinate file
    if (feof(vcf)) {
      break;
    }
    int words = SplitLine(split, line, "\t ");
    // break loop if the line isn't a coordinate line
    char error[LINE];
    if (!CheckVtfCoordinateLine(words, split, error, indexed)) {
      if (words == 0 || split[0][0] == '#' ||
          strcasecmp(split[0], "pbc") == 0) {
        break;
      } else { // error if not the start of the next step's preamble
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", vcf_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - %s\n", error);
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
  } //}}}
  // error - if ordered timestep, all beads must be present //{{{
  if (!indexed && lines != (*Counts).Beads) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", vcf_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - ordered timestep needs all beads in the vcf file\n");
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
  // copy beads present in the vcf file to a new BEAD struct & Index array //{{{
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
  // count numbers of beads of different types present in the vcf file //{{{
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
   * ii) copy the molecule types present in the coordinate file to a new
   *     MOLECULETYPE struct
   * iii) determine what bead types are present in the 'new' system (only parts
   *      of a molecule may be in the vcf)
   *      iii)a. copy only bonds between beads present in the vcf file
   *      iii)b. copy only angles between beads present in the vcf file
   *             TODO: test that it workd
   * iv) determine BTypes, charge & mass of the 'new' molecule types (taking
   *     into account that only parts of the original molecules may be present)
   */
  // i) //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int old_type = (*MoleculeType)[i].BType[j];
      if (FindBeadType2((*BeadType)[old_type].Name,
                        c_new.TypesOfBeads, bt_new) != -1) {
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
      if (FindBeadType2((*BeadType)[old_type].Name,
                        c_new.TypesOfBeads, bt_new) != -1) {
        strcpy(mt_new[count].Name, (*MoleculeType)[i].Name);
        mt_new[count].Number = (*MoleculeType)[i].Number;
        // copy beads present in the vcf to the new MOLECULETYPE.Bead array
        mt_new[count].nBeads = 0;
        mt_new[count].Bead = calloc(1, sizeof(int));
        // go through all original types, testing each one
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          int btype = (*MoleculeType)[i].Bead[k];
          btype = FindBeadType2((*BeadType)[btype].Name,
                                c_new.TypesOfBeads, bt_new);
          if (btype != -1) {
            int id = mt_new[count].nBeads++;
            mt_new[count].Bead = realloc(mt_new[count].Bead,
                                         mt_new[count].nBeads*sizeof(int));
            mt_new[count].Bead[id] = btype;
          }
        }
        // iii)a. //{{{
        if ((*MoleculeType)[i].nBonds > 0) {
          mt_new[count].nBonds = 0;
          mt_new[count].Bond = calloc(1, sizeof(int *));
          // go through all original types, testing each pair of beads
          for (int k = 0; k < (*MoleculeType)[i].nBonds; k++) {
            int id1 = (*MoleculeType)[i].Bond[k][0];
            int id2 = (*MoleculeType)[i].Bond[k][1];
            int old_type1 = (*MoleculeType)[i].Bead[id1];
            int old_type2 = (*MoleculeType)[i].Bead[id2];
            int btype1 = FindBeadType2((*BeadType)[old_type1].Name,
                                       c_new.TypesOfBeads, bt_new);
            int btype2 = FindBeadType2((*BeadType)[old_type2].Name,
                                       c_new.TypesOfBeads, bt_new);
            if (btype1 != -1 && btype2 != -1) {
              int id = mt_new[count].nBonds++;
              mt_new[count].Bond = realloc(mt_new[count].Bond,
                                           mt_new[count].nBonds*sizeof(int *));
              mt_new[count].Bond[id] = calloc(3,sizeof(int));
              /* these indices correspond to the old MOLECULETYPE that contains
               * all beads - will be adjusted later */
              mt_new[count].Bond[id][0] = (*MoleculeType)[i].Bond[k][0];
              mt_new[count].Bond[id][1] = (*MoleculeType)[i].Bond[k][1];
              mt_new[count].Bond[id][2] = (*MoleculeType)[i].Bond[k][2];
            }
          }
          // if there are no bonds, free the array
          if (mt_new[count].nBonds == 0) {
            free(mt_new[count].Bond);
          }
        } //}}}
        // iii).b //{{{
        if ((*MoleculeType)[i].nAngles > 0) {
          mt_new[count].nAngles = 0;
          mt_new[count].Angle = calloc(1, sizeof(int *));
          // go through all original types, testing each pair of beads
          for (int k = 0; k < (*MoleculeType)[i].nAngles; k++) {
            int id1 = (*MoleculeType)[i].Angle[k][0];
            int id2 = (*MoleculeType)[i].Angle[k][1];
            int id3 = (*MoleculeType)[i].Angle[k][2];
            int old_type1 = (*MoleculeType)[i].Bead[id1];
            int old_type2 = (*MoleculeType)[i].Bead[id2];
            int old_type3 = (*MoleculeType)[i].Bead[id3];
            int btype1 = FindBeadType2((*BeadType)[old_type1].Name,
                                       c_new.TypesOfBeads, bt_new);
            int btype2 = FindBeadType2((*BeadType)[old_type2].Name,
                                       c_new.TypesOfBeads, bt_new);
            int btype3 = FindBeadType2((*BeadType)[old_type3].Name,
                                       c_new.TypesOfBeads, bt_new);
            if (btype1 != -1 && btype2 != -1 && btype3 != -1) {
              int id = mt_new[count].nBonds++;
              mt_new[count].Bond = realloc(mt_new[count].Bond,
                                           mt_new[count].nBonds*sizeof(int *));
              mt_new[count].Bond[id] = calloc(3,sizeof(int));
              /* these indices correspond to the old MOLECULETYPE that contains
               * all beads - will be adjusted later */
              mt_new[count].Angle[id][0] = (*MoleculeType)[i].Angle[k][0];
              mt_new[count].Angle[id][1] = (*MoleculeType)[i].Angle[k][1];
              mt_new[count].Angle[id][2] = (*MoleculeType)[i].Angle[k][2];
              mt_new[count].Angle[id][3] = (*MoleculeType)[i].Angle[k][3];
            }
          }
          // if there are no angles, free the array
          if (mt_new[count].nAngles == 0) {
            free(mt_new[count].Angle);
          }
        } //}}}
        // rebase bonded beads' indices according to the new MOLECULETYPE //{{{
        /*
         * i) count how much to subtract from old indices in the new
         *    MOLECULETYPE's Bond array by going through the old MOLECULETYPE's
         *    index of bead types and testing if they are in the coordinate file
         * ii) subtract the required amount from each index in every bond
         * iii) subtract the required amount from each index in every angle
         */
        int subtract[(*MoleculeType)[i].nBeads];
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          subtract[k] = 0;
        }
        // i)
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          int old_type = (*MoleculeType)[i].Bead[k];
          int btype = FindBeadType2((*BeadType)[old_type].Name,
                                    c_new.TypesOfBeads, bt_new);
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
        }
        // iii)
        for (int k = 0; k < mt_new[count].nAngles; k++) {
          int id1 = mt_new[count].Angle[k][0];
          int id2 = mt_new[count].Angle[k][1];
          int id3 = mt_new[count].Angle[k][2];
          mt_new[count].Angle[k][0] -= subtract[id1];
          mt_new[count].Angle[k][1] -= subtract[id2];
          mt_new[count].Angle[k][2] -= subtract[id3];
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
  // copy molecules present in the vcf file to a new MOLECULE struct //{{{
  MOLECULE *m_new = calloc(c_new.Molecules, sizeof(MOLECULE));
  count = 0; // count molecules
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int old_mtype = (*Molecule)[i].Type;
    int new_mtype = FindMoleculeType2((*MoleculeType)[old_mtype].Name,
                                      c_new.TypesOfMolecules, mt_new);
    if (new_mtype != -1) {
      m_new[count].Type = new_mtype;
      m_new[count].Bead = calloc(mt_new[new_mtype].nBeads, sizeof(int));
      int count2 = 0; // counts beads in count-th molecule
      for (int j = 0; j < (*MoleculeType)[old_mtype].nBeads; j++) {
        int old_btype = (*MoleculeType)[old_mtype].Bead[j];
        if (FindBeadType2((*BeadType)[old_btype].Name,
                          c_new.TypesOfBeads, bt_new) != -1) {
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
  // BEAD struct - realloc and copy; redo Index array
  *Bead = realloc(*Bead, c_new.Beads*sizeof(BEAD));
  for (int i = 0; i < c_new.Beads; i++) {
    (*Bead)[i] = b_new[i];
    (*Index)[(*Bead)[i].Index] = i;
  }
  // BEADTYPE struct - realloc and copy
  *BeadType = realloc(*BeadType, c_new.TypesOfBeads*sizeof(BEADTYPE));
  for (int i = 0; i < c_new.TypesOfBeads; i++) {
    (*BeadType)[i] = bt_new[i];
  }
  // MOLECULETYPE - free memory because CopyMoleculeType calloc's it
  FreeMoleculeType((*Counts).TypesOfMolecules, MoleculeType);
  CopyMoleculeType(c_new.TypesOfMolecules, MoleculeType, mt_new);
  // MOLECULE - free memory because CopyMolecule calloc's it
  FreeMolecule((*Counts).Molecules, Molecule);
  CopyMolecule(c_new.Molecules, *MoleculeType, Molecule, m_new);
  // copy the new Counts struct back to the original one
  *Counts = c_new; //}}}
  // free memory
  FreeSystemInfo(*Counts, &mt_new, &m_new, &bt_new, &b_new, &index_new);
  return indexed;
} //}}}

// ReadVtfStructure() //{{{
// TODO: split into more functions?
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
      default_atom_line = -1, // line number of the first 'atom default' line
      count_bond_lines = 0, // number of bond lines
      count_comment_lines = 0, // number of comments (#) or blank lines
      atom_names = 0, res_names = 0, // number of unieque bead & molecule names
      last_struct = 0; // number of the last structure line (i.e., bond or atom)
  char **atom_name = calloc(1, sizeof(char *)); // names of different atoms
  char **res_name = calloc(1, sizeof(char *)); // names of different molecules
  while(true) {
    char error[LINE] = {'\0'}, // error message; no error for strlen(error) == 0
         split[30][100], line[LINE];
    fgets(line, sizeof line, vsf);
    if (feof(vsf)) { // exit while loop on vsf file finish
      break;
    }
    int words = SplitLine(split, line, "\t ");
    // in case of a vtf file, stop reading at the 'timestep' line
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
          last_struct = count_atom_lines +
                        count_bond_lines +
                        count_comment_lines;
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
            YellowText(STDERR_FILENO);
            fprintf(stderr, "\nWarning: ");
            CyanText(STDERR_FILENO);
            fprintf(stderr, "%s", vsf_file);
            YellowText(STDERR_FILENO);
            fprintf(stderr, " - multiple 'atom default lines'\n");
            fprintf(stderr, "         Using line number ");
            CyanText(STDERR_FILENO);
            fprintf(stderr, "%d\n", default_atom_line+1);
            ResetColour(STDERR_FILENO);
          } else {
            default_atom_line = count_atom_lines +
                                count_comment_lines +
                                count_bond_lines - 1; // so it starts from 0
          }
        } else {
          // check for highest index (i.e., the number of beads in vsf)
          if (atoi(split[1]) >= (*Counts).BeadsInVsf) {
            (*Counts).BeadsInVsf = atoi(split[1]) + 1; // vtf indices start at 0
          }
        }
        // save bead & molecule names (if not saved already)
        for (int i = 0; i < words; i+= 2) {
          if (split[i][0] == 'n') { // atom name
            // is the name saved already?
            int in = -1;
            for (int j = 0; j < atom_names; j++) {
              if (strncmp(split[i+1], atom_name[j], BEAD_NAME) == 0) {
                in = j;
                break;
              }
            }
            // add new name if not saved yet
            if (in == -1) {
              atom_name = realloc(atom_name, (atom_names+1)*sizeof(char *));
              atom_name[atom_names] = calloc(BEAD_NAME+1, sizeof(char));
              strncpy(atom_name[atom_names], split[i+1], BEAD_NAME);
              atom_names++;
            }
          } else if (strncmp(split[i], "res", 3) == 0 &&
                     strcmp(split[i], "resid") != 0) { // molecule name
            // is the name saved already?
            int in = -1;
            for (int j = 0; j < res_names; j++) {
              if (strncmp(split[i+1], res_name[j], MOL_NAME) == 0) {
                in = j;
                break;
              }
            }
            // add new name if not saved yet
            if (in == -1) {
              res_name = realloc(res_name, (res_names+1)*sizeof(char *));
              res_name[res_names] = calloc(MOL_NAME+1, sizeof(char));
              strncpy(res_name[res_names], split[i+1], MOL_NAME);
              res_names++;
            }
          }
        }
        break; //}}}
      case 'b': // bond line //{{{
        if (CheckVtfBondLine(words, split, error)) {
          count_bond_lines++;
          last_struct = count_atom_lines +
                        count_bond_lines +
                        count_comment_lines;
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
    fprintf(stderr, "       (i.e., when 'atom default' is omitted,");
    fprintf(stderr, " there must be a line for each atom)\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // total number of lines in the structure section (file type dependent) //{{{
  int total_lines = last_struct;
  if (!vtf) {
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
    // array to connect atom ids in vsf with 'struct atom' numbering
    int *atom_id = calloc((*Counts).BeadsInVsf, sizeof(int));
    for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
      atom_id[i] = -1; // if it stays -1, it's 'atom default'
    }
    // structure to save info from bond lines
    struct bond {
      int index1, index2;
    } *bond = calloc(count_bond_lines, sizeof(struct bond)); //}}}
  // 2) save all vsf lines and count number of molecules //{{{
  fsetpos(vsf, &pos); // restore file pointer to the beginning of vsf_file
  int count_atoms = 0, count_bonds = 0;
  for (int line_no = 0; line_no < total_lines; line_no++) {
    char split[30][100], line[LINE];
    fgets(line, sizeof line, vsf);
    int words = SplitLine(split, line, "\t ");
    // in case of a vtf file, finish reading at the 'timestep' line
    if (vtf && CheckVtfTimestepLine(words, split)) {
      break;
    }
    // skip blank, comment, and pbc lines, (i.e., count them as comment lines)
    if (words == 0 || split[0][0] == '#' || strcasecmp(split[0], "pbc") == 0) {
      continue;
    }
    switch (split[0][0]) {
      case 'a': // save atom line into atom struct //{{{
        // ignore multiple 'atom default' lines (use only the first one)
        if (strcmp(split[1], "default") != 0 || line_no == default_atom_line) {
          int id; // linetype-based bead index
          if (line_no == default_atom_line) { // 'atom default' line
            id = count_atom_lines - 1;
            atom[id].index = -1;
          } else { // 'atom <id>' line
            id = count_atoms++;
            atom_id[atoi(split[1])] = id;
            atom[id].index = atoi(split[1]);
          }
          // go through the line
          for (int i = 0; i < words; i+=2) {
            // if the bead is in a molecule, count it as bonded
            if (strncmp(split[i], "resid", 3) == 0 &&
                strcmp(split[i], "resname") != 0) {
              if (atoi(split[i+1]) > (*Counts).Molecules) {
                (*Counts).Molecules = atoi(split[i+1]); // molecule index
              }
              atom[id].resid = atoi(split[i+1]) - 1; // vsf resid starts with 1
            } else if (split[i][0] == 'n') { // bead name
              for (int j = 0; j < atom_names; j++) {
                if (strncmp(split[i+1], atom_name[j], 16) == 0) {
                  atom[id].name = j;
                  break;
                }
              }
            } else if (strcmp(split[i], "resname") == 0) { // molecule name
              for (int j = 0; j < res_names; j++) {
                if (strncmp(split[i+1], res_name[j], 8) == 0) {
                  atom[id].resname = j;
                  break;
                }
              }
            } else if (strcmp(split[i], "charge") == 0 ||
                       strcmp(split[i], "q") == 0) { // bead charge
              atom[id].charge = atof(split[i+1]);
            } else if (split[i][0] == 'm') { // bead mass
              atom[id].mass = atof(split[i+1]);
            } else if (split[i][0] == 'r' &&
                       strncmp(split[i], "res", 3) != 0) { // bead radius
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
  fclose(vsf); // everything is in memory now
  // error - bond with bead id higher than the number of beads //{{{
  for (int i = 0; i < count_bond_lines; i++) {
    if (bond[i].index1 >= (*Counts).BeadsInVsf ||
        bond[i].index2 >= (*Counts).BeadsInVsf) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", vsf_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - bead index too high in ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "bond %d:%d\n\n", bond[i].index1, bond[i].index2);
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  } //}}}
  // fill atom_id for 'atom default' beads //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if (atom_id[i] == -1) {
      atom_id[i] = count_atom_lines - 1;
    }
  } //}}}
  // 3) fill BEADTYPE struct //{{{
  BEADTYPE *bt = calloc(1, sizeof(BEADTYPE));
  // save default type - if there is one
  if (atom[count_atom_lines-1].index == -1) {
    (*Counts).TypesOfBeads = 1;
    strcpy(bt[0].Name, atom_name[atom[count_atom_lines-1].name]);
    bt[0].Charge = atom[count_atom_lines-1].charge;
    bt[0].Mass = atom[count_atom_lines-1].mass;
    bt[0].Radius = atom[count_atom_lines-1].radius;
  }
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
     *      should be three different types
     * iv) however, only some lines can be ambiguous
     *       atom 0 n x q 1 m 1
     *       atom 0 n x     m 1
     *       atom 0 n x q 0 m 1
     *       atom 0 n x q 0
     *     should be still be three different types (the last two should be
     *     considered the same) TODO: for now, they're counted as four types
     * TODO: test third type against a pair to be joined completely; i.e.,
     * test charge/mass/radius against both; although maybe just somehow
     * check for undefined? If everything is defined, then those that
     * aren't to be merged simply aren't to be merged. So maybe at first
     * just create those that differ, disregarding the undefined-ity and
     * then decide which ones containing some undefined-ity to join...
     */
    // create bead types according to bead properties //{{{
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
      if (btype == -1) { // new bead type?
        NewBeadType(&bt, &(*Counts).TypesOfBeads, atom_name[atom[i].name],
                    atom[i].charge, atom[i].mass, atom[i].radius);
        btype = (*Counts).TypesOfBeads - 1;
      }
      // count number of beads of given type (excluding 'atom default'
      // type, if there is one)
      if (default_atom_line == -1 || atom[i].index != (count_atom_lines-1)) {
        bt[btype].Number++;
      }
    } //}}}
    // count number of beads of default type if 'atom default' is present //{{{
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
     *    are to be merged
     * ii) lines 'atom <id> name x mass 1'
     *           'atom <id> name x'
     *           'atom <id> name x mass 2 radius 1'
     *     are not to be merged, as it's uncertain if the second atom should
     *     have mass 1 or 2
     *
     * helper matrix: merge[i][j]
     * i) i!=j: 0 means do not merge i and j, 1 means merge i and j
     * ii) i=j: 0 means i is to be merged with some other bead type; 1 means
     *     i doesn't need merging (or was already merged)
     */
    // 1) assume all bead types are to be merged //{{{
    // TODO: it should be easier to assume the opposite, shouldn't it?
    bool merge[(*Counts).TypesOfBeads][(*Counts).TypesOfBeads];
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        merge[i][j] = true;
      }
    } //}}}
    // 2) exclude from merging bead type pairs that have: //{{{
    //    i) different names
    //    ii) different and well defined values of charge/radius/mass
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[i].Name, bt[j].Name) != 0) { // i)
          merge[i][j] = false;
        } else if ((bt[i].Charge != bt[j].Charge && // ii)
                    bt[i].Charge != CHARGE &&
                    bt[j].Charge != CHARGE) ||
                   (bt[i].Radius != bt[j].Radius &&
                    bt[i].Radius != RADIUS &&
                    bt[j].Radius != RADIUS) ||
                   (bt[i].Mass != bt[j].Mass &&
                    bt[i].Mass != MASS &&
                    bt[j].Mass != MASS)) {
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
     * i) assuming, the first two types (i and j) share the name with the third
     *    type (k), exclude all three pairs from merging if k isn't to be
     *    merged with either i or j
     * ii) if i and j aren't to be merged but share the same name nothing
     *     of the same name is to be merged
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
    // 5) merge the bead types to be merged //{{{
    /*
     * Example of the merge array (see Examples/TODO: where? for the
     * corresponding system):
     *
     *       0 1 2 3 4 5
     *       A B A A B B
     *      ------------
     * 0 A | 0 0 1 1 0 0
     * 1 B | 0 1 0 0 0 0
     * 2 A | 1 0 0 1 0 0
     * 3 A | 1 0 1 1 0 0
     * 4 B | 0 0 0 0 1 0
     * 5 B | 0 0 0 0 0 1
     *
     *  Here, merge[0][0]=0 means that type 0 (the first A) will be merged
     *  with types that have have merge[0][j]=1, namely types 2 and 3 (the
     *  other types called A). Conversely, merge[1][1]=1 meeans that type 1
     *  (the first B) will not be merged with anything (even thought type 4
     *  and 5 share the name B).
     */
    // copy all bead types to a temporary struct //{{{
    BEADTYPE *temp = calloc((*Counts).TypesOfBeads, sizeof(BEADTYPE));
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      temp[i] = bt[i];
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
    // reorder the types, placing same-named ones next to each other //{{{
    count = 0;
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      // copy all bead types to a temporary struct
      temp[i] = bt[i];
      // was 'i' already copied back to bt?
      temp[i].Use = false;
    }
    // copy the bead types back in a proper order
    count = 0;
    for (int i = 0; i < ((*Counts).TypesOfBeads-1); i++) {
      if (!temp[i].Use) {
        bt[count] = temp[i];
        temp[i].Use = true;
        for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
          if (strcmp(temp[i].Name, temp[j].Name) == 0 && !temp[j].Use) {
            count++;
            bt[count] = temp[j];
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
      bt[btype].Number++;
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
        char name[BEAD_NAME+1];
        // shorten name if necessary
        if (count < 10) {
          strncpy(name, bt[j].Name, BEAD_NAME-2);
        } else if (count < 100) {
          strncpy(name, bt[j].Name, BEAD_NAME-3);
        }
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
        // shorten name if necessary to append '_<int>'
        char name[MOL_NAME+1];
        strcpy(name, mt[j].Name);
        if (count < 10) {
          name[MOL_NAME-2] = '\0';
        } else if (count < 100) {
          name[MOL_NAME-3] = '\0';
        }
        snprintf(mt[j].Name, MOL_NAME+1, "%s_%d", name, count);
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
  FreeBead((*Counts).BeadsInVsf, &bead_all);
  FreeMolecule((*Counts).Molecules, &mol);
  FreeMoleculeType((*Counts).TypesOfMolecules, &mt);
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
} //}}}

// FullVtfRead() //{{{
void FullVtfRead(char *vsf_file, char *vcf_file, bool detailed, bool vtf, bool *indexed, int *struct_lines,
                 VECTOR *BoxLength, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {
  // read the whole structure section
  ReadVtfStructure(vsf_file, detailed, Counts, BeadType, Bead, Index, MoleculeType, Molecule);

  // count structure lines if vtf as the coordinate file (-1 otherwise)
  if (vcf_file[0] != '\0') {
    // get box size
    *BoxLength = GetPBC2(vcf_file);
    // number of structure lines (useful if vcf_file is vtf)
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
  }

  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, vsf_file);

  // allocate memory for aggregates
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(1, sizeof(int));
  }
} //}}}

// ReadTimestepPreamble() //{{{
/**
 * Function to read timestep preamble
 */
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
    strcpy(line2, line); // save line string to strcat to stuff array
    words = SplitLine(split, line, " \t");
    // comment line - copy to stuff array
    if (split[0][0] == '#') {
      SafeStrcat(stuff, line2, LINE);
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
    // continue reading the same line if it's too long
    while (strlen(line2) == (LINE-1)) {
      fgets(line2, sizeof(line2), vcf_file);
      if (feof(vcf_file)) {
        return -1;
      }
      SafeStrcat(stuff, line2, LINE);
      count_lines++;
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
  char line2[LINE];
  fgets(line2, sizeof(line2), vcf_file);
  fsetpos(vcf_file, &position); // restore pointer position
  fgets(line2, sizeof(line2), vcf_file);
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

// SkipVcfCoor() //{{{
/**
 * Function to skip one timestep in a vcf coordinate file.
 */
void SkipVcfCoor(FILE *vcf_file, char *input_coor, COUNTS Counts, char **stuff) {
  bool rubbish; // testing timestep type - not used here
  // read and skip timestep preamble
  int count_lines = ReadTimestepPreamble(&rubbish, input_coor, vcf_file, stuff, false);
  for (int i = 0; i < count_lines; i++) {
    char line[LINE];
    fgets(line, sizeof(line), vcf_file);
  }
  // skip coordinate lines
  for (int i = 0; i < Counts.Beads; i++) {
    char line[LINE];
    fgets(line, sizeof(line), vcf_file);
    // premature end of file?
    if (feof(vcf_file) == EOF) {
      RedText(STDERR_FILENO);
      fprintf(stderr, "\nError: premature end of ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " file\n\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  }
} //}}}

// ReadAggregates() //{{{
/**
 * Function reading information about aggregates from `.agg` file (\ref AggregateFile) generated by Aggregates utility.
 */
void ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts, AGGREGATE **Aggregate,
                    BEADTYPE *BeadType, BEAD **Bead,
                    MOLECULETYPE *MoleculeType, MOLECULE **Molecule, int *Index) {

  char line[LINE], split[30][100];
  // read (Last) Step line
  fgets(line, sizeof(line), fr);
  int words = SplitLine(split, line, " \t");
  // error if the first line is 'L[ast Step]' or isn't 'Step: <int>'//{{{
  if (split[0][0] == 'L') {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: premature end of ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " file");
    ResetColour(STDERR_FILENO);
    exit(1);
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
  } //}}}

  // initialize array of number of aggregates per bead //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].nAggregates = 0;
  } //}}}

  // get number of aggregates //{{{
  fgets(line, sizeof(line), fr);
  words = SplitLine(split, line, " \t");
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
} //}}}

// SkipAgg() //{{{
/**
 * Function to skip one timestep in aggregate file.
 */
void SkipAgg(FILE *agg, char *agg_file) {
  char line[LINE], split[30][100];
  // read (Last) Step line
  fgets(line, sizeof(line), agg);
  int words = SplitLine(split, line, " \t");
  // error if the first line is 'L[ast Step]' or isn't 'Step: <int>'//{{{
  if (split[0][0] == 'L') {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: premature end of ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " file");
    ResetColour(STDERR_FILENO);
    exit(1);
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
  } //}}}

  // get number of aggregates
  fgets(line, sizeof(line), agg);
  words = SplitLine(split, line, " \t");
  // Error - number of aggregates must be <int> //{{{
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
  } //}}}
  int aggs = atoi(split[0]);
  int test;
  for (int i = 0; i < (2*aggs); i++) {
    while ((test = getc(agg)) != '\n') {
      if (test == EOF) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: premature end of ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", agg_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " file");
        ResetColour(STDERR_FILENO);
        exit(1);
      }
    }
  }
  // skip empty line at the end
  for (int i = 0; i < 2; i++) {
    while ((test = getc(agg)) != '\n') {
      if (test == EOF) {
        RedText(STDERR_FILENO);
        fprintf(stderr, "\nError: premature end of ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", agg_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " file");
        ResetColour(STDERR_FILENO);
        exit(1);
      }
    }
  }
} //}}}

// ReadFieldPbc() //{{{
/**
 * Function reading box size from the first line of the FIELD-like file.
 */
bool ReadFieldPbc(char *field, VECTOR *BoxLength) {
  bool pbc = false; // assume the first line doesn't contain box size
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}
  // read first line
  char line[LINE], split[30][100];
  fgets(line, sizeof(line), fr);
  int words = SplitLine(split, line, " \t");
  // box size must be 'pbc <double> <double> <double>', so test:
  // 1) enough strings
  // 2) 'pbc' string
  // 3) positive doubles
  if (words >= 4 && // 1
      strcasecmp(split[0], "pbc") == 0 && // 3) pbc keyword
      IsPosDouble(split[1]) && // 2) first <double>
      IsPosDouble(split[2]) && // 2) second <double>
      IsPosDouble(split[3])) { // 2) third <double>
    (*BoxLength).x = atof(split[1]);
    (*BoxLength).y = atof(split[2]);
    (*BoxLength).z = atof(split[3]);
    pbc = true;
  }
  fclose(fr);
  return pbc;
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

  char line[LINE], split[30][100];

  // read number of bead types //{{{
  bool missing = true; // is 'species' keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    int words = SplitLine(split, line, " \t");
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
    int words = SplitLine(split, line, " \t");

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

    snprintf((*BeadType)[i].Name, BEAD_NAME+1, "%s", split[0]);
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

  char line[LINE], split[30][100];

  // read number of molecule types //{{{
  bool missing = true; // is molecule keyword missing?
  while(fgets(line, sizeof(line), fr)) {
    int words = SplitLine(split, line, " \t");
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
  } //}}}

  // test proper number of 'finish' keywords //{{{
  fpos_t position;
  fgetpos(fr, &position); // save pointer position
  int count = 0;
  while(fgets(line, sizeof(line), fr)) {
    SplitLine(split, line, " \t");
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

  // read info about molecule types
  fgetpos(fr, &position); // save pointer position
  (*Counts).TypesOfBonds = 0;
  (*Counts).TypesOfAngles = 0;
  (*bond_type) = calloc(1, sizeof(struct Params));
  (*angle_type) = calloc(1, sizeof(struct Params));
  // stored bond & angle types - temporary //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // read i-th molecule name
    fgets(line, sizeof(line), fr);
    SplitLine(split, line, " \t");
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
    snprintf((*MoleculeType)[i].Name, MOL_NAME+1, "%s", split[0]);
    // read number of the molecules //{{{
    fgets(line, sizeof(line), fr);
    int words = SplitLine(split, line, " \t");
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
    words = SplitLine(split, line, " \t");
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
      words = SplitLine(split, line, " \t");
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
    words = SplitLine(split, line, " \t");
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
      words = SplitLine(split, line, " \t");
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
    words = SplitLine(split, line, " \t");
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
        words = SplitLine(split, line, " \t");
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
      SplitLine(split, line, " \t");
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
      SplitLine(split, line, " \t");
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
      SplitLine(split, line, " \t");
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
      SplitLine(split, line, " \t");
      for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
        if ((*angle_type)[k].a == atof(split[4]) && (*angle_type)[k].b == atof(split[5])) {
          (*MoleculeType)[i].Angle[j][3] = k;
          break;
        }
      }
    } //}}}
    // skip till 'finish' //{{{
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, " \t");
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
               PARAMS **bond_type, PARAMS **angle_type,
               PARAMS **dihedral_type) {

  // read pbc if required //{{{
  if (BoxLength != NULL && !ReadFieldPbc(field, BoxLength)) {
    RedText(STDERR_FILENO);
    fprintf(stderr, "\nError: ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - first line must start with box size, ");
    fprintf(stderr, "i.e., 'pbc <double> <double> <double>'\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  ReadFieldBeadType(field, Counts, BeadType, Bead);
  ReadFieldMolecules(field, Counts, BeadType, Bead, MoleculeType, Molecule,
                     bond_type, angle_type);
  // allocate Bead[].Aggregate array - needed only to free() //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(10, sizeof(int));
  } //}}}
  FillMolBTypes((*Counts).TypesOfMolecules, MoleculeType);
  // fill Molecule & Bead structs //{{{
  *Molecule = calloc((*Counts).Molecules,sizeof(struct Molecule));
  int count_mol = 0, count_bead = (*Counts).Unbonded;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].Number; j++) {
      (*Molecule)[count_mol].Type = i;
      (*Molecule)[count_mol].Bead =
          malloc((*MoleculeType)[i].nBeads*sizeof(int));
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
  char split[30][100];
  int words;
  do {
    fgets(line, sizeof(line), fr);
    words = SplitLine(split, line, " \t");
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
        words = SplitLine(split, line, " \t");
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
            snprintf((*BeadType)[i].Name, BEAD_NAME+1, "%s", split[3]);
          } else if (strlen(split[2]) > 1) { // comment form '#name'
            for (int j = 0; j < strlen(split[2]); j++) {
              split[2][j] = split[2][j+1];
            }
            strncpy((*BeadType)[i].Name, split[2], BEAD_NAME);
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
        words = SplitLine(split, line, " \t");
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
        words = SplitLine(split, line, " \t");
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
        words = SplitLine(split, line, " \t");
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
      free(mols);
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
        words = SplitLine(split, line, " \t");
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
      int mtype_bead_order = 0;
      for (int i = 0; i < (*Counts).Molecules; i++) {
        // is molecule 'i' of known type?
        bool exists = false;
        // do 'i' and given molecule type share connectivity and bead order?
        bool same_conn = true, same_bead = true;
        // molecule type with which 'i' shares bead order - naming purposes
        int type_bead_order = -1;
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
            if (same_bead && type_bead_order == -1) {
              type_bead_order = j;
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
          if (same_bead && !same_conn) { // same beads - mol<type_bead_order>-b#
            diff_conn[type_bead_order]++;
            // shorten name if necessary to append '-b<int>'
            char name[MOL_NAME+1];
            strcpy(name, (*MoleculeType)[type_bead_order].Name);
            if (diff_conn[type_bead_order] < 10) {
              name[MOL_NAME-3] = '\0';
            } else if (diff_conn[type_bead_order] < 100) {
              name[MOL_NAME-4] = '\0';
            }
            snprintf((*MoleculeType)[mtype].Name, MOL_NAME+1, "%s-b%d", name, diff_conn[type_bead_order]);
          } else { // same connectivity or both different - mol#
            mtype_bead_order++;
            snprintf((*MoleculeType)[mtype].Name, MOL_NAME+1, "mol%d", mtype_bead_order);
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
        words = SplitLine(split, line, " \t");
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
            // shorten name if necessary to append '-a<int>'
            char name[MOL_NAME+1];
            strcpy(name, (*MoleculeType)[mtype].Name);
            if (extra[mtype] < 10) {
              name[MOL_NAME-3] = '\0';
            } else if (extra[mtype] < 100) {
              name[MOL_NAME-4] = '\0';
            }
            snprintf((*MoleculeType)[mtype].Name, MOL_NAME+1, "%s-a%d", name, extra[mtype]);
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
    words = SplitLine(split, line, " \t");
  } //}}}
  free(mols);

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

// SkipCoorSteps() { //{{{
/**
 * Function to skip timesteps (coordinate file only) from the beginning (-st
 * option).
 */
int SkipCoorSteps(FILE *vcf, char *input_coor, COUNTS Counts, int start, bool silent) {
  int test;
  int count = 0;
  char *stuff = calloc(LINE, sizeof(char)); // just for SkipVcfCoor
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}
    SkipVcfCoor(vcf, input_coor, Counts, &stuff);
  }
  // print starting step? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Starting step: %d\n", start);
  } //}}}
  // is the vcf file continuing?
  if (ErrorDiscard(start, count, input_coor, vcf)) {
    exit(1);
  }
  free(stuff);
  return count;
} //}}}

// SkipCoorAggSteps() { //{{{
/**
 * Function to skip timesteps (coordinate and aggregate files) from the
 * beginning (-st option).
 */
int SkipCoorAggSteps(FILE *vcf, char *input_coor, FILE *agg, char *input_agg,
                     COUNTS Counts, int start, bool silent) {
  int test, count = 0;
  char *stuff = calloc(LINE, sizeof(char)); // just for SkipVcfCoor
  // TODO: use feof()
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}
    SkipAgg(agg, input_agg);
    SkipVcfCoor(vcf, input_coor, Counts, &stuff);
  }
  // print starting step? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Starting step: %d\n", start);
  } //}}}
  // is the vcf file continuing?
  if (ErrorDiscard(start, count, input_coor, vcf)) {
    exit(1);
  }
  return count;
} //}}}
