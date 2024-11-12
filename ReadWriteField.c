#include "ReadWriteField.h"
#include "Errors.h"
#include "General.h"

static void FieldReadSpecies(const char *file, SYSTEM *System);
static void FieldReadMolecules(const char *file, SYSTEM *System);
static void SkipTilKeyword(FILE *f, const char *file,
                           const char *keyword, int *line_count);
static bool GetBeadIds(const int num, long id[num], MOLECULETYPE mt_i);
static PARAMS GetParams(const char *file, const int num,
                        const MOLECULETYPE mt, bool *warned);
static bool ReadStuff(const char *file, FILE *fr, int *line_count,
                      const int type, const int mtype, SYSTEM *System);

// GetBeadIds() //{{{
static bool GetBeadIds(const int num, long id[num], MOLECULETYPE mt_i) {
  bool err = false;
  for (int aa = 0; aa < (num - 1); aa++) {
    if (words < num ||
        !IsNaturalNumber(split[aa+1], &id[aa])) {
      err = true;
      break;
    }
  }
  if (err) {
    snprintf(ERROR_MSG, LINE, "incorrect line for molecule %s%s",
             ErrYellow(), mt_i.Name);
    return false;
  }
  // error - bead index is too high
  err = false;
  for (int aa = 0; aa < (num - 1); aa++) {
    if (id[aa] > mt_i.nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("bead index in a improper is too high");
    snprintf(ERROR_MSG, LINE, "bead index in molecule %s%s is too high",
             ErrYellow(), mt_i.Name);
    return false;
  }
  return true;
} //}}}
// GetParams() //{{{
static PARAMS GetParams(const char *file, const int num,
                        const MOLECULETYPE mt, bool *warned) {
  PARAMS values = InitParams;
  for (int k = num; k < (num + 3); k++) {
    double *ptr = NULL;
    if (k == num) {
      ptr = &values.a;
    } else if (k == (num + 1)) {
      ptr = &values.b;
    } else {
      ptr = &values.c;
    }
    if (words > k) {
      if (strcmp(split[k], "???") == 0 && !(*warned)) {
        snprintf(ERROR_MSG, LINE, "undefined parameter (\'???\') "
                 "in molecule %s%s", ErrYellow(), mt.Name);
        PrintWarnFile(file, "\0", "\0");
        *warned = true;
      } else if (!IsRealNumber(split[k], ptr)) {
        break;
      }
    }
  }
  return values;
} //}}}
// ReadStuff() //{{{
static bool ReadStuff(const char *file, FILE *fr, int *line_count,
                      const int type, const int mtype, SYSTEM *System) {
  // assign arrays to pointers based on bond/angle/dihedral/improper //{{{
  PARAMS **StuffType;
  COUNT *Count = &System->Count;
  int *CountStuff;
  int *n_stuff;
  MOLECULETYPE *mt = &System->MoleculeType[mtype];
  int (**Stuff)[5];
  int num = 0;
  if (type == 0) { // bonds
    num = 3;
    n_stuff = &mt->nBonds;
    StuffType = &System->BondType;
    CountStuff = &Count->BondType;
    Stuff = &mt->Bond;
  } else if (type == 1) { // angles
    num = 4;
    n_stuff = &mt->nAngles;
    StuffType = &System->AngleType;
    CountStuff = &Count->AngleType;
    Stuff = &mt->Angle;
  } else if (type == 2) { // dihedreals
    num = 5;
    n_stuff = &mt->nDihedrals;
    StuffType = &System->DihedralType;
    CountStuff = &Count->DihedralType;
    Stuff = &mt->Dihedral;
  } else if (type == 3) { // impropers
    num = 5;
    n_stuff = &mt->nImpropers;
    StuffType = &System->ImproperType;
    CountStuff = &Count->ImproperType;
    Stuff = &mt->Improper;
  } else { // error
    err_msg("ReadStuff(): type must be 0 to 3");
    PrintError();
    exit(1);
  } //}}}
  *Stuff = malloc(sizeof **Stuff * *n_stuff);
  bool warned = false;
  for (int j = 0; j < *n_stuff; j++) {
    (*line_count)++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Molecules' section");
      exit(1);
    } //}}}
    long beads[num];
    if (!GetBeadIds(num, beads, *mt)) {
      return false;
    }
    // assign bead ids to proper 'stuff'
    for (int aa = 0; aa < (num - 1); aa++) {
      (*Stuff)[j][aa] = beads[aa] - 1;
    }
    PARAMS values = GetParams(file, num, *mt, &warned);
    // assign stuff type
    int n_type = -1;
    if (values.a != 0 || values.b != 0 || values.c != 0) {
      // find if this bond type already exists
      for (int k = 0; k < *CountStuff; k++) {
        if ((*StuffType)[k].a == values.a &&
            (*StuffType)[k].b == values.b &&
            (*StuffType)[k].c == values.c) {
          n_type = k;
          break;
        }
      }
      // create a new bond type if necessary
      if (n_type == -1) {
        n_type = *CountStuff;
        (*CountStuff)++;
        *StuffType = s_realloc(*StuffType, *CountStuff * sizeof **StuffType);
        (*StuffType)[n_type] = values;
      }
    }
    (*Stuff)[j][num-1] = n_type;
  }
  return true;
} //}}}

// static void SkipTilKeyword() //{{{
static void SkipTilKeyword(FILE *f, const char *file,
                           const char *keyword, int *line_count) {
  bool test;
  while ((test = ReadAndSplitLine(f, SPL_STR, " \t\n"))) {
    (*line_count)++;
    if (words > 0 && strncasecmp(split[0], keyword, 6) == 0) {
      break;
    }
  }
  // error - missing 'Species' line
  if (!test) {
    snprintf(ERROR_MSG, LINE, "missing '%s' line", keyword);
    ErrorEOF(file, ERROR_MSG);
    exit(1);
  }
} //}}}
static void FieldReadSpecies(const char *file, SYSTEM *System) { //{{{
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  SkipTilKeyword(fr, file, "Species", &line_count);
  COUNT *Count = &System->Count;
  // read number of bead types //{{{
  long types;
  if (words < 2 || !IsWholeNumber(split[1], &types)) {
    err_msg("incorrect 'Species' keyword line");
    PrintErrorFileLine(file, line_count);
    exit(1);
  } //}}}
  // read bead types //{{{
  bool warned = false; // warn about undefined mass/charge only once
  for (int i = 0; i < types; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Species' section");
      exit(1);
    }
    double mass, charge;
    long number;
    // error - illegal 'species' line //{{{
    if (words < 4 ||
        (!IsPosRealNumber(split[1], &mass) && strcmp(split[1], "???") != 0) ||
        (!IsRealNumber(split[2], &charge) && strcmp(split[2], "???") != 0) ||
        !IsWholeNumber(split[3], &number)) {
      err_msg("incorrect 'Species' line");
      PrintErrorFileLine(file, line_count);
      exit(1);
    } //}}}
    if (strcmp(split[1], "???") == 0) {
      mass = MASS;
      if (!warned) {
        err_msg("undefined charge/mass (\'???\' in species section)");
        PrintWarnFile(file, "\0", "\0");
        putc('\n', stderr);
        warned = true;
      }
    }
    if (strcmp(split[2], "???") == 0) {
      charge = CHARGE;
      if (!warned) {
        err_msg("undefined charge/mass (\'???\' in species section)");
        PrintWarnFile(file, "\0", "\0");
        putc('\n', stderr);
        warned = true;
      }
    }
    NewBeadType(&System->BeadType, &Count->BeadType, split[0], charge, mass,
                HIGHNUM);
    System->BeadType[Count->BeadType-1].Number = number;
    Count->Bead += number;
    Count->Unbonded += number;
  } //}}}
  fclose(fr);
} //}}}
static void FieldReadMolecules(const char *file, SYSTEM *System) { //{{{
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  SkipTilKeyword(fr, file, "Molecules", &line_count);
  // read number of types //{{{
  long val;
  if (words < 2 || !IsWholeNumber(split[1], &val)) {
    err_msg("incorrect 'Molecules' keyword line");
    goto error;
  } //}}}
  COUNT *Count = &System->Count;
  Count->MoleculeType = val;
  if (Count->MoleculeType > 0) {
    System->MoleculeType = s_realloc(System->MoleculeType,
                                     sizeof *System->MoleculeType *
                                     Count->MoleculeType);
    // read molecule types
    /*
     * Order of entries:
     *   1) <molecule name>
     *   2) nummol[s] <int>
     *   3) bead[s] <int>
     *      <int> lines: <bead name> <x> <y> <z>
     *   4) bond[s] <int> (optional)
     *      <int> lines: harm <id1> <id2> <k> <r_0>
     *   5) angle[s] <int> (optional)
     *      <int> lines: harm <id1> <id2> <id3> <k> <theta_0>
     *   6) dihed[rals] <int> (optional)
     *      <int> lines: harm <id1> <id2> <id3> <id4> <k>
     *   7) improp[ers] <int> (optional)
     *      <int> lines: harm <id1> <id2> <id3> <id4> <k>
     *   8) finish keyword
     */
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      InitMoleculeType(mt_i);
      // 1) name //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words == 0) {
        err_msg("missing molecule name");
        goto error;
      }
      s_strcpy(mt_i->Name, split[0], MOL_NAME); //}}}
      // 2) number of molecules //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words < 2 || strcasecmp(split[0], "nummols") != 0 ||
          // !IsNaturalNumber(split[1], &val)) {
          !IsIntegerNumber(split[1], &val) || val < 0) {
        err_msg("incorrect 'nummols' line in a molecule entry");
        goto error;
      } else if (val == 0) {
        do {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
        } while(strcmp(split[0], "finish") != 0);
        i--;
        Count->MoleculeType--;
        continue;
      }
      mt_i->Number = val; //}}}
      // 3) beads in the molecule //{{{
      // a) number of beads //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words < 2 || strncasecmp(split[0], "beads", 4) != 0 ||
          !IsNaturalNumber(split[1], &val)) {
        err_msg("incorrect 'beads' line in a molecule entry");
        goto error;
      }
      mt_i->nBeads = val; //}}}
      mt_i->Bead = malloc(sizeof *mt_i->Bead * mt_i->nBeads);
      double (*coor)[3] = malloc(sizeof *coor * mt_i->nBeads);
      // b) beads themselves //{{{
      for (int j = 0; j < mt_i->nBeads; j++) {
        line_count++;
        // read a line //{{{
        if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
          ErrorEOF(file, "incomplete 'Molecules' section");
          exit(1);
        } //}}}
        // error - incorrect line //{{{
        if (words < 4 ||
            !IsRealNumber(split[1], &coor[j][0]) ||
            !IsRealNumber(split[2], &coor[j][1]) ||
            !IsRealNumber(split[3], &coor[j][2])) {
          err_msg("incorrect bead coordinate line in a molecule entry");
          goto error;
        } //}}}
        int btype = FindBeadType(split[0], *System);
        // error - unknown type //{{{
        if (btype == -1) {
          err_msg("unknown bead in a molecule entry");
          PrintErrorFileLine(file, line_count);
          ErrorBeadType(split[0], *System);
          exit(1);
        } //}}}
        mt_i->Bead[j] = btype;
      }                                //}}}
      int count = Count->Bead,         // for filling Molecule[] & Bead[]
          mol_count = Count->Molecule; //
      Count->Molecule += mt_i->Number;
      Count->Bead += mt_i->Number * mt_i->nBeads;
      Count->Bonded += mt_i->Number * mt_i->nBeads;
      System->Molecule = s_realloc(System->Molecule,
                                   sizeof *System->Molecule * Count->Molecule);
      System->Bead = s_realloc(System->Bead,
                               sizeof *System->Bead * Count->Bead);
      System->Bonded = s_realloc(System->Bonded,
                                 sizeof System->Bonded * Count->Bonded);
      for (int j = count; j < Count->Bead; j++) {
        InitBead(&System->Bead[j]);
      }
      // c) fill Molecule[] & Bead[] //{{{
      for (int j = 0; j < mt_i->Number; j++) {
        MOLECULE *mol = &System->Molecule[mol_count];
        mol->Type = i;
        mol->Index = mol_count;
        mol->Bead = malloc(sizeof *mol->Bead * mt_i->nBeads);
        mol->InTimestep = true;
        for (int k = 0; k < mt_i->nBeads; k++) {
          BEAD *bead = &System->Bead[count];
          bead->Type = mt_i->Bead[k];
          bead->Molecule = mol_count;
          for (int dd = 0; dd < 3; dd++) {
            bead->Position[dd] = coor[k][dd];
          }
          System->BeadType[bead->Type].Number++;
          System->Bonded[count - Count->Unbonded] = count;
          mol->Bead[k] = count;
          count++;
        }
        mol_count++;
      }
      // pro forma checks //{{{
      if (count != Count->Bead) {
        err_msg("wrong number of beads; should never happen!");
        PrintErrorFile(file, "\0", "\0");
        exit(1);
      }
      if (mol_count != Count->Molecule) {
        err_msg("wrong number of molecules; should never happen!");
        PrintErrorFile(file, "\0", "\0");
        exit(1);
      } //}}}
      //}}}
      free(coor); //}}}
      // 4) bonds in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "bonds", 4) == 0) {
        // a) number of bonds //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          err_msg("incorrect 'bond' line in a molecule entry");
          goto error;
        }
        mt_i->nBonds = val; //}}}
        // b) bonds themselves & bond types
        // TODO: for now, only harmonic bonds are considered
        if (!ReadStuff(file, fr, &line_count, 0, i, System)) {
          goto error;
        }
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE, "unrecognised line in an entry for "
                 "molecule %s%s", ErrYellow(), mt_i->Name);
        goto error;
      } //}}}
      // 5) angles in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "angles", 4) == 0) {
        // a) number of angles //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          err_msg("incorrect 'angle' line in a molecule entry");
          goto error;
        }
        mt_i->nAngles = val; //}}}
        // b) angles themselves & angle types
        // TODO: for now, only harmonic angles are considered
        if (!ReadStuff(file, fr, &line_count, 1, i, System)) {
          goto error;
        }
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE, "unrecognised line in an entry for molecule "
                 "%s%s", ErrYellow(), mt_i->Name);
        goto error;
      } //}}}
      // 6) dihedrals in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "dihedrals", 5) == 0) {
        // a) number of dihedrals //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          err_msg("incorrect 'dihedral' line in a molecule entry");
          goto error;
        }
        mt_i->nDihedrals = val; //}}}
        // b) dihedrals themselves & dihedral types
        // TODO: for now, only harmonic dihedrals are considered
        if (!ReadStuff(file, fr, &line_count, 2, i, System)) {
          goto error;
        }
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE, "unrecognised line in an entry for molecule "
                 "%s%s", ErrYellow(), mt_i->Name);
        goto error;
      } //}}}
      // 7) impropers in the molecule (if present) //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words > 1 && strncasecmp(split[0], "impropers", 6) == 0) {
        // a) number of impropers //{{{
        if (!IsNaturalNumber(split[1], &val)) {
          err_msg("incorrect 'improper' line in a molecule entry");
          goto error;
        }
        mt_i->nImpropers = val; //}}}
        // b) impropers themselves & improper types
        // TODO: for now, only lammps' cvff impropers are considered
        if (!ReadStuff(file, fr, &line_count, 3, i, System)) {
          goto error;
        }
      } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
        continue;
      } else {
        snprintf(ERROR_MSG, LINE, "unrecognised line in an entry for molecule "
                 "%s%s", ErrYellow(), mt_i->Name);
        PrintErrorFileLine(file, line_count);
        exit(1);
      } //}}}
      // finish keyword //{{{
      line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
        ErrorEOF(file, "incomplete 'Molecules' section");
        exit(1);
      } //}}}
      if (words == 0 || strcasecmp(split[0], "finish") != 0) {
        snprintf(ERROR_MSG, LINE,
                 "missing 'finish' at the end of an entry for "
                 "molecule %s%s",
                 ErrYellow(), mt_i->Name);
        PrintErrorFileLine(file, line_count);
        exit(1);
      } //}}}
    }
    Count->HighestResid = Count->Molecule;
  }
  fclose(fr);
  return;
error: // unrecognised line
  PrintErrorFileLine(file, line_count);
  exit(1);
} //}}}

SYSTEM FieldRead(const char *file) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  // read box dimensions, if present //{{{
  FILE *fr = OpenFile(file, "r");
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "empty file");
    exit(1);
  }
  double box[3];
  if (words > 2 && IsRealNumber(split[0], &box[0]) &&
                   IsRealNumber(split[1], &box[1]) &&
                   IsRealNumber(split[2], &box[2])) {
    System.Box.Length[0] = box[0];
    System.Box.Length[1] = box[1];
    System.Box.Length[2] = box[2];
    CalculateBoxData(&System.Box, 0);
  }
  fclose(fr); //}}}
  FieldReadSpecies(file, &System);
  // fill System.Bead & System.Unbonded //{{{
  if (Count->Unbonded > 0) {
    System.Bead = s_realloc(System.Bead, Count->Bead * sizeof *System.Bead);
    System.Unbonded = s_realloc(System.Unbonded,
                                sizeof *System.Unbonded * Count->Unbonded);
    System.BeadCoor = s_realloc(System.BeadCoor,
                                Count->Bead * sizeof *System.BeadCoor);
    int count = 0;
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = 0; j < System.BeadType[i].Number; j++) {
        InitBead(&System.Bead[count]);
        System.Bead[count].Type = i;
        System.Bead[count].Molecule = -1;
        System.Unbonded[count] = count;
        count++;
      }
    }
    if (count != Count->Unbonded) {
      err_msg("wrong number of unbonded beads; should never happen!");
      PrintErrorFile(file, "\0", "\0");
      exit(1);
    }
  } //}}}
  FieldReadMolecules(file, &System);
  RemoveExtraTypes(&System);
  MergeBeadTypes(&System, true);
  MergeMoleculeTypes(&System);
  if (Count->Bead == 0) {
    err_msg("No beads in the file");
    PrintWarnFile(file, "\0", "\0");
  } else {
    System.BeadCoor = s_realloc(System.BeadCoor,
                                Count->Bead * sizeof *System.BeadCoor);
  }
  if (Count->Molecule > 0) {
    System.MoleculeCoor = s_realloc(System.MoleculeCoor, Count->Molecule *
                                    sizeof *System.MoleculeCoor);
  }
  FillSystemNonessentials(&System);
  CheckSystem(System, file);
  return System;
} //}}}
void WriteField(SYSTEM System, char *file_field, int argc, char **argv) { //{{{
  FILE *fw = OpenFile(file_field, "w");
  BOX *box = &System.Box;
  if (box->Volume != -1) {
    fprintf(fw, "%.3f %.3f %.3f ",
            box->Length[0], box->Length[1], box->Length[2]);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 0) {
      fprintf(fw, "%lf %lf %lf ", box->alpha, box->beta, box->gamma);
    }
  }
  fprintf(fw, "Created via AnalysisTools v%s", VERSION);
  fprintf(fw, "(https://github.com/KaGaSi/AnalysisTools); Command: ");
  PrintCommand(fw, argc, argv);
  putc('\n', fw);
  COUNT *Count = &System.Count;
  // print species section //{{{
  fprintf(fw, "species %d <name> <m> <q> <# of unbonded>\n", Count->BeadType);
  // count unbonded beads of each type
  int *unbonded = calloc(Count->BeadType, sizeof *unbonded);
  for (int i = 0; i < Count->Unbonded; i++) {
    int id = System.Unbonded[i];
    int btype = System.Bead[id].Type;
    unbonded[btype]++;
  }
  // print the lines
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt_i = &System.BeadType[i];
    fprintf(fw, "%16s", bt_i->Name);
    if (bt_i->Mass == MASS) {
      fprintf(fw, " %8s", "???");
    } else {
      fprintf(fw, " %8.5f", bt_i->Mass);
    }
    if (bt_i->Charge == CHARGE) {
      fprintf(fw, " %8s", "???");
    } else {
      fprintf(fw, " %8.5f", bt_i->Charge);
    }
    fprintf(fw, " %5d\n", unbonded[i]);
  }
  free(unbonded); //}}}
  // print molecules section //{{{
  fprintf(fw, "molecule %d\n", Count->MoleculeType);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    fprintf(fw, "%s\n", mt_i->Name);
    fprintf(fw, "nummols %d\n", mt_i->Number);
    fprintf(fw, "beads %d\n", mt_i->nBeads);
    int mol = mt_i->Index[0];
    // beads
    for (int j = 0; j < mt_i->nBeads; j++) {
      int id = System.Molecule[mol].Bead[j];
      int bt = mt_i->Bead[j];
      double (*pos)[3] = &System.Bead[id].Position;
      fprintf(fw, "%16s %8.5f %8.5f %8.5f\n", System.BeadType[bt].Name,
              (*pos)[0], (*pos)[1], (*pos)[2]);
    }
    // bonds (if present)
    if (mt_i->nBonds > 0) {
      fprintf(fw, "bonds %d\n", mt_i->nBonds);
      for (int j = 0; j < mt_i->nBonds; j++) {
        // TODO harm only for now
        fprintf(fw, "harm %5d %5d", mt_i->Bond[j][0] + 1, mt_i->Bond[j][1] + 1);
        int type = mt_i->Bond[j][2];
        if (type != -1) {
          fprintf(fw, " %lf %lf\n", System.BondType[type].a,
                  System.BondType[type].b);
        } else {
          fprintf(fw, " ???   ???\n");
        }
      }
    }
    // angles (if present)
    if (mt_i->nAngles > 0) {
      fprintf(fw, "angles %d\n", mt_i->nAngles);
      for (int j = 0; j < mt_i->nAngles; j++) {
        // TODO harm only for now
        fprintf(fw, "harm %5d %5d %5d", mt_i->Angle[j][0] + 1,
                mt_i->Angle[j][1] + 1,
                mt_i->Angle[j][2] + 1);
        int type = mt_i->Angle[j][3];
        if (type != -1) {
          fprintf(fw, " %lf %lf\n", System.AngleType[type].a,
                  System.AngleType[type].b);
        } else {
          fprintf(fw, " ???   ???\n");
        }
      }
    }
    // dihedrals (if present)
    if (mt_i->nDihedrals > 0) {
      fprintf(fw, "dihedrals %d ", mt_i->nDihedrals);
      fprintf(fw, "# lammps' harmonic style\n");
      for (int j = 0; j < mt_i->nDihedrals; j++) {
        // TODO harm (lammps) only for now
        fprintf(fw, "harm %5d %5d %5d %5d", mt_i->Dihedral[j][0] + 1,
                mt_i->Dihedral[j][1] + 1,
                mt_i->Dihedral[j][2] + 1,
                mt_i->Dihedral[j][3] + 1);
        int type = mt_i->Dihedral[j][4];
        if (type != -1) {
          fprintf(fw, " %lf %lf %lf\n", System.DihedralType[type].a,
                  System.DihedralType[type].b,
                  System.DihedralType[type].c);
        } else {
          fprintf(fw, " ???\n");
        }
      }
    }
    // impropers (if present)
    if (mt_i->nImpropers > 0) {
      fprintf(fw, "impropers %d ", mt_i->nImpropers);
      fprintf(fw, "# lammps' cvff style\n");
      for (int j = 0; j < mt_i->nImpropers; j++) {
        // TODO harm only for now
        fprintf(fw, "cvff %5d %5d %5d %5d", mt_i->Improper[j][0] + 1,
                mt_i->Improper[j][1] + 1,
                mt_i->Improper[j][2] + 1,
                mt_i->Improper[j][3] + 1);
        int type = mt_i->Improper[j][4];
        if (type != -1) {
          fprintf(fw, " %lf %lf %lf\n", System.ImproperType[type].a,
                  System.ImproperType[type].b,
                  System.ImproperType[type].c);
        } else {
          fprintf(fw, " ???\n");
        }
      }
    }
    fprintf(fw, "finish\n");
  } //}}}
  fclose(fw);
} //}}}
