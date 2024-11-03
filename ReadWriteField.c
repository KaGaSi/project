#include "ReadWriteField.h"
#include "General.h"

/*
 * Function to read dl_meso FIELD-like file as a structure file
 */
// Helper functions for FIELD-like file
// read Species section
static void FieldReadSpecies(const char *file, SYSTEM *System);
// read Molecules section
static void FieldReadMolecules(const char *file, SYSTEM *System);

static void FieldReadSpecies(const char *file, SYSTEM *System) { //{{{
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // skip till line starting with 'Species' keyword //{{{
  bool test;
  while ((test = ReadAndSplitLine(fr, SPL_STR, " \t\n"))) {
    line_count++;
    if (words > 0 && strncasecmp(split[0], "species", 6) == 0) {
      break;
    }
  }
  // error - missing 'Species' line
  if (!test) {
    ErrorEOF(file, "missing 'Species' line");
    exit(1);
  } //}}}
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
  // skip till line starting with 'Molecule' keyword //{{{
  bool test;
  while ((test = ReadAndSplitLine(fr, SPL_STR, " \t\n"))) {
    line_count++;
    if (words > 0 && strncasecmp(split[0], "molecules", 7) == 0) {
      break;
    }
  }
  // error - missing 'Species' line
  if (!test) {
    ErrorEOF(file, "missing 'Molecules' line");
    exit(1);
  } //}}}
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
          err_msg("incorrect 'bonds' line in a molecule entry");
          PrintErrorFileLine(file, line_count);
          exit(1);
        }
        mt_i->nBonds = val; //}}}
        mt_i->Bond = malloc(sizeof *mt_i->Bond * mt_i->nBonds);
        // b) bonds themselves & bond types //{{{
        // TODO: for now, only harmonic bonds are considered
        bool warned = false;
        for (int j = 0; j < mt_i->nBonds; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[2];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 3 ||
              !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              beads[0] == beads[1]) {
            err_msg("incorrect bond line in a molecule entry");
            goto error;
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads) {
            err_msg("bead index in a bond is too high");
            goto error;
          }                                //}}}
          mt_i->Bond[j][0] = beads[0] - 1; // in FIELD, indices start from 1
          mt_i->Bond[j][1] = beads[1] - 1; //
          // read up to three values for bond type
          for (int k = 3; k <= 5; k++) {
            double *ptr = NULL;
            if (k == 3) {
              ptr = &values.a;
            } else if (k == 4) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined bond parameter (\'???\') "
                         "in molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          }
          // assign bond type //{{{
          int bond_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this bond type already exists
            for (int k = 0; k < Count->BondType; k++) {
              if (System->BondType[k].a == values.a &&
                  System->BondType[k].b == values.b &&
                  System->BondType[k].c == values.c) {
                bond_type = k;
                break;
              }
            }
            // create a new bond type if necessary
            if (bond_type == -1) {
              bond_type = Count->BondType;
              Count->BondType++;
              System->BondType = s_realloc(System->BondType, Count->BondType *
                                           sizeof *System->BondType);
              System->BondType[bond_type] = InitParams;
              System->BondType[bond_type].a = values.a;
              System->BondType[bond_type].b = values.b;
              System->BondType[bond_type].c = values.c;
            }
          } //}}}
          mt_i->Bond[j][2] = bond_type;
        } //}}}
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
          snprintf(ERROR_MSG, LINE, "incorrect line in 'angles' line for "
                   "molecule %s%s", ErrYellow(), mt_i->Name);
          goto error;
        }
        mt_i->nAngles = val; //}}}
        mt_i->Angle = malloc(sizeof *mt_i->Angle * mt_i->nAngles);
        // b) angles themselves & angle types //{{{
        // TODO: for now, only harmonic angles are considered
        bool warned = false; // to warn of undefined angle type only once
        for (int j = 0; j < mt_i->nAngles; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[3];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 4 || !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              !IsNaturalNumber(split[3], &beads[2]) || beads[0] == beads[1] ||
              beads[0] == beads[2] || beads[1] == beads[2]) {
            snprintf(ERROR_MSG, LINE, "incorrect line in an entry for "
                     "molecule %s%s", ErrYellow(), mt_i->Name);
            goto error;
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads ||
              beads[2] > mt_i->nBeads) {
            snprintf(ERROR_MSG, LINE, "bead index in an angle is too high in "
                     "molecule %s%s", ErrYellow(), mt_i->Name);
            goto error;
          } //}}}
          mt_i->Angle[j][0] = beads[0] - 1; // in FIELD, bead indices go from 1
          mt_i->Angle[j][1] = beads[1] - 1; //
          mt_i->Angle[j][2] = beads[2] - 1; //
          // read up to three values for angle type
          for (int k = 4; k <= 6; k++) {
            double *ptr = NULL;
            if (k == 4) {
              ptr = &values.a;
            } else if (k == 5) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined angle parameter (\'???\') "
                         " for molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          }
          // assign angle type //{{{
          int angle_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this angle type already exists
            for (int k = 0; k < Count->AngleType; k++) {
              if (System->AngleType[k].a == values.a &&
                  System->AngleType[k].b == values.b &&
                  System->AngleType[k].c == values.c) {
                angle_type = k;
                break;
              }
            }
            // create a new angle type if necessary
            if (angle_type == -1) {
              angle_type = Count->AngleType;
              Count->AngleType++;
              System->AngleType = s_realloc(System->AngleType,
                                            Count->AngleType *
                                            sizeof *System->AngleType);
              System->AngleType[angle_type] = InitParams;
              System->AngleType[angle_type].a = values.a;
              System->AngleType[angle_type].b = values.b;
              System->AngleType[angle_type].c = values.c;
            }
          } //}}}
          mt_i->Angle[j][3] = angle_type;
        } //}}}
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
          snprintf(ERROR_MSG, LINE, "incorrect 'dihedrals' line for molecule "
                   "%s%s", ErrYellow(), mt_i->Name);
          PrintErrorFileLine(file, line_count);
          exit(1);
        }
        mt_i->nDihedrals = val; //}}}
        mt_i->Dihedral = malloc(sizeof *mt_i->Dihedral * mt_i->nDihedrals);
        // b) dihedrals themselves & dihedral types //{{{
        // TODO: for now, only harmonic dihedrals are considered
        bool warned = false; // to warn of undefined type only once
        for (int j = 0; j < mt_i->nDihedrals; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[4];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 5 ||
              !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              !IsNaturalNumber(split[3], &beads[2]) ||
              !IsNaturalNumber(split[4], &beads[3]) ||
              beads[0] == beads[1] ||
              beads[0] == beads[2] ||
              beads[0] == beads[3] ||
              beads[1] == beads[2] ||
              beads[1] == beads[3] ||
              beads[2] == beads[3]) {
            snprintf(ERROR_MSG, LINE, "incorrect dihedral line for molecule "
                     "%s%s", ErrYellow(), mt_i->Name);
            PrintErrorFileLine(file, line_count);
            exit(1);
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads ||
              beads[1] > mt_i->nBeads ||
              beads[2] > mt_i->nBeads ||
              beads[3] > mt_i->nBeads) {
            err_msg("bead index in a improper is too high");
            goto error;
          }                                    //}}}
          mt_i->Dihedral[j][0] = beads[0] - 1; // in FIELD, indices start from 1
          mt_i->Dihedral[j][1] = beads[1] - 1; //
          mt_i->Dihedral[j][2] = beads[2] - 1; //
          mt_i->Dihedral[j][3] = beads[3] - 1; //
          // read up to three values for dihedral type //{{{
          for (int k = 5; k <= 7; k++) {
            double *ptr = NULL;
            if (k == 5) {
              ptr = &values.a;
            } else if (k == 6) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined dihedral paramter "
                         "(\'???\') in molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          } //}}}
          // assign dihedral type //{{{
          int dihedral_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this dihedral type already exists
            for (int k = 0; k < Count->DihedralType; k++) {
              if (System->DihedralType[k].a == values.a &&
                  System->DihedralType[k].b == values.b &&
                  System->DihedralType[k].c == values.c) {
                dihedral_type = k;
                break;
              }
            }
            // create a new dihedral type if necessary
            if (dihedral_type == -1) {
              dihedral_type = Count->DihedralType;
              Count->DihedralType++;
              System->DihedralType = s_realloc(System->DihedralType,
                                               Count->DihedralType *
                                               sizeof *System->DihedralType);
              System->DihedralType[dihedral_type] = InitParams;
              System->DihedralType[dihedral_type] = values;
            }
          }
          mt_i->Dihedral[j][4] = dihedral_type; //}}}
        }                                       //}}}
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
          snprintf(ERROR_MSG, LINE,
                   "incorrect 'impropers' line in an entry "
                   "for molecule %s%s",
                   ErrYellow(), mt_i->Name);
          PrintErrorFileLine(file, line_count);
          exit(1);
        }
        mt_i->nImpropers = val; //}}}
        mt_i->Improper = malloc(sizeof *mt_i->Improper * mt_i->nImpropers);
        // b) impropers themselves & improper types //{{{
        bool warned = false; // to warn of undefined type only once
        for (int j = 0; j < mt_i->nImpropers; j++) {
          line_count++;
          // read a line //{{{
          if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
            ErrorEOF(file, "incomplete 'Molecules' section");
            exit(1);
          } //}}}
          long beads[4];
          PARAMS values = InitParams;
          // error - incorrect line //{{{
          if (words < 5 ||
              !IsNaturalNumber(split[1], &beads[0]) ||
              !IsNaturalNumber(split[2], &beads[1]) ||
              !IsNaturalNumber(split[3], &beads[2]) ||
              !IsNaturalNumber(split[4], &beads[3]) ||
              beads[0] == beads[1] ||
              beads[0] == beads[2] ||
              beads[0] == beads[3] ||
              beads[1] == beads[2] ||
              beads[1] == beads[3] ||
              beads[2] == beads[3]) {
            snprintf(ERROR_MSG, LINE, "incorrect improper line for molecule "
                     "%s%s", ErrYellow(), mt_i->Name);
            PrintErrorFileLine(file, line_count);
            exit(1);
          } //}}}
          // error - bead index is too high //{{{
          if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads ||
              beads[2] > mt_i->nBeads || beads[3] > mt_i->nBeads) {
            err_msg("bead index in a improper is too high");
            goto error;
          }                                    //}}}
          mt_i->Improper[j][0] = beads[0] - 1; // in FIELD, indices start from 1
          mt_i->Improper[j][1] = beads[1] - 1; //
          mt_i->Improper[j][2] = beads[2] - 1; //
          mt_i->Improper[j][3] = beads[3] - 1; //
          // read up to three values for improper type //{{{
          for (int k = 5; k <= 7; k++) {
            double *ptr = NULL;
            if (k == 5) {
              ptr = &values.a;
            } else if (k == 6) {
              ptr = &values.b;
            } else {
              ptr = &values.c;
            }
            if (words > k) {
              if (strcmp(split[k], "???") == 0 && !warned) {
                snprintf(ERROR_MSG, LINE, "undefined improper parameter "
                         "(\'???\') in molecule %s%s", ErrYellow(), mt_i->Name);
                PrintWarnFile(file, "\0", "\0");
                warned = true;
              } else if (!IsRealNumber(split[k], ptr)) {
                break;
              }
            }
          } //}}}
          // assign improper type //{{{
          int improper_type = -1;
          if (values.a != 0 || values.b != 0 || values.c != 0) {
            // find if this improper type already exists
            for (int k = 0; k < Count->ImproperType; k++) {
              if (System->ImproperType[k].a == values.a &&
                  System->ImproperType[k].b == values.b &&
                  System->ImproperType[k].c == values.c) {
                improper_type = k;
                break;
              }
            }
            // create a new improper type if necessary
            if (improper_type == -1) {
              improper_type = Count->ImproperType;
              Count->ImproperType++;
              System->ImproperType = s_realloc(System->ImproperType,
                                               Count->ImproperType *
                                               sizeof *System->ImproperType);
              System->ImproperType[improper_type] = InitParams;
              System->ImproperType[improper_type] = values;
            }
          }
          mt_i->Improper[j][4] = improper_type; //}}}
        }                                       //}}}
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
error: // unrecognised line //{{{
  PrintErrorFileLine(file, line_count);
  exit(1); //}}}
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
