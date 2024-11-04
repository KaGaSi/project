#include "ReadWriteLdata.h"
#include "General.h"
#include "ReadWrite.h"
#include "System.h"

// TODO: LmpDataReadDihedralCoeffs() and LmpDataReadImproperCoeffs() should read
//       up to three numbers, not assuming any format of the potential

// Helper functions for lmpdata file
// read header of the lammps data file
static int LmpDataReadHeader(FILE *fr, const char *file,
                             SYSTEM *System, int *line_count);
// read body of the lammps data file
static void LmpDataReadBody(FILE *fr, const char *file, SYSTEM *System,
                            const int atom_types, int *line_count);
// functions to read the various sections in the lammps data file
static void LmpDataReadMasses(FILE *fr, const char *file, BEADTYPE *name_mass,
                              const int atom_types, int *line_count);
// TODO: non-harmonic bonds
static void LmpDataReadBondCoeffs(FILE *fr, const char *file,
                                  SYSTEM *System, int *line_count);
// TODO: non-harmonic angles
static void LmpDataReadAngleCoeffs(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count);
// TODO: non-harmonic (lammps type) angles
static void LmpDataReadDihedralCoeffs(FILE *fr, const char *file,
                                      SYSTEM *System, int *line_count);
// TODO: non-cvff (lammps type) impropers
static void LmpDataReadImproperCoeffs(FILE *fr, const char *file,
                                      SYSTEM *System, int *line_count);
static void LmpDataReadAtoms(FILE *fr, const char *file, SYSTEM *System,
                             const BEADTYPE *name_mass, const int atom_types,
                             int *line_count, const int mode);
static void LmpDataReadVelocities(FILE *fr, const char *file,
                                  SYSTEM *System, int *line_count);
static void LmpDataReadBonds(FILE *fr, const char *file, const COUNT Count,
                             int (*bond)[3], int *line_count);
static void LmpDataReadAngles(FILE *fr, const char *file, const COUNT Count,
                              int (*angle)[4], int *line_count);
static void LmpDataReadDihedrals(FILE *fr, const char *file, const COUNT Count,
                                 int (*dihedral)[5], int *line_count);
static void LmpDataReadImpropers(FILE *fr, const char *file, const COUNT Count,
                                 int (*improper)[5], int *line_count);
static bool ReadPbc(BOX *box);

// read header //{{{
static int LmpDataReadHeader(FILE *fr, const char *file,
                             SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // ignore the first line (comment) //{{{
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "\0");
    exit(1);
  } //}}}
  *line_count = 1;
  // read until a line starting with capital letter //{{{
  double atom_types = 0;
  fpos_t position;
  do {
    (*line_count)++;
    fgetpos(fr, &position);
    // read a line
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete lammps data file header");
      exit(1);
    }
    // evaluate the line
    long val;
    // <int> atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      Count->Bead = val;
      System->Bead = s_realloc(System->Bead,
                               Count->Bead * sizeof *System->Bead);
      /* Count->Bead = Count->BeadType as each bead can have different charge,
       * so at first, each bead will be its own type
       */
      Count->BeadType = val;
      System->BeadType = s_realloc(System->BeadType,
                                   Count->BeadType * sizeof *System->BeadType);
      for (int i = 0; i < Count->Bead; i++) {
        InitBead(&System->Bead[i]);
        InitBeadType(&System->BeadType[i]);
      } //}}}
    // <int> bonds //{{{
    } else if (words > 1 && strcmp(split[1], "bonds") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Bond = val; //}}}
    // <int> angles //{{{
    } else if (words > 1 && strcmp(split[1], "angles") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Angle = val; //}}}
    // <int> dihedrals //{{{
    } else if (words > 1 && strcmp(split[1], "dihedrals") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Dihedral = val; //}}}
    // <int> impropers //{{{
    } else if (words > 1 && strcmp(split[1], "impropers") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->Improper = val; //}}}
    // <int> atom types //{{{
    } else if (words > 2 && strcmp(split[1], "atom") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      atom_types = val; //}}}
    // <int> bond types //{{{
    } else if (words > 2 && strcmp(split[1], "bond") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (strcmp(split[0], "???") != 0) {
        if (!IsWholeNumber(split[0], &val)) {
          goto error;
        }
        Count->BondType = val;
        if (Count->BondType > 0) {
          System->BondType = s_realloc(System->BondType, Count->BondType *
                                       sizeof *System->BondType);
          for (int i = 0; i < Count->BondType; i++) {
            System->BondType[i] = InitParams;
          }
        }
      }
      //}}}
    // <int> angle types //{{{
    } else if (words > 2 && strcmp(split[1], "angle") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->AngleType = val;
      if (Count->AngleType > 0) {
        System->AngleType = s_realloc(System->AngleType, Count->AngleType *
                                      sizeof *System->AngleType);
        for (int i = 0; i < Count->AngleType; i++) {
          System->AngleType[i] = InitParams;
        }
      }
      //}}}
    // <int> dihedral types //{{{
    } else if (words > 2 && strcmp(split[1], "dihedral") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->DihedralType = val;
      if (Count->DihedralType > 0) {
        System->DihedralType = s_realloc(System->DihedralType,
                                         Count->DihedralType *
                                         sizeof *System->DihedralType);
        for (int i = 0; i < Count->DihedralType; i++) {
          System->DihedralType[i] = InitParams;
        }
      }
      //}}}
    // <int> improper types //{{{
    } else if (words > 2 && strcmp(split[1], "improper") == 0 &&
               strcmp(split[2], "types") == 0) {
      if (!IsWholeNumber(split[0], &val)) {
        goto error;
      }
      Count->ImproperType = val;
      if (Count->ImproperType > 0) {
        System->ImproperType = s_realloc(System->ImproperType,
                                         Count->ImproperType *
                                         sizeof *System->ImproperType);
        for (int i = 0; i < Count->ImproperType; i++) {
          System->ImproperType[i] = InitParams;
        }
      }
    }
      //}}}
    if (!ReadPbc(&System->Box)) {
      goto error;
    }
  } while (words == 0 || split[0][0] < 'A' || split[0][0] > 'Z'); //}}}
  // return file pointer to before the first capital-letter-starting line
  fsetpos(fr, &position);
  (*line_count)--;
  // check all necessary parts were read //{{{
  if (Count->Bead == 0 || Count->BeadType == 0) {
    err_msg("missing atom or atom types in the file header");
    PrintErrorFile(file, "\0", "\0");
    putc('\n', stderr);
    exit(1);
  }
  if (System->Box.OrthoLength[0] == -1 ||
      System->Box.OrthoLength[1] == -1 ||
      System->Box.OrthoLength[2] == -1) {
    err_msg("missing box size in the file header");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
  }
  //}}}
  return atom_types;
error: // unrecognised line //{{{
  err_msg("unrecognised line in the file header");
  PrintErrorFileLine(file, *line_count);
  exit(1); //}}}
} //}}}
SYSTEM LmpDataReadStruct(const char *file) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  FILE *fr = OpenFile(file, "r");
  int line_count = 0;
  int lmp_types = LmpDataReadHeader(fr, file, &System, &line_count);
  LmpDataReadBody(fr, file, &System, lmp_types, &line_count);
  fclose(fr);
  CalculateBoxData(&System.Box, 1);
  RemoveExtraTypes(&System);
  MergeBeadTypes(&System, true);
  MergeMoleculeTypes(&System);
  FillSystemNonessentials(&System);
  // take the data file as a coordinate file as well
  Count->BeadCoor = Count->Bead;
  Count->UnbondedCoor = Count->Unbonded;
  Count->BondedCoor = Count->Bonded;
  System.BeadCoor = s_realloc(System.BeadCoor,
                              Count->Bead * sizeof *System.BeadCoor);
  int c_unbonded = 0, c_bonded = 0;
  for (int i = 0; i < Count->Bead; i++) {
    System.BeadCoor[i] = i;
    if (System.Bead[i].Molecule == -1) {
      System.UnbondedCoor[c_unbonded] = i;
      c_unbonded++;
    } else {
      System.BondedCoor[c_bonded] = i;
      c_bonded++;
    }
  }
  // make molecule indices go from 0 (just to avoid large numbers) //{{{
  // 1) find the lowest molecule index
  int min_id = 1e7;
  for (int i = 0; i < System.Count.Molecule; i++) {
    if (System.Molecule[i].Index < min_id) {
      min_id = System.Molecule[i].Index;
    }
  }
  // 2) subtract the lowest molecule index from molecule indices
  for (int i = 0; i < System.Count.Molecule; i++) {
    System.Molecule[i].Index -= min_id;
  } //}}}
  if (Count->Molecule) {
    System.MoleculeCoor = s_realloc(System.MoleculeCoor, Count->Molecule *
                                    sizeof *System.MoleculeCoor);
  }
  Count->MoleculeCoor = Count->Molecule;
  for (int i = 0; i < Count->MoleculeCoor; i++) {
    System.MoleculeCoor[i] = i;
  }
  CheckSystem(System, file);
  return System;
} //}}}
// TODO: make Atoms section type detection (i.e., # <type>) into a function
// TODO: Atoms & Velocities section can be switched
// TODO: add Forces section reading
// read lammps data file as a coordinate file //{{{
int LmpDataReadTimestep(FILE *fr, const char *file,
                        SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // set 'in timestep' to all beads //{{{
  for (int i = 0; i < Count->Bead; i++) {
    System->Bead[i].InTimestep = true;
  } //}}}
  // set 'in timestep' to all molecules //{{{
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].InTimestep = false;
  } //}}}
  Count->MoleculeCoor = Count->Molecule;
  for (int i = 0; i < Count->MoleculeCoor; i++) {
    System->MoleculeCoor[i] = i;
  }

  // ignore the first line (comment) //{{{
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "\0");
    exit(1);
  } //}}}
  *line_count = 1;
  // read numer of atoms & box size //{{{
  // read until the first capital-letter-starting line
  fpos_t position;
  do {
    (*line_count)++;
    fgetpos(fr, &position);
    // read a line
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete lammps data file header");
      exit(1);
    }
    long val;
    // <int> atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsNaturalNumber(split[0], &val) || val == 0) {
        goto error;
      }
      if (Count->Bead != val) {
        err_msg("coordinate lammps data file must contain all beads");
        PrintErrorFile(file, "\0", "\0");
        exit(1);
      } //}}}
    }
    if (!ReadPbc(&System->Box)) {
      goto error;
    }
  } while (words == 0 || split[0][0] < 'A' || split[0][0] > 'Z');
  //  return file pointer to before the first capital-letter-starting line
  fsetpos(fr, &position);
  (*line_count)--;
  CalculateBoxData(&System->Box, 1); //}}}
  // find 'Atoms' section (and skip the next blank line) //{{{
  do {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
  } while (words == 0 || strcmp(split[0], "Atoms") != 0);
  int mode = 0;
  if (words > 2 && split[1][0] == '#') {
    if (strcmp(split[2], "full") == 0) {
      mode = 0; // for 'Atoms # bond|angle|molecular|charge'
    } else if (strcmp(split[2], "bond") == 0 ||
               strcmp(split[2], "angle") == 0 ||
               strcmp(split[2], "molecular") == 0 ||
               strcmp(split[2], "charge") == 0) {
      mode = 1; // for 'Atoms # bond|angle|molecular'
    } else if (strcmp(split[2], "atomic") == 0) {
      mode = 2; // for 'Atoms # atomic'
    } else {    // for 'Atoms # <something else>'
      snprintf(ERROR_MSG, LINE, "unsupported atom_style (Atoms section) "
               "%s%s; assuming 'full'", ErrYellow(), split[2]);
      PrintWarnFile(file, "\0", "\0");
    }
  }
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  } //}}}
  // read atom lines //{{{
  Count->BeadCoor = Count->Bead;
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
    long id;
    double pos[3];
    if (mode == 0) {
      // 'Atoms # full' or just 'Atoms':
      // <bead id> <mol id> <bead type id> <charge> <coordinates>
      if (words < 7 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsRealNumber(split[4], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[5], &pos[1]) || //
          !IsRealNumber(split[6], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    } else if (mode == 1 || mode == 3) {
      // 'Atoms # bond|angle|molecular':
      // <bead id> <mol id> <bead type id> <coordinates>
      // or 'Atoms # charge':
      // <bead id> <bead type id> <charge> <coordinates>
      if (words < 6 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsRealNumber(split[3], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[4], &pos[1]) || //
          !IsRealNumber(split[5], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    } else { // mode == 2
      // 'Atoms # atomic': <bead id> <bead type id> <coordinates>
      if (words < 5 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsRealNumber(split[2], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[3], &pos[1]) || //
          !IsRealNumber(split[4], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    }
    id--; // in lammps data file, ids start from 1
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = pos[dd] - System->Box.Low[dd];
    }
    System->BeadCoor[i] = id;
  } //}}}
  // find 'Velocities' section (and skip the next blank line) //{{{
  do {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ChangeBoxByLow(System, -1);
      FillInCoor(System);
      return 1; // Velocities section is not mandatory
    }
  } while (words == 0 || strcmp(split[0], "Velocities") != 0);
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ChangeBoxByLow(System, -1);
    FillInCoor(System);
    return -2;
  } //}}}
  // read velocity lines //{{{
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      return -2;
    }
    long id;
    double vel[3];
    // error - wrong line //{{{
    if (words < 4 || !IsNaturalNumber(split[0], &id) ||
        id > Count->Bead ||                // bead index
        !IsRealNumber(split[1], &vel[0]) || // bead velocities
        !IsRealNumber(split[2], &vel[1]) || //
        !IsRealNumber(split[3], &vel[2])) { //
      err_msg("wrong line in Velocities section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }     //}}}
    id--; // in lammps data file, ids start from 1
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Velocity[dd] = vel[dd];
    }
  } //}}}
  ChangeBoxByLow(System, -1);
  FillInCoor(System);
  return 1;
  error: // unrecognised line //{{{
    err_msg("unrecognised line in the file header");
    PrintErrorFileLine(file, *line_count);
    exit(1); //}}}
} //}}}
// read body //{{{
static void LmpDataReadBody(FILE *fr, const char *file, SYSTEM *System,
                            const int atom_types, int *line_count) {
  COUNT *Count = &System->Count;
  BEADTYPE *name_mass = calloc(atom_types, sizeof *name_mass);
  // create arrays for bonds/angles/dihedrals/impropers //{{{
  int(*bond)[3], (*angle)[4], (*dihedral)[5], (*improper)[5];
  if (Count->Bond > 0) {
    bond = calloc(Count->Bond, sizeof *bond);
  } else {
    bond = calloc(1, sizeof *bond);
  }
  if (Count->Angle > 0) {
    angle = calloc(Count->Angle, sizeof *angle);
  } else {
    angle = calloc(1, sizeof *angle);
  }
  if (Count->Dihedral > 0) {
    dihedral = calloc(Count->Dihedral, sizeof *dihedral);
  } else {
    dihedral = calloc(1, sizeof *dihedral);
  }
  if (Count->Improper > 0) {
    improper = calloc(Count->Improper, sizeof *improper);
  } else {
    improper = calloc(1, sizeof *improper);
  } //}}}
  // variables to check all sections are present
  bool atoms = false,                // Atoms section
       bonds[2] = {false, false},     // Bond Coeffs & Bonds sections
       angles[2] = {false, false},    // Angle Coeffs & Angles sections
       dihedrals[2] = {false, false}, // Dihedral Coeffs & Dihedrals sections
       impropers[2] = {false, false}; // Improper Coeffs & Impropers sections
  // read the rest of the file //{{{
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    (*line_count)++;
    // evaluate the line and read the appropriate section
    if (words > 0 && strcmp(split[0], "Masses") == 0) {
      LmpDataReadMasses(fr, file, name_mass, atom_types, line_count);
    } else if (words > 1 && strcmp(split[0], "Bond") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      bonds[0] = true;
      LmpDataReadBondCoeffs(fr, file, System, line_count);
    } else if (words > 1 && strcmp(split[0], "Angle") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      angles[0] = true;
      LmpDataReadAngleCoeffs(fr, file, System, line_count);
    } else if (words > 1 && strcmp(split[0], "Dihedral") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      dihedrals[0] = true;
      LmpDataReadDihedralCoeffs(fr, file, System, line_count);
    } else if (words > 1 && strcmp(split[0], "Improper") == 0 &&
               strcmp(split[1], "Coeffs") == 0) {
      impropers[0] = true;
      LmpDataReadImproperCoeffs(fr, file, System, line_count);
    } else if (words > 0 && strcmp(split[0], "Atoms") == 0) {
      atoms = true;
      int mode = 0; // for 'Atoms' or 'Atoms # full'
      if (words > 2 && split[1][0] == '#') {
        if (strcmp(split[2], "full") == 0) {
          mode = 0; // for 'Atoms # bond|angle|molecular'
        } else if (strcmp(split[2], "bond") == 0 ||
                   strcmp(split[2], "angle") == 0 ||
                   strcmp(split[2], "molecular") == 0) {
          mode = 1; // for 'Atoms # bond|angle|molecular'
        } else if (strcmp(split[2], "atomic") == 0) {
          mode = 2; // for 'Atoms # atomic'
        } else if (strcmp(split[2], "charge") == 0) {
          mode = 3; // for 'Atoms # charge'
        } else {    // for 'Atoms # <something else>'
          snprintf(ERROR_MSG, LINE, "unsupported atom_style (Atoms section) "
                   "%s%s; assuming 'full'", ErrYellow(), split[2]);
          PrintWarnFile(file, "\0", "\0");
        }
      }
      LmpDataReadAtoms(fr, file, System, name_mass, atom_types, line_count,
                       mode);
    } else if (words > 0 && strcmp(split[0], "Velocities") == 0) {
      LmpDataReadVelocities(fr, file, System, line_count);
    } else if (words > 0 && strcmp(split[0], "Bonds") == 0) {
      bonds[1] = true;
      LmpDataReadBonds(fr, file, *Count, bond, line_count);
    } else if (words > 0 && strcmp(split[0], "Angles") == 0) {
      angles[1] = true;
      LmpDataReadAngles(fr, file, *Count, angle, line_count);
    } else if (words > 0 && strcmp(split[0], "Dihedrals") == 0) {
      dihedrals[1] = true;
      LmpDataReadDihedrals(fr, file, *Count, dihedral, line_count);
    } else if (words > 0 && strcmp(split[0], "Impropers") == 0) {
      impropers[1] = true;
      LmpDataReadImpropers(fr, file, *Count, improper, line_count);
    }
  } //}}}
  free(name_mass);
  // errors/warnings - missing sections //{{{
  if (!atoms) {
    err_msg("Missing Atoms section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Bond > 0 && !bonds[0]) {
    err_msg("Missing Bond Coeffs; setting all bond parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->BondType; i++) {
      System->BondType[i].a = 0;
      System->BondType[i].b = 0;
      System->BondType[i].c = 0;
    }
  }
  if (Count->Bond > 0 && !bonds[1]) {
    err_msg("Missing Bonds section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Angle > 0 && !angles[0]) {
    err_msg("Missing Angle Coeffs; setting all angle parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->AngleType; i++) {
      System->AngleType[i].a = 0;
      System->AngleType[i].b = 0;
      System->AngleType[i].c = 0;
    }
  }
  if (Count->Angle > 0 && !angles[1]) {
    err_msg("Missing Angles section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Dihedral > 0 && !dihedrals[0]) {
    err_msg("Missing Dihedral Coeffs; setting all dihedral parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->DihedralType; i++) {
      System->DihedralType[i].a = 0;
      System->DihedralType[i].b = 0;
      System->DihedralType[i].c = 0;
    }
  }
  if (Count->Dihedral > 0 && !dihedrals[1]) {
    err_msg("Missing Dihedrals section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  }
  if (Count->Improper > 0 && !impropers[0]) {
    err_msg("Missing Improper Coeffs; setting all improper parameters to 0");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    for (int i = 0; i < Count->ImproperType; i++) {
      System->ImproperType[i].a = 0;
      System->ImproperType[i].b = 0;
      System->ImproperType[i].c = 0;
    }
  }
  if (Count->Improper > 0 && !impropers[1]) {
    err_msg("Missing Impropers section");
    PrintErrorFile(file, "\0", "\0");
    exit(1);
  } //}}}
  Count->MoleculeType = Count->Molecule;
  /*
   * Sort beads ascendingly in each MoleculeType (i.e., in each molecule);
   * necessary because bonds (etc.) are added according to their ascending bead
   * indexes, and no connction to the MoleculeType[].Bead array is made
   * TODO: I should rather fill Molecule[].Bead array and only later fill
   *       MoleculeType[].Bead based on connectivity of the molecule in
   *       question. Somehow...
   */
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    SortArrayInt(mt->Bead, mt->nBeads, 0);
  }
  CopyMoleculeTypeBeadsToMoleculeBeads(System);
  FillMTypeStuff(System, 0, 3, bond, Count->Bond);
  FillMTypeStuff(System, 1, 4, angle, Count->Angle);
  FillMTypeStuff(System, 2, 5, dihedral, Count->Dihedral);
  FillMTypeStuff(System, 3, 5, improper, Count->Improper);
  MinimizeMTypeStuffIds(System);
  free(bond);
  free(angle);
  free(dihedral);
  free(improper);
} //}}}
// read Masses section //{{{
static void LmpDataReadMasses(FILE *fr, const char *file, BEADTYPE *name_mass,
                              const int atom_types, int *line_count) {
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Masses' section");
    exit(1);
  }
  for (int i = 0; i < atom_types; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Masses' section");
      exit(1);
    }
    long type;
    double mass;
    if (words < 2 || !IsNaturalNumber(split[0], &type) || type > atom_types ||
        (!IsPosRealNumber(split[1], &mass) && strcmp(split[1], "???") != 0)) {
      err_msg("wrong line in Masses section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }
    type--;
    if (strcmp(split[1], "???") == 0) {
      err_msg("undefined mass (\'???\' in Masses section)");
      PrintWarnFileLine(file, *line_count);
      name_mass[type].Mass = MASS;
    } else {
      name_mass[type].Mass = mass;
    }
    // check for bead type name: # <name>
    if (words > 3 && split[2][0] == '#') {
      s_strcpy(name_mass[type].Name, split[3], BEAD_NAME);
    } else {
      snprintf(name_mass[type].Name, BEAD_NAME, "b%d", i);
    }
  }
} //}}}
// read Bond Coeffs section //{{{
// TODO: non-harmonic bonds
static void LmpDataReadBondCoeffs(FILE *fr, const char *file,
                                  SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // error - no bond types //{{{
  if (Count->BondType == 0) {
    err_msg("Bond Coeffs in a file with no bond types");
    PrintWarnFile(file, "\0", "\0");
    putc('\n', stderr);
    return;
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Bond Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->BondType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Bond Coeffs' section");
      exit(1);
    }
    long type;
    double a, b;
    // error - wrong line //{{{
    if (words < 3 ||
        !IsNaturalNumber(split[0], &type) ||
        type > Count->BondType ||
        !IsRealNumber(split[1], &a) || a < 0 ||
        !IsRealNumber(split[2], &b) || b < 0) {
      err_msg("wrong line in Bond Coeffs section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }                                     //}}}
    System->BondType[type - 1].a = a * 2; // lammps uses k/2 as spring constant
    System->BondType[type - 1].b = b;
  }
} //}}}
// read Angle Coeffs section //{{{
// TODO: non-harmonic angles
static void LmpDataReadAngleCoeffs(FILE *fr, const char *file,
                                   SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // error - no angle types //{{{
  if (Count->AngleType == 0) {
    err_msg("Angle Coeffs in a file with no angle types");
    PrintWarnFile(file, "\0", "\0");
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Angle Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->AngleType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Angle Coeffs' section");
      exit(1);
    }
    long type;
    double a, b;
    // error - wrong line //{{{
    if (words < 3 ||
        !IsNaturalNumber(split[0], &type) ||
        type > Count->AngleType ||
        !IsRealNumber(split[1], &a) || a < 0 ||
        !IsRealNumber(split[2], &b) || b < 0) {
      err_msg("wrong line in Angle Coeffs section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    }                                      //}}}
    System->AngleType[type - 1].a = a * 2; // lammps uses k/2 as spring constant
    System->AngleType[type - 1].b = b;
  }
} //}}}
// read Dihedral Coeffs section //{{{
// TODO: non-harmonic (lammps type) angles
static void LmpDataReadDihedralCoeffs(FILE *fr, const char *file,
                                      SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // error - no dihedral types //{{{
  if (Count->DihedralType == 0) {
    err_msg("Dihedral Coeffs in a file with no dihedral types");
    PrintWarnFile(file, "\0", "\0");
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Dihedral Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->DihedralType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Dihedral Coeffs' section");
      exit(1);
    }
    long type;
    double a;
    long b, c;
    // error - wrong line //{{{
    // see https://docs.lammps.org/dihedral_harmonic.html for details
    if (words < 4 || !IsNaturalNumber(split[0], &type) ||
        type > Count->DihedralType ||
        !IsRealNumber(split[1], &a) ||
        (!IsIntegerNumber(split[2], &b) && labs(b) != 1) ||
        !IsWholeNumber(split[3], &c)) {
      err_msg("wrong line in Dihedral Coeffs section "
                        "(for now, only harmonic type accepted)");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    } //}}}
    System->DihedralType[type-1].a = a * 2;
    System->DihedralType[type-1].b = b;
    System->DihedralType[type-1].c = c;
  }
} //}}}
// read Improper Coeffs section //{{{
// TODO: non-cvff (lammps type) impropers
static void LmpDataReadImproperCoeffs(FILE *fr, const char *file,
                                      SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // error - no improper types //{{{
  if (Count->ImproperType == 0) {
    err_msg("Improper Coeffs in a file with no improper types");
    PrintWarnFile(file, "\0", "\0");
  } //}}}
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Improper Coeffs' section");
    exit(1);
  }
  for (int i = 0; i < Count->ImproperType; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Improper Coeffs' section");
      exit(1);
    }
    long type;
    double a;
    long b, c;
    // error - wrong line //{{{
    // see https://docs.lammps.org/improper_cvff.html for details
    if (words < 4 || !IsNaturalNumber(split[0], &type) ||
        type > Count->DihedralType ||
        !IsRealNumber(split[1], &a) ||
        (!IsIntegerNumber(split[2], &b) && labs(b) != 1) ||
        !IsWholeNumber(split[3], &c)) {
      err_msg("wrong line in Improper Coeffs (only cvff type accepted)");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    } //}}}
    System->ImproperType[type - 1].a = a * 2;
    System->ImproperType[type - 1].b = b;
    System->ImproperType[type - 1].c = c;
  }
} //}}}
// read Atoms section //{{{
// mode: 0..full; 1..angle/bond; 2..atomic; 3..molecular; 4..charge
static void LmpDataReadAtoms(FILE *fr, const char *file, SYSTEM *System,
                             const BEADTYPE *name_mass, const int atom_types,
                             int *line_count, const int mode) {
  COUNT *Count = &System->Count;
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Atoms' section");
    exit(1);
  }
  bool warned = false; // to warn of undefined charge only once
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Atoms' section");
      exit(1);
    }
    long id, resid, type;
    double q;
    double pos[3];
    // check the line is correct & read/assign values //{{{
    // 'Atoms # full': <bead id> <mol id> <bead type id> <charge> <coordinates>
    if (mode == 0) {
      if (words < 7 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                   // bead index
          !IsIntegerNumber(split[1], &resid) || // molecule index
          !IsWholeNumber(split[2], &type) ||    // bead type (mass-defined)
          (!IsRealNumber(split[3], &q) &&       // bead charge
           strcmp(split[3], "???") != 0) ||     //
          !IsRealNumber(split[4], &pos[0]) ||    // Cartesean coordinates
          !IsRealNumber(split[5], &pos[1]) ||    //
          !IsRealNumber(split[6], &pos[2])) {    //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
    } else if (mode == 1) {
      // 'Atoms # bond|angle|molecular':
      // <bead id> <mol id> <bead type id> <coordinates>
      if (words < 6 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                   // bead index
          !IsIntegerNumber(split[1], &resid) || // molecule index
          !IsWholeNumber(split[2], &type) ||    // bead type (mass-defined)
          !IsRealNumber(split[3], &pos[0]) ||    // Cartesean coordinates
          !IsRealNumber(split[4], &pos[1]) ||    //
          !IsRealNumber(split[5], &pos[2])) {    //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
      q = CHARGE;
    } else if (mode == 2) {
      // 'Atoms # atomic': <bead id> <bead type id> <coordinates>
      if (words < 5 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsWholeNumber(split[1], &type) || // bead type (mass-defined)
          !IsRealNumber(split[2], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[3], &pos[1]) || //
          !IsRealNumber(split[4], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
      resid = -1;
      q = CHARGE;
    } else {
      // 'Atoms # charge': <bead id> <bead type id> <charge> <coordinates>
      if (words < 6 || !IsNaturalNumber(split[0], &id) ||
          id > Count->Bead ||                // bead index
          !IsWholeNumber(split[1], &type) || // bead type (mass-defined)
          (!IsRealNumber(split[2], &q) &&    // bead charge
           strcmp(split[3], "???") != 0) ||  //
          !IsRealNumber(split[3], &pos[0]) || // Cartesean coordinates
          !IsRealNumber(split[4], &pos[1]) || //
          !IsRealNumber(split[5], &pos[2])) { //
        err_msg("wrong line in Atoms section");
        PrintErrorFileLine(file, *line_count);
        exit(1);
      }
      resid = -1;
    }
    //}}}
    id--;   // in lammps data file, bead ids start from 1
    type--; // in lammps data file, bead type ids start from 1
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = pos[dd];
    }
    b->InTimestep = true;
    b->Type = id;
    /*
     * bead type is same as bead index but with lmp-defined mass, i.e., there
     * is the same number of beads and bead types; the single-bead bead types
     * will be merged later
     */
    BEADTYPE *bt = &System->BeadType[id];
    if (type > atom_types) {
      snprintf(ERROR_MSG, LINE, "bead type is too high; id is %s%ld%s in "
               "a file with %s%d%s atom types", ErrRed(), type, ErrYellow(),
               ErrRed(), atom_types, ErrYellow());
      PrintErrorFile(file, "\0", "\0");
      exit(1);
    }
    bt->Mass = name_mass[type].Mass;
    if (strcmp(split[3], "???") == 0) {
      if (!warned) {
        err_msg("undefined charge (\'???\' in Atoms section)");
        PrintWarnFile(file, "\0", "\0");
        putc('\n', stderr);
        warned = true;
      }
      bt->Charge = CHARGE;
    } else {
      bt->Charge = q;
    }
    bt->Number = 1;
    s_strcpy(bt->Name, name_mass[type].Name, BEAD_NAME);
    if (resid >= 0) { // bead in a molecule //{{{
      Count->Bonded++;
      b->Molecule = resid;
      // resid ids may be discontinuous, so define Molecule for all possible ids
      if (resid > Count->HighestResid) {
        System->MoleculeType = s_realloc(System->MoleculeType, (resid + 1) *
                                         sizeof *System->MoleculeType);
        System->Molecule = s_realloc(System->Molecule, (resid + 1) *
                                     sizeof *System->Molecule);
        for (int j = (Count->HighestResid + 1); j <= resid; j++) {
          InitMoleculeType(&System->MoleculeType[j]);
          InitMolecule(&System->Molecule[j]);
        }
        Count->HighestResid = resid; // goes from 0
        Count->Molecule = resid + 1;
      }
      // does the molecule type (i.e., that molecule id) exist?
      MOLECULETYPE *mt_resid = &System->MoleculeType[resid];
      if (mt_resid->Number == 0) { // no, create new type
        char name[MOL_NAME] = "m";
        for (int j = 0; j < words; j++) {
          if (split[j][0] == '#' && words > j) {
            snprintf(name, MOL_NAME, "%s", split[j+1]);
            mt_resid->Flag = true;
            break;
          }
        }
        s_strcpy(mt_resid->Name, name, MOL_NAME);
        mt_resid->Number = 1;
        mt_resid->nBeads = 1;
        mt_resid->Bead = malloc(sizeof *mt_resid->Bead);
        mt_resid->Bead[0] = id;
      } else { // yes, assign the atom to that type
        int bead = mt_resid->nBeads;
        mt_resid->nBeads++;
        mt_resid->Bead = s_realloc(mt_resid->Bead, mt_resid->nBeads *
                                   sizeof *mt_resid->Bead);
        mt_resid->Bead[bead] = id;
      }
    } else {
      Count->Unbonded++;
    } //}}}
  }
} //}}}
// read Velocities section //{{{
static void LmpDataReadVelocities(FILE *fr, const char *file,
                                  SYSTEM *System, int *line_count) {
  COUNT *Count = &System->Count;
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Velocities' section");
    exit(1);
  }
  for (int i = 0; i < Count->Bead; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Velocities' section");
      exit(1);
    }
    long id;
    double vel[3];
    // error - wrong line //{{{
    if (words < 4 || !IsNaturalNumber(split[0], &id) ||
        id > Count->Bead ||                // bead index
        !IsRealNumber(split[1], &vel[0]) || // bead velocities
        !IsRealNumber(split[2], &vel[1]) || //
        !IsRealNumber(split[3], &vel[2])) { //
      err_msg("wrong line in Velocities section");
      PrintErrorFileLine(file, *line_count);
      exit(1);
    } //}}}
    id--;
    BEAD *b = &System->Bead[id];
    for (int dd = 0; dd < 3; dd++) {
      b->Velocity[dd] = vel[dd];
    }
  }
} //}}}
// read Bonds section //{{{
static void LmpDataReadBonds(FILE *fr, const char *file, const COUNT Count,
                             int (*bond)[3], int *line_count) {
  bool *found = calloc(Count.Bond, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Bonds' section");
    exit(1);
  }
  bool warned = false;
  for (int i = 0; i < Count.Bond; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Bonds' section");
      exit(1);
    }
    long id, type, b_id[2];
    // errors //{{{
    // not <int> <int> <int> <int> line
    if (words < 4 || !IsNaturalNumber(split[0], &id) ||
        (!IsNaturalNumber(split[1], &type) && strcmp(split[1], "???") != 0) ||
        !IsNaturalNumber(split[2], &b_id[0]) ||
        !IsNaturalNumber(split[3], &b_id[1])) {
      err_msg("wrong line in Bonds section");
      goto error;
    }
    if (id > Count.Bond) {
      err_msg("Bond index is too high");
      goto error;
    }
    if (type > Count.BondType) {
      err_msg("Bond type is too high");
      goto error;
    }
    if (b_id[0] > Count.Bead || b_id[1] > Count.Bead) {
      err_msg("Bead index in a bond is too high");
      goto error;
    }
    if (b_id[0] == b_id[1]) {
      err_msg("Same beads in a bond");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated bond index");
      goto error;
    } //}}}
    bond[id-1][0] = b_id[0] - 1;
    bond[id-1][1] = b_id[1] - 1;
    if (strcmp(split[1], "???") == 0) {
      if (!warned) {
        err_msg("undefined bond type (\'???\' in Bonds section)");
        PrintWarnFile(file, "\0", "\0");
        warned = true;
      }
      bond[id-1][2] = -1;
    } else {
      bond[id-1][2] = type - 1;
    }
    found[id - 1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
// read Angles section //{{{
static void LmpDataReadAngles(FILE *fr, const char *file, const COUNT Count,
                              int (*angle)[4], int *line_count) {
  bool *found = calloc(Count.Angle, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Angles' section");
    exit(1);
  }
  for (int i = 0; i < Count.Angle; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Angles' section");
      exit(1);
    }
    long id, type, a_id[3];
    // errors //{{{
    // not <int> <int> <int> <int> line
    if (words < 5 || !IsNaturalNumber(split[0], &id) ||
        (!IsNaturalNumber(split[1], &type) && strcmp(split[1], "???") != 0) ||
        !IsNaturalNumber(split[2], &a_id[0]) ||
        !IsNaturalNumber(split[3], &a_id[1]) ||
        !IsNaturalNumber(split[4], &a_id[2])) { //
      err_msg("wrong line in Angles section");
      goto error;
    }
    if (id > Count.Angle) {
      err_msg("Angle index is too high");
      goto error;
    }
    if (type > Count.AngleType) {
      err_msg("Angle type is too high");
      goto error;
    }
    if (a_id[0] > Count.Bead || a_id[1] > Count.Bead || a_id[2] > Count.Bead) {
      err_msg("Bead index in an angle is too high");
      goto error;
    }
    if (a_id[0] == a_id[1] || a_id[0] == a_id[2] || a_id[1] == a_id[2]) {
      err_msg("Same beads in an angle");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated angle index");
      goto error;
    } //}}}
    angle[id-1][0] = a_id[0] - 1;
    angle[id-1][1] = a_id[1] - 1;
    angle[id-1][2] = a_id[2] - 1;
    angle[id-1][3] = type - 1;
    found[id-1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
// read Dihedrals section //{{{
static void LmpDataReadDihedrals(FILE *fr, const char *file, const COUNT Count,
                                 int (*dihedral)[5], int *line_count) {
  bool *found = calloc(Count.Dihedral, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Dihedrals' section");
    exit(1);
  }
  for (int i = 0; i < Count.Dihedral; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Dihedrals' section");
      exit(1);
    }
    long id, type, d_id[4];
    // errors //{{{
    // not <int> <int> <int> <int> line
    if (words < 6 || !IsNaturalNumber(split[0], &id) ||
        !IsNaturalNumber(split[1], &type) ||
        !IsNaturalNumber(split[2], &d_id[0]) ||
        !IsNaturalNumber(split[3], &d_id[1]) ||
        !IsNaturalNumber(split[4], &d_id[2]) ||
        !IsNaturalNumber(split[5], &d_id[3])) {
      err_msg("wrong line in Dihedrals section");
      goto error;
    }
    if (id > Count.Dihedral) {
      err_msg("Dihedral index is too high");
      goto error;
    }
    if (type > Count.DihedralType) {
      err_msg("Dihedral type is too high");
      goto error;
    }
    if (d_id[0] > Count.Bead ||
        d_id[1] > Count.Bead ||
        d_id[2] > Count.Bead ||
        d_id[3] > Count.Bead) {
      err_msg("Bead index in a dihedral is too high");
      goto error;
    }
    if (d_id[0] == d_id[1] || d_id[0] == d_id[2] || d_id[0] == d_id[3] ||
        d_id[1] == d_id[2] || d_id[1] == d_id[3] || d_id[2] == d_id[3]) {
      err_msg("Same beads in a single dihedral");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated dihedral index");
      goto error;
    } //}}}
    dihedral[id-1][0] = d_id[0] - 1;
    dihedral[id-1][1] = d_id[1] - 1;
    dihedral[id-1][2] = d_id[2] - 1;
    dihedral[id-1][3] = d_id[3] - 1;
    dihedral[id-1][4] = type - 1;
    found[id-1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
// read Impropers section //{{{
static void LmpDataReadImpropers(FILE *fr, const char *file, const COUNT Count,
                                 int (*improper)[5], int *line_count) {
  bool *found = calloc(Count.Angle, sizeof *found);
  // skip one line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "incomplete 'Impropers' section");
    exit(1);
  }
  for (int i = 0; i < Count.Improper; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete 'Impropers' section");
      exit(1);
    }
    long id, type, i_id[4];
    // errors //{{{
    // not 6x<int>
    if (words < 6 || !IsNaturalNumber(split[0], &id) ||
        !IsNaturalNumber(split[1], &type) ||
        !IsNaturalNumber(split[2], &i_id[0]) ||
        !IsNaturalNumber(split[3], &i_id[1]) ||
        !IsNaturalNumber(split[4], &i_id[2]) ||
        !IsNaturalNumber(split[5], &i_id[3])) {
      err_msg("wrong line in Impropers section");
      goto error;
    }
    if (id > Count.Improper) {
      err_msg("Improper index is too high");
      goto error;
    }
    if (type > Count.ImproperType) {
      err_msg("Improper type is too high");
      goto error;
    }
    if (i_id[0] > Count.Bead || i_id[1] > Count.Bead || i_id[2] > Count.Bead ||
        i_id[3] > Count.Bead) {
      err_msg("Bead index in a improper is too high");
      goto error;
    }
    if (i_id[0] == i_id[1] || i_id[0] == i_id[2] || i_id[0] == i_id[3] ||
        i_id[1] == i_id[2] || i_id[1] == i_id[3] || i_id[2] == i_id[3]) {
      err_msg("Same beads in a single improper");
      goto error;
    }
    if (found[id - 1]) {
      err_msg("Repeated improper index");
      goto error;
    } //}}}
    improper[id-1][0] = i_id[0] - 1;
    improper[id-1][1] = i_id[1] - 1;
    improper[id-1][2] = i_id[2] - 1;
    improper[id-1][3] = i_id[3] - 1;
    improper[id-1][4] = type - 1;
    found[id-1] = true;
  }
  free(found);
  return;
error:
  PrintErrorFileLine(file, *line_count);
  exit(1);
} //}}}
static bool ReadPbc(BOX *box) { //{{{
  double lo, hi;
  char *a[3] = {"xlo", "ylo", "zlo"};
  char *b[3] = {"xhi", "yhi", "zhi"};
  for (int i = 0; i < 3; i++) {
    if (words > 3 && strcmp(split[2], a[i]) == 0 &&
                     strcmp(split[3], b[i]) == 0) {
      if (!IsRealNumber(split[0], &lo) || !IsRealNumber(split[1], &hi)) {
        return false;
      }
      box->Low[i] = lo;
      box->OrthoLength[i] = hi - lo;
    }
  }
  // <double> <double> <double> xy xz yz
  if (words > 5 && strcmp(split[3], "xy") == 0 &&
             strcmp(split[4], "xz") == 0 && strcmp(split[5], "yz") == 0) {
    double xy, xz, yz;
    if (!IsRealNumber(split[0], &xy) || !IsRealNumber(split[1], &xz) ||
        !IsRealNumber(split[2], &yz)) {
      return false;
    }
    box->transform[0][1] = xy;
    box->transform[0][2] = xz;
    box->transform[1][2] = yz;
  }
  return true;
} //}}}

// write single bond/angle/dihedral/improper 'stuff' //{{{
void WriteStuff(FILE *fw, const SYSTEM System, const int mol,
                int *count, const int num, int (**arr)[num], const int n) {
  (*count)++;
  bool in = true;
  for (int aa = 0; aa < (num - 1); aa++) {
    int id = (*arr)[n][aa];
    id = System.Molecule[mol].Bead[id];
    if (!System.Bead[id].InTimestep) {
      in = false;
      break;
    }
  }
  if (in) {
    fprintf(fw, "%7d", *count);
    if ((*arr)[n][num-1] != -1) {
      fprintf(fw, " %6d", (*arr)[n][num-1] + 1);
    } else {
      fprintf(fw, "   ???");
    }
    for (int aa = 0; aa < (num - 1); aa++) {
      int id = (*arr)[n][aa];
      id = System.Molecule[mol].Bead[id];
      fprintf(fw, " %5d", id + 1);
    }
    putc('\n', fw);
  }
} //}}}
// write all bond/angle/dihedral/improper 'stuff'
void WriteAllStuff(FILE *fw, const SYSTEM System, const int type) { //{{{
  int num = 0;
  if (type == 0) { // bond
    num = 3;
    fprintf(fw, "\nBonds\n\n");
  } else if (type == 1) { // angle
    num = 4;
    fprintf(fw, "\nAngles\n\n");
  } else if (type == 2) { // dihedral
    num = 5;
    fprintf(fw, "\nDihedrals\n\n");
  } else /* if (type == 0) */ { // improper
    num = 5;
    fprintf(fw, "\nImpropers\n\n");
  }
  int count = 0;
  for (int i = 0; i < System.Count.Molecule; i++) {
    int mtype = System.Molecule[i].Type;
    MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
    int (**arr)[num];
    int *n;
    if (type == 0) { // bond
      arr = &mt_i->Bond;
      n = &mt_i->nBonds;
    } else if (type == 1) { // angle
      arr = &mt_i->Angle;
      n = &mt_i->nAngles;
    } else if (type == 2) { // dihedral
      arr = &mt_i->Dihedral;
      n = &mt_i->nDihedrals;
    } else /* if (type == 0) */ { // improper
      arr = &mt_i->Improper;
      n = &mt_i->nImpropers;
    }
    for (int j = 0; j < *n; j++) {
      WriteStuff(fw, System, i, &count, num, arr, j);
    }
  }
} //}}}
// WriteLmpData() //{{{
void WriteLmpData(const SYSTEM System, const char *file, const bool mass,
                  const int argc, char **argv) {
  FILE *fw = OpenFile(file, "w");
  fprintf(fw, "Created via AnalysisTools v%s ", VERSION);
  fprintf(fw, "(https://github.com/KaGaSi/AnalysisTools); Command: ");
  PrintCommand(fw, argc, argv);
  putc('\n', fw);
  const COUNT *Count = &System.Count;
  // create new SYSTEM structure to figure out bead types if mass == true //{{{
  int mass_types = 0;
  int *bt_masstype_to_old = calloc(Count->BeadType, sizeof *bt_masstype_to_old);
  int *bt_old_to_masstype = calloc(Count->BeadType, sizeof *bt_old_to_masstype);
  if (mass) {
    for (int i = 0; i < Count->BeadType; i++) {
      BEADTYPE *bt_i = &System.BeadType[i];
      bool found = false;
      for (int j = 0; j < mass_types; j++) {
        BEADTYPE *bt_j = &System.BeadType[bt_masstype_to_old[j]];
        if (bt_i->Mass == bt_j->Mass) {
          // bt_masstype_to_old[j] = i;
          bt_old_to_masstype[i] = j;
          found = true;
          break;
        }
      }
      if (!found) {
        bt_masstype_to_old[mass_types] = i;
        bt_old_to_masstype[i] = mass_types;
        mass_types++;
      }
    }
  } else {
    mass_types = Count->BeadType;
    for (int i = 0; i < Count->BeadType; i++) {
      bt_old_to_masstype[i] = i;
      bt_masstype_to_old[i] = i;
    }
  } //}}}
  // print counts //{{{
  fprintf(fw, "%10d atoms\n", Count->Bead);
  fprintf(fw, "%10d bonds\n", Count->Bond);
  fprintf(fw, "%10d angles\n", Count->Angle);
  fprintf(fw, "%10d dihedrals\n", Count->Dihedral);
  fprintf(fw, "%10d impropers\n", Count->Improper);
  putc('\n', fw);
  // // add one atom type for extra (possibly srp)
  // fprintf(fw, "%10d atom types\n", Count->BeadType + 1);
  fprintf(fw, "%10d atom types\n", mass_types);
  if (Count->Bond > 0 && Count->BondType > 0) {
    fprintf(fw, "%10d bond types\n", Count->BondType);
  }
  if (Count->Angle > 0 && Count->AngleType > 0) {
    fprintf(fw, "%10d angle types\n", Count->AngleType);
  }
  if (Count->Dihedral > 0 && Count->DihedralType > 0) {
    fprintf(fw, "%10d dihedral types\n", Count->DihedralType);
  }
  if (Count->Improper > 0 && Count->ImproperType > 0) {
    fprintf(fw, "%10d improper types\n", Count->ImproperType);
  }
  putc('\n', fw); //}}}
  // print box size //{{{
  fprintf(fw, "%.3f", System.Box.Low[0]);
  fprintf(fw, " %.3f xlo xhi\n", System.Box.Low[0] + System.Box.OrthoLength[0]);
  fprintf(fw, "%.3f", System.Box.Low[1]);
  fprintf(fw, " %.3f ylo yhi\n", System.Box.Low[1] + System.Box.OrthoLength[1]);
  fprintf(fw, "%.3f", System.Box.Low[2]);
  fprintf(fw, " %.3f zlo zhi\n", System.Box.Low[2] + System.Box.OrthoLength[2]);
  if (System.Box.alpha != 90 ||
      System.Box.beta != 90 ||
      System.Box.gamma != 90) {
    fprintf(fw, "%.3f %.3f %.3f xy xz yz\n", System.Box.transform[0][1],
                                             System.Box.transform[0][2],
                                             System.Box.transform[1][2]);
  }
  putc('\n', fw); //}}}
  // print bead type masses //{{{
  fprintf(fw, "Masses\n\n");
  for (int i = 0; i < mass_types; i++) {
    BEADTYPE *bt = &System.BeadType[bt_masstype_to_old[i]];
    fprintf(fw, "%5d", i + 1);
    if (bt->Mass == MASS) {
      fprintf(fw, " ???");
    } else {
      fprintf(fw, " %lf", bt->Mass);
    }
    fprintf(fw, " # %s\n", bt->Name);
  } //}}}
  // print various coeffs //{{{
  if (Count->BondType > 0) {
    fprintf(fw, "\nBond Coeffs\n\n");
    for (int i = 0; i < Count->BondType; i++) {
      fprintf(fw, "%5d %lf %lf\n", i + 1, System.BondType[i].a / 2,
                                          System.BondType[i].b);
    }
  }
  if (Count->AngleType > 0) {
    fprintf(fw, "\nAngle Coeffs\n\n");
    for (int i = 0; i < Count->AngleType; i++) {
      fprintf(fw, "%5d %lf %lf\n", i + 1, System.AngleType[i].a / 2,
                                          System.AngleType[i].b);
    }
  }
  if (Count->DihedralType > 0) {
    fprintf(fw, "\nDihedral Coeffs\n\n");
    for (int i = 0; i < Count->DihedralType; i++) {
      fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.DihedralType[i].a / 2,
                                              System.DihedralType[i].b,
                                              System.DihedralType[i].c);
    }
  }
  if (Count->ImproperType > 0) {
    fprintf(fw, "\nImproper Coeffs\n\n");
    for (int i = 0; i < Count->ImproperType; i++) {
      fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.ImproperType[i].a / 2,
                                              System.ImproperType[i].b,
                                              System.ImproperType[i].c);
    }
  } //}}}
  // print atoms //{{{
  // if there is 0 molecule index, saved indices will get +1
  // TODO: why would I need to go from 1?
  bool zero = false;
  // for (int i = 0; i < Count->Molecule; i++) {
  //   if (System.Molecule[i].Index == 0) {
  //     zero = true;
  //     break;
  //   }
  // }
  fprintf(fw, "\nAtoms # full\n\n");
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    // <bead id>
    fprintf(fw, "%7d", id + 1);
    // <molecule id (-1 for no molecule)>
    int mol = bead->Molecule;
    if (mol != -1) {
      int id = System.Molecule[mol].Index;
      if (zero) {
        id++;
      }
      fprintf(fw, " %5d", id);
    } else {
      fprintf(fw, " %5d", -1);
    }
    // <bead type id>
    fprintf(fw, " %5d", bt_old_to_masstype[bead->Type] + 1);
    // <charge> from original System as the charge can differ in lmp data file
    int type = bead->Type;
    double q = System.BeadType[type].Charge;
    if (q == CHARGE) {
      fprintf(fw, " %15s", "???");
    } else {
      fprintf(fw, " %15f", q);
    }
    // coordinates
    for (int dd = 0; dd < 3; dd++) {
      fprintf(fw, " %15f", bead->Position[dd] + System.Box.Low[dd]);
    }
    // molecule name
    if (mol != -1) {
      int type = System.Molecule[mol].Type;
      fprintf(fw, " # %s", System.MoleculeType[type].Name);
    }
    putc('\n', fw);
  } //}}}
  // print velocities (if at least one non-zero) //{{{
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    double (*vel)[3] = &System.Bead[id].Velocity;
    if (fabs((*vel)[0]) > 1e-5 ||
        fabs((*vel)[1]) > 1e-5 ||
        fabs((*vel)[2]) > 1e-5) {
      fprintf(fw, "\nVelocities\n\n");
      for (int j = 0; j < Count->BeadCoor; j++) {
        id = System.BeadCoor[j];
        vel = &System.Bead[id].Velocity;
        fprintf(fw, "%7d", id + 1);
        for (int dd = 0; dd < 3; dd++) {
          fprintf(fw, " %15f", (*vel)[dd]);
        }
        putc('\n', fw);
      }
      break;
    }
  } //}}}
  // print bonds/angles/dihedrals/impropers
  if (Count->Bond > 0) {
    WriteAllStuff(fw, System, 0);
  }
  if (Count->Angle > 0) {
    WriteAllStuff(fw, System, 1);
  }
  if (Count->Dihedral > 0) {
    WriteAllStuff(fw, System, 2);
  }
  if (Count->Improper > 0) {
    WriteAllStuff(fw, System, 3);
  }
  free(bt_masstype_to_old);
  free(bt_old_to_masstype);
  fclose(fw);
} //}}}
