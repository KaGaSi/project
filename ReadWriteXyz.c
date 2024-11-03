#include "ReadWriteXyz.h"

/*
 * Functions to read xyz file as a coordinate file via XyzReadTimestep() and
 * XyzSkipTimestep() and as a structure file via XyzReadStruct()
 */
static bool XyzCheckCoorLine(double coor[3]);

SYSTEM XyzReadStruct(const char *file) { //{{{
  SYSTEM Sys;
  InitSystem(&Sys);
  COUNT *Count = &Sys.Count;
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  // read number of beads //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing number of beads");
    exit(1);
  }
  long val;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("wrong first line of an xyz file");
    PrintErrorFileLine(file, line_count);
    exit(1);
  } //}}}
  Count->Bead = val;
  Count->BeadCoor = val;
  Count->Unbonded = val;
  Count->UnbondedCoor = val;
  Sys.Bead = s_realloc(Sys.Bead, sizeof *Sys.Bead * Count->Bead);
  Sys.BeadCoor = s_realloc(Sys.BeadCoor, sizeof *Sys.BeadCoor * Count->Bead);
  // read next line //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing comment line");
    exit(1);
  } //}}}
  // read atoms //{{{
  for (int i = 0; i < Count->Bead; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete timestep");
      exit(1);
    }
    double coor[3];
    if (!XyzCheckCoorLine(coor)) {
      err_msg("wrong coordinate line");
      PrintErrorFileLine(file, line_count);
      exit(1);
    }
    BEAD *b = &Sys.Bead[i];
    InitBead(b);
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = coor[dd];
    }
    b->InTimestep = true;
    Sys.BeadCoor[i] = i;
    // assign bead type or create a new one
    bool new = true;
    for (int j = 0; j < Count->BeadType; j++) {
      BEADTYPE *bt = &Sys.BeadType[j];
      if (strcmp(split[0], bt->Name) == 0) {
        bt->Number++;
        b->Type = j;
        new = false;
        break;
      }
    }
    if (new) {
      int type = Count->BeadType;
      NewBeadType(&Sys.BeadType, &Count->BeadType, split[0], HIGHNUM, HIGHNUM, HIGHNUM);
      BEADTYPE *bt_new = &Sys.BeadType[type];
      bt_new->Number = 1;
      b->Type = type;
    }
  } //}}}
  fclose(fr);
  FillSystemNonessentials(&Sys);
  for (int i = 0; i < Count->Bead; i++) {
    Sys.UnbondedCoor[i] = i;
  }
  CheckSystem(Sys, file);
  return Sys;
} //}}}
// XYZReadTimestep() //{{{
int XyzReadTimestep(FILE *fr, const char *file,
                    SYSTEM *System, int *line_count) {
  System->Count.MoleculeCoor = 0;
  // read number of beads //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  long val;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("wrong first line of an xyz timestep");
    PrintWarnFileLine(file, *line_count);
    return -1;
  } else if (val > System->Count.Bead) {
    snprintf(ERROR_MSG, LINE,
             "too many beads in the timestep (maximum number is %s%d%s)",
             ErrYellow(), System->Count.Bead, ErrRed());
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  System->Count.BeadCoor = val; //}}}
  // ignore next line //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing comment line");
    return -2;
  }
  (*line_count)++; //}}}
  // read coordinates //{{{
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      snprintf(ERROR_MSG, LINE, "premature end of file "
               "(%s%d%s coordinate lines instead of %s%d%s)", ErrYellow(),
               i, ErrRed(), ErrYellow(), System->Count.BeadCoor, ErrRed());
      PrintErrorFileLine(file, *line_count);
      return -2;
    }
    double coor[3];
    if (XyzCheckCoorLine(coor)) {
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[i].Position[dd] = coor[dd];
      }
      System->BeadCoor[i] = i;
      (*line_count)++;
      double vel[3];
      if (words > 6 &&
          IsRealNumber(split[4], &vel[0]) &&
          IsRealNumber(split[5], &vel[1]) &&
          IsRealNumber(split[6], &vel[2])) {
        for (int dd = 0; dd < 3; dd++) {
          System->Bead[i].Velocity[dd] = vel[dd];
        }
      } else {
        for (int dd = 0; dd < 3; dd++) {
          System->Bead[i].Velocity[dd] = 0;
        }
      }
    } else {
      err_msg("unrecognized line (not enough valid coordinate lines)");
      PrintErrorFileLine(file, *line_count);
      return -1;
    }
  } //}}}
  FillInCoor(System);
  return true;
} //}}}
bool XyzSkipTimestep(FILE *fr, const char *file, int *line_count) { //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  long val;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("wrong first line of an xyz timestep");
    PrintWarnFileLine(file, *line_count);
    return -1;
  }
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  for (int i = 0; i < val; i++) {
    (*line_count)++;
    char a;
    while ((a = getc(fr) != '\n') && a != EOF)
      ;
  }
  return true;
} //}}}
static bool XyzCheckCoorLine(double coor[3]) { //{{{
  if (words > 3 &&
      IsRealNumber(split[1], &coor[0]) &&
      IsRealNumber(split[2], &coor[1]) &&
      IsRealNumber(split[3], &coor[2])) {
    return true;
  } else {
    return false;
  }
} //}}}

void XyzWriteCoor(FILE *fw, const bool *write, const SYSTEM System) { //{{{
  // find out number of beads to save
  int count = 0;
  bool none = true; // to make sure there are beads to save
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    if (write[id]) {
      none = false;
      count++;
    }
  }
  if (none) {
    err_msg("no beads to save");
    PrintWarning();
    return;
  }
  // write pbc on the second line
  fprintf(fw, "%d\n", count);
  const BOX *box = &System.Box;
  if (box->Volume != -1) {
    fprintf(fw, "%.3f %.3f %.3f", box->Length[0],
                                  box->Length[1],
                                  box->Length[2]);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 90) {
      fprintf(fw, " %lf %lf %lf", box->alpha, box->beta, box->gamma);
    }
  }
  putc('\n', fw);
  // write the coodinates
  // for (int i = 0; i < System.Count.BeadCoor; i++) {
  //   int id = System.BeadCoor[i];
  for (int id = 0; id < System.Count.Bead; id++) {
    BEAD *bead = &System.Bead[id];
    if (write[id] && bead->InTimestep) {
      int type = bead->Type;
      fprintf(fw, "%8s %8.4f %8.4f %8.4f\n", System.BeadType[type].Name,
              bead->Position[0], bead->Position[1], bead->Position[2]);
    }
  }
} //}}}
