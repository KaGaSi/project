#include "General.h"

// IsDouble() //{{{
/**
 * Function to test if provided string is a real number.
 */
bool IsDouble(char *a) {
  // wrong first character - can be minus, dot, or number
  if (a[0] != '-' && a[0] != '.' && (a[0] < '0' || a[0] > '9')) {
    return false;
  }
  // only one dot can be present
  bool dot = false;
  if (a[0] == '.') {
    dot = true;
  }
  // test the remaining characters - either digit, or dot (but only 1 in total)
  for (int i = 1; i < strlen(a); i++) {
    if (a[i] == '.') {
      if (dot) { // has there been a dot already?
        return false;
      } else {
        dot = true;
      }
    } else if (a[i] < '0' || a[i] > '9') {
      return false;
    }
  }
  return true;
} //}}}

// IsPosDouble() //{{{
/**
 * Function to test if provided string is a non-negative real number.
 */
bool IsPosDouble(char *a) {
  // wrong first character - can be minus, dot, or number
  if (a[0] != '.' && (a[0] < '0' || a[0] > '9')) {
    return false;
  }
  // only one dot can be present
  bool dot = false;
  if (a[0] == '.') {
    dot = true;
  }
  // test the remaining characters - either digit, or dot (but only 1 in total)
  for (int i = 1; i < strlen(a); i++) {
    if (a[i] == '.') {
      if (dot) { // has there been a dot already?
        return false;
      } else {
        dot = true;
      }
    } else if (a[i] < '0' || a[i] > '9') {
      return false;
    }
  }
  return true;
} //}}}

// IsInteger() //{{{
/**
 * Function to test if provided string is a non-negative whole number.
 */
bool IsInteger(char *a) {
  // test the remaining characters - either digit, or dot (but only 1 in total)
  for (int i = 0; i < strlen(a); i++) {
    if (a[i] < '0' || a[i] > '9') {
      return false;
    }
  }
  return true;
} //}}}

// Length() //{{{
/**
 * Function to calculate vector length.
 */
double Length(Vector a) {
  double length = sqrt(SQR(a.x) + SQR(a.y) + SQR(a.z));
  return length;
} //}}}

// Min3() //{{{
/**
 * Function returning the lowest number from three floats.
 */
double Min3(double x, double y, double z) {

  double min;
  if (x > y) {
    if (y > z) {
      min = z;
    } else {
      min = y;
    }
  } else if (x > z) {
    min = z;
  } else {
    min = x;
  }

  return min;
} //}}}

// Max3() //{{{
/**
 * Function returning the highest number from three floats.
 */
double Max3(double x, double y, double z) {

  double max;
  if (x < y) {
    if (y < z) {
      max = z;
    } else {
      max = y;
    }
  } else if (x < z) {
    max = z;
  } else {
    max = x;
  }

  return max;
} //}}}

// Sort3() //{{{
/**
 * Function returning sorted numbers x < y < z.
 */
Vector Sort3(Vector in) {

  Vector out;

  if (in.x < in.y) {
    if (in.y < in.z) {
      out.x = in.x;
      out.y = in.y;
      out.z = in.z;
    } else if (in.x < in.z) {
      out.x = in.x;
      out.y = in.z;
      out.z = in.y;
    } else {
      out.x = in.z;
      out.y = in.x;
      out.z = in.y;
    }
  } else {
    if (in.x < in.z) {
      out.x = in.y;
      out.y = in.x;
      out.z = in.z;
    } else if (in.y < in.z) {
      out.x = in.y;
      out.y = in.z;
      out.z = in.x;
    } else {
      out.x = in.z;
      out.y = in.y;
      out.z = in.x;
    }
  }

  return out;
} //}}}

// Swap() //{{{
/**
 * Swap two integers.
 */
void Swap(int *a, int *b) {
  int swap = *a;
  *a = *b;
  *b = swap;
}
// }}}

// SwapDouble() //{{{
/**
 * Swap two doubles.
 */
void SwapDouble(double *a, double *b) {
  double swap = *a;
  *a = *b;
  *b = swap;
}
// }}}

// SwapBool() //{{{
/**
 * Swap two booleans.
 */
void SwapBool(bool *a, bool *b) {
  bool swap = *a;
  *a = *b;
  *b = swap;
}
// }}}

// SortArray() //{{{
/**
 * Sort an array using the bubble sort algorithm. If mode = 0, sort
 * ascendingly; if mode = 1, sort descendingly.
 */
void SortArray(int **array, int length, int mode) {
  if (mode != 0 && mode != 1) {
    fprintf(stderr, "\nError - SortArray(): use 0 or 1 for sorting mode\n");
    exit(1);
  }
  for (int i = 0 ; i < (length-1); i++) {
    bool done = true;
    for (int j = 0 ; j < (length-i-1); j++) {
      if (mode == 0 && (*array)[j] > (*array)[j+1]) {
        Swap(&(*array)[j], &(*array)[j+1]);
        done = false;
      }
      if (mode == 1 && (*array)[j] < (*array)[j+1]) {
        Swap(&(*array)[j], &(*array)[j+1]);
        done = false;
      }
    }
    if (done)
      break;
  }
} //}}}

// Distance() //{{{
/**
 * Function calculating distance vector between two beads. It removes
 * periodic boundary conditions and returns x, y, and z distances in the
 * range <0, BoxLength/2).
 */
Vector Distance(Vector id1, Vector id2, Vector BoxLength) {

  Vector rij;

  // distance vector
  rij.x = id1.x - id2.x;
  rij.y = id1.y - id2.y;
  rij.z = id1.z - id2.z;

  // remove periodic boundary conditions in x-direction
  while (rij.x >= (BoxLength.x/2))
    rij.x = rij.x - BoxLength.x;
  while (rij.x < -(BoxLength.x/2))
    rij.x = rij.x + BoxLength.x;
  // in y-direction
  while (rij.y >= (BoxLength.y/2))
    rij.y = rij.y - BoxLength.y;
  while (rij.y < -(BoxLength.y/2))
    rij.y = rij.y + BoxLength.y;
  // in z-direction
  while (rij.z >= (BoxLength.z/2))
    rij.z = rij.z - BoxLength.z;
  while (rij.z < -(BoxLength.z/2))
    rij.z = rij.z + BoxLength.z;

  return rij;
} //}}}

// SplitLine() //{{{
/**
 * Function that splits the provided line into individual strings (using tab,
 * space, and colon as a delimiter) and removes newline character from the end
 * of the last string.
 */
int SplitLine(char out[30][100], char *line) {
  // trim whitespaces at the beginning and end of line
  strcpy(line, TrimLine(line));
  // split into words separated by " ", tab, or colon
  char *split[30];
  split[0] = strtok(line, " \t:"); // first word
  int words = 0;
  while (split[words] != NULL && words < 29) {
    words++; // start from 1, as the first split is already done
    split[words] = strtok(NULL, " \t:");
  }
  // if the last word ends with newline, make it into '\0'
  if (split[words-1][strlen(split[words-1])-1] == '\n') {
    split[words-1][strlen(split[words-1])-1] = '\0';
  }
  if (split[words-1][0] == '\n') {
    words--;
  }
  // if the last words is just '\0', disregard it
  for (int i = 0; i < words; i++) {
    strcpy(out[i], split[i]);
  }
 return words;
} //}}}

// TrimLine //{{{
/**
 * Function to trim whitespace from the
 * beginning and end of a string.
 */
char * TrimLine(char *line) {
  int length = strlen(line);
  static char trimmed[LINE];
  strcpy(trimmed, line);
  // 1) trailing whitespace
  while (length > 1 &&
         (trimmed[length-1] == ' ' ||
          trimmed[length-1] == '\n' ||
          trimmed[length-1] == '\t')) {
    trimmed[length-1] = '\0';
    length--;
  }
  // 2) preceding whitespace
  while (length > 1 &&
         (trimmed[0] == ' ' ||
          trimmed[0] == '\n' ||
          trimmed[0] == '\t')) {
    for (int i = 0; i < length; i++) { // line[length] contains '\0'
      trimmed[i] = trimmed[i+1];
    }
    length--;
  }
  // 3) if only 1 character remains, make it into newline (otherwise segfault happens)
  if (length == 1 && (trimmed[0] == ' ' || trimmed[0] == '\t')) {
    trimmed[0] = '\n';
  }

  return trimmed;
} //}}}
