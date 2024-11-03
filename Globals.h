#ifndef CONSTANTS_H
#define CONSTANTS_H

extern const double PI;

enum { SPL_STR = 32, SPL_LEN = 64, LINE = 1024 };

extern char line[LINE];
extern char *split[SPL_STR];
extern int words;

extern char ERROR_MSG[LINE];

extern const char *BLACK;
extern const char *RED;
extern const char *GREEN;
extern const char *YELLOW;
extern const char *BLUE;
extern const char *MAGENTA;
extern const char *CYAN;
extern const char *WHITE;
extern const char *C_RESET;
#endif
