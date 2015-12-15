#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "CStructs.h"
#include "Structure.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   SelectedVcf <input.vsf> <start> <skip> ");
  fprintf(stderr, "<output.vcf> <type names> <options>\n\n");

  fprintf(stderr, "   <input.vcf>     input filename (vcf format)\n");
  fprintf(stderr, "   <start>         number of timestep to start from\n");
  fprintf(stderr, "   <skip>          leave out every 'skip' steps\n");
  fprintf(stderr, "   <output.vcf>    output filename (vcf format)\n");
  fprintf(stderr, "   <type names>    name of bead types to save\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -v   verbose output\n");
  fprintf(stderr, "      -h   print this help and exit\n");
  fprintf(stderr, "      -i   use input .vsf file different from dl_meso.vsf\n\n");
} //}}}

int main(int argc, char *argv[]) {

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  printf("\n\n"); //}}}

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("SelectedVcf creates new <output.vcf> file from <input.vcf>\n");
      printf("containing only selected bead types. Also <start> timesteps\n");
      printf("can be omitted and every <skip> timesteps can be left out.\n");
      printf("The program uses dl_meso.vsf (or other input structure file)\n");
      printf("and FIELD files to determine all information about the system.\n\n");

      printf("Usage:\n");
      printf("   SelectedVcf <input.vcf> <start> <skip> ");
      printf("<output.vcf> <type names> <options>\n\n");

      printf("   <input.vcf>     input filename (vcf format)\n");
      printf("   <start>         number of timestep to start from\n");
      printf("   <skip>          leave out every 'skip' steps\n");
      printf("   <output.vcf>    output filename (vcf format)\n");
      printf("   <type names>    name of bead types to save\n");
      printf("   <options>\n");
      printf("      -v   verbose output\n");
      printf("      -h   print this help and exit\n");
      printf("      -i   use input .vsf file different from dl_meso.vsf\n\n");
      exit(0);
    }
  } //}}}

  // check if correct number of arguments //{{{
  if (argc < 6) {
    fprintf(stderr, "Too little arguments!\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  int count = 0; // count arguments

  // <input.vcf> - filename of input vcf file (must end with .vcf) //{{{
  char input_vcf[32];
  strcpy(input_vcf, argv[++count]);

  // test if <input.vcf> filename ends with '.vsf' (required by VMD)
  char *dot = strrchr(input_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<input.vcf> '%s' does not have .vcf ending!\n", input_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // <start> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <start>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int start = atoi(argv[count]); //}}}

  // <skip> - number of steps to skip per one used //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <skip>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int skip = atoi(argv[count]); //}}}

  // <output.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char output_vcf[32];
  strcpy(output_vcf, argv[++count]);

  // test if <output.vcf> filename ends with '.vcf' (required by VMD)
  dot = strrchr(output_vcf, '.');
  if (!dot || strcmp(dot, ".vcf")) {
    fprintf(stderr, "<output.vcf> '%s' does not have .vcf ending!\n", output_vcf);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  //TODO: first read beadtype from .vcf (or take all if All.vcf); read into struct BeadType Selected
  //TODO: read into struct Counts whatever
  // initial malloc for array holding ids of selected beadtypes
  int *selected;
  selected = malloc(Counts.TypesOfBeads*sizeof(int));

  // selected beadtypes to save //{{{
  while (++count < argc || argv[count][0] != '-') { // go through remaing arguments not starting with '-'

    int type = FindType(argv[count], Counts, Selected);
    selected[] = type;

    // Error - beadtype name doesn't exist //{{{
    if (j == 0) {
      fprintf(stderr, "Beadtype %s does not exist!\n", argv[i]);
      exit(1);
    } //}}}
  }
  if (selected[0] == TypesOfBeads ) {
    fprintf(stderr, "All beadtypes are selected which would only copy All.vcf,\n");
    fprintf(stderr, "therefore this is not allowed!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }//}}}

  // -v option - verbose output //{{{
  count++;
  bool verbose = false;
  for (int i = count; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      verbose = true;

      break;
    }
  } //}}}

  // -i option - file name of input structure file //{{{
  char input_vsf[32];
  input_vsf[0] = '\0'; // check if -i option is used
  for (int i = count; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "Missing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");
        exit(1);
      }

      // check if .vsf ending is present
      dot = strrchr(argv[i+1], '.');
      if (!dot || strcmp(dot, ".vsf")) {
        fprintf(stderr, "'-i' argument '%s' does not have .vsf ending!\n", argv[i+1]);
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(input_vsf, argv[i+1]);
    }
  }

  // -i option is not used
  if (input_vsf[0] == '\0') {
    strcpy(input_vsf, "dl_meso.vsf");
  } //}}}

  if (verbose) {
    printf("%d %d\n", skip, start);
  }

  return 0;
}
