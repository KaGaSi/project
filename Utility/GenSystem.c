#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stderr, "\033[1;33m");
    fprintf(stdout, "\
Warning: This utility was not extensively tested and is in fact not a \
very good generator of initial configuration.\n\n");
    fprintf(stderr, "\033[0m");
    fprintf(stdout, "\
GenSystem reads information from a FIELD-like file to create vsf \
structure file and generate coordinates for all beads (used, e.g., as \
initial configuration for a simulation). This utility only creates linear \
molecules no matter the connectivity in the provided FIELD-like file. \
GenSystem also lacks checking for errors in the FIELD-like file so if an \
incorrect file is provided, the utility will exhibit undefined behaviour \
(it will freeze, crash, or produce wrong results).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <out.vsf> <out.vcf> <options>\n\n", cmd);
  fprintf(ptr, "   <out.vsf>     output structure file (vsf format)\n");
  fprintf(ptr, "   <out.vcf>     output coordinate file (vcf format)\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -f <name>  FIELD-like file (default: FIELD)\n");
  fprintf(ptr, "      -v         verbose output\n");
  fprintf(ptr, "      -h         print this help and exit\n");
  fprintf(ptr, "      --version  print version number and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 2; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // options before reading system data //{{{
  // output verbosity //{{{
  bool verbose = BoolOption(argc, argv, "-v"); // verbose output
  // }}}

  // FIELD-like file //{{{
  char *input = calloc(LINE, sizeof(char *));
  if (FileOption(argc, argv, "-f", &input)) {
    exit(1);
  }
  if (input[0] == '\0') {
    strcpy(input, "FIELD");
  } //}}}
  //}}}

  count = 0; // count arguments

  // <out.vsf> - output structure file (must end with .vsf) //{{{
  char output[LINE];
  strcpy(output, argv[++count]);

  // test if <output.vsf> filename ends with '.vsf' (required by VMD)
  int ext = 1;
  char extension[1][5];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output vcf file //{{{
  char *output_vcf = calloc(LINE, sizeof(char));
  char stuff[LINE];
  strcpy(output_vcf, argv[++count]);

  // test if outpuf_vcf has '.vcf' extension - required by vmd //{{{
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}
  //}}}

  // variables - structures //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  VECTOR BoxLength;
  int *Index;
  PARAMS *bond_type; // information about bond types
  PARAMS *angle_type; // information about angle types //}}}

  ReadField(input, &BoxLength, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule, &bond_type, &angle_type);

  // set all molecules & beads to 'write' //{{{
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = true;
  }
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = true;
  } //}}}

  // create & fill output vsf file
  WriteVsf(output, Counts, BeadType, Bead, MoleculeType, Molecule, false);

  // print information - verbose option //{{{
  if (verbose) {
    char null[1] = {'\0'};
    fprintf(stdout, "\n\n");
    VerboseOutput(null, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // Generate coordinates //{{{
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(input, "r")) == NULL) {
    ErrorFileOpen(input, 'r');
    exit(1);
  } //}}}

  char line[LINE], split[30][100], delim[8];
  strcpy(delim, " \t");
  // skip to molecule section //{{{
  while(fgets(line, sizeof(line), fr)) {
    SplitLine(split, line, delim);
    if (strncmp(split[0], "molecules", 8) == 0 ||
        strncmp(split[0], "Molecules", 8) == 0 ||
        strncmp(split[0], "MOLECULES", 8) == 0 ) {
      Counts.TypesOfMolecules = atoi(split[1]);
      break;
    }
  } //}}}

  // read bond lenghts for all molecule types and create prototype molecules //{{{
  printf("\nbox = (%lf, %lf, %lf)\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  double *prototype_z[Counts.TypesOfMolecules]; // prototype molecule z coordinate
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    // allocate array for bond lengths
    prototype_z[i] = calloc(MoleculeType[i].nBeads,sizeof(double));
    prototype_z[i][0] = 0;

    // skip to bonds section //{{{
    while(fgets(line, sizeof(line), fr)) {
      SplitLine(split, line, delim);
      if (strncmp(split[0], "bonds", 4) == 0 ||
          strncmp(split[0], "Bonds", 4) == 0 ||
          strncmp(split[0], "BONDS", 4) == 0 ) {
        break;
      }
    } //}}}

//  for (int j = 0; j < MoleculeType[i].nBonds; j++) {
    // if we assume linear molecule, there must be nBeads-1 bonds
    for (int j = 0; j < (MoleculeType[i].nBeads-1); j++) {
      fgets(line, sizeof(line), fr);
      SplitLine(split, line, delim);
      double length = atof(split[4]); // bond length is in fifth columne
      if (fabs(length) < 0.01) {
        length = 0.7;
      }
      prototype_z[i][j+1] = prototype_z[i][j] + length;
    }
  } //}}}
  fclose(fr);

  // longest molecule type //{{{
  int longest = -1;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (MoleculeType[i].nBeads > longest) {
      longest = i;
    }
  }
  double length = prototype_z[longest][MoleculeType[longest].nBeads-1];
  printf("Length of the molecule: %lf\n", length); //}}}

  double dist = 0.7; // distance between layers (and beads)
  int x_n = BoxLength.x / dist, y_n = BoxLength.y / dist, z_n = BoxLength.z / dist; // number of positions per axes
  int x = 0, y = 0; // coordinates that are incrementally increased
  double z = 0.1; // coordinates that are incrementally increased
  int count_free = 0; // count number of unbonded beads that are placed

  int layer[Counts.TypesOfMolecules];
  int total_layers = 0;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    layer[i] = x_n * (int)(BoxLength.z / (prototype_z[i][MoleculeType[i].nBeads-1] + dist));
    printf("Maximum number of %s molecules per layer: %d\n", MoleculeType[i].Name, layer[i]);
    layer[i] = MoleculeType[i].Number / layer[i] + 1;
    total_layers += layer[i];
  }
  int sol_layers = Counts.Unbonded / (x_n * z_n);
  int free_layer_per_mol_layer = sol_layers / total_layers;
  printf("%d molecule layers and %d solvent layers (%d solvent layers per one molecule layer).\n", total_layers, sol_layers, free_layer_per_mol_layer);

  // molecules layers with solvent layers in between //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int id = Molecule[i].Bead[j];
      Bead[id].Position.x = (x % x_n) * dist + 0.1;
      Bead[id].Position.y = (y % y_n) * dist + 0.1;
      Bead[id].Position.z = z + prototype_z[type][j];
    }

    z += prototype_z[type][MoleculeType[type].nBeads-1] + dist;
    // add solvent when molecule do not fit till the end of z direction
    if ((z+prototype_z[type][MoleculeType[type].nBeads-1]) >= BoxLength.z) {
      for (; count_free < Counts.Unbonded; count_free++) {
        Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
        Bead[count_free].Position.y = (y % y_n) * dist + 0.1;
        Bead[count_free].Position.z = z;
        z += dist;
        if (z >= BoxLength.z) {
          break;
        }
      }
      z = 0.1;
    }
    // if z is filled, change also x and y
    if (z == 0.1) {
      x++;
      int y_old = y;
      y = x / x_n;
      // when y changes, whole xz layer of molecules was added -- add xz layers of unbonded beads
      if ((y-y_old) > 0) {
        for (int j = 0; j < free_layer_per_mol_layer; j++) {
          for (; count_free < Counts.Unbonded; count_free++) {
            Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
            Bead[count_free].Position.y = (y % y_n) * dist + 0.1;
            Bead[count_free].Position.z = z;
            z += dist;
            if (z >= BoxLength.z) {
              z = 0.1;
            }
            if (z == 0.1) {
              x++;
              y_old = y;
              y = x / x_n;
              // when y changes, whole xz layer of unbonded beads was added
              if ((y-y_old) > 0) {
                break;
              }
            }
          }
        }
      }
    }
  } //}}}

  // remaining unbonded beads //{{{
  for (; count_free < Counts.Unbonded; count_free++) {
    Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
    Bead[count_free].Position.y = (y % y_n) * dist + 0.1;
    Bead[count_free].Position.z = z;
    z += dist;
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z == 0.1) {
      x++;
      y = x / x_n;
    }
  } //}}}

  // tweak y_dist and do it again //{{{
  double y_dist = ((BoxLength.y - 0.5) / (y * dist)) * dist;
  x = 0; // coordinates that are incrementally increased
  y = 0; // coordinates that are incrementally increased
  z = 0.1; // coordinates that are incrementally increased
  double y_coor = 0.1;
  count_free = 0; // count number of unbonded beads that are placed

  // molecules layers with solvent layers in between //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int id = Molecule[i].Bead[j];
      Bead[id].Position.x = (x % x_n) * dist + 0.1;
      Bead[id].Position.y = y_coor;
      Bead[id].Position.z = z + prototype_z[type][j];
    }

    z += prototype_z[type][MoleculeType[type].nBeads-1] + dist;
    // add solvent when molecule do not fit till the end of z direction
    if ((z+prototype_z[type][MoleculeType[type].nBeads-1]) >= BoxLength.z) {
      for (; count_free < Counts.Unbonded; count_free++) {
        Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
        Bead[count_free].Position.y = y_coor;
        Bead[count_free].Position.z = z;
        z += dist;
        if (z >= BoxLength.z) {
          break;
        }
      }
      z = 0.1;
    }
    // if z is filled, change also x and y
    if (z == 0.1) {
      x++;
      int y_old = y;
      y = x / x_n;
      // when y changes, whole xz layer of molecules was added -- add xz layers of unbonded beads
      if ((y-y_old) > 0) {
        y_coor += y_dist;
        for (int j = 0; j < free_layer_per_mol_layer; j++) {
          for (; count_free < Counts.Unbonded; count_free++) {
            Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
            Bead[count_free].Position.y = y_coor;
            Bead[count_free].Position.z = z;
            z += dist;
            if (z >= BoxLength.z) {
              z = 0.1;
            }
            if (z == 0.1) {
              x++;
              y_old = y;
              y = x / x_n;
              // when y changes, whole xz layer of unbonded beads was added
              if ((y-y_old) > 0) {
                y_coor += y_dist;
                break;
              }
            }
          }
        }
      }
    }
  } //}}}

  // remaining unbonded beads //{{{
  for (; count_free < Counts.Unbonded; count_free++) {
    Bead[count_free].Position.x = (x % x_n) * dist + 0.1;
    Bead[count_free].Position.y = y_coor;
    Bead[count_free].Position.z = z;
    z += dist;
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z >= BoxLength.z) {
      z = 0.1;
    }
    if (z == 0.1) {
      x++;
      int y_old = y;
      y = x / x_n;
      if ((y-y_old) > 0) {
        y_coor += y_dist;
      }
    }
  } //}}}
  //}}}
  //}}}

  // open output .vcf file for writing //{{{
  FILE *out;
  if ((out = fopen(output_vcf, "w")) == NULL) {
    ErrorFileOpen(output_vcf, 'w');
    exit(1);
  } //}}}

  strcpy(stuff, "# by GenSystem\n");

  // write pbc
  fprintf(out, "pbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);
  // write coordinates
  WriteCoorIndexed(out, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);
  fclose(out);

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    free(prototype_z[i]);
  }
  free(input);
  free(output_vcf);
  free(angle_type);
  free(bond_type); //}}}

  return 0;
}
