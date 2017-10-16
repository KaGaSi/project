#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "../AnalysisTools.h"
#include "../Options.h"
#include "Aggregates.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <distance> <contacts> ", cmd);
  fprintf(stderr, "<output.agg> <type names> <options>\n\n");

  fprintf(stderr, "   <input.vcf>       input filename (vcf format)\n");
  fprintf(stderr, "   <distance>        minimum distance for contact for aggregate check\n");
  fprintf(stderr, "   <contacts>        minimum number of contacts for aggregate check\n");
  fprintf(stderr, "   <output.agg>      output filename (agg format)\n");
  fprintf(stderr, "   <type names>      names of bead types to use for closeness calculation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "     -x <name(s)>    exclude specified molecule(s)\n");
  fprintf(stderr, "     -j <joined.vcf> output vcf file with joined coordinates\n");
  CommonHelp(1);
} //}}}

// CalculateAggregates() //{{{
/**
 * Function to determine distribution of molecules in aggregates.
 */
void CalculateAggregates(Aggregate **Aggregate, Counts *Counts, int sqdist, int contacts,
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule **Molecule) {

  // zeroize //{{{
  (*Counts).Aggregates = 0;
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Aggregate)[i].nMolecules = 0;
    (*Aggregate)[i].nBeads = 0;
    (*Aggregate)[i].nMonomers = 0;
  }

  for (int i = 0; i < ((*Counts).Bonded+(*Counts).Unbonded); i++) {
    (*Bead)[i].nAggregates = 0;
  } //}}}

  // allocate & zeroize contact[][] (triangular matrix) and moved array //{{{
  int **contact = malloc((*Counts).Molecules*sizeof(int *));
  int *moved = malloc((*Counts).Molecules*sizeof(int));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    contact[i] = malloc((i+1)*sizeof(int));
  }

  // zeroize arrays
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Molecule)[i].Aggregate = -1;
    moved[i] = 0;
    for (int j = 0; j < i; j++) /* j == i is really not needed */
      contact[i][j] = 0;
  } //}}}

  // cell size for cell linked list
  double cell_size = sqrt(sqdist);

  // number of cells in all three dimensions //{{{
  IntVector n_cells;
  n_cells.x = ceil(BoxLength.x/cell_size),
  n_cells.y = ceil(BoxLength.y/cell_size),
  n_cells.z = ceil(BoxLength.z/cell_size); //}}}

  // allocate memory for arrays for cell linked list
  int *Head = malloc((n_cells.x*n_cells.y*n_cells.z)*sizeof(int));
  int *Link = malloc(((*Counts).Unbonded+(*Counts).Bonded)*sizeof(int));

  // initialize Head array //{{{
  for (int i = 0; i < (n_cells.x*n_cells.y*n_cells.z); i++) {
    Head[i] = -1;
  } //}}}

  // sort beads into cells //{{{
  for (int i = 0; i < ((*Counts).Unbonded+(*Counts).Bonded); i++) {
    int cell = (int)((*Bead)[i].Position.x / cell_size)
             + (int)((*Bead)[i].Position.y / cell_size) * n_cells.x
             + (int)((*Bead)[i].Position.z / cell_size) * n_cells.x * n_cells.y;
    Link[i] = Head[cell];
    Head[cell] = i;
  } //}}}

  // coordinates of adjoining cells //{{{
  int Dcx[14] = {0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
  int Dcy[14] = {0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
  int Dcz[14] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}; //}}}

  // count contacts between all molecules pairs (using cell linked list) //{{{
  for (int c1z = 0; c1z < n_cells.z; c1z++) {
    for (int c1y = 0; c1y < n_cells.y; c1y++) {
      for (int c1x = 0; c1x < n_cells.x; c1x++) {

        // select first cell
        int cell1 = c1x + c1y * n_cells.x + c1z * n_cells.x * n_cells.y;

        // select first bead in the cell 'cell1'
        int i = Head[cell1];

        while (i != -1) {
          for (int k = 0; k < 14; k++) {
            int c2x = c1x + Dcx[k]; //{{{
            int c2y = c1y + Dcy[k];
            int c2z = c1z + Dcz[k]; //}}}

            // periodic boundary conditions for cells //{{{
            if (c2x >= n_cells.x)
              c2x -= n_cells.x;
            else if (c2x < 0)
              c2x += n_cells.x;

            if (c2y >= n_cells.y)
              c2y -= n_cells.y;
            else if (c2y < 0)
              c2y += n_cells.y;

            if (c2z >= n_cells.z)
              c2z -= n_cells.z; //}}}

            // select second cell
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y;

            // select bead in the cell 'cell2' //{{{
            int j;
            if (cell1 == cell2) { // next bead in 'cell1'
              j = Link[i];
            } else { // first bead in 'cell2'
              j = Head[cell2];
            } //}}}

            while (j != -1) {
              int type_i = (*Molecule)[(*Bead)[i].Molecule].Type;
              int type_j = (*Molecule)[(*Bead)[j].Molecule].Type;

              if (BeadType[(*Bead)[i].Type].Use && // bead i must be of specified type
                  BeadType[(*Bead)[j].Type].Use && // bead j as well
                  MoleculeType[type_i].Use && // molecule with i cannot be excluded
                  MoleculeType[type_j].Use) { // molecule with j the same

                // calculate distance between i and j beads
                Vector rij = Distance((*Bead)[i].Position, (*Bead)[j].Position, BoxLength);

                // are 'i' and 'j' close enough?
                if ((*Bead)[i].Molecule != (*Bead)[j].Molecule &&
                    (SQR(rij.x) + SQR(rij.y) + SQR(rij.z)) <= sqdist) {
                  if (i > j) {
                    contact[(*Bead)[i].Molecule][(*Bead)[j].Molecule]++;
                  } else {
                    contact[(*Bead)[j].Molecule][(*Bead)[i].Molecule]++;
                  }
                }
              }

              j = Link[j];
            }
          }
          i = Link[i];
        }
      }
    }
  } //}}}

  // evaluate the contacts //{{{
  // first molecule
  for (int i = 1; i < (*Counts).Molecules; i++) {

    // second molecule
    for (int j = 0; j < i; j++) {

      int agg_i = (*Molecule)[i].Aggregate;
      int agg_j = (*Molecule)[j].Aggregate;

      // molecules 'i' and 'j' are in contact //{{{
      if (contact[i][j] >= contacts) {
        // create new aggregate if 'j' isn'it in any //{{{
        if (agg_j == -1) {
          agg_j = (*Counts).Aggregates;

          (*Aggregate)[agg_j].nMolecules = 1;
          (*Aggregate)[agg_j].Molecule[0] = j;

          (*Molecule)[j].Aggregate = agg_j;

          (*Counts).Aggregates++;
        } //}}}

        // add 'i' to aggregate if 'i' isn't in any //{{{
        if (agg_i == -1) {
          (*Aggregate)[agg_j].Molecule[(*Aggregate)[agg_j].nMolecules] = i;

          (*Molecule)[i].Aggregate = agg_j;

          (*Aggregate)[agg_j].nMolecules++;
        } //}}}

        // each residue in different aggregate => unite aggregates
        if (agg_i != -1 && agg_j != -1 && agg_i != agg_j) {

          // add molecules from aggregate 'agg_i' to 'agg_j' //{{{
          int mols = (*Aggregate)[agg_j].nMolecules;

          (*Aggregate)[agg_j].nMolecules += (*Aggregate)[agg_i].nMolecules;
          int id1 = 0;
          // copy molecule ids from Aggregate[agg_i-1] to Aggregate[agg_j-1] */
          for (int k = mols; k < (*Aggregate)[agg_j].nMolecules; k++) {
            (*Aggregate)[agg_j].Molecule[k] = (*Aggregate)[agg_i].Molecule[id1];

            (*Molecule)[(*Aggregate)[agg_j].Molecule[k]].Aggregate = agg_j;

            id1++;
          } //}}}

          // move aggregates with id greater then agg_i to id-1 //{{{
          for (int k = (agg_i+1); k < (*Counts).Aggregates; k++) {

            (*Aggregate)[k-1].nMolecules = (*Aggregate)[k].nMolecules;

            // move every molecule from aggregate 'k' to aggregate 'k-1'
            for (int l = 0; l < (*Aggregate)[k].nMolecules; l++) {
              (*Aggregate)[k-1].Molecule[l] = (*Aggregate)[k].Molecule[l];

              (*Molecule)[(*Aggregate)[k-1].Molecule[l]].Aggregate = k - 1;
            }
          } //}}}

          // reduce number of aggregates (two aggregates were merged)
          (*Counts).Aggregates--;
        } //}}}
      // or 'i' and 'j' aren't in contact and 'j' isn't in any aggregate =>  new aggregate for 'j' */ //{{{
      } else if (agg_j == -1) {
        (*Aggregate)[(*Counts).Aggregates].nMolecules = 1;
        (*Aggregate)[(*Counts).Aggregates].Molecule[0] = j;

        (*Molecule)[j].Aggregate = (*Counts).Aggregates;

        (*Counts).Aggregates++;
      } //}}}
    }
  } //}}}

  // if residue with highest id is in no aggregate, create it //{{{
  // check if highest id residue is in aggregate //{{{
  bool test = false;
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    for (int j = 1; j < (*Aggregate)[i].nMolecules; j++) {
      if ((*Aggregate)[i].Molecule[j] == ((*Counts).Molecules-1)) {
        test = 1;
      }
    }
  } //}}}

  /* highest id residue not in any aggregate => create separate one */ //{{{
  if (!test) {
    (*Aggregate)[(*Counts).Aggregates].nMolecules = 1;
    (*Aggregate)[(*Counts).Aggregates].Molecule[0] = (*Counts).Molecules - 1;

    (*Counts).Aggregates++;
  } //}}}
  //}}}

  // bubble sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {

    for (int j = 0 ; j < ((*Aggregate)[i].nMolecules-1); j++) {
      bool done = true;

      for (int k = 0 ; k < ((*Aggregate)[i].nMolecules-j-1); k++) {

        if ((*Aggregate)[i].Molecule[k] > (*Aggregate)[i].Molecule[k+1]) {

          int swap = (*Aggregate)[i].Molecule[k];
          (*Aggregate)[i].Molecule[k] = (*Aggregate)[i].Molecule[k+1];
          (*Aggregate)[i].Molecule[k+1] = swap;

          done = false;
        }
      }
      if (done)
        break;
    }
  } //}}}

  // bubble sort aggregates according to ascending ids of first molecules //{{{
  for (int i = 0; i < ((*Counts).Aggregates-1); i++) {
    bool done = true;

    for (int j = 0; j < ((*Counts).Aggregates-i-1); j++) {

      if ((*Aggregate)[j].Molecule[0] > (*Aggregate)[j+1].Molecule[0]) {
        // swtich numbers of molecules
        int swap = (*Aggregate)[j].nMolecules;
        (*Aggregate)[j].nMolecules = (*Aggregate)[j+1].nMolecules;
        (*Aggregate)[j+1].nMolecules = swap;

        // switch the whole Aggregate[].Molecule array (no idea which aggregate contains more molecules)
        for (int k = 0; k < (*Counts).Molecules; k++) {
          swap = (*Aggregate)[j].Molecule[k];
          (*Aggregate)[j].Molecule[k] = (*Aggregate)[j+1].Molecule[k];
          (*Aggregate)[j+1].Molecule[k] = swap;
        }

        done = false;
      }
    }
    if (done)
      break;
  } //}}}

  // assign bonded beads to Aggregate struct //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {

    // go through all molecules in aggregate 'i'
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];

      // copy all bead in molecule 'mol' to Aggregate struct
      for (int k = 0; k < MoleculeType[(*Molecule)[mol].Type].nBeads; k++) {
        (*Aggregate)[i].Bead[(*Aggregate)[i].nBeads] = (*Molecule)[mol].Bead[k];
        (*Aggregate)[i].nBeads++;

        // every bead from molecule is only in one aggregate
        (*Bead)[(*Molecule)[mol].Bead[k]].nAggregates = 1;
        (*Bead)[(*Molecule)[mol].Bead[k]].Aggregate[0] = i;
      }
    }
  } //}}}

  // find monomeric beads close to aggregates (using cell linked list) //{{{
  for (int c1z = 0; c1z < n_cells.z; c1z++) {
    for (int c1y = 0; c1y < n_cells.y; c1y++) {
      for (int c1x = 0; c1x < n_cells.x; c1x++) {

        // select first cell
        int cell1 = c1x + c1y * n_cells.x + c1z * n_cells.x * n_cells.y;

        int i = Head[cell1];

        while (i != -1) {
          for (int k = 0; k < 14; k++) {
            int c2x = c1x + Dcx[k];
            int c2y = c1y + Dcy[k];
            int c2z = c1z + Dcz[k];

            // cell periodic boundary condition //{{{
            if (c2x >= n_cells.x)
              c2x -= n_cells.x;
            else if (c2x < 0)
              c2x += n_cells.x;

            if (c2y >= n_cells.y)
              c2y -= n_cells.y;
            else if (c2y < 0)
              c2y += n_cells.y;

            if (c2z >= n_cells.z)
              c2z -= n_cells.z; //}}}

            // select second cell
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y;

            int j;
            if (cell1 == cell2) {
              j = Link[i];
            } else {
              j = Head[cell2];
            }

            while (j != -1) {
              // test if the monmeric bead is near aggregate
              if ((*Bead)[i].Molecule == -1 && // monomeric 'i'
                  (*Bead)[j].Molecule != -1) { // 'j' in molecule //{{{

                int agg_j = (*Bead)[j].Aggregate[0];

                // test if 'i' is already in 'j''s aggregate //{{{
                bool in_agg = false;
                for (int k = 0; k < (*Bead)[i].nAggregates; k++) {
                  if ((*Bead)[i].Aggregate[k] == agg_j) {
                    in_agg = true;
                    break;
                  }
                } //}}}

                if (!in_agg) {
                  // calculate distance between i and j beads
                  Vector rij = Distance((*Bead)[i].Position, (*Bead)[j].Position, BoxLength);

                  // test if 'i' is near 'j''s aggregate
                  if ((SQR(rij.x)+SQR(rij.y)+SQR(rij.z)) < sqdist) {
                    (*Aggregate)[agg_j].Monomer[(*Aggregate)[agg_j].nMonomers] = (*Bead)[i].Index;
                    (*Aggregate)[agg_j].nMonomers++;

                    (*Bead)[i].Aggregate[(*Bead)[i].nAggregates] = agg_j;
                    (*Bead)[i].nAggregates++;
                  }
                } //}}}
              } else if ((*Bead)[j].Molecule == -1 && // monomeric 'j'
                         (*Bead)[i].Molecule != -1) { // 'i' in molecule //{{{

                int agg_i = (*Bead)[i].Aggregate[0];

                // test if 'j' is already in 'i''s aggregate //{{{
                bool in_agg = false;
                for (int k = 0; k < (*Bead)[j].nAggregates; k++) {
                  if ((*Bead)[j].Aggregate[k] == agg_i) {
                    in_agg = true;
                    break;
                  }
                } //}}}

                if (!in_agg) {
                  // calculate distance between i and j beads
                  Vector rij = Distance((*Bead)[i].Position, (*Bead)[j].Position, BoxLength);

                  // test if 'j' is near 'i''s aggregate
                  if ((SQR(rij.x)+SQR(rij.y)+SQR(rij.z)) < sqdist) {
                    (*Aggregate)[agg_i].Monomer[(*Aggregate)[agg_i].nMonomers] = (*Bead)[j].Index;
                    (*Aggregate)[agg_i].nMonomers++;

                    (*Bead)[j].Aggregate[(*Bead)[j].nAggregates] = agg_i;
                    (*Bead)[j].nAggregates++;
                  }
                }
              } //}}}

              j = Link[j];
            }
          }
          i = Link[i];
        }
      }
    }
  } //}}}

  // bubble sort monomers in aggregates according to ascending ids //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {

    for (int j = 0 ; j < ((*Aggregate)[i].nMonomers-1); j++) {
      bool done = true;

      for (int k = 0 ; k < ((*Aggregate)[i].nMonomers-j-1); k++) {

        if ((*Aggregate)[i].Monomer[k] > (*Aggregate)[i].Monomer[k+1]) {

          int swap = (*Aggregate)[i].Monomer[k];
          (*Aggregate)[i].Monomer[k] = (*Aggregate)[i].Monomer[k+1];
          (*Aggregate)[i].Monomer[k+1] = swap;

          done = false;
        }
      }
      if (done)
        break;
    }
  } //}}}

  // free memory //{{{
  free(Head);
  free(Link);
  for (int i = 0; i < (*Counts).Molecules; i++)
    free(contact[i]);
  free(contact);
  free(moved); //}}}
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("Aggregates utility determines which molecules belong to which aggregate on  \n");
      printf("the basis of given parameters - the minimum distance at which a pair of     \n");
      printf("beads from different molecules is considered a contact and the minimum      \n");
      printf("number of such contacts between two molecules to consider them as belonging \n");
      printf("to the same aggregate. Only distances between specified bead types are      \n");
      printf("considered. Information about aggregates in each timestep is written to     \n");
      printf(".agg file. Also joined coordinates can be written to an output .vcf file.   \n\n");

      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <distance> <contacts> ", argv[0]);
      printf("<output.agg> <type names> <options>\n\n");

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <distance>          minimum distance for contact for aggregate check\n");
      printf("   <contacts>          minimum number of contacts for aggregate check\n");
      printf("   <output.agg>        output filename (agg format)\n");
      printf("   <type names>        names of bead types for closeness calculation\n");
      printf("   <options>\n");
      printf("      -x <name(s)>     exclude specified molecule(s)\n");
      printf("      -j <joined.vcf>  output vcf file with joined coordinates\n");
      CommonHelp(0);
      exit(0);
    }
  }

  int options = 5; //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-b") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "-V") != 0 &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-j") != 0) {

      fprintf(stderr, "Non-existent option '%s'!\n", argv[i]);
      ErrorHelp(argv[0]);
      exit(1);
    }
  } //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  for (int i = 1; i < argc && argv[count+1][0] != '-'; i++) {
    count++;
  }

  if (count < options) {
    fprintf(stderr, "Too little mandatory arguments (%d instead of at least %d)!\n\n", count, options);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

//// -j <out.vcf> - filename of output vcf file (must end with .vcf) //{{{
//char joined_vcf[32];
//joined_vcf[0] = '\0'; // no -j option
//for (int i = 1; i < argc; i++) {
//  if (strcmp(argv[i], "-j") == 0) {

//    // wrong argument to -i option
//    if ((i+1) >= argc || argv[i+1][0] == '-') {
//      fprintf(stderr, "Missing argument to '-j' option ");
//      fprintf(stderr, "(or filename beginning with a dash)!\n");
//      exit(1);
//    }

//    strcpy(joined_vcf, argv[i+1]);

//    // test if <joined.vcf> filename ends with '.vcf' (required by VMD)
//    char *dot = strrchr(joined_vcf, '.');
//    if (!dot || strcmp(dot, ".vcf")) {
//      fprintf(stderr, "<joined.vcf> '%s' does not have .vcf ending!\n", joined_vcf);
//      ErrorHelp(argv[0]);
//      exit(1);
//    }
//  }
//} //}}}
  char joined_vcf[32];
  bool error = JoinCoorOption(argc, argv, joined_vcf); //{{{
  if (error) {
    exit(1);
  } //}}}

  // standard options //{{{
  char *vsf_file = calloc(32,sizeof(char *));
  char *bonds_file = calloc(32,sizeof(char *));
  bool verbose, verbose2, silent, script;
  error = CommonOptions(argc, argv, &vsf_file, &bonds_file, &verbose, &verbose2, &silent, &script);

  // was there error during CommonOptions()?
  if (error) {
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      printf(" %s", argv[i]);
    printf("\n\n");
  } //}}}

  count = 0; // count mandatory arguments

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

  // <distance> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <distance>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double distance = atof(argv[count]); //}}}

  // <contacts> - number of steps to skip per one used //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <contacts>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  int contacts = atoi(argv[count]); //}}}

  // <output.agg> - filename of output agg file (must end with .agg) //{{{
  char output_agg[32];
  strcpy(output_agg, argv[++count]);

  // test if <output.agg> filename ends with '.agg' (required by VMD)
  dot = strrchr(output_agg, '.');
  if (!dot || strcmp(dot, ".agg")) {
    fprintf(stderr, "<output.agg> '%s' does not have .agg ending!\n", output_agg);
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(vsf_file, input_vcf, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(vsf_file);

  // <type names> - names of bead types to use for closeness calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", argv[count], input_vcf);
      exit(1);
    }

    BeadType[type].Use = true;
  } //}}}

  // '-x' option //{{{
  error = ExcludeOption(argc, argv, Counts, &MoleculeType);
  if (error) {
    exit(1);
  } //}}}

  // print command to output .agg file //{{{
  FILE *out;
  if ((out = fopen(output_agg, "w")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_agg);
    exit(1);
  }

  // print command to stdout
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out);

  fclose(out); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf);
    exit(1);
  } //}}}

  // get pbc from coordinate file //{{{
  char str[32];
  // skip till 'pbc' keyword
  do {
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read string from '%s' file!\n", input_vcf);
    }
  } while (strcmp(str, "pbc") != 0);

  // read pbc
  Vector BoxLength;
  if (fscanf(vcf, "%lf %lf %lf", &BoxLength.x, &BoxLength.y, &BoxLength.z) != 3) {
    fprintf(stderr, "Cannot read pbc from %s!\n", input_vcf);
    exit(1);
  }

  // skip remainder of pbc line
  while (getc(vcf) != '\n')
    ;
  // skip blank line
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    printf("   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // write bead type names and pbc to <joined.vcf> if '-j' option was used //{{{
  if (joined_vcf[0] != '\0') {

    // bead types are to be written in joined.vcf -- probably will never change,
    // since the beadtypes should correspond to those in .agg file
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      BeadType[i].Write = true;
    }

    // open <joined.vcf>
    FILE *joined;
    if ((joined = fopen(joined_vcf, "w")) == NULL) {
      fprintf(stderr, "Cannot open output %s vcf file!\n", joined_vcf);
      exit(1);
    }

    // write bead type names
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      // only those bead types that are to be used
      if (BeadType[i].Write) {
        fprintf(joined, "# %s\n", BeadType[i].Name);
      }
    }

    fprintf(joined, "\npbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);

    fclose(joined);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = calloc(Counts.Molecules,sizeof(*Aggregate));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded,sizeof(int));
    // assumes all bonded beads can be in one aggregate - memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded,sizeof(int));
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules,sizeof(int));
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(verbose2, input_vcf, bonds_file, Counts, BeadType, Bead, MoleculeType, Molecule);

    printf("\n   Distance for closeness check: %lf\n", distance);
    printf("   Number of needed contacts for aggregate check: %d\n", contacts);
  }

  // bonds file is not needed anymore
  free(bonds_file); //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int test;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;
    if (!silent) {
      if (script) {
        printf("Step: %6d\n", count);
      } else {
        fflush(stdout);
        printf("\rStep: %6d", count);
      }
    }

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf, Counts, &Bead, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
        // print newline to stdout if Step... doesn't end with one
        if (!script && !silent) {
          putchar('\n');
        }
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      }
    } //}}}

    RestorePBC(Counts, BoxLength, &Bead);

    CalculateAggregates(&Aggregate, &Counts, SQR(distance), contacts, BoxLength, BeadType, &Bead, MoleculeType, &Molecule);

    // calculate & write joined coordinatest to <joined.vcf> if '-j' option is used //{{{
    if (joined_vcf[0] != '\0') {

      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

      // open <joined.vcf> file //{{{
      FILE *joined;
      if ((joined = fopen(joined_vcf, "a")) == NULL) {
        fprintf(stderr, "Cannot open output %s vcf file!\n", joined_vcf);
        exit(1);
      } //}}}

      WriteCoorIndexed(joined, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

      fclose(joined);
    } //}}}

    // open output .agg file for appending //{{{
    if ((out = fopen(output_agg, "a")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", output_agg);
      exit(1);
    } //}}}

    // write data to output .agg file //{{{

    // find the number of aggregates - remove aggregates only of excluded mols //{{{
    int no_excluded_aggs = 0;
    int test_count = 0; // to test that every molecule is in an aggregate
    for (int i = 0; i < Counts.Aggregates; i++) {
      Aggregate[i].Use = false;

      test_count += Aggregate[i].nMolecules;

      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
        if (MoleculeType[moltype].Use) {
          Aggregate[i].Use = true;
          no_excluded_aggs++;
          break;
        }
      }
    } //}}}

    // are all molecules accounted for? //{{{
    if (test_count != Counts.Molecules) {
      fprintf(stderr, "Not all molecules were assigned to aggregates!\n");
      fprintf(stderr, "Counts.Molecules = %5d; Molecules in aggregates: %d\n", Counts.Molecules, test_count);
      exit(1);
    } //}}}

    // print number of aggregates to agg file
    fprintf(out, "\nStep: %d\n%d\n\n", count, no_excluded_aggs);

    // go through all aggregates
    for (int i = 0; i < Counts.Aggregates; i++) {

      if (Aggregate[i].Use) {

        // go through all molecules in aggregate 'i'
        fprintf(out, "%d :", Aggregate[i].nMolecules);
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          fprintf(out, " %d", Aggregate[i].Molecule[j]+1);
        }
        putc('\n', out);

        // go through all monomeric beads in aggregate 'i'
        fprintf(out, "   %d :", Aggregate[i].nMonomers);
        for (int j = 0; j < Aggregate[i].nMonomers; j++) {
          fprintf(out, " %d", Aggregate[i].Monomer[j]);
        }
        putc('\n', out);
      }
    }

    fclose(out); //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}
  }

  fclose(vcf);

  if (!silent) {
    if (script) {
      printf("Last Step: %6d\n", count);
    } else {
      fflush(stdout);
      printf("\rLast Step: %6d\n", count);
    }
  } //}}}

  // print last step number to <output.agg> //{{{
  // open output .agg file for appending
  if ((out = fopen(output_agg, "a")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", output_agg);
    exit(1);
  }

  fprintf(out, "\nLast Step: %d\n", count);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(BeadType);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  //}}}

  return 0;
}
