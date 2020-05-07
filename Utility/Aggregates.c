#include "../AnalysisTools.h"
#include "Aggregates.h"

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Aggregates utility determines which molecules belong to which aggregate on \
the basis of given parameters - the minimum distance at which a pair of beads \
from different molecules is considered in contact and the minimum number of \
such contacts between two molecules to consider them as belonging to the same \
aggregate. Only distances between specified bead types are considered. \
Information about aggregates in each timestep is written to '.agg' file (see \
documentation for the format of this file). Coordinates of joined aggregates \
can be written to an output '.vcf' file (with indexed timesteps).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <distance> <contacts> ", cmd);
  fprintf(ptr, "<output.agg> <bead name(s)> <options>\n\n");

  fprintf(ptr, "   <input>               input coordinate file (either vcf or vtf format)\n");
  fprintf(ptr, "   <distance>            minimum distance for beads to be considered in contact\n");
  fprintf(ptr, "   <contacts>            minimum number of contacts for aggregate check\n");
  fprintf(ptr, "   <output.agg>          output filename with '.agg' ending\n");
  fprintf(ptr, "   <bead name(s)>        names of bead types for closeness calculation\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      -x <mol name(s)>   exclude specified molecule(s)\n");
  fprintf(ptr, "      -xm <mol name(s)>  exclude molecules close to specified molecule(s)\n");
  fprintf(ptr, "      -j <output.vcf>    output vcf file with joined coordinates\n");
  CommonHelp(error);
} //}}}

// CalculateAggregates() //{{{
/**
 * Function to determine distribution of molecules in aggregates.
 */
void CalculateAggregates(Aggregate **Aggregate, Counts *Counts, int sqdist, int contacts,
                         int *xm_mols, bool **xm_use_mol,
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule **Molecule) {

  // zeroize //{{{
  (*Counts).Aggregates = 0;
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Aggregate)[i].nMolecules = 0;
    (*Aggregate)[i].nBeads = 0;
    (*Aggregate)[i].nMonomers = 0;
  }

  for (int i = 0; i < (*Counts).Beads; i++) {
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
  int *Link = malloc((*Counts).Beads*sizeof(int));

  // initialize Head array //{{{
  for (int i = 0; i < (n_cells.x*n_cells.y*n_cells.z); i++) {
    Head[i] = -1;
  } //}}}

  // sort beads into cells //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
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

  // disqualify '-xm'ed molecules //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    if (xm_mols[(*Molecule)[i].Type]) {
      (*xm_use_mol)[i] = false;
    }
  }

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
            int c2z = c1z + Dcz[k];

            // periodic boundary conditions for cells
            if (c2x >= n_cells.x)
              c2x -= n_cells.x;
            else if (c2x < 0)
              c2x += n_cells.x;

            if (c2y >= n_cells.y)
              c2y -= n_cells.y;
            else if (c2y < 0)
              c2y += n_cells.y;

            if (c2z >= n_cells.z)
              c2z -= n_cells.z;

            // select second cell
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y; //}}}

            // select bead in the cell 'cell2' //{{{
            int j;
            if (cell1 == cell2) { // next bead in 'cell1'
              j = Link[i];
            } else { // first bead in 'cell2'
              j = Head[cell2];
            } //}}}

            while (j != -1) {
              if ((*Bead)[i].Molecule != -1 && (*Bead)[j].Molecule != -1) { // both i and j must be in molecule)
                int btype_i = (*Bead)[i].Type;
                int btype_j = (*Bead)[j].Type;
                int mol_i = (*Bead)[i].Molecule;
                int mol_j = (*Bead)[j].Molecule;
                int mtype_i = (*Molecule)[mol_i].Type;
                int mtype_j = (*Molecule)[mol_j].Type;

                // one must be used and the other not
                if ((BeadType[btype_i].Use && (*xm_use_mol)[mol_i] && xm_mols[mtype_j]) ||
                    (BeadType[btype_j].Use && (*xm_use_mol)[mol_j] && xm_mols[mtype_i])) {

                  Vector rij = Distance((*Bead)[i].Position, (*Bead)[j].Position, BoxLength);
                  rij.x = SQR(rij.x) + SQR(rij.y) + SQR(rij.z);

                  if (rij.x <= sqdist) {
                    if ((*xm_use_mol)[mol_i]) {
                      (*xm_use_mol)[mol_i] = false;
                    } else {
                      (*xm_use_mol)[mol_j] = false;
                    }
                    break;
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
            int c2z = c1z + Dcz[k];

            // periodic boundary conditions for cells
            if (c2x >= n_cells.x)
              c2x -= n_cells.x;
            else if (c2x < 0)
              c2x += n_cells.x;

            if (c2y >= n_cells.y)
              c2y -= n_cells.y;
            else if (c2y < 0)
              c2y += n_cells.y;

            if (c2z >= n_cells.z)
              c2z -= n_cells.z;

            // select second cell
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y; //}}}

            // select bead in the cell 'cell2' //{{{
            int j;
            if (cell1 == cell2) { // next bead in 'cell1'
              j = Link[i];
            } else { // first bead in 'cell2'
              j = Head[cell2];
            } //}}}

            while (j != -1) {
              int mol_i = (*Bead)[i].Molecule;
              int mol_j = (*Bead)[j].Molecule;
              if (mol_i != -1 && mol_j != -1) { // both i and j must be in molecule
                int btype_i = (*Bead)[i].Type;
                int btype_j = (*Bead)[j].Type;
                int mtype_i = (*Molecule)[mol_i].Type;
                int mtype_j = (*Molecule)[mol_j].Type;

                if (BeadType[btype_i].Use && BeadType[btype_j].Use && // beads must be of specified type
                    MoleculeType[mtype_i].Use && MoleculeType[mtype_j].Use && // molecules can't be excluded via -x option
                    (*xm_use_mol)[mol_i] && (*xm_use_mol)[mol_j]) { // molecules can't be excluded via -xm option

                  // calculate distance between i and j beads
                  Vector rij = Distance((*Bead)[i].Position, (*Bead)[j].Position, BoxLength);
                  rij.x = SQR(rij.x) + SQR(rij.y) + SQR(rij.z);

                  // are 'i' and 'j' close enough?
                  if (mol_i != mol_j && rij.x <= sqdist) {
                    // xm option
                    if (i > j) {
                      contact[mol_i][mol_j]++;
                    } else {
                      contact[mol_j][mol_i]++;
                    }
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
          (*Molecule)[j].Aggregate = agg_j;

          (*Aggregate)[agg_j].nMolecules = 1;
          (*Aggregate)[agg_j].Molecule[0] = j;

          (*Counts).Aggregates++;
        } //}}}

        // add 'mol_i' to 'agg_j' aggregate (that contains 'mol_j' molecule) if 'i' isn't in an agg //{{{
        if (agg_i == -1) {
          int mols = (*Aggregate)[agg_j].nMolecules;
          (*Aggregate)[agg_j].Molecule[mols] = i;
          (*Aggregate)[agg_j].nMolecules++;

          (*Molecule)[i].Aggregate = agg_j;
        } //}}}

        // 'mol_i' and 'mol_j' molecules are in different aggregate => unite aggregates
        if (agg_i != -1 && agg_j != -1 && agg_i != agg_j) {

          // add molecules from aggregate 'agg_i' to 'agg_j' //{{{
          int mols = (*Aggregate)[agg_j].nMolecules;
          (*Aggregate)[agg_j].nMolecules += (*Aggregate)[agg_i].nMolecules;

          // copy molecule ids from Aggregate[agg_i-1] to Aggregate[agg_j-1]
          int id1 = 0;
          for (int k = mols; k < (*Aggregate)[agg_j].nMolecules; k++) {
            int mol = (*Aggregate)[agg_i].Molecule[id1];
            (*Aggregate)[agg_j].Molecule[k] = mol;
            (*Molecule)[mol].Aggregate = agg_j;
            id1++;
          } //}}}

          // move aggregates with id greater then agg_i to id-1 //{{{
          for (int k = (agg_i+1); k < (*Counts).Aggregates; k++) {

            (*Aggregate)[k-1].nMolecules = (*Aggregate)[k].nMolecules;

            // move every molecule from aggregate 'k' to aggregate 'k-1'
            for (int l = 0; l < (*Aggregate)[k].nMolecules; l++) {
              int mol = (*Aggregate)[k].Molecule[l];
              (*Aggregate)[k-1].Molecule[l] = mol;
              (*Molecule)[mol].Aggregate = k - 1;
            }
          } //}}}

          // reduce number of aggregates (two aggregates were merged)
          (*Counts).Aggregates--;
        } //}}}
      } else if (agg_j == -1) { // or 'i' and 'j' aren't in contact and 'j' isn't in any aggregate =>  new aggregate for 'j' */ //{{{
        int agg_j = (*Counts).Aggregates;
        (*Molecule)[j].Aggregate = agg_j;

        (*Aggregate)[agg_j].nMolecules = 1;
        (*Aggregate)[agg_j].Molecule[0] = j;

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
    int aggs = (*Counts).Aggregates;
    (*Aggregate)[aggs].nMolecules = 1;
    (*Aggregate)[aggs].Molecule[0] = (*Counts).Molecules - 1;

    (*Counts).Aggregates++;
  } //}}}
  //}}}

  // sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    SortArray(&(*Aggregate)[i].Molecule, (*Aggregate)[i].nMolecules, 0);
  } //}}}

  // sort aggregates according to ascending ids of first molecules
  SortAggStruct(&(*Aggregate), *Counts);

  // assign bonded beads to Aggregate struct //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    // go through all molecules in aggregate 'i'
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];

      // copy all bead in molecule 'mol' to Aggregate struct
      int mtype = (*Molecule)[mol].Type;
      for (int k = 0; k < MoleculeType[mtype].nBeads; k++) {
        int beads = (*Aggregate)[i].nBeads;
        (*Aggregate)[i].Bead[beads] = (*Molecule)[mol].Bead[k];
        (*Aggregate)[i].nBeads++;

        // every bead from molecule 'mol' is only in one aggregate
        int id = (*Molecule)[mol].Bead[k];
        (*Bead)[id].nAggregates = 1;
        (*Bead)[id].Aggregate[0] = i;
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
                int beads_j = (*Aggregate)[agg_j].nMonomers;

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
                  if ((SQR(rij.x)+SQR(rij.y)+SQR(rij.z)) <= sqdist) {
                    (*Aggregate)[agg_j].Monomer[beads_j] = i;
                    (*Aggregate)[agg_j].nMonomers++;

                    int aggs = (*Bead)[i].nAggregates;
                    (*Bead)[i].nAggregates++;
                    (*Bead)[i].Aggregate = realloc((*Bead)[i].Aggregate, (*Bead)[i].nAggregates*sizeof(int));
                    (*Bead)[i].Aggregate[aggs] = agg_j;
                  }
                } //}}}
              } else if ((*Bead)[j].Molecule == -1 && // monomeric 'j'
                         (*Bead)[i].Molecule != -1) { // 'i' in molecule //{{{

                int agg_i = (*Bead)[i].Aggregate[0];
                int mono_i = (*Aggregate)[agg_i].nMonomers;

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
                  if ((SQR(rij.x)+SQR(rij.y)+SQR(rij.z)) <= sqdist) {
                    (*Aggregate)[agg_i].Monomer[mono_i] = j;
                    (*Aggregate)[agg_i].nMonomers++;

                    int aggs = (*Bead)[j].nAggregates;
                    (*Bead)[j].nAggregates++;
                    (*Bead)[j].Aggregate = realloc((*Bead)[j].Aggregate, (*Bead)[j].nAggregates*sizeof(int));
                    (*Bead)[j].Aggregate[aggs] = agg_i;
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

  // sort monomers in aggregates according to ascending ids //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    SortArray(&(*Aggregate)[i].Monomer, (*Aggregate)[i].nMonomers, 0);
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
  int req_args = 5; //}}}

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
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--script") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-xm") != 0 &&
        strcmp(argv[i], "-j") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - filename of input coordinate file //{{{
  char input_coor[LINE];
  strcpy(input_coor, argv[++count]);

  // test if <input> ends with '.vcf' or '.vtf' (required by VMD)
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_coor, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  }
  // if vtf, copy to input_vsf
  char *input_vsf = calloc(LINE,sizeof(char));
  if (strcmp(strrchr(input_coor, '.'),".vtf") == 0) {
    strcpy(input_vsf, input_coor);
  } else {
    strcpy(input_vsf, "traject.vsf");
  } //}}}

  // <distance> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<distance>");
    Help(argv[0], true);
    exit(1);
  }
  double distance = atof(argv[count]); //}}}

  // <contacts> - number of steps to skip per one used //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    ErrorNaN("<contacts>");
    Help(argv[0], true);
    exit(1);
  }
  int contacts = atoi(argv[count]); //}}}

  // <output.agg> - filename of output agg file (must end with .agg) //{{{
  char output_agg[LINE];
  strcpy(output_agg, argv[++count]);

  // test if <output.agg> ends with '.agg'
  ext = 1;
  strcpy(extension[0], ".agg");
  if (ErrorExtension(output_agg, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  bool script;
  CommonOptions(argc, argv, &input_vsf, &verbose, &silent, &script);

  // save coordinates of joined aggregates //{{{
  char joined_vcf[LINE];
  if (JoinCoorOption(argc, argv, joined_vcf)) {
    exit(1);
  }

  // test if <joined.vcf> filename ends with '.vcf' (required by VMD)
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (joined_vcf[0] != '\0') {
    if (ErrorExtension(joined_vcf, ext, extension)) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    for (int i = 0; i < argc; i++)
      fprintf(stdout, " %s", argv[i]);
    fprintf(stdout, "\n\n");
  } //}}}

  // variables - structures //{{{
  BeadType *BeadType; // structure with info about all bead types
  MoleculeType *MoleculeType; // structure with info about all molecule types
  Bead *Bead; // structure with info about every bead
  int *Index; // reverse of Bead[].Index
  Molecule *Molecule; // structure with info about every molecule
  Counts Counts; // structure with number of beads, molecules, etc. //}}}

  // read system information
  bool indexed = ReadStructure(input_vsf, input_coor, &Counts, &BeadType, &Bead, &Index, &MoleculeType, &Molecule);

  // vsf file is not needed anymore
  free(input_vsf);

  // <type names> - names of bead types to use for closeness calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      ErrorBeadType(input_coor, argv[count], Counts, BeadType);
      exit(1);
    }
    if (BeadType[type].Use) {
      fprintf(stderr, "\nError: bead type %s specified more than once\n\n", argv[count]);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }

  // used molecule type = write molecule type -- for now
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  } //}}}

  // '-xm' option //{{{
  int *xm_mols = calloc(Counts.TypesOfMolecules,sizeof(int *));
  if (MoleculeTypeOption2(argc, argv, "-xm", &xm_mols, Counts, &MoleculeType)) {
    exit(1);
  }

  // set all individual molecules to be used - changes in CalculateAggregates()
  bool *xm_use_mol = calloc(Counts.Molecules, sizeof(bool));

  // is -xm in use?
  bool xm = false;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (xm_mols[i]) {
      xm = true;
      break;
    }
  } //}}}

  // print command to output .agg file //{{{
  FILE *out;
  if ((out = fopen(output_agg, "w")) == NULL) {
    ErrorFileOpen(output_agg, 'w');
    exit(1);
  }

  // print command to stdout
  for (int i = 0; i < argc; i++)
    fprintf(out, " %s", argv[i]);
  putc('\n', out);

  fclose(out); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  } //}}}

  Vector BoxLength = GetPBC(vcf, input_coor);

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
      ErrorFileOpen(joined_vcf, 'w');
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

  // create array for the first line of a timestep ('# <number and/or other comment>')
  char *stuff = calloc(LINE, sizeof(char));

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
    VerboseOutput(input_coor, Counts, BoxLength, BeadType, Bead, MoleculeType, Molecule);

    fprintf(stdout, "\n   Distance for closeness check: %lf\n", distance);
    fprintf(stdout, "   Number of needed contacts for aggregate check: %d\n", contacts);
    if (xm) {
      fprintf(stdout, "   Ignore molecules close to:");
      for (int i = 0; i < Counts.TypesOfMolecules; i++) {
        if (xm_mols[i]){
          fprintf(stdout, " %s", MoleculeType[i].Name);
        }
        putchar('\n');
      }
    }
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  int test;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    count++;

    // print step? //{{{
    if (!silent && !script) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}

    for (int i = 0; i < Counts.Molecules; i++) {
      xm_use_mol[i] = true;
    }

    // read coordinates //{{{
    if ((test = ReadCoordinates(indexed, vcf, Counts, Index, &Bead, &stuff)) != 0) {
      // print newline to stdout if Step... doesn't end with one
      ErrorCoorRead(input_coor, test, count, stuff, input_vsf);
      exit(1);
    } //}}}

    RestorePBC(Counts, BoxLength, &Bead);

    CalculateAggregates(&Aggregate, &Counts, SQR(distance), contacts, xm_mols, &xm_use_mol, BoxLength, BeadType, &Bead, MoleculeType, &Molecule);

    // calculate & write joined coordinatest to <joined.vcf> if '-j' option is used //{{{
    if (joined_vcf[0] != '\0') {
//    PrintAggregate(Counts, Index, MoleculeType, Molecule, Bead, BeadType, Aggregate);
      RemovePBCMolecules(Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

      // open <joined.vcf> file //{{{
      FILE *joined;
      if ((joined = fopen(joined_vcf, "a")) == NULL) {
        ErrorFileOpen(joined_vcf, 'a');
        exit(1);
      } //}}}

      WriteCoorIndexed(joined, Counts, BeadType, Bead, MoleculeType, Molecule, stuff);

      fclose(joined);
    } //}}}

    // open output .agg file for appending //{{{
    if ((out = fopen(output_agg, "a")) == NULL) {
      ErrorFileOpen(output_agg, 'a');
      exit(1);
    } //}}}

    // write data to output .agg file //{{{
    // find the number of aggregates - remove aggregates only of excluded mols //{{{
    int no_excluded_aggs = 0;
    int test_count = 0; // to test that every molecule is in an aggregate
    for (int i = 0; i < Counts.Aggregates; i++) {
      Aggregate[i].Use = false;

      test_count += Aggregate[i].nMolecules;

      if (Aggregate[i].nMolecules != 1 || xm_use_mol[Aggregate[i].Molecule[0]]) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
          if (MoleculeType[moltype].Use) {
            Aggregate[i].Use = true;
            no_excluded_aggs++;
            break;
          }
        }
      }
    } //}}}

    // are all molecules accounted for? //{{{
    if (test_count != Counts.Molecules) {
      fprintf(stderr, "\nError: not all molecules were assigned to aggregates\n");
      fprintf(stderr, "       Counts.Molecules = %5d; Molecules in aggregates: %d\n\n", Counts.Molecules, test_count);
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
          fprintf(out, " %d", Bead[Aggregate[i].Monomer[j]].Index);
        }
        putc('\n', out);
      }
    }

    fclose(out); //}}}
  }

  fclose(vcf);

  if (!silent) {
    if (script) {
      fprintf(stdout, "Last Step: %d\n", count);
    } else {
      fflush(stdout);
      fprintf(stdout, "\r                       ");
      fprintf(stdout, "\rLast Step: %d\n", count);
    }
  } //}}}

  // print last step number to <output.agg> //{{{
  // open output .agg file for appending
  if ((out = fopen(output_agg, "a")) == NULL) {
    ErrorFileOpen(output_agg, 'a');
    exit(1);
  }

  fprintf(out, "\nLast Step: %d\n", count);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(xm_use_mol);
  free(xm_mols);
  free(BeadType);
  free(Index);
  FreeAggregate(Counts, &Aggregate);
  FreeMoleculeType(Counts, &MoleculeType);
  FreeMolecule(Counts, &Molecule);
  FreeBead(Counts, &Bead);
  free(stuff);
  //}}}

  return 0;
}
