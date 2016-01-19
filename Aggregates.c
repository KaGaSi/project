#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "AnalysisTools.h"

void ErrorHelp(char cmd[50]) { //{{{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "   %s <input.vcf> <distance> <contacts> ", cmd);
  fprintf(stderr, "<output.agg> <type names> <options>\n\n");

  fprintf(stderr, "   <input.vcf>         input filename (vcf format)\n");
  fprintf(stderr, "   <distance>          minimum distance for contact for aggregate check\n");
  fprintf(stderr, "   <contacts>          minimum number of contacts for aggregate check\n");
  fprintf(stderr, "   <output.agg>        output filename (agg format)\n");
  fprintf(stderr, "   <type names>        names of bead types to use for closeness calculation\n");
  fprintf(stderr, "   <options>\n");
  fprintf(stderr, "      -i <name>        use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>        file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -j <joined.vcf>  output vcf file with joined coordinates\n");
  fprintf(stderr, "      -v               verbose output\n");
  fprintf(stderr, "      -V               verbose output with comments from input .vcf file\n");
  fprintf(stderr, "      -h               print this help and exit\n");
} //}}}

void CalculateAggregates(Aggregate **Aggregate, Counts *Counts, int sqdist, int contacts, //{{{
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule) {

  // zeroize //{{{
  (*Counts).Aggregates = 0;
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Aggregate)[i].nMolecules = 0;
    (*Aggregate)[i].nBeads = 0;
    (*Aggregate)[i].nMonomers = 0;
  } //}}}

  // allocate & zeroize contact[][] (triangular matrix) and moved array //{{{
  int **contact = malloc((*Counts).Molecules*sizeof(int *));
  int *moved = malloc((*Counts).Molecules*sizeof(int));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    contact[i] = malloc((i+1)*sizeof(int));
  }

  // zeroize arrays
  for (int i = 0; i < (*Counts).Molecules; i++) {
    moved[i] = 0;
    for (int j = 0; j < i; j++) /* j == i is really not needed */
      contact[i][j] = 0;
  } //}}}

  // count contacts between all molecules pairs //{{{
  // go over all pairs of molecules
  for (int i = 1; i < (*Counts).Molecules; i++) { // first molecule
    for (int j = 0; j < i; j++) { // second molecule

      // go over all beads in first molecule
      for (int k = 0; k < MoleculeType[Molecule[i].Type].nBeads; k++) {
        int id1 = Molecule[i].Bead[k];

        // go over all beads in second molecule if the bead type of id1 is in use
        for (int l = 0; BeadType[(*Bead)[id1].Type].Use && l < MoleculeType[Molecule[j].Type].nBeads; l++) {
          int id2 = Molecule[j].Bead[l];

          // should bead of this type be used to calculate aggregates?
          if (BeadType[(*Bead)[id2].Type].Use) {

            // calculate distance between k-th bead in molecule i and l-th bead in molecule j
            Vector rij = DistanceBetweenBeads(id1, id2, *Bead, BoxLength);

            // are 'id1' and 'id2' close enough?
            if ((SQR(rij.x) + SQR(rij.y) + SQR(rij.z)) <= sqdist)
              contact[i][j]++;
          }
        }
      }
    }
  } //}}}

  // evaluate the contacts //{{{
  // first molecule
  for (int i = 1; i < (*Counts).Molecules; i++) {

    // second molecule
    for (int j = 0; j < i; j++) {

      // find out if Molecule[i] is in any aggregate //{{{
      int testi = -1;
      for (int k = 0; k < (*Counts).Aggregates; k++) {
        for (int l = 0; l < (*Aggregate)[k].nMolecules; l++) {
          if (i == (*Aggregate)[k].Molecule[l]) {
            testi = k;
            break;
          }
        }

        // if 'i' is in aggregate, no need to go through the rest of aggregates
        if (testi != -1)
          break;
      } //}}}

      // find out if Molecule[j] is in any aggregate //{{{
      int testj = -1;
      for (int k = 0; k < (*Counts).Aggregates; k++) {
        for (int l = 0; l < (*Aggregate)[k].nMolecules; l++) {
          if (j == (*Aggregate)[k].Molecule[l]) {
            testj = k;
            break;
          }
        }

        // if 'j' is in aggregate, no need to go through the rest of aggregates
        if (testj != -1)
          break;
      } //}}}

 //   printf("%d\n", (*Counts).Aggregates);
 //         printf("testi=%d testj=%d\n", testi, testj);

      // molecules 'i' and 'j' are in contact //{{{
      if (contact[i][j] >= contacts) {
        // create new aggregate if 'j' isn'it in any //{{{
        if (testj == -1) {
          testj = (*Counts).Aggregates;

          (*Aggregate)[testj].nMolecules = 1;
          (*Aggregate)[testj].Molecule[0] = j;

          (*Counts).Aggregates++;
        } //}}}

        // add 'i' to aggregate if 'i' isn't in any //{{{
        if (testi == -1) {
          (*Aggregate)[testj].Molecule[(*Aggregate)[testj].nMolecules] = i;

          (*Aggregate)[testj].nMolecules++;
        } //}}}

        // each residue in different aggregate => unite aggregates
        if (testi != -1 && testj != -1 && testi != testj) {

          // add molecules from aggregate 'testi' to 'testj' //{{{
          int mols = (*Aggregate)[testj].nMolecules;

          (*Aggregate)[testj].nMolecules += (*Aggregate)[testi].nMolecules;
          int id1 = 0;
          // copy molecule ids from Aggregate[testi-1] to Aggregate[testj-1] */
          for (int k = mols; k < (*Aggregate)[testj].nMolecules; k++) {
            (*Aggregate)[testj].Molecule[k] = (*Aggregate)[testi].Molecule[id1];
            id1++;
          } //}}}

          // move aggregates with id greater then testi to id-1 //{{{
          for (int k = (testi+1); k < (*Counts).Aggregates; k++) {

            (*Aggregate)[k-1].nMolecules = (*Aggregate)[k].nMolecules;

            // move every molecule from aggregate 'k' to aggregate 'k-1'
            for (int l = 0; l < (*Aggregate)[k].nMolecules; l++) {
              (*Aggregate)[k-1].Molecule[l] = (*Aggregate)[k].Molecule[l];
            }
          } //}}}

          // reduce number of aggregates (two aggregates were merged)
          (*Counts).Aggregates--;
        } //}}}
      // or 'i' and 'j' aren't in contact and 'j' isn't in any aggregate =>  new aggregate for 'j' */ //{{{
      } else if (testj == -1) {
        (*Aggregate)[(*Counts).Aggregates].nMolecules = 1;
        (*Aggregate)[(*Counts).Aggregates].Molecule[0] = j;

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

  FillAggregateBeads(Aggregate, *Counts, MoleculeType, Molecule);

  // calculate the number of monomeric beads in aggregate //{{{
  // go through all unbonded beads
  for (int i = 0; i < ((*Counts).Unbonded+(*Counts).Bonded); i++) {
    if ((*Bead)[i].Molecule == -1) { // -1 means 'in no molecule'

      // go through all aggregates
      for (int j = 0; j < (*Counts).Aggregates; j++) {
        bool in_agg = false;

        // go through all molecules in aggregate 'j'
        for (int k = 0; k < (*Aggregate)[j].nMolecules; k++) {
          int id = (*Aggregate)[j].Molecule[k];

          // go through all beads in molecule 'id'
          for (int l = 0; l < MoleculeType[Molecule[id].Type].nBeads; l++) {

            // calculate distance between monomeric bead 'i' and bead 'Molecule[id].Bead[l]'
            Vector dist = DistanceBetweenBeads(i, Molecule[id].Bead[l], *Bead, BoxLength);

            if ((SQR(dist.x)+SQR(dist.y)+SQR(dist.z)) < sqdist) {
              (*Aggregate)[j].Monomer[(*Aggregate)[j].nMonomers] = (*Bead)[i].Index;
              (*Aggregate)[j].nMonomers++;

              in_agg = true;

              break;
            }
          }

          if (in_agg)
            break;
        }
      }
    }
  } //}}}

//// remove aggregate PBC
//bool *moved = malloc(Counts.Molecules*sizeof(bool));

//// go through all aggregates larger than unimers
//for (int i = 0; i < Counts.Aggregates; i++) {

//  // skip removing pbc for unimers
//  if (Aggregate[i].nMolecules == 1)
//    continue;

//  // negate moved array, while first bead is not to move //{{{
//  for (int j = 1; j < Counts.Molecules; j++) {
//    moved[j] = false;
//  }
//  moved[0] = true; //}}}

//  bool done = false;
//  while (!done) {

//    // go through all molecule pairs
//    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
//      for (int k = 0; moved[j] && k < Aggregate[i].nMolecules; k++) {

//        // use only moved molecule 'mol_j' and unmoved molecule 'mol_k'
//        if (moved[j] && !moved[k]) { // automatically follows that j!=k
//          int mol1 = Aggregate[i].Molecule[j];
//          int mol2 = Aggregate[i].Molecule[k];

//          // go through all bead pairs in the two molecules
//          for (int l = 0; l < MoleculeType[Molecule[mol1].Type].nBeads; l++) {
//            for (int m = 0; m < MoleculeType[Molecule[mol2].Type].nBeads; m++) {
//            }
//            // if molekule 'k' (or 'mol2') has been moved, skip also remainder of 'mol1'
//            if (moved[k]) {
//              break;
//            }
//          }
//        }
//      }
//    }

//    done = true;
//    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
//      if (!moved[j]) {
//        done = false;
//        break;
//      }
//    }
//  }
//}

  // free memory //{{{
  for (int i = 0; i < (*Counts).Molecules; i++)
    free(contact[i]);
  free(contact);
  free(moved); //}}}
} //}}}

void RemovePBCAggregates(double distance, Aggregate *Aggregate, Counts Counts, //{{{
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule) {

  bool *moved = malloc(Counts.Molecules*sizeof(bool));

  // go through all aggregates larger than unimers
  for (int i = 0; i < Counts.Aggregates; i++) {

    // negate moved array, while first bead is not to move //{{{
    for (int j = 1; j < Counts.Molecules; j++) {
      moved[j] = false;
    }
    moved[0] = true; //}}}

    bool done = false;
    while (!done) {

      // go through all molecule pairs
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        for (int k = 0; k < Aggregate[i].nMolecules; k++) {

          // use only moved molecule 'mol1' and unmoved molecule 'mol2'
          if (moved[j] && !moved[k]) { // automatically follows that j != k
            int mol1 = Aggregate[i].Molecule[j];
            int mol2 = Aggregate[i].Molecule[k];

            // go through all bead pairs in the two molecules
            for (int l = 0; l < MoleculeType[Molecule[mol1].Type].nBeads; l++) {
              for (int m = 0; m < MoleculeType[Molecule[mol2].Type].nBeads; m++) {
                int bead1 = Molecule[mol1].Bead[l];
                int bead2 = Molecule[mol2].Bead[m];

                // use only bead types that were used to assign molecules to aggregates
                if (BeadType[(*Bead)[bead1].Type].Use &&
                    BeadType[(*Bead)[bead2].Type].Use) {

                  // calculate distance between 'bead1' and 'bead2'
                  Vector dist = DistanceBetweenBeads(bead1, bead2, *Bead, BoxLength);
                  dist.x = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));

                  // move 'mol2' if 'bead1' and 'bead2' are in contact
                  if (dist.x < distance) {

                    // distance vector between 'bead1' and 'bead2'
                    dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;

                    // if 'bead1' and 'bead2' are too far in x-direction, move 'mol2' in x-direction //{{{
                    if (dist.x > (BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x += BoxLength.x;
                      }
                    } else if (dist.x <= -(BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x -= BoxLength.x;
                      }
                    } //}}}

                    // if 'bead1' and 'bead2' are too far in y-direction, move 'mol2' in y-direction //{{{
                    if (dist.y > (BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y += BoxLength.y;
                      }
                    } else if (dist.y <= -(BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y -= BoxLength.y;
                      }
                    } //}}}

                    // if 'bead1' and 'bead2' are too far in z-direction, move 'mol2' in x-direction //{{{
                    if (dist.z > (BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z += BoxLength.z;
                      }
                    } else if (dist.z <= -(BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[Molecule[mol2].Type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z -= BoxLength.z;
                      }
                    } //}}}

                    moved[k] = true;

                    // skip remainder of 'mol2' (or 'k')
                    break;
                  }
                }
              }
              // if molekule 'k' (or 'mol2') has been moved, skip also remainder of 'mol1'
              if (moved[k]) {
                break;
              }
            }
          }
        }
      }

      done = true;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (!moved[j]) {
          done = false;
          break;
        }
      }
    }
  }

  free(moved);
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
      printf(".agg file. Also joined coordinates can be written to an output .vcf file.   \n");
      printf("The utility uses dl_meso.vsf (or other input structure file) and FIELD      \n");
      printf("(along with optional bond file) files to determine all information about    \n");
      printf("the system.                                                                 \n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <distance> <contacts> ", argv[0]);
      printf("<output.agg> <type names> <options>\n\n");

      printf("   <input.vcf>         input filename (vcf format)\n");
      printf("   <output.agg>        output filename (agg format)\n");
      printf("   <distance>          minimum distance for contact for aggregate check\n");
      printf("   <contacts>          minimum number of contacts for aggregate check\n");
      printf("   <type names>        names of bead types for closeness calculation\n");
      printf("   <options>\n");
      printf("      -i <name>        use input .vsf file different from dl_meso.vsf\n");
      printf("      -b <name>        file containing bond alternatives to FIELD\n");
      printf("      -j <joined.vcf>  output vcf file with joined coordinates\n");
      printf("      -v               verbose output\n");
      printf("      -V               verbose output with comments from input .vcf file\n");
      printf("      -h               print this help and exit\n");
      exit(0);
    }
  } //}}}

  // print command to stdout //{{{
  for (int i = 0; i < argc; i++)
    printf(" %s", argv[i]);
  printf("\n\n"); //}}}

  // check if correct number of arguments //{{{
  if (argc < 6) {
    fprintf(stderr, "Too little arguments!\n\n");
    ErrorHelp(argv[0]);
    exit(1);
  } //}}}

  // -i <name> option - filename of input structure file //{{{
  char vsf_file[32];
  vsf_file[0] = '\0'; // check if -i option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-i' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");
        exit(1);
      }

      // check if .vsf ending is present
      char *vsf = strrchr(argv[i+1], '.');
      if (!vsf || strcmp(vsf, ".vsf")) {
        fprintf(stderr, "'-i' arguments does not have .vsf ending!\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(vsf_file, argv[i+1]);
    }
  }

  // -i option is not used
  if (vsf_file[0] == '\0') {
    strcpy(vsf_file, "dl_meso.vsf");
  } //}}}

  // -b <name> option - filename of input bond file //{{{
  char bonds_file[32];
  bonds_file[0] = '\0'; // check if -b option is used
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-b") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-b' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n\n");
        ErrorHelp(argv[0]);
        exit(1);
      }

      strcpy(bonds_file, argv[i+1]);
    }
  } //}}}

  // -j <joined.vcf> - filename of output vcf file (must end with .vcf) //{{{
  char joined_vcf[32];
  joined_vcf[0] = '\0'; // no -j option
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-j") == 0) {

      // wrong argument to -i option
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        fprintf(stderr, "\nMissing argument to '-j' option ");
        fprintf(stderr, "(or filename beginning with a dash)!\n");
        exit(1);
      }

      strcpy(joined_vcf, argv[i+1]);

      // test if <joined.vcf> filename ends with '.vcf' (required by VMD)
      char *dot = strrchr(joined_vcf, '.');
      if (!dot || strcmp(dot, ".vcf")) {
        fprintf(stderr, "<joined.vcf> '%s' does not have .vcf ending!\n", joined_vcf);
        ErrorHelp(argv[0]);
        exit(1);
      }
    }
  } //}}}

  // -v option - verbose output //{{{
  bool verbose = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0) {
      verbose = true;

      break;
    }
  } //}}}

  // -V option - verbose output with comments from input .vcf file //{{{
  bool verbose2 = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-V") == 0) {
      verbose = true;
      verbose2 = true;

      break;
    }
  } //}}}

  int count = 0; // count mandatory arguments

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

  // <type names> - names of bead types to use for closeness calculation //{{{
  while (++count < argc && argv[count][0] != '-') {
    int type = FindType(argv[count], Counts, BeadType);

    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      fprintf(stderr, "Bead type '%s' is not in %s coordinate file!\n", argv[count], input_vcf);
      exit(1);
    }

    BeadType[type].Use = true;
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    printf("\n   Read from FIELD\n");
    printf("Counts.{");
    printf("TypesOfBeads =%3d, ", Counts.TypesOfBeads);
    printf("Bonded =%7d, ", Counts.Bonded);
    printf("Unboded =%7d, ", Counts.Unbonded);
    printf("TypesOfMolecules =%3d, ", Counts.TypesOfMolecules);
    printf("Molecules =%4d}\n", Counts.Molecules);
    printf("\ntotal number of beads: %d\n\n", Counts.Bonded+Counts.Unbonded);

    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      printf("BeadType[%2d].{", i);
      printf("Name =%10s, ", BeadType[i].Name);
      printf("Number =%7d, ", BeadType[i].Number);
      printf("Charge =%6.2f, ", BeadType[i].Charge);
      printf("Mass =%5.2f, ", BeadType[i].Mass);
      printf("Use = %s}\n", BeadType[i].Use ? "True" : "False");
    }
    putchar('\n');

    for (int i = 0; i < Counts.TypesOfMolecules; i++) {
      printf("MoleculeType[%d].{", i);
      printf("Name =%10s", MoleculeType[i].Name);
      printf(", Number =%4d", MoleculeType[i].Number);
      printf(", nBeads =%3d", MoleculeType[i].nBeads);
      printf(", nBonds =%3d", MoleculeType[i].nBonds);
      if (bonds_file[0] == '\0') { // all bonds taken from FIELD
        printf(", Bonds from 'FIELD'}\n");
      } else {
        // go through bond file to find out if molecule type 'i' is there
        FILE *bond;
        if ((bond = fopen(bonds_file, "r")) == NULL) {
          fprintf(stderr, "Cannot open file %s with '-v' option!\n", bonds_file);
          exit(1);
        }

        int test;
        char str[32];
        while ((test = getc(bond)) != EOF) {
          ungetc(test, bond);

          if ((fscanf(bond, "%s %d", str, &test)) != 2) {
            fprintf(stderr, "Cannot read string or number of bonds from %s with '-v' option!\n", bonds_file);
            exit(1);
          }

          if (strcmp(str, MoleculeType[i].Name) == 0) {
            printf(", Bonds from '%s'}\n", bonds_file);
            break;
          }

          while (getc(bond) != '\n')
            ;
        }

        // if not in bonds_file, then bonds taken from FIELD
        if (test == EOF) {
          printf(", Bonds from 'FIELD'}\n");
        }

        fclose(bond);
      }
    }

    printf("\n   Distance for closeness check: %lf\n", distance);
    printf("   Number of needed contacts for aggregate check: %d\n", contacts);
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

    // for now all bead types are to be written in joined.vcf
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

  // main loop //{{{
  count = 0; // count timesteps
  int test;
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    fflush(stdout);
    printf("\rStep: %6d", ++count);

    // read indexed timestep from input .vcf file //{{{
    if (indexed) {
      if ((test = ReadCoorIndexed(vcf, Counts, &Bead, &stuff)) != 0) {
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      } //}}}
    // or read ordered timestep from input .vcf file //{{{
    } else {
      if ((test = ReadCoorOrdered(vcf, Counts, &Bead, &stuff)) != 0) {
        fprintf(stderr, "Cannot read coordinates from %s! (%d. step; %d. bead)\n", input_vcf, count, test);
        exit(1);
      }
    } //}}}

    CalculateAggregates(&Aggregate, &Counts, SQR(distance), contacts, BoxLength, BeadType, &Bead, MoleculeType, Molecule);

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

      WriteCoorIndexed(joined, Counts, BeadType, Bead, stuff);

      fclose(joined);
    } //}}}

    // open output .agg file for appending //{{{
    if ((out = fopen(output_agg, "a")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", output_agg);
      exit(1);
    } //}}}

    // write data to output .agg file //{{{
    fprintf(out, "\nStep: %d\n%d\n\n", count, Counts.Aggregates);
    int test_count = 0; // to test that all molecules are in aggregate

    // go through all aggregates
    for (int i = 0; i < Counts.Aggregates; i++) {

      test_count += Aggregate[i].nMolecules;

      // go through all molecules in aggregate 'i'
      fprintf(out, "%d :", Aggregate[i].nMolecules);
      for (int j = 0; j < Aggregate[i].nMolecules; j++ ) {
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

    fclose(out); //}}}

    // making sure all molecules are in aggregates //{{{
    if (test_count != Counts.Molecules) {
      fprintf(stderr, "Not all molecules were assigned to aggregates!\n");
      fprintf(stderr, "Counts.Molecules = %5d; Molecules in aggregates: %d\n", Counts.Molecules, test_count);
      exit(1);
    } //}}}

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}
  }

  fclose(vcf);

  fflush(stdout);
  printf("\rLast Step: %6d\n", count); //}}}

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
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      free(MoleculeType[i].Bond[j]);
    }
    free(MoleculeType[i].Bond);
  }
  free(MoleculeType);
  free(Bead);
  for (int i = 0; i < Counts.Molecules; i++) {
    free(Molecule[i].Bead);

    free(Aggregate[i].Molecule);
    free(Aggregate[i].Bead);
    free(Aggregate[i].Monomer);
  }
  free(Molecule);
  free(Aggregate);
  free(stuff);
  //}}}

  return 0;
}
