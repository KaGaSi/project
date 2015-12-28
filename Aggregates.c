#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "AnalysisTools.h"

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
  fprintf(stderr, "      -i <name>      use input .vsf file different from dl_meso.vsf\n");
  fprintf(stderr, "      -b <name>      file containing bond alternatives to FIELD\n");
  fprintf(stderr, "      -v             verbose output\n");
  fprintf(stderr, "      -V             verbose output with comments from input .vcf file\n");
  fprintf(stderr, "      -h             print this help and exit\n");
} //}}}

// CalculateAggregates() //{{{
void CalculateAggregates(Aggregate **Aggregate, Counts *Counts, Vector BoxLength, int sqdist,
                         BeadType *BeadType, Bead *Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule) {

  (*Counts).Aggregates = 0;

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
  for (int i = 0; i < (*Counts).Molecules; i++) { // first molecule
    for (int j = 0; j < i; j++) { // second molecule

      // go over all beads in first molecule
      for (int k = 0; k < MoleculeType[Molecule[i].Type].nBeads; k++) {
        int id1 = Molecule[i].Bead[k];

        // should bead of this type be used to calculate aggregates?
        if (BeadType[Bead[id1].Type].Use == 2) {

          // go over all beads in second molecule
          for (int l = 0; l < MoleculeType[Molecule[j].Type].nBeads; l++) {
            int id2 = Molecule[j].Bead[l];

            // should bead of this type be used to calculate aggregates?
            if (BeadType[Bead[id2].Type].Use == 2) {

              // calculate distance between k-th bead in molecule i and l-th bead in molecule j
              Vector rij = DistanceBetweenBeads(id1, id2, Bead, BoxLength);

              // distance squared
              rij.x = SQR(rij.x) + SQR(rij.y) + SQR(rij.z);

              if (rij.x <= sqdist)
                contact[i][j]++;
            }
          }
        }
      }
    }
  } //}}}

///* calculate number of residues in vcf file by testing
// * if bead type of the first bead in residue is in vcf file */ //{{{
//residues = 0;
//for (i = 0; i < (*Counts).Molecules; i++) { /* starts from 0 to correctly calculate
//                                            the number of residues in vcf file */
//  /* test if residue bead types are in vcf */
//  if (ResidInVCF(i, beads)) { /* use only residues with beads in vcf file */
//    residues++;
//  }
//} //}}}

///* evaluate the contacts */ //{{{
//NumberOfAggregates = 0;
//for (i = 1; i < (*Counts).Molecules; i++) {
//  /* test if residue bead types are in vcf */
//  if (ResidInVCF(i,beads)) { /* use only residues with beads in vcf file */
//    res = i;
//    for (j = 0; j < i; j++) { /* second residue */
//      /* test if residue bead types are in vcf */
//      if (ResidInVCF(i,beads)) { /* use only residues with beads in vcf file */
//        /* find out if Molecule[j] is in any aggregate */ //{{{
//        testj = -1;
//        for (k = 0; k < NumberOfAggregates; k++) {
//          for (l = 1; l <= Aggregate[k].Residues[0]; l++) {
//            if (j == Aggregate[k].Residues[l]) {
//              testj = k + 1;
//              break;
//            }
//          }
//          if (testj != -1) /* j is in aggregate? */
//            break; /* no need to go through the rest of aggregates */
//        } //}}}

//        /* find out if Molecule[i] is in any aggregate */ //{{{
//        testi = -1;
//        for (k = 0; k < NumberOfAggregates; k++) {
//          for (l = 1; l <= Aggregate[k].Residues[0]; l++) {
//            if (i == Aggregate[k].Residues[l]) {
//              testi = k + 1;
//              break;
//            }
//          }
//          if (testi != -1) /* i is in aggregate? */
//            break; /* no need to go through the rest of aggregates */
//        } //}}}

//        /* test whether Molecule[i] and Molecule[j] are in contact */
//        if (contact[i][j] >= contacts) { /* i and j are in contact */
//          /* create new aggregate if j isn'it in any */ //{{{
//          if (testj == -1) {
//            NumberOfAggregates++;
//            testj = NumberOfAggregates;
//            Aggregate[testj-1].Residues[0] = 1;
//            Aggregate[testj-1].Residues[1] = j;

//            Molecule[j].Aggregate = testj - 1;

//            Aggregate[testj-1].Beads[0] = Molecule[j].Beads[0];
//            for (k = 1; k <= Molecule[j].Beads[0]; k++)
//              Aggregate[testj-1].Beads[k] = Molecule[j].Beads[k];
//          } //}}}

//          /* add i to j's aggregate */ //{{{
//          if (testi == -1) {
//            /* one more residue in j's aggregate */ //{{{
//            Aggregate[testj-1].Residues[0]++;
//            /* add i to j's aggregate */
//            Aggregate[testj-1].Residues[Aggregate[testj-1].Residues[0]] = i;
//            /* is in aggregate #(testj-1) */
//            Molecule[i].Aggregate = testj - 1; //}}}

//            /* add i's beads to j's aggregate */ //{{{
//            for (k = 1; k <= Molecule[i].Beads[0]; k++) {
//              Aggregate[testj-1].Beads[Aggregate[testj-1].Beads[0]+k] = Molecule[i].Beads[k];
//            }
//            /* increment number of beads in j's aggregate */
//            Aggregate[testj-1].Beads[0] += Molecule[i].Beads[0]; //}}}
//          } //}}}

//          /* each residue in different aggregate => unite aggregates */ //{{{
//          if (testi != -1 && testj != -1 && testi != testj) {
//            /* make sure testj<testi and i & j correspond to testi & testj */ //{{{
//            if (testj > testi) {
//              id2 = testj;
//              testj = testi;
//              testi = id2;

//              id2 = i;
//              i = j;
//              j = id2;
//            } //}}}

//            /* add beads from aggregate (testi-1) to (testj-1) */ //{{{
//            l = Aggregate[testj-1].Beads[0]; /* initial number of beads in Aggregate[testj-1] */
//            l++;
//            Aggregate[testj-1].Beads[0] += Aggregate[testi-1].Beads[0]; /* new number of beads in Aggregate[testj-1] */
//            /* copy bead ids from Aggregate[testi-1] to Aggregate[testj-1] */
//            id1 = 1;
//            for (k = l; k <= Aggregate[testj-1].Beads[0]; k++) {
//              Aggregate[testj-1].Beads[k] = Aggregate[testi-1].Beads[id1];
//              id1++;
//            } //}}}

//            /* add residues from aggregate (testi-1) to (testj-1) */ //{{{
//            l = Aggregate[testj-1].Residues[0]; /* initial number of residues in Aggregate[testj-1] */
//            l++;
//            Aggregate[testj-1].Residues[0] += Aggregate[testi-1].Residues[0]; /* new number of residues in Aggregate[testj-1] */
//            id1 = 1;
//            /* copy residue ids from Aggregate[testi-1] to Aggregate[testj-1] */
//            for (k = l; k <= Aggregate[testj-1].Residues[0]; k++) {
//              Aggregate[testj-1].Residues[k] = Aggregate[testi-1].Residues[id1];
//              Molecule[Aggregate[testi-1].Residues[id1]].Aggregate = testj - 1;
//              id1++;
//            } //}}}

//            /* move aggregates with id greater then testi-i to id-1 */ //{{{
//            for (k = testi; k < NumberOfAggregates; k++) {
//              /* move k's residues to (k-1) */ //{{{

//              Aggregate[k-1].Residues[0] = Aggregate[k].Residues[0];
//              for (l = 1; l <= Aggregate[k].Residues[0]; l++) {
//                Aggregate[k-1].Residues[l] = Aggregate[k].Residues[l];
//                Molecule[Aggregate[k].Residues[l]].Aggregate = k - 1;
//              } //}}}

//              /* move k's beads to (k-1) */ //{{{
//              for (l = 0; l < Aggregate[k].Beads[0]; l++) {
//                Aggregate[k-1].Beads[l] = Aggregate[k].Beads[l];
//              } //}}}
//            } //}}}

//            /* reduce number of aggregates (two aggregates were merged) */
//            NumberOfAggregates--;
//          } //}}}

//          /* residues i and j aren't in contact and j isn't in any aggregate,
//           * so create a new aggregate with just one residue */ //{{{
//        } else if (testj == -1) {
//          NumberOfAggregates++;
//          Aggregate[NumberOfAggregates-1].Residues[0] = 1;
//          Aggregate[NumberOfAggregates-1].Residues[1] = j;

//          Molecule[j].Aggregate = NumberOfAggregates - 1;

//          Aggregate[NumberOfAggregates-1].Beads[0] = Molecule[j].Beads[0];
//          for (k = 1; k <= Molecule[j].Beads[0]; k++)
//            Aggregate[NumberOfAggregates-1].Beads[k] = Molecule[j].Beads[k];
//        } //}}}
//      }
//    }
//  }
//} //}}}

///* if residue with highest id is in no aggregate,
// * a separate aggregate must to be created */ //{{{

///* check if highest id residue is in aggregate */ //{{{
//testi = -1;
//l = 0; /* total number in all aggregates for later check */
//for (i = 0; i < NumberOfAggregates; i++) {
//  l += Aggregate[i].Residues[0];
//  for (j = 1; j <= Aggregate[i].Residues[0]; j++) {
//    if (Aggregate[i].Residues[j] == res) {
//      testi = 1;
//    }
//  }
//} //}}}

///* highest id residue not in any aggregate => create separate one */ //{{{
//if (testi == -1) {
//  l++;
//  NumberOfAggregates++;
//  Aggregate[NumberOfAggregates-1].Residues[0] = 1;
//  Aggregate[NumberOfAggregates-1].Residues[1] = res;

//  Molecule[(*Counts).Molecules-1].Aggregate = NumberOfAggregates - 1;

//  Aggregate[NumberOfAggregates-1].Beads[0] = Molecule[res].Beads[0];
//  for (k = 1; k <= Molecule[res].Beads[0]; k++) {
//    Aggregate[NumberOfAggregates-1].Beads[k] = Molecule[res].Beads[k];
//  }
//} //}}}
////}}}

///* test if correct number of residues in all aggregates */ //{{{
//if (l != residues) {
//  printf("ERROR:     (*Counts).Molecules = %d\n", (*Counts).Molecules);
//  printf("  # of residues in vcf file = %d\n", residues);
//  printf("# of residues in aggregates = %d\n\n", l);
//  exit(1);
//} //}}}

///* bubble sort for ascending residue ids */ //{{{
//for (i = 0; i < NumberOfAggregates; i++) {
//  for (j = 1 ; j <= (Aggregate[i].Residues[0]-1); j++) {
//    id1 = 0; /* swapped? */
//    for (k = 1 ; k <= (Aggregate[i].Residues[0]-j); k++) {
//      if (Aggregate[i].Residues[k] > Aggregate[i].Residues[k+1]) {
//        id1 = Aggregate[i].Residues[k];
//        Aggregate[i].Residues[k] = Aggregate[i].Residues[k+1];
//        Aggregate[i].Residues[k+1] = id1;
//        id1 = 1;
//      }
//    }
//    if (id1 == 0)
//      break;
//  }
//} //}}}

///* bubble sort for ascending bead ids */ //{{{
//for (i = 0; i < NumberOfAggregates; i++) {
//  for (j = 1 ; j <= (Aggregate[i].Beads[0]-1); j++) {
//    id1 = 0; /* swapped? */
//    for (k = 1 ; k <= (Aggregate[i].Beads[0]-j); k++) {
//      if (Aggregate[i].Beads[k] > Aggregate[i].Beads[k+1]) {
//        id1 = Aggregate[i].Beads[k];
//        Aggregate[i].Beads[k] = Aggregate[i].Beads[k+1];
//        Aggregate[i].Beads[k+1] = id1;
//        id1 = 1;
//      }
//    }
//    if (id1 == 0)
//      break;
//  }
//} //}}}

//// join all residues; aka remove pbc //{{{
//for (i = 0; i < (*Counts).Molecules; i++)
//  RemovePBCResidue(i); //}}}

///* put together all aggregates */ //{{{
//for (i = 0; i < NumberOfAggregates; i++) {
//  testj = 0; /* just to make sure do{}while isn't infinite */

//  moved[Aggregate[i].Residues[1]] = 1; /* first residue of aggregate isn't to move */

//  /* do stuff until all residues have moved[]=1 */
//  do {
//    for (j = 1; j <= Aggregate[i].Residues[0]; j++) {
//      id1 = Aggregate[i].Residues[j];
//      // if id1 was moved, look at other residues to move
//      if (moved[id1] == 1) { // always true for first residue of aggregate
//        for (k = 1; k <= Aggregate[i].Residues[0]; k++) {
//          if (k != j) { /* look only at different residues */
//            id2 = Aggregate[i].Residues[k];

//            // number of contact between id1 and id2 //{{{
//            if (id1 > id2) // contact[][] is triangular matrix
//              l = contact[id1][id2];
//            else
//              l = contact[id2][id1]; //}}}

//            /* if id1 and id2 are supposed to be close and id2 wasn't moved,
//             * check if they're close or separated by boxlength */ //{{{
//            if (l >= contacts && moved[id2] == 0) {
//              moveid2 = Molecule[id2].Beads[1];

//              // find bead in id1 closest to the first bead in id2 //{{{
//              min_dist = 1e5; // just take any high number
//              for (l = 1; l <= Molecule[id1].Beads[0]; l++) {
//                rij.x = Bead[moveid2].Position.x - Bead[Molecule[id1].Beads[l]].Position.x; //{{{
//                rij.y = Bead[moveid2].Position.y - Bead[Molecule[id1].Beads[l]].Position.y;
//                rij.z = Bead[moveid2].Position.z - Bead[Molecule[id1].Beads[l]].Position.z; //}}}
//                rij = RijPeriodicBoundary(rij);
//                rij.x = sqrt(SQR(rij.x) + SQR(rij.y) + SQR(rij.z));
//                if (rij.x < min_dist) {
//                  min_dist = rij.x;
//                  moveid1 = Molecule[id1].Beads[l];
//                }
//              } //}}}

//              rij.x = Bead[moveid2].Position.x - Bead[moveid1].Position.x; //{{{
//              rij.y = Bead[moveid2].Position.y - Bead[moveid1].Position.y;
//              rij.z = Bead[moveid2].Position.z - Bead[moveid1].Position.z; //}}}

//              rij = RijPeriodicBoundary(rij);
//              Bead[moveid2].Position.x = Bead[moveid1].Position.x + rij.x; //{{{
//              Bead[moveid2].Position.y = Bead[moveid1].Position.y + rij.y;
//              Bead[moveid2].Position.z = Bead[moveid1].Position.z + rij.z; //}}}

//              // move the rest of the beads of Molecule[id2]
//              RemovePBCResidue(id2);

//              // residue id2 is not to move again
//              moved[id2] = 1;
//            } //}}}
//          }
//        }
//      }
//    }
//    // test if all aggregates are moved //{{{
//    testi = -1;
//    whatever = 0;
//    for (j = 1; j <= Aggregate[i].Residues[0]; j++) {
//      if (moved[Aggregate[i].Residues[j]] == 0) {
//        testi = 0;
//        whatever++;
//      }
//    }
//    if (testi == 0)
//      testj++; //}}}
//  } while (whatever != 0);

//  /* move center of aggregate's mass to the box */ //{{{
//  rij = AggCenterOfMass(i, selected);
//  if (rij.x > BoxLength.x) {
//    for (j = 1; j <= Aggregate[i].Beads[0]; j++)
//      Bead[Aggregate[i].Beads[j]].Position.x -= BoxLength.x;
//  } else if (rij.x < 0) {
//    for (j = 1; j <= Aggregate[i].Beads[0]; j++)
//      Bead[Aggregate[i].Beads[j]].Position.x += BoxLength.x;
//  }
//  if (rij.y > BoxLength.y) { //{{{
//    for (j = 1; j <= Aggregate[i].Beads[0]; j++)
//      Bead[Aggregate[i].Beads[j]].Position.y -= BoxLength.y;
//  } else if (rij.y < 0) {
//    for (j = 1; j <= Aggregate[i].Beads[0]; j++)
//      Bead[Aggregate[i].Beads[j]].Position.y += BoxLength.y;
//  } //}}}
//  if (rij.z > BoxLength.z) { //{{{
//    for (j = 1; j <= Aggregate[i].Beads[0]; j++)
//      Bead[Aggregate[i].Beads[j]].Position.z -= BoxLength.z;
//  } else if (rij.z < 0) {
//    for (j = 1; j <= Aggregate[i].Beads[0]; j++)
//      Bead[Aggregate[i].Beads[j]].Position.z += BoxLength.z;
//  } //}}}
//  //}}}
//} //}}}

  /* free memory */ //{{{
  for (int i = 0; i < (*Counts).Molecules; i++)
    free(contact[i]);
  free(contact);
  free(moved); //}}}
} //}}}

int main(int argc, char *argv[]) {

  // -h option - print help and exit //{{{
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printf("Aggregates determines which molecules belong to which aggregate\n");
      printf("on the basis of given parameters - the minimum distance at which\n");
      printf("a pair of beads from different molecules is considered a contact;\n");
      printf("the minimum number of such contacts between two molecules to consider\n");
      printf("them as belonging to the same aggregate. Only distances between\n");
      printf("specified bead types are considered. Information about aggregates\n");
      printf("in each timestep is written to .agg file. The program uses dl_meso.vsf\n");
      printf("(or other input structure file) and FIELD (along with optional\n");
      printf("bond file) files to determine all information about the system.\n\n");

      printf("Usage:\n");
      printf("   %s <input.vcf> <distance> <contacts> ", argv[0]);
      printf("<output.agg> <type names> <options>\n\n");

      printf("   <input.vcf>       input filename (vcf format)\n");
      printf("   <output.agg>      output filename (agg format)\n");
      printf("   <distance>        minimum distance for contact for aggregate check\n");
      printf("   <contacts>        minimum number of contacts for aggregate check\n");
      printf("   <type names>      names of bead types for closeness calculation\n");
      printf("   <options>\n");
      printf("      -i <name>      use input .vsf file different from dl_meso.vsf\n");
      printf("      -b <name>      file containing bond alternatives to FIELD\n");
      printf("      -v             verbose output\n");
      printf("      -V             verbose output with comments from input .vcf file\n");
      printf("      -h             print this help and exit\n");
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
    fprintf(stderr, "Non-numeric argement for <start>!\n");
    ErrorHelp(argv[0]);
    exit(1);
  }
  double distance = atof(argv[count]); //}}}

  // <contacts> - number of steps to skip per one used //{{{
  // Error - non-numeric argument
  if (argv[++count][0] < '0' || argv[count][0] > '9') {
    fprintf(stderr, "Non-numeric argement for <skip>!\n");
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
  ReadStructure(vsf_file, bonds_file, &Counts, &BeadType, &Bead, &MoleculeType, &Molecule);

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_vcf, "r")) == NULL) {
    fprintf(stderr, "Cannot open file %s!\n", input_vcf);
    exit(1);
  } //}}}

  // vcf possibly contains indexed timesteps - find which bead types are in vcf file //{{{
  int beadcount = 0; // number of beads in every indexed timestep
  int test;
  while ((test = getc(vcf)) == '#') {

    // read bead type name from vcf //{{{
    char str[32];
    if (fscanf(vcf, "%s", str) != 1) {
      fprintf(stderr, "Cannot read bead type name from the beginning\n");
      fprintf(stderr, "of %s file (should contain indexed timesteps)!\n", input_vcf);
      exit(1);
    } //}}}

    // find type id of the bead type name
    int type = FindType(str, Counts, BeadType);

    // bead type is in vcf file
    BeadType[type].Use = 1;

    beadcount += BeadType[type].Number;

    // skip the rest of line //{{{
    while (getc(vcf) != '\n')
      ; //}}}

    // verbose output //{{{
    if (verbose) {
      printf("   in vcf file: %s\n", str);
    } //}}}
  }
  ungetc(test, vcf);
  //}}}

  // vcf contains ordered timesteps? If so, all bead types are in vcf file //{{{
  if (beadcount == 0) {
    beadcount = Counts.Bonded + Counts.Unbonded;

    // all bead types are in vcf file
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      BeadType[i].Use = 1;
    }

    if (verbose) {
      printf("   in vcf: all bead types (ordered timesteps)\n");
    }
  } //}}}

  // <type names> - names of bead types to use for closeness calculation //{{{
  // first run through to test that every specified bead type is in vcf file as well
  for (int i = (count+1); i < argc && argv[i][0] != '-'; i++) {
    int type = FindType(argv[i], Counts, BeadType);

    if (BeadType[type].Use == 0) {
      fprintf(stderr, "Bead type %s required for calculation is not in %s coordinate file!\n", BeadType[type].Name, input_vcf);
      exit(1);
    }

    BeadType[type].Use = 2;
  } //}}}

  // print information - verbose output //{{{
  if (verbose || verbose2) {
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
      printf("Use = %d}\n", BeadType[i].Use);
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

  // read pbc from coordinate file //{{{
  char str[32];
  Vector BoxLength;
  if (fscanf(vcf, "%s %lf %lf %lf", str, &BoxLength.x, &BoxLength.y, &BoxLength.z) != 4 ||
      strcmp(str, "pbc") != 0) {
    fprintf(stderr, "Cannot read pbc from %s (should be first line)!\n", input_vcf);
    exit(1);
  }

  while (getc(vcf) != '\n')
    ;
  while (getc(vcf) != '\n')
    ;

  // print pbc if verbose output
  if (verbose) {
    printf("   box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } //}}}

  // create array for the first line of a timestep ('# <number and/or other comment>') //{{{
  char *stuff;
  stuff = malloc(128*sizeof(int));

  // initialize the array
  for (int i = 0; i < 128; i++) {
    stuff[i] = '\0';
  } //}}}

  // allocate Aggregate struct //{{{
  Aggregate *Aggregate = malloc(Counts.Molecules*sizeof(*Aggregate));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate - memory-heavy, but reliable
    Aggregate[i].Monomer = malloc(Counts.Unbonded*sizeof(int));
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = malloc(Counts.Molecules*sizeof(int));
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  while ((test = getc(vcf)) != EOF) {
    ungetc(test, vcf);

    fflush(stdout);
    printf("\rStep: %6d", ++count);

    // either read indexed timestep from input .vcf file //{{{
    if (beadcount != (Counts.Bonded+Counts.Unbonded)) {
      if ((test = ReadCoorIndexed(vcf, beadcount, &Bead, &stuff)) != 0) {
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

    CalculateAggregates(&Aggregate, &Counts, BoxLength, SQR(distance), BeadType, Bead, MoleculeType, Molecule);

    // open output .agg file for appending //{{{
    if ((out = fopen(output_agg, "a")) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", output_agg);
      exit(1);
    } //}}}

    // do stuff

    fclose(out);

    // print comment at the beginning of a timestep - detailed verbose output //{{{
    if (verbose2) {
      printf("\n%s", stuff);
    } //}}}
  }

  fflush(stdout);
  printf("\rLast Step: %6d\n", count);

  fclose(vcf); //}}}

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
    free(Aggregate[i].Monomer);
    free(Aggregate[i].Molecule);
  }
  free(Molecule);
  free(Aggregate);
  free(stuff);
  //}}}

  return 0;
}
