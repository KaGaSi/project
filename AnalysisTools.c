#include "AnalysisTools.h"

// VerboseOutput() //{{{
/**
 * Function providing standard verbose output (for cases when verbose
 * option is used). It prints most of the information about used system.
 */
void VerboseOutput(char *input_vcf, Counts Counts, Vector BoxLength,
                   BeadType *BeadType, Bead *Bead,
                   MoleculeType *MoleculeType, Molecule *Molecule) {

  putchar('\n');
  if (BoxLength.x != -1) {
    fprintf(stdout, "Box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } else {
    fprintf(stdout, "Unknown box size because no coordinate file is provided.\n\n");
  }
  PrintCounts(Counts);
  PrintBeadType(Counts, BeadType);
  PrintMoleculeType(Counts, BeadType, MoleculeType);
  putchar('\n');
} //}}}

// PrintCounts()  //{{{
/**
 * Function printing Counts structure.
 */
void PrintCounts(Counts Counts) {
  fprintf(stdout, "Counts.{");
  fprintf(stdout, "TypesOfBeads = %d, ", Counts.TypesOfBeads);
  fprintf(stdout, "Bonded = %d, ", Counts.Bonded);
  fprintf(stdout, "Unbonded = %d, ", Counts.Unbonded);
  fprintf(stdout, "Beads = %d, ", Counts.Beads);
  fprintf(stdout, "BeadsInVsf = %d, ", Counts.BeadsInVsf);
  fprintf(stdout, "TypesOfMolecules = %d, ", Counts.TypesOfMolecules);
  fprintf(stdout, "Molecules = %d}\n\n", Counts.Molecules);
} //}}}

// PrintBeadType()  //{{{
/**
 * Function printing BeadType structure.
 */
void PrintBeadType(Counts Counts, BeadType *BeadType) {
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(stdout, "BeadType[%2d].{", i);
    fprintf(stdout, "Name =%10s, ", BeadType[i].Name);
    fprintf(stdout, "Number =%7d, ", BeadType[i].Number);
    fprintf(stdout, "Charge =%9.5f, ", BeadType[i].Charge);
    fprintf(stdout, "Mass = %.5f}\n", BeadType[i].Mass);
//  fprintf(stdout, "Use = %3s, ", BeadType[i].Use? "Yes":"No");
//  fprintf(stdout, "Write = %3s}\n", BeadType[i].Write? "Yes":"No");
  }
  putchar('\n');
} //}}}

// PrintMoleculeType()  //{{{
/**
 * Function printing MoleculeType structure.
 */
void PrintMoleculeType(Counts Counts, BeadType *BeadType, MoleculeType *MoleculeType) {
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(stdout, "MoleculeType[%2d].{\n", i);
    fprintf(stdout, "  Name    = %s\n", MoleculeType[i].Name);
    fprintf(stdout, "  Number  = %d\n", MoleculeType[i].Number);
    fprintf(stdout, "  nBeads  = %d\n  Bead    = {", MoleculeType[i].nBeads);
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].Bead[j]].Name);
    }
    fprintf(stdout, "}\n  nBonds  = %d\n  Bond    = {", MoleculeType[i].nBonds);
    for (int j = 0; j < MoleculeType[i].nBonds; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%d-%d", MoleculeType[i].Bond[j][0]+1, MoleculeType[i].Bond[j][1]+1);
    }
    fprintf(stdout, "}\n  nBTypes = %d\n  BType   = {", MoleculeType[0].nBTypes);
    for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].BType[j]].Name);
    }
    fprintf(stdout, "}\n  Mass    = %.5f}\n", MoleculeType[i].Mass);
  }
} //}}}

// PrintBead() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead(Counts Counts, int *Index, BeadType *BeadType, Bead *Bead) {
  fprintf(stdout, "Beads - <i> (<Bead[i].Index>; <Index[i]>)\n");
  for (int i = 0; i < Counts.Beads; i++) {
    int type = Bead[i].Type;
    fprintf(stdout, "   %6d (%6d; %6d) %8s molecule: ", i, Bead[i].Index, Index[i], BeadType[type].Name);
    if (Bead[i].Molecule == -1) {
      fprintf(stdout, "None\n");
    } else {
      fprintf(stdout, "%6d\n", Bead[i].Molecule+1);
    }
  }
} //}}}

// PrintMolecule() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule(Counts Counts, int *Index,
                   MoleculeType *MoleculeType, Molecule *Molecule,
                   BeadType *BeadType, Bead *Bead) {
  fprintf(stdout, "Molecules\n");
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    fprintf(stdout, "Molecule %3d (%s):\n", i+1, MoleculeType[type].Name);
    fprintf(stdout, " BEAD INDICES (%d): internal (in vsf)\n", MoleculeType[type].nBeads);
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      fprintf(stdout, "   %d (%d)\n", Molecule[i].Bead[j], Bead[Molecule[i].Bead[j]].Index);
    }
    fprintf(stdout, " BONDS (%d): internal (in vsf)\n", MoleculeType[type].nBonds);
    for (int j = 0; j < MoleculeType[type].nBonds; j++) {
      int bead1 = MoleculeType[type].Bond[j][0];
      int bead1_1 = Molecule[i].Bead[bead1];
      int bead2 = MoleculeType[type].Bond[j][1];
      int bead2_1 = Molecule[i].Bead[bead2];
      fprintf(stdout, "   %d-%d (%d-%d)\n", bead1_1, bead2_1, Bead[bead1_1].Index, Bead[bead2_1].Index);
    }
  }
  fprintf(stdout, "\n");
} //}}}

// PrintAggregate() //{{{
/**
 * Function printing Aggregate structure.
 */
void PrintAggregate(Counts Counts, int *Index, MoleculeType *MoleculeType, Molecule *Molecule, Bead *Bead, BeadType *BeadType, Aggregate *Aggregate) {
  fprintf(stdout, "Aggregates: %d\n", Counts.Aggregates);
  for (int i = 0; i < Counts.Aggregates; i++) {
    // print molecules
    fprintf(stdout, " %d mols:", Aggregate[i].nMolecules);
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol = Aggregate[i].Molecule[j];
      int type = Molecule[mol].Type;
      fprintf(stdout, " %d (%d)", mol, type);
      if (j != (Aggregate[i].nMolecules-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
    // print bonded beads
    fprintf(stdout, " %d bonded beads:", Aggregate[i].nBeads);
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      fprintf(stdout, " %d (%d)", bead, Bead[bead].Index);
      if (j != (Aggregate[i].nBeads-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
    // print monomeric beads
    fprintf(stdout, " %d free beads:", Aggregate[i].nMonomers);
    for (int j = 0; j < Aggregate[i].nMonomers; j++) {
      int bead = Aggregate[i].Monomer[j];
      fprintf(stdout, " %d (%d)", bead, Bead[bead].Index);
      if (j != (Aggregate[i].nMonomers-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
  }
} //}}}

// FindBeadType() //{{{
/* Function to identify type of bead from its name; returns -1 on non-existent
 * bead name.
 */
int FindBeadType(char *name, Counts Counts, BeadType *BeadType) {
  int type;

  // compare give 'name' with all known bead types & return bead type id
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (strcmp(name, BeadType[i].Name) == 0) {
      type = i;
      return type;
    }
  }

  // name isn't in BeadType struct
  return -1;
} //}}}

// FindMoleculeType() //{{{
int FindMoleculeType(char *name, Counts Counts, MoleculeType *MoleculeType) {
  int type;

  // compare give 'name' with all known bead types & return bead type id
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (strcmp(name, MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }

  // name isn't in MoleculeType struct
  return(-1);
} //}}}

// RemovePBCMolecules() //{{{
/**
 * Function to remove periodic boundary conditions from all individual
 * molecules, thus joining them
 */
void RemovePBCMolecules(Counts Counts, Vector BoxLength,
                        BeadType *BeadType, Bead **Bead,
                        MoleculeType *MoleculeType, Molecule *Molecule) {

  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      (*Bead)[Molecule[i].Bead[j]].Flag = false; // no beads moved yet
    }
    (*Bead)[Molecule[i].Bead[MoleculeType[type].Bond[0][0]]].Flag = true;
    bool done = false;
    while (!done) {
      for (int j = 0; j < MoleculeType[type].nBonds; j++) {
        int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];

        // move id1, if id2 is moved already
        if (!(*Bead)[id1].Flag && (*Bead)[id2].Flag) {
          Vector dist = Distance((*Bead)[id2].Position, (*Bead)[id1].Position, BoxLength);

          (*Bead)[id1].Position.x = (*Bead)[id2].Position.x - dist.x;
          (*Bead)[id1].Position.y = (*Bead)[id2].Position.y - dist.y;
          (*Bead)[id1].Position.z = (*Bead)[id2].Position.z - dist.z;
          (*Bead)[id1].Flag = true;
        // move id2, if id1 was moved already
        } else if ((*Bead)[id1].Flag && !(*Bead)[id2].Flag) {
          Vector dist = Distance((*Bead)[id1].Position, (*Bead)[id2].Position, BoxLength);

          (*Bead)[id2].Position.x = (*Bead)[id1].Position.x - dist.x;
          (*Bead)[id2].Position.y = (*Bead)[id1].Position.y - dist.y;
          (*Bead)[id2].Position.z = (*Bead)[id1].Position.z - dist.z;
          (*Bead)[id2].Flag = true;
        }
      }

      done = true;
      for (int j = 1; j < MoleculeType[type].nBeads; j++) {
        if (!(*Bead)[Molecule[i].Bead[j]].Flag) {
          done = false;
          break;
        }
      }
    }
  }
} //}}}

// RemovePBCAggregates() //{{{
/**
 * Function to remove periodic boundary conditions from all aggregates,
 * thus joining them.
 */
void RemovePBCAggregates(double distance, Aggregate *Aggregate, Counts Counts,
                         Vector BoxLength, BeadType *BeadType, Bead **Bead,
                         MoleculeType *MoleculeType, Molecule *Molecule) {

  bool *moved = malloc(Counts.Molecules*sizeof(bool));

  // go through all aggregates larger than unimers and put all molecules together //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {

    // negate moved array, but the first molecule is not to move //{{{
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
            int mol1 = Aggregate[i].Molecule[j],
                mol2 = Aggregate[i].Molecule[k],
                mol1_type = Molecule[mol1].Type,
                mol2_type = Molecule[mol2].Type;

            // go through all bead pairs in the two molecules
            for (int l = 0; l < MoleculeType[mol1_type].nBeads; l++) {
              for (int m = 0; m < MoleculeType[mol2_type].nBeads; m++) {
                int bead1 = Molecule[mol1].Bead[l];
                int bead2 = Molecule[mol2].Bead[m];

                // use only bead types that were used to assign molecules to aggregates
                if (BeadType[(*Bead)[bead1].Type].Use &&
                    BeadType[(*Bead)[bead2].Type].Use) {

                  // calculate distance between 'bead1' and 'bead2'
                  Vector dist = Distance((*Bead)[bead1].Position, (*Bead)[bead2].Position, BoxLength);
                  dist.x = Length(dist);

                  // move 'mol2' (or 'k') if 'bead1' and 'bead2' are in contact
                  if (dist.x <= distance) {

                    // distance vector between 'bead1' and 'bead2' //{{{
                    dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z; //}}}
                    // if 'bead1' and 'bead2' are too far in x-direction, move 'mol2' in x-direction //{{{
                    while (dist.x > (BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x += BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    }
                    while (dist.x <= -(BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x -= BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    } //}}}
                    // if 'bead1' and 'bead2' are too far in y-direction, move 'mol2' in y-direction //{{{
                    while (dist.y > (BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y += BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    }
                    while (dist.y <= -(BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y -= BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    } //}}}
                    // if 'bead1' and 'bead2' are too far in z-direction, move 'mol2' in x-direction //{{{
                    while (dist.z > (BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z += BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    }
                    while (dist.z <= -(BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z -= BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    } //}}}

                    moved[k] = true;

                    // skip remainder of 'mol2' (or 'k')
                    break;
                  }
                }
              }
              // if molekule 'k' (or 'mol2') has been moved, skip also remainder of molecules 'mol1'
              if (moved[k]) {
                break;
              }
            }
          }
        }
      }

      // check if all molecules have moved //{{{
      done = true;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (!moved[j]) {
          done = false;
          break;
        }
      } //}}}
    }
  }
  free(moved); //}}}

  // put aggregates' centre of mass into the simulation box //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    Vector com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead,
                              *Bead, BeadType);

    // by how many BoxLength's should com by moved?
    // for distant aggregates - it shouldn't happen, but better safe than sorry
    IntVector move;
    move.x = com.x / BoxLength.x;
    move.y = com.y / BoxLength.y;
    move.z = com.z / BoxLength.z;
    if (com.x < 0) {
      move.x--;
    }
    if (com.y < 0) {
      move.y--;
    }
    if (com.z < 0) {
      move.z--;
    }

    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      (*Bead)[bead].Position.x -= move.x * BoxLength.x;
      (*Bead)[bead].Position.y -= move.y * BoxLength.y;
      (*Bead)[bead].Position.z -= move.z * BoxLength.z;
    }
  } //}}}
} //}}}

// RestorePBC() //{{{
/**
 * Function to restore removed periodic boundary conditions. Used also in case
 * of cell linked lists, because they need coordinates <0, BoxLength>.
 */
void RestorePBC(Counts Counts, Vector BoxLength, Bead **Bead) {

  for (int i = 0; i < (Counts.Unbonded+Counts.Bonded); i++) {
    // x direction
    while ((*Bead)[i].Position.x >= BoxLength.x) {
      (*Bead)[i].Position.x -= BoxLength.x;
    }
    while ((*Bead)[i].Position.x < 0) {
      (*Bead)[i].Position.x += BoxLength.x;
    }
    // y direction
    while ((*Bead)[i].Position.y >= BoxLength.y) {
      (*Bead)[i].Position.y -= BoxLength.y;
    }
    while ((*Bead)[i].Position.y < 0) {
      (*Bead)[i].Position.y += BoxLength.y;
    }
    // z direction
    while ((*Bead)[i].Position.z >= BoxLength.z) {
      (*Bead)[i].Position.z -= BoxLength.z;
    }
    while ((*Bead)[i].Position.z < 0) {
      (*Bead)[i].Position.z += BoxLength.z;
    }
  }
} //}}}

// CentreOfMass() //{{{
/**
 * Function to calculate centre of mass for a given list of beads.
 */
Vector CentreOfMass(int n, int *list, Bead *Bead, BeadType *BeadType) {

  Vector com;
  com.x = 0;
  com.y = 0;
  com.z = 0;

  for (int i = 0; i < n; i++) {
    com.x += Bead[list[i]].Position.x * BeadType[Bead[list[i]].Type].Mass;
    com.y += Bead[list[i]].Position.y * BeadType[Bead[list[i]].Type].Mass;
    com.z += Bead[list[i]].Position.z * BeadType[Bead[list[i]].Type].Mass;
  }
  com.x /= n;
  com.y /= n;
  com.z /= n;

  return com;
} //}}}

// Gyration() //{{{
/**
 * Function to calculate the principle moments of the gyration tensor.
 */
Vector Gyration(int n, int *list, Counts Counts, Vector BoxLength, BeadType *BeadType, Bead **Bead) {
  // gyration tensor (3x3 array) //{{{
  // use long double to ensure precision -- previous problem with truncation in short chains
  struct Tensor {
    LongVector x, y, z;
  } GyrationTensor;

  GyrationTensor.x.x = 0;
  GyrationTensor.x.y = 0;
  GyrationTensor.x.z = 0;
  GyrationTensor.y.x = 0;
  GyrationTensor.y.y = 0;
  GyrationTensor.y.z = 0;
  GyrationTensor.z.x = 0;
  GyrationTensor.z.y = 0;
  GyrationTensor.z.z = 0; //}}}

// test print of given coordinates -- uncomment if need be //{{{
//for (int i = 0; i < n; i++) {
//  fprintf(stderr, " %10.5f %10.5f %10.5f \n", (*Bead)[list[i]].Position.x,
//                                              (*Bead)[list[i]].Position.y,
//                                              (*Bead)[list[i]].Position.z);
//} //}}}

  Vector com = CentreOfMass(n, list, *Bead, BeadType);
//fprintf(stderr, "%lf %lf %lf\n", com.x, com.y, com.z);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= com.x;
    (*Bead)[list[i]].Position.y -= com.y;
    (*Bead)[list[i]].Position.z -= com.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    GyrationTensor.x.x += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.x;
    GyrationTensor.x.y += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.y;
    GyrationTensor.x.z += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.z;
    GyrationTensor.y.y += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.y;
    GyrationTensor.y.z += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.z;
    GyrationTensor.z.z += (*Bead)[list[i]].Position.z * (*Bead)[list[i]].Position.z;
  }
  GyrationTensor.x.x /= n;
  GyrationTensor.x.y /= n;
  GyrationTensor.x.z /= n;
  GyrationTensor.y.y /= n;
  GyrationTensor.y.z /= n;
  GyrationTensor.z.z /= n; //}}}
//fprintf(stderr, "Tensor: (%lf, %lf, %lf)\n", GyrationTensor.x.x, GyrationTensor.x.y, GyrationTensor.x.z);
//fprintf(stderr, "        (%lf, %lf, %lf)\n", GyrationTensor.x.y, GyrationTensor.y.y, GyrationTensor.y.z);
//fprintf(stderr, "        (%lf, %lf, %lf)\n", GyrationTensor.x.z, GyrationTensor.y.z, GyrationTensor.z.z);

  // char polynomial: a_cube * x^3 + b_cube * x^2 + c_cube * x + d_cube = 0 //{{{
  long double a_cube = -1;
  long double b_cube = GyrationTensor.x.x + GyrationTensor.y.y + GyrationTensor.z.z;
  long double c_cube = - GyrationTensor.x.x * GyrationTensor.y.y
                  - GyrationTensor.x.x * GyrationTensor.z.z
                  - GyrationTensor.y.y * GyrationTensor.z.z
                  + SQR(GyrationTensor.y.z)
                  + SQR(GyrationTensor.x.y)
                  + SQR(GyrationTensor.x.z);
  long double d_cube = + GyrationTensor.x.x * GyrationTensor.y.y * GyrationTensor.z.z
                  + 2 * GyrationTensor.x.y * GyrationTensor.y.z * GyrationTensor.x.z
                  - SQR(GyrationTensor.x.z) * GyrationTensor.y.y
                  - SQR(GyrationTensor.x.y) * GyrationTensor.z.z
                  - SQR(GyrationTensor.y.z) * GyrationTensor.x.x; //}}}
//fprintf(stderr, "character: %lfx^3 + %lfx^2 + %lfx^1 + %lfx^0;\n", a_cube, b_cube, c_cube, d_cube);

  // first root: either 0 or Newton's iterative method to get it //{{{
  long double root0 = 0;
  if (fabs(d_cube) > 0.0000000001L) {
    // derivative of char. polynomial: a_deriv * x^2 + b_deriv * x + c_deriv
    long double a_deriv = 3 * a_cube;
    long double b_deriv = 2 * b_cube;
    long double c_deriv = c_cube;

    long double root1 = 1;

    while (fabs(root0-root1) > 0.0000000001L) {
      long double f_root0 = (a_cube * CUBE(root0) + b_cube * SQR(root0) + c_cube * root0 + d_cube);
      long double f_deriv_root0 = (a_deriv * SQR(root0) + b_deriv * root0 + c_deriv);
      root1 = root0 - f_root0 / f_deriv_root0;

      // swap root0 and root1 for the next iteration
      long double tmp = root0;
      root0 = root1;
      root1 = tmp;
    }
  } //}}}
//fprintf(stderr, "root0=%lf; ", root0);

  // determine paremeters of quadratic equation a_quad * x^2 + b_quad * x + c_quad = 0 //{{{
  // derived by division: (x^3 + (b_cube/a_cube) * x^2 + (c_cube/a_cube) * x + (d_cube/a_cube)):(x - root0)
  long double a_quad = 1;
  long double b_quad = b_cube / a_cube + root0;
  long double c_quad = SQR(root0) + b_cube / a_cube * root0 + c_cube/a_cube; //}}}
//fprintf(stderr, "quad: %lfx^2 + %lfx + %lf; ", a_quad, b_quad, c_quad);

  // calculate & sort eigenvalues //{{{
  LongVector eigen;
  eigen.x = root0; // found out by Newton's method
  // roots of the quadratic equation
  eigen.y = (-b_quad + sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);
  eigen.z = (-b_quad - sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);

  Vector eigen2;
  eigen2.x = eigen.x;
  eigen2.y = eigen.y;
  eigen2.z = eigen.z;
  eigen2 = Sort3(eigen2); //}}}
//fprintf(stderr, "eigen=(%lf, %lf, %lf)\n", eigen.x, eigen.y, eigen.z);

  return eigen2;
} //}}}

// EvaluateContacts() //{{{
/**
 * Function evaluating contacts for aggregate detection for Aggregate and
 * Aggregate-NotSameBeads utilities.
 */
void EvaluateContacts(Counts *Counts, Aggregate **Aggregate,
                      Molecule **Molecule,
                      int contacts, int **contact) {
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
  }

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
} //}}}

// SortAggStruct() //{{{
/**
 * Sort an Aggregate struct using the bubble sort algorithm. The resulting
 * struct is arranged so that aggregates with the first molecule's lower id
 * come first.
 *
 * TODO: why doesn't it sort Aggregate[].Monomers? ...because there aren't any yet
 */
void SortAggStruct(Aggregate **Aggregate, Counts Counts,
                   Molecule *Molecule, MoleculeType *MoleculeType,
                   Bead **Bead, BeadType *BeadType) {
  for (int i = 0; i < (Counts.Aggregates-1); i++) {
    bool done = true;
    for (int j = 0; j < (Counts.Aggregates-i-1); j++) {
      if ((*Aggregate)[j].Molecule[0] > (*Aggregate)[j+1].Molecule[0]) {
        Swap(&(*Aggregate)[j].nMolecules, &(*Aggregate)[j+1].nMolecules);
        // switch the whole Aggregate[].Molecule array
        int mols; // number of molecules in the larger aggregate
        if ((*Aggregate)[j].nMolecules > (*Aggregate)[j+1].nMolecules) {
          mols = (*Aggregate)[j].nMolecules;
        } else {
          mols = (*Aggregate)[j+1].nMolecules;
        }
        for (int k = 0; k < mols; k++) {
          Swap(&(*Aggregate)[j].Molecule[k], &(*Aggregate)[j+1].Molecule[k]);
        }
        // switch bonded beads array
        Swap(&(*Aggregate)[j].nBeads, &(*Aggregate)[j+1].nBeads);
        int beads; // number of beads in the larger aggregate
        if ((*Aggregate)[j].nBeads > (*Aggregate)[j+1].nBeads) {
          beads = (*Aggregate)[j].nBeads;
        } else {
          beads = (*Aggregate)[j+1].nBeads;
        }
        for (int k = 0; k < beads; k++) {
          Swap(&(*Aggregate)[j].Bead[k], &(*Aggregate)[j+1].Bead[k]);
        }
        // switch monomer beads array
        Swap(&(*Aggregate)[j].nMonomers, &(*Aggregate)[j+1].nMonomers);
        int mons; // larger number of monomers of the two aggregates
        if ((*Aggregate)[j].nMonomers > (*Aggregate)[j+1].nMonomers) {
          mons = (*Aggregate)[j].nMonomers;
        } else {
          mons = (*Aggregate)[j+1].nMonomers;
        }
        for (int k = 0; k < mons; k++) {
          Swap(&(*Aggregate)[j].Monomer[k], &(*Aggregate)[j+1].Monomer[k]);
        }
        done = false;
      }
    }
    if (done)
      break;
  }

  // TODO: test that it's working correctly be prints in Aggregate's main()
  // re-assign aggregate id to every bonded bead in the aggregate, correcting after sorting //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];
      int mtype = Molecule[mol].Type;
      for (int k = 0; k < MoleculeType[mtype].nBeads; k++) {
        int id = Molecule[mol].Bead[k];
        (*Bead)[id].nAggregates = 1;
        (*Bead)[id].Aggregate[0] = i;
      }
    }
  } //}}}
} //}}}

// LinkedList() //{{{
void LinkedList(Vector BoxLength, Counts Counts, Bead *Bead,
                int **Head, int **Link, double cell_size, IntVector *n_cells,
                int *Dcx, int *Dcy, int *Dcz) {

  (*n_cells).x = ceil(BoxLength.x/cell_size),
  (*n_cells).y = ceil(BoxLength.y/cell_size),
  (*n_cells).z = ceil(BoxLength.z/cell_size);

  // allocate arrays
  *Head = malloc(((*n_cells).x*(*n_cells).y*(*n_cells).z)*sizeof(int));
  *Link = malloc(Counts.Beads*sizeof(int));
  for (int i = 0; i < ((*n_cells).x*(*n_cells).y*(*n_cells).z); i++) {
    (*Head)[i] = -1;
  }

  // sort beads into cells //{{{
  for (int i = 0; i < Counts.Beads; i++) {
    int cell = (int)(Bead[i].Position.x / cell_size)
             + (int)(Bead[i].Position.y / cell_size) * (*n_cells).x
             + (int)(Bead[i].Position.z / cell_size) * (*n_cells).x * (*n_cells).y;
    (*Link)[i] = (*Head)[cell];
    (*Head)[cell] = i;
  } //}}}

  // coordinates of adjoining cells //{{{
  int x[14] = {0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
  int y[14] = {0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
  int z[14] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (int i = 0; i < 14; i++) {
    Dcx[i] = x[i];
    Dcy[i] = y[i];
    Dcz[i] = z[i];
  }
  //}}}
} //}}}

// SortBonds() //{{{
/**
 * Function to sort a bond array. As each bond contains a 2-member array, first
 * check that first bead id is lower than the second. Then sort the bonds
 * according to the index of the first bead in each bond.
 */
void SortBonds(int **bond, int length) {
  // first, check order in every bond
  for (int j = 0; j < length; j++) {
    if (bond[j][0] > bond[j][1]) {
      Swap(&bond[j][0], &bond[j][1]);
    }
  }
  // second, bubble sort bonds
  for (int j = 0; j < (length-1); j++) {
    bool swap = false;
    for (int k = 0; k < (length-j-1); k++) {
      if ((bond[k][0] > bond[k+1][0]) || // swap if first beads are in wrong order
          (bond[k][0] == bond[k+1][0] && // or if they're the same, but second ones are in wrong order
          bond[k][1] > bond[k+1][1])) {
        Swap(&bond[k][0], &bond[k+1][0]);
        Swap(&bond[k][1], &bond[k+1][1]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}

// ZeroCounts() //{{{
/**
 * Zeroize Counts struct
 */
void ZeroCounts(Counts *Counts) {
  (*Counts).TypesOfBeads = 0;
  (*Counts).TypesOfMolecules = 0;
  (*Counts).Beads = 0;
  (*Counts).Bonded = 0;
  (*Counts).Unbonded = 0;
  (*Counts).BeadsInVsf = 0;
  (*Counts).Molecules = 0;
  (*Counts).Aggregates = 0;
} //}}}

// FreeBead() //{{{
/**
 * Free memory allocated for Bead struct array. This function makes it
 * easier to add other arrays to the Bead struct in the future
 */
void FreeBead(Counts Counts, Bead **Bead) {
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    free((*Bead)[i].Aggregate);
  }
  free(*Bead);
} //}}}

// FreeMolecule() //{{{
/**
 * Free memory allocated for Molecule struct array. This function makes it
 * easier other arrays to the Molecule struct in the future
 */
void FreeMolecule(Counts Counts, Molecule **Molecule) {
  for (int i = 0; i < Counts.Molecules; i++) {
    free((*Molecule)[i].Bead);
  }
  free(*Molecule);
} //}}}

// FreeMoleculeType() //{{{
/**
 * Free memory allocated for MoleculeType struct array. This function makes
 * it easier other arrays to the MoleculeType struct in the future
 */
void FreeMoleculeType(Counts Counts, MoleculeType **MoleculeType) {
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      free((*MoleculeType)[i].Bond[j]);
    }
    free((*MoleculeType)[i].Bead);
    free((*MoleculeType)[i].Bond);
    free((*MoleculeType)[i].BType);
  }
  free(*MoleculeType);
} //}}}

// FreeAggregate() //{{{
/**
 * Free memory allocated for Aggregate struct array. This function makes it
 * easier other arrays to the Aggregate struct in the future
 */
void FreeAggregate(Counts Counts, Aggregate **Aggregate) {
  for (int i = 0; i < Counts.Molecules; i++) {
    free((*Aggregate)[i].Molecule);
    free((*Aggregate)[i].Bead);
    free((*Aggregate)[i].Monomer);
  }
  free(*Aggregate);
} //}}}
