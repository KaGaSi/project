/* This program creates a nanoparticle of a given radius. The inner part is
 * formed by a C60 fullerene with diameter of 0.5 of the nanoparticle
 * radius. The shell is formed by a C60 fullerene with added beads in the
 * middle of every C60 face. There 153 beads in the nanoparticle.
 *
 * Bonds are created between all beads, that is 153 * 152 / 2 = 11628 bonds
 * in all.
 *
 * Two files are created: C60.vtf to check how the nanoparticle looks and
 * Nano-FIELD with coordinates and bonds in the format required by the
 * FIELD file in DL_MESO simulation package.
 *
 * Usage:
 * C60 <radius>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

int main(int argc, char *argv[]) {

  double C60[60][3],
         E[90][3],
         H_temp[60][3], H[20][3],
         P_temp[60][3], P[12][3],
         Shell[92][3],
         bond = 5;

  int beads = 1 + 60 + 60 + 20 + 12;

  // C60 //{{{
  // coordinates //{{{
  C60[  0][0] =  2.16650;
  C60[  0][1] =  0.59060;
  C60[  0][2] =  2.58740;
  C60[  1][0] =  3.03780;
  C60[  1][1] =  0.17660;
  C60[  1][2] =  1.59180;
  C60[  2][0] =  1.27860;
  C60[  2][1] = -0.30980;
  C60[  2][2] =  3.16790;
  C60[  3][0] =  3.01180;
  C60[  3][1] = -1.14340;
  C60[  3][2] =  1.16540;
  C60[  4][0] =  3.10340;
  C60[  4][1] = -1.43350;
  C60[  4][2] = -0.19300;
  C60[  5][0] =  3.15030;
  C60[  5][1] =  1.21060;
  C60[  5][2] =  0.66820;
  C60[  6][0] =  3.24280;
  C60[  6][1] =  0.91490;
  C60[  6][2] = -0.68590;
  C60[  7][0] =  3.21920;
  C60[  7][1] = -0.40230;
  C60[  7][2] = -1.12070;
  C60[  8][0] = -0.43930;
  C60[  8][1] =  1.35270;
  C60[  8][2] =  3.12710;
  C60[  9][0] =  0.43630;
  C60[  9][1] =  2.26180;
  C60[  9][2] =  2.55420;
  C60[ 10][0] = -0.02960;
  C60[ 10][1] =  0.06330;
  C60[ 10][2] =  3.43790;
  C60[ 11][0] =  1.74420;
  C60[ 11][1] =  1.87900;
  C60[ 11][2] =  2.28300;
  C60[ 12][0] =  2.35190;
  C60[ 12][1] =  2.26760;
  C60[ 12][2] =  1.09900;
  C60[ 13][0] = -0.26330;
  C60[ 13][1] =  3.02680;
  C60[ 13][2] =  1.63260;
  C60[ 14][0] =  0.33740;
  C60[ 14][1] =  3.40540;
  C60[ 14][2] =  0.43730;
  C60[ 15][0] =  1.65160;
  C60[ 15][1] =  3.02780;
  C60[ 15][2] =  0.17070;
  C60[ 16][0] = -2.09030;
  C60[ 16][1] = -0.82250;
  C60[ 16][2] =  2.59550;
  C60[ 17][0] = -2.51110;
  C60[ 17][1] =  0.46640;
  C60[ 17][2] =  2.28540;
  C60[ 18][0] = -0.84490;
  C60[ 18][1] = -1.02520;
  C60[ 18][2] =  3.17380;
  C60[ 19][0] = -1.68740;
  C60[ 19][1] =  1.55330;
  C60[ 19][2] =  2.55120;
  C60[ 20][0] = -1.58430;
  C60[ 20][1] =  2.58580;
  C60[ 20][2] =  1.63190;
  C60[ 21][0] = -3.23140;
  C60[ 21][1] =  0.40610;
  C60[ 21][2] =  1.10070;
  C60[ 22][0] = -3.12270;
  C60[ 22][1] =  1.44100;
  C60[ 22][2] =  0.17460;
  C60[ 23][0] = -2.29470;
  C60[ 23][1] =  2.52910;
  C60[ 23][2] =  0.43990;
  C60[ 24][0] = -0.49080;
  C60[ 24][1] = -2.91330;
  C60[ 24][2] =  1.73650;
  C60[ 25][0] = -1.74300;
  C60[ 25][1] = -2.71240;
  C60[ 25][2] =  1.16370;
  C60[ 26][0] = -0.03930;
  C60[ 26][1] = -2.06840;
  C60[ 26][2] =  2.74530;
  C60[ 27][0] = -2.54860;
  C60[ 27][1] = -1.66500;
  C60[ 27][2] =  1.59420;
  C60[ 28][0] = -3.26020;
  C60[ 28][1] = -0.91410;
  C60[ 28][2] =  0.67010;
  C60[ 29][0] = -1.65430;
  C60[ 29][1] = -3.00610;
  C60[ 29][2] = -0.18970;
  C60[ 30][0] = -2.35420;
  C60[ 30][1] = -2.24390;
  C60[ 30][2] = -1.11700;
  C60[ 31][0] = -3.16430;
  C60[ 31][1] = -1.19490;
  C60[ 31][2] = -0.68780;
  C60[ 32][0] =  2.13640;
  C60[ 32][1] = -2.05530;
  C60[ 32][2] =  1.73580;
  C60[ 33][0] =  1.68950;
  C60[ 33][1] = -2.90090;
  C60[ 33][2] =  0.72930;
  C60[ 34][0] =  1.27850;
  C60[ 34][1] = -1.63660;
  C60[ 34][2] =  2.74350;
  C60[ 35][0] =  0.36780;
  C60[ 35][1] = -3.33270;
  C60[ 35][2] =  0.73020;
  C60[ 36][0] = -0.34400;
  C60[ 36][1] = -3.39040;
  C60[ 36][2] = -0.45940;
  C60[ 37][0] =  2.28890;
  C60[ 37][1] = -2.52500;
  C60[ 37][2] = -0.46400;
  C60[ 38][0] =  1.57900;
  C60[ 38][1] = -2.57180;
  C60[ 38][2] = -1.65800;
  C60[ 39][0] =  0.25600;
  C60[ 39][1] = -3.00540;
  C60[ 39][2] = -1.65310;
  C60[ 40][0] = -2.18280;
  C60[ 40][1] = -0.57830;
  C60[ 40][2] = -2.59790;
  C60[ 41][0] = -1.74800;
  C60[ 41][1] = -1.86940;
  C60[ 41][2] = -2.30830;
  C60[ 42][0] = -0.43850;
  C60[ 42][1] = -2.24690;
  C60[ 42][2] = -2.58450;
  C60[ 43][0] = -1.28150;
  C60[ 43][1] =  0.31890;
  C60[ 43][2] = -3.16710;
  C60[ 44][0] = -2.15260;
  C60[ 44][1] =  2.05450;
  C60[ 44][2] = -1.73780;
  C60[ 45][0] = -3.04850;
  C60[ 45][1] =  1.15350;
  C60[ 45][2] = -1.18110;
  C60[ 46][0] = -3.06560;
  C60[ 46][1] = -0.16290;
  C60[ 46][2] = -1.61070;
  C60[ 47][0] = -1.26610;
  C60[ 47][1] =  1.64070;
  C60[ 47][2] = -2.72710;
  C60[ 48][0] =  0.50390;
  C60[ 48][1] =  2.93610;
  C60[ 48][2] = -1.74180;
  C60[ 49][0] = -0.37880;
  C60[ 49][1] =  3.35610;
  C60[ 49][2] = -0.75130;
  C60[ 50][0] = -1.69430;
  C60[ 50][1] =  2.91860;
  C60[ 50][2] = -0.74910;
  C60[ 51][0] =  0.05210;
  C60[ 51][1] =  2.07300;
  C60[ 51][2] = -2.73550;
  C60[ 52][0] =  2.09760;
  C60[ 52][1] =  0.83400;
  C60[ 52][2] = -2.60510;
  C60[ 53][0] =  2.55170;
  C60[ 53][1] =  1.69230;
  C60[ 53][2] = -1.61070;
  C60[ 54][0] =  1.75890;
  C60[ 54][1] =  2.74520;
  C60[ 54][2] = -1.18240;
  C60[ 55][0] =  0.84200;
  C60[ 55][1] =  1.02060;
  C60[ 55][2] = -3.17860;
  C60[ 56][0] =  0.44610;
  C60[ 56][1] = -1.34950;
  C60[ 56][2] = -3.16610;
  C60[ 57][0] =  1.69830;
  C60[ 57][1] = -1.54850;
  C60[ 57][2] = -2.59080;
  C60[ 58][0] =  2.51840;
  C60[ 58][1] = -0.46230;
  C60[ 58][2] = -2.31710;
  C60[ 59][0] =  0.02180;
  C60[ 59][1] = -0.06450;
  C60[ 59][2] = -3.45850; //}}}

  // distnce from C60 center //{{{
  double C60_radius = sqrt(C60[0][0]*C60[0][0]
                         + C60[0][1]*C60[0][1]
                         + C60[0][2]*C60[0][2]); //}}}
  //}}}

  // edge distance //{{{
  double edge = sqrt((C60[0][0] - C60[2][0])*(C60[0][0] - C60[2][0]) +
                     (C60[0][1] - C60[2][1])*(C60[0][1] - C60[2][1]) +
                     (C60[0][2] - C60[2][2])*(C60[0][2] - C60[2][2])); //}}}

  // coordinates of every edge middle //{{{
  int edge_centers = 0;
  for (int i = 0; i < 60; i++) {
    for (int j = (i+1); j < 60; j++) {
      double dist = sqrt((C60[i][0] - C60[j][0])*(C60[i][0] - C60[j][0]) +
                         (C60[i][1] - C60[j][1])*(C60[i][1] - C60[j][1]) +
                         (C60[i][2] - C60[j][2])*(C60[i][2] - C60[j][2]));

      if ((dist-edge) < 0.01 && (dist-edge) > -0.01) {
        for (int k = 0; k < 3; k++) {
          E[edge_centers][k] = (C60[i][k] + C60[j][k]) / 2;
        }

        edge_centers++;
      }
    }
  } //}}}
  printf("Edge centers: %d\n", edge_centers);

  // coordinates of every hexagonal face //{{{
  int hex_centers = 0;
  double hex = 2.779;
  for (int i = 0; i < 60; i++) {
    for (int j = (i+1); j < 60; j++) {
      double dist = sqrt((C60[i][0] - C60[j][0])*(C60[i][0] - C60[j][0]) +
                         (C60[i][1] - C60[j][1])*(C60[i][1] - C60[j][1]) +
                         (C60[i][2] - C60[j][2])*(C60[i][2] - C60[j][2]));

      if ((dist-hex) < 0.05 && (dist-hex) > -0.05) {
        for (int k = 0; k < 3; k++) {
          H_temp[hex_centers][k] = (C60[i][k] + C60[j][k]) / 2;
        }

        dist = sqrt(H_temp[hex_centers][0]*H_temp[hex_centers][0] +
                    H_temp[hex_centers][1]*H_temp[hex_centers][1] +
                    H_temp[hex_centers][2]*H_temp[hex_centers][2]);

        // make distance to center the same as for C60 vertices
        for (int k = 0; k < 3; k++) {
          H_temp[hex_centers][k] *= C60_radius / dist;
        }

        hex_centers++;
      }
    }
  }
  printf("Hexagon centers: %d\n", hex_centers);

  // reduce the 60 hex center to proper 20 hex centers //{{{
  hex_centers = 0;
  for (int i = 0; i < 60; i++) {
    for (int j = 0; j < i; j++) {
      double dist = sqrt((H_temp[i][0] - H_temp[j][0])*(H_temp[i][0] - H_temp[j][0]) +
                         (H_temp[i][1] - H_temp[j][1])*(H_temp[i][1] - H_temp[j][1]) +
                         (H_temp[i][2] - H_temp[j][2])*(H_temp[i][2] - H_temp[j][2]));

      if (dist < 1) {
        for (int k = 0; k < j; k++) {
          double dist = sqrt((H_temp[i][0] - H_temp[k][0])*(H_temp[i][0] - H_temp[k][0]) +
                             (H_temp[i][1] - H_temp[k][1])*(H_temp[i][1] - H_temp[k][1]) +
                             (H_temp[i][2] - H_temp[k][2])*(H_temp[i][2] - H_temp[k][2]));

          if (dist < 1) {
            for (int l = 0; l < 3; l++) {
              H[hex_centers][l] = (H_temp[i][l] +
                                   H_temp[j][l] +
                                   H_temp[k][l]) / 3;
            }

            hex_centers++;
          }
        }
      }
    }
  } //}}}
  //}}}
  printf("Hexagon centers: %d\n", hex_centers);

  // coordinates of every pentagonal face //{{{
  int pent_centers = 0;
  double pent = 2.141;
  for (int i = 0; i < 60; i++) {
    for (int j = 0; j < 90; j++) {
      double dist = sqrt((C60[i][0] - E[j][0])*(C60[i][0] - E[j][0]) +
                         (C60[i][1] - E[j][1])*(C60[i][1] - E[j][1]) +
                         (C60[i][2] - E[j][2])*(C60[i][2] - E[j][2]));

      if ((dist-pent) < 0.05 && (dist-pent) > -0.05) {
        for (int k = 0; k < 3; k++) {
          P_temp[pent_centers][k] = (C60[i][k] + E[j][k]) / 2;
        }

        dist = sqrt(P_temp[pent_centers][0]*P_temp[pent_centers][0] +
                    P_temp[pent_centers][1]*P_temp[pent_centers][1] +
                    P_temp[pent_centers][2]*P_temp[pent_centers][2]);

        // make distance to center the same as for C60 vertices
        for (int k = 0; k < 3; k++) {
          P_temp[pent_centers][k] *= C60_radius / dist;
        }

        pent_centers++;
      }
    }
  }
  printf("Pentagon centers: %d\n", pent_centers);

  // reduce the 60 pent center to proper 12 pent centers //{{{
  pent_centers = 0;
  for (int i = 0; i < 60; i++) {
    for (int j = 0; j < i; j++) {
      double dist = sqrt((P_temp[i][0] - P_temp[j][0])*(P_temp[i][0] - P_temp[j][0]) +
                         (P_temp[i][1] - P_temp[j][1])*(P_temp[i][1] - P_temp[j][1]) +
                         (P_temp[i][2] - P_temp[j][2])*(P_temp[i][2] - P_temp[j][2]));

      if (dist < 1) {
        for (int k = 0; k < j; k++) {
          double dist = sqrt((P_temp[i][0] - P_temp[k][0])*(P_temp[i][0] - P_temp[k][0]) +
                             (P_temp[i][1] - P_temp[k][1])*(P_temp[i][1] - P_temp[k][1]) +
                             (P_temp[i][2] - P_temp[k][2])*(P_temp[i][2] - P_temp[k][2]));

          if (dist < 1) {
            for (int l = 0; l < k; l++) {
              double dist = sqrt((P_temp[i][0] - P_temp[l][0])*(P_temp[i][0] - P_temp[l][0]) +
                                 (P_temp[i][1] - P_temp[l][1])*(P_temp[i][1] - P_temp[l][1]) +
                                 (P_temp[i][2] - P_temp[l][2])*(P_temp[i][2] - P_temp[l][2]));

              if (dist < 1) {
                for (int m = 0; m < l; m++) {
                  double dist = sqrt((P_temp[i][0] - P_temp[m][0])*(P_temp[i][0] - P_temp[m][0]) +
                                     (P_temp[i][1] - P_temp[m][1])*(P_temp[i][1] - P_temp[m][1]) +
                                     (P_temp[i][2] - P_temp[m][2])*(P_temp[i][2] - P_temp[m][2]));
                  if (dist < 1) {
                    for (int n = 0; n < 3; n++) {
                      P[pent_centers][n] = (P_temp[i][n] +
                                            P_temp[j][n] +
                                            P_temp[k][n] +
                                            P_temp[l][n] +
                                            P_temp[m][n]) / 5;
                    }

                    pent_centers++;
                  }
                }
              }
            }
          }
        }
      }
    }
  } //}}}
  //}}}
  printf("Pentagon centers: %d\n", pent_centers);

  // Shell coordinates //{{{
  int count = 0;
  for (int i = 0; i < 60; i++) {
    double dist = sqrt(C60[count][0]*C60[count][0] +
                       C60[count][1]*C60[count][1] +
                       C60[count][2]*C60[count][2]);
    for (int j = 0; j < 3; j++) {
      Shell[count][j] = C60[i][j] / dist * atof(argv[1]);
    }
    count++;
  }
  for (int i = 0; i < 20; i++) {
    double dist = sqrt(H[count][0]*H[count][0] +
                       H[count][1]*H[count][1] +
                       H[count][2]*H[count][2]);
    for (int j = 0; j < 3; j++) {
      Shell[count][j] = H[i][j] / dist * atof(argv[1]);
    }
    count++;
  }
  for (int i = 0; i < 12; i++) {
    double dist = sqrt(P[count][0]*P[count][0] +
                       P[count][1]*P[count][1] +
                       P[count][2]*P[count][2]);
    for (int j = 0; j < 3; j++) {
      Shell[count][j] = P[i][j] / dist * atof(argv[1]);
    }
    count++;
  } //}}}

  // Inner C60 //{{{
//// distnce from C60 center //{{{
//C60_radius = sqrt(C60[0][0]*C60[0][0]
//                + C60[0][1]*C60[0][1]
//                + C60[0][2]*C60[0][2]); //}}}

  // coefficient for multiplying of coordinates
  double RequiredRadiusCoefficient = atof(argv[1]) * 1 / (2 * C60_radius);

  for (int i = 0; i < 60; i++) {
    for (int j = 0; j < 3; j++) {
      C60[i][j] *= RequiredRadiusCoefficient;
    }
  } //}}}

  FILE *fw1 = fopen("C60.vtf", "w");
  FILE *fw2 = fopen("Nano-FIELD", "w");

  // write coordinates to FIELD //{{{
  fprintf(fw2, "beads %d\n", beads);
  // center bead
  fprintf(fw2, "Inner   0.0        0.0        0.0\n");

//// octahedron //{{{
//for (int i = 0; i < 6; i++) {
//  fprintf(fw2, "Inner");
//  for (int j = 0; j < 3; j++) {
//    fprintf(fw2, " %10f", Octahedron[i][j]);
//  }
//  putc('\n', fw2);
//} //}}}

//// icosahedron //{{{
//for (int i = 0; i < 12; i++) {
//  fprintf(fw2, "Inner");
//  for (int j = 0; j < 3; j++) {
//    fprintf(fw2, " %10f", Icosahedron[i][j]);
//  }
//  putc('\n', fw2);
//} //}}}

  // inner C60 //{{{
  for (int i = 0; i < 60; i++) {
    fprintf(fw2, "Inner");
    for (int j = 0; j < 3; j++) {
      fprintf(fw2, " %10f", C60[i][j]);
    }
    putc('\n', fw2);
  } //}}}

  // Shell //{{{
  for (int i = 0; i < 92; i++) {
    fprintf(fw2, "Shell");
    for (int j = 0; j < 3; j++) {
      fprintf(fw2, " %10f", Shell[i][j]);
    }
    putc('\n', fw2);
  } //}}}
  //}}}

  fprintf(fw2, "bond 11628\n");

  // vtf structure - beads //{{{

  // center bead
  count = 0;
  fprintf(fw1, "atom %d n Inner m 1 q 0\n", count);

//// inner beads - octahedron
//for (int i = 0; i < 6; i++) {
//  fprintf(fw1, "atom %d n Inner\n", ++count);
//}

//// inner beads - icosahedron
//for (int i = 0; i < 12; i++) {
//  fprintf(fw1, "atom %d n Inner\n", ++count);
//}

  // inner beads - C60
  for (int i = 0; i < 60; i++) {
    fprintf(fw1, "atom %d n Inner\n", ++count);
  }

  // shell beads - C60
  fprintf(fw1, "atom %d n Shell m 1 q 0\n", ++count);
  for (int i = 1; i < 92; i++) {
    fprintf(fw1, "atom %d n Shell\n", ++count);
  }
  //}}}

  // bonds between all beads //{{{
  double All[beads][3]; //{{{
  count = 0;
  All[count][0] = 0;
  All[count][1] = 0;
  All[count][2] = 0;
  for (int i = 0; i < 60; i++) {
    count++;
    for (int j = 0; j < 3; j++) {
      All[count][j] = C60[i][j];
    }
  }
  for (int i = 0; i < 92; i++) {
    count++;
    for (int j = 0; j < 3; j++) {
      All[count][j] = Shell[i][j];
    }
  } //}}}

  for (int i = 1; i < beads; i++) {
    for (int j = 0; j < i; j++) {
      double dist = sqrt((All[i][0] - All[j][0])*(All[i][0] - All[j][0]) +
                         (All[i][1] - All[j][1])*(All[i][1] - All[j][1]) +
                         (All[i][2] - All[j][2])*(All[i][2] - All[j][2]));

      fprintf(fw1, "bond %d:%d\n", j, i);
      fprintf(fw2, "harm %3d %3d %lf %2.2f\n", j+1, i+1, bond, dist);
    }
  } //}}}

  // vtf coordinates //{{{
  fprintf(fw1, "\ntimestep ordered\n");

  // center bead
  fprintf(fw1, "   0.0        0.0        0.0\n");

//// octahedron //{{{
//for (int i = 0; i < 6; i++) {
//  fprintf(fw1, "%lf %lf %lf\n", Octahedron[i][0], Octahedron[i][1], Octahedron[i][2]);
//}

//// icosahedron
//for (int i = 0; i < 12; i++) {
//  fprintf(fw1, "%lf %lf %lf\n", Icosahedron[i][0], Icosahedron[i][1], Icosahedron[i][2]);
//} //}}}

  // inner C60
  for (int i = 0; i < 60; i++) {
    fprintf(fw1, "%lf %lf %lf\n", C60[i][0], C60[i][1], C60[i][2]);
  }

  // shell
  for (int i = 0; i < 92; i++) {
    fprintf(fw1, "%lf %lf %lf\n", Shell[i][0], Shell[i][1], Shell[i][2]);
  }
  //}}}

  fclose(fw1);
  fclose(fw2);

  return 0;
}
