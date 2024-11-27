#include "../AnalysisTools.h"
#include <stdio.h>

// name for the monovalent cation
static char *name = "C";
// name for the divalent cation (molecule with two connected beads)
static char *name_mol = "CA2";
// calculate distance between two points, accounting for pbc //{{{
inline static double DistLength(const double v1[3], const double v2[3],
                                const double box[3]) {
  double dist[3];
  Distance(v1, v2, box, dist);
  return VectLength(dist);
} //}}}
// get two molecule/bead types ordered so the first one < second one //{{{
typedef struct {
  int a, b;
} Types;
inline static Types SortTypes(const int type1, const int type2) {
  Types tp;
  tp.a = type1;
  tp.b = type2;
  if (tp.a > tp.b) {
    SwapInt(&tp.a, &tp.b);
  }
  return tp;
} //}}}

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Count contacts intra- and intermolecular contacts between specified bead types \
in specified molecule types. Also counts 3-body contacts, specifically, \
two of the specified bead types with 'C' bead type or 'CA2' molecule type's \
geometric centre (these names are hardcoded).\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <output> <double> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output file\n");
  fprintf(ptr, "<skip>              number of in-between beads to skip\n");
  fprintf(ptr, "<dist>              minimum distance contact check\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -mt <name(s)>     use specified molecule type(s)\n");
  fprintf(ptr, "  -bt <name(s)>     use specified bead type(s)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  // here com option variables
  bool bool_opt;     // --bool
  bool *mt;          // -mt
  bool *bt;          // -bt
  int int1;          // -int1
  int int2[2];       // -int2
  char f_file[LINE]; // -f (filename)
  int f_list[100],   // -f (list of numbers)
      f_num;         // -f (number of those numbers)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity //{{{
  int common = 8, all = common + 2, count = 0,
      req_arg = 4;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version", "-mt", "-bt"); //}}}

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // mandatory options //{{{
  // <input> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <output> - output file name
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // <skip> - how many beads to skip at least betweem contact-able beads
  long skip = 0;
  if (!IsWholeNumber(argv[++count], &skip)) {
    ErrorNaN("<skip>");
    Help(argv[0], true, common, option);
    exit(1);
  }
  double dist_check = 0;
  if (!IsPosRealNumber(argv[++count], &dist_check)) {
    ErrorNaN("<dist>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  double *boxlength = System.Box.Length;

  // molecule/bead type options //{{{
  // molecule types to calculate contacts for
  opt->mt = calloc(System.Count.MoleculeType, sizeof *opt->mt);
  if (!MoleculeTypeOption(argc, argv, "-mt", true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }
  // bead types to calculate contacts for
  opt->bt = calloc(System.Count.BeadType, sizeof *opt->bt);
  if (!BeadTypeOption(argc, argv, "-bt", true, opt->bt, System)) {
    InitBoolArray(opt->bt, Count->BeadType, true);
  } //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // arrays for all the necessary stuff //{{{
  // count molecules of each type
  int *c_mtype = calloc(Count->MoleculeType, sizeof *c_mtype);
  // count per moltype intramolecular contacts
  long int ***intra_mol = calloc(Count->MoleculeType, sizeof *intra_mol);
  long int ****intra_3body = calloc(Count->MoleculeType, sizeof *intra_3body);
  for (int i = 0; i < Count->MoleculeType; i++) {
    intra_mol[i] = calloc(Count->BeadType, sizeof *intra_mol[i]);
    intra_3body[i] = calloc(Count->BeadType, sizeof *intra_3body[i]);
    for (int j = 0; j < Count->BeadType; j++) {
      intra_mol[i][j] = calloc(Count->BeadType, sizeof *intra_mol[i][j]);
      intra_3body[i][j] = calloc(Count->BeadType, sizeof *intra_3body[i][j]);
      for (int k = 0; k < Count->BeadType; k++) {
        intra_3body[i][j][k] = calloc(Count->BeadType,
                                      sizeof *intra_3body[i][j][k]);
      }
    }
  }
  // count molecules of each type
  int **c_mtype_mtype = calloc(Count->MoleculeType, sizeof *c_mtype_mtype);
  // count per moltype-moltype pair intermolecular contacts
  long int ****inter_mol = calloc(Count->MoleculeType, sizeof *inter_mol);
  long int *****inter_3body = calloc(Count->MoleculeType, sizeof *inter_3body);
  for (int i = 0; i < Count->MoleculeType; i++) {
    c_mtype_mtype[i] = calloc(Count->MoleculeType, sizeof *c_mtype_mtype[i]);
    inter_mol[i] = calloc(Count->MoleculeType, sizeof *inter_mol[i]);
    inter_3body[i] = calloc(Count->MoleculeType, sizeof *inter_3body[i]);
    for (int j = 0; j < Count->MoleculeType; j++) {
      inter_mol[i][j] = calloc(Count->BeadType, sizeof *inter_mol[i][j]);
      inter_3body[i][j] = calloc(Count->BeadType, sizeof *inter_3body[i][j]);
      for (int k = 0; k < Count->BeadType; k++) {
        inter_mol[i][j][k] = calloc(Count->BeadType,
                                    sizeof *inter_mol[i][j][k]);
        inter_3body[i][j][k] = calloc(Count->BeadType,
                                      sizeof *inter_3body[i][j][k]);
        for (int l = 0; l < Count->BeadType; l++) {
          inter_3body[i][j][k][l] = calloc(Count->BeadType,
                                           sizeof *inter_3body[i][j][k][l]);
        }
      }
    }
  } //}}}
  int **contacts = calloc(Count->Bonded * Count->Bonded, sizeof *contacts);
  for (int i = 0; i < (Count->Bonded * Count->Bonded); i++) {
    contacts[i] = calloc(2, sizeof *contacts[i]);
  }

  PrintByline(fout, argc, argv);
  FILE *fw = OpenFile(fout, "a");
  count = 1;
  fprintf(fw, "# (%d) step\n", count++);
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      fprintf(fw, "# molecule %s: ", System.MoleculeType[i].Name);
      for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
        for (int k = j; k < System.MoleculeType[i].nBTypes; k++) {
          int btj = System.MoleculeType[i].BType[j];
          int btk = System.MoleculeType[i].BType[k];
          if (opt->bt[btj] && opt->bt[btk]) {
            fprintf(fw, " (%d) %s-%s;", count++,
                    System.BeadType[btj].Name, System.BeadType[btk].Name);
            if (FindBeadType(name, System) != -1) {
              fprintf(fw, " (%d) %s-%s-%s;", count++, System.BeadType[btj].Name,
                      System.BeadType[btk].Name, name);
            }
            if (FindMoleculeName(name_mol, System) != -1) {
              fprintf(fw, " (%d) %s-%s-%s", count++, System.BeadType[btj].Name,
                      System.BeadType[btk].Name, name_mol);
            }
          }
        }
      }
      putc('\n', fw);
    }
  }

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      // go over all molecules in the coordinate file //{{{
      int count_contacts = 0;
      for (int m1 = 0; m1 < Count->MoleculeCoor; m1++) {
        int m1_id = System.MoleculeCoor[m1];
        MOLECULE *mol1 = &System.Molecule[m1_id];
        int m1_type = mol1->Type;
        MOLECULETYPE *mt1 = &System.MoleculeType[m1_type];
        if (opt->mt[m1_type]) { // use only specified molecule types
          c_mtype[m1_type]++;
          // go over all beads in the molecule
          for (int b1 = 0; b1 < mt1->nBeads; b1++) {
            int b1_id = mol1->Bead[b1];
            BEAD *bead1 = &System.Bead[b1_id];
            if (opt->bt[bead1->Type]) { // use only specified bead types
              // 1) intramolecular contacts
              // go over all other beads that are far enough from the first bead
              for (int b2 = (b1 + 1 + skip); b2 < mt1->nBeads; b2++) {
                int b2_id = mol1->Bead[b2];
                BEAD *bead2 = &System.Bead[b2_id];
                if (opt->bt[bead2->Type]) { // use only specified bead types
                  double d = DistLength(bead1->Position,
                                    bead2->Position, boxlength);
                  // are the two beads are close enough?
                  if (d < dist_check) {
                    contacts[count_contacts][0] = b1_id;
                    contacts[count_contacts][1] = b2_id;
                    count_contacts++;
                  }
                }
              }
              // 2) intermolecular contacts
              // go over other molecules
              for (int m2 = (m1 + 1); m2 < Count->MoleculeCoor; m2++) {
                int m2_id = System.MoleculeCoor[m2];
                MOLECULE *mol2 = &System.Molecule[m2_id];
                int m2_type = mol2->Type;
                MOLECULETYPE *mt2 = &System.MoleculeType[m2_type];
                if (opt->mt[m2_type]) { // use only specified molecule types
                  Types mt = SortTypes(m1_type, m2_type);
                  // printf("YYY %d %d .. %d %d\n", m1_type, mt[0], m2_type, mt[1]);
                  c_mtype_mtype[mt.a][mt.b]++;
                  for (int b2 = 0; b2 < mt2->nBeads; b2++) {
                    int b2_id = mol2->Bead[b2];
                    BEAD *bead2 = &System.Bead[b2_id];
                    if (opt->bt[bead2->Type]) { // use only specified bead types
                      double d = DistLength(bead1->Position, bead2->Position,
                                            boxlength);
                      // count contact if the two beads are close enough
                      if (d < dist_check) {
                        contacts[count_contacts][0] = b1_id;
                        contacts[count_contacts][1] = b2_id;
                        count_contacts++;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      int ****per_step = calloc(Count->MoleculeType, sizeof *per_step);
      for (int i = 0; i < Count->MoleculeType; i++) {
        per_step[i] = calloc(Count->BeadType, sizeof *per_step[i]);
        for (int j = 0; j < Count->BeadType; j++) {
          per_step[i][j] = calloc(Count->BeadType, sizeof *per_step[i][j]);
          for (int k = 0; k < Count->BeadType; k++) {
            per_step[i][j][k] = calloc(3, sizeof *per_step[i][j][k]);
          }
        }
      }
      bool *three_body = calloc(count_contacts, sizeof *three_body);
      char tcl[LINE] = "";
      FILE *out_vmd = NULL;
      if (count_contacts > 0) {
        snprintf(tcl, LINE, "contacts-%04d.tcl", count_used);
        out_vmd = OpenFile(tcl, "w");
      }
      count_used++;
      // i) is monovalent counterion near (globally defined name)?
      int bt3 = FindBeadType(name, System);
      if (bt3 != -1) {
        for (int i = 0; i < System.BeadType[bt3].Number; i++) {
          int b3_id = System.BeadType[bt3].Index[i];
          BEAD *b3 = &System.Bead[b3_id];
          for (int j = 0; j < count_contacts; j++) {
            int b1_id = contacts[j][0];
            int b2_id = contacts[j][1];
            BEAD *b1 = &System.Bead[b1_id];
            BEAD *b2 = &System.Bead[b2_id];
            Types bt = SortTypes(b1->Type, b2->Type);
            int m1_type = System.Molecule[b1->Molecule].Type;
            int m2_type = System.Molecule[b2->Molecule].Type;
            Types mt = SortTypes(m1_type, m2_type);
            double d[2];
            d[0] = DistLength(b1->Position, b3->Position, boxlength);
            d[1] = DistLength(b2->Position, b3->Position, boxlength);
            if (d[0] < dist_check && d[1] < dist_check) {
              if (mt.a == mt.b) {
                intra_3body[mt.a][bt.a][bt.b][bt3]++;
              } else {
                inter_3body[mt.a][mt.b][bt.a][bt.b][bt3]++;
              }
              per_step[mt.a][bt.a][bt.b][1]++;
              three_body[j] = true;
              fprintf(out_vmd, "set rep [expr $rep + 1]\n");
              fprintf(out_vmd, "mol addrep ${mol}\n");
              fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d %d\n",
                      b1_id, b2_id, b3_id);
              fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
              fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
              // break;
            }
          }
        }
      }
      int m3_type = FindMoleculeName(name_mol, System);
      if (m3_type != -1) {
        MOLECULETYPE *mt3 = &System.MoleculeType[m3_type];
        for (int i = 0; i < mt3->Number; i++) {
          MOLECULE *mol = &System.Molecule[mt3->Index[i]];
          double gc[3];
          GeomCentre(mt3->nBeads, mol->Bead, System.Bead, gc);
          int btype = mt3->BType[0];
          for (int j = 0; j < count_contacts; j++) {
            int b1_id = contacts[j][0];
            int b2_id = contacts[j][1];
            BEAD *b1 = &System.Bead[b1_id];
            BEAD *b2 = &System.Bead[b2_id];
            Types bt = SortTypes(b1->Type, b2->Type);
            int m1_type = System.Molecule[b1->Molecule].Type;
            int m2_type = System.Molecule[b2->Molecule].Type;
            Types mt = SortTypes(m1_type, m2_type);
            double d[2];
            d[0] = DistLength(b1->Position, gc, boxlength);
            d[1] = DistLength(b2->Position, gc, boxlength);
            if (d[0] < dist_check && d[1] < dist_check) {
              if (mt.a == mt.b) {
                intra_3body[mt.a][bt.a][bt.b][btype]++;
              } else {
                inter_3body[mt.a][mt.b][bt.a][bt.b][btype]++;
              }
              per_step[mt.a][bt.a][bt.b][2]++;
              three_body[j] = true;
              fprintf(out_vmd, "set rep [expr $rep + 1]\n");
              fprintf(out_vmd, "mol addrep ${mol}\n");
              fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d ",
                      b1_id, b2_id);
              fprintf(out_vmd, "or resid %d\n", mol->Index);
              fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
              fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
              // break;
            }
          }
        }
      }
      for (int i = 0; i < count_contacts; i++) {
        if (!three_body[i]) {
          int b1_id = contacts[i][0];
          int b2_id = contacts[i][1];
          BEAD *b1 = &System.Bead[b1_id];
          BEAD *b2 = &System.Bead[b2_id];
          Types bt = SortTypes(b1->Type, b2->Type);
          int m1_type = System.Molecule[b1->Molecule].Type;
          int m2_type = System.Molecule[b2->Molecule].Type;
          Types mt = SortTypes(m1_type, m2_type);
          per_step[mt.a][bt.a][bt.b][0]++;
          if (mt.a == mt.b) {
            intra_mol[mt.a][bt.a][bt.b]++;
          } else {
            inter_mol[mt.a][mt.b][bt.a][bt.b]++;
          }
          fprintf(out_vmd, "set rep [expr $rep + 1]\n");
          fprintf(out_vmd, "mol addrep ${mol}\n");
          fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d\n",
                  b1_id, b2_id);
          fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
          fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
        }
      }

      fprintf(fw, "%5d", count_used);
      for (int i = 0; i < Count->MoleculeType; i++) {
        if (opt->mt[i]) {
          for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
            for (int k = j; k < System.MoleculeType[i].nBTypes; k++) {
              int btj = System.MoleculeType[i].BType[j];
              int btk = System.MoleculeType[i].BType[k];
              if (opt->bt[btj] && opt->bt[btk]) {
                double avg = (double)(per_step[i][btj][btk][0]) /
                             System.MoleculeType[i].Number;
                fprintf(fw, " %lf", avg);
                avg = (double)(per_step[i][btj][btk][1]) /
                      System.MoleculeType[i].Number;
                fprintf(fw, " %lf", avg);
                avg = (double)(per_step[i][btj][btk][2]) /
                      System.MoleculeType[i].Number;
                fprintf(fw, " %lf", avg);
              }
            }
          }
        }
      }
      putc('\n', fw);

      free(three_body);
      for (int i = 0; i < Count->MoleculeType; i++) {
        for (int j = 0; j < Count->BeadType; j++) {
          for (int k = 0; k < Count->BeadType; k++) {
            free(per_step[i][j][k]);
          }
          free(per_step[i][j]);
        }
        free(per_step[i]);
      }
      free(per_step);
      //}}}
      if (count_contacts > 0) {
        fclose(out_vmd);
      }
      //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // print results //{{{
  for (int m1 = 0; m1 < Count->MoleculeType; m1++) {
    MOLECULETYPE *mt1 = &System.MoleculeType[m1];
    if (opt->mt[m1]) { // use only specified molecule types
      for (int b1 = 0; b1 < mt1->nBTypes; b1++) {
        int bt1 = mt1->BType[b1];
        if (opt->bt[bt1]) { // only specified bead types
          // 1) intramolecular contacts
          for (int b2 = b1; b2 < mt1->nBTypes; b2++) {
            int bt2 = mt1->BType[b2];
            if (opt->bt[bt2] && c_mtype[m1] > 0) { // only specified bead types
              double avg = (double)(intra_mol[m1][bt1][bt2]) / c_mtype[m1];
              fprintf(fw, "# mol %s, %s-%s contacts: %lf\n", mt1->Name,
                      System.BeadType[bt1].Name, System.BeadType[bt2].Name,
                      avg);
              int bt3 = FindBeadType(name, System);
              if (bt3 != -1) {
                avg = (double)(intra_3body[m1][bt1][bt2][bt3]) / c_mtype[m1];
                fprintf(fw, "#   3body (%s): %lf\n",
                        System.BeadType[bt3].Name, avg);
              }
              bt3 = FindMoleculeName(name_mol, System);
              if (bt3 != -1) {
                bt3 = System.MoleculeType[bt3].BType[0];
                avg = (double)(intra_3body[m1][bt1][bt2][bt3]) / c_mtype[m1];
                fprintf(fw, "#   3body (%s): %lf\n",
                        System.BeadType[bt3].Name, avg);
              }
            }
          }
          // 2) intermolecular contacts
          // go over other molecules
          for (int m2 = m1; m2 < Count->MoleculeType; m2++) {
            MOLECULETYPE *mt2 = &System.MoleculeType[m2];
            if (opt->mt[m2]) { // use only specified molecule types
              for (int b2 = 0; b2 < mt2->nBTypes; b2++) {
                int bt2 = mt2->BType[b2];
                if (opt->bt[bt2] && c_mtype_mtype[m1][m2] > 0) {
                  double avg = (double)(inter_mol[m1][m2][bt1][bt2]) /
                               c_mtype_mtype[m1][m2];
                  fprintf(fw, "# mol %s-%s, %s-%s contacts: %lf\n",
                          mt1->Name, mt2->Name,
                          System.BeadType[bt1].Name,
                          System.BeadType[bt2].Name, avg);
                  int bt3 = FindBeadType(name, System);
                  if (bt3 != -1) {
                    avg = (double)(inter_3body[m1][m2][bt1][bt2][bt3]) /
                          c_mtype_mtype[m1][m2];
                    fprintf(fw, "#   3body (%s): %lf\n",
                            System.BeadType[bt3].Name, avg);
                  }
                  if (bt3 != -1) {
                    bt3 = FindMoleculeName(name_mol, System);
                    bt3 = System.MoleculeType[bt3].BType[0];
                    avg = (double)(inter_3body[m1][m2][bt1][bt2][bt3]) /
                          c_mtype_mtype[m1][m2];
                    fprintf(fw, "#   3body (%s): %lf\n",
                            System.BeadType[bt3].Name, avg);
                  }
                }
              }
            }
          }
        }
      }
    }
  } //}}}
  fclose(fw);

  // free memory - to make valgrind happy //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        free(intra_3body[i][j][k]);
      }
      free(intra_mol[i][j]);
      free(intra_3body[i][j]);
    }
    free(intra_mol[i]);
    free(intra_3body[i]);
  }
  free(intra_mol);
  free(intra_3body);
  free(c_mtype);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->MoleculeType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        for (int l = 0; l < Count->BeadType; l++) {
          free(inter_3body[i][j][k][l]);
        }
        free(inter_mol[i][j][k]);
        free(inter_3body[i][j][k]);
      }
      free(inter_mol[i][j]);
      free(inter_3body[i][j]);
    }
    free(inter_mol[i]);
    free(inter_3body[i]);
    free(c_mtype_mtype[i]);
  }
  free(inter_mol);
  free(inter_3body);
  free(c_mtype_mtype);
  for (int i = 0; i < (Count->Bonded * Count->Bonded); i++) {
    free(contacts[i]);
  }
  free(contacts);
  free(opt->mt);
  free(opt->bt);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
