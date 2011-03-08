/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 03/04/2011
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(backbone,FixBackbone)

#else

#ifndef LMP_FIX_BACKBONE_H
#define LMP_FIX_BACKBONE_H

#include "fix.h"

#include "smart_matrix_lib.h"
#include "fragment_memory.h"

namespace LAMMPS_NS {

class FixBackbone : public Fix {
 public:
  FixBackbone(class LAMMPS *, int, char **);
  ~FixBackbone();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  void init_list(int, class NeighList *);

// private:
public:
  double epsilon;
  double k_chain[3], k_shake, k_chi, k_rama;
  double k_excluded_C, k_excluded_O;
  double r_ncb0, r_cpcb0, r_ncp0, chi0;
  double r_sh1, r_sh2, r_sh3;
  double rC_ex0, rO_ex0;
  int p; // Excluded volume: (r-r0)^p
  int n_rama_par, n_rama_p_par;
  static const int i_rp = 6; // Where Proline rama parametors start
  bool ssweight[12];
  double w[12], sigma[12], phiw[12], phi0[12], psiw[12], psi0[12], *aps[12];
  double hbscl[4][9], sigma_NO, sigma_HO, NO_zero, HO_zero;
  double dssp_hdrgn_cut, pref[2], d_nu0;
  double k_P_AP[3], P_AP_pref, P_AP_cut;
  int i_diff_P_AP, i_med_max, i_med_min;
  double m_anti_HB[20][20][2], m_anti_NHB[20][20][2], m_para_HB[20][20][2];
  double m_para_one[20], m_anti_one[20];
  double k_water, k_burial;
  double well_r_min[5], well_r_max[5], treshold, water_kappa, water_kappa_sigma, burial_kappa;
  double water_gamma[5][20][20][2], burial_gamma[20][3];
  double burial_ro_min[3], burial_ro_max[3];
  double k_helix, helix_gamma_p, helix_gamma_w, h4prob[20];
  double helix_kappa, helix_kappa_sigma, helix_treshold, helix_cutoff;
  int well_flag[5], n_wells, contact_cutoff, helix_i_diff;
  double k_amh_go, amh_go_rc;
  int amh_go_p;
  Gamma_Array *amh_go_gamma;
  double **amh_go_force;
  int *amh_go_force_map;
  double k_frag_mem[2];
  Fragment_Memory *m_amh_go;
  int i_fm_min, i_fm_med_min, i_fm_med_max;
  char fmem_file[100];
  int igroup2, group2bit;
  int igroup3, group3bit;
  double foriginal[4],foriginal_all[4];
  int force_flag;
  int nlevels_respa;
  bool allocated;
  class NeighList *list;         // standard neighbor list used by most pairs
  
  int ntimestep;
  int n, nn;
  double an, bn, cn, ap, bp, cp, ah, bh, ch;
  int *alpha_carbons;
  int *beta_atoms;
  int *oxygens;
  int *res_no, *res_info;
  double **xca, **xcb, **xo, **xn, **xcp, **xh;
  double **x, **f;
  int *image;
  double prd[3], half_prd[3];
  int *periodicity;
  bool abc_flag, chain_flag, shake_flag, chi_flag, rama_flag, rama_p_flag, excluded_flag, p_excluded_flag, r6_excluded_flag;
  bool ssweight_flag, dssp_hdrgn_flag, p_ap_flag, water_flag, burial_flag, helix_flag, amh_go_flag, frag_mem_flag;
  
  enum Atoms{CA0 = 0, CA1, CA2, O0, O1, nAtoms};
  enum Angles{PHI = 0, PSI, nAngles};
  enum ResInfo{NONE=0, LOCAL, GHOST, OFF};
  
  char se[1000]; // Protein sequance
  
 private:
  void compute_backbond();
  void compute_chain_potential(int i);
  void compute_shake(int i);
  void compute_chi_potential(int i);
  void compute_rama_potential(int i);
  void compute_excluded_volume();
  void compute_p_degree_excluded_volume();
  void compute_r6_excluded_volume();
  void compute_dssp_hdrgn(int i, int j);
  void compute_P_AP_potential(int i, int j);
  void compute_water_potential(int i, int j);
  void compute_burial_potential(int i);
  void compute_helix_potential(int i, int j);
  void compute_amh_go_model();
  void compute_fragment_memory_potential(int i);

  void allocate();
  inline void Construct_Computational_Arrays();
  int Tag(int index);

  void calcDihedralAndSlopes(int, double& angle, int iAng);
  double y_slope[nAngles][nAtoms][3], x_slope[nAngles][nAtoms][3];
  inline double PeriodicityCorrection(double d, int i);
  inline bool isFirst(int index);
  inline bool isLast(int index);
  inline double anti_HB(int res1, int res2, int k);
  inline double anti_NHB(int res1, int res2, int k);
  inline double para_HB(int res1, int res2, int k);
  inline double para_one(int res);
  inline double anti_one(int res);

  cP_AP<double, FixBackbone> *p_ap;
  cR<double, FixBackbone> *R;
  cWell<double, FixBackbone> *well;
  cWell<double, FixBackbone> *helix_well;
  
  WPV water_par;
  WPV helix_par;

  void out_xyz_and_force(int coord=0);
  void initforce(double **tmpf);
  void copyforce(double **tmpf);
  double findBiggestForce(double **tmpf, int iatom);
  void printforce(double **tmpf);
  void sumForces(double *fsum);
  double maxsumf;
  int stepmaxsumf;
  int sStep, eStep;
  int Step;
  FILE *fout;
};

}

#endif
#endif
