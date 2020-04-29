/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Solvent Separated Barrier Potential was contributed by Nick Schafer

Last Update: 03/23/2011
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
  // Standart lammaps interface
  FixBackbone(class LAMMPS *, int, char **);
  ~FixBackbone();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void pre_force_respa(int, int, int);
  void min_pre_force(int);
  double compute_scalar();
  double compute_vector(int);
  void init_list(int, class NeighList *);

// private:
public:
  // Global energy scale
  double epsilon;
  
  // Backbone parameters
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
  
  // Hydrogen bonding parameters
  double hbscl[4][9], sigma_NO, sigma_HO, NO_zero, HO_zero;
  double k_dssp, dssp_hdrgn_cut, pref[2], d_nu0;
  
  // P_AP Liquid Crystal potential parameters
  double k_global_P_AP, k_betapred_P_AP, k_P_AP[3], P_AP_pref, P_AP_cut;
  int i_diff_P_AP, i_med_max, i_med_min;
  
  // Water mediated interactions parameters
  double m_anti_HB[20][20][2], m_anti_NHB[20][20][2], m_para_HB[20][20][2];
  double m_para_one[20], m_anti_one[20];
  double k_water;
  double well_r_min[5], well_r_max[5], treshold, water_kappa, water_kappa_sigma, burial_kappa;
  double water_gamma[5][20][20][2];
  int well_flag[5], n_wells, contact_cutoff;
  
  // Phosphorylation via hypercharged glutamate interaction
  double phosph_water_gamma[5][20][20][2];
  double k_hypercharge;
  int n_phosph_res;
  double phosph_res[20];
  int *phosph_map;
  // Burial potential parameters
  double k_burial;
  double burial_ro_min[3], burial_ro_max[3];
  double burial_gamma[20][3];
  
  // Helical hydrogen bonding parameters
  double k_helix, helix_gamma_p, helix_gamma_w, h4prob[20];
  double helix_kappa, helix_kappa_sigma, helix_treshold, helix_cutoff;
  double helix_well_r_min[5], helix_well_r_max[5];
  double helix_sigma_HO, helix_sigma_NO, helix_HO_zero, helix_NO_zero;
  int helix_well_flag[5], n_helix_wells, helix_i_diff;
  int pro_accepter_flag;
  double h4prob_pro_accepter;
  
  // Non-additive AMH-Go parameters
  double k_amh_go, amh_go_rc;
  double amh_go_p;
  Fragment_Memory *m_amh_go;
  Gamma_Array *amh_go_gamma;
  double **amh_go_force;
  int *amh_go_force_map;
  double *amh_go_norm;
  
  // Fragment Memory parameters
  double k_frag_mem;
  int n_frag_mems, **frag_mem_map, *ilen_fm_map;
  Fragment_Memory **frag_mems;
  Gamma_Array *fm_gamma;
  char frag_mems_file[100];
  char fm_gamma_file[100];
  double fm_sigma_exp;
  
  // Table Fragment Memory parameters
  TBV **fm_table;
  int tb_size, tb_nbrs;
  double tb_rmin, tb_rmax, tb_dr;
  
  // Vector Fragment Memory
  double k_vec_frag_mem;
  double vfm_sigma, vfm_sigma_sq;
  double frag_table_well_width;
  int fm_energy_debug_flag;

  
  // Table Vector Fragment Memory
  TBV **vfm_table;
  int vfm_tb_size;
  double vfm_tb_vmin, vfm_tb_vmax, vfm_tb_dv;

  // Solvent separated barrier
  double k_solventb1, k_solventb2;
  double ssb_kappa, ssb_rmin1, ssb_rmax1, ssb_rmin2, ssb_rmax2;
  int ssb_ij_sep;
  bool ssb_rad_cor;
  double ssb_rshift[20];

  // Standart lammaps interface
  int igroup2, group2bit;
  int igroup3, group3bit;
  int force_flag;
  int nlevels_respa;
  bool allocated;
  class NeighList *list;         // standard neighbor list used by most pairs
  
  int ntimestep;
  int n, nn; // n is the total number of residues, nn is the local number of residues
  double an, bn, cn, ap, bp, cp, ah, bh, ch;
  int *alpha_carbons;
  int *beta_atoms;
  int *oxygens;
  int *res_no, *res_info, *chain_no;
  int *res_no_l;
  double **xca, **xcb, **xo, **xn, **xcp, **xh;
  double **x, **f;
  int *image;
  double prd[3], half_prd[3];
  int *periodicity;
  bool abc_flag, chain_flag, shake_flag, chi_flag, rama_flag, rama_p_flag, excluded_flag, p_excluded_flag, r6_excluded_flag;
  bool ssweight_flag, dssp_hdrgn_flag, p_ap_flag, water_flag, burial_flag, helix_flag, amh_go_flag, frag_mem_flag, ssb_flag;
  bool phosph_flag;
  bool frag_mem_tb_flag, vec_frag_mem_flag, vec_frag_mem_tb_flag;
  
  enum Atoms{CA0 = 0, CA1, CA2, O0, O1, nAtoms};
  enum Angles{PHI = 0, PSI, nAngles};
  enum ResInfo{NONE=0, LOCAL, GHOST, OFF};
  
  char *se; // Protein sequance
  int nch, ch_len[100], ch_pos[100];
  
  double energy[15], energy_all[15];
  enum EnergyTerms{ET_TOTAL=0, ET_CHAIN, ET_SHAKE, ET_CHI, ET_RAMA, ET_VEXCLUDED, ET_DSSP, ET_PAP, 
                    ET_WATER, ET_BURIAL, ET_HELIX, ET_AMHGO, ET_FRAGMEM, ET_VFRAGMEM, ET_SSB, nEnergyTerms};
  
  double ctime[15], previous_time;
  enum ComputeTime{TIME_CHAIN=0, TIME_SHAKE, TIME_CHI, TIME_RAMA, TIME_VEXCLUDED, TIME_DSSP, TIME_PAP, 
  					TIME_WATER, TIME_BURIAL, TIME_HELIX, TIME_AMHGO, TIME_FRAGMEM, TIME_VFRAGMEM, TIME_SSB, TIME_N};
  
 private:
  void compute_backbone();
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
  void compute_solvent_barrier(int i, int j);
  void compute_fragment_memory_table();
  void table_fragment_memory(int i, int j);
  void compute_amhgo_normalization();
  void compute_vector_fragment_memory_potential(int i);
  void compute_vector_fragment_memory_table();
  void table_vector_fragment_memory(int i, int j);

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
  inline double get_water_gamma(int i_resno, int j_resno, int i_well, int ires_type, int jres_type, int local_dens);
  inline double get_burial_gamma(int i_resno, int irestype, int local_dens);
  
  inline void print_log(char *line);
  void final_log_output();
  Fragment_Memory **read_mems(char *mems_file, int &n_mems);
  bool isEmptyString(char *str);
  char *ltrim(char *s);
  char *rtrim(char *s);
  char *trim(char *s);
  
  void timerBegin();
  void timerEnd(int which);

  cP_AP<double, FixBackbone> *p_ap;
  cR<double, FixBackbone> *R;
  cWell<double, FixBackbone> *well;
  cWell<double, FixBackbone> *helix_well;
  
  WPV water_par;
  WPV helix_par;
  
  FILE *efile;
  
  FILE *dout;
  int sStep, eStep;
  void print_forces(int coord=0);
  
/*  double tmpforce1[1000][3];
  double tmpforce2[1000][3];
  double tmpmax;
  double tmpmax2;
  int iresmax, imax, jmax, steptmp;
  int iresmax2, steptmp2, jresmax2;*/
};

}

#endif
#endif
