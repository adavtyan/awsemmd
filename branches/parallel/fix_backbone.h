/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Solvent Separated Barrier Potential was contributed by Nick Schafer

Membrane Potential was contributed by Leonardo Boechi and Bobby Kim

Last Update: 03/9/2012
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
  int z_res[2000];
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
  int **frustration_censoring_map;
  
  // Fragment Memory parameters
  double k_frag_mem;
  int n_frag_mems, **frag_mem_map, *ilen_fm_map;
  Fragment_Memory **frag_mems;
  Gamma_Array *fm_gamma;
  char frag_mems_file[100];
  char fm_gamma_file[100];
  double fm_sigma_exp;

  // Fragment Frustratometer parameters
  int n_decoy_mems, **decoy_mem_map, *ilen_decoy_map;
  Fragment_Memory **decoy_mems;
  char decoy_mems_file[100];
  int num_decoy_calcs;
  double **decoy_energy;
  int frag_frust_output_freq;
  char frag_frust_mode[100];
  int frag_frust_read_flag, frag_frust_shuffle_flag;
  double *frag_frust_read_mean, *frag_frust_read_variance;
  double frag_frust_well_width, frag_frust_seqsep_gamma;
  int frag_frust_seqsep_flag;
  bool frag_frust_normalizeInteraction; 

  // Tertiary Frustratometer parameters
  double tert_frust_cutoff;
  int tert_frust_ndecoys, tert_frust_output_freq;
  char tert_frust_mode[100];
  double *tert_frust_decoy_energies;
  double *decoy_ixn_stats;
  bool already_computed_configurational_decoys;

  // nmer frustratometer parameters
  int nmer_frust_size, nmer_frust_ndecoys, nmer_frust_output_freq;
  int nmer_contacts_cutoff;
  double nmer_frust_cutoff, nmer_frust_min_frust_threshold, nmer_frust_high_frust_threshold;
  double *nmer_frust_decoy_energies;
  double *nmer_decoy_ixn_stats;  
  char *nmer_seq_i, *nmer_seq_j, *nmer_seq_k;
  char *nmer_ss_i, *nmer_ss_j, *nmer_ss_k;
  bool nmer_output_neutral_flag;
  bool nmer_frust_trap_flag, nmer_frust_draw_trap_flag;
  double nmer_frust_trap_num_sigma, nmer_frust_ss_frac;
  char nmer_frust_mode[100];

  // Table Fragment Memory parameters
    TBV **fm_table;
//  TBV ****fm_table;
  int tb_size, tb_nbrs;
  double tb_rmin, tb_rmax, tb_dr;
  
  // Vector Fragment Memory
  double k_vec_frag_mem;
  double vfm_sigma, vfm_sigma_sq;
  double frag_table_well_width;
  int fm_energy_debug_flag;

  // Solvent separated barrier
  double k_solventb1, k_solventb2;
  double ssb_kappa, ssb_rmin1, ssb_rmax1, ssb_rmin2, ssb_rmax2;
  int ssb_ij_sep;
  bool ssb_rad_cor;
  double ssb_rshift[20];

  // Electrostatic Interaction from Huckel Method 
  double k_PlusPlus, k_MinusMinus, k_PlusMinus;
  double k_screening;
  double screening_length;
  double dielectric_constant, ion_concentration;
  double *charge_on_residue; 
  bool huckel_flag; //flag to turn on DebyeHuckel
  int debye_huckel_min_sep; // minimum sequence separation for DH interaction

  // Amylometer variables
  char amylometer_sequence_file[100];  
  int amylometer_nmer_size;
  int** nmer_array;
  int amylometer_mode;
  int number_of_nmers;
  char amylometer_structure_file[100];
  double amylometer_contact_cutoff;

  // Membrane potential
  double k_overall_memb;
  double memb_dens_offset;
  double k_bin;
  double memb_xo[3];
  int    memb_pore_type;
  double memb_len;
  double rho0_max;
  double rho0_distor;
  double g_memb[3][4];

  // Selection Temperature
  char selection_temperature_file_name[100];
  int selection_temperature_output_frequency;
  char selection_temperature_sequences_file_name[100];
  char selection_temperature_residues_file_name[100];
  char selection_temperature_sequence_energies_output_file_name[100];
  char selection_temperature_output_contact_list_file_name[100];
  char **selection_temperature_sequences;
  int num_selection_temperature_sequences;
  int num_selection_temperature_residues;
  int *selection_temperature_residues;
  bool selection_temperature_output_interaction_energies_flag;
  bool selection_temperature_output_contact_list_flag;
  bool selection_temperature_evaluate_sequence_energies_flag;
  double selection_temperature_rij_cutoff;
  int selection_temperature_min_seq_sep;

  // Monte Carlo Sequence Optimization
  double mcso_start_temp, mcso_end_temp;
  int mcso_num_steps;
  char mcso_se[1000];
  char mcso_seq_output_file_name[100];
  char mcso_energy_output_file_name[100];
  static const double k_b = 0.001987;

  // Optimization block parameters
  int optimization_output_freq;

  // Burial Optimization block parameters
  int burial_optimization_output_freq;

  // DebyeHuckel Optimization block parameters
  int debyehuckel_optimization_output_freq;
  char shuffler_mode[100];

  // Average Sequence Optimization parameters
  double **average_sequence;
  int average_sequence_optimization_output_freq;
  char average_sequence_input_file_name[100];

  // Mutate Sequence parameters
  char mutate_sequence_sequences_file_name[100];
  int mutate_sequence_number_of_sequences;
  char **mutate_sequence_sequences;
  int mutate_sequence_sequence_index;

  // Output_Per_Residue_Contacts parameters
  int output_per_residue_contacts_frequency;
  char output_per_residue_contacts_mode[100];
  char output_per_residue_contacts_native_structure_file_name[100];
  int output_per_residue_contacts_min_seq_sep;
  double output_per_residue_contacts_rij_threshold;
  char output_per_residue_contacts_file_name[100];
  Fragment_Memory *output_per_residue_contacts_structure;

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
  int *alpha_carbons, *beta_atoms, *oxygens;
  int *res_no,   *res_info, *chain_no;
  double **xca, **xcb, **xo, **xn, **xcp, **xh;
  int *res_no_l;

  double **x, **f;
  int *image;
  double prd[3], half_prd[3];
  int *periodicity;
  bool abc_flag, chain_flag, shake_flag, chi_flag, rama_flag, rama_p_flag, excluded_flag, p_excluded_flag, r6_excluded_flag;
  bool ssweight_flag, dssp_hdrgn_flag, p_ap_flag, water_flag, burial_flag, helix_flag, amh_go_flag, frag_mem_flag, ssb_flag;
  bool phosph_flag;
  bool vec_frag_mem_flag;
  bool amylometer_flag;
  bool frag_mem_tb_flag;
  bool memb_flag;
  bool frag_frust_flag;
  bool tert_frust_flag;
  bool nmer_frust_flag;
  bool selection_temperature_flag;
  bool frustration_censoring_flag;
  bool optimization_flag;
  bool burial_optimization_flag;
  bool debyehuckel_optimization_flag;
  bool average_sequence_optimization_flag;
  bool shuffler_flag;
  bool mutate_sequence_flag;
  bool monte_carlo_seq_opt_flag;
  bool output_per_residue_contacts_flag;

  enum Atoms{CA0 = 0, CA1, CA2, O0, O1, nAtoms};
  enum Angles{PHI = 0, PSI, nAngles};
  enum ResInfo{NONE=0, LOCAL, GHOST, OFF};
  
  char se[10000]; // Protein sequance
  int nch, ch_len[100], ch_pos[100];
  
  double energy[17], energy_all[17];
  enum EnergyTerms{ET_TOTAL=0, ET_CHAIN, ET_SHAKE, ET_CHI, ET_RAMA, ET_VEXCLUDED, ET_DSSP, ET_PAP, 
		   ET_WATER, ET_BURIAL, ET_HELIX, ET_AMHGO, ET_FRAGMEM, ET_VFRAGMEM, ET_MEMB, ET_SSB, ET_DH, nEnergyTerms};
  
  double ctime[17], previous_time;
  enum ComputeTime{TIME_CHAIN=0, TIME_SHAKE, TIME_CHI, TIME_RAMA, TIME_VEXCLUDED, TIME_DSSP, TIME_PAP, 
		   TIME_WATER, TIME_BURIAL, TIME_HELIX, TIME_AMHGO, TIME_FRAGMEM, TIME_VFRAGMEM, TIME_MEMB, TIME_SSB, TIME_DH, TIME_N};
  
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
  void compute_helix_potential(int i);
  void compute_amh_go_model();
  void compute_fragment_memory_potential(int i);
  void compute_decoy_memory_potential(int i, int decoy_calc);
  void randomize_decoys();
  void compute_generated_decoy_energies();
  void compute_fragment_frustration();
  void compute_solvent_barrier(int i, int j);
  void compute_fragment_memory_table();
  void output_fragment_memory_table();
  void table_fragment_memory(int i, int j);
  void compute_amhgo_normalization();
  void compute_vector_fragment_memory_potential(int i);
  void compute_amylometer();
  void read_amylometer_sequences(char *amylometer_sequence_file, int amylometer_nmer_size, int amylometer_mode);
  void compute_membrane_potential(int i);
  void compute_DebyeHuckel_Interaction(int i, int j);

  // Tertiary Frustratometer Functions
  void compute_tert_frust();
  double compute_native_ixn(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j);
  void compute_decoy_ixns(int i_resno_orig, int j_resno_orig, double rij_orig, double rho_i_orig, double rho_j_orig);
  double compute_water_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j);
  double compute_burial_energy(int i_resno, int ires_type, double rho_i);
  int get_random_residue_index();
  double get_residue_distance(int i_resno, int j_resno);
  double get_residue_density(int i);
  int get_residue_type(int i);
  double compute_frustration_index(double native_energy, double *decoy_stats);
  double compute_array_mean(double *array, int arraysize);
  double compute_array_std(double *array, int arraysize);
  void compute_tert_frust_singleresidue();
  double compute_singleresidue_native_ixn(int i_resno, int ires_type, double rho_i, int i_chno, double cutoff, bool nmercalc);
  void compute_singleresidue_decoy_ixns(int i_resno, double rho_i, int i_chno);
  // nmer frustratometer functions
  void compute_nmer_frust();
  void compute_singlenmer_frust();
  int compute_nmer_contacts(int i, int j);
  void get_nmer_seq(int i, char *nmer_seq, int backward);
  void get_nmer_secondary_structure(int i, char *nmer_secondary_structure);
  double compute_nmer_native_ixn(int i, int j);
  double compute_singlenmer_native_ixn(int i);
  void compute_nmer_decoy_ixns(int i, int j);
  void compute_singlenmer_decoy_ixns(int i);
  int compute_nmer_traps(int i, int j, int atomselect, double threshold_energy, char *nmer_seq_1, char *nmer_seq_2);
  int get_nmer_ss_dist(char *nmer_ss_j, char *nmer_ss_k);

  // Selection Temperature functions
  void output_selection_temperature_data();

  // Monte Carlo Sequence Optimization functions
  void compute_mcso();
  double compute_total_burial_energy();
  double compute_total_contact_energy();

  // Optimziation functions
  void compute_optimization();
  void shuffler();
  double compute_direct_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j);
  double compute_proteinmed_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j);
  double compute_watermed_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j);

  // Mutate_Sequence functions
  void mutate_sequence();

  // Average Sequence Optimization functions and variables
  void compute_average_sequence_optimization();
  
  // Burial Optimization functions
  void compute_burial_optimization();

  // DebyeHuckel Optimization functions
  void compute_debyehuckel_optimization();

  // Output_Per_Residue_Contacts functions
  void output_per_residue_contacts();

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
  FILE *fmenergiesfile;

  FILE *dout;
  int sStep, eStep;
  void print_forces(int coord=0);

  // Fragment Frustration files
  FILE *fragment_frustration_file;
  FILE *fragment_frustration_gap_file;
  FILE *fragment_frustration_variance_file;
  FILE *fragment_frustration_decoy_data;
  FILE *fragment_frustration_native_data;

  // Tertiary Frustration files
  FILE *tert_frust_output_file;
  FILE *tert_frust_vmd_script;

  // nmer frustration files
  FILE *nmer_frust_output_file;
  FILE *nmer_frust_vmd_script;
  FILE *nmer_frust_trap_file;

  // Selection temperature file
  FILE *selection_temperature_file;
  FILE *selection_temperature_sequence_energies_output_file;
  FILE *selection_temperature_contact_list_file;
  
  // Monte Carlo Sequence Optimization files
  FILE *mcso_seq_output_file;
  FILE *mcso_energy_output_file;

  // Optimization file
  FILE *optimization_file;
  FILE *native_optimization_file;
  FILE *optimization_norm_file;
  FILE *native_optimization_norm_file;
  
  // Burial Optimization file
  FILE *burial_optimization_file;
  FILE *native_burial_optimization_file;
  FILE *burial_optimization_norm_file;

  // DebyeHuckel Optimization file
  FILE *debyehuckel_optimization_file;
  FILE *debyehuckel_native_optimization_file;
  FILE *debyehuckel_optimization_norm_file;
  FILE *debyehuckel_native_optimization_norm_file;

  // Average Sequence Optimization files
  FILE *average_sequence_optimization_file;
  FILE *average_sequence_optimization_norm_file;

  // Mutate Sequences files
  FILE *mutate_sequences_file;

  // Output_Per_Residue_Contacts files
  FILE *output_per_residue_contacts_file;

  };
}

#endif
#endif
