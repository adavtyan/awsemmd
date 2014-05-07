/* ----------------------------------------------------------------------
   Copyright (2010) Aram Davtyan and Garegin Papoian

   Papoian's Group, University of Maryland at Collage Park
   http://papoian.chem.umd.edu/

   Solvent Separated Barrier Potential was contributed by Nick Schafer

   Membrane Potential was contributed by Leonardo Boechi and Bobby Kim

   Last Update: 03/9/2012
   ------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_backbone.h"
#include "atom.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "group.h"
#include "domain.h"
#include "memory.h"
#include "comm.h"
#include "timer.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using std::ifstream;

#define delta 0.00001

using namespace LAMMPS_NS;
using namespace FixConst;

//double fm_f[100][2][3], tfm_f[100][2][3];
//double err=0.0, err_max=0.0, err_max2=0.0;

/* ---------------------------------------------------------------------- */

// {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
// {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
int se_map[] = {0, 0, 4, 3, 6, 13, 7, 8, 9, 0, 11, 10, 12, 2, 0, 14, 5, 1, 15, 16, 0, 19, 17, 0, 18, 0};

// Four letter classes
// 1) SHL: Small Hydrophilic (ALA, GLY, PRO, SER THR) or (A, G, P, S, T) or {0, 7, 14, 15, 16}
// 2) AHL: Acidic Hydrophilic (ASN, ASP, GLN, GLU) or (N, D, Q, E) or {2, 3, 5, 6}
// 3) BAS: Basic (ARG HIS LYS) or (R, H, K) or {1, 8, 11}
// 4) HPB: Hydrophobic (CYS, ILE, LEU, MET, PHE, TRP, TYR, VAL) or (C, I, L, M, F, W, Y, V)  or {4, 9, 10, 12, 13, 17, 18, 19}
int bb_four_letter_map[] = {1, 3, 2, 2, 4, 2, 2, 1, 3, 4, 4, 3, 4, 4, 1, 1, 1, 4, 4, 4};

void itoa(int a, char *buf, int s)
{
  int b = abs(a);
  int c, i;
  i=0;
  while (b>0) {
    c = b - int(b/10)*10;
    b = b/10;
    buf[i] = c + '0';
    i++;
  }
  buf[i]='\0';
}

inline void FixBackbone::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

FixBackbone::FixBackbone(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix backbone command");
	
  efile = fopen("energy.log", "w");
  fmenergiesfile = fopen("fmenergies.log", "w");

  char buff[5];
  char forcefile[20]="";
  itoa(comm->me+1,buff,10);
  strcpy(forcefile,"forces\0");
  if (comm->nprocs>1) strcat(forcefile, buff);
  strcat(forcefile, ".dat");
  dout = fopen(forcefile, "w");
	
  char eheader[] = "Step   \tChain   \tShake   \tChi     \tRama    \tExcluded\tDSSP    \tP_AP    \tWater   \tBurial  \tHelix   \tAMH-Go  \tFrag_Mem\tVec_FM  \tMembrane\tSSB     \tElectro.\tVTotal\n";
  fprintf(efile, "%s", eheader);

  scalar_flag = 1;
  vector_flag = 1;
  thermo_energy = 1;
  size_vector = nEnergyTerms-1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
	
  abc_flag = chain_flag = shake_flag = chi_flag = rama_flag = rama_p_flag = excluded_flag = p_excluded_flag = r6_excluded_flag = 0;
  ssweight_flag = dssp_hdrgn_flag = p_ap_flag = water_flag = burial_flag = helix_flag = amh_go_flag = frag_mem_flag = vec_frag_mem_flag = 0;
  ssb_flag = frag_mem_tb_flag = phosph_flag = amylometer_flag = memb_flag = selection_temperature_flag = 0;
  frag_frust_flag = tert_frust_flag = nmer_frust_flag = optimization_flag = burial_optimization_flag = 0;
  huckel_flag = debyehuckel_optimization_flag = 0;
  average_sequence_optimization_flag = 0;
  shuffler_flag = 0;
  mutate_sequence_flag = 0;
  monte_carlo_seq_opt_flag = 0;
  output_per_residue_contacts_flag = 0;

  epsilon = 1.0; // general energy scale
  p = 2; // for excluded volume
	
  int i, j;
	
  for (i=0;i<12;i++) ssweight[i] = false;
	
  for (i=0;i<TIME_N;i++) ctime[i] = 0.0;

  // backbone geometry coefficients
  an = 0.4831806; bn = 0.7032820; cn = -0.1864262;
  ap = 0.4436538; bp = 0.2352006; cp = 0.3211455;
  ah = 0.8409657; bh = 0.8929599; ch = -0.7338894;
	
  igroup2 = group->find(arg[3]);
  if (igroup2 == -1) 
    error->all(FLERR,"Could not find fix backbone beta atoms group ID"); 
  igroup3 = group->find(arg[4]);
  if (igroup3 == -1) 
    error->all(FLERR,"Could not find fix backbone oxygen atoms group ID"); 
  if (igroup2 == igroup || igroup3 == igroup || igroup2 == igroup3) 
    error->all(FLERR,"Two groups cannot be the same in fix backbone"); 
  if (group->count(igroup)!=group->count(igroup2) || group->count(igroup2)!=group->count(igroup3))
    error->all(FLERR,"All groups must contain the same # of atoms in fix backbone");
  group2bit = group->bitmask[igroup2];
  group3bit = group->bitmask[igroup3];
	
  char varsection[30];
  ifstream in(arg[5]);
  if (!in) error->all(FLERR,"Coefficient file was not found!");
  while (!in.eof()) {
    in >> varsection;
    if (strcmp(varsection, "[ABC]")==0) {
      abc_flag = 1;
      if (comm->me==0) print_log("ABC flag on\n");
      in >> an >> bn >> cn;
      in >> ap >> bp >> cp;
      in >> ah >> bh >> ch;
    } else if (strcmp(varsection, "[Chain]")==0) {
      chain_flag = 1;
      if (comm->me==0) print_log("Chain flag on\n");
      in >> k_chain[0] >> k_chain[1] >> k_chain[2]; 
      in >> r_ncb0 >> r_cpcb0 >> r_ncp0;
    } else if (strcmp(varsection, "[Shake]")==0) {
      shake_flag = 1;
      if (comm->me==0) print_log("Shake flag on\n");
      in >> k_shake >> r_sh1 >> r_sh2 >> r_sh3;
    } else if (strcmp(varsection, "[Chi]")==0) {
      chi_flag = 1;
      if (comm->me==0) print_log("Chi flag on\n");
      in >> k_chi >> chi0;
    } else if (strcmp(varsection, "[Excluded]")==0) {
      excluded_flag = 1;
      if (comm->me==0) print_log("Excluded flag on\n");
      in >> k_excluded_C >> rC_ex0;
      in >> k_excluded_O >> rO_ex0;
    } else if (strcmp(varsection, "[Excluded_P]")==0) {
      p_excluded_flag = 1;
      if (comm->me==0) print_log("Excluded_P flag on\n");
      in >> p;
      in >> k_excluded_C >> rC_ex0;
      in >> k_excluded_O >> rO_ex0;
    } else if (strcmp(varsection, "[Excluded_R6]")==0) {
      r6_excluded_flag = 1;
      if (comm->me==0) print_log("Excluded_R6 flag on\n");
      in >> k_excluded_C >> rC_ex0;
      in >> k_excluded_O >> rO_ex0;
    } else if (strcmp(varsection, "[Rama]")==0) {
      rama_flag = 1;
      if (comm->me==0) print_log("Rama flag on\n");
      in >> k_rama;
      in >> n_rama_par;
      for (int j=0;j<n_rama_par;j++) {
	in >> w[j] >> sigma[j] >> phiw[j] >> phi0[j] >> psiw[j] >> psi0[j];
      }
    } else if (strcmp(varsection, "[Rama_P]")==0) {
      rama_p_flag = 1;
      if (comm->me==0) print_log("Rama_P flag on\n");
      in >> n_rama_p_par;
      for (int j=0;j<n_rama_p_par;j++) {
	in >> w[j+i_rp] >> sigma[j+i_rp] >> phiw[j+i_rp] >> phi0[j+i_rp] >> psiw[j+i_rp] >> psi0[j+i_rp];
      }
    } else if (strcmp(varsection, "[SSWeight]")==0) {
      ssweight_flag = 1;
      if (comm->me==0) print_log("SSWeight flag on\n");
      for (int j=0;j<12;++j)
	in >> ssweight[j];
    } else if (strcmp(varsection, "[Dssp_Hdrgn]")==0) {
      dssp_hdrgn_flag = 1;
      if (comm->me==0) print_log("Dssp_Hdrgn flag on\n");
      in >> k_dssp;
      in >> hbscl[0][0] >> hbscl[0][1];
      for (int j=0;j<7;++j) in >> hbscl[1][j];
      for (int j=0;j<9;++j) in >> hbscl[2][j];
      for (int j=0;j<9;++j) in >> hbscl[3][j];
      in >> sigma_HO >> sigma_NO;
      in >> HO_zero >> NO_zero;
      in >> dssp_hdrgn_cut;
      in >> pref[0] >> pref[1];
      in >> d_nu0;
    } else if (strcmp(varsection, "[P_AP]")==0) {
      p_ap_flag = 1;
      if (comm->me==0) print_log("P_AP flag on\n");
      in >> k_global_P_AP;
      in >> k_betapred_P_AP;
      in >> k_P_AP[0] >> k_P_AP[1] >> k_P_AP[2];
      in >> P_AP_cut;
      in >> P_AP_pref;
      in >> i_med_min >> i_med_max;
      in >> i_diff_P_AP;
    } else if (strcmp(varsection, "[Water]")==0) {
      water_flag = 1;
      if (comm->me==0) print_log("Water flag on\n");
      in >> k_water;
      in >> water_kappa >> water_kappa_sigma;
      in >> treshold;
      in >> contact_cutoff;
      in >> n_wells;
      for (int j=0;j<n_wells;++j)
	in >> well_r_min[j] >> well_r_max[j] >> well_flag[j];
    } else if (strcmp(varsection, "[Burial]")==0) {
      burial_flag = 1;
      if (comm->me==0) print_log("Burial flag on\n");
      in >> k_burial;
      in >> burial_kappa;
      in >> burial_ro_min[0] >> burial_ro_max[0];
      in >> burial_ro_min[1] >> burial_ro_max[1];
      in >> burial_ro_min[2] >> burial_ro_max[2];
    } else if (strcmp(varsection, "[Helix]")==0) {
      helix_flag = 1;
      if (comm->me==0) print_log("Helix flag on\n");
      in >> k_helix;
      in >> helix_gamma_p >> helix_gamma_w;
      in >> helix_kappa >> helix_kappa_sigma;
      in >> helix_treshold;
      in >> helix_i_diff;
      in >> helix_cutoff;
      n_helix_wells = 1;
      helix_well_flag[0] = 1;
      in >> helix_well_r_min[0] >> helix_well_r_max[0];
      for (int j=0;j<20;++j)
	in >> h4prob[j];
      // h4prob coefficent for proline if it is aceptor
      // It will be used only if pro_accepter_flag=1
      in >> pro_accepter_flag >> h4prob_pro_accepter;
      in >> helix_sigma_HO >> helix_sigma_NO;
      in >> helix_HO_zero >> helix_NO_zero;
    } else if (strcmp(varsection, "[AMH-Go]")==0) {
      amh_go_flag = 1;
      if (comm->me==0) print_log("AMH-Go flag on\n");
      in >> k_amh_go;
      in >> amh_go_p;
      in >> amh_go_rc;
      in >> frustration_censoring_flag;
    } else if (strcmp(varsection, "[Fragment_Memory]")==0) {
      frag_mem_flag = 1;
      if (comm->me==0) print_log("Fragment_Memory flag on\n");
      in >> k_frag_mem;
      in >> frag_mems_file;
      in >> fm_gamma_file;
    } else if (strcmp(varsection, "[Fragment_Memory_Table]")==0) {
      frag_mem_tb_flag = 1;
      if (comm->me==0) print_log("Fragment_Memory_Table flag on\n");
      in >> k_frag_mem;
      in >> frag_mems_file;
      in >> fm_gamma_file;
      in >> tb_rmin >> tb_rmax >> tb_dr;
      tb_size = (int)((tb_rmax-tb_rmin)/tb_dr)+2;
      in >> frag_table_well_width;
      in >> fm_energy_debug_flag;
      in >> fm_sigma_exp;      
    } else if (strcmp(varsection, "[Vector_Fragment_Memory]")==0) {
      vec_frag_mem_flag = 1;
      if (comm->me==0) print_log("Vector Fragment_Memory flag on\n");
      in >> k_vec_frag_mem;
      in >> vfm_sigma;
      vfm_sigma_sq = vfm_sigma*vfm_sigma;
    } else if (strcmp(varsection, "[Solvent_Barrier]")==0) {
      ssb_flag = 1;
      if (comm->me==0) print_log("Solvent separated barrier flag on\n");
      in >> k_solventb1;
      in >> ssb_rmin1 >> ssb_rmax1;
      in >> k_solventb2;
      in >> ssb_rmin2 >> ssb_rmax2;
      in >> ssb_kappa;
      in >> ssb_ij_sep;
      in >> ssb_rad_cor;
      for (int j=0;j<20;++j)
        in >> ssb_rshift[j];
    } else if (strcmp(varsection, "[Membrane]")==0) {
      memb_flag = 1;
      print_log("Membrane flag on\n");
      in >> k_overall_memb;
      in >> memb_dens_offset;
      in >> k_bin;
      in >> memb_xo[0] >> memb_xo[1] >> memb_xo[2];
      in >> memb_pore_type;
      in >> memb_len;
      in >> rho0_max;
      in >> rho0_distor;
      for (int i=0;i<3;++i)
        for (int j=0;j<4;++j) 
          in >> g_memb[i][j]; 
    } else if (strcmp(varsection, "[Fragment_Frustratometer]")==0) {
      // The fragment frustratometer requires the fragment memory potential to be active
      if (!frag_mem_flag && !frag_mem_tb_flag) error->all(FLERR,"Cannot run Fragment_Frustratometer without Fragment_Memory or Fragment_Memory_Table.");
      frag_frust_flag = 1; // activate flag for fragment frustratometer
      if (comm->me==0) print_log("Fragment_Frustratometer flag on\n");
      in >> frag_frust_mode; // the possible modes are "read" and "shuffle"
      if (strcmp(frag_frust_mode, "shuffle")==0) {
	if (comm->me==0) print_log("Fragment_Frustratometer in shuffle mode\n");
	frag_frust_shuffle_flag=1; // activate "shuffle" specific flag
	in >> decoy_mems_file; // read in the decoy fragments that will be shuffled to generated decoy energies
	in >> num_decoy_calcs; // this is the number of times that the decoy fragments will be shuffled
	in >> frag_frust_output_freq; // this is the number of steps between frustration calculations
      }
      else if (strcmp(frag_frust_mode, "read")==0) {
	if (comm->me==0) print_log("Fragment_Frustratometer in read mode\n");
	frag_frust_read_flag=1; // activate "read" specific flag
	in >> decoy_mems_file; // read in the decoy structures that will be used to generate the decoy energies
	in >> frag_frust_output_freq; // this is the number of steps between frustration calculations
	in >> frag_frust_well_width; // parameter to tune well width, default is 1.0
	in >> frag_frust_seqsep_flag >> frag_frust_seqsep_gamma; // flag and parameter to tune sequence separation dependent gamma
	in >> frag_frust_normalizeInteraction; // flag that determines whether or not the fragment interaction is normalized by the width of the interaction
      }
      else {
	// throw an error if the "mode" is anything but "read" or "shuffle"
	error->all(FLERR,"Only \"shuffle\" and \"read\" are acceptable modes for the Fragment_Frustratometer.");
      }
    } else if (strcmp(varsection, "[Tertiary_Frustratometer]")==0) {
      tert_frust_flag = 1;
      if (comm->me==0) print_log("Tertiary_Frustratometer flag on\n");
      in >> tert_frust_cutoff;
      in >> tert_frust_ndecoys;
      in >> tert_frust_output_freq;
      in >> tert_frust_mode;
      // Set the value of this flag to 0 so that the configurational decoy statistics will be computed at least once
      already_computed_configurational_decoys = 0;
      if (strcmp(tert_frust_mode, "configurational")!=0 && strcmp(tert_frust_mode, "mutational")!=0 && strcmp(tert_frust_mode, "singleresidue")!=0) {
	// throw an error if the "mode" is anything but "configurational" or "mutational"
	error->all(FLERR,"Only \"configurational\", \"mutational\", \"singleresidue\" are acceptable modes for the Tertiary_Frustratometer.");
      }
    } else if (strcmp(varsection, "[Nmer_Frustratometer]")==0) {
      nmer_frust_flag = 1;
      if (comm->me==0) print_log("Nmer_Frustratometer flag on\n");
      in >> nmer_frust_size;
      in >> nmer_frust_cutoff;
      in >> nmer_contacts_cutoff;
      in >> nmer_frust_ndecoys;
      in >> nmer_frust_output_freq;
      in >> nmer_frust_min_frust_threshold >> nmer_frust_high_frust_threshold >> nmer_output_neutral_flag;
      in >> nmer_frust_trap_flag >> nmer_frust_draw_trap_flag >> nmer_frust_trap_num_sigma >> nmer_frust_ss_frac;
      in >> nmer_frust_mode;
      if (strcmp(nmer_frust_mode, "pairwise")!=0 && strcmp(nmer_frust_mode, "singlenmer")!=0) {
	// throw an error if the "mode" is anything but "configurational" or "mutational"
	error->all(FLERR,"Only \"pairwise\", \"singlenmer\" are acceptable modes for the Nmer_Frustratometer.");
      }	    
    } else if (strcmp(varsection, "[Phosphorylation]")==0) {
      if (!water_flag) error->all(FLERR,"Cannot run phosphorylation without water potential");
      phosph_flag = 1;
      if (comm->me==0) print_log("Phosphorylation flag on\n");
      in >> k_hypercharge;
      in >> n_phosph_res;
      if (n_phosph_res > 20) error->all(FLERR,"Number of phosphorylated residues may not exceed 20");
      for (int i=0;i<n_phosph_res;++i)
	in >> phosph_res[i]; 
    } else if (strcmp(varsection, "[Epsilon]")==0) {
      in >> epsilon;
    } else if (strcmp(varsection, "[Amylometer]")==0) {
      amylometer_flag = 1;
      if (comm->me==0) print_log("Amylometer flag on\n");
      in >> amylometer_sequence_file;
      in >> amylometer_nmer_size;
      // 1 == self-only, 2 == heterogeneous
      in >> amylometer_mode;
      if (amylometer_mode == 2) {
	in >> amylometer_structure_file;
	in >> amylometer_contact_cutoff;
      }
      read_amylometer_sequences(amylometer_sequence_file, amylometer_nmer_size, amylometer_mode);
    } else if (strcmp(varsection, "[Selection_Temperature]")==0) {
      selection_temperature_flag = 1;
      if (comm->me==0) print_log("Selection_Temperature flag on \n");
      // outputting interaction energies
      in >> selection_temperature_output_frequency;
      in >> selection_temperature_output_interaction_energies_flag;
      in >> selection_temperature_file_name;
      // evaluating multiple sequence energies
      in >> selection_temperature_evaluate_sequence_energies_flag;
      in >> selection_temperature_sequences_file_name;
      in >> selection_temperature_residues_file_name;
      in >> selection_temperature_sequence_energies_output_file_name;
      // outputting contact lists
      in >> selection_temperature_output_contact_list_flag;
      in >> selection_temperature_rij_cutoff;
      in >> selection_temperature_min_seq_sep;
      in >> selection_temperature_output_contact_list_file_name;
    } else if (strcmp(varsection, "[Monte_Carlo_Seq_Opt]")==0) {
      monte_carlo_seq_opt_flag = 1;
      if (comm->me==0) print_log("Monte_Carlo_Seq_Opt flag on \n");
      in >> mcso_start_temp >> mcso_end_temp >> mcso_num_steps;
      in >> mcso_seq_output_file_name;
      in >> mcso_energy_output_file_name;
    } else if (strcmp(varsection, "[Optimization]")==0) {
      optimization_flag = 1;
      if (comm->me==0) print_log("Optimization flag on\n");
      in >> optimization_output_freq;
    } else if (strcmp(varsection, "[Average_Sequence_Optimization]")==0) {
      average_sequence_optimization_flag = 1;
      if (comm->me==0) print_log("Average Sequence Optimization flag on\n");
      in >> average_sequence_optimization_output_freq;
      in >> average_sequence_input_file_name;
    }
    else if (strcmp(varsection, "[Burial_Optimization]")==0) {
      burial_optimization_flag = 1;
      if (comm->me==0) print_log("Burial Optimization flag on\n");
      in >> burial_optimization_output_freq;
    } else if (strcmp(varsection, "[DebyeHuckel]")==0) {
      huckel_flag = 1;
      if (comm->me==0) print_log("DebyeHuckel on\n");
      in >> k_PlusPlus >> k_MinusMinus >> k_PlusMinus;
      in >> k_screening;
      in >> screening_length;
      fprintf(screen, "Debye-Huckel Screening Length = %8.6f Angstroms\n", screening_length);  
      in >> debye_huckel_min_sep;
    } else if (strcmp(varsection, "[DebyeHuckel_Optimization]")==0) {
      debyehuckel_optimization_flag = 1;
      if (comm->me==0) print_log("DebyeHuckel_Optimization flag on\n");
      in >> debyehuckel_optimization_output_freq;
    } else if (strcmp(varsection, "[Shuffler]")==0) {
      in >> shuffler_flag;
      in >> shuffler_mode;
      if ( shuffler_flag == 1 ) {
	if (comm->me==0) print_log("Shuffler flag on\n");
      }
    } else if (strcmp(varsection, "[Mutate_Sequence]")==0) {
      in >> mutate_sequence_flag;
      in >> mutate_sequence_sequences_file_name;
      if ( mutate_sequence_flag == 1 ) {
	if (comm->me==0) print_log("Mutate_Sequence flag on\n");
      }
    } else if (strcmp(varsection, "[Output_Per_Residue_Contacts]")==0) {
      output_per_residue_contacts_flag = 1;
      if (comm->me==0) print_log("Output_Per_Residue_Contacts flag on\n");
      in >> output_per_residue_contacts_frequency;
      in >> output_per_residue_contacts_mode;
      if (strcmp(output_per_residue_contacts_mode, "native")!=0 && strcmp(output_per_residue_contacts_mode, "all")!=0) {
	error->all(FLERR,"Only \"native\" or \"all\" are acceptable modes for Output_Per_Residue_Contacts.");
      }
      if (strcmp(output_per_residue_contacts_mode, "native")==0) {
	in >> output_per_residue_contacts_native_structure_file_name;
      }
      in >> output_per_residue_contacts_min_seq_sep >> output_per_residue_contacts_rij_threshold;
      in >> output_per_residue_contacts_file_name;
    } 
      
    varsection[0]='\0'; // Clear buffer
  }
  in.close();
  if (comm->me==0) print_log("\n");

  // Read sequance file
  ifstream ins(arg[6]);
  if (!ins) error->all(FLERR,"Sequence file was not found");
  char buf[10000];
  se[0]='\0';
  nch = 0;
  while (!ins.eof()) {
    ins >> buf;
    if (buf[0]=='#' || isEmptyString(buf)) continue;
    ch_pos[nch] = strlen(se)+1;
    strcat(se, buf);
    ch_len[nch] = strlen(buf);
    nch++;
    buf[0]='\0';
  }
  ins.close();
	
  force_flag = 0;
  n = (int)(group->count(igroup)+1e-12);
  for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;
  x = atom->x;
  f = atom->f;
  image = atom->image;
  prd[0] = domain->xprd;
  prd[1] = domain->yprd;
  prd[2] = domain->zprd;
  half_prd[0] = prd[0]/2;
  half_prd[1] = prd[1]/2;
  half_prd[2] = prd[2]/2; 
  periodicity = domain->periodicity;
  allocated = false;

  allocate();

  if (dssp_hdrgn_flag) {
    ifstream in_anti_HB("anti_HB");
    ifstream in_anti_NHB("anti_NHB");
    ifstream in_para_HB("para_HB");
    ifstream in_para_one("para_one");
    ifstream in_anti_one("anti_one");
    
    if (!in_anti_HB) error->all(FLERR,"File anti_HB doesn't exist");
    if (!in_anti_NHB) error->all(FLERR,"File anti_NHB doesn't exist");
    if (!in_para_HB) error->all(FLERR,"File para_HB doesn't exist");
    if (!in_para_one) error->all(FLERR,"File para_one doesn't exist");
    if (!in_anti_one) error->all(FLERR,"File anti_one doesn't exist");
    
    for (i=0;i<20;++i) {
      in_para_one >> m_para_one[i];
      in_anti_one >> m_anti_one[i];
      for (j=0;j<20;++j) {
        in_anti_HB >> m_anti_HB[i][j][0];
        in_anti_NHB >> m_anti_NHB[i][j][0];
        in_para_HB >> m_para_HB[i][j][0];
      }
    }
    for (i=0;i<20;++i) {
      for (j=0;j<20;++j) {
        in_anti_HB >> m_anti_HB[i][j][1];
        in_anti_NHB >> m_anti_NHB[i][j][1];
        in_para_HB >> m_para_HB[i][j][1];
      }
    }
    in_anti_HB.close();
    in_anti_NHB.close();
    in_para_HB.close();
    in_para_one.close();
    in_anti_one.close();
  }

  if (ssweight_flag) {
    ifstream in_ssw("ssweight");
    if (!in_ssw) error->all(FLERR,"File ssweight doesn't exist");
    for (j=0;j<n;++j) {
      for (i=0;i<12;++i) {
	if (ssweight[i]) in_ssw >> aps[i][j]; else aps[i][j] = 0.0;
      }
    }
    in_ssw.close();
  }

  if (memb_flag) {
    ifstream in_memb_zim("zim");
    if (!in_memb_zim) error->all(FLERR,"File zim doesn't exist");
    // what's happen if zim file is not correct
    for (i=0;i<n;++i) {
      in_memb_zim >> z_res[i];
    }
    in_memb_zim.close();
  }


  if (water_flag) {
    ifstream in_wg("gamma.dat");
    if (!in_wg) error->all(FLERR,"File gamma.dat doesn't exist");
    for (int i_well=0;i_well<n_wells;++i_well) {
      for (i=0;i<20;++i) {
	for (j=i;j<20;++j) {
	  in_wg >> water_gamma[i_well][i][j][0] >> water_gamma[i_well][i][j][1];
	  water_gamma[i_well][j][i][0] = water_gamma[i_well][i][j][0];
	  water_gamma[i_well][j][i][1] = water_gamma[i_well][i][j][1];
	}
      }
    }
    in_wg.close();
  }
	
  if (phosph_flag) {
    for (int i_well=0;i_well<n_wells;++i_well) {
      for (i=0;i<20;++i) {
	for (j=i;j<20;++j) {
	  phosph_water_gamma[i_well][i][j][0] = phosph_water_gamma[i_well][j][i][0] = water_gamma[i_well][i][j][0];
	  phosph_water_gamma[i_well][i][j][1] = phosph_water_gamma[i_well][j][i][1] = water_gamma[i_well][i][j][1];
	}
      }
    }
	  

    //replacing serine interaction gammas with hypercharged glutamate interaction gammas
    for (int i_well=0;i_well<n_wells;++i_well) {
      for (i=0;i<20;++i) {
	if (bb_four_letter_map[i]==1) {
	  phosph_water_gamma[i_well][i][15][0] = phosph_water_gamma[i_well][15][i][0] = phosph_water_gamma[i_well][i][6][0]*k_hypercharge;
	  phosph_water_gamma[i_well][i][15][1] = phosph_water_gamma[i_well][15][i][1] = phosph_water_gamma[i_well][i][6][1]*k_hypercharge;
	}
	else if (bb_four_letter_map[i]==2 || bb_four_letter_map[i]==3) {
	  phosph_water_gamma[i_well][i][15][0] = phosph_water_gamma[i_well][15][i][0] = phosph_water_gamma[i_well][i][6][0]*pow(k_hypercharge,2);
	  phosph_water_gamma[i_well][i][15][1] = phosph_water_gamma[i_well][15][i][1] = phosph_water_gamma[i_well][i][6][1]*pow(k_hypercharge,2);
	}
	else {
	  phosph_water_gamma[i_well][i][15][0] = phosph_water_gamma[i_well][15][i][0] = phosph_water_gamma[i_well][i][6][0];
	  phosph_water_gamma[i_well][i][15][1] = phosph_water_gamma[i_well][15][i][1] = phosph_water_gamma[i_well][i][6][1];
	}
      }
    }
    //create map of phosphorylated residues
    phosph_map = new int[n];
    for (int i=0;i<n;++i) {
      phosph_map[i]=0;
    }  
    for (int j=0;j<n_phosph_res;++j) {
      if (phosph_res[j]!=0) {
	int dummy = phosph_res[j]-1;
	phosph_map[dummy]=1;
      }
    }	  
  }
	
  if (burial_flag) {
    ifstream in_brg("burial_gamma.dat");
    if (!in_brg) error->all(FLERR,"File burial_gamma.dat doesn't exist");
    for (i=0;i<20;++i) {
      in_brg >> burial_gamma[i][0] >> burial_gamma[i][1] >> burial_gamma[i][2];
    }
    in_brg.close();
  }
	
  if (amh_go_flag) {
    char amhgo_gamma_file[] = "amh-go.gamma";
    amh_go_gamma = new Gamma_Array(amhgo_gamma_file);
    if (amh_go_gamma->error==amh_go_gamma->ERR_FILE) error->all(FLERR,"Cannot read file amh-go.gamma");
    if (amh_go_gamma->error==amh_go_gamma->ERR_CLASS_DEF) error->all(FLERR,"AMH_Go: Wrong definition of sequance separation classes");
    if (amh_go_gamma->error==amh_go_gamma->ERR_GAMMA) error->all(FLERR,"AMH_Go: Incorrect entery in gamma file");
    if (amh_go_gamma->error==amh_go_gamma->ERR_G_CLASS) error->all(FLERR,"AMH_Go: Wrong sequance separation class tag");
    if (amh_go_gamma->error==amh_go_gamma->ERR_ASSIGN) error->all(FLERR,"AMH_Go: Cannot build gamma array");
    
    char amhgo_mem_file[] = "amh-go.gro";
    m_amh_go = new Fragment_Memory(0, 0, n, 1.0, amhgo_mem_file);
    if (m_amh_go->error==m_amh_go->ERR_FILE) error->all(FLERR,"Cannot read file amh-go.gro");
    if (m_amh_go->error==m_amh_go->ERR_ATOM_COUNT) error->all(FLERR,"AMH_Go: Wrong atom count in memory structure file");
    if (m_amh_go->error==m_amh_go->ERR_RES) error->all(FLERR,"AMH_Go: Unknown residue");

    // if frustration censoring is on, read in frustration censored interactions
    if (frustration_censoring_flag == 1) {
      std::ifstream infile("frustration_censored_contacts.dat");
      int i, j;
      while(infile >> i >> j) {
	frustration_censoring_map[i-1][j-1] = 1;
      }
    }

    // Calculate normalization factor for AMH-GO potential
    compute_amhgo_normalization();
  }
	
	
  if (frag_mem_flag || frag_mem_tb_flag) {
    if (comm->me==0) print_log("Reading fragments...\n");
		
    fm_gamma = new Gamma_Array(fm_gamma_file);
    if (fm_gamma->error==fm_gamma->ERR_FILE) error->all(FLERR,"Fragment_Memory: Cannot read gamma file");
    if (fm_gamma->error==fm_gamma->ERR_CLASS_DEF) error->all(FLERR,"Fragment_Memory: Wrong definition of sequance separation classes");
    if (fm_gamma->error==fm_gamma->ERR_GAMMA) error->all(FLERR,"Fragment_Memory: Incorrect entery in gamma file");
    if (fm_gamma->error==fm_gamma->ERR_G_CLASS) error->all(FLERR,"Fragment_Memory: Wrong sequance separation class tag");
    if (fm_gamma->error==fm_gamma->ERR_ASSIGN) error->all(FLERR,"Fragment_Memory: Cannot build gamma array");
		
    // read frag_mems_file and create a list of the fragments
    frag_mems = read_mems(frag_mems_file, n_frag_mems);
		
    // allocate frag_mem_map and ilen_fm_map
    ilen_fm_map = new int[n]; // Number of fragments for residue i
    frag_mem_map = new int*[n]; // Memory Fragments map
    for (i=0;i<n;++i) {
      ilen_fm_map[i] = 0;
      frag_mem_map[i] = NULL;
    }
		
    // Fill Fragment Memory map
    int k, pos, len, min_sep;
    min_sep = fm_gamma->minSep();
    for (k=0;k<n_frag_mems;++k) {
      pos = frag_mems[k]->pos;
      len = frag_mems[k]->len;
      
      if (pos+len>n) {
	fprintf(stderr, "pos %d len %d n %d\n", pos, len, n); 
	error->all(FLERR,"Fragment_Memory: Incorrectly defined memory fragment");
      }
      
      for (i=pos; i<pos+len-min_sep; ++i) {
	ilen_fm_map[i]++;
	frag_mem_map[i] = (int *) memory->srealloc(frag_mem_map[i],ilen_fm_map[i]*sizeof(int),"modify:frag_mem_map");
	frag_mem_map[i][ilen_fm_map[i]-1] = k;
      }
    }
  }
  
  // if the fragment frustratometer flag is on, perform appropriate initializations
  if (frag_frust_flag) {
    // open fragment frustration file for writing
    fragment_frustration_file = fopen("fragment_frustration.dat","w");
    fragment_frustration_gap_file = fopen("fragment_frustration_gap.dat","w");
    fragment_frustration_variance_file = fopen("fragment_frustration_variance.dat","w");
    fragment_frustration_decoy_data = fopen("fragment_frustration_decoy.dat","w");
    fragment_frustration_native_data = fopen("fragment_frustration_native.dat","w");

    if (comm->me==0) print_log("Reading decoy fragments...\n");
    // create a decoy memory array by reading in the appropriate file
    decoy_mems = read_mems(decoy_mems_file, n_decoy_mems); // n_decoy_mems is set equal to the number of decoys in the read_mems function
    // because the number of decoy calculations is set by the size of the decoy list in "read" mode, we need to initialize the variable here
    if (frag_frust_read_flag) {
      num_decoy_calcs = n_decoy_mems+1; // add one so that the "native" energy can occupy the 0 index
    }

    // allocate decoy_mem_map and ilen_decoy_map
    ilen_decoy_map = new int[n]; // Number of decoys for residue i
    decoy_mem_map = new int*[n]; // decoy Memory Fragments map
    for (i=0;i<n;++i) {
      ilen_decoy_map[i] = 0;
      decoy_mem_map[i] = NULL;
    }
	
    // Fill Decoy Memory map
    int k, pos, len, min_sep;
    min_sep = fm_gamma->minSep();
    for (k=0;k<n_decoy_mems;++k) {
      pos = decoy_mems[k]->pos;
      len = decoy_mems[k]->len;
      
      if (pos+len>n) {
	fprintf(stderr, "pos %d len %d n %d\n", pos, len, n); 
	error->all(FLERR,"Fragment_Frustratometer: Incorrectly defined memory fragment");
      }
      
      for (i=pos; i<pos+len-min_sep; ++i) {
	ilen_decoy_map[i]++;
	decoy_mem_map[i] = (int *) memory->srealloc(decoy_mem_map[i],ilen_decoy_map[i]*sizeof(int),"modify:decoy_mem_map");
	decoy_mem_map[i][ilen_decoy_map[i]-1] = k;
      }
    }

    // Allocate decoy_energy array
    decoy_energy = new double*[n];
    for (i=0;i<n;i++)
      {
	decoy_energy[i] = new double[num_decoy_calcs];
	for(int decoyindex=0; decoyindex<num_decoy_calcs; decoyindex++)
	  {
	    decoy_energy[i][decoyindex] = 0.0;
	  }
      }
    // if in "read" mode, allocate per residue mean and variance arrays and compute generated decoy energies
    if (frag_frust_read_flag)
      {
	frag_frust_read_mean = new double[n];
	frag_frust_read_variance = new double[n];
      }
  }

  // if tert_frust_flag is on, perform appropriate initializations
  if(tert_frust_flag) {
    tert_frust_decoy_energies = new double[tert_frust_ndecoys];
    decoy_ixn_stats = new double[2];
    tert_frust_output_file = fopen("tertiary_frustration.dat","w");
    tert_frust_vmd_script = fopen("tertiary_frustration.tcl","w");
    if (strcmp(tert_frust_mode, "configurational")==0 || strcmp(tert_frust_mode, "mutational")==0) {
      fprintf(tert_frust_output_file,"# i j i_chain j_chain xi yi zi xj yj zj r_ij rho_i rho_j a_i a_j native_energy <decoy_energies> std(decoy_energies) f_ij\n");
    }
    else if (strcmp(tert_frust_mode, "singleresidue")==0) {
      fprintf(tert_frust_output_file,"# i i_chain xi yi zi rho_i a_i native_energy <decoy_energies> std(decoy_energies) f_i\n");
    }
  }

  // if nmer_frust_flag is on, perform appropriate initializations
  if(nmer_frust_flag) {
    nmer_frust_decoy_energies = new double[nmer_frust_ndecoys];
    nmer_decoy_ixn_stats = new double[2];
    nmer_seq_i = new char[nmer_frust_size+1]; // extend the array
    nmer_seq_i[nmer_frust_size] = '\0';       // and null terminate it so that it can be printed properly
    nmer_seq_j = new char[nmer_frust_size+1]; // extend the array
    nmer_seq_j[nmer_frust_size] = '\0';       // and null terminate it so that it can be printed properly
    nmer_seq_k = new char[nmer_frust_size+1]; // extend the array
    nmer_seq_k[nmer_frust_size] = '\0';       // and null terminate it so that it can be printed properly
    nmer_ss_i = new char[nmer_frust_size+1]; // extend the array
    nmer_ss_i[nmer_frust_size] = '\0';       // and null terminate it so that it can be printed properly
    nmer_ss_j = new char[nmer_frust_size+1]; // extend the array
    nmer_ss_j[nmer_frust_size] = '\0';       // and null terminate it so that it can be printed properly
    nmer_ss_k = new char[nmer_frust_size+1]; // extend the array
    nmer_ss_k[nmer_frust_size] = '\0';       // and null terminate it so that it can be printed properly
    nmer_frust_output_file = fopen("nmer_frustration.dat","w");
    nmer_frust_vmd_script = fopen("nmer_frustration.tcl","w");
    if (strcmp(nmer_frust_mode, "pairwise")==0) {
      fprintf(nmer_frust_output_file,"# i j ncontacts a_i a_j native_energy <decoy_energies> std(decoy_energies) f_ij\n");
    }
    else if (strcmp(nmer_frust_mode, "singlenmer")==0) {
      fprintf(nmer_frust_output_file,"# i a_i native_energy <decoy_energies> std(decoy_energies) f_ij\n");
    }
    if(nmer_frust_trap_flag) {
      nmer_frust_trap_file = fopen("nmer_traps.dat", "w");
      fprintf(nmer_frust_trap_file,"# i a_i ss_i j a_j ss_j threshold_energy k a_k ss_k direction trap_energy\n");
    }
  }

  // Selection temperature file
  if (selection_temperature_flag) {
    if (selection_temperature_output_interaction_energies_flag) {
      selection_temperature_file = fopen(selection_temperature_file_name, "w");
    }
    if (selection_temperature_evaluate_sequence_energies_flag) {
      selection_temperature_sequence_energies_output_file = fopen(selection_temperature_sequence_energies_output_file_name, "w");
      fprintf(selection_temperature_file, "# i j a_i a_j rij rho_i rho_j water burial_i burial_j\n");
      // read in sequences in selection temperature sequences file
      char temp_sequence[1000];
      ifstream selection_temperature_sequences_file(selection_temperature_sequences_file_name);
      selection_temperature_sequences_file >> num_selection_temperature_sequences;
      selection_temperature_sequences = new char*[num_selection_temperature_sequences];
      for (int i=0;i<num_selection_temperature_sequences;i++) {
	selection_temperature_sequences[i] = new char[n];
      }
      for(int i_sequence = 0; i_sequence < num_selection_temperature_sequences; i_sequence++) {
	selection_temperature_sequences_file >> temp_sequence;
	strcpy(selection_temperature_sequences[i_sequence],temp_sequence);
      }
      selection_temperature_sequences_file.close();

      // read in residues in selection temperature residues file
      int temp_res_index;
      ifstream selection_temperature_residues_file(selection_temperature_residues_file_name);
      selection_temperature_residues_file >> num_selection_temperature_residues;
      selection_temperature_residues = new int[num_selection_temperature_residues];
      for(int i=0; i<num_selection_temperature_residues;i++) {
	selection_temperature_residues_file >> temp_res_index;
	selection_temperature_residues[i] = temp_res_index;
      }
      selection_temperature_residues_file.close();
    }
    if (selection_temperature_output_contact_list_flag) {
      selection_temperature_contact_list_file = fopen(selection_temperature_output_contact_list_file_name, "w");
    }
  }

  if (monte_carlo_seq_opt_flag) {
    mcso_seq_output_file = fopen(mcso_seq_output_file_name, "w");
    mcso_energy_output_file = fopen(mcso_energy_output_file_name, "w");
  }
  
  // if optimization_flag is on, perform appropriate initializations
  if(optimization_flag) {
    optimization_file = fopen("optimization_energies.dat","w");    
    
    native_optimization_file = fopen("native_optimization_energies.dat","w");
    
    optimization_norm_file = fopen("optimization_norms.dat","w");
    native_optimization_norm_file = fopen("native_optimization_norms.dat","w");
  }
  
  if (burial_optimization_flag) { 
    burial_optimization_file = fopen("burial_optimization_energies.dat","w");
    native_burial_optimization_file = fopen("native_burial_optimization_energies.dat","w");	
    burial_optimization_norm_file = fopen("burial_optimization_norm.dat","w");
  }

  if (debyehuckel_optimization_flag) { 
    debyehuckel_optimization_file = fopen("debyehuckel_optimization_energies.dat","w");
    debyehuckel_native_optimization_file = fopen("debyehuckel_native_optimization_energies.dat","w");	
    debyehuckel_optimization_norm_file = fopen("debyehuckel_optimization_norm.dat","w");
    debyehuckel_native_optimization_norm_file = fopen("debyehuckel_native_optimization_norm.dat","w");
  }

  // if optimization_flag is on, perform appropriate initializations
  if(average_sequence_optimization_flag) {
    average_sequence_optimization_file = fopen("average_sequence_optimization_energies.dat","w");    
    average_sequence_optimization_norm_file = fopen("average_sequence_optimization_norms.dat","w");
    ifstream in_average_sequence(average_sequence_input_file_name);
    for(i=0;i<n;i++) {
      for(j=0;j<20;j++) {
	in_average_sequence >> average_sequence[i][j];
      }
    }
    in_average_sequence.close();
  }

  // if Mutate_Sequence flag is on, perform the appropriate initializations
  if (mutate_sequence_flag) {
    // read in sequences in selection temperature sequences file
    char temp_sequence[1000];
    ifstream mutate_sequence_sequences_file(mutate_sequence_sequences_file_name);
    mutate_sequence_sequences_file >> mutate_sequence_number_of_sequences;
    mutate_sequence_sequences = new char*[mutate_sequence_number_of_sequences];
    for (int i=0;i<mutate_sequence_number_of_sequences;i++) {
      mutate_sequence_sequences[i] = new char[n];
    }
    for(int i_sequence = 0; i_sequence < mutate_sequence_number_of_sequences; i_sequence++) {
      mutate_sequence_sequences_file >> temp_sequence;
      strcpy(mutate_sequence_sequences[i_sequence],temp_sequence);
    }
    mutate_sequence_sequences_file.close();
    
    mutate_sequence_sequence_index = 0;
  }

  
  // Allocate the table
  if (frag_mem_tb_flag) {
    if (fm_gamma->maxSep()!=-1)
      tb_nbrs = fm_gamma->maxSep()-fm_gamma->minSep()+1;
    else
      tb_nbrs = n - fm_gamma->minSep();
		
    fm_table = new TBV*[4*n*tb_nbrs];
		
    for (i=0; i<4*n*tb_nbrs; ++i) {
      fm_table[i] = NULL;
    }
		
    //		Construct_Computational_Arrays();
    if (comm->me==0) print_log("Computing FM table...\n");
    compute_fragment_memory_table();
  }

  // If using Debye_Huckel potential, read charges from file
  // Skip if DebyeHuckel optimization is on, because it uses a residue type based potential
  // instead of residue index based potential (so that sequence shuffling can be used)
  if (huckel_flag && !debyehuckel_optimization_flag) {
    int residue_number, total_residues; 
    double charge_value;
    double total_charge =0;
    ifstream input_charge("charge_on_residues.dat");
    if (!input_charge) error->all(FLERR,"File charge_on_residues.dat doesn't exist");
    input_charge >> total_residues;
    
    //fprintf(screen, "check charge data \n");
    fprintf(screen, "Number of Charge input = %5d \n", total_residues);
    for(int ires = 0; ires<total_residues; ires++)
      {
	input_charge >> residue_number >> charge_value;
	int res_min_one = residue_number -1;
	charge_on_residue[res_min_one] = charge_value;
	total_charge = total_charge + charge_value;
	//fprintf(screen, "residue=%5d, charge on residue =%8.6f\n", residue_number, charge_value);
	//fprintf(screen, "residue=%5d, charge on residue =%8.6f\n", res_min_one, charge_on_residue[res_min_one]);
      }
    input_charge.close();
    fprintf(screen, "Total Charge on the System = %8.4f\n", total_charge ); 
  }

  if (output_per_residue_contacts_flag) {
    output_per_residue_contacts_file = fopen(output_per_residue_contacts_file_name, "w");
    if (strcmp(output_per_residue_contacts_mode, "native")==0) {
      output_per_residue_contacts_structure = new Fragment_Memory(0, 0, n, 1.0, output_per_residue_contacts_native_structure_file_name);
      if (output_per_residue_contacts_structure->error==output_per_residue_contacts_structure->ERR_FILE) error->all(FLERR,"Cannot read file for Output_Per_Residue_Contacts block");
      if (output_per_residue_contacts_structure->error==output_per_residue_contacts_structure->ERR_ATOM_COUNT) error->all(FLERR,"Output_Per_Residue_Contacts: Wrong atom count in Output_Per_Residue_Contacts block structure file");
      if (output_per_residue_contacts_structure->error==output_per_residue_contacts_structure->ERR_RES) error->all(FLERR,"Output_Per_Residue_Contacts: Unknown residue in Output_Per_Residue_Contacts block structure file");
    }
  }
  
  sStep=0, eStep=0;
  ifstream in_rs("record_steps");
  in_rs >> sStep >> eStep;
  in_rs.close();

}

void FixBackbone::final_log_output()
{
  double time, tmp;
  char txt_timer[][11] = {"Chain", "Shake", "Chi", "Rama", "Vexcluded", "DSSP", "PAP", "Water", "Burial", "Helix", "AHM-Go", "Frag_Mem", "Vec_FM", "Membrane", "SSB", "DH"};
  int me,nprocs;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  fprintf(dout, "\n");
  for (int i=0;i<TIME_N;++i) {
    time = ctime[i];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      fprintf(dout, "%s time = %g\n", txt_timer[i], time); 
    }
  }
  fprintf(dout, "\n");
}

/* ---------------------------------------------------------------------- */
FixBackbone::~FixBackbone()
{
  final_log_output();

  if (allocated) {
    for (int i=0;i<n;i++) {
      delete [] xca[i];
      delete [] xcb[i];
      delete [] xo[i];
      delete [] xn[i];
      delete [] xcp[i];
      delete [] xh[i];
    }

    for (int i=0;i<12;i++) delete [] aps[i];

    delete [] alpha_carbons;
    delete [] beta_atoms;
    delete [] oxygens;
    delete [] res_no;
    delete [] res_no_l;
    delete [] res_info;
    delete [] chain_no;
    delete [] xca;
    delete [] xcb;
    delete [] xo;
    delete [] xn;
    delete [] xcp;
    delete [] xh;

    if (p_ap_flag) {
      delete p_ap;
    }
    delete R;
		
    if (amh_go_flag) {
      for (int i=0;i<3*n;i++) {
	delete [] amh_go_force[i];
      }
      
      delete [] amh_go_force;
      delete [] amh_go_force_map;
			
      delete m_amh_go;
      delete amh_go_gamma;
      if (frustration_censoring_flag == 1) {
	for (int i=0;i<n;i++) {
	  delete [] frustration_censoring_map[i];
	}
	delete [] frustration_censoring_map;
      }
    }
    
    if (frag_mem_flag || frag_mem_tb_flag) {
      delete fm_gamma;
			
      for (int i=0;i<n_frag_mems;i++) delete frag_mems[i];
      if (n_frag_mems>0) memory->sfree(frag_mems);
			
      for (int i=0;i<n;++i) memory->sfree(frag_mem_map[i]);
      delete [] frag_mem_map;
      delete [] ilen_fm_map;
    }
  }
	
  if (frag_mem_tb_flag) {
    for (int i=0; i<4*n*tb_nbrs; ++i) {
      if (fm_table[i])
	delete [] fm_table[i];
    }
    delete [] fm_table;
	
    /*		for (int i=0; i<n; ++i) {
		for (int j=0; j<tb_nbrs; ++j) {
		for (int k=0; j<4; ++k) {
		delete [] fm_table[i][j][k];
		}
		delete [] fm_table[i][j];
		}
		delete [] fm_table[i];
		}		
		delete [] fm_table;*/
  }

  // if the fragment frustratometer was on, close the fragment frustration file
  if (frag_frust_flag) {
    fclose(fragment_frustration_file);
    fclose(fragment_frustration_gap_file);
    fclose(fragment_frustration_variance_file);
    fclose(fragment_frustration_decoy_data);
    fclose(fragment_frustration_native_data);
    for (int i=0;i<n;++i) memory->sfree(decoy_mem_map[i]);
    delete [] decoy_mem_map;
    delete [] ilen_decoy_map;
    for (int i=0;i<n_decoy_mems;i++) delete decoy_mems[i];
    if (n_decoy_mems>0) memory->sfree(decoy_mems);
    if (frag_frust_read_flag) {
      delete [] frag_frust_read_mean;
      delete [] frag_frust_read_variance;
    }
  }

  // if the tertiary frustratometer was on, write the end of the vmd script and close the files
  if (tert_frust_flag) {
    fclose(tert_frust_output_file);
    fclose(tert_frust_vmd_script);
  }

  // if the nmer frustratometer was on, write the end of the vmd script and close the files
  if (nmer_frust_flag) {
    fclose(nmer_frust_output_file);
    fclose(nmer_frust_vmd_script);
    if(nmer_frust_trap_flag) {
      fclose(nmer_frust_trap_file);
    }
  }

  // if the selection temperature was being output, close the file
  if (selection_temperature_flag) {
    if (selection_temperature_output_interaction_energies_flag) {
      fclose(selection_temperature_file);
    }
    if (selection_temperature_evaluate_sequence_energies_flag) {
      fclose(selection_temperature_sequence_energies_output_file);
    }
    if (selection_temperature_output_contact_list_flag) {
      fclose(selection_temperature_contact_list_file);
    }
  }
  
  // if the mcso sequences and energies were being output, close the files
  if (monte_carlo_seq_opt_flag) {
    fclose(mcso_seq_output_file);
    fclose(mcso_energy_output_file);
  }

  // if the optimization block was on, close the files
  if (optimization_flag) {
    fclose(optimization_file);
    fclose(native_optimization_file);	
    fclose(optimization_norm_file);
    fclose(native_optimization_norm_file);
  }
  if (burial_optimization_flag) {
    fclose(burial_optimization_file);
    fclose(native_burial_optimization_file);
    fclose(burial_optimization_norm_file);
  }
  if (debyehuckel_optimization_flag) {
    fclose(debyehuckel_optimization_file);
    fclose(debyehuckel_native_optimization_file);
    fclose(debyehuckel_optimization_norm_file);
    fclose(debyehuckel_native_optimization_norm_file);
  }
  // if the Average Sequence optimization block was on, close the files
  if (average_sequence_optimization_flag) {
    fclose(average_sequence_optimization_file);
    fclose(average_sequence_optimization_norm_file);
  }

  if (huckel_flag) {
    delete[] charge_on_residue;
  }

  if (output_per_residue_contacts_flag) {
    fclose(output_per_residue_contacts_file);
  }
  
  fclose(efile);
}

void FixBackbone::allocate()
{
  int i, j, k;

  alpha_carbons = new int[n];
  beta_atoms = new int[n];
  oxygens = new int[n];
  res_no = new int[n];
  res_no_l = new int[n];
  for (i=0; i<n; ++i){
	res_no_l[i]=-1;
  }
  res_info = new int[n];
  chain_no = new int[n];
	
  xca = new double*[n];
  xcb = new double*[n];
  xo = new double*[n];
  xn = new double*[n];
  xcp = new double*[n];
  xh = new double*[n];

  if (huckel_flag) {
    charge_on_residue = new double[n];
  }

  if (water_flag) {
    water_par = WPV(water_kappa, water_kappa_sigma, treshold, n_wells, well_flag, well_r_min, well_r_max);
    well = new cWell<double, FixBackbone>(n, n, n_wells, water_par, &ntimestep, this);
  }

  if (helix_flag) {
    helix_par = WPV(helix_kappa, helix_kappa_sigma, helix_treshold, n_helix_wells, helix_well_flag, helix_well_r_min, helix_well_r_max);
    helix_well = new cWell<double, FixBackbone>(n, n, n_helix_wells, helix_par, &ntimestep, this);
  }

  if (p_ap_flag) {
    p_ap = new cP_AP<double, FixBackbone>(n, n, &ntimestep, this);
  }

  R = new cR<double, FixBackbone>(n, n, &ntimestep, this);

  for (i = 0; i < n; ++i) {
    // Ca, Cb and O coordinates
    xca[i] = new double [3];
    xcb[i] = new double [3];
    xo[i] = new double [3];
		
    // Nitrogen and C prime coordinates
    xn[i] = new double [3];
    xcp[i] = new double [3];
    xh[i] = new double [3];
    
    if (huckel_flag) {
      charge_on_residue[i] = 0.0;
    }
  }

  for (i = 0; i < 12; ++i) {
    aps[i] = new double[n];
  }

  xn[0][0] = 0;
  xn[0][1] = 0;
  xn[0][2] = 0;
  xcp[n-1][0] = 0; 
  xcp[n-1][1] = 0; 
  xcp[n-1][2] = 0; 
  xh[0][0] = 0;
  xh[0][1] = 0;
  xh[0][2] = 0;
	
  if (amh_go_flag) {
    amh_go_force = new double*[3*n];
    amh_go_force_map = new int[3*n];
    for (int i=0;i<3*n;i++) {
      amh_go_force[i] = new double[3];
    }
    amh_go_norm = new double[nch];

    // if frustration censoring is on, allocate the frustration_censoring_map table
    if (frustration_censoring_flag == 1){
      frustration_censoring_map = new int*[n];
      for (int i=0;i<n;i++) {
	frustration_censoring_map[i] = new int[n];
	for (int j=0;j<n;j++) {
	  frustration_censoring_map[i][j] = 0;
	}
      }
    }
  }
  
  if (average_sequence_optimization_flag) {
    average_sequence = new double*[n];
    for (int i=0;i<n;i++) {
      average_sequence[i] = new double[20];
      for (int j=0;j<20;j++) {
	average_sequence[i][j] = 0.0;
      }
    }
  }
	
  allocated = true;
}

/* ---------------------------------------------------------------------- */
inline bool FixBackbone::isFirst(int index)
{
  if (ch_pos[chain_no[index]-1]==res_no[index]) return true;

  return false;
}

inline bool FixBackbone::isLast(int index)
{
  int ch_no = chain_no[index]-1;
  if (ch_pos[ch_no] + ch_len[ch_no]-1==res_no[index]) return true;

  return false;
}

int FixBackbone::Tag(int index) {
  if (index==-1) return -1;
  return atom->tag[index];
}

inline void FixBackbone::Construct_Computational_Arrays()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int *mol_tag = atom->molecule;
  int *res_tag = atom->residue;

  int i;
  for (i=0; i<n; ++i){
	res_no_l[i]=-1;
  }

  // Creating index arrays for Alpha_Carbons, Beta_Atoms and Oxygens
  nn = 0;
  int last = 0;
  for (i = 0; i < n; ++i) {
    int min[3] = {-1, -1, -1}, jm[3] = {-1, -1, -1}, amin = -1;
    for (int j = 0; j < nall; ++j) {
      if (i==0 && res_tag[j]<=0 && (mask[j] & groupbit || mask[j] & group2bit || mask[j] & group3bit) )
	error->all(FLERR,"Molecular tag must be positive in fix backbone");
			
      if ( (mask[j] & groupbit) && res_tag[j]>last ) {
	if (res_tag[j]<min[0] || min[0]==-1) {
	  min[0] = res_tag[j];
	  jm[0] = j;
	}
      }
      if ( (mask[j] & group2bit) && res_tag[j]>last ) {
	if (res_tag[j]<min[1] || min[1]==-1) {
	  min[1] = res_tag[j];
	  jm[1] = j;
	}
      }
      if ( (mask[j] & group3bit) && res_tag[j]>last ) {
	if (res_tag[j]<min[2] || min[2]==-1) {
	  min[2] = res_tag[j];
	  jm[2] = j;
	}
      }
    }
		
    amin = MIN(min[0], MIN(min[1], min[2]));
    if (amin==-1) break;

    if (min[0]!=amin) jm[0] = -1;
    if (min[1]!=amin) jm[1] = -1;
    if (min[2]!=amin) jm[2] = -1;

    alpha_carbons[nn] = jm[0];
    beta_atoms[nn] = jm[1];
    oxygens[nn] = jm[2];
    res_no[nn] = amin; //amin range [1,n]
    res_no_l[res_no[nn]-1]=nn; //local i=res_no_l[i_resno];
    last = amin;
    nn++;
  }

  for (i = 0; i < nn; ++i) {
    chain_no[i] = -1;
	
    // Checking sequance and marking residues
    if (alpha_carbons[i]!=-1) {
			
      // Making sure chain tags match for same residue atoms
      if ( (beta_atoms[i]!=-1 && mol_tag[alpha_carbons[i]]!=mol_tag[beta_atoms[i]]) ||
	   (oxygens[i]!=-1 && mol_tag[alpha_carbons[i]]!=mol_tag[oxygens[i]]) ) {
	error->all(FLERR,"Atoms in a residue have different chain tag");
      }			
      chain_no[i] = mol_tag[alpha_carbons[i]];
			
      if (chain_no[i]<=0 || chain_no[i]>nch)
	error->all(FLERR,"Chain tag is out of range");
			
      // Checking for correct residue numbering
      if ( res_no[i]<ch_pos[chain_no[i]-1] || res_no[i]>ch_pos[chain_no[i]-1]+ch_len[chain_no[i]-1]-1 )
	error->all(FLERR,"Residue tag is out of range");

      if (alpha_carbons[i]<nlocal) {
	//if (beta_atoms[i]==-1 || oxygens[i]==-1) {
	//  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 001)");
	//}
	//if ( !isFirst(i) && (i==0 || res_info[i-1]==OFF) ) {
	//  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 002)");
	//}
	res_info[i] = LOCAL;
      } else {
	if ( i>0 && !isFirst(i) && res_info[i-1]==LOCAL ) res_info[i] = GHOST;
	else if (i<nn-1 && !isLast(i) && alpha_carbons[i+1]<nlocal && alpha_carbons[i+1]!=-1) {
	  //if (oxygens[i]==-1) {
	  //  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 003)");
	  //}
	  res_info[i] = GHOST;
	} else if (oxygens[i]==-1 || beta_atoms[i]==-1) {
	  res_info[i] = OFF;
	} else res_info[i] = GHOST;
      }
			
    } else res_info[i] = OFF;
		
    //if (i>0 && res_info[i-1]==LOCAL && res_info[i]==OFF) {
    //  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 004)");
    //}
  }
	
  /*	if (ntimestep==0) {
	for (i = 0; i < nn; ++i) {
	fprintf(dout, "%d %d %d %d %d %d\n", i, res_no[i], res_info[i], alpha_carbons[i], beta_atoms[i], oxygens[i]);
	}
	}*/
}

/*inline void FixBackbone::Construct_Computational_Arrays()
  {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  int *mol_tag = atom->molecule;

  int i;

  // Creating index arrays for Alpha_Carbons, Beta_Atoms and Oxygens
  nn = 0;
  int last = 0;
  for (i = 0; i < n; ++i) {
  int min[3] = {-1, -1, -1}, jm[3] = {-1, -1, -1}, amin = -1;
  for (int j = 0; j < nall; ++j) {
  if (i==0 && mol_tag[j]<=0)
  error->all(FLERR,"Molecular tag must be positive in fix backbone");
			
  if ( (mask[j] & groupbit) && mol_tag[j]>last ) {
  if (mol_tag[j]<min[0] || min[0]==-1) {
  min[0] = mol_tag[j];
  jm[0] = j;
  }
  }
  if ( (mask[j] & group2bit) && mol_tag[j]>last ) {
  if (mol_tag[j]<min[1] || min[1]==-1) {
  min[1] = mol_tag[j];
  jm[1] = j;
  }
  }
  if ( (mask[j] & group3bit) && mol_tag[j]>last ) {
  if (mol_tag[j]<min[2] || min[2]==-1) {
  min[2] = mol_tag[j];
  jm[2] = j;
  }
  }
  }
		
  amin = MIN(min[0], MIN(min[1], min[2]));
  if (amin==-1) break;

  if (min[0]!=amin) jm[0] = -1;
  if (min[1]!=amin) jm[1] = -1;
  if (min[2]!=amin) jm[2] = -1;

  alpha_carbons[nn] = jm[0];
  beta_atoms[nn] = jm[1];
  oxygens[nn] = jm[2];
  res_no[nn] = amin;
  last = amin;
  nn++;
  }

  for (i = 0; i < nn; ++i) {
  // Checking sequance and marking residues
  if (alpha_carbons[i]!=-1) {
  if (alpha_carbons[i]<nlocal) {
  if (beta_atoms[i]==-1 || oxygens[i]==-1) {
  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 001)");
  }
  if ( !isFirst(i) && (i==0 || res_info[i-1]==OFF) ) {
  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 002)");
  }
  res_info[i] = LOCAL;
  } else {
  if ( i>0 && !isFirst(i) && res_info[i-1]==LOCAL ) res_info[i] = GHOST;
  else if (i<nn-1 && !isLast(i) && alpha_carbons[i+1]<nlocal && alpha_carbons[i+1]!=-1) {
  if (oxygens[i]==-1) {
  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 003)");
  }
  res_info[i] = GHOST;
  } else res_info[i] = OFF;
  }
			
  } else res_info[i] = OFF;
		
  if (i>0 && res_info[i-1]==LOCAL && res_info[i]==OFF) {
  error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 004)");
  }
  }
  }*/

/* ---------------------------------------------------------------------- */

int FixBackbone::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= THERMO_ENERGY;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackbone::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
	
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  //  neighbor->requests[irequest]->occasional = 0;

}

/* ---------------------------------------------------------------------- */

void FixBackbone::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBackbone::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    pre_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    pre_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixBackbone::min_setup(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

inline double adotb(double *a, double *b)
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

inline double cross(double *a, double *b, int index)
{
  switch (index) {
  case 0:
    return a[1]*b[2] - a[2]*b[1];
  case 1:
    return a[2]*b[0] - a[0]*b[2];
  case 2:
    return a[0]*b[1] - a[1]*b[0];
  }
	
  return 0;
}

inline double atan2(double y, double x)
{
  if (x==0) {
    if (y>0) return M_PI_2;
    else if (y<0) return -M_PI_2;
    else return NULL;
  } else {
    return atan(y/x) + (x>0 ? 0 : (y>=0 ? M_PI : -M_PI) );
  }
}

inline double FixBackbone::PeriodicityCorrection(double d, int i)
{
  if (!periodicity[i]) return d;
  else return ( d > half_prd[i] ? d - prd[i] : (d < -half_prd[i] ? d + prd[i] : d) );
}

bool FixBackbone::isEmptyString(char *str)
{
  int len = strlen(str);
  
  if (len==0) return true;
  
  for (int i=0;i<len;++i) {
    if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') return false;
  }

  return true;
}

char *FixBackbone::ltrim(char *s) 
{     
  while(isspace(*s)) s++;     
  return s; 
}  

char *FixBackbone::rtrim(char *s) 
{
  char* back;
  int len = strlen(s);

  if(len == 0)
    return(s); 

  back = s + len;     
  while(isspace(*--back));     
  *(back+1) = '\0';     
  return s; 
}  

char *FixBackbone::trim(char *s) 
{     
  return rtrim(ltrim(s));  
} 

Fragment_Memory **FixBackbone::read_mems(char *mems_file, int &n_mems)
{
  int file_state, nstr;
  int tpos, fpos, len;
  double weight;
  char ln[500], *line, *str[10];
  FILE *file;
  Fragment_Memory **mems_array = NULL;
  
  enum File_States{FS_NONE=0, FS_TARGET, FS_MEMS};
  
  file = fopen(mems_file,"r");
  if (!file) error->all(FLERR,"Fragment_Memory: Error opening mem file");
  
  n_mems = 0;
  file_state = FS_NONE;
  while ( fgets ( ln, sizeof ln, file ) != NULL ) {
    line = trim(ln);
    
    if (line[0]=='#') continue;
    if (line[0]=='[') file_state = FS_NONE;
    if (isEmptyString(line)) { file_state = FS_NONE; continue; }
        
    switch (file_state) {
    case FS_MEMS:
      nstr = 0;
      str[nstr] = strtok(line," \t\n");
      while ( str[nstr]!=NULL ) {
        nstr++;
        if (nstr>5) break;
        str[nstr] = strtok(NULL," \t\n");
      }
      if (nstr!=5) error->all(FLERR,"Fragment_Memory: Error reading mem file");
        
      tpos = atoi(str[1])-1;
      fpos = atoi(str[2])-1;
      len = atoi(str[3]);
      weight = atof(str[4]);
        
      n_mems++;
      mems_array = (Fragment_Memory **) memory->srealloc(mems_array,n_mems*sizeof(Fragment_Memory *),"modify:mems_array");
      mems_array[n_mems-1] = new Fragment_Memory(tpos, fpos, len, weight, str[0], vec_frag_mem_flag);
      
      if (mems_array[n_mems-1]->error!=Fragment_Memory::ERR_NONE) {
        if (screen) fprintf(screen, "Error reading %s file!\n", str[0]);
        if (logfile) fprintf(logfile, "Error reading %s file!\n", str[0]);
      } 
      if (mems_array[n_mems-1]->error==Fragment_Memory::ERR_FILE) error->all(FLERR,"Fragment_Memory: Cannot read the file");
      if (mems_array[n_mems-1]->error==Fragment_Memory::ERR_ATOM_COUNT) error->all(FLERR,"Fragment_Memory: Wrong atom count in memory structure file");
      if (mems_array[n_mems-1]->error==Fragment_Memory::ERR_RES) error->all(FLERR,"Fragment_Memory: Unknown residue");
      
      if (mems_array[n_mems-1]->pos+mems_array[n_mems-1]->len>n) {
      	if (screen) fprintf(screen, "Error reading %s file!\n", str[0]);
        if (logfile) fprintf(logfile, "Error reading %s file!\n", str[0]);
        fprintf(stderr, "pos %d len %d n %d\n", mems_array[n_mems-1]->pos, mems_array[n_mems-1]->len, n); 
      	error->all(FLERR,"read_mems: Fragment_Memory: Incorrectly defined memory fragment");
      }
        
      break;
    case FS_NONE:
      if (strcmp(line, "[Target]")==0)
        file_state = FS_TARGET;
      else if (strcmp(line, "[Memories]")==0)
        file_state = FS_MEMS;
      break;
    }
  }
  
  fclose(file);
  
  return mems_array;
}

void FixBackbone::timerBegin()
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  previous_time = MPI_Wtime();
}

void FixBackbone::timerEnd(int which)
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  double current_time = MPI_Wtime();
  ctime[which] += current_time - previous_time;
  previous_time = current_time;
}

/* ---------------------------------------------------------------------- */

void FixBackbone::compute_chain_potential(int i)
{	
  double dx[3], r, dr, force;
  int i_resno = res_no[i]-1;
  int ip1, im1; 
  // N(i) - Cb(i)
  if (!isFirst(i) && se[i_resno]!='G') {
    im1 = res_no_l[i_resno-1];
    if (im1!=-1 && (res_info[im1]==LOCAL || res_info[im1]==GHOST)){
		dx[0] = xn[i][0] - xcb[i][0];
		dx[1] = xn[i][1] - xcb[i][1];
		dx[2] = xn[i][2] - xcb[i][2];
		r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
		dr = r - r_ncb0;
		force = 2*epsilon*k_chain[0]*dr/r;
		
		energy[ET_CHAIN] += epsilon*k_chain[0]*dr*dr;
	
		f[alpha_carbons[im1]][0] -= an*dx[0]*force;
		f[alpha_carbons[im1]][1] -= an*dx[1]*force;
		f[alpha_carbons[im1]][2] -= an*dx[2]*force;
				
		f[oxygens[im1]][0] -= cn*dx[0]*force;
		f[oxygens[im1]][1] -= cn*dx[1]*force;
		f[oxygens[im1]][2] -= cn*dx[2]*force;
		
		f[alpha_carbons[i]][0] -= bn*dx[0]*force;
		f[alpha_carbons[i]][1] -= bn*dx[1]*force;
		f[alpha_carbons[i]][2] -= bn*dx[2]*force;	
		
		f[beta_atoms[i]][0] -= -dx[0]*force;
		f[beta_atoms[i]][1] -= -dx[1]*force;
		f[beta_atoms[i]][2] -= -dx[2]*force;
    }
  }
	
  // Cp(i) - Cb(i)
  if (!isLast(i) && se[i_resno]!='G') {
	ip1 = res_no_l[i_resno+1];
	if (ip1!=-1 && (res_info[ip1]==LOCAL || res_info[ip1]==GHOST)){
		dx[0] = xcp[i][0] - xcb[i][0];
		dx[1] = xcp[i][1] - xcb[i][1];
		dx[2] = xcp[i][2] - xcb[i][2];
		r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
		dr = r - r_cpcb0;
		force = 2*epsilon*k_chain[1]*dr/r;
		
		energy[ET_CHAIN] += epsilon*k_chain[1]*dr*dr;  
		
		f[alpha_carbons[ip1]][0] -= bp*dx[0]*force;
		f[alpha_carbons[ip1]][1] -= bp*dx[1]*force;
		f[alpha_carbons[ip1]][2] -= bp*dx[2]*force;
		
		f[alpha_carbons[i]][0] -= ap*dx[0]*force;
		f[alpha_carbons[i]][1] -= ap*dx[1]*force;
		f[alpha_carbons[i]][2] -= ap*dx[2]*force;
		
		f[oxygens[i]][0] -= cp*dx[0]*force;
		f[oxygens[i]][1] -= cp*dx[1]*force;
		f[oxygens[i]][2] -= cp*dx[2]*force;
		
		f[beta_atoms[i]][0] -= -dx[0]*force;
		f[beta_atoms[i]][1] -= -dx[1]*force;
		f[beta_atoms[i]][2] -= -dx[2]*force;
	}
  }


  // N(i) - Cp(i)
  if (!isFirst(i) && !isLast(i)) {
	im1 = res_no_l[i_resno-1];
	ip1 = res_no_l[i_resno+1];
	if( im1!=-1 && ip1!=-1 && (res_info[im1]==LOCAL || res_info[im1]==GHOST) && (res_info[ip1]==LOCAL || res_info[ip1]==GHOST)){
		dx[0] = xn[i][0] - xcp[i][0];
		dx[1] = xn[i][1] - xcp[i][1];
		dx[2] = xn[i][2] - xcp[i][2];
		r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
		dr = r - r_ncp0;
		force = 2*epsilon*k_chain[2]*dr/r;
		
		energy[ET_CHAIN] += epsilon*k_chain[2]*dr*dr;
		
		f[alpha_carbons[im1]][0] -= an*dx[0]*force;
		f[alpha_carbons[im1]][1] -= an*dx[1]*force;
		f[alpha_carbons[im1]][2] -= an*dx[2]*force;
				
		f[oxygens[im1]][0] -= cn*dx[0]*force;
		f[oxygens[im1]][1] -= cn*dx[1]*force;
		f[oxygens[im1]][2] -= cn*dx[2]*force;
			
		f[alpha_carbons[ip1]][0] -= -bp*dx[0]*force;
		f[alpha_carbons[ip1]][1] -= -bp*dx[1]*force;
		f[alpha_carbons[ip1]][2] -= -bp*dx[2]*force;
		
		f[alpha_carbons[i]][0] -= (bn-ap)*dx[0]*force;
		f[alpha_carbons[i]][1] -= (bn-ap)*dx[1]*force;
		f[alpha_carbons[i]][2] -= (bn-ap)*dx[2]*force;
		
		f[oxygens[i]][0] -= -cp*dx[0]*force;
		f[oxygens[i]][1] -= -cp*dx[1]*force;
		f[oxygens[i]][2] -= -cp*dx[2]*force;
	}
  }
}

void FixBackbone::compute_shake(int i)
{	
  double dx[3], r, dr, force;
	
  // r_sh1 = r_Ca(i) - rCa(i+1)
  // r_sh2 = r_Ca(i) - r_O(i)
  // r_sh3 = r_O(i) - r_Ca(i+1)
	
  // Ca(i) - Ca(i+1)
  if (!isLast(i)) {
    dx[0] = xca[i][0] - xca[i+1][0];
    dx[1] = xca[i][1] - xca[i+1][1];
    dx[2] = xca[i][2] - xca[i+1][2];
    r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    dr = r - r_sh1;
    force = 2*epsilon*k_shake*dr/r;
		
    energy[ET_SHAKE] += epsilon*k_shake*dr*dr;
		
    f[alpha_carbons[i]][0] -= dx[0]*force;
    f[alpha_carbons[i]][1] -= dx[1]*force;
    f[alpha_carbons[i]][2] -= dx[2]*force;
	
    f[alpha_carbons[i+1]][0] -= -dx[0]*force;
    f[alpha_carbons[i+1]][1] -= -dx[1]*force;
    f[alpha_carbons[i+1]][2] -= -dx[2]*force;
  }
	
  // Ca(i) - O(i)
  dx[0] = xca[i][0] - xo[i][0];
  dx[1] = xca[i][1] - xo[i][1];
  dx[2] = xca[i][2] - xo[i][2];
  r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  dr = r - r_sh2;
  force = 2*epsilon*k_shake*dr/r;
	
  energy[ET_SHAKE] += epsilon*k_shake*dr*dr;
	
  f[alpha_carbons[i]][0] -= dx[0]*force;
  f[alpha_carbons[i]][1] -= dx[1]*force;
  f[alpha_carbons[i]][2] -= dx[2]*force;
	
  f[oxygens[i]][0] -= -dx[0]*force;
  f[oxygens[i]][1] -= -dx[1]*force;
  f[oxygens[i]][2] -= -dx[2]*force;
	
  // O(i) - Ca(i+1)
  if (!isLast(i)) {
    dx[0] = xo[i][0] - xca[i+1][0];
    dx[1] = xo[i][1] - xca[i+1][1];
    dx[2] = xo[i][2] - xca[i+1][2];
    r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    dr = r - r_sh3;
    force = 2*epsilon*k_shake*dr/r;
		
    energy[ET_SHAKE] += epsilon*k_shake*dr*dr;
		
    f[oxygens[i]][0] -= dx[0]*force;
    f[oxygens[i]][1] -= dx[1]*force;
    f[oxygens[i]][2] -= dx[2]*force;
	
    f[alpha_carbons[i+1]][0] -= -dx[0]*force;
    f[alpha_carbons[i+1]][1] -= -dx[1]*force;
    f[alpha_carbons[i+1]][2] -= -dx[2]*force;
  }
}

void FixBackbone::compute_chi_potential(int i)
{
  double dx[3], r, dr, force;
  double a[3], b[3], c[3], arvsq, brvsq, crvsq;
  double axb[3], cxa[3], bxc[3], aprl[3], bprl[3], cprl[3];
  double norm, chi, dchi;
  int i_resno = res_no[i]-1;
  int im1, ip1;
	
  a[0] = xcp[i][0] - xca[i][0];
  a[1] = xcp[i][1] - xca[i][1];
  a[2] = xcp[i][2] - xca[i][2];

  b[0] = xca[i][0] - xn[i][0];
  b[1] = xca[i][1] - xn[i][1];
  b[2] = xca[i][2] - xn[i][2];

  c[0] = xca[i][0] - xcb[i][0];
  c[1] = xca[i][1] - xcb[i][1];
  c[2] = xca[i][2] - xcb[i][2];

  arvsq = 1.0/adotb(a, a);
  brvsq = 1.0/adotb(b, b);
  crvsq = 1.0/adotb(c, c);

  norm = sqrt(arvsq*brvsq*crvsq);

  axb[0] = cross(a, b, 0);
  axb[1] = cross(a, b, 1);
  axb[2] = cross(a, b, 2);

  cxa[0] = cross(c, a, 0);
  cxa[1] = cross(c, a, 1);
  cxa[2] = cross(c, a, 2);

  bxc[0] =  cross(b, c, 0);
  bxc[1] =  cross(b, c, 1);
  bxc[2] =  cross(b, c, 2);

  chi = adotb(axb, c)*norm;

  // partial derivitives
  aprl[0] = norm*bxc[0] - arvsq*chi*a[0];
  aprl[1] = norm*bxc[1] - arvsq*chi*a[1];
  aprl[2] = norm*bxc[2] - arvsq*chi*a[2];

  bprl[0] = norm*cxa[0] - brvsq*chi*b[0];
  bprl[1] = norm*cxa[1] - brvsq*chi*b[1];
  bprl[2] = norm*cxa[2] - brvsq*chi*b[2];

  cprl[0] = norm*axb[0] - crvsq*chi*c[0];
  cprl[1] = norm*axb[1] - crvsq*chi*c[1];
  cprl[2] = norm*axb[2] - crvsq*chi*c[2];

  dchi = chi - chi0;
  force = 2*epsilon*k_chi*dchi;
	
  energy[ET_CHI] += epsilon*k_chi*dchi*dchi;
	
  if (!isFirst(i)) {
    im1 = res_no_l[i_resno-1];
    if (im1!=-1 && (res_info[im1]==LOCAL || res_info[im1]==GHOST)){ 
		f[alpha_carbons[im1]][0] -= -an*bprl[0]*force;
		f[alpha_carbons[im1]][1] -= -an*bprl[1]*force;
		f[alpha_carbons[im1]][2] -= -an*bprl[2]*force;
			
		f[oxygens[im1]][0] -= -cn*bprl[0]*force;
		f[oxygens[im1]][1] -= -cn*bprl[1]*force;
		f[oxygens[im1]][2] -= -cn*bprl[2]*force;
    }
  }
	
  if (!isLast(i)) {
    ip1 = res_no_l[i_resno+1];
    if (ip1!=-1 && (res_info[ip1]==LOCAL || res_info[ip1]==GHOST)){
		f[alpha_carbons[ip1]][0] -= bp*aprl[0]*force;
		f[alpha_carbons[ip1]][1] -= bp*aprl[1]*force;
		f[alpha_carbons[ip1]][2] -= bp*aprl[2]*force;
    }
  }

  f[alpha_carbons[i]][0] -= (cprl[0] + (1-bn)*bprl[0] + (ap-1)*aprl[0])*force;
  f[alpha_carbons[i]][1] -= (cprl[1] + (1-bn)*bprl[1] + (ap-1)*aprl[1])*force;
  f[alpha_carbons[i]][2] -= (cprl[2] + (1-bn)*bprl[2] + (ap-1)*aprl[2])*force;

  f[oxygens[i]][0] -= cp*aprl[0]*force;
  f[oxygens[i]][1] -= cp*aprl[1]*force;
  f[oxygens[i]][2] -= cp*aprl[2]*force;

  f[beta_atoms[i]][0] -= -cprl[0]*force;
  f[beta_atoms[i]][1] -= -cprl[1]*force;
  f[beta_atoms[i]][2] -= -cprl[2]*force;
}

void FixBackbone::calcDihedralAndSlopes(int i, double& angle, int iAng)
{
  double a[3], b[3], c[3];
  double bxa[3], cxa[3], cxb[3];
  double adb, bdc, adc, b2, bm, cdbxa;
  double X, Y, X2Y2;
  double dAngle_y, dAngle_x;
  double h1, h2, h3;
  int i_resno = res_no[i]-1;	
  int im1 = res_no_l[i_resno-1];
  int ip1 = res_no_l[i_resno+1];

  if (im1!=-1 && ip1!=-1){
   if (iAng==PHI) {
    a[0] = xcp[i][0] - xca[i][0];
    a[1] = xcp[i][1] - xca[i][1];
    a[2] = xcp[i][2] - xca[i][2];
		
    b[0] = xca[i][0] - xn[i][0];
    b[1] = xca[i][1] - xn[i][1];
    b[2] = xca[i][2] - xn[i][2];
		
    c[0] = xn[i][0] - xcp[im1][0];
    c[1] = xn[i][1] - xcp[im1][1];
    c[2] = xn[i][2] - xcp[im1][2];
   } else {
    a[0] = xn[ip1][0] - xcp[i][0];
    a[1] = xn[ip1][1] - xcp[i][1];
    a[2] = xn[ip1][2] - xcp[i][2];
		
    b[0] = xcp[i][0] - xca[i][0];
    b[1] = xcp[i][1] - xca[i][1];
    b[2] = xcp[i][2] - xca[i][2];
		
    c[0] = xca[i][0] - xn[i][0];
    c[1] = xca[i][1] - xn[i][1];
    c[2] = xca[i][2] - xn[i][2];
   }
  }
  adb = adotb(a, b);
  bdc = adotb(b, c);
  adc = adotb(a, c);
  b2 = adotb(b, b);
  bm = sqrt(b2);
	
  bxa[0] = cross(b, a, 0);
  bxa[1] = cross(b, a, 1);
  bxa[2] = cross(b, a, 2);
	
  cxa[0] = cross(c, a, 0);
  cxa[1] = cross(c, a, 1);
  cxa[2] = cross(c, a, 2);
	
  cxb[0] = cross(c, b, 0);
  cxb[1] = cross(c, b, 1);
  cxb[2] = cross(c, b, 2);
	
  cdbxa = adotb(c, bxa);
	
  Y = bm*cdbxa;
  X = adb*bdc - b2*adc;
	
  angle = atan2(Y, X);
	
  X2Y2 = (X*X + Y*Y);
  dAngle_y = X/X2Y2;	
  dAngle_x = -Y/X2Y2;
	
  for (int l=0;l<3;l++) {
    if (iAng==PHI) {
      y_slope[iAng][CA0][l] = dAngle_y*( -an*b[l]*cdbxa/bm + bm*( (an-ap)*bxa[l] + an*cxa[l] ) );	
      y_slope[iAng][CA1][l] = dAngle_y*( (1-bn)*b[l]*cdbxa/bm + bm*( (bn-bp)*bxa[l] + (ap-1)*cxb[l] - (1-bn)*cxa[l] ) );
      y_slope[iAng][CA2][l] = dAngle_y*bm*bp*cxb[l];
      y_slope[iAng][O0][l] = dAngle_y*( -cn*b[l]*cdbxa/bm + bm*( (cn-cp)*bxa[l] + cn*cxa[l] ) );	
      y_slope[iAng][O1][l] = dAngle_y*bm*cp*cxb[l];
			
      h1 = b[l]*bdc - c[l]*b2;
      h2 = a[l]*bdc - 2*b[l]*adc + c[l]*adb;
      h3 = b[l]*adb - a[l]*b2;
      x_slope[iAng][CA0][l] = dAngle_x*( -an*h2 + (an-ap)*h3 );
      x_slope[iAng][CA1][l] = dAngle_x*( (ap-1)*h1 + (1-bn)*h2 + (bn-bp)*h3 );
      x_slope[iAng][CA2][l] = dAngle_x*bp*h1;
      x_slope[iAng][O0][l] = dAngle_x*( -cn*h2 + (cn-cp)*h3 );
      x_slope[iAng][O1][l] = dAngle_x*cp*h1;
    } else {
      y_slope[iAng][CA0][l] = -dAngle_y*bm*an*bxa[l];
      y_slope[iAng][CA1][l] = dAngle_y*( (ap-1)*b[l]*cdbxa/bm + bm*( (1-bn)*bxa[l] + (an-ap)*cxb[l] - (ap-1)*cxa[l] ) );
      y_slope[iAng][CA2][l] = dAngle_y*( bp*b[l]*cdbxa/bm + bm*( (bn-bp)*cxb[l] - bp*cxa[l] ) );
      y_slope[iAng][O0][l] = -dAngle_y*bm*cn*bxa[l];
      y_slope[iAng][O1][l] = dAngle_y*( cp*b[l]*cdbxa/bm + bm*( (cn-cp)*cxb[l] - cp*cxa[l] ) );
			
      h1 = b[l]*bdc - c[l]*b2;
      h2 = a[l]*bdc - 2*b[l]*adc + c[l]*adb;
      h3 = b[l]*adb - a[l]*b2;
      x_slope[iAng][CA0][l] = -dAngle_x*an*h3;
      x_slope[iAng][CA1][l] = dAngle_x*( (an-ap)*h1 + (ap-1)*h2 + (1-bn)*h3 );
      x_slope[iAng][CA2][l] = dAngle_x*( (bn-bp)*h1 + bp*h2 );
      x_slope[iAng][O0][l] = -dAngle_x*cn*h3;
      x_slope[iAng][O1][l] = dAngle_x*( (cn-cp)*h1 + cp*h2 );
    }
  }
}

void FixBackbone::compute_rama_potential(int i)
{
  double V, phi, psi;
  double force, force1[nAngles];
  int jStart, nEnd;
  int j, ia, l;

  int i_resno   = res_no[i]-1;
  int im1 = res_no_l[i_resno-1];
  int ip1 = res_no_l[i_resno+1];		
	
	
  calcDihedralAndSlopes(i, phi, PHI);
  calcDihedralAndSlopes(i, psi, PSI);
	
  jStart = 0;
  nEnd = n_rama_par;
  if (se[i_resno]=='P' && rama_p_flag) {
    jStart = i_rp;
    nEnd = i_rp + n_rama_p_par;
  }

  for (j=jStart;j<nEnd;j++) {
    V = epsilon*k_rama*w[j]*exp( -sigma[j]*( phiw[j]*pow(cos(phi + phi0[j]) - 1, 2) + psiw[j]*pow(cos(psi + psi0[j]) - 1, 2) ) );

    if (ssweight[j]) {
      if (aps[j][i_resno]==0.0) continue;
      V *= aps[j][i_resno];
    }

    force = 2*V*sigma[j];
    force1[PHI] = force*phiw[j]*(cos(phi + phi0[j]) - 1)*sin(phi + phi0[j]);
    force1[PSI] = force*psiw[j]*(cos(psi + psi0[j]) - 1)*sin(psi + psi0[j]);
			
    energy[ET_RAMA] += -V;
    if (im1!=-1 && ip1!=-1 && (res_info[im1]==LOCAL || res_info[im1]==GHOST) && (res_info[ip1]==LOCAL || res_info[ip1]==GHOST)){
		for (ia=0; ia<nAngles; ia++) {
		  for (l=0; l<3; l++) {				
			f[alpha_carbons[im1]][l] += force1[ia]*(y_slope[ia][CA0][l] + x_slope[ia][CA0][l]);
			f[alpha_carbons[i]][l] += force1[ia]*(y_slope[ia][CA1][l] + x_slope[ia][CA1][l]);
			f[alpha_carbons[ip1]][l] += force1[ia]*(y_slope[ia][CA2][l] + x_slope[ia][CA2][l]);
						
			f[oxygens[im1]][l] += force1[ia]*(y_slope[ia][O0][l] + x_slope[ia][O0][l]);
			f[oxygens[i]][l] += force1[ia]*(y_slope[ia][O1][l] + x_slope[ia][O1][l]);
		  }
		}
    }
  }
} 

void FixBackbone::compute_excluded_volume()
{ 
  double dx[3], r, dr, force;
	
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      // Ca(i) - Cb(j)
      dx[0] = xca[i][0] - xcb[j][0];
      dx[1] = xca[i][1] - xcb[j][1];
      dx[2] = xca[i][2] - xcb[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rC_ex0) {
	dr = r - rC_ex0;
	force = 2*epsilon*k_excluded_C*dr/r;
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_C*dr*dr;
				
	f[alpha_carbons[i]][0] -= dx[0]*force;
	f[alpha_carbons[i]][1] -= dx[1]*force;
	f[alpha_carbons[i]][2] -= dx[2]*force;
				
	f[beta_atoms[j]][0] -= -dx[0]*force;
	f[beta_atoms[j]][1] -= -dx[1]*force;
	f[beta_atoms[j]][2] -= -dx[2]*force;
      }
			
      if (j<=i) continue;

      // Ca(i) - Ca(j)
      dx[0] = xca[i][0] - xca[j][0];
      dx[1] = xca[i][1] - xca[j][1];
      dx[2] = xca[i][2] - xca[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rC_ex0) {
	dr = r - rC_ex0;
	force = 2*epsilon*k_excluded_C*dr/r;
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_C*dr*dr;
				
	f[alpha_carbons[i]][0] -= dx[0]*force;
	f[alpha_carbons[i]][1] -= dx[1]*force;
	f[alpha_carbons[i]][2] -= dx[2]*force;
				
	f[alpha_carbons[j]][0] -= -dx[0]*force;
	f[alpha_carbons[j]][1] -= -dx[1]*force;
	f[alpha_carbons[j]][2] -= -dx[2]*force;
      }

      // Cb(i) - Cb(j)
      dx[0] = xcb[i][0] - xcb[j][0];
      dx[1] = xcb[i][1] - xcb[j][1];
      dx[2] = xcb[i][2] - xcb[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rC_ex0) {
	dr = r - rC_ex0;
	force = 2*epsilon*k_excluded_C*dr/r;
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_C*dr*dr;
				
	f[beta_atoms[i]][0] -= dx[0]*force;
	f[beta_atoms[i]][1] -= dx[1]*force;
	f[beta_atoms[i]][2] -= dx[2]*force;
				
	f[beta_atoms[j]][0] -= -dx[0]*force;
	f[beta_atoms[j]][1] -= -dx[1]*force;
	f[beta_atoms[j]][2] -= -dx[2]*force;
      }
			
      // O(i) - O(j)
      dx[0] = xo[i][0] - xo[j][0];
      dx[1] = xo[i][1] - xo[j][1];
      dx[2] = xo[i][2] - xo[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rO_ex0) {
	dr = r - rO_ex0;
	force = 2*epsilon*k_excluded_O*dr/r;
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_O*dr*dr;
				
	f[oxygens[i]][0] -= dx[0]*force;
	f[oxygens[i]][1] -= dx[1]*force;
	f[oxygens[i]][2] -= dx[2]*force;
				
	f[oxygens[j]][0] -= -dx[0]*force;
	f[oxygens[j]][1] -= -dx[1]*force;
	f[oxygens[j]][2] -= -dx[2]*force;
      }
    }
  }
}

void FixBackbone::compute_p_degree_excluded_volume()
{ 
  int sign = (p%2==0 ? 1 : -1);
  double factorC = sign/pow(rC_ex0, p-2);
  double factorO = sign/pow(rO_ex0, p-2);

  double dx[3], r, dr, force;
	
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      // Ca(i) - Cb(j)
      dx[0] = xca[i][0] - xcb[j][0];
      dx[1] = xca[i][1] - xcb[j][1];
      dx[2] = xca[i][2] - xcb[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rC_ex0) {
	dr = r - rC_ex0;
	force = factorC*p*epsilon*k_excluded_C*pow(dr, p-1)/r;
			
	energy[ET_VEXCLUDED] += factorC*epsilon*k_excluded_C*pow(dr, p);
				
	f[alpha_carbons[i]][0] -= dx[0]*force;
	f[alpha_carbons[i]][1] -= dx[1]*force;
	f[alpha_carbons[i]][2] -= dx[2]*force;
				
	f[beta_atoms[j]][0] -= -dx[0]*force;
	f[beta_atoms[j]][1] -= -dx[1]*force;
	f[beta_atoms[j]][2] -= -dx[2]*force;
      }
			
      if (j<=i) continue;

      // Ca(i) - Ca(j)
      dx[0] = xca[i][0] - xca[j][0];
      dx[1] = xca[i][1] - xca[j][1];
      dx[2] = xca[i][2] - xca[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rC_ex0) {
	dr = r - rC_ex0;
	force = factorC*p*epsilon*k_excluded_C*pow(dr, p-1)/r;
			
	energy[ET_VEXCLUDED] += factorC*epsilon*k_excluded_C*pow(dr, p);
				
	f[alpha_carbons[i]][0] -= dx[0]*force;
	f[alpha_carbons[i]][1] -= dx[1]*force;
	f[alpha_carbons[i]][2] -= dx[2]*force;
				
	f[alpha_carbons[j]][0] -= -dx[0]*force;
	f[alpha_carbons[j]][1] -= -dx[1]*force;
	f[alpha_carbons[j]][2] -= -dx[2]*force;
      }

      // Cb(i) - Cb(j)
      dx[0] = xcb[i][0] - xcb[j][0];
      dx[1] = xcb[i][1] - xcb[j][1];
      dx[2] = xcb[i][2] - xcb[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rC_ex0) {
	dr = r - rC_ex0;
	force = factorC*p*epsilon*k_excluded_C*pow(dr, p-1)/r;
			
	energy[ET_VEXCLUDED] += factorC*epsilon*k_excluded_C*pow(dr, p);
				
	f[beta_atoms[i]][0] -= dx[0]*force;
	f[beta_atoms[i]][1] -= dx[1]*force;
	f[beta_atoms[i]][2] -= dx[2]*force;
				
	f[beta_atoms[j]][0] -= -dx[0]*force;
	f[beta_atoms[j]][1] -= -dx[1]*force;
	f[beta_atoms[j]][2] -= -dx[2]*force;
      }
			
      // O(i) - O(j)
      dx[0] = xo[i][0] - xo[j][0];
      dx[1] = xo[i][1] - xo[j][1];
      dx[2] = xo[i][2] - xo[j][2];
      r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      if (r<rO_ex0) {
	dr = r - rO_ex0;
	force = factorO*p*epsilon*k_excluded_O*pow(dr, p-1)/r;
			
	energy[ET_VEXCLUDED] += factorO*epsilon*k_excluded_O*pow(dr, p);
				
	f[oxygens[i]][0] -= dx[0]*force;
	f[oxygens[i]][1] -= dx[1]*force;
	f[oxygens[i]][2] -= dx[2]*force;
				
	f[oxygens[j]][0] -= -dx[0]*force;
	f[oxygens[j]][1] -= -dx[1]*force;
	f[oxygens[j]][2] -= -dx[2]*force;
      }
    }
  }
}

void FixBackbone::compute_r6_excluded_volume()
{
  double dx[3], r, rsq, force;
	
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) {
      // Ca(i) - Cb(j)
      dx[0] = xca[i][0] - xcb[j][0];
      dx[1] = xca[i][1] - xcb[j][1];
      dx[2] = xca[i][2] - xcb[j][2];
      rsq = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
      r = sqrt(rsq);
      if (r<rC_ex0) {
	force = -6*epsilon*k_excluded_C/pow(rsq, 4);
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_C/pow(rsq, 3);
				
	f[alpha_carbons[i]][0] -= dx[0]*force;
	f[alpha_carbons[i]][1] -= dx[1]*force;
	f[alpha_carbons[i]][2] -= dx[2]*force;
				
	f[beta_atoms[j]][0] -= -dx[0]*force;
	f[beta_atoms[j]][1] -= -dx[1]*force;
	f[beta_atoms[j]][2] -= -dx[2]*force;
      }
			
      if (j<=i) continue;

      // Ca(i) - Ca(j)
      dx[0] = xca[i][0] - xca[j][0];
      dx[1] = xca[i][1] - xca[j][1];
      dx[2] = xca[i][2] - xca[j][2];
      rsq = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
      r = sqrt(rsq);
      if (r<rC_ex0) {
	force = -6*epsilon*k_excluded_C/pow(rsq, 4);
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_C/pow(rsq, 3);
				
	f[alpha_carbons[i]][0] -= dx[0]*force;
	f[alpha_carbons[i]][1] -= dx[1]*force;
	f[alpha_carbons[i]][2] -= dx[2]*force;
				
	f[alpha_carbons[j]][0] -= -dx[0]*force;
	f[alpha_carbons[j]][1] -= -dx[1]*force;
	f[alpha_carbons[j]][2] -= -dx[2]*force;
      }

      // Cb(i) - Cb(j)
      dx[0] = xcb[i][0] - xcb[j][0];
      dx[1] = xcb[i][1] - xcb[j][1];
      dx[2] = xcb[i][2] - xcb[j][2];
      rsq = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
      r = sqrt(rsq);
      if (r<rC_ex0) {
	force = -6*epsilon*k_excluded_C/pow(rsq, 4);
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_C/pow(rsq, 3);
				
	f[beta_atoms[i]][0] -= dx[0]*force;
	f[beta_atoms[i]][1] -= dx[1]*force;
	f[beta_atoms[i]][2] -= dx[2]*force;
				
	f[beta_atoms[j]][0] -= -dx[0]*force;
	f[beta_atoms[j]][1] -= -dx[1]*force;
	f[beta_atoms[j]][2] -= -dx[2]*force;
      }
			
      // O(i) - O(j)
      dx[0] = xo[i][0] - xo[j][0];
      dx[1] = xo[i][1] - xo[j][1];
      dx[2] = xo[i][2] - xo[j][2];
      rsq = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
      r = sqrt(rsq);
      if (r<rO_ex0) {
	force = -6*epsilon*k_excluded_O/pow(rsq, 4);
			
	energy[ET_VEXCLUDED] += epsilon*k_excluded_O/pow(rsq, 3);
				
	f[oxygens[i]][0] -= dx[0]*force;
	f[oxygens[i]][1] -= dx[1]*force;
	f[oxygens[i]][2] -= dx[2]*force;
				
	f[oxygens[j]][0] -= -dx[0]*force;
	f[oxygens[j]][1] -= -dx[1]*force;
	f[oxygens[j]][2] -= -dx[2]*force;
      }
    }
  }
}

/*
  hbscl(1) short range additive
  hbscl(2) short range repulsive
  hbscl(3) med range additive
  hbscl(4) med range repulsive
  hbscl(5) med range anti-parallel non-add
  hbscl(6) med range parallel non-add
  hbscl(7) long range additive
  hbscl(8) long range repulsive
  hbscl(9) long range anti-parallel non-add
  hbscl(10) long range parallel non-add
  seq-dep terms (that depend non-additively on identities of residues in H-bond pair)
  hbscl(11) med range seq-dep anti-P for H bonded pair
  hbscl(12) med range seq-dep anti-P for non-H bonded pair
  hbscl(13) long range seq-dep anti-P for H bonded pair
  hbscl(14) long range seq-dep anti-P for non-H bonded pair
  hbscl(15) long range seq-dep parallel (all cross strand pairs same for parallel)
  seq-dep terms (that depend *additively* on identities of residues in H-bond pair)
  hbscl(16) seq-dep anti-P weight  
  hbscl(17) seq-dep parallel weight
*/

inline double FixBackbone::anti_HB(int res1, int res2, int k)
{
  return m_anti_HB[se_map[res1-'A']][se_map[res2-'A']][k];
}

inline double FixBackbone::anti_NHB(int res1, int res2, int k)
{
  return m_anti_NHB[se_map[res1-'A']][se_map[res2-'A']][k];
}

inline double FixBackbone::para_HB(int res1, int res2, int k)
{
  return m_para_HB[se_map[res1-'A']][se_map[res2-'A']][k];
}

inline double FixBackbone::para_one(int res)
{
  return m_para_one[se_map[res-'A']];
}

inline double FixBackbone::anti_one(int res)
{
  return m_anti_one[se_map[res-'A']];
}

inline double FixBackbone::get_water_gamma(int i_resno, int j_resno, int i_well, int ires_type, int jres_type, int water_prot_flag)
{
   
  if (!phosph_flag) {
    return water_gamma[i_well][ires_type][jres_type][water_prot_flag];
  } 
  else {
    if (phosph_map[i_resno] || phosph_map[j_resno]) {
      return phosph_water_gamma[i_well][ires_type][jres_type][water_prot_flag];
    }
    else {
      return water_gamma[i_well][ires_type][jres_type][water_prot_flag];
    }
  }
}

inline double FixBackbone::get_burial_gamma(int i_resno, int ires_type, int local_dens)
{
  if (!phosph_flag) {
    return burial_gamma[ires_type][local_dens];
  }
  else {
    if (phosph_map[i_resno]) {
      return burial_gamma[6][local_dens];  // returns burial gamma for glutamine instead of serine
    }
    else {
      return burial_gamma[ires_type][local_dens];
    }
  }
}

void FixBackbone::compute_dssp_hdrgn(int i, int j)
{
  if (R->rNO(i,j)>dssp_hdrgn_cut) return;

  bool i_repulsive = true, i_AP = true, i_P = true, i_theta[4] = {true, true, true, true};
  double lambda[4], R_NO[4], R_HO[4], theta[4], nu[2], th;
  double r_nu[2], prdnu[2], prd_theta[4][2], V[4], VTotal;
  double dxnu[2][3], xNO[4][3], xHO[4][3];
  double force, force1, force2, theta_sum;
  double theta_seq_anti_HB[2], theta_seq_anti_NHB[2], theta_seq_para_HB[2];
  int hb_class;
  int k;

  int i_resno = res_no[i]-1;
  int j_resno = res_no[j]-1;
  int i_chno = chain_no[i]-1;
  int j_chno = chain_no[j]-1;

  int im1 = res_no_l[i_resno-1];
  int im2 = res_no_l[i_resno-2];
  int ip1 = res_no_l[i_resno+1];
  int ip2 = res_no_l[i_resno+2];

  int jm1 = res_no_l[j_resno-1];
  int jm2 = res_no_l[j_resno-2];
  int jp1 = res_no_l[j_resno+1];
  int jp2 = res_no_l[j_resno+2];

  if (im1==-1 || im2==-1 || ip1==-1 || ip2==-1 || jm1==-1 || jm2==-1 || jp1==-1 || jp2==-1)return;
  
  if (res_info[im1]!=LOCAL && res_info[im1]!=GHOST) return;
  if (res_info[im2]!=LOCAL && res_info[im2]!=GHOST) return;
  if (res_info[ip1]!=LOCAL && res_info[ip1]!=GHOST) return;
  if (res_info[ip2]!=LOCAL && res_info[ip2]!=GHOST) return;
  if (res_info[jm1]!=LOCAL && res_info[jm1]!=GHOST) return;
  if (res_info[jm2]!=LOCAL && res_info[jm2]!=GHOST) return;
  if (res_info[jp1]!=LOCAL && res_info[jp1]!=GHOST) return;
  if (res_info[jp2]!=LOCAL && res_info[jp2]!=GHOST) return;
  
  if (chain_no[im2] != chain_no[ip2]) return;
  if (chain_no[jm2] != chain_no[jp2]) return;
  
  if ( isLast(j) || se[j_resno+1]=='P' ) i_repulsive = false;
  if ( isFirst(i) || isLast(j) || se[i_resno]=='P' ) i_AP = false;
  if ( i_resno>=(ch_pos[i_chno]+ch_len[i_chno]-1)-2 || isLast(j) || se[i_resno+2]=='P' ) i_P = false;
	
  /*if ( j_resno==n-1 || se[j_resno+1]=='P' ) i_repulsive = false;
	if ( i_resno==0 || j_resno==n-1 || se[i_resno]=='P' ) i_AP = false;
	if ( i_resno==n-2 || j_resno==n-1 || se[i_resno+2]=='P' ) i_P = false;*/

  for (k=0;k<2;++k) {
	  if (i_AP) {
		  theta_seq_anti_HB[k]=0.5*anti_HB(se[i_resno], se[j_resno], k);
		  theta_seq_anti_NHB[k]=0.25*( anti_NHB( se[i_resno+1], se[j_resno-1], k) + anti_NHB(se[i_resno-1], se[j_resno+1], k) );
	  } else {
		  theta_seq_anti_HB[k] = 0.0;
		  theta_seq_anti_NHB[k] = 0.0;
	  }
	  if (i_P) {
		  theta_seq_para_HB[k]=para_HB(se[i_resno+1], se[j_resno], k);
	  } else {
		  theta_seq_para_HB[k] = 0.0;
	  }
  }

  if ( i_chno==j_chno && abs(j_resno - i_resno) < 45 ) {
	  if (abs(j_resno - i_resno) < 4) { // |j-i|<4, class 1
		  lambda[0] = -hbscl[0][0];
		  lambda[1] = -hbscl[0][1];
		  lambda[2] = 0.0;
		  lambda[3] = 0.0;
		  hb_class = 1;
	  } else if (abs(j_resno - i_resno) < 18) { // 4 =< |j-i| < 18, class 2
		  if (aps[n_rama_par-1][i_resno]==1.0 && aps[n_rama_par-1][j_resno]==1.0) {
			  lambda[0] = -hbscl[1][0];
			  lambda[1] = -hbscl[1][1];
			  lambda[2] = -hbscl[1][2]
			                        -hbscl[1][3]*theta_seq_anti_HB[0]
			                        -hbscl[1][4]*theta_seq_anti_NHB[0]
			                        -hbscl[1][5]*(anti_one(se[i_resno])+anti_one(se[j_resno]));
			  lambda[3] = -hbscl[1][6];
		  } else { // 4 =< |j-i| < 18 and ssweight predicts not beta, class 2
			  lambda[0] = 0.0;
			  lambda[1] = -hbscl[1][1];
			  lambda[2] = 0.0;
			  lambda[3] = 0.0;
		  }
		  hb_class = 2;
	  } else if (abs(j_resno - i_resno) < 45) { // 18 =< |j-i| < 45, class 3
		  lambda[0] = -hbscl[2][0];
		  lambda[1] = -hbscl[2][1];
		  lambda[2] = -hbscl[2][2]
		                        -hbscl[2][3]*theta_seq_anti_HB[1]
		                        -hbscl[2][4]*theta_seq_anti_NHB[1]
		                        -hbscl[2][5]*(anti_one(se[i_resno])+anti_one(se[j_resno]));
		  lambda[3] = -hbscl[2][6]
		                        -hbscl[2][7]*theta_seq_para_HB[1]
		                        -hbscl[2][8]*(para_one(se[i_resno+1])+para_one(se[j_resno]));
		  hb_class = 3;
	  }
  } else { // (|j-i| >= 45) OR (j,i are on different chain), class 4
	  lambda[0] = -hbscl[3][0];
	  lambda[1] = -hbscl[3][1];
	  lambda[2] = -hbscl[3][2]
	                        -hbscl[3][3]*theta_seq_anti_HB[1]
	                        -hbscl[3][4]*theta_seq_anti_NHB[1]
	                        -hbscl[3][5]*(anti_one(se[i_resno])+anti_one(se[j_resno]));
	  lambda[3] = -hbscl[3][6]
	                        -hbscl[3][7]*theta_seq_para_HB[1]
	                        -hbscl[3][8]*(para_one(se[i_resno+1])+para_one(se[j_resno]));
	  hb_class = 4;
  }

  i_theta[1] = i_repulsive;
  i_theta[2] = i_AP;
  i_theta[3] = i_P;

  R_NO[0]=R->rNO(i, j);
  R_HO[0]=R->rHO(i, j);

  xNO[0][0] = xo[i][0] - xn[j][0];
  xNO[0][1] = xo[i][1] - xn[j][1];
  xNO[0][2] = xo[i][2] - xn[j][2];

  xHO[0][0] = xo[i][0] - xh[j][0];
  xHO[0][1] = xo[i][1] - xh[j][1];
  xHO[0][2] = xo[i][2] - xh[j][2];

  if (i_repulsive) {
	  R_NO[1]=R->rNO(i, jp1);
	  R_HO[1]=R->rHO(i, jp1);

	  xNO[1][0] = xo[i][0] - xn[jp1][0];
	  xNO[1][1] = xo[i][1] - xn[jp1][1];
	  xNO[1][2] = xo[i][2] - xn[jp1][2];

	  xHO[1][0] = xo[i][0] - xh[jp1][0];
	  xHO[1][1] = xo[i][1] - xh[jp1][1];
	  xHO[1][2] = xo[i][2] - xh[jp1][2];
  }

  if (i_AP) {
	  R_NO[2]=R->rNO(j, i);
	  R_HO[2]=R->rHO(j, i);

	  xNO[2][0] = xo[j][0] - xn[i][0];
	  xNO[2][1] = xo[j][1] - xn[i][1];
	  xNO[2][2] = xo[j][2] - xn[i][2];

	  xHO[2][0] = xo[j][0] - xh[i][0];
	  xHO[2][1] = xo[j][1] - xh[i][1];
	  xHO[2][2] = xo[j][2] - xh[i][2];
  }

  if (i_P) {
	  R_NO[3]=R->rNO(j, ip2);
	  R_HO[3]=R->rHO(j, ip2);

	  xNO[3][0] = xo[j][0] - xn[ip2][0];
	  xNO[3][1] = xo[j][1] - xn[ip2][1];
	  xNO[3][2] = xo[j][2] - xn[ip2][2];

	  xHO[3][0] = xo[j][0] - xh[ip2][0];
	  xHO[3][1] = xo[j][1] - xh[ip2][1];
	  xHO[3][2] = xo[j][2] - xh[ip2][2];
  }


  for (k=0;k<4;++k) {
	  if (i_theta[k]) {
		  theta[k] = exp( - pow(R_NO[k] - NO_zero, 2)/(2.0*pow(sigma_NO, 2)) - pow(R_HO[k] - HO_zero, 2)/(2.0*pow(sigma_HO, 2)) );

		  prd_theta[k][0] = - (R_NO[k] - NO_zero)/(pow(sigma_NO, 2)*R_NO[k]);
		  prd_theta[k][1] = - (R_HO[k] - HO_zero)/(pow(sigma_HO, 2)*R_HO[k]);
	  } else {
		  theta[k] = 0.0;
		  prd_theta[k][0] = prd_theta[k][1] = 0.0;
	  }
  }

  if (im2 > 0 && !isFirst(im1) && !isFirst(im2) && ip2 < nn && !isLast(ip1) && hb_class!=2) {
	  dxnu[0][0] = xca[ip2][0]-xca[im2][0];
	  dxnu[0][1] = xca[ip2][1]-xca[im2][1];
	  dxnu[0][2] = xca[ip2][2]-xca[im2][2];

	  r_nu[0] = sqrt (pow(dxnu[0][0], 2) + pow(dxnu[0][1], 2) + pow(dxnu[0][2], 2) );

	  th = tanh(pref[0]*(r_nu[0] - d_nu0));
	  nu[0] = 0.5*(1+th);

	  prdnu[0] = 0.5*pref[0]*(1-pow(th,2))/r_nu[0];
  } else nu[0] = 1.0;

  if (jm2 > 0 && !isFirst(jm1) && !isFirst(jm2) && jp2 < nn && !isLast(jp1) && hb_class!=2) {
	  dxnu[1][0] = xca[jp2][0]-xca[jm2][0];
	  dxnu[1][1] = xca[jp2][1]-xca[jm2][1];
	  dxnu[1][2] = xca[jp2][2]-xca[jm2][2];

	  r_nu[1] = sqrt (pow(dxnu[1][0], 2) + pow(dxnu[1][1], 2) + pow(dxnu[1][2], 2) );

	  th = tanh(pref[1]*(r_nu[1] - d_nu0));
	  nu[1] = 0.5*(1+th);

	  prdnu[1] = 0.5*pref[1]*(1-pow(th,2))/r_nu[1];
  } else nu[1] = 1.0;


  theta_sum = lambda[0]*theta[0] 
                              + lambda[1]*theta[0]*theta[1] 
                              + lambda[2]*theta[0]*theta[2] 
                              + lambda[3]*theta[0]*theta[3];

  V[0] = k_dssp*epsilon*lambda[0]*theta[0]*nu[0]*nu[1];
  V[1] = k_dssp*epsilon*lambda[1]*theta[0]*theta[1]*nu[0]*nu[1];
  V[2] = k_dssp*epsilon*lambda[2]*theta[0]*theta[2]*nu[0]*nu[1];
  V[3] = k_dssp*epsilon*lambda[3]*theta[0]*theta[3]*nu[0]*nu[1];

  VTotal = V[0] + V[1] + V[2] + V[3];

  energy[ET_DSSP] +=  epsilon*VTotal;

  if (im2 > 0 && !isFirst(im1) && !isFirst(im2) && ip2 < nn && !isLast(ip1) && hb_class!=2) {
	  force = k_dssp*epsilon*theta_sum*prdnu[0]*nu[1];
	  f[alpha_carbons[im2]][0] -= -force*dxnu[0][0];
	  f[alpha_carbons[im2]][1] -= -force*dxnu[0][1];
	  f[alpha_carbons[im2]][2] -= -force*dxnu[0][2];

	  f[alpha_carbons[ip2]][0] -= force*dxnu[0][0];
	  f[alpha_carbons[ip2]][1] -= force*dxnu[0][1];
	  f[alpha_carbons[ip2]][2] -= force*dxnu[0][2];
  }

  if (jm2 > 0 && !isFirst(jm1) && !isFirst(jm2) && jp2 < nn && !isLast(jp1) && hb_class!=2) {
	  force = k_dssp*epsilon*theta_sum*nu[0]*prdnu[1];
	  f[alpha_carbons[jm2]][0] -= -force*dxnu[1][0];
	  f[alpha_carbons[jm2]][1] -= -force*dxnu[1][1];
	  f[alpha_carbons[jm2]][2] -= -force*dxnu[1][2];

	  f[alpha_carbons[jp2]][0] -= force*dxnu[1][0];
	  f[alpha_carbons[jp2]][1] -= force*dxnu[1][1];
	  f[alpha_carbons[jp2]][2] -= force*dxnu[1][2];
  }

  f[alpha_carbons[jm1]][0] -= -VTotal*(an*prd_theta[0][0]*xNO[0][0] + ah*prd_theta[0][1]*xHO[0][0]);
  f[alpha_carbons[jm1]][1] -= -VTotal*(an*prd_theta[0][0]*xNO[0][1] + ah*prd_theta[0][1]*xHO[0][1]);
  f[alpha_carbons[jm1]][2] -= -VTotal*(an*prd_theta[0][0]*xNO[0][2] + ah*prd_theta[0][1]*xHO[0][2]);

  f[alpha_carbons[j]][0] -= -VTotal*(bn*prd_theta[0][0]*xNO[0][0] + bh*prd_theta[0][1]*xHO[0][0]);
  f[alpha_carbons[j]][1] -= -VTotal*(bn*prd_theta[0][0]*xNO[0][1] + bh*prd_theta[0][1]*xHO[0][1]);
  f[alpha_carbons[j]][2] -= -VTotal*(bn*prd_theta[0][0]*xNO[0][2] + bh*prd_theta[0][1]*xHO[0][2]);

  f[oxygens[jm1]][0] -= -VTotal*(cn*prd_theta[0][0]*xNO[0][0] + ch*prd_theta[0][1]*xHO[0][0]);
  f[oxygens[jm1]][1] -= -VTotal*(cn*prd_theta[0][0]*xNO[0][1] + ch*prd_theta[0][1]*xHO[0][1]);
  f[oxygens[jm1]][2] -= -VTotal*(cn*prd_theta[0][0]*xNO[0][2] + ch*prd_theta[0][1]*xHO[0][2]);

  f[oxygens[i]][0] -= VTotal*(prd_theta[0][0]*xNO[0][0] + prd_theta[0][1]*xHO[0][0]);
  f[oxygens[i]][1] -= VTotal*(prd_theta[0][0]*xNO[0][1] + prd_theta[0][1]*xHO[0][1]);
  f[oxygens[i]][2] -= VTotal*(prd_theta[0][0]*xNO[0][2] + prd_theta[0][1]*xHO[0][2]);

  if (i_repulsive) {
	  f[alpha_carbons[j]][0] -= -V[1]*(an*prd_theta[1][0]*xNO[1][0] + ah*prd_theta[1][1]*xHO[1][0]);
	  f[alpha_carbons[j]][1] -= -V[1]*(an*prd_theta[1][0]*xNO[1][1] + ah*prd_theta[1][1]*xHO[1][1]);
	  f[alpha_carbons[j]][2] -= -V[1]*(an*prd_theta[1][0]*xNO[1][2] + ah*prd_theta[1][1]*xHO[1][2]);

	  f[alpha_carbons[jp1]][0] -= -V[1]*(bn*prd_theta[1][0]*xNO[1][0] + bh*prd_theta[1][1]*xHO[1][0]);
	  f[alpha_carbons[jp1]][1] -= -V[1]*(bn*prd_theta[1][0]*xNO[1][1] + bh*prd_theta[1][1]*xHO[1][1]);
	  f[alpha_carbons[jp1]][2] -= -V[1]*(bn*prd_theta[1][0]*xNO[1][2] + bh*prd_theta[1][1]*xHO[1][2]);

	  f[oxygens[j]][0] -= -V[1]*(cn*prd_theta[1][0]*xNO[1][0] + ch*prd_theta[1][1]*xHO[1][0]);
	  f[oxygens[j]][1] -= -V[1]*(cn*prd_theta[1][0]*xNO[1][1] + ch*prd_theta[1][1]*xHO[1][1]);
	  f[oxygens[j]][2] -= -V[1]*(cn*prd_theta[1][0]*xNO[1][2] + ch*prd_theta[1][1]*xHO[1][2]);

	  f[oxygens[i]][0] -= V[1]*(prd_theta[1][0]*xNO[1][0] + prd_theta[1][1]*xHO[1][0]);
	  f[oxygens[i]][1] -= V[1]*(prd_theta[1][0]*xNO[1][1] + prd_theta[1][1]*xHO[1][1]);
	  f[oxygens[i]][2] -= V[1]*(prd_theta[1][0]*xNO[1][2] + prd_theta[1][1]*xHO[1][2]);
  }


  if (i_AP) {
	  f[alpha_carbons[im1]][0] -= -V[2]*(an*prd_theta[2][0]*xNO[2][0] + ah*prd_theta[2][1]*xHO[2][0]);
	  f[alpha_carbons[im1]][1] -= -V[2]*(an*prd_theta[2][0]*xNO[2][1] + ah*prd_theta[2][1]*xHO[2][1]);
	  f[alpha_carbons[im1]][2] -= -V[2]*(an*prd_theta[2][0]*xNO[2][2] + ah*prd_theta[2][1]*xHO[2][2]);

	  f[alpha_carbons[i]][0] -= -V[2]*(bn*prd_theta[2][0]*xNO[2][0] + bh*prd_theta[2][1]*xHO[2][0]);
	  f[alpha_carbons[i]][1] -= -V[2]*(bn*prd_theta[2][0]*xNO[2][1] + bh*prd_theta[2][1]*xHO[2][1]);
	  f[alpha_carbons[i]][2] -= -V[2]*(bn*prd_theta[2][0]*xNO[2][2] + bh*prd_theta[2][1]*xHO[2][2]);

	  f[oxygens[im1]][0] -= -V[2]*(cn*prd_theta[2][0]*xNO[2][0] + ch*prd_theta[2][1]*xHO[2][0]);
	  f[oxygens[im1]][1] -= -V[2]*(cn*prd_theta[2][0]*xNO[2][1] + ch*prd_theta[2][1]*xHO[2][1]);
	  f[oxygens[im1]][2] -= -V[2]*(cn*prd_theta[2][0]*xNO[2][2] + ch*prd_theta[2][1]*xHO[2][2]);

	  f[oxygens[j]][0] -= V[2]*(prd_theta[2][0]*xNO[2][0] + prd_theta[2][1]*xHO[2][0]);
	  f[oxygens[j]][1] -= V[2]*(prd_theta[2][0]*xNO[2][1] + prd_theta[2][1]*xHO[2][1]);
	  f[oxygens[j]][2] -= V[2]*(prd_theta[2][0]*xNO[2][2] + prd_theta[2][1]*xHO[2][2]);
  }


  if (i_P) {
	  f[alpha_carbons[ip1]][0] -= -V[3]*(an*prd_theta[3][0]*xNO[3][0] + ah*prd_theta[3][1]*xHO[3][0]);
	  f[alpha_carbons[ip1]][1] -= -V[3]*(an*prd_theta[3][0]*xNO[3][1] + ah*prd_theta[3][1]*xHO[3][1]);
	  f[alpha_carbons[ip1]][2] -= -V[3]*(an*prd_theta[3][0]*xNO[3][2] + ah*prd_theta[3][1]*xHO[3][2]);

	  f[alpha_carbons[ip2]][0] -= -V[3]*(bn*prd_theta[3][0]*xNO[3][0] + bh*prd_theta[3][1]*xHO[3][0]);
	  f[alpha_carbons[ip2]][1] -= -V[3]*(bn*prd_theta[3][0]*xNO[3][1] + bh*prd_theta[3][1]*xHO[3][1]);
	  f[alpha_carbons[ip2]][2] -= -V[3]*(bn*prd_theta[3][0]*xNO[3][2] + bh*prd_theta[3][1]*xHO[3][2]);

	  f[oxygens[ip1]][0] -= -V[3]*(cn*prd_theta[3][0]*xNO[3][0] + ch*prd_theta[3][1]*xHO[3][0]);
	  f[oxygens[ip1]][1] -= -V[3]*(cn*prd_theta[3][0]*xNO[3][1] + ch*prd_theta[3][1]*xHO[3][1]);
	  f[oxygens[ip1]][2] -= -V[3]*(cn*prd_theta[3][0]*xNO[3][2] + ch*prd_theta[3][1]*xHO[3][2]);

	  f[oxygens[j]][0] -= V[3]*(prd_theta[3][0]*xNO[3][0] + prd_theta[3][1]*xHO[3][0]);
	  f[oxygens[j]][1] -= V[3]*(prd_theta[3][0]*xNO[3][1] + prd_theta[3][1]*xHO[3][1]);
	  f[oxygens[j]][2] -= V[3]*(prd_theta[3][0]*xNO[3][2] + prd_theta[3][1]*xHO[3][2]);
  }
}

void FixBackbone::compute_P_AP_potential(int i, int j)
{
  double K, force[2], dx[2][3];
  bool i_AP_med, i_AP_long, i_P;

  int i_resno = res_no[i]-1;
  int j_resno = res_no[j]-1;
  int i_chno = chain_no[i]-1;
  int j_chno = chain_no[j]-1;

  // anti-parallel hairpin
  i_AP_med = i_chno==j_chno && i_resno<n-(i_med_min+2*i_diff_P_AP) && j_resno>=i_resno+(i_med_min+2*i_diff_P_AP) && j_resno<=MIN(i_resno+i_med_max+2*i_diff_P_AP,n-1);

  // anti-parallel beta strands
  int ip4 = res_no_l[i_resno+i_diff_P_AP];
  int jm4 = res_no_l[j_resno-i_diff_P_AP];
  int ip4_chno, jm4_chno;
  if (ip4==-1 || jm4==-1){ //can't access chain_no[] properly
	i_AP_long = false;
  }
  else {
	  ip4_chno = chain_no[ip4]-1;
	  jm4_chno = chain_no[jm4]-1;
	  i_AP_long = (i_chno==j_chno && i_resno<n-(i_med_max+2*i_diff_P_AP+1) && j_resno>=i_resno+(i_med_max+2*i_diff_P_AP+1) && j_resno<n) || (i_chno!=j_chno && ip4_chno==i_chno && jm4_chno==j_chno);
  }

  //parallel beta strands
  int jp4, jp4_chno;
  if (j_resno+i_diff_P_AP >= n){
	  i_P = false ;
  }
  else{
	jp4 = res_no_l[j_resno+i_diff_P_AP];
	if (jp4==-1){
		i_P = false;
	}
	else {
		jp4_chno = chain_no[jp4]-1;
		i_P = (i_chno==jp4_chno && i_resno<n-(i_med_max+1+i_diff_P_AP) && j_resno>=i_resno+(i_med_max+1) && j_resno<n-i_diff_P_AP) || (i_chno!=j_chno && ip4_chno==i_chno && jp4_chno==j_chno);
	}
  }

  if (i_AP_med || i_AP_long) {
    if ((res_info[ip4]!=LOCAL && res_info[ip4]!=GHOST) || (res_info[jm4]!=LOCAL && res_info[jm4]!=GHOST)) return;
    if (aps[n_rama_par-1][i_resno]==1.0 && aps[n_rama_par-1][j_resno]==1.0) {
      K = (i_AP_med ? k_P_AP[0] : 0.0) + (i_AP_long ? k_P_AP[1]*k_betapred_P_AP : 0.0);
    } else {
      K = (i_AP_med ? k_P_AP[0] : 0.0) + (i_AP_long ? k_P_AP[1] : 0.0);
    }

    energy[ET_PAP] += -k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->nu(ip4, jm4);

    dx[0][0] = xca[i][0] - xca[j][0];
    dx[0][1] = xca[i][1] - xca[j][1];
    dx[0][2] = xca[i][2] - xca[j][2];

    dx[1][0] = xca[ip4][0] - xca[jm4][0];
    dx[1][1] = xca[ip4][1] - xca[jm4][1];
    dx[1][2] = xca[ip4][2] - xca[jm4][2];

    force[0] = k_global_P_AP*epsilon*K*p_ap->prd_nu(i, j)*p_ap->nu(ip4, jm4);
    force[1] = k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->prd_nu(ip4, jm4);
	
    f[alpha_carbons[i]][0] -= force[0]*dx[0][0];
    f[alpha_carbons[i]][1] -= force[0]*dx[0][1];
    f[alpha_carbons[i]][2] -= force[0]*dx[0][2];
	
    f[alpha_carbons[j]][0] -= -force[0]*dx[0][0];
    f[alpha_carbons[j]][1] -= -force[0]*dx[0][1];
    f[alpha_carbons[j]][2] -= -force[0]*dx[0][2];

    f[alpha_carbons[ip4]][0] -= force[1]*dx[1][0];
    f[alpha_carbons[ip4]][1] -= force[1]*dx[1][1];
    f[alpha_carbons[ip4]][2] -= force[1]*dx[1][2];

    f[alpha_carbons[jm4]][0] -= -force[1]*dx[1][0];
    f[alpha_carbons[jm4]][1] -= -force[1]*dx[1][1];
    f[alpha_carbons[jm4]][2] -= -force[1]*dx[1][2];
  }
    
  if (i_P) {
    if ((res_info[ip4]!=LOCAL && res_info[ip4]!=GHOST) || (res_info[jp4]!=LOCAL && res_info[jp4]!=GHOST)) return;
    if (aps[n_rama_par-1][i_resno]==1.0 && aps[n_rama_par-1][j_resno]==1.0) {
      K = k_P_AP[2]*k_betapred_P_AP;
    } else {
      K = k_P_AP[2];
    }

    energy[ET_PAP] += -k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->nu(ip4, jp4);

    dx[0][0] = xca[i][0] - xca[j][0];
    dx[0][1] = xca[i][1] - xca[j][1];
    dx[0][2] = xca[i][2] - xca[j][2];

    dx[1][0] = xca[ip4][0] - xca[jp4][0];
    dx[1][1] = xca[ip4][1] - xca[jp4][1];
    dx[1][2] = xca[ip4][2] - xca[jp4][2];

    force[0] = k_global_P_AP*epsilon*K*p_ap->prd_nu(i, j)*p_ap->nu(ip4, jp4);
    force[1] = k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->prd_nu(ip4, jp4);

    f[alpha_carbons[i]][0] -= force[0]*dx[0][0];
    f[alpha_carbons[i]][1] -= force[0]*dx[0][1];
    f[alpha_carbons[i]][2] -= force[0]*dx[0][2];

    f[alpha_carbons[j]][0] -= -force[0]*dx[0][0];
    f[alpha_carbons[j]][1] -= -force[0]*dx[0][1];
    f[alpha_carbons[j]][2] -= -force[0]*dx[0][2];

    f[alpha_carbons[ip4]][0] -= force[1]*dx[1][0];
    f[alpha_carbons[ip4]][1] -= force[1]*dx[1][1];
    f[alpha_carbons[ip4]][2] -= force[1]*dx[1][2];

    f[alpha_carbons[jp4]][0] -= -force[1]*dx[1][0];
    f[alpha_carbons[jp4]][1] -= -force[1]*dx[1][1];
    f[alpha_carbons[jp4]][2] -= -force[1]*dx[1][2];
  }
}

void FixBackbone::compute_water_potential(int i, int j)
{	
  if (chain_no[i]==chain_no[j] && res_no[j]-res_no[i]<contact_cutoff) return;
  if (chain_no[i]!=chain_no[j] && res_no[j] < res_no[i]) return;

  double dx[3], sigma_gamma, theta_gamma, force;
  double *xi, *xj, *xk;
  double water_gamma_0, water_gamma_1;
  int iatom, jatom, katom, i_well, k;
  int k_resno, k_chno;
  bool direct_contact;

  int i_resno = res_no[i]-1;
  int j_resno = res_no[j]-1;
	
  int i_chno = chain_no[i]-1;
  int j_chno = chain_no[j]-1;
	
  int ires_type = se_map[se[i_resno]-'A'];
  int jres_type = se_map[se[j_resno]-'A'];
  
  if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
  else { xi = xcb[i]; iatom  = beta_atoms[i]; }
  if (se[j_resno]=='G') { xj = xca[j]; jatom = alpha_carbons[j]; }
  else { xj = xcb[j]; jatom  = beta_atoms[j]; }
	
  dx[0] = xi[0] - xj[0];
  dx[1] = xi[1] - xj[1];
  dx[2] = xi[2] - xj[2];

  for (i_well=0;i_well<n_wells;++i_well) {
    if (!well_flag[i_well]) continue;
    if (fabs(well->theta(i, j, i_well))<delta) continue;

    direct_contact = false;
		
    water_gamma_0 = get_water_gamma(i_resno, j_resno, i_well, ires_type, jres_type, 0);
    water_gamma_1 = get_water_gamma(i_resno, j_resno, i_well, ires_type, jres_type, 1);
		
    // Optimization for gamma[0]==gamma[1]
    if (fabs(water_gamma_0 - water_gamma_1)<delta) direct_contact = true;
		
    if (direct_contact) {
      sigma_gamma = (water_gamma_0 + water_gamma_1)/2;
      theta_gamma = 0;
    } else {	
      sigma_gamma = (1.0 - well->sigma(i, j))*water_gamma_0 + well->sigma(i, j)*water_gamma_1;
      theta_gamma = (water_gamma_1 - water_gamma_0)*well->theta(i, j, i_well);
    }		
	    
    energy[ET_WATER] += -epsilon*k_water*sigma_gamma*well->theta(i, j, i_well);
		  
    force = epsilon*k_water*sigma_gamma*well->prd_theta(i, j, i_well);
		
    f[iatom][0] += force*dx[0];
    f[iatom][1] += force*dx[1];
    f[iatom][2] += force*dx[2];
		
    f[jatom][0] += -force*dx[0];
    f[jatom][1] += -force*dx[1];
    f[jatom][2] += -force*dx[2];
		
    if (!direct_contact) {
      for (k=0;k<nn;++k) {
	if (res_info[k]==OFF) continue;

	if (se[res_no[k]-1]=='G') { xk = xca[k]; katom = alpha_carbons[k]; }
	else { xk = xcb[k]; katom  = beta_atoms[k]; }
				
	k_resno = res_no[k]-1;
	k_chno = chain_no[k]-1;
				
	if (abs(k_resno-i_resno)>1 || k_chno!=i_chno) {
	  dx[0] = xi[0] - xk[0];
	  dx[1] = xi[1] - xk[1];
	  dx[2] = xi[2] - xk[2];
					
	  force = epsilon*k_water*theta_gamma*well->prd_H(i)*well->H(j)*well->prd_theta(i, k, 0);
					
	  f[iatom][0] += force*dx[0];
	  f[iatom][1] += force*dx[1];
	  f[iatom][2] += force*dx[2];
					
	  f[katom][0] += -force*dx[0];
	  f[katom][1] += -force*dx[1];
	  f[katom][2] += -force*dx[2];
	}
	if (abs(k_resno-j_resno)>1 || k_chno!=j_chno) {
	  dx[0] = xj[0] - xk[0];
	  dx[1] = xj[1] - xk[1];
	  dx[2] = xj[2] - xk[2];
	
	  force = epsilon*k_water*theta_gamma*well->H(i)*well->prd_H(j)*well->prd_theta(j, k, 0);
				
	  f[jatom][0] += force*dx[0];
	  f[jatom][1] += force*dx[1];
	  f[jatom][2] += force*dx[2];
					
	  f[katom][0] += -force*dx[0];
	  f[katom][1] += -force*dx[1];
	  f[katom][2] += -force*dx[2];
	}
      }
    }
  }
}

void FixBackbone::compute_burial_potential(int i)
{
  double t[3][2], dx[3], force[3], force2, *xi, *xk;
  double burial_gamma_0, burial_gamma_1, burial_gamma_2;
  int iatom, katom, k, k_resno, k_chno;
  int i_resno = res_no[i]-1;
  int i_chno = chain_no[i]-1;
  
  int ires_type = se_map[se[i_resno]-'A'];
  
  if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
  else { xi = xcb[i]; iatom  = beta_atoms[i]; }
  
  t[0][0] = tanh( burial_kappa*(well->ro(i) - burial_ro_min[0]) );
  t[0][1] = tanh( burial_kappa*(burial_ro_max[0] - well->ro(i)) );
  t[1][0] = tanh( burial_kappa*(well->ro(i) - burial_ro_min[1]) );
  t[1][1] = tanh( burial_kappa*(burial_ro_max[1] - well->ro(i)) );
  t[2][0] = tanh( burial_kappa*(well->ro(i) - burial_ro_min[2]) );
  t[2][1] = tanh( burial_kappa*(burial_ro_max[2] - well->ro(i)) );
  
  burial_gamma_0 = get_burial_gamma(i_resno, ires_type, 0);
  burial_gamma_1 = get_burial_gamma(i_resno, ires_type, 1);
  burial_gamma_2 = get_burial_gamma(i_resno, ires_type, 2);

  energy[ET_BURIAL] += -0.5*epsilon*k_burial*burial_gamma_0*(t[0][0] + t[0][1]);
  energy[ET_BURIAL] += -0.5*epsilon*k_burial*burial_gamma_1*(t[1][0] + t[1][1]);
  energy[ET_BURIAL] += -0.5*epsilon*k_burial*burial_gamma_2*(t[2][0] + t[2][1]);
  
  force[0] = 0.5*epsilon*k_burial*burial_gamma_0*burial_kappa*( t[0][1]*t[0][1] - t[0][0]*t[0][0] );
  force[1] = 0.5*epsilon*k_burial*burial_gamma_1*burial_kappa*( t[1][1]*t[1][1] - t[1][0]*t[1][0] );
  force[2] = 0.5*epsilon*k_burial*burial_gamma_2*burial_kappa*( t[2][1]*t[2][1] - t[2][0]*t[2][0] );
  
  for (k=0;k<nn;++k) {
    if (res_info[k]==OFF) continue;
  
    k_resno = res_no[k]-1;
    k_chno = chain_no[k]-1;
    
    if (abs(k_resno-i_resno)>1 || i_chno!=k_chno) {
      if (se[res_no[k]-1]=='G') { xk = xca[k]; katom = alpha_carbons[k]; }
      else { xk = xcb[k]; katom  = beta_atoms[k]; }

      dx[0] = xi[0] - xk[0];
      dx[1] = xi[1] - xk[1];
      dx[2] = xi[2] - xk[2];
      
      force2 = (force[0] + force[1] + force[2])*well->prd_theta(i,k,0);
      
      f[iatom][0] += force2*dx[0];
      f[iatom][1] += force2*dx[1];
      f[iatom][2] += force2*dx[2];
      
      f[katom][0] += -force2*dx[0];
      f[katom][1] += -force2*dx[1];
      f[katom][2] += -force2*dx[2];
    }
  }
}

void FixBackbone::compute_helix_potential(int i)
{
  int i_resno = res_no[i]-1;
  int j = res_no_l[i_resno+helix_i_diff];
  int j_resno = res_no[j]-1;
  int jm1 = res_no_l[j_resno-1];
  int i_chno = chain_no[i]-1;
  int j_chno = chain_no[j]-1;

  if (j==-1   || (res_info[j]  !=LOCAL && res_info[j]  !=GHOST)) return;
  if (R->rNO(i,j)>helix_cutoff) return;
  if (jm1==-1 || (res_info[jm1]!=LOCAL && res_info[jm1]!=GHOST)) return;
  if (i_chno!=j_chno) return;

  double R_NO, R_HO, xNO[3], xHO[3], dx[3];
  double pair_theta, prd_pair_theta[2], prob_sum;
  double pair_theta_gamma, sigmma_gamma, V;
  double force;
  double *xi, *xj, *xk;
  double hp1, hp2;
  int iatom, jatom, katom, k;
  int k_resno, k_chno;

  int ires_type = se_map[se[i_resno]-'A'];
  int jres_type = se_map[se[j_resno]-'A'];

  R_NO=R->rNO(i, j);
  R_HO=R->rHO(i, j);

  xNO[0] = xo[i][0] - xn[j][0];
  xNO[1] = xo[i][1] - xn[j][1];
  xNO[2] = xo[i][2] - xn[j][2];

  xHO[0] = xo[i][0] - xh[j][0];
  xHO[1] = xo[i][1] - xh[j][1];
  xHO[2] = xo[i][2] - xh[j][2];
	
  double h4probi = h4prob[ires_type];
  if (se[i_resno]=='P' && pro_accepter_flag) h4probi = h4prob_pro_accepter;
	
  prob_sum = h4probi + h4prob[jres_type]; // sequence-identity weight
	
  pair_theta = prob_sum*exp( - pow(R_NO - helix_NO_zero, 2)/(2.0*pow(helix_sigma_NO, 2)) - pow(R_HO - helix_HO_zero, 2)/(2.0*pow(helix_sigma_HO, 2)) );

  prd_pair_theta[0] = - (R_NO - helix_NO_zero)/(pow(helix_sigma_NO, 2)*R_NO);
  prd_pair_theta[1] = - (R_HO - helix_HO_zero)/(pow(helix_sigma_HO, 2)*R_HO);
	
  if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
  else { xi = xcb[i]; iatom  = beta_atoms[i]; if (iatom==-1) return;}
  if (se[j_resno]=='G') { xj = xca[j]; jatom = alpha_carbons[j]; }
  else { xj = xcb[j]; jatom  = beta_atoms[j]; if (jatom==-1) return;}
	
  dx[0] = xi[0] - xj[0];
  dx[1] = xi[1] - xj[1];
  dx[2] = xi[2] - xj[2];
	
  sigmma_gamma = helix_gamma_p*(1.0-helix_well->sigma(i, j)) + helix_gamma_w*helix_well->sigma(i, j);

  pair_theta_gamma = -epsilon*k_helix*(helix_gamma_w - helix_gamma_p)*pair_theta;
	
  V = -epsilon*k_helix*sigmma_gamma*pair_theta;

  energy[ET_HELIX] += V;
	
  f[alpha_carbons[jm1]][0] -= -V*(an*prd_pair_theta[0]*xNO[0] + ah*prd_pair_theta[1]*xHO[0]);
  f[alpha_carbons[jm1]][1] -= -V*(an*prd_pair_theta[0]*xNO[1] + ah*prd_pair_theta[1]*xHO[1]);
  f[alpha_carbons[jm1]][2] -= -V*(an*prd_pair_theta[0]*xNO[2] + ah*prd_pair_theta[1]*xHO[2]);

  f[alpha_carbons[j]][0] -= -V*(bn*prd_pair_theta[0]*xNO[0] + bh*prd_pair_theta[1]*xHO[0]);
  f[alpha_carbons[j]][1] -= -V*(bn*prd_pair_theta[0]*xNO[1] + bh*prd_pair_theta[1]*xHO[1]);
  f[alpha_carbons[j]][2] -= -V*(bn*prd_pair_theta[0]*xNO[2] + bh*prd_pair_theta[1]*xHO[2]);

  f[oxygens[jm1]][0] -= -V*(cn*prd_pair_theta[0]*xNO[0] + ch*prd_pair_theta[1]*xHO[0]);
  f[oxygens[jm1]][1] -= -V*(cn*prd_pair_theta[0]*xNO[1] + ch*prd_pair_theta[1]*xHO[1]);
  f[oxygens[jm1]][2] -= -V*(cn*prd_pair_theta[0]*xNO[2] + ch*prd_pair_theta[1]*xHO[2]);

  f[oxygens[i]][0] -= V*(prd_pair_theta[0]*xNO[0] + prd_pair_theta[1]*xHO[0]);
  f[oxygens[i]][1] -= V*(prd_pair_theta[0]*xNO[1] + prd_pair_theta[1]*xHO[1]);
  f[oxygens[i]][2] -= V*(prd_pair_theta[0]*xNO[2] + prd_pair_theta[1]*xHO[2]);

  for (k=0;k<nn;++k) {
    if (res_info[k]!=LOCAL || res_info[k]!=GHOST) continue;
	
    k_resno = res_no[k]-1;
    k_chno = chain_no[k]-1;
    
    if (se[res_no[k]-1]=='G') { xk = xca[k]; katom = alpha_carbons[k]; }
    else { xk = xcb[k]; katom  = beta_atoms[k]; if(katom==-1) continue; }
		
    if (abs(k_resno-i_resno)>1 || k_chno!=i_chno) {
      dx[0] = xi[0] - xk[0];
      dx[1] = xi[1] - xk[1];
      dx[2] = xi[2] - xk[2];
			
      force = pair_theta_gamma*helix_well->prd_H(i)*helix_well->H(j)*helix_well->prd_theta(i, k, 0);
			
      f[iatom][0] -= force*dx[0];
      f[iatom][1] -= force*dx[1];
      f[iatom][2] -= force*dx[2];
			
      f[katom][0] -= -force*dx[0];
      f[katom][1] -= -force*dx[1];
      f[katom][2] -= -force*dx[2];
    }
    if (abs(k_resno-j_resno)>1 || k_chno!=j_chno) {
      dx[0] = xj[0] - xk[0];
      dx[1] = xj[1] - xk[1];
      dx[2] = xj[2] - xk[2];

      force = pair_theta_gamma*helix_well->H(i)*helix_well->prd_H(j)*helix_well->prd_theta(j, k, 0);
			
      f[jatom][0] -= force*dx[0];
      f[jatom][1] -= force*dx[1];
      f[jatom][2] -= force*dx[2];
			
      f[katom][0] -= -force*dx[0];
      f[katom][1] -= -force*dx[1];
      f[katom][2] -= -force*dx[2];
    }
  }
}


void FixBackbone::compute_amhgo_normalization()
{
  // compute normalization constant for the amhgo potential
  // the constant is called "a" and is given in Eqn. 8 in 
  // Eastwood and Wolynes 2000 "Role of explicitly..."
  // a = 1/(8N) \sum_i abs(\sum_(j in native contact) gamma_ij)^p
	
  int i, j, ich, res0, resn, iatom, jatom;
  int ires_type, jres_type;
  double amhgo_gamma, rnative;
  double normi;
	
  // Loop over chains
  for (ich=0;ich<nch;++ich) {
    res0 = ch_pos[ich]-1;
    resn = ch_pos[ich]+ch_len[ich]-1;
	
    // Double loop over all residue pairs
    for (i=res0;i<resn;++i) {
      ires_type = se_map[se[i]-'A'];
			
      for (iatom=Fragment_Memory::FM_CA; iatom<=Fragment_Memory::FM_CB - (se[i]=='G' ? 1 : 0); ++iatom) {
	normi = 0.0;
				
	for (j=res0;j<resn;++j) {
	  jres_type = se_map[se[j]-'A'];
					
	  for (jatom=Fragment_Memory::FM_CA; jatom<=Fragment_Memory::FM_CB - (se[j]=='G' ? 1 : 0); ++jatom) {
					
	    if (abs(i-j)<amh_go_gamma->minSep()) continue;

	    // if frustration censoring is on, check to see if interaction is censored
	    // if (frustration_censoring_flag == 1){
	    //   if (frustration_censoring_map[i][j] == 1 || frustration_censoring_map[j][i] == 1) continue;
	    // }
						
	    rnative = m_amh_go->Rf(i, iatom, j, jatom);
	    if (rnative<amh_go_rc) {
	      amhgo_gamma = amh_go_gamma->getGamma(ires_type, jres_type, i, j);
	      normi +=amhgo_gamma;
	    }
	  }
	}
	amh_go_norm[ich] += pow(fabs(normi), amh_go_p);
      }
    }
    amh_go_norm[ich] /= 8*resn;
  }	
}

void FixBackbone::compute_amh_go_model()
{
  int i, j, k, ii, jj, inum, jnum, ires, jres, iatom, jatom, ires_type, jres_type;
  int imol, jmol;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double xi[3], xj[3], dx[3], r, dr, drsq, rnative, amhgo_sigma_sq, amhgo_gamma;
  double Eij, Ei=0.0, E=0.0, force, factor;
  
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *mask = atom->mask;
  
  int nforces; // Number of atoms' forces in the buffer
    
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    ires = atom->residue[i];
    imol = atom->molecule[i];
    ires_type = se_map[se[ires-1]-'A'];
    
    // atom i is either C-Alpha or C-Bata and is LOCAL
    if ( (mask[i]&groupbit || (mask[i]&group2bit && se[ires-1]!='G') ) && i<nlocal ) {
      xi[0] = x[i][0];
      xi[1] = x[i][1];
      xi[2] = x[i][2];
      
      if (domain->xperiodic) xi[0] += prd[0]*((image[i] & 1023) - 512);
      if (domain->yperiodic) xi[1] += prd[1]*((image[i] >> 10 & 1023) - 512);
      if (domain->zperiodic) xi[2] += prd[2]*((image[i] >> 20) - 512);
      
      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      nforces = 1;
      amh_go_force_map[0] = i;
      amh_go_force[0][0] = amh_go_force[0][1] = amh_go_force[0][2] = 0.0;
      
      Ei = 0.0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        jres = atom->residue[j];
        jmol = atom->molecule[j];
        jres_type = se_map[se[jres-1]-'A'];

	// test to see if the interactions between ires and jres are censored
	if (frustration_censoring_flag == 1){
	  if (frustration_censoring_map[ires-1][jres-1] == 1 || frustration_censoring_map[jres-1][ires-1] == 1){
	    continue;
	  }
	}
        
        // atom j is either C-Alpha or C-Bata
        if ( (mask[j]&groupbit || (mask[j]&group2bit && se[jres-1]!='G') ) && abs(ires-jres)>=amh_go_gamma->minSep() && imol==jmol ) {
          xj[0] = x[j][0];
          xj[1] = x[j][1];
          xj[2] = x[j][2];
          
          if (domain->xperiodic) xj[0] += prd[0]*((image[j] & 1023) - 512);
          if (domain->yperiodic) xj[1] += prd[1]*((image[j] >> 10 & 1023) - 512);
          if (domain->zperiodic) xj[2] += prd[2]*((image[j] >> 20) - 512);

          if (mask[i]&groupbit) iatom = Fragment_Memory::FM_CA; else iatom = Fragment_Memory::FM_CB;
          if (mask[j]&groupbit) jatom = Fragment_Memory::FM_CA; else jatom = Fragment_Memory::FM_CB;
          rnative = m_amh_go->Rf(ires-1, iatom, jres-1, jatom);          

          if (rnative<amh_go_rc) {
	    dx[0] = xi[0] - xj[0];
            dx[1] = xi[1] - xj[1];
            dx[2] = xi[2] - xj[2];
          
            r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
            dr = r - rnative;
            drsq = dr*dr;
            
            amhgo_sigma_sq = pow(abs(ires-jres), 0.3);
            
            // drsq < 12*Log[10]*sigma_sq ~= 27.6*sigma_sq
            // this equivalent to having exp[-drsq/2*sigma_sq]=10^-6
            if (drsq<27.6*amhgo_sigma_sq) {
            
              amhgo_gamma = amh_go_gamma->getGamma(ires_type, jres_type, ires-1, jres-1);
              if (amh_go_gamma->error==amh_go_gamma->ERR_CALL) error->all(FLERR,"AMH-Go: Wrong call of getGamma() function");
			  
              Eij = amhgo_gamma*exp(-drsq/(2*amhgo_sigma_sq));
			  
              force = Eij*dr/(amhgo_sigma_sq*r);
			  
              amh_go_force[0][0] += force*dx[0];
              amh_go_force[0][1] += force*dx[1];
              amh_go_force[0][2] += force*dx[2];
			  
              amh_go_force[nforces][0] = -force*dx[0];
              amh_go_force[nforces][1] = -force*dx[1];
              amh_go_force[nforces][2] = -force*dx[2];
			  
              amh_go_force_map[nforces] = j;
              nforces++;
			  
              Ei += Eij;
            }
          }
        }
      }
      
      factor = -0.5*epsilon*k_amh_go*amh_go_p*pow(Ei, amh_go_p-1)/amh_go_norm[imol-1];
      for (k=0;k<nforces;k++) {
        f[amh_go_force_map[k]][0] += factor*amh_go_force[k][0];
        f[amh_go_force_map[k]][1] += factor*amh_go_force[k][1];
        f[amh_go_force_map[k]][2] += factor*amh_go_force[k][2];
      }
      
      E += -0.5*epsilon*k_amh_go*pow(Ei, amh_go_p)/amh_go_norm[imol-1];
    }
  }
  
  energy[ET_AMHGO] += E;
}

void FixBackbone::compute_vector_fragment_memory_potential(int i)
{
  int j, js, je, i_fm;
  int i_resno, j_resno, ires_type, jres_type;
  double vi[3], vj[3], vmi, vmj, vmsqi, vmsqj, vp, vpn, gc, gf, dg;
  double V, epsilon_k_weight, force, forcei[3], forcej[3];
  Fragment_Memory *frag;
  
  i_resno = res_no[i]-1;
  ires_type = se_map[se[i_resno]-'A'];
  
  for (i_fm=0; i_fm<ilen_fm_map[i_resno]; ++i_fm) {
    frag = frag_mems[ frag_mem_map[i_resno][i_fm] ];
    
    epsilon_k_weight = epsilon*k_vec_frag_mem;
    
    js = i+fm_gamma->minSep();
    je = MIN(frag->pos+frag->len-1, i+fm_gamma->maxSep());
    if (je>=n || res_no[je]-res_no[i]!=je-i) error->all(FLERR,"Missing residues in memory potential");
    
    for (j=js;j<=je;++j) {
      j_resno = res_no[j]-1;
      jres_type = se_map[se[j_resno]-'A'];
      
      if (chain_no[i]!=chain_no[j]) error->all(FLERR,"Fragment Memory: Interaction between residues of different chains");
      
      if (se[i_resno]!='G' && se[j_resno]!='G' && frag->getSe(i_resno)!='G' && frag->getSe(j_resno)!='G') {
	    vi[0] = xcb[i][0] - xca[i][0];
	    vi[1] = xcb[i][1] - xca[i][1];
	    vi[2] = xcb[i][2] - xca[i][2];
	    
	    vj[0] = xcb[j][0] - xca[j][0];
	    vj[1] = xcb[j][1] - xca[j][1];
	    vj[2] = xcb[j][2] - xca[j][2];
	    
	    vmsqi = vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2];
	    vmsqj = vj[0]*vj[0]+vj[1]*vj[1]+vj[2]*vj[2];
	    vmi = sqrt(vmsqi);
	    vmj = sqrt(vmsqj);
	    vp = vi[0]*vj[0]+vi[1]*vj[1]+vi[2]*vj[2];
	    
	    vpn = vp/(vmi*vmj);
	    gc = acos(vpn);
	    
	    gf = frag->VMf(i_resno, j_resno);
	    if (frag->error==frag->ERR_CALL || frag->error==frag->ERR_VFM_GLY)
	      error->all(FLERR,"Vector_Fragment_Memory: Wrong call of VMf() function");
	    
	    dg = gc - gf;
	    
	    V = -epsilon_k_weight*exp(-dg*dg/(2*vfm_sigma_sq));
	    
	    energy[ET_VFRAGMEM] += V;
	    
	    force = -V*dg/(vfm_sigma_sq*vmi*vmj*sqrt(1-vpn*vpn));
	    
	    forcei[0] = force*(vj[0]-vi[0]*vp/vmsqi);
	    forcei[1] = force*(vj[1]-vi[1]*vp/vmsqi);
	    forcei[2] = force*(vj[2]-vi[2]*vp/vmsqi);
	    
	    forcej[0] = force*(vi[0]-vj[0]*vp/vmsqj);
	    forcej[1] = force*(vi[1]-vj[1]*vp/vmsqj);
	    forcej[2] = force*(vi[2]-vj[2]*vp/vmsqj);
	    
	    f[alpha_carbons[i]][0] += -forcei[0];
	    f[alpha_carbons[i]][1] += -forcei[1];
	    f[alpha_carbons[i]][2] += -forcei[2];
	    
	    f[beta_atoms[i]][0] += forcei[0];
	    f[beta_atoms[i]][1] += forcei[1];
	    f[beta_atoms[i]][2] += forcei[2];
	    
	    f[alpha_carbons[j]][0] += -forcej[0];
	    f[alpha_carbons[j]][1] += -forcej[1];
	    f[alpha_carbons[j]][2] += -forcej[2];
	    
	    f[beta_atoms[j]][0] += forcej[0];
	    f[beta_atoms[j]][1] += forcej[1];
	    f[beta_atoms[j]][2] += forcej[2];
	  }
    }
  }
}

/*void FixBackbone::compute_vector_fragment_memory_table()
{
  int i, j, js, je, ir, i_fm, k, itb, iatom[4], jatom[4], iatom_type[4], jatom_type[4];
  int i_first_res, i_last_res, i_resno, j_resno, ires_type, jres_type;
  double r, rf, dr, drsq, V, force;
  double fm_sigma_sq, frag_mem_gamma, epsilon_k_weight, epsilon_k_weight_gamma;
  Fragment_Memory *frag;
  
  int j, js, je, i_fm;
  int i_resno, j_resno, ires_type, jres_type;
  double vi[3], vj[3], vmi, vmj, vmsqi, vmsqj, vp, vpn, gc, gf, dg;
  double V, epsilon_k_weight, force, forcei[3], forcej[3];
  Fragment_Memory *frag;
  
  for (i=0; i<n; ++i) {	  
    //	  i_resno = res_no[i]-1;
    i_resno = i;
    ires_type = se_map[se[i_resno]-'A'];
	  
    for (i_fm=0; i_fm<ilen_fm_map[i_resno]; ++i_fm) {
      frag = frag_mems[ frag_mem_map[i_resno][i_fm] ];
		
      epsilon_k_weight = epsilon*k_vec_frag_mem;
		
      js = i+fm_gamma->minSep();
      je = MIN(frag->pos+frag->len-1, i+fm_gamma->maxSep());
      //		if (je>=n || res_no[je]-res_no[i]!=je-i) error->all(FLERR,"Missing residues in Vector Fragment Memory potential");
      if (je>=n) error->all(FLERR,"Missing residues in memory potential");
		
      for (j=js;j<=je;++j) {
	//		  j_resno = res_no[j]-1;
	j_resno = j;
	jres_type = se_map[se[j_resno]-'A'];
		  
	//		  if (chain_no[i]!=chain_no[j]) error->all(FLERR,"Vector Fragment Memory: Interaction between residues of different chains");
	
	if (se[i_resno]!='G' && se[j_resno]!='G' && frag->getSe(i_resno)!='G' && frag->getSe(j_resno)!='G') {
		vi[0] = xcb[i][0] - xca[i][0];
	    vi[1] = xcb[i][1] - xca[i][1];
	    vi[2] = xcb[i][2] - xca[i][2];
	    
	    vj[0] = xcb[j][0] - xca[j][0];
	    vj[1] = xcb[j][1] - xca[j][1];
	    vj[2] = xcb[j][2] - xca[j][2];
	    
	    vmsqi = vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2];
	    vmsqj = vj[0]*vj[0]+vj[1]*vj[1]+vj[2]*vj[2];
	    vmi = sqrt(vmsqi);
	    vmj = sqrt(vmsqj);
	    vp = vi[0]*vj[0]+vi[1]*vj[1]+vi[2]*vj[2];
	    
	    vpn = vp/(vmi*vmj);
	    gc = acos(vpn);
	    
	    gf = frag->VMf(i_resno, j_resno);
	    if (frag->error==frag->ERR_CALL || frag->error==frag->ERR_VFM_GLY)
	      error->all(FLERR,"Vector_Fragment_Memory: Wrong call of VMf() function");
	    
	    dg = gc - gf;
	    
	    V = -epsilon_k_weight*exp(-dg*dg/(2*vfm_sigma_sq));
	    
	    force = -V*dg/(vfm_sigma_sq*vmi*vmj*sqrt(1-vpn*vpn));
	    
	    forcei[0] = force*(vj[0]-vi[0]*vp/vmsqi);
	    forcei[1] = force*(vj[1]-vi[1]*vp/vmsqi);
	    forcei[2] = force*(vj[2]-vi[2]*vp/vmsqi);
	    
	    forcej[0] = force*(vi[0]-vj[0]*vp/vmsqj);
	    forcej[1] = force*(vi[1]-vj[1]*vp/vmsqj);
	    forcej[2] = force*(vi[2]-vj[2]*vp/vmsqj);
		  	
	  	itb = 4*tb_nbrs*i + 4*(j-js) + k;
	  	if (!fm_table[itb])
	    	fm_table[itb] = new TBV[tb_size];
			
	  	for (ir=0;ir<tb_size;++ir) {
	    	r = tb_rmin + ir*tb_dr;
				
	    dr = r - rf;
	    drsq = dr*dr;
				
	    V = -epsilon_k_weight_gamma*exp(-drsq/(2*fm_sigma_sq));
				
	    fm_table[itb][ir].energy += V;
				
	    fm_table[itb][ir].force += V*dr/(fm_sigma_sq*r);
				
	    //				fm_table[i][j-js][k][ir].energy = -epsilon_k_weight_gamma*exp(-drsq/(2*fm_sigma_sq));
				
	    //				fm_table[i][j-js][k][ir].force = V*dr/(fm_sigma_sq*r);
	  }
	}
      }
    }
  }
}*/

void FixBackbone::compute_fragment_memory_potential(int i)
{
  int j, js, je, i_fm, k, iatom[4], jatom[4], iatom_type[4], jatom_type[4];
  int i_first_res, i_last_res, i_resno, j_resno, ires_type, jres_type;
  double *xi[4], *xj[4], dx[3], r, rf, dr, drsq, V, force;
  double fm_sigma_sq, frag_mem_gamma, epsilon_k_weight, epsilon_k_weight_gamma;
  Fragment_Memory *frag;
  
  iatom_type[0] = Fragment_Memory::FM_CA;
  iatom_type[1] = Fragment_Memory::FM_CA;
  iatom_type[2] = Fragment_Memory::FM_CB;
  iatom_type[3] = Fragment_Memory::FM_CB;
  
  jatom_type[0] = Fragment_Memory::FM_CA;
  jatom_type[1] = Fragment_Memory::FM_CB;
  jatom_type[2] = Fragment_Memory::FM_CA;
  jatom_type[3] = Fragment_Memory::FM_CB;
  
  xi[0] = xca[i];
  xi[1] = xca[i];
  xi[2] = xcb[i];
  xi[3] = xcb[i];
  
  iatom[0] = alpha_carbons[i];
  iatom[1] = alpha_carbons[i];
  iatom[2] = beta_atoms[i];
  iatom[3] = beta_atoms[i];
  
  i_resno = res_no[i]-1;
  ires_type = se_map[se[i_resno]-'A'];
  
  for (i_fm=0; i_fm<ilen_fm_map[i_resno]; ++i_fm) {
    frag = frag_mems[ frag_mem_map[i_resno][i_fm] ];
    
    epsilon_k_weight = epsilon*k_frag_mem*frag->weight;
    
    js = i+fm_gamma->minSep();
    je = MIN(frag->pos+frag->len-1, i+fm_gamma->maxSep());
    if (je>=n || res_no[je]-res_no[i]!=je-i) error->all(FLERR,"Missing residues in memory potential");
    
    for (j=js;j<=je;++j) {
      j_resno = res_no[j]-1;
      jres_type = se_map[se[j_resno]-'A'];
      
      if (chain_no[i]!=chain_no[j]) error->all(FLERR,"Fragment Memory: Interaction between residues of different chains");
      
      fm_sigma_sq = pow(abs(i_resno-j_resno), 2*fm_sigma_exp);
      
      if (!fm_gamma->fourResTypes()) {
	frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, i_resno, j_resno);
      } else {
	frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, frag->resType(i_resno), frag->resType(j_resno), i_resno, j_resno);
      }
      if (fm_gamma->error==fm_gamma->ERR_CALL) error->all(FLERR,"Fragment_Memory: Wrong call of getGamma() function");
      
      epsilon_k_weight_gamma = epsilon_k_weight*frag_mem_gamma;
      
      xj[0] = xca[j];
      xj[1] = xcb[j];
      xj[2] = xca[j];
      xj[3] = xcb[j];
      
      jatom[0] = alpha_carbons[j];
      jatom[1] = beta_atoms[j];
      jatom[2] = alpha_carbons[j];
      jatom[3] = beta_atoms[j];
      
      for (k=0;k<4;++k) {
	if ( iatom_type[k]==frag->FM_CB && (se[i_resno]=='G' || frag->getSe(i_resno)=='G') ) continue;
	if ( jatom_type[k]==frag->FM_CB && (se[j_resno]=='G' || frag->getSe(j_resno)=='G') ) continue;
        
        dx[0] = xi[k][0] - xj[k][0];
        dx[1] = xi[k][1] - xj[k][1];
        dx[2] = xi[k][2] - xj[k][2];

        r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
        rf = frag->Rf(i_resno, iatom_type[k], j_resno, jatom_type[k]);
        if (frag->error==frag->ERR_CALL) error->all(FLERR,"Fragment_Memory: Wrong call of Rf() function");
        dr = r - rf;
        drsq = dr*dr;
        
        V = -epsilon_k_weight_gamma*exp(-drsq/(2*fm_sigma_sq));
        
        energy[ET_FRAGMEM] += V;
        
        force = V*dr/(fm_sigma_sq*r);
        
        f[iatom[k]][0] += force*dx[0];
        f[iatom[k]][1] += force*dx[1];
        f[iatom[k]][2] += force*dx[2];
        
        f[jatom[k]][0] += -force*dx[0];
        f[jatom[k]][1] += -force*dx[1];
        f[jatom[k]][2] += -force*dx[2];
      }
    }
  }
}

void FixBackbone::table_fragment_memory(int i, int j)
{
  static double tmax=0.0;

  int k, i_resno, j_resno, tb_i, tb_j, itb, iatom_type[4], jatom_type[4], iatom[4], jatom[4], ir;
  double *xi[4], *xj[4], dx[3], r, r1, r2;
  double V, ff, v1, v2, f1, f2;

  i_resno = res_no[i]-1;
  j_resno = res_no[j]-1;
  
  if ( j_resno-i_resno<fm_gamma->minSep() ) return;
  if ( fm_gamma->maxSep()!=-1 && j_resno-i_resno>fm_gamma->maxSep() ) return;
  
  tb_i = i_resno;
  tb_j = j_resno - i_resno - fm_gamma->minSep();

  itb = 4*tb_nbrs*tb_i + 4*tb_j;
  if (!fm_table[itb]) return;
  
  iatom_type[0] = Fragment_Memory::FM_CA;
  iatom_type[1] = Fragment_Memory::FM_CA;
  iatom_type[2] = Fragment_Memory::FM_CB;
  iatom_type[3] = Fragment_Memory::FM_CB;
  
  jatom_type[0] = Fragment_Memory::FM_CA; 
  jatom_type[1] = Fragment_Memory::FM_CB; 
  jatom_type[2] = Fragment_Memory::FM_CA; 
  jatom_type[3] = Fragment_Memory::FM_CB;
  
  iatom[0] = alpha_carbons[i];
  iatom[1] = alpha_carbons[i];
  iatom[2] = beta_atoms[i];
  iatom[3] = beta_atoms[i];
  
  jatom[0] = alpha_carbons[j];
  jatom[1] = beta_atoms[j];
  jatom[2] = alpha_carbons[j];
  jatom[3] = beta_atoms[j];
  
  xi[0] = xca[i];
  xi[1] = xca[i];
  xi[2] = xcb[i];
  xi[3] = xcb[i];
  
  xj[0] = xca[j];
  xj[1] = xcb[j];
  xj[2] = xca[j];
  xj[3] = xcb[j];
  
  for (k=0;k<4;++k) {
  
    if (se[i_resno]=='G' && iatom_type[k]==Fragment_Memory::FM_CB) continue;
    if (se[j_resno]=='G' && jatom_type[k]==Fragment_Memory::FM_CB) continue;
    
    dx[0] = xi[k][0] - xj[k][0];
    dx[1] = xi[k][1] - xj[k][1];
    dx[2] = xi[k][2] - xj[k][2];

    r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    
    if (r>=tb_rmin && r<=tb_rmax) {
      ir = int((r-tb_rmin)/tb_dr);
    	
      itb = 4*tb_nbrs*tb_i + 4*tb_j + k;
    	
      if (!fm_table[itb]) return;
    	
      if (ir<0 || ir>=tb_size) error->all(FLERR,"Table Fragment Memory: ir is out of range.");
    	
      // Energy and force values are obtained from trangle interpolation
      r1 = tb_rmin + (double)ir*tb_dr;
      r2 = tb_rmin + (double)(ir+1)*tb_dr;
    	
      //    	v1 = fm_table[tb_i][tb_j][k][ir].energy;
      //    	v2 = fm_table[tb_i][tb_j][k][ir+1].energy;
      v1 = fm_table[itb][ir].energy;
      v2 = fm_table[itb][ir+1].energy;
    	
      V = ((v2-v1)*r + v1*r2 - v2*r1)/(r2-r1);
    	
      //    	f1 = fm_table[tb_i][tb_j][k][ir].force;
      //    	f2 = fm_table[tb_i][tb_j][k][ir+1].force;
      f1 = fm_table[itb][ir].force;
      f2 = fm_table[itb][ir+1].force;
    	
      ff = ((f2-f1)*r + f1*r2 - f2*r1)/(r2-r1);
    	   	
      energy[ET_FRAGMEM] += V;
    	
      f[iatom[k]][0] += ff*dx[0];
      f[iatom[k]][1] += ff*dx[1];
      f[iatom[k]][2] += ff*dx[2];
        
      f[jatom[k]][0] += -ff*dx[0];
      f[jatom[k]][1] += -ff*dx[1];
      f[jatom[k]][2] += -ff*dx[2];
    } else {
      fprintf(screen,  "r=%f, i=%d, j=%d\n", r, i, j);
      fprintf(logfile, "r=%f, i=%d, j=%d\n", r, i, j);
      error->all(FLERR,"Table Fragment Memory: r is out of computed range.");
    }	    
  }
}

void FixBackbone::compute_decoy_memory_potential(int i, int decoy_calc)
{
  // adapted from compute_fragment_memory_potential
  // given a residue i and a decoy_calc index, calculates a decoy memory energy
  // the energy is stored for residues i and j in the decoy_energy array
  // for decoy_calc > 0, the fragment memory energy is calculated using a shuffled fragment library (only used in "shuffle" mode)
  // for decoy_calc = 0, the fragment memory energy is calculated using the default fragment library 
  // the decoy_calc = 0 energy is used as the native energy in compute_fragment_frustration

  int j, js, je, i_fm, k, iatom[4], jatom[4], iatom_type[4], jatom_type[4];
  int i_resno, j_resno, ires_type, jres_type;
  double *xi[4], *xj[4], dx[3], r, rf, dr, drsq, V;
  double fm_sigma_sq, frag_mem_gamma, epsilon_k_weight, epsilon_k_weight_gamma, k_seqsep;
  Fragment_Memory *frag;
  int num_frags;
  
  iatom_type[0] = Fragment_Memory::FM_CA;
  iatom_type[1] = Fragment_Memory::FM_CA;
  iatom_type[2] = Fragment_Memory::FM_CB;
  iatom_type[3] = Fragment_Memory::FM_CB;
  
  jatom_type[0] = Fragment_Memory::FM_CA; 
  jatom_type[1] = Fragment_Memory::FM_CB; 
  jatom_type[2] = Fragment_Memory::FM_CA; 
  jatom_type[3] = Fragment_Memory::FM_CB;
  
  xi[0] = xca[i];
  xi[1] = xca[i];
  xi[2] = xcb[i];
  xi[3] = xcb[i];
  
  i_resno = res_no[i]-1;
  ires_type = se_map[se[i_resno]-'A'];
  
  // initialize num_frags as the number of fragments for a residue i
  if (decoy_calc == 0) 
    {
      num_frags=ilen_fm_map[i_resno];
    }
  else 
    {
      num_frags=ilen_decoy_map[i_resno];
    }
  
  // loop over fragments associated with residue i
  for (i_fm=0; i_fm<num_frags; ++i_fm) 
    {      
      
      // declare frag to be a fragment memory object
      if (decoy_calc == 0)
	{
	  frag = frag_mems[ frag_mem_map[i_resno][i_fm] ];
	}
      else
	{
	  frag = decoy_mems[ decoy_mem_map[i_resno][i_fm] ];
	}
      

      
      // loop over all residues j associated with residue i for fragment i_fm 
      js = i+fm_gamma->minSep();
      je = MIN(frag->pos+frag->len-1, i+fm_gamma->maxSep());
      
      epsilon_k_weight = epsilon*k_frag_mem*frag->weight;
      
      if (je>=n || res_no[je]-res_no[i]!=je-i) error->all(FLERR,"Missing residues in decoy memory potential");
      
      for (j=js;j<=je;++j) 
	{
	  j_resno = res_no[j]-1;
	  jres_type = se_map[se[j_resno]-'A'];
	  
	  if (chain_no[i]!=chain_no[j]) error->all(FLERR,"Decoy Memory: Interaction between residues of different chains");
	  
	  fm_sigma_sq = pow(abs(i_resno-j_resno), 2*fm_sigma_exp);
	  fm_sigma_sq = fm_sigma_sq*frag_frust_well_width*frag_frust_well_width;

	  if (!fm_gamma->fourResTypes()) 
	    {
	      frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, i_resno, j_resno);
	    } 
	  else 
	    {
	      frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, frag->resType(i_resno), frag->resType(j_resno), i_resno, j_resno);
	    }
	  if (fm_gamma->error==fm_gamma->ERR_CALL) error->all(FLERR,"Decoy_Memory: Wrong call of getGamma() function");
	  
	  // sequence distance dependent gamma
	  if (frag_frust_seqsep_flag)
	    {
	      k_seqsep = pow((abs(i_resno - j_resno)-fm_gamma->minSep()+1),-frag_frust_seqsep_gamma);
	      frag_mem_gamma *=k_seqsep;
	    }
	  
	  epsilon_k_weight_gamma = epsilon_k_weight*frag_mem_gamma;
	  
	  xj[0] = xca[j];
	  xj[1] = xcb[j];
	  xj[2] = xca[j];
	  xj[3] = xcb[j];
	  
	  // loop over combinations of CA, CB pairs
	  for (k=0;k<4;++k) {
	    if ( iatom_type[k]==frag->FM_CB && (se[i_resno]=='G' || frag->getSe(i_resno)=='G') ) continue;
	    if ( jatom_type[k]==frag->FM_CB && (se[j_resno]=='G' || frag->getSe(j_resno)=='G') ) continue;
	    
	    dx[0] = xi[k][0] - xj[k][0];
	    dx[1] = xi[k][1] - xj[k][1];
	    dx[2] = xi[k][2] - xj[k][2];
	    
	    r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
	    rf = frag->Rf(i_resno, iatom_type[k], j_resno, jatom_type[k]);
	    if (frag->error==frag->ERR_CALL) error->all(FLERR,"Fragment_Frustratometer: Wrong call of Rf() function");
	    dr = r - rf;
	    drsq = dr*dr;

	    if(frag_frust_normalizeInteraction)
	      {
		V *= 1/sqrt(fm_sigma_sq);
	      }
	    
	    V = -epsilon_k_weight_gamma*exp(-drsq/(2*fm_sigma_sq));
	    
	    // add decoy memory energy to both residue i and j 
	    decoy_energy[i_resno][decoy_calc] += V;
	    decoy_energy[j_resno][decoy_calc] += V;
	  }
	}
    }
}

// This function will shuffle the positions of the "decoy_mems" array
// This is used in "shuffle" mode to generate the decoy energies
void FixBackbone::randomize_decoys()
{
  //loops over decoy_mems and randomizes the starting position of each fragment object

  int i, k, pos, len, min_sep, random_position;

  for(i=0; i<n_decoy_mems; i++) 
    {
      // randomize the position of each decoy memory object such that the end of the fragment does not exceed the length of the protein
      random_position = rand() % (n-decoy_mems[i]->len+1);
      decoy_mems[i]->pos = random_position;
    }
  
  // repopulate decoy memory map
  for (i=0;i<n;++i) 
    {
      ilen_decoy_map[i] = 0;
      decoy_mem_map[i] = NULL;
    }
  
  min_sep = fm_gamma->minSep();
  
  for (k=0;k<n_decoy_mems;++k) 
    {
      pos = decoy_mems[k]->pos;
      len = decoy_mems[k]->len;
      
      if (pos+len>n) 
	{
	  fprintf(stderr, "pos %d len %d n %d\n", pos, len, n); 
	  error->all(FLERR,"Fragment_Frustratometer: Incorrectly defined memory fragment");
	}
      
      for (i=pos; i<pos+len-min_sep; ++i) 
	{
	  ilen_decoy_map[i]++;
	  decoy_mem_map[i] = (int *) memory->srealloc(decoy_mem_map[i],ilen_decoy_map[i]*sizeof(int),"modify:decoy_mem_map");
	  decoy_mem_map[i][ilen_decoy_map[i]-1] = k;
	}
    }
}

// This routine is used in both "shuffle" and "read" modes to compute the per residue
// fragment frustration. Because the number of decoy fragments can be larger or smaller
// than the number of "native" fragments in "shuffle" mode, a normalization is applied.
// In "read" mode, all of the decoy structures should be the same size as the target
// structure and uses the same number of fragments to calculate the energy, so no normalization
// is necessary.
void FixBackbone::compute_fragment_frustration()
{
  int residueindex, decoyindex;
  double averagedecoyenergy, variancedecoyenergy, frustrationindex;
  double nativeenergy;
  
  // if using shuffle method, normalize energies
  if (frag_frust_shuffle_flag)
    {
      // normalize native and decoy energies by number of respective memories
      for (residueindex=0; residueindex<n; residueindex++)
      	{
      	  decoy_energy[residueindex][0] /= n_frag_mems;
	  
      	  for (decoyindex=1;decoyindex<num_decoy_calcs;decoyindex++)
      	    {
      	      decoy_energy[residueindex][decoyindex] /= n_decoy_mems;
      	    }
      	}
    }
  
  // perform per residue statistics and calculate frustration index
  for (residueindex=0; residueindex<n; residueindex++)
    {
      // in "shuffle" mode, the mean and variance must be calculated from the decoy energy
      // distribution during every frustration calculation
      if (frag_frust_shuffle_flag)
	{
	  // zero per-residue variables
	  averagedecoyenergy=0.0;
	  variancedecoyenergy=0.0;
	  // i=0 is native, start at i=1
	  // compute average decoy energy
	  for (decoyindex=1;decoyindex<num_decoy_calcs;decoyindex++)
	    {
	      averagedecoyenergy += decoy_energy[residueindex][decoyindex];
	    }
	  // divide sum over decoys by num_decoy_calcs-1 (because native is excluded from sum)
	  averagedecoyenergy /= (num_decoy_calcs-1);
	  
	  // compute variance decoy energy
	  for (decoyindex=1;decoyindex<num_decoy_calcs;decoyindex++)
	    {
	      variancedecoyenergy += pow(decoy_energy[residueindex][decoyindex]-averagedecoyenergy,2);
	    }
	  variancedecoyenergy /= (num_decoy_calcs-1); 
	}
      // in "read" mode, the decoy energy distribution is computed once
      // and reused for every frustration calculation
      else if (frag_frust_read_flag)
	{
	  averagedecoyenergy = frag_frust_read_mean[residueindex];
	  variancedecoyenergy = frag_frust_read_variance[residueindex];
	}
      else
	{
	  // throw an error because only shuffle and read are valid modes
	  error->all(FLERR,"Fragment_Frustratometer: only shuffle and read are valid modes.");
	}
      
      // compute frustration index
      nativeenergy = decoy_energy[residueindex][0]; // the "native" energy is stored in index 0 of decoy_energy
      // the frustration index: f_i = E_i-<E_d>/(sqrt(/\E_d^2))
      frustrationindex = (nativeenergy-averagedecoyenergy)/(sqrt(variancedecoyenergy));
      // print out the frustration index to the fragment_frustration file
      fprintf(fragment_frustration_file, "%f ",frustrationindex);
      fprintf(fragment_frustration_gap_file, "%f ",(nativeenergy-averagedecoyenergy));
      fprintf(fragment_frustration_variance_file, "%f ",sqrt(variancedecoyenergy));
    }
  // print a new line for each new frustration calculation
  fprintf(fragment_frustration_file, "\n");
  fprintf(fragment_frustration_gap_file, "\n");
  fprintf(fragment_frustration_variance_file, "\n");
}

// This routine is used in "read" mode to calculate the decoy energy distribution
// It treats the decoy structures as if they were memories and uses the Rf() function
// to return the appropriate pairwise distances.
void FixBackbone::compute_generated_decoy_energies()
{
  int i, j, js, je, i_fm, k, iatom_type[4], jatom_type[4];
  int i_resno, j_resno, ires_type, jres_type;
  double r, rf, dr, drsq, V;
  double fm_sigma_sq, frag_mem_gamma, epsilon_k_weight, epsilon_k_weight_gamma, k_seqsep;
  Fragment_Memory *frag, *decoy;
  int num_frags;
	  
  int idecoy;
  int decoyindex;
  // for each generated decoy in the mem file
  for (idecoy=0; idecoy<n_decoy_mems; idecoy++)
    {
      // set decoy memory object
      decoy = decoy_mems[idecoy];

      // for each residue in the protein
      for (i=0;i<n;i++)
	{
	  i_resno = res_no[i]-1;
	  	  
	  iatom_type[0] = Fragment_Memory::FM_CA;
	  iatom_type[1] = Fragment_Memory::FM_CA;
	  iatom_type[2] = Fragment_Memory::FM_CB;
	  iatom_type[3] = Fragment_Memory::FM_CB;
	  
	  jatom_type[0] = Fragment_Memory::FM_CA; 
	  jatom_type[1] = Fragment_Memory::FM_CB; 
	  jatom_type[2] = Fragment_Memory::FM_CA; 
	  jatom_type[3] = Fragment_Memory::FM_CB;

	  i_resno = res_no[i]-1;
	  ires_type = se_map[se[i_resno]-'A'];

	  // set number of frags for residue i_resno
	  num_frags = ilen_fm_map[i_resno];

	  // loop over fragments associated with residue i
	  for (i_fm=0; i_fm<num_frags; ++i_fm) 
	    {      
	      frag = frag_mems[ frag_mem_map[i_resno][i_fm] ];
	      
	      epsilon_k_weight = epsilon*k_frag_mem*frag->weight;
	      
	      // loop over all residues j associated with residue i for fragment i_fm 
	      js = i+fm_gamma->minSep();
	      je = MIN(frag->pos+frag->len-1, i+fm_gamma->maxSep());
	      
	      if (je>=n || res_no[je]-res_no[i]!=je-i) error->all(FLERR,"Missing residues in decoy memory potential");
	      
	      for (j=js;j<=je;++j) 
		{
		  j_resno = res_no[j]-1;
		  jres_type = se_map[se[j_resno]-'A'];
		  
		  fm_sigma_sq = pow(abs(i_resno-j_resno), 2*fm_sigma_exp);
		  fm_sigma_sq = fm_sigma_sq*frag_frust_well_width*frag_frust_well_width;
		  if (!fm_gamma->fourResTypes()) 
		    {
		      frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, i_resno, j_resno);
		    } 
		  else 
		    {
		      frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, frag->resType(i_resno), frag->resType(j_resno), i_resno, j_resno);
		    }
		  if (fm_gamma->error==fm_gamma->ERR_CALL) error->all(FLERR,"Decoy_Memory: Wrong call of getGamma() function");
		  
		  // sequence distance dependent gamma
		  if (frag_frust_seqsep_flag)
		    {
		      k_seqsep = pow((abs(i_resno - j_resno)-fm_gamma->minSep()+1),-frag_frust_seqsep_gamma);
		      frag_mem_gamma *=k_seqsep;
		    }

		  epsilon_k_weight_gamma = epsilon_k_weight*frag_mem_gamma;
		  
		  // loop over combinations of CA, CB pairs
		  for (k=0;k<4;++k) 
		    {
		      if ( iatom_type[k]==frag->FM_CB && (se[i_resno]=='G' || frag->getSe(i_resno)=='G') ) continue;
		      if ( jatom_type[k]==frag->FM_CB && (se[j_resno]=='G' || frag->getSe(j_resno)=='G') ) continue;

		      r = decoy->Rf(i_resno, iatom_type[k], j_resno, jatom_type[k]);
		      rf = frag->Rf(i_resno, iatom_type[k], j_resno, jatom_type[k]);

		      if (frag->error==frag->ERR_CALL) error->all(FLERR,"Fragment_Frustratometer: Wrong call of Rf() function");
		      dr = r - rf;
		      drsq = dr*dr;
		      
		      if(frag_frust_normalizeInteraction)
			{
			  V *= 1/sqrt(fm_sigma_sq);
			}

		      V = -epsilon_k_weight_gamma*exp(-drsq/(2*fm_sigma_sq));
		      
		      // add decoy memory energy to both residue i and j 
		      decoy_energy[i_resno][idecoy+1] += V; // shift all decoy energies by one
		      decoy_energy[j_resno][idecoy+1] += V; // so that the native energy is in index 0
		    }
		}
	    }
	}
    }
  // compute statistics on the decoy_energy array
  for (i=0; i<n; i++)
    {
      fprintf(fragment_frustration_native_data,"%f ",decoy_energy[i][0]);
    }

  // for every residue, average over idecoy to get frag_frust_read_mean and frag_frust_read_variance
  for (i=0; i<n; i++)
    {
      // zero per-residue variables
      frag_frust_read_mean[i]=0.0;
      frag_frust_read_variance[i]=0.0;
      // i=0 is native, start at i=1
      // compute average decoy energy
      // note that num_decoy_calcs is the size of the decoy_energy array, which in this case
      // is one larger than n_decoy_mems because the native energy is stored in index 0
      // which is why the loop starts with index 1
      for (decoyindex=1;decoyindex<num_decoy_calcs;decoyindex++)
	{
	  fprintf(fragment_frustration_decoy_data," %f\n",decoy_energy[i][decoyindex]);
	  frag_frust_read_mean[i] += decoy_energy[i][decoyindex];
	}
      // divide sum over decoys by n_decoy_mems
      frag_frust_read_mean[i] /= n_decoy_mems;
      
      // compute variance decoy energy
      // note that num_decoy_calcs is the size of the decoy_energy array, which in this case
      // is one larger than n_decoy_mems because the native energy is stored in index 0
      // which is why the loop starts with index 1
      for (decoyindex=1;decoyindex<num_decoy_calcs;decoyindex++)
	{
	  frag_frust_read_variance[i] += pow(decoy_energy[i][decoyindex]-frag_frust_read_mean[i],2);
	}
      frag_frust_read_variance[i] /= n_decoy_mems;      
    }
}

void FixBackbone::output_selection_temperature_data()
{
  double *xi, *xj, dx[3];
  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rij, rho_i, rho_j;
  double water_energy, burial_energy_i, burial_energy_j;
  
  if (selection_temperature_output_interaction_energies_flag) {
    // Loop over original sequence and output detailed information
    // Double loop over all residue pairs
    for (i=0;i<n;++i) {
      // get information about residue i
      i_resno = res_no[i]-1;
      ires_type = get_residue_type(i_resno);
      i_chno = chain_no[i]-1;
      for (j=i+1;j<n;++j) {
	// get information about residue j
	j_resno = res_no[j]-1;
	jres_type = get_residue_type(j_resno);
	j_chno = chain_no[j]-1;
	// get the distance between i and j
	rij = get_residue_distance(i_resno, j_resno);
	rho_i = get_residue_density(i_resno);
	rho_j = get_residue_density(j_resno);
	// compute the energies for the (i,j) pair
	water_energy = 0.0;
	if ((abs(i-j)>=contact_cutoff || i_chno != j_chno)) {
	  water_energy = compute_water_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
	}
	
	burial_energy_i = compute_burial_energy(i_resno, ires_type, rho_i);
	burial_energy_j = compute_burial_energy(j_resno, jres_type, rho_j);
	
	fprintf(selection_temperature_file,"%d %d %c %c %f %f %f %f %f %f\n", i+1, j+1, se[i], se[j], rij, rho_i, rho_j, water_energy, burial_energy_i, burial_energy_j);
      }
    }
  }

  if (selection_temperature_output_contact_list_flag) {
    fprintf(selection_temperature_contact_list_file,"# timestep: %d\n", ntimestep);
    // Loop over all pairs of residues, output those in contact
    for (i=0;i<n;++i) {
      // get information about residue i
      i_resno = res_no[i]-1;
      ires_type = get_residue_type(i_resno);
      i_chno = chain_no[i]-1;
      for (j=i+1;j<n;++j) {
	// get information about residue j
	j_resno = res_no[j]-1;
	jres_type = get_residue_type(j_resno);
	j_chno = chain_no[j]-1;
	// get the distance between i and j
	rij = get_residue_distance(i_resno, j_resno);
	if ((abs(i-j)>=selection_temperature_min_seq_sep || i_chno != j_chno)) {
	  if (rij < selection_temperature_rij_cutoff) {
	    fprintf(selection_temperature_contact_list_file,"%d %d\n", i+1, j+1);
	  }
	}
      }
    }
  }

  if (selection_temperature_evaluate_sequence_energies_flag) {
    // Loop over sequences in selection temperature sequences file and output energy
    // Sum the energies only from those residues in the selection temperature residues file
    int i_sel_temp = 0;
    int j_sel_temp = 0;
    double temp_sequence_energy = 0.0;
    for (int i_sequence = 0; i_sequence<num_selection_temperature_sequences; i_sequence++) {
      //printf("%d\n", i_sequence);
      temp_sequence_energy = 0.0;
      i_sel_temp = 0;
      for (i=0;i<n;++i) {
	j_sel_temp = i_sel_temp+1;
	// get information about residue i
	i_resno = res_no[i]-1;
	if (!(i_resno == selection_temperature_residues[i_sel_temp]-1) || selection_temperature_sequences[i_sequence][i] == '*') {
	  continue;
	}
	i_sel_temp++;
	ires_type = se_map[selection_temperature_sequences[i_sequence][i]-'A'];
	i_chno = chain_no[i]-1;
	rho_i = get_residue_density(i_resno);
	temp_sequence_energy += compute_burial_energy(i_resno, ires_type, rho_i);
	for (j=i+1;j<n;++j) {
	  // get information about residue j
	  j_resno = res_no[j]-1;
	  if (!(j_resno == selection_temperature_residues[j_sel_temp]-1) || selection_temperature_sequences[i_sequence][j] == '*') {
	    continue;
	  }
	  j_sel_temp++;
	  jres_type = se_map[selection_temperature_sequences[i_sequence][j]-'A'];
	  j_chno = chain_no[j]-1;
	  // get the distance between i and j
	  rij = get_residue_distance(i_resno, j_resno);
	  rho_j = get_residue_density(j_resno);
	  // compute the energies for the (i,j) pair
	  water_energy = 0.0;
	  if ((abs(i-j)>=contact_cutoff || i_chno != j_chno)) {
	    water_energy = compute_water_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
	  }
	  temp_sequence_energy += water_energy;
	}
      }
      fprintf(selection_temperature_sequence_energies_output_file, "%f\n", temp_sequence_energy);
    }
  }
}

void FixBackbone::compute_mcso()
{
  int i;
  int rand_res_1, rand_res_2;
  char temp_res;

  double total_energy = 0.0;
  double new_total_energy = 0.0;
  double energy_difference = 0.0;
  double mcso_temp = mcso_start_temp;
  int mcso_temp_step;
  double random_probability = 0.0;
  double mcso_increment = ((mcso_end_temp-mcso_start_temp)/(double)mcso_num_steps);

  for (mcso_temp_step=0; mcso_temp_step<mcso_num_steps;mcso_temp_step++) {
    // compute total energy
    total_energy = compute_total_burial_energy() + compute_total_contact_energy();
    // printf("\n inital total energy %f\n", total_energy);
    // save old sequence
    for (i=0;i<n;i++) {
      mcso_se[i] = se[i];
    }
    // permute two residues
    rand_res_1 = rand() % n;
    rand_res_2 = rand() % n;
    temp_res = se[rand_res_1];
    se[rand_res_1] = se[rand_res_2];
    se[rand_res_2] = temp_res;
    // compute new total energy
    new_total_energy = compute_total_burial_energy() + compute_total_contact_energy();
    // printf("\n new total energy %f\n", new_total_energy);
    // accept or reject based on temperature
    energy_difference = new_total_energy - total_energy;
    // printf("\n energy difference %f\n", energy_difference);
    mcso_temp += mcso_increment;
    if (energy_difference > 0) {
      random_probability = (double)rand()/RAND_MAX;
      // printf("\n random probability %f\n", random_probability);
      // printf("\n exponential factor %f\n", exp(-energy_difference/mcso_temp));
      // printf("\n temperature %f\n", mcso_temp);
      if (random_probability > exp(-energy_difference/(k_b*mcso_temp))) {
	// printf("\n rejected!\n");
	  // if reject, put the old sequence back
	  for (i=0;i<n;i++) {
	    se[i] = mcso_se[i];
	  }	
	}
      else {
	// printf("\n accepted!\n");
	total_energy = new_total_energy;
      }
    }
    // output the sequence and energy
    for (i=0;i<n;i++) {
      fprintf(mcso_seq_output_file,"%c", se[i]);
    }	
    fprintf(mcso_seq_output_file,"\n");
    fprintf(mcso_energy_output_file,"%f\n",total_energy);
  }
}

double FixBackbone::compute_total_burial_energy()
{
  int i;
  int ires_type, i_resno, i_chno;
  double rho_i;

  double total_burial_energy = 0.0;

  for (i=0;i<n;++i) {
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = get_residue_type(i_resno);
    i_chno = chain_no[i]-1;
    rho_i = get_residue_density(i_resno);
    total_burial_energy += compute_burial_energy(i_resno, ires_type, rho_i);
  }

  return total_burial_energy;
}

double FixBackbone::compute_total_contact_energy()
{
  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rij, rho_i, rho_j;
  double total_water_energy = 0.0;
  
  // Double loop over all residue pairs
  for (i=0;i<n;++i) {
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = get_residue_type(i_resno);
    i_chno = chain_no[i]-1;
    for (j=i+1;j<n;++j) {
      // get information about residue j
      j_resno = res_no[j]-1;
      jres_type = get_residue_type(j_resno);
      j_chno = chain_no[j]-1;
      // get the distance between i and j
      rij = get_residue_distance(i_resno, j_resno);
      rho_i = get_residue_density(i_resno);
      rho_j = get_residue_density(j_resno);
      // compute the energies for the (i,j) pair
      if ((abs(i-j)>=contact_cutoff || i_chno != j_chno)) {
	total_water_energy += compute_water_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
      }
    }
  }

  return total_water_energy;
}

void FixBackbone::compute_tert_frust()
{
  double *xi, *xj, dx[3];
  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rij, rho_i, rho_j;
  double native_energy;
  double frustration_index;
  int atomselect;

  atomselect = 0; // for the vmd script output

  // Double loop over all residue pairs
  for (i=0;i<n;++i) {
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = get_residue_type(i_resno);
    i_chno = chain_no[i]-1;

    for (j=i+1;j<n;++j) {
      // get information about residue j
      j_resno = res_no[j]-1;
      jres_type = get_residue_type(j_resno);
      j_chno = chain_no[j]-1;
      
      // get the distance between i and j
      rij = get_residue_distance(i_resno, j_resno);

      // if the atoms are within the threshold, compute the frustration
      if (rij < tert_frust_cutoff && (abs(i-j)>=contact_cutoff || i_chno != j_chno)) {
	// Get coordinates of relevant atoms
	if (se[i_resno]=='G') { xi = xca[i]; }
	else { xi = xcb[i]; }
	if (se[j_resno]=='G') { xj = xca[j]; }
	else { xj = xcb[j]; }
	rho_i = get_residue_density(i_resno);
	rho_j = get_residue_density(j_resno);
	native_energy = compute_native_ixn(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
	// if the mode is not configurational or if the mode is configurational and we have
	// not already computed the decoy energies, compute the decoy energies
	// the configurational decoy statistics only need to be computed once because they are
	// the same for every contact
	if (!strcmp(tert_frust_mode, "configurational")==0 || (strcmp(tert_frust_mode, "configurational")==0 && !already_computed_configurational_decoys)) {
	  compute_decoy_ixns(i_resno, j_resno, rij, rho_i, rho_j);
	}
	frustration_index = compute_frustration_index(native_energy, decoy_ixn_stats);
	// write information out to output file
 	fprintf(tert_frust_output_file,"%5d %5d %3d %3d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %c %c %8.3f %8.3f %8.3f %8.3f\n", i_resno+1, j_resno+1, i_chno+1, j_chno+1, xi[0], xi[1], xi[2], xj[0], xj[1], xj[2], rij, rho_i, rho_j, se[i_resno], se[j_resno], native_energy, decoy_ixn_stats[0], decoy_ixn_stats[1], frustration_index);
	if(frustration_index > 0.78 || frustration_index < -1) {
	  // write information out to vmd script
	  fprintf(tert_frust_vmd_script,"set sel%d [atomselect top \"resid %d and name CA\"]\n", i_resno, i_resno+1);
	  fprintf(tert_frust_vmd_script,"set sel%d [atomselect top \"resid %d and name CA\"]\n", j_resno, j_resno+1);
	  fprintf(tert_frust_vmd_script,"lassign [atomselect%d get {x y z}] pos1\n",atomselect);
	  atomselect += 1;
	  fprintf(tert_frust_vmd_script,"lassign [atomselect%d get {x y z}] pos2\n",atomselect);
	  atomselect += 1;
	  if(frustration_index > 0.78) {
	    fprintf(tert_frust_vmd_script,"draw color green\n");
	  }
	  else {
	    fprintf(tert_frust_vmd_script,"draw color red\n");
	  }
	  if(rij < well->par.well_r_max[0]) {
	    fprintf(tert_frust_vmd_script,"draw line $pos1 $pos2 style solid width 1\n");
	  }
	  else {
	    fprintf(tert_frust_vmd_script,"draw line $pos1 $pos2 style dashed width 2\n");
	  }
	}
      }
    }
  }
  
  // after looping over all pairs, write out the end of the vmd script
  fprintf(tert_frust_vmd_script, "mol modselect 0 top \"all\"\n");
  fprintf(tert_frust_vmd_script, "mol modstyle 0 top newcartoon\n");
  fprintf(tert_frust_vmd_script, "mol modcolor 0 top colorid 15\n");
}

void FixBackbone::compute_tert_frust_singleresidue()
{
  int i;
  int ires_type, i_resno, i_chno;
  double rho_i;
  double native_energy;
  double frustration_index;
  int atomselect;
  double *xi;

  atomselect = 0; // for the vmd script output

  // Loop over all residues 
  for (i=0;i<n;++i) {
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = get_residue_type(i_resno);
    rho_i = get_residue_density(i_resno);
    i_chno = chain_no[i]-1;
    if (se[i_resno]=='G') { xi = xca[i]; }
    else { xi = xcb[i]; }

    // compute native energy
    native_energy = compute_singleresidue_native_ixn(i_resno, ires_type, rho_i, i_chno, tert_frust_cutoff, 0);

    // compute decoy energies
    compute_singleresidue_decoy_ixns(i_resno, rho_i, i_chno);

    // compute frustration index
    frustration_index = compute_frustration_index(native_energy, decoy_ixn_stats);

    // write information out to output file
    fprintf(tert_frust_output_file,"%5d %5d %8.3f %8.3f %8.3f %8.3f %c %8.3f %8.3f %8.3f %8.3f\n", i_resno+1, i_chno+1, xi[0], xi[1], xi[2], rho_i, se[i_resno], native_energy, decoy_ixn_stats[0], decoy_ixn_stats[1], frustration_index);
    if(frustration_index > 0.78 || frustration_index < -1) {
      // write information out to vmd script
      atomselect += 1;
      fprintf(tert_frust_vmd_script,"mol addrep 0\n");
      fprintf(tert_frust_vmd_script,"mol modselect %d 0 resid %d\n", atomselect, i_resno+1);
      fprintf(tert_frust_vmd_script,"mol modstyle %d 0 VDW %f 12.000000\n", atomselect, 0.5*abs(frustration_index));
      fprintf(tert_frust_vmd_script,"mol modmaterial %d 0 Transparent\n", atomselect);
      if(frustration_index > 0.78) {
    	// color the residue green\n
        fprintf(tert_frust_vmd_script,"mol modcolor %d 0 ColorID 7\n", atomselect);
      }
      else {
    	// color the residue red\n
	fprintf(tert_frust_vmd_script,"mol modcolor %d 0 ColorID 1\n", atomselect);
      }
    }
  }
  
  // after looping over all pairs, write out the end of the vmd script
  fprintf(tert_frust_vmd_script, "mol modselect 0 top \"all\"\n");
  fprintf(tert_frust_vmd_script, "mol modstyle 0 top newcartoon\n");
  fprintf(tert_frust_vmd_script, "mol modcolor 0 top colorid 15\n");
}

double FixBackbone::compute_native_ixn(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j)
{
  double water_energy, burial_energy_i, burial_energy_j, rik, rjk, rho_k;
  int k, kres_type;

  // compute the energies for the (i,j) pair
  water_energy = compute_water_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
  burial_energy_i = compute_burial_energy(i_resno, ires_type, rho_i);
  burial_energy_j = compute_burial_energy(j_resno, jres_type, rho_j);

  // in configurational mode, only the (i,j) contact contributes to the native energy
  // so we are ready to return the sum of water_ij+burial_i+burial_j
  if (strcmp(tert_frust_mode, "configurational")==0) {
    return water_energy + burial_energy_i + burial_energy_j;
  }
  // in mutational mode, all (i,k) and (j,k) pairs also contribute to the native energy
  // so we have to compute those first, then return the sum
  else if (strcmp(tert_frust_mode, "mutational")==0) {
    for (k=0; k<n; k++) {
      // make sure to exclude the native interactions because you have already counted those
      if (k==i_resno || k==j_resno) {
	continue;
      }
      // get interaction parameters for resdiue k
      rho_k = get_residue_density(k);
      kres_type = get_residue_type(k);

      // check to see if i and k are in contact; if so, add the energy
      rik = get_residue_distance(i_resno, k);
      if (rik < tert_frust_cutoff) {
	// add (i,k) contribution
	water_energy += compute_water_energy(rik, i_resno, k, ires_type, kres_type, rho_i, rho_k);
      }
      // check to see if j and k are in contact; if so, add the energy
      rjk = get_residue_distance(j_resno, k);
      if (rjk < tert_frust_cutoff) {	
	// add (j,k) contribution
	water_energy += compute_water_energy(rjk, j_resno, k, jres_type, kres_type, rho_j, rho_k);
      }
    }
    return water_energy+burial_energy_i+burial_energy_j;
  }
}

void FixBackbone::compute_decoy_ixns(int i_resno, int j_resno, double rij_orig, double rho_i_orig, double rho_j_orig)
{
  int decoy_i, rand_i_resno, rand_j_resno, ires_type, jres_type, k, kres_type;
  double rij, rho_i, rho_j, water_energy, burial_energy_i, burial_energy_j, rik, rjk, rho_k;
  
  for (decoy_i=0; decoy_i<tert_frust_ndecoys; decoy_i++) {
    if (strcmp(tert_frust_mode, "configurational")==0) {
      // choose random rij, rho_i, rho_j 
      rand_i_resno = get_random_residue_index();
      rand_j_resno = get_random_residue_index();
      rij = get_residue_distance(rand_i_resno, rand_j_resno);
      // make sure that the randomly chosen residues are in contact
      while(rij > tert_frust_cutoff) { 
	rand_i_resno = get_random_residue_index();
	rand_j_resno = get_random_residue_index();
	rij = get_residue_distance(rand_i_resno, rand_j_resno);
      }
      // get new pair of random residues for burial term
      rand_i_resno = get_random_residue_index();
      rand_j_resno = get_random_residue_index();
      rho_i = get_residue_density(rand_i_resno);
      rho_j = get_residue_density(rand_j_resno);
    }
    else {
      // if in mutational mode, use configurational parameters passed into the function
      rij = rij_orig;
      rho_i = rho_i_orig;
      rho_j = rho_j_orig;
    }

    // choose random ires_type, jres_type
    rand_i_resno = get_random_residue_index();
    rand_j_resno = get_random_residue_index();
    ires_type = get_residue_type(rand_i_resno);
    jres_type = get_residue_type(rand_j_resno);

    // compute energy terms for the (i,j) pair
    water_energy = compute_water_energy(rij, rand_i_resno, rand_j_resno, ires_type, jres_type, rho_i, rho_j);
    burial_energy_i = compute_burial_energy(rand_i_resno, ires_type, rho_i);
    burial_energy_j = compute_burial_energy(rand_j_resno, jres_type, rho_j);

    // in mutational mode, all (i,k) and (j,k) pairs also contribute to the decoy energy
    // so we have to compute those first, then return the sum
    if (strcmp(tert_frust_mode, "mutational")==0) {
      for (k=0; k<n; k++) {
	// make sure to exclude the native interactions because you have already counted those
	if (k==i_resno || k==j_resno) {
	  continue;
	}
	// get interaction parameters for resdiue k
	rho_k = get_residue_density(k);
	kres_type = get_residue_type(k);
	
	// check to see if i and k are in contact; if so, add the energy
	rik = get_residue_distance(i_resno, k);
	if (rik < tert_frust_cutoff) {
	  // add (i,k) contribution
	  water_energy += compute_water_energy(rik, rand_i_resno, k, ires_type, kres_type, rho_i, rho_k);
	}
	// check to see if j and k are in contact; if so, add the energy
	rjk = get_residue_distance(j_resno, k);
	if (rjk < tert_frust_cutoff) {	
	  // add (j,k) contribution
	  water_energy += compute_water_energy(rjk, rand_j_resno, k, jres_type, kres_type, rho_j, rho_k);
	}
      }
    }

    // sum the energy terms, store in array
    tert_frust_decoy_energies[decoy_i] = water_energy + burial_energy_i + burial_energy_j;
  }

  // save the mean and standard deviation into the decoy_ixn_stats array
  decoy_ixn_stats[0] = compute_array_mean(tert_frust_decoy_energies, tert_frust_ndecoys);
  decoy_ixn_stats[1] = compute_array_std(tert_frust_decoy_energies, tert_frust_ndecoys);

  // if we get here and the mode is configurational, we have already computed the statistics
  // for the decoys and we don't need to repeat the calculation again
  // this is to save time when using the configurational mode of the frustratometer
  if (strcmp(tert_frust_mode, "configurational")==0) {
    already_computed_configurational_decoys = 1;
  }
}

double FixBackbone::compute_singleresidue_native_ixn(int i_resno, int ires_type, double rho_i, int i_chno, double cutoff, bool nmercalc)
{
  double water_energy, burial_energy_i, rij, rho_j;
  int j, j_resno, jres_type, j_chno;

  // find burial energy for residue i
  burial_energy_i = compute_burial_energy(i_resno, ires_type, rho_i);

  // initialize water energy
  water_energy = 0.0;

  // loop over all possible interacting residues
  for (j=0; j<n; j++) {
    // get information about residue j
    j_resno = res_no[j]-1;
    jres_type = get_residue_type(j_resno);
    j_chno = chain_no[j]-1;
    rho_j = get_residue_density(j_resno);

    // don't interact with self
    if (i_resno == j_resno) {
      continue;
    }

    // don't double count water energy for nmer calculations
    if (j_resno > i_resno && nmercalc) {
      continue;
    }

    // find distance between residues i and j
    rij = get_residue_distance(i_resno, j_resno);
    
    // if within the interaction distance, compute energy
    if (rij < cutoff && (abs(i_resno-j_resno)>=contact_cutoff || i_chno != j_chno)) {
      // compute the energies for the (i,j) pair
      water_energy += compute_water_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
    }
  }

  return water_energy + burial_energy_i;
}

void FixBackbone::compute_singleresidue_decoy_ixns(int i_resno, double rho_i, int i_chno)
{
  int decoy_i, rand_i_resno, ires_type;
  
  for (decoy_i=0; decoy_i<tert_frust_ndecoys; decoy_i++) {
    // randomize ires_type
    rand_i_resno = get_random_residue_index();
    ires_type = get_residue_type(rand_i_resno);

    // compute the decoy energy
    tert_frust_decoy_energies[decoy_i] = compute_singleresidue_native_ixn(i_resno, ires_type, rho_i, i_chno, tert_frust_cutoff, 0);
  }

  // save the mean and standard deviation into the decoy_ixn_stats array
  decoy_ixn_stats[0] = compute_array_mean(tert_frust_decoy_energies, tert_frust_ndecoys);
  decoy_ixn_stats[1] = compute_array_std(tert_frust_decoy_energies, tert_frust_ndecoys);
}

double FixBackbone::compute_array_mean(double *array, int arraysize)
{
  double mean;
  int i;

  mean = 0;
  for(i = 0; i < arraysize; i++) {
    mean += array[i];
  }
  mean /= (double)arraysize;

  return mean;
}

double FixBackbone::compute_array_std(double *array, int arraysize)
{
  double mean, std;
  int i;

  mean = compute_array_mean(array, arraysize);

  std = 0.0;
  for(i = 0; i < arraysize; i++) {
    std += (array[i] - mean)*(array[i] - mean);
  }
  std /= (double)arraysize;
  std = sqrt(std);

  return std;
}

double FixBackbone::compute_water_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j)
{
  double water_gamma_0_direct, water_gamma_1_direct, water_gamma_prot_mediated, water_gamma_wat_mediated;
  double sigma_wat, sigma_prot, sigma_gamma_direct, sigma_gamma_mediated;
  double t_min_direct, t_max_direct, theta_direct, t_min_mediated, t_max_mediated, theta_mediated;
  double water_energy;

  // if(abs(i_resno-j_resno)<contact_cutoff) return 0.0;

  water_gamma_0_direct = get_water_gamma(i_resno, j_resno, 0, ires_type, jres_type, 0);
  water_gamma_1_direct = get_water_gamma(i_resno, j_resno, 0, ires_type, jres_type, 1);

  water_gamma_prot_mediated = get_water_gamma(i_resno, j_resno, 1, ires_type, jres_type, 0);
  water_gamma_wat_mediated = get_water_gamma(i_resno, j_resno, 1, ires_type, jres_type, 1);

  sigma_wat = 0.25*(1.0 - tanh(well->par.kappa_sigma*(rho_i-well->par.treshold)))*(1.0 - tanh(well->par.kappa_sigma*(rho_j-well->par.treshold)));
  sigma_prot = 1.0 - sigma_wat;
  
  sigma_gamma_direct = (water_gamma_0_direct + water_gamma_1_direct)/2;
  sigma_gamma_mediated = sigma_prot*water_gamma_prot_mediated + sigma_wat*water_gamma_wat_mediated;

  t_min_direct = tanh( well->par.kappa*(rij - well->par.well_r_min[0]) );
  t_max_direct = tanh( well->par.kappa*(well->par.well_r_max[0] - rij) );
  theta_direct = 0.25*(1.0 + t_min_direct)*(1.0 + t_max_direct);

  t_min_mediated = tanh( well->par.kappa*(rij - well->par.well_r_min[1]) );
  t_max_mediated = tanh( well->par.kappa*(well->par.well_r_max[1] - rij) );
  theta_mediated = 0.25*(1.0 + t_min_mediated)*(1.0 + t_max_mediated);

  water_energy = -epsilon*k_water*(sigma_gamma_direct*theta_direct+sigma_gamma_mediated*theta_mediated);

  return water_energy;
}

double FixBackbone::compute_burial_energy(int i_resno, int ires_type, double rho_i)
{
  double t[3][2];
  double burial_gamma_0, burial_gamma_1, burial_gamma_2, burial_energy;

  t[0][0] = tanh( burial_kappa*(rho_i - burial_ro_min[0]) );
  t[0][1] = tanh( burial_kappa*(burial_ro_max[0] - rho_i) );
  t[1][0] = tanh( burial_kappa*(rho_i - burial_ro_min[1]) );
  t[1][1] = tanh( burial_kappa*(burial_ro_max[1] - rho_i) );
  t[2][0] = tanh( burial_kappa*(rho_i - burial_ro_min[2]) );
  t[2][1] = tanh( burial_kappa*(burial_ro_max[2] - rho_i) );
  
  burial_gamma_0 = get_burial_gamma(i_resno, ires_type, 0);
  burial_gamma_1 = get_burial_gamma(i_resno, ires_type, 1);
  burial_gamma_2 = get_burial_gamma(i_resno, ires_type, 2);

  burial_energy = 0.0;
  burial_energy += -0.5*epsilon*k_burial*burial_gamma_0*(t[0][0] + t[0][1]);
  burial_energy += -0.5*epsilon*k_burial*burial_gamma_1*(t[1][0] + t[1][1]);
  burial_energy += -0.5*epsilon*k_burial*burial_gamma_2*(t[2][0] + t[2][1]);

  return burial_energy;
}

// generates a random but valid residue index
int FixBackbone::get_random_residue_index()
{
  int index;
  index = rand() % n; 
  return index;
}

// returns the CB-CB distance between two residues (CA for GLY)
double FixBackbone::get_residue_distance(int i_resno, int j_resno)
{
  double dx[3];
  double *xi, *xj;
  double r;
	
  if (se[i_resno]=='G') { xi = xca[i_resno]; }
  else { xi = xcb[i_resno]; }
  if (se[j_resno]=='G') { xj = xca[j_resno]; }
  else { xj = xcb[j_resno]; }
	
  dx[0] = xi[0] - xj[0];
  dx[1] = xi[1] - xj[1];
  dx[2] = xi[2] - xj[2];
  
  r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

  return r;
}

// returns the local density around residue i
double FixBackbone::get_residue_density(int i)
{
  return well->ro(i);
}

// returns the residue type of residue i_resno
int FixBackbone::get_residue_type(int i_resno)
{
  return se_map[se[i_resno]-'A'];
}

// given the native energy and the mean and variance of the decoy energy distribution, computes the frustration index
double FixBackbone::compute_frustration_index(double native_energy, double *decoy_stats)
{
  double frustration_index;

  frustration_index = (decoy_stats[0] - native_energy)/decoy_stats[1];

  return frustration_index;
}

void FixBackbone::compute_nmer_frust()
{
  int i, j, atomselect, nmer_contacts;
  double native_energy, frustration_index;

  atomselect = 0; // for the vmd script output

  // Double loop over all nmers
  for (i=0;i<n-nmer_frust_size;++i) {
    // get sequence of nmer starting at i
    get_nmer_seq(i, nmer_seq_i, 0);
    for (j=i+1;j<n-nmer_frust_size;++j) {
      // get sequence of nmer starting at j
      get_nmer_seq(j, nmer_seq_j, 0);
      // compute the number of contacts between the two nmers
      nmer_contacts = compute_nmer_contacts(i, j);
      // if the nmers have at least the number of required contacts and are non-overlapping, compute the frustration
      // this will need to be fixed to take multiple chains into account
      if (nmer_contacts >= nmer_contacts_cutoff && j-i >= nmer_frust_size) {
	native_energy = compute_nmer_native_ixn(i, j);
	compute_nmer_decoy_ixns(i, j);
	if (nmer_frust_trap_flag) {
	  // compute traps and keep track of how many times you wrote to tcl file
	  atomselect = compute_nmer_traps(i, j, atomselect, native_energy-nmer_frust_trap_num_sigma*nmer_decoy_ixn_stats[1], nmer_seq_i, nmer_seq_j); 
	  atomselect = compute_nmer_traps(j, i, atomselect, native_energy-nmer_frust_trap_num_sigma*nmer_decoy_ixn_stats[1], nmer_seq_j, nmer_seq_i);
	}
	frustration_index = compute_frustration_index(native_energy, nmer_decoy_ixn_stats);
	// write information out to output file
	fprintf(nmer_frust_output_file,"%d %d %d %s %s %f %f %f %f\n", i+1, j+1, nmer_contacts, nmer_seq_i, nmer_seq_j, native_energy, nmer_decoy_ixn_stats[0], nmer_decoy_ixn_stats[1], frustration_index);
	
	if(frustration_index > nmer_frust_min_frust_threshold || frustration_index < nmer_frust_high_frust_threshold || nmer_output_neutral_flag) {
	  // write information out to vmd script
	  fprintf(nmer_frust_vmd_script,"set sel%d [atomselect top \"resid %d and name CA\"]\n", i+nmer_frust_size/2, i+1+nmer_frust_size/2);
	  fprintf(nmer_frust_vmd_script,"set sel%d [atomselect top \"resid %d and name CA\"]\n", j+nmer_frust_size/2, j+1+nmer_frust_size/2);
	  fprintf(nmer_frust_vmd_script,"lassign [atomselect%d get {x y z}] pos1\n",atomselect);
	  atomselect += 1;
	  fprintf(nmer_frust_vmd_script,"lassign [atomselect%d get {x y z}] pos2\n",atomselect);
	  atomselect += 1;
	  if(frustration_index > nmer_frust_min_frust_threshold) {
	    fprintf(nmer_frust_vmd_script,"draw color green\n");
	  }
	  else if (frustration_index < nmer_frust_high_frust_threshold) {
	    fprintf(nmer_frust_vmd_script,"draw color red\n");
	  }
	  else {
	    fprintf(nmer_frust_vmd_script,"draw color blue\n");
	  }
	  fprintf(nmer_frust_vmd_script,"draw line $pos1 $pos2 style solid width 1\n");
	}
      }
    }
  }
  
  // after looping over all pairs, write out the end of the vmd script
  fprintf(nmer_frust_vmd_script, "mol modselect 0 top \"all\"\n");
  fprintf(nmer_frust_vmd_script, "mol modstyle 0 top newcartoon\n");
  fprintf(nmer_frust_vmd_script, "mol modcolor 0 top colorid 15\n");
}

void FixBackbone::compute_singlenmer_frust()
{
  double native_energy, frustration_index;
  int i, i_resno, atomselect;

  // write out style information to the vmd script
  fprintf(nmer_frust_vmd_script, "mol modselect 0 top \"all\"\n");
  fprintf(nmer_frust_vmd_script, "mol modstyle 0 top newcartoon\n");
  fprintf(nmer_frust_vmd_script, "mol modcolor 0 top colorid 15\n");

  // initialize representation counter
  atomselect = 0;

  // Loop over each nmer
  for (i=0;i<n-nmer_frust_size+1;++i) {
    // get the nmer sequence
    get_nmer_seq(i, nmer_seq_i, 0);
    // get the reside number
    i_resno = res_no[i]-1;
    // calculate the native energy
    native_energy = compute_singlenmer_native_ixn(i_resno);
    // calculate the decoy energies
    compute_singlenmer_decoy_ixns(i_resno);
    // Calculate frustration index
    frustration_index = compute_frustration_index(native_energy, nmer_decoy_ixn_stats);
    // write information out to output file
    fprintf(nmer_frust_output_file,"%d %s %f %f %f %f\n", i+1, nmer_seq_i, native_energy, nmer_decoy_ixn_stats[0], nmer_decoy_ixn_stats[1], frustration_index);
    
    if(frustration_index > nmer_frust_min_frust_threshold || frustration_index < nmer_frust_high_frust_threshold) {
      // write information out to vmd script
      atomselect += 1;
      fprintf(nmer_frust_vmd_script,"mol addrep 0\n");
      fprintf(nmer_frust_vmd_script,"mol modselect %d 0 resid %d to %d\n", atomselect, i_resno+1, i_resno+nmer_frust_size);
      fprintf(nmer_frust_vmd_script,"mol modstyle %d 0 VDW %f 12.000000\n", atomselect, 0.5*abs(frustration_index));
      // fprintf(nmer_frust_vmd_script,"mol modstyle %d 0 NewCartoon %f 10.000000 4.100000 0\n", atomselect, 0.5*abs(frustration_index));
      // fprintf(nmer_frust_vmd_script,"mol modstyle %d 0 QuickSurf 0.500000 0.500000 1.000000 1.000000\n", atomselect);
      fprintf(nmer_frust_vmd_script,"mol modmaterial %d 0 Transparent\n", atomselect);
      if(frustration_index > nmer_frust_min_frust_threshold) {
    	// color the residue green\n
        fprintf(nmer_frust_vmd_script,"mol modcolor %d 0 ColorID 7\n", atomselect);
      }
      else if(frustration_index < nmer_frust_high_frust_threshold) {
    	// color the residue red\n
	fprintf(nmer_frust_vmd_script,"mol modcolor %d 0 ColorID 1\n", atomselect);
      }
    }
  }  
}

double FixBackbone::compute_singlenmer_native_ixn(int i_resno)
{
  int j, j_resno, jres_type, j_chno;
  double native_energy, rho_j;

  // initialize native energy
  native_energy = 0.0;

  // Loop over each residue in the nmer
  for (j=i_resno;j<i_resno+nmer_frust_size;j++) {
    // get information about residue i
    j_resno = res_no[j]-1;
    jres_type = get_residue_type(j_resno);
    j_chno = chain_no[j]-1;
    rho_j = get_residue_density(j_resno);
    
    // Calculate native energy
    native_energy += compute_singleresidue_native_ixn(j_resno, jres_type, rho_j, j_chno, nmer_frust_cutoff, 1);
  }

  return native_energy;
}

void FixBackbone::compute_singlenmer_decoy_ixns(int i_resno)
{
  double rij, rho_i, rho_j;
  int jres_type, j_rand, decoy_i, j_resno, j_chno;
  int j;

  // do the decoy calculation nmer_frust_ndecoys times
  for (decoy_i=0; decoy_i<nmer_frust_ndecoys; decoy_i++) {
    // zero out this spot in the decoy energy array
    nmer_frust_decoy_energies[decoy_i] = 0.0;
    // get a random index to define the sequence of the decoy nmer
    j_rand = get_random_residue_index();
    // if the nmer would go off of the end of the sequence
    // then choose new j_rand and j_rand and repeat
    while(j_rand + nmer_frust_size > n) {
      j_rand = get_random_residue_index();
    }
    // Mutate the residues
    get_nmer_seq(j_rand, nmer_seq_j, 0);
    // Loop over each residue in the nmer
    for (j=i_resno;j<i_resno+nmer_frust_size;j++) {
      // get information about residue j
      j_resno = res_no[j]-1;
      // choose the next element in the nmer_seq_j array
      jres_type = se_map[nmer_seq_j[j-i_resno]-'A'];
      // choose a random residue type
      // jres_type = get_residue_type(get_random_residue_index());
      j_chno = chain_no[j]-1;
      rho_j = get_residue_density(j_resno);
      
      // compute decoy interaction energy, add to total
      nmer_frust_decoy_energies[decoy_i] += compute_singleresidue_native_ixn(j_resno, jres_type, rho_j, j_chno, nmer_frust_cutoff, 1);
    }
  }

  nmer_decoy_ixn_stats[0] = compute_array_mean(nmer_frust_decoy_energies, nmer_frust_ndecoys);
  nmer_decoy_ixn_stats[1] = compute_array_std(nmer_frust_decoy_energies, nmer_frust_ndecoys);
}

int FixBackbone::get_nmer_ss_dist(char *nmer_ss_j, char *nmer_ss_k)
{
  int i;
  int ss_dist;

  ss_dist=0;

  for(i=0; i<nmer_frust_size; i++) {
    if(nmer_ss_j[i]!=nmer_ss_k[i]) {
      ss_dist++;
    }
  }

  return ss_dist;
}

int FixBackbone::compute_nmer_traps(int i_start, int j_start, int atomselect, double threshold_energy, char *nmer_seq_1, char *nmer_seq_2)
{
  int k_start;
  double total_trap_energy;
  int tcl_index;
  int i, j, k;
  int ires_type, jres_type;
  double rho_i, rho_j, rij;
  int ss_dist;
  int backward;

  get_nmer_secondary_structure(i_start, nmer_ss_i);
  get_nmer_secondary_structure(j_start, nmer_ss_j);

  // start writing to tcl script at the appropriate index
  tcl_index = atomselect;
  static int rep_index=1;

  // loop over possible forward and backward sequences
  for (backward=0;backward<2;backward++) {
    // loop over all possible nmer traps
    for (k_start=0;k_start<n-nmer_frust_size;k_start++) {
      // if the nmer at position k_start doesn't overlap the nmers starting at position i_start, 
      // swap the sequence at k_start into j_start and calculate the energy
      if (abs(k_start-j_start)<=nmer_frust_size || abs(k_start-i_start)<=nmer_frust_size && i_start!=k_start) {
	continue;
      }
      get_nmer_seq(k_start, nmer_seq_k, backward);
      get_nmer_secondary_structure(k_start, nmer_ss_k);
      ss_dist = get_nmer_ss_dist(nmer_ss_j, nmer_ss_k);
      // check to make sure that the secondary structure distance constraint is satisfied
      if (ss_dist > nmer_frust_size*(1-nmer_frust_ss_frac)) {
	continue;
      }
      total_trap_energy = 0.0;
    
      // loop over all residues individually, compute burial energies
      for (i = i_start; i < i_start+nmer_frust_size; i++) {
	ires_type = get_residue_type(i);
	rho_i = get_residue_density(i);
	total_trap_energy += compute_burial_energy(i, ires_type, rho_i);
      }
    
      for (j = j_start; j < j_start+nmer_frust_size; j++) {
	// get the sequence starting from k rather than j
	jres_type = get_residue_type(((1-backward)*(j-j_start))+backward*(nmer_frust_size-(j-j_start)));
	rho_j = get_residue_density(j);
	total_trap_energy += compute_burial_energy(j, jres_type, rho_j);
      }
    
      // loop over all pairs of residues between the two nmers, compute water interaction
      for (i = i_start; i < i_start+nmer_frust_size; i++) {
	// get information about residue i
	ires_type = get_residue_type(i);
      
	for (j = j_start; j < j_start+nmer_frust_size; j++) {
	  // get the sequence starting from k rather than j
	  jres_type = get_residue_type(((1-backward)*(j-j_start))+backward*(nmer_frust_size-(j-j_start)));
	
	  // get interaction parameters
	  rij = get_residue_distance(i, j);
	  rho_i = get_residue_density(i);
	  rho_j = get_residue_density(j);
	
	  // compute water interaction energy, add to total
	  total_trap_energy += compute_water_energy(rij, i, j, ires_type, jres_type, rho_i, rho_j);
	}
      }
      if(total_trap_energy<threshold_energy) {
	// write out to trap file and trap tcl script
	if (backward) {
	  fprintf(nmer_frust_trap_file,"%d %s %s %d %s %s %f %d %s <-- %s %f \n", i_start+1, nmer_seq_1, nmer_ss_i, j_start+1, nmer_seq_2, nmer_ss_j, threshold_energy, k_start+1, nmer_seq_k, nmer_ss_k, total_trap_energy);
	}
	else {
	  fprintf(nmer_frust_trap_file,"%d %s %s %d %s %s %f %d %s --> %s %f \n", i_start+1, nmer_seq_1, nmer_ss_i, j_start+1, nmer_seq_2, nmer_ss_j, threshold_energy, k_start+1, nmer_seq_k, nmer_ss_k, total_trap_energy);
	}
	if(nmer_frust_draw_trap_flag) {
	  // if this is a case of "self-recognition", make that part of the sequence purple
	  if (i_start == k_start) {
	    fprintf(nmer_frust_vmd_script,"mol addrep 0\n",rep_index);
	    fprintf(nmer_frust_vmd_script,"mol modselect %d 0 resid %d to %d\n",rep_index,i_start+1,i_start+nmer_frust_size);
	    fprintf(nmer_frust_vmd_script,"mol modcolor %d 0 ColorID 11\n",rep_index);
	    fprintf(nmer_frust_vmd_script,"mol modstyle %d 0 NewCartoon 0.350000 10.000000 4.100000 0\n",rep_index);
	    rep_index++;
	  }
	  // if not self-recognition, draw a purple line
	  else {
	    fprintf(nmer_frust_vmd_script,"set sel%d [atomselect top \"resid %d and name CA\"]\n", i_start+nmer_frust_size/2, i_start+1+nmer_frust_size/2);
	    fprintf(nmer_frust_vmd_script,"set sel%d [atomselect top \"resid %d and name CA\"]\n", k_start+nmer_frust_size/2, k_start+1+nmer_frust_size/2);
	    fprintf(nmer_frust_vmd_script,"lassign [atomselect%d get {x y z}] pos1\n",tcl_index);
	    tcl_index += 1;
	    fprintf(nmer_frust_vmd_script,"lassign [atomselect%d get {x y z}] pos2\n",tcl_index);
	    tcl_index += 1;
	    fprintf(nmer_frust_vmd_script,"draw color purple\n");
	    if(backward) {
	      fprintf(nmer_frust_vmd_script,"draw line $pos1 $pos2 style dashed width 1\n");  
	    }
	    else {
	      fprintf(nmer_frust_vmd_script,"draw line $pos1 $pos2 style solid width 1\n");  
	    }
	  }
	}
      }
    }
  }
  return tcl_index; // return the new tcl_index (atomselect) so that you can continue writing to the tcl file in the other functions
}

// computes the number of contacts between two nmers of size nmer_frust_size starting at positions i_resno and j_resno
int FixBackbone::compute_nmer_contacts(int i_start, int j_start)
{
  int numcontacts;
  int i,j;
  double distance;

  numcontacts = 0;

  // loop over all pairs of residues between the two nmers
  for (i = i_start; i < i_start+nmer_frust_size; i++) {
    for (j = j_start; j < j_start+nmer_frust_size; j++) {
      // if you are not beyond minimum sequence separation for the water potential, then you are not a contact
      if(abs(i-j) < contact_cutoff) continue;
      // compute distance between the two residues
      distance = get_residue_distance(i, j);
      // if less than the threshold, incrememnt number of contacts
      if (distance < nmer_frust_cutoff) {
	numcontacts++;
      }
    }
  }

  return numcontacts;
}

// returns a string that is the sequence of the nmer of size nmer_frust_size starting at position i
void FixBackbone::get_nmer_seq(int i_start, char *nmer_seq, int backward)
{
  int i;

  for(i=0; i<nmer_frust_size; i++) {
    nmer_seq[i] = se[((1-backward)*(i_start+i))+backward*(i_start+nmer_frust_size-i)];
  }
  //printf("%s\n",nmer_seq);
}

void FixBackbone::get_nmer_secondary_structure(int i_start, char *nmer_secondary_structure)
{
  int i;

  for(i=0; i<nmer_frust_size; i++) {
    // if assigned to be beta, show an E
    if(aps[4][i+i_start]==1.0) nmer_secondary_structure[i] = 'E';
    // if assigned to be helix, show an H
    if(aps[3][i+i_start]==1.0) nmer_secondary_structure[i] = 'H';
    // if assigned to both, show an !
    if(aps[4][i+i_start]==1.0 && aps[3][i+i_start]==1.0) nmer_secondary_structure[i] = '!';
    // if anything else, show a -
    if(aps[4][i+i_start]!=1.0 && aps[3][i+i_start]!=1.0) nmer_secondary_structure[i] = '-';
  }
}

double FixBackbone::compute_nmer_native_ixn(int i_start, int j_start)
{
  double total_native_energy, rij, rho_i, rho_j;
  int ires_type, jres_type;
  int i, j;

  total_native_energy = 0.0;

  // loop over all residues individually, compute burial energies
  for (i = i_start; i < i_start+nmer_frust_size; i++) {
    ires_type = get_residue_type(i);
    rho_i = get_residue_density(i);
    total_native_energy += compute_burial_energy(i, ires_type, rho_i);
  }

  for (j = j_start; j < j_start+nmer_frust_size; j++) {
    jres_type = get_residue_type(j);
    rho_j = get_residue_density(j);
    total_native_energy += compute_burial_energy(j, jres_type, rho_j);
  }

  // loop over all pairs of residues between the two nmers, compute water interaction
  for (i = i_start; i < i_start+nmer_frust_size; i++) {
    // get information about residue i
    ires_type = get_residue_type(i);
    
    for (j = j_start; j < j_start+nmer_frust_size; j++) {
      // get information about residue j
      jres_type = get_residue_type(j);

      // get interaction parameters
      rij = get_residue_distance(i, j);
      rho_i = get_residue_density(i);
      rho_j = get_residue_density(j);

      // compute water interaction energy, add to total
      total_native_energy += compute_water_energy(rij, i, j, ires_type, jres_type, rho_i, rho_j);
    }
  }

  return total_native_energy;
}

void FixBackbone::compute_nmer_decoy_ixns(int i_start, int j_start)
{
  double total_decoy_energy, rij, rho_i, rho_j, decoy_energy;
  int ires_type, jres_type, i_rand, j_rand, decoy_i;
  int i, j, itemp, jtemp;
  int atomselect;

  atomselect = 0;

  // do the decoy calculation nmer_frust_ndecoys times
  for (decoy_i=0; decoy_i<nmer_frust_ndecoys; decoy_i++) {
    // zero out this spot in the decoy energy array
    nmer_frust_decoy_energies[decoy_i] = 0.0;
    // get two random indices to define the sequence of the decoy nmers
    i_rand = get_random_residue_index();
    j_rand = get_random_residue_index();
    // make sure that j > i so that we can test for overlap
    if(i_rand > j_rand) {
      itemp = i_rand;
      jtemp = j_rand;
      j_rand = itemp;
      i_rand = jtemp;
    }
    // if either nmer would go off of the end of the sequence or they are overlapping
    // then choose new i_rand and j_rand and repeat
    while(i_rand + nmer_frust_size > n || j_rand + nmer_frust_size > n || j_rand-i_rand < nmer_frust_size) {
      i_rand = get_random_residue_index();
      j_rand = get_random_residue_index();
      if(i_rand > j_rand) {
	itemp = i_rand;
	jtemp = j_rand;
	j_rand = itemp;
	i_rand = jtemp;
      }
    }
    // loop over all residues individually, compute burial energies
    for (i = i_start; i < i_start+nmer_frust_size; i++) {
      ires_type = get_residue_type(i_rand+i-i_start);
      rho_i = get_residue_density(i);
      nmer_frust_decoy_energies[decoy_i] += compute_burial_energy(i, ires_type, rho_i);
    }
    for (j = j_start; j < j_start+nmer_frust_size; j++) {
      jres_type = get_residue_type(j_rand+j-j_start);
      rho_j = get_residue_density(j);
      nmer_frust_decoy_energies[decoy_i] += compute_burial_energy(j, jres_type, rho_j);
    }

    // loop over all pairs of residues between the two nmers
    for (i = i_start; i < i_start+nmer_frust_size; i++) {
      // assign random residue type to residue i
      ires_type = get_residue_type(i_rand+i-i_start);
      
      for (j = j_start; j < j_start+nmer_frust_size; j++) {
	// assign random residue type to residue j
	jres_type = get_residue_type(j_rand+j-j_start);
	
	// get interaction parameters
	rij = get_residue_distance(i, j);
	rho_i = get_residue_density(i);
	rho_j = get_residue_density(j);

	// compute decoy interaction energy, add to total
	nmer_frust_decoy_energies[decoy_i] += compute_water_energy(rij, i, j, ires_type, jres_type, rho_i, rho_j);
      }
    }
  }

  nmer_decoy_ixn_stats[0] = compute_array_mean(nmer_frust_decoy_energies, nmer_frust_ndecoys);
  nmer_decoy_ixn_stats[1] = compute_array_std(nmer_frust_decoy_energies, nmer_frust_ndecoys);
}

void FixBackbone::compute_fragment_memory_table()
{
  int i, j, js, je, ir, i_fm, k, itb, iatom[4], jatom[4], iatom_type[4], jatom_type[4];
  int i_first_res, i_last_res, i_resno, j_resno, ires_type, jres_type;
  double r, rf, dr, drsq, V, force;
  double fm_sigma_sq, frag_mem_gamma, epsilon_k_weight, epsilon_k_weight_gamma;
  Fragment_Memory *frag;
  
  iatom_type[0] = Fragment_Memory::FM_CA;
  iatom_type[1] = Fragment_Memory::FM_CA;
  iatom_type[2] = Fragment_Memory::FM_CB;
  iatom_type[3] = Fragment_Memory::FM_CB;
  
  jatom_type[0] = Fragment_Memory::FM_CA; 
  jatom_type[1] = Fragment_Memory::FM_CB; 
  jatom_type[2] = Fragment_Memory::FM_CA; 
  jatom_type[3] = Fragment_Memory::FM_CB;
  
  for (i=0; i<n; ++i) {  
    iatom[0] = alpha_carbons[i];
    iatom[1] = alpha_carbons[i];
    iatom[2] = beta_atoms[i];
    iatom[3] = beta_atoms[i];
	  
    //	  i_resno = res_no[i]-1;
    i_resno = i;
    ires_type = se_map[se[i_resno]-'A'];
	  
    for (i_fm=0; i_fm<ilen_fm_map[i_resno]; ++i_fm) {
      frag = frag_mems[ frag_mem_map[i_resno][i_fm] ];
		
      epsilon_k_weight = epsilon*k_frag_mem*frag->weight;
		
      js = i+fm_gamma->minSep();
      je = MIN(frag->pos+frag->len-1, i+fm_gamma->maxSep());
      //		if (je>=n || res_no[je]-res_no[i]!=je-i) error->all(FLERR,"Missing residues in memory potential");
      if (je>=n) error->all(FLERR,"Missing residues in memory potential");
		
      for (j=js;j<=je;++j) {
	//		  j_resno = res_no[j]-1;
	j_resno = j;
	jres_type = se_map[se[j_resno]-'A'];
		  
	//		  if (chain_no[i]!=chain_no[j]) error->all(FLERR,"Fragment Memory: Interaction between residues of different chains");
		  
	fm_sigma_sq = pow(abs(i_resno-j_resno), 2*fm_sigma_exp);
	fm_sigma_sq = fm_sigma_sq*frag_table_well_width*frag_table_well_width;
			  
	if (!fm_gamma->fourResTypes()) {
	  frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, i_resno, j_resno);
	} else {
	  frag_mem_gamma = fm_gamma->getGamma(ires_type, jres_type, frag->resType(i_resno), frag->resType(j_resno), i_resno, j_resno);
	}
	if (fm_gamma->error==fm_gamma->ERR_CALL) error->all(FLERR,"Fragment_Memory: Wrong call of getGamma() function");
		  
	epsilon_k_weight_gamma = epsilon_k_weight*frag_mem_gamma;
		  
	jatom[0] = alpha_carbons[j];
	jatom[1] = beta_atoms[j];
	jatom[2] = alpha_carbons[j];
	jatom[3] = beta_atoms[j];
		  
	for (k=0;k<4;++k) {
	  if ( iatom_type[k]==frag->FM_CB && (se[i_resno]=='G' || frag->getSe(i_resno)=='G') ) continue;
	  if ( jatom_type[k]==frag->FM_CB && (se[j_resno]=='G' || frag->getSe(j_resno)=='G') ) continue;
			
	  itb = 4*tb_nbrs*i + 4*(j-js) + k;
	  if (!fm_table[itb])
	    fm_table[itb] = new TBV[tb_size];
			
	  rf = frag->Rf(i_resno, iatom_type[k], j_resno, jatom_type[k]);
	  if (frag->error==frag->ERR_CALL) error->all(FLERR,"Fragment_Memory: Wrong call of Rf() function");
	  for (ir=0;ir<tb_size;++ir) {
	    r = tb_rmin + ir*tb_dr;
				
	    dr = r - rf;
	    drsq = dr*dr;
				
	    V = -epsilon_k_weight_gamma*exp(-drsq/(2*fm_sigma_sq));
				
	    fm_table[itb][ir].energy += V;
			
	    fm_table[itb][ir].force += V*dr/(fm_sigma_sq*r);
				
	    //				fm_table[i][j-js][k][ir].energy = -epsilon_k_weight_gamma*exp(-drsq/(2*fm_sigma_sq));
				
	    //				fm_table[i][j-js][k][ir].force = V*dr/(fm_sigma_sq*r);
	  }
	}
      }
    }
  }
  if(fm_energy_debug_flag) {
    output_fragment_memory_table();
  }
}

// this routine outputs the contents of fm_table (only the energies) to a file called fmenergies.log
void FixBackbone::output_fragment_memory_table()
{
  int itb,ir; // loop variables
  double energyvalue;

  // loop over all interactions
  for (itb=0; itb<4*n*tb_nbrs; itb++) {
    // loop over all distances
    for (ir=0; ir<tb_size; ir++) {
      // don't try to read out the energies if it was never allocated because of the exception for glycines
      // instead, output a row of zeroes so that the row indices remain meaningful
      if (fm_table[itb] == NULL)
	{
	  energyvalue = 0.0;
	}
      else
	{
	  energyvalue = fm_table[itb][ir].energy;
	}
      // output fm_table energies to file for debugging/visualization
      fprintf(fmenergiesfile,"%.4f ",energyvalue);
    }
    // put every interaction on a new line of the file
    fprintf(fmenergiesfile,"\n");
  }
}

void FixBackbone::compute_membrane_potential(int i)
{

//  k_bin is coming from the input
//  gamma[0][0] is an array coming from input
//  k_overall_memb coming from input
//  rho0_distor is coming from the input 
//  rho0_max= is coming from the input
//  memb_len= is coming from the input
//  memb_pore_type = is coming from the input


  double V, x_actual, y_actual, z_actual;
  double dx, dy, dz, *xi;
  double memb_a, memb_b;
  double rho_actual, rho0;
  double s_per, s_mem, s_cyt, s_por, s_nopor;
  double dz_per, dz_mem, dz_cyt, dr1_dz;
  double dz_s_por, dx_s_por, dy_s_por;
  double dz_s_nopor, dx_s_nopor, dy_s_nopor;
  double dz_s_por_smem, dz_s_nopor_smem;
  double dV_dx, dV_dy, dV_dz;
  double memb_force_x, memb_force_y, memb_force_z;

  int iatom;

//  int z_res[i] = is coming from zim file

  int i_resno = res_no[i]-1;

  if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
  else { xi = xcb[i]; iatom  = beta_atoms[i]; }


  x_actual=xi[0];
  y_actual=xi[1];
  z_actual=xi[2];

  dx=x_actual-memb_xo[0];
  dy=y_actual-memb_xo[1];
  dz=z_actual-memb_xo[2];

  memb_a = rho0_distor*rho0_max;
  memb_b = memb_len/2;

  rho_actual = sqrt(dx*dx+dy*dy);

// if (memb_pore_type == 0){
  rho0=(rho0_max-memb_a) + ((memb_a)/(memb_len))*(dz+memb_b);
// }
// else if (memb_pore_type == 1) {
//  rho0=rho0_max-(dz<=memb_b ? sqrt(1-pow(dz/memb_b,2)):0)*memb_a;
// }



//definition of swithing functions
  s_per=0.5*(1+tanh(k_bin*(dz-memb_b)));
  s_mem=0.5*((tanh(k_bin*(dz+memb_b)))+(tanh(k_bin*(memb_b-dz))));
  s_cyt=0.5*(1+tanh(k_bin*(-memb_b-dz)));
  s_por=0.5*(1-(tanh(k_bin*(rho_actual-rho0))));
  s_nopor=(1-s_por);

  if (z_res[i] == 1) {
    V=(-g_memb[0][0]*s_per)+(-g_memb[0][1]*s_cyt)+(g_memb[0][2]*s_mem*s_nopor)+(-g_memb[0][3]*s_mem*s_por);
  }
  else if (z_res[i] == 2){
    V=(g_memb[1][0]*s_per)+(g_memb[1][1]*s_cyt)+(-g_memb[1][2]*s_mem*s_nopor)+(g_memb[1][3]*s_mem*s_por);
  }
  else if (z_res[i] == 3) {
    V=(-g_memb[2][0]*s_per)+(-g_memb[2][1]*s_cyt)+(g_memb[2][2]*s_mem*s_nopor)+(-g_memb[2][3]*s_mem*s_por);
  }

// modify energy matrix to accomodate membrane potential
  energy[ET_MEMB] += epsilon*k_overall_memb*V;

// parcial derivatives
dz_per=0.5*k_bin*(1-pow((tanh(k_bin*(dz-memb_b))),2));
dz_mem=-0.5*k_bin*pow(tanh(k_bin*(dz+memb_b)),2)+0.5*k_bin*pow(tanh(k_bin*(memb_b-dz)),2);
dz_cyt=-0.5*k_bin*(1-pow((tanh(k_bin*(-memb_b-dz))),2));

// if (memb_pore_type == 0){
  dr1_dz=((memb_a)/(memb_len));
//  }
// else if (memb_pore_type == 1){
//  dr1_dz=(1/(dz<=memb_b ? (sqrt(1-pow(dz/memb_b,2))):0))*(dz/pow(memb_b,2))*memb_a;
//  }

dz_s_por=0.5*k_bin*(1-pow(tanh(k_bin*(rho_actual-rho0)),2))*dr1_dz;

dx_s_por=((-0.5*k_bin*dx)/rho_actual)*(1-pow(tanh(k_bin*(rho_actual-rho0)),2));
dy_s_por=((-0.5*k_bin*dy)/rho_actual)*(1-pow(tanh(k_bin*(rho_actual-rho0)),2));

//some definitions to be used in "calculate general derivatives"
dx_s_nopor=-dx_s_por;
dy_s_nopor=-dy_s_por;
dz_s_nopor=-dz_s_por;

dz_s_por_smem=s_mem*dz_s_por+dz_mem*s_por;
dz_s_nopor_smem=s_mem*dz_s_nopor+dz_mem*s_nopor;

//calculate general derivatives
 if (z_res[i] == 1) {
   dV_dx=g_memb[0][2]*s_mem*dx_s_nopor+(-g_memb[0][3])*s_mem*dx_s_por;
   dV_dy=g_memb[0][2]*s_mem*dy_s_nopor+(-g_memb[0][3])*s_mem*dy_s_por;
   dV_dz=-g_memb[0][0]*dz_per+(-g_memb[0][1])*dz_cyt+g_memb[0][2]*dz_s_nopor_smem+(-g_memb[0][3])*dz_s_por_smem;
   }
 else if (z_res[i] == 2){
   dV_dx=-g_memb[1][2]*s_mem*dx_s_nopor+g_memb[1][3]*s_mem*dx_s_por;
   dV_dy=-g_memb[1][2]*s_mem*dy_s_nopor+g_memb[1][3]*s_mem*dy_s_por;
   dV_dz=g_memb[1][0]*dz_per+g_memb[1][1]*dz_cyt+(-g_memb[1][2])*dz_s_nopor_smem+g_memb[1][3]*dz_s_por_smem;
   }
 else if (z_res[i] == 3){
   dV_dx=g_memb[2][2]*s_mem*dx_s_nopor+(-g_memb[2][3])*s_mem*dx_s_por;
   dV_dy=g_memb[2][2]*s_mem*dy_s_nopor+(-g_memb[2][3])*s_mem*dy_s_por;
   dV_dz=-g_memb[2][0]*dz_per+(-g_memb[2][1])*dz_cyt+g_memb[2][2]*dz_s_nopor_smem+(-g_memb[2][3])*dz_s_por_smem;
   }

// calculate forces
  memb_force_x = -epsilon*k_overall_memb*dV_dx;
  memb_force_y = -epsilon*k_overall_memb*dV_dy;
  memb_force_z = -epsilon*k_overall_memb*dV_dz;

// add forces
  f[iatom][0] += memb_force_x;
  f[iatom][1] += memb_force_y;
  f[iatom][2] += memb_force_z;

}

void FixBackbone::compute_solvent_barrier(int i, int j)
{
  if (chain_no[i]==chain_no[j] && res_no[j]-res_no[i]<ssb_ij_sep) return;

  double dx[3], force1, force2;
  double *xi, *xj, r, rmin1, rmax1, rmin2, rmax2, rshift;
  double t_min1, t_max1, theta1, t_min2, t_max2, theta2;
  int iatom, jatom;
  
  int i_resno = res_no[i]-1;
  int j_resno = res_no[j]-1;
  
  int ires_type = se_map[se[i_resno]-'A'];
  int jres_type = se_map[se[j_resno]-'A'];
  
  if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
  else { xi = xcb[i]; iatom  = beta_atoms[i]; }
  if (se[j_resno]=='G') { xj = xca[j]; jatom = alpha_carbons[j]; }
  else { xj = xcb[j]; jatom  = beta_atoms[j]; }
  
  dx[0] = xi[0] - xj[0];
  dx[1] = xi[1] - xj[1];
  dx[2] = xi[2] - xj[2];

  r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

  rmin1 = ssb_rmin1;
  rmax1 = ssb_rmax1;
  rmin2 = ssb_rmin2;
  rmax2 = ssb_rmax2;
  if(ssb_rad_cor){
    rshift = ssb_rshift[ires_type]+ssb_rshift[jres_type];
    rmin1 += rshift;
    rmax1 += rshift;    
    rmin2 += rshift;
    rmax2 += rshift;
  }

  // apply a distance cutoff criterion, cutoff = rmax + 10/kappa
  if(r>rmax1+10/ssb_kappa && r>rmax2+10/ssb_kappa) return;

  t_min1 = tanh(ssb_kappa*(r-rmin1));
  t_max1 = tanh(ssb_kappa*(rmax1-r));
  t_min2 = tanh(ssb_kappa*(r-rmin2));
  t_max2 = tanh(ssb_kappa*(rmax2-r));

  theta1 = 0.5*(t_min1+t_max1);
  theta2 = 0.5*(t_min2+t_max2);

  energy[ET_SSB] += epsilon*k_solventb1*theta1;
  energy[ET_SSB] += epsilon*k_solventb2*theta2;

  force1 = -epsilon*k_solventb1*ssb_kappa*theta1*(t_max1-t_min1)/r;
  force2 = -epsilon*k_solventb2*ssb_kappa*theta2*(t_max2-t_min2)/r;

  f[iatom][0] += force1*dx[0];
  f[iatom][1] += force1*dx[1];
  f[iatom][2] += force1*dx[2];
  f[iatom][0] += force2*dx[0];
  f[iatom][1] += force2*dx[1];
  f[iatom][2] += force2*dx[2];

  f[jatom][0] += -force1*dx[0];
  f[jatom][1] += -force1*dx[1];
  f[jatom][2] += -force1*dx[2];
  f[jatom][0] += -force2*dx[0];
  f[jatom][1] += -force2*dx[1];
  f[jatom][2] += -force2*dx[2];
}

void FixBackbone::compute_DebyeHuckel_Interaction(int i, int j)
{
  if (abs(i-j)<debye_huckel_min_sep) return;
  
  double dx[3];
  double *xi, *xj, r;
  int iatom, jatom;
  double charge_i = 0.0;
  double charge_j = 0.0;
  double term_qq_by_r = 0.0;
  double force_term = 0.0;
    
  charge_i = charge_on_residue[i]; 
  charge_j = charge_on_residue[j]; 
  
  if (charge_i == 0 && charge_j == 0) return;
  
  int i_resno = res_no[i]-1;
  int j_resno = res_no[j]-1;
  
  int ires_type = se_map[se[i_resno]-'A'];
  int jres_type = se_map[se[j_resno]-'A'];
  
  if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
  else { xi = xcb[i]; iatom  = beta_atoms[i]; }
  if (se[j_resno]=='G') { xj = xca[j]; jatom = alpha_carbons[j]; }
  else { xj = xcb[j]; jatom  = beta_atoms[j]; }
  
  dx[0] = xi[0] - xj[0];
  dx[1] = xi[1] - xj[1];
  dx[2] = xi[2] - xj[2];
  
  r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  
  if( (charge_i > 0.0) && (charge_j > 0.0) ) {
    term_qq_by_r = k_PlusPlus*charge_i*charge_j/r;
  }
  else if(charge_i < 0.0 && charge_j < 0.0) {
    term_qq_by_r = k_MinusMinus*charge_i*charge_j/r;
  }
  else if( (charge_i < 0.0 && charge_j > 0.0) || (charge_i > 0.0 && charge_j < 0.0)) {
    term_qq_by_r = k_PlusMinus*charge_i*charge_j/r;
  }
  
  double term_energy = epsilon*term_qq_by_r*exp(-k_screening*r/screening_length);
  energy[ET_DH] += term_energy;
  
  force_term = (term_energy/r)*(1.0/r + k_screening/screening_length);
  
  f[iatom][0] += force_term*dx[0];
  f[iatom][1] += force_term*dx[1];
  f[iatom][2] += force_term*dx[2];
  
  f[jatom][0] += -force_term*dx[0];
  f[jatom][1] += -force_term*dx[1];
  f[jatom][2] += -force_term*dx[2];
}

void FixBackbone::compute_debyehuckel_optimization()
{
  // computes and writes out the energies for debyehuckel

  double *xi, *xj, dx[3];
  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rij;
  double debyehuckel_energy;
  double debyehuckel_energies[2][2];
  double contact_norm[2][2];
  int i_charge_type, j_charge_type;
  double charge_i, charge_j;

  debyehuckel_energy=0.0;

  // array initialization
  for (i=0;i<2;++i) {
    for (j=i;j<2;++j) {
      debyehuckel_energies[i][j] = debyehuckel_energies[j][i] = 0.0;
      contact_norm[i][j] = contact_norm[j][i] = 0.0;
    }
  }
  
  // Double loop over all residue pairs
  for (i=0;i<n;++i) {
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = se_map[se[i_resno]-'A'];
    i_chno = chain_no[i]-1;

    // check if ires_type is D, E, R or K; if not, skip; if so, assign charge type
    if (se[i_resno]=='R' || se[i_resno]=='K') {
      i_charge_type = 0;
      charge_i = 1.0;
    }
    else if (se[i_resno]=='D' || se[i_resno]=='E') {
      i_charge_type = 1;
      charge_i = -1.0;
    }
    else {
      continue;
    }
    
    for (j=i+1;j<n;++j) {
      // get information about residue j
      j_resno = res_no[j]-1;
      jres_type = se_map[se[j_resno]-'A'];
      j_chno = chain_no[j]-1;

      // check if ires_type is D, E, R or K; if not, skip; if so, assign charge type
      if (se[j_resno]=='R' || se[j_resno]=='K') {
	j_charge_type = 0;
	charge_j = 1.0;
      }
      else if (se[j_resno]=='D' || se[j_resno]=='E') {
	j_charge_type = 1;
	charge_j = -1.0;
      }
      else {
	continue;
      }
   
      // Select beta atom unless the residue type is GLY, then select alpha carbon
      if (se[i_resno]=='G') { xi = xca[i]; }
      else { xi = xcb[i]; }
      if (se[j_resno]=='G') { xj = xca[j]; }
      else { xj = xcb[j]; }
	  
      // compute distance between the two atoms
      dx[0] = xi[0] - xj[0];
      dx[1] = xi[1] - xj[1];
      dx[2] = xi[2] - xj[2];
      rij = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      
      // if the atoms are within the threshold, compute the energies
      if (abs(i-j)>=debye_huckel_min_sep || i_chno != j_chno) {
	// calculate debyehuckel energies
	double term_qq_by_r = 0.0;

	if( (charge_i > 0.0) && (charge_j > 0.0) ) {
	  term_qq_by_r = k_PlusPlus*charge_i*charge_j/rij;
	}
	else if(charge_i < 0.0 && charge_j < 0.0) {
	  term_qq_by_r = k_MinusMinus*charge_i*charge_j/rij;
	}
	else if( (charge_i < 0.0 && charge_j > 0.0) || (charge_i > 0.0 && charge_j < 0.0)) {
	  term_qq_by_r = k_PlusMinus*charge_i*charge_j/rij;
	}
	
	double term_energy = epsilon*term_qq_by_r*exp(-k_screening*rij/screening_length);
	debyehuckel_energies[i_charge_type][j_charge_type] += term_energy;
	contact_norm[i_charge_type][j_charge_type] += 1.0;
      }
    }
  }

  for (i=0;i<2;++i) {
    for (j=i;j<2;++j) {
      debyehuckel_energies[i][j] = debyehuckel_energies[i][j] + debyehuckel_energies[j][i];
      contact_norm[i][j] = contact_norm[i][j] + contact_norm[j][i];
    }
  }

  for (i=0;i<2;++i) {
    debyehuckel_energies[i][i] /= 2.0;
    contact_norm[i][i] /= 2.0;
  }

  // write output calculated using native sequence on step 0
  // if step !=0 then write output calculated with shuffled sequence
  if (ntimestep == 0){
    fprintf(debyehuckel_native_optimization_file,"%f %f %f \n", debyehuckel_energies[0][0], debyehuckel_energies[1][1], debyehuckel_energies[1][0]);
    fprintf(debyehuckel_native_optimization_norm_file,"%f %f %f \n", contact_norm[0][0], contact_norm[1][1], contact_norm[1][0]);
  }
  else {
    fprintf(debyehuckel_optimization_file,"%f %f %f \n", debyehuckel_energies[0][0], debyehuckel_energies[1][1], debyehuckel_energies[1][0]);
    fprintf(debyehuckel_optimization_norm_file,"%f %f %f \n", contact_norm[0][0], contact_norm[1][1], contact_norm[1][0]);
  }
}

void FixBackbone::read_amylometer_sequences(char *amylometer_sequence_file, int amylometer_nmer_size, int amylometer_mode)
{
  // Read in sequences, split into n-mers
  FILE *file;  
  char ln[10000], *line;
  size_t number_of_aminoacids;
  number_of_nmers = 0;
  
  file = fopen(amylometer_sequence_file,"r");

  if (!file) error->all(FLERR,"Amylometer: Error opening amylometer sequences file");
  while ( fgets ( ln, sizeof ln, file ) != NULL ) {
    line = trim(ln);
    number_of_aminoacids = strlen(line);
    if (line[0]=='#') continue;
    for (int i=0; i<number_of_aminoacids-amylometer_nmer_size+1; i++) {
      for (int j=0; j<amylometer_nmer_size; j++) {
	//fprintf(output_file, "%c", line[i+j]);
      }
      number_of_nmers++;
      //fprintf(output_file, "\n");
    }
  }
  fclose(file);

  // allocate nmer_array
  nmer_array = (int**)malloc(number_of_nmers * sizeof(int*));
  for (int i = 0; i < number_of_nmers; i++) {
    nmer_array[i] = (int*)malloc(amylometer_nmer_size * sizeof(int));
  }

  file = fopen(amylometer_sequence_file,"r");
  if (!file) error->all(FLERR,"Amylometer: Error opening amylometer sequences file\n");
  int nmer_index = 0;
  while ( fgets ( ln, sizeof ln, file ) != NULL ) {
    line = trim(ln);
    number_of_aminoacids = strlen(line);
    if (line[0]=='#') continue;
    for (int i=0; i<number_of_aminoacids-amylometer_nmer_size+1; i++) {
      for (int j=0; j<amylometer_nmer_size; j++) {
	nmer_array[nmer_index][j] = line[i+j];
      }
      nmer_index++;
    }
  }
  fclose(file);

  // for (int i=0; i < number_of_nmers; i++) {
  //   for (int j=0; j < amylometer_nmer_size; j++) {
  //     printf("%d",se_map[nmer_array[i][j]-'A']);
  //   }
  //   printf("\n");
  // }

  return;
}

void FixBackbone::compute_amylometer()
{
  if (comm->me==0) print_log("Running amylometer...\n");
  FILE *amylometer_energy_file;
  amylometer_energy_file = fopen("amylometer_energy.log", "w");
  char eheader[] = "\tChain   \tShake   \tChi     \tRama    \tExcluded\tDSSP    \tP_AP    \tWater   \tBurial  \tHelix   \tAMH-Go  \tFrag_Mem\tVec_FM  \tSSB     \tVTotal\n";
  fprintf(amylometer_energy_file, "%s", eheader);
  
  FILE *nmer_output_file;
  nmer_output_file = fopen("nmer_output","w");

  // homogeneous mode
  if (amylometer_mode == 1) {
    fprintf(nmer_output_file, "nmer\n");
    // loop over all nmers
    for (int i=0; i<number_of_nmers; i++) { 
      // mutate sequence
      for (int j=0; j<n; j++) {
	se[j] = nmer_array[i][j % amylometer_nmer_size];
	if (j < amylometer_nmer_size) {
	  fprintf(nmer_output_file, "%c",se[j]);
	}
      }
      fprintf(nmer_output_file, "\n");
      // calculate energy, output
      for (int k=0;k<nEnergyTerms;++k) energy_all[k] = 0.0; // clear energy values
      compute_backbone();                                   // compute energies
      // output energies to file
      for (int k=1;k<nEnergyTerms;++k) fprintf(amylometer_energy_file, "\t%8.6f", energy_all[k]);
      fprintf(amylometer_energy_file, "\t%8.6f\n", energy_all[ET_TOTAL]);
    }
  }
  // heterogeneous mode
  else if (amylometer_mode == 2) {
    Fragment_Memory *native_structure = new Fragment_Memory(0, 0, amylometer_nmer_size+number_of_nmers-1, 0.0, amylometer_structure_file, 0);
    int iatom_type = Fragment_Memory::FM_CA;
    int jatom_type = Fragment_Memory::FM_CA;
    int index1 = 0;
    int index2 = 0;
    double native_distance = 0.0;
    int native_contacts = 0;
    int resindex1 = 0;
    int resindex2 = 0;
    double average_distance = 0.0;

    fprintf(nmer_output_file, "nmer1  nmer2 \tss \tnc \t<r>\n");
    for (int i=0; i<2*number_of_nmers; i++) {
      for (int j=0; j<number_of_nmers; j++) {
	index1 = i % number_of_nmers;
	index2 = j;
	// count native contacts
	native_contacts = 0;
	average_distance = 0.0;
	for (int q=0; q<amylometer_nmer_size; q++) {
	  if (i < number_of_nmers) {
	    resindex1 = index1 + q;
	    resindex2 = index2 + q;
	  }
	  else if (i >= number_of_nmers) {
	    resindex1 = index1 + q;
	    resindex2 = index2 + (amylometer_nmer_size - q - 1);
	  }
	  native_distance = native_structure->Rf(resindex1, iatom_type, resindex2, jatom_type);
	  average_distance += native_distance;
	  if (native_distance < amylometer_contact_cutoff && abs(resindex1-resindex2) > amylometer_nmer_size) {
	    native_contacts++;
	  }
	}
	average_distance /= double(amylometer_nmer_size);
	// mutate sequence
	// loop over each pair of nmers
	for (int k=0; k<n/(amylometer_nmer_size*2); k++) {
	  // loop within each pair of nmers
	  for (int l=0; l<amylometer_nmer_size*2; l++) {
	    // put a space between the nmers
	    if (l == amylometer_nmer_size) {
	      fprintf(nmer_output_file, " ");
	    }
	    // for the first nmer...
	    if (l < amylometer_nmer_size) {
	      // always put it going forward
	      se[k*2*amylometer_nmer_size+l]= nmer_array[i % number_of_nmers][l % amylometer_nmer_size];
	      // if this is your first time through the pair of nmers, write to nmer_output file
	      if (k == 0) {
		fprintf(nmer_output_file, "%c",se[k*2*amylometer_nmer_size+l]);
	      }
	    }
	    // for the second nmer...
	    else if (l >= amylometer_nmer_size) {
	      // if i is in the first half, write the second nmer forwards
	      if (i < number_of_nmers) {
		se[k*2*amylometer_nmer_size+l] = nmer_array[j][l % amylometer_nmer_size];
	      }
	      // if i is in the second half, write the second nmer backwards
	      else if (i >= number_of_nmers) {
		se[k*2*amylometer_nmer_size+l] = nmer_array[j][amylometer_nmer_size - 1 - (l % amylometer_nmer_size)];
	      }
	      // if this is your first time through the pair of nmers, write to nmer_output file
	      if (k == 0) {
		fprintf(nmer_output_file, "%c",se[k*2*amylometer_nmer_size+l]);
	      }
	    }
	  }
	}
	// done looping through the pair of nmers
	fprintf(nmer_output_file, "\t%3d \t%3d \t%3.1f\n", abs(index1-index2), native_contacts, average_distance);
	// print out whole sequence (for debugging)
	// for (int i=0; i<n; i++) {
	//   printf("%c", se[i]);
	// }
	// printf("\n");
	// calculate energy, output
	for (int m=0;m<nEnergyTerms;++m) energy_all[m] = 0.0; // clear energy values
	compute_backbone();                                   // compute energies
	// output energies to file
	for (int m=1;m<nEnergyTerms;++m) fprintf(amylometer_energy_file, "\t%8.6f", energy_all[m]);
	fprintf(amylometer_energy_file, "\t%8.6f\n", energy_all[ET_TOTAL]);
      }
    }
  }
  // give an error if an incorrect mode was used in fix_backbone_coeff.dat
  else {
    error->all(FLERR,"Amylometer: invalid amylometer mode\n");
  }

  fclose(amylometer_energy_file);
  return;
}

void FixBackbone::compute_optimization()
{
  // computes and writes out the energies for all interaction types 
  // for direct, protein-mediated, and water-mediated potentials

  double *xi, *xj, dx[3];
  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rij, rho_i, rho_j;
  double direct_energy, proteinmed_energy, watermed_energy;
  double direct_energies[20][20], protein_energies[20][20], water_energies[20][20];
  double contact_norm[20][20];

  direct_energy=0.0;
  proteinmed_energy=0.0;
  watermed_energy=0.0;

  // array initialization
  for (i=0;i<20;++i) {
    for (j=i;j<20;++j) {
      direct_energies[i][j] = direct_energies[j][i] = 0.0;
      protein_energies[i][j] = protein_energies[j][i] = 0.0;
      water_energies[i][j] = water_energies[j][i] = 0.0;      
      contact_norm[i][j] = contact_norm[j][i] = 0.0;
    }
  }
  
  // Double loop over all residue pairs
  for (i=0;i<n;++i) {
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = se_map[se[i_resno]-'A'];
    i_chno = chain_no[i]-1;
    
    rho_i = get_residue_density(i);
    
    for (j=i+1;j<n;++j) {
      // get information about residue j
      j_resno = res_no[j]-1;
      jres_type = se_map[se[j_resno]-'A'];
      j_chno = chain_no[j]-1;
      
      // Select beta atom unless the residue type is GLY, then select alpha carbon
      if (se[i_resno]=='G') { xi = xca[i]; }
      else { xi = xcb[i]; }
      if (se[j_resno]=='G') { xj = xca[j]; }
      else { xj = xcb[j]; }
	  
      // compute distance between the two atoms
      dx[0] = xi[0] - xj[0];
      dx[1] = xi[1] - xj[1];
      dx[2] = xi[2] - xj[2];
      rij = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      
      // if the atoms are within the threshold, compute the energies
      if (abs(i-j)>=contact_cutoff || i_chno != j_chno) {
	// rho_i = get_residue_density(i);
	rho_j = get_residue_density(j);
	// calculate direct, protein-mediated, water-mediated energies
	direct_energies[ires_type][jres_type] += compute_direct_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j); 
	protein_energies[ires_type][jres_type] += compute_proteinmed_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
	water_energies[ires_type][jres_type] += compute_watermed_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
	contact_norm[ires_type][jres_type] += 1.0;
      }
    }
  }

  for (i=0;i<20;++i) {
    for (j=i;j<20;++j) {
      direct_energies[i][j] = direct_energies[i][j] + direct_energies[j][i];
      protein_energies[i][j] = protein_energies[i][j] + protein_energies[j][i];
      water_energies[i][j] = water_energies[i][j] + water_energies[j][i];
      contact_norm[i][j] = contact_norm[i][j] + contact_norm[j][i];
    }
  }

  for (i=0;i<20;++i) {
    direct_energies[i][i] /= 2.0;
    protein_energies[i][i] /= 2.0;
    water_energies[i][i] /= 2.0;
    contact_norm[i][i] /= 2.0;
  }


  // write output calculated using native sequence on step 0
  // if step !=0 then write output calculated with shuffled sequence
  if (ntimestep == 0){
    for (i=0;i<20;++i) {
      for (j=i;j<20;++j) {
	fprintf(native_optimization_file,"%f %f %f \n", direct_energies[i][j], protein_energies[i][j], water_energies[i][j]);
	fprintf(native_optimization_norm_file,"%f \n", contact_norm[i][j]);
      }
    }
  }
  else {
    for (i=0;i<20;++i) {
      for (j=i;j<20;++j) {
	fprintf(optimization_file,"%f %f %f \n", direct_energies[i][j], protein_energies[i][j], water_energies[i][j]);
	fprintf(optimization_norm_file,"%f \n", contact_norm[i][j]);
      }
    }
  }
}

void FixBackbone::compute_average_sequence_optimization()
{
  // computes and writes out the energies for all interaction types 
  // for direct, protein-mediated, and water-mediated potentials
  // uses an "average sequence" computed from a multiple sequence alignment
  // to spread the contribution to the counts across multiple residue pairs

  double *xi, *xj, dx[3];
  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rij, rho_i, rho_j;
  double direct_energy, proteinmed_energy, watermed_energy;
  double direct_energies[20][20], protein_energies[20][20], water_energies[20][20];
  double contact_norm[20][20];
  double average_sequence_factor;
  int itype, jtype;

  direct_energy=0.0;
  proteinmed_energy=0.0;
  watermed_energy=0.0;

  // array initialization
  for (i=0;i<20;++i) {
    for (j=i;j<20;++j) {
      direct_energies[i][j] = direct_energies[j][i] = 0.0;
      protein_energies[i][j] = protein_energies[j][i] = 0.0;
      water_energies[i][j] = water_energies[j][i] = 0.0;      
      contact_norm[i][j] = contact_norm[j][i] = 0.0;
    }
  }
  
  // Double loop over all residue pairs
  for (i=0;i<n;++i) {
    // Loop over all 20 residue types for residue i
    for (itype=0;itype<20;++itype) {
      // get information about residue i
      i_resno = res_no[i]-1;
      if (average_sequence[i_resno][itype] == 0.0) {
	continue;
      }
      ires_type = itype;
      i_chno = chain_no[i]-1;
      
      rho_i = get_residue_density(i);
    
      for (j=i+1;j<n;++j) {
	// Loop over all 20 residue types for residue j
	for (jtype=0;jtype<20;++jtype) {
	  // get information about residue j
	  j_resno = res_no[j]-1;
	  if (average_sequence[j_resno][jtype] == 0.0) {
	    continue;
	  }
	  jres_type = jtype;
	  j_chno = chain_no[j]-1;
	  
	  // Select beta atom unless the residue type is GLY, then select alpha carbon
	  if (itype==9) { xi = xca[i]; }
	  else { xi = xcb[i]; }
	  if (jtype==9) { xj = xca[j]; }
	  else { xj = xcb[j]; }
	  
	  // compute distance between the two atoms
	  dx[0] = xi[0] - xj[0];
	  dx[1] = xi[1] - xj[1];
	  dx[2] = xi[2] - xj[2];
	  rij = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	  
	  // if the atoms are within the threshold, compute the energies
	  if (abs(i-j)>=contact_cutoff || i_chno != j_chno) {
	    // rho_i = get_residue_density(i);
	    rho_j = get_residue_density(j);
	    average_sequence_factor = average_sequence[i_resno][itype]*average_sequence[j_resno][jtype];
	    // calculate direct, protein-mediated, water-mediated energies
	    direct_energies[ires_type][jres_type] += average_sequence_factor*compute_direct_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j); 
	    protein_energies[ires_type][jres_type] += average_sequence_factor*compute_proteinmed_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
	    water_energies[ires_type][jres_type] += average_sequence_factor*compute_watermed_energy(rij, i_resno, j_resno, ires_type, jres_type, rho_i, rho_j);
	    contact_norm[ires_type][jres_type] += average_sequence_factor*1.0;
	  }
	}
      }
    }
  }
  
  for (i=0;i<20;++i) {
    for (j=i;j<20;++j) {
      direct_energies[i][j] = direct_energies[i][j] + direct_energies[j][i];
      protein_energies[i][j] = protein_energies[i][j] + protein_energies[j][i];
      water_energies[i][j] = water_energies[i][j] + water_energies[j][i];
      contact_norm[i][j] = contact_norm[i][j] + contact_norm[j][i];
    }
  }

  for (i=0;i<20;++i) {
    direct_energies[i][i] /= 2.0;
    protein_energies[i][i] /= 2.0;
    water_energies[i][i] /= 2.0;
    contact_norm[i][i] /= 2.0;
  }


  // write output calculated energies
  for (i=0;i<20;++i) {
    for (j=i;j<20;++j) {
      fprintf(average_sequence_optimization_file,"%f %f %f \n", direct_energies[i][j], protein_energies[i][j], water_energies[i][j]);
      fprintf(average_sequence_optimization_norm_file,"%f \n", contact_norm[i][j]);
    }
  }
}

void FixBackbone::shuffler()
{
  double residue_density_i;
  double residue_density_j;

  // If doing a normal, full shuffling
  if (strcmp(shuffler_mode, "normal")==0) {
    // sequence shuffler
    for (int i=0; i<n; i++) {
      int r = i + (rand() % (n-i)); // Random remaining position.
      int temp = se[i]; se[i] = se[r]; se[r] = temp;
    }  
  }
  else if (strcmp(shuffler_mode, "burial")==0) {
    // burial-constrained sequence shuffler
    for (int shuffle_iteration=0; shuffle_iteration<1000; shuffle_iteration++) {
      for (int i=0; i<n; i++) {
	residue_density_i = get_residue_density(i);
	int j = i + (rand() % (n-i)); // Random remaining position.
	residue_density_j = get_residue_density(j);
	if ((residue_density_i > burial_ro_min[0] && residue_density_i < burial_ro_max[0] && residue_density_j > burial_ro_min[0] && residue_density_j < burial_ro_max[0]) ||
	    (residue_density_i > burial_ro_min[1] && residue_density_i < burial_ro_max[1] && residue_density_j > burial_ro_min[1] && residue_density_j < burial_ro_max[1]) ||
	    (residue_density_i > burial_ro_min[2] && residue_density_i < burial_ro_max[2] && residue_density_j > burial_ro_min[2] && residue_density_j < burial_ro_max[2])) {
	  // if criterion is met, swap residues
	  int temp = se[i]; se[i] = se[j]; se[j] = temp;
	}  
      }
    }
  }
  else {
    printf("Unrecognized shuffler mode %s\n", shuffler_mode);
  }
}

double FixBackbone::compute_direct_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j)
{
  double water_gamma_0_direct, water_gamma_1_direct;
  double sigma_wat, sigma_prot, sigma_gamma_direct, sigma_gamma_mediated;
  double t_min_direct, t_max_direct, theta_direct, t_min_mediated, t_max_mediated, theta_mediated;
  double direct_energy;

  if(abs(i_resno-j_resno)<contact_cutoff) return 0.0;

  // grab direct contact gamma for ires_type and jres_type
  // water_gamma_0_direct and water_gamma_1_direct should be equivalent (relic of compute_water_potential)
  water_gamma_0_direct = get_water_gamma(i_resno, j_resno, 0, ires_type, jres_type, 0);
  water_gamma_1_direct = get_water_gamma(i_resno, j_resno, 0, ires_type, jres_type, 1);
  sigma_gamma_direct = (water_gamma_0_direct + water_gamma_1_direct)/2;

  // compute theta function (eq. 9 in AWSEM SI)
  t_min_direct = tanh( well->par.kappa*(rij - well->par.well_r_min[0]) );
  t_max_direct = tanh( well->par.kappa*(well->par.well_r_max[0] - rij) );
  theta_direct = 0.25*(1.0 + t_min_direct)*(1.0 + t_max_direct);

  direct_energy = -epsilon*k_water*(sigma_gamma_direct*theta_direct);

  return direct_energy;
}

double FixBackbone::compute_proteinmed_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j)
{
  double water_gamma_prot_mediated;
  double sigma_wat, sigma_prot;
  double t_min_mediated, t_max_mediated, theta_mediated;
  double proteinmed_energy;

  if(abs(i_resno-j_resno)<contact_cutoff) return 0.0;

  water_gamma_prot_mediated = get_water_gamma(i_resno, j_resno, 1, ires_type, jres_type, 0);

  // compute sigma_prot
  sigma_wat = 0.25*(1.0 - tanh(well->par.kappa_sigma*(rho_i-well->par.treshold)))*(1.0 - tanh(well->par.kappa_sigma*(rho_j-well->par.treshold)));
  sigma_prot = 1.0 - sigma_wat;
  
  // compute theta function
  t_min_mediated = tanh( well->par.kappa*(rij - well->par.well_r_min[1]) );
  t_max_mediated = tanh( well->par.kappa*(well->par.well_r_max[1] - rij) );
  theta_mediated = 0.25*(1.0 + t_min_mediated)*(1.0 + t_max_mediated);

  proteinmed_energy = -epsilon*k_water*sigma_prot*water_gamma_prot_mediated*theta_mediated;

  return proteinmed_energy;
}

double FixBackbone::compute_watermed_energy(double rij, int i_resno, int j_resno, int ires_type, int jres_type, double rho_i, double rho_j)
{
  double water_gamma_wat_mediated;
  double sigma_wat;
  double t_min_mediated, t_max_mediated, theta_mediated;
  double watermed_energy;  

  if(abs(i_resno-j_resno)<contact_cutoff) return 0.0;

  water_gamma_wat_mediated = get_water_gamma(i_resno, j_resno, 1, ires_type, jres_type, 1);

  // compute sigma_wat
  sigma_wat = 0.25*(1.0 - tanh(well->par.kappa_sigma*(rho_i-well->par.treshold)))*(1.0 - tanh(well->par.kappa_sigma*(rho_j-well->par.treshold)));

  // compute theta function
  t_min_mediated = tanh( well->par.kappa*(rij - well->par.well_r_min[1]) );
  t_max_mediated = tanh( well->par.kappa*(well->par.well_r_max[1] - rij) );
  theta_mediated = 0.25*(1.0 + t_min_mediated)*(1.0 + t_max_mediated);

  watermed_energy = -epsilon*k_water*sigma_wat*water_gamma_wat_mediated*theta_mediated;

  return watermed_energy;
}

void FixBackbone::compute_burial_optimization()
{
  // computes and writes out the energies for all interaction types 
  // for direct, protein-mediated, and water-mediated potentials

  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rho_i, rho_j;
  double burial_energy;
 
  double t[3][2];
  double burial_gamma_0, burial_gamma_1, burial_gamma_2;
  double norm_array[20], burial_array[3][20];

  // array initialization
  for (i=0;i<20;++i) {
    norm_array[i]=0.0;
    for (j=0;j<3;++j) {
      burial_array[j][i]=0.0;
    }
  }
    
  for (i=0;i<n;++i) {
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = se_map[se[i_resno]-'A'];
    i_chno = chain_no[i]-1;
    
    rho_i = get_residue_density(i);
    
    t[0][0] = tanh( burial_kappa*(rho_i - burial_ro_min[0]) );
    t[0][1] = tanh( burial_kappa*(burial_ro_max[0] - rho_i) );
    t[1][0] = tanh( burial_kappa*(rho_i - burial_ro_min[1]) );
    t[1][1] = tanh( burial_kappa*(burial_ro_max[1] - rho_i) );
    t[2][0] = tanh( burial_kappa*(rho_i - burial_ro_min[2]) );
    t[2][1] = tanh( burial_kappa*(burial_ro_max[2] - rho_i) );
	
    burial_gamma_0 = get_burial_gamma(i_resno, ires_type, 0);
    burial_gamma_1 = get_burial_gamma(i_resno, ires_type, 1);
    burial_gamma_2 = get_burial_gamma(i_resno, ires_type, 2);
    
    burial_array[0][ires_type] +=-0.5*epsilon*k_burial*burial_gamma_0*(t[0][0] + t[0][1]);
    burial_array[1][ires_type] += -0.5*epsilon*k_burial*burial_gamma_1*(t[1][0] + t[1][1]);
    burial_array[2][ires_type] += -0.5*epsilon*k_burial*burial_gamma_2*(t[2][0] + t[2][1]);
    
    norm_array[ires_type] += 1.0;
  }
  
  // write output calculated using native sequence on step 0
  // if step !=0 then write output calculated with shuffled sequence
  if (ntimestep == 0){
    for (i=0;i<20;++i) {
      fprintf(native_burial_optimization_file,"%f %f %f \n", burial_array[0][i], burial_array[1][i], burial_array[2][i]);
      fprintf(burial_optimization_norm_file,"%f \n",norm_array[i]);
    }
  }
  else {
    for (i=0;i<20;++i) {
      fprintf(burial_optimization_file,"%f %f %f \n", burial_array[0][i], burial_array[1][i], burial_array[2][i]);
    }
  }
}

// Mutate the sequence
void FixBackbone::mutate_sequence()
{
  // Copy the new sequence into the se array
  strcpy(se,mutate_sequence_sequences[mutate_sequence_sequence_index]);
  // Increment the sequence index
  mutate_sequence_sequence_index++;
}

// Output the number of contacts that each residue in the protein is currently making
void FixBackbone::output_per_residue_contacts()
{
  int i, j;
  int ires_type, jres_type, i_resno, j_resno, i_chno, j_chno;
  double rij;
  int num_i_contacts;
  
  // output the timestep in the first column
  fprintf(output_per_residue_contacts_file, "%8d ", ntimestep);
  // Double loop over all residue pairs
  for (i=0;i<n;++i) {
    // Zero out the number of contacts
    num_i_contacts = 0;
    // get information about residue i
    i_resno = res_no[i]-1;
    ires_type = get_residue_type(i_resno);
    i_chno = chain_no[i]-1;
    
    // We do want to double count contacts in this case, so start at j=0
    for (j=0;j<n;++j) {
      // get information about residue j
      j_resno = res_no[j]-1;
      jres_type = get_residue_type(j_resno);
      j_chno = chain_no[j]-1;
      
      // get the distance between i and j
      rij = get_residue_distance(i_resno, j_resno);

      // if the atoms are within the threshold, add one to the number of contacts for residue i
      if (rij < output_per_residue_contacts_rij_threshold && (abs(i-j)>=output_per_residue_contacts_min_seq_sep || i_chno != j_chno)) {
	num_i_contacts++;
      }
    }
    // output the number of contacts for residue i
    fprintf(output_per_residue_contacts_file, "%3d ", num_i_contacts);
  }
  // print out a new line after looping over all the residues
  fprintf(output_per_residue_contacts_file, "\n", num_i_contacts);
  
}


void FixBackbone::print_forces(int coord)
{
  int index;

  if (coord==1) {
    fprintf(dout, "rca = {");
    for (int i=0;i<nn;i++) {
      index = alpha_carbons[i];
      if (index!=-1) {
	fprintf(dout, "{%.8f, %.8f, %.8f}", x[index][0], x[index][1], x[index][2]);
	if (i!=nn-1) fprintf(dout, ",\n");
      }
    }
    fprintf(dout, "};\n\n");

    fprintf(dout, "rcb = {");
    for (int i=0;i<nn;i++) {
      index = beta_atoms[i];
      if (index!=-1) {
	fprintf(dout, "{%.8f, %.8f, %.8f}", x[index][0], x[index][1], x[index][2]);
	if (i!=nn-1) fprintf(dout, ",\n");
      }
    }
    fprintf(dout, "};\n\n");

    fprintf(dout, "ro = {");
    for (int i=0;i<nn;i++) {
      index = oxygens[i];
      if (index!=-1) {
	fprintf(dout, "{%.8f, %.8f, %.8f}", x[index][0], x[index][1], x[index][2]);
	if (i!=nn-1) fprintf(dout, ",\n");
      }
    }
    fprintf(dout, "};\n\n");

    fprintf(dout, "rn = {");
    for (int i=0;i<nn;i++) {
      fprintf(dout, "{%.8f, %.8f, %.8f}", xn[i][0], xn[i][1], xn[i][2]);
      if (i!=nn-1) fprintf(dout, ",\n");
    }
    fprintf(dout, "};\n\n");

    fprintf(dout, "rcp = {");
    for (int i=0;i<nn;i++) {
      fprintf(dout, "{%.8f, %.8f, %.8f}", xcp[i][0], xcp[i][1], xcp[i][2]);
      if (i!=nn-1) fprintf(dout, ",\n");
    }
    fprintf(dout, "};\n\n");
        
    fprintf(dout, "rh = {");
    for (int i=0;i<nn;i++) {
      fprintf(dout, "{%.8f, %.8f, %.8f}", xh[i][0], xh[i][1], xh[i][2]);
      if (i!=nn-1) fprintf(dout, ",\n");
    }
    fprintf(dout, "};\n\n\n");
  }
	
  fprintf(dout, "fca = {");
  for (int i=0;i<nn;i++) {
    index = alpha_carbons[i];
    if (index!=-1) {
      fprintf(dout, "{%.8f, %.8f, %.8f}", f[index][0], f[index][1], f[index][2]);
      if (i!=nn-1) fprintf(dout, ",\n");
    }
  }
  fprintf(dout, "};\n\n");

  fprintf(dout, "fcb = {");
  for (int i=0;i<nn;i++) {
    index = beta_atoms[i];
    if (index!=-1) {
      fprintf(dout, "{%.8f, %.8f, %.8f}", f[index][0], f[index][1], f[index][2]);
      if (i!=nn-1) fprintf(dout, ",\n");
    }
  }
  fprintf(dout, "};\n\n");

  fprintf(dout, "fo = {");
  for (int i=0;i<nn;i++) {
    index = oxygens[i];
    if (index!=-1) {
      fprintf(dout, "{%.8f, %.8f, %.8f}", f[index][0], f[index][1], f[index][2]);
      if (i!=nn-1) fprintf(dout, ",\n");
    }
  }
  fprintf(dout, "};\n\n\n\n");
}

void FixBackbone::compute_backbone()
{     
  ntimestep = update->ntimestep;

  if(atom->nlocal==0) return;
	
  if (comm->nprocs>1 || ntimestep==0)
    Construct_Computational_Arrays();

  x = atom->x;
  f = atom->f;
  image = atom->image;

  int i, j, xbox, ybox, zbox;
  int i_resno, j_resno;
  int i_chno, j_chno;
	
  for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;
  force_flag = 0;

  for (i=0;i<nn;++i) {		
    if ( (res_info[i]==LOCAL || res_info[i]==GHOST) ) {
      if (domain->xperiodic) {
	xbox = (image[alpha_carbons[i]] & 1023) - 512;
	xca[i][0] = x[alpha_carbons[i]][0] + xbox*prd[0];
      } else xca[i][0] = x[alpha_carbons[i]][0];
      if (domain->yperiodic) {
	ybox = (image[alpha_carbons[i]] >> 10 & 1023) - 512;
	xca[i][1] = x[alpha_carbons[i]][1] + ybox*prd[1];
      } else xca[i][1] = x[alpha_carbons[i]][1];
      if (domain->zperiodic) {
	zbox = (image[alpha_carbons[i]] >> 20) - 512;
	xca[i][2] = x[alpha_carbons[i]][2] + zbox*prd[2];
      } else xca[i][2] = x[alpha_carbons[i]][2];
					
      if (beta_atoms[i]!=-1) {
	if (domain->xperiodic) {
	  xbox = (image[beta_atoms[i]] & 1023) - 512;
	  xcb[i][0] = x[beta_atoms[i]][0] + xbox*prd[0];
	} else xcb[i][0] = x[beta_atoms[i]][0];
	if (domain->yperiodic) {
	  ybox = (image[beta_atoms[i]] >> 10 & 1023) - 512;
	  xcb[i][1] = x[beta_atoms[i]][1] + ybox*prd[1];
	} else xcb[i][1] = x[beta_atoms[i]][1];
	if (domain->zperiodic) {
	  zbox = (image[beta_atoms[i]] >> 20) - 512;
	  xcb[i][2] = x[beta_atoms[i]][2] + zbox*prd[2];
	} else xcb[i][2] = x[beta_atoms[i]][2];	
      }
		
      if (oxygens[i]!=-1) {
	if (domain->xperiodic) {
	  xbox = (image[oxygens[i]] & 1023) - 512;
	  xo[i][0] = x[oxygens[i]][0] + xbox*prd[0];
	} else xo[i][0] = x[oxygens[i]][0];
	if (domain->yperiodic) {
	  ybox = (image[oxygens[i]] >> 10 & 1023) - 512;
	  xo[i][1] = x[oxygens[i]][1] + ybox*prd[1];
	} else xo[i][1] = x[oxygens[i]][1];
	if (domain->zperiodic) {
	  zbox = (image[oxygens[i]] >> 20) - 512;
	  xo[i][2] = x[oxygens[i]][2] + zbox*prd[2];
	} else xo[i][2] = x[oxygens[i]][2];
      }
    }

    i_resno=res_no[i]-1;
    int im1 = res_no_l[i_resno-1];	
    if ( im1!=-1 && i_resno>0 && !isFirst(i) && (res_info[i]==LOCAL || res_info[i]==GHOST) && (res_info[im1]==LOCAL || res_info[im1]==GHOST) ) {
      xn[i][0] = an*xca[im1][0] + bn*xca[i][0] + cn*xo[im1][0];
      xn[i][1] = an*xca[im1][1] + bn*xca[i][1] + cn*xo[im1][1];
      xn[i][2] = an*xca[im1][2] + bn*xca[i][2] + cn*xo[im1][2];

      xh[i][0] = ah*xca[im1][0] + bh*xca[i][0] + ch*xo[im1][0];
      xh[i][1] = ah*xca[im1][1] + bh*xca[i][1] + ch*xo[im1][1];
      xh[i][2] = ah*xca[im1][2] + bh*xca[i][2] + ch*xo[im1][2];
    } else {
      xn[i][0] = xn[i][1] = xn[i][2] = 0.0;

      xh[i][0] = xh[i][1] = xh[i][2] = 0.0;
    }
		
    if ( im1!=-1 && i_resno>0) {
      if (!isFirst(i) && (res_info[i]==LOCAL || res_info[i]==GHOST) && (res_info[im1]==LOCAL || res_info[im1]==GHOST) ) {
	xcp[im1][0] = ap*xca[im1][0] + bp*xca[i][0] + cp*xo[im1][0];
	xcp[im1][1] = ap*xca[im1][1] + bp*xca[i][1] + cp*xo[im1][1];
	xcp[im1][2] = ap*xca[im1][2] + bp*xca[i][2] + cp*xo[im1][2];
      } else if (im1!=-1){
	xcp[im1][0] = xcp[im1][1] = xcp[im1][2] = 0.0;
      }
    }

  }
  xcp[nn-1][0] = xcp[nn-1][1] = xcp[nn-1][2] = 0.0;

#ifdef DEBUGFORCES

  if (ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "AtStart: %d\n", ntimestep);
    fprintf(dout, "Number of residues %d\n", n);
    fprintf(dout, "Local Number of residues %d\n\n", nn);
    print_forces(1);
  }
	
  timerBegin();
  
  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    if (chain_flag && res_info[i]==LOCAL)
      compute_chain_potential(i);
  }
	
  if (chain_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Chain: %d\n", ntimestep);
    fprintf(dout, "Chain_Energy: %f\n\n", energy[ET_CHAIN]);
    print_forces();
  }
	
  timerEnd(TIME_CHAIN);
	
  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    if (!isFirst(i) && !isLast(i) && chi_flag && res_info[i]==LOCAL && se[i_resno]!='G')
      compute_chi_potential(i);
  }
	
  if (chi_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Chi: %d\n", ntimestep);
    fprintf(dout, "Chi_Energy: %f\n\n", energy[ET_CHI]);
    print_forces();
  }
	
  timerEnd(TIME_CHI);
	
  for (i=0;i<nn;i++) {
    if (shake_flag && res_info[i]==LOCAL)
      compute_shake(i);
  }
	
  timerEnd(TIME_SHAKE);

  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    if (!isFirst(i) && !isLast(i) && rama_flag && res_info[i]==LOCAL && se[i_resno]!='G')
      compute_rama_potential(i);			
  }
	
  if (rama_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Rama: %d\n", ntimestep);
    fprintf(dout, "Rama_Energy: %f\n\n", energy[ET_RAMA]);
    print_forces();
  }
	
  timerEnd(TIME_RAMA);

  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      if (i_resno-2>=0 && i_resno+2<n && j_resno-2>=0 && j_resno+2<n && !isLast(i) && !isFirst(j) && ( i_chno!=j_chno || abs(j_resno-i_resno)>2 ) && dssp_hdrgn_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) && se[j_resno]!='P')
	//			if (!isLast(i) && !isFirst(j) && abs(j_resno-i_resno)>2 && dssp_hdrgn_flag && res_info[i]==LOCAL && res_info[j]==LOCAL && se[j_resno]!='P')
	compute_dssp_hdrgn(i, j);
    }
  }
	
  if (dssp_hdrgn_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "DSSP: %d\n", ntimestep);
    fprintf(dout, "DSSP_Energy: %f\n\n", energy[ET_DSSP]);
    print_forces();
  }
	
  timerEnd(TIME_DSSP);

  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      //if (i<n-i_med_min && j>=i+i_med_min && p_ap_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j] == GHOST))
      if (j_resno>i_resno+i_med_min && p_ap_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_P_AP_potential(i, j);
    }
  }
	
  if (p_ap_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "P_AP: %d\n", ntimestep);
    fprintf(dout, "P_AP_Energy: %f\n\n", energy[ET_PAP]);
    print_forces();
  }
	
  timerEnd(TIME_PAP);

  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      //if (water_flag && ( i_chno!=j_chno || j_resno-i_resno>=contact_cutoff ) && res_info[i]==LOCAL)
      //if (water_flag && j_resno-i_resno>=contact_cutoff && res_info[i]==LOCAL)
      if (water_flag && ( (i_chno!=j_chno && j_resno > i_resno ) || ( i_chno == j_chno && j_resno-i_resno>=contact_cutoff) ) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_water_potential(i, j);
    }
  }
	
  if (water_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Water: %d\n", ntimestep);
    fprintf(dout, "Water_Energy: %f\n\n", energy[ET_WATER]);
    print_forces();
  }
	
  timerEnd(TIME_WATER);
	
  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      if (frag_mem_tb_flag && j_resno-i_resno>=fm_gamma->minSep() && (fm_gamma->maxSep()==-1 || j_resno-i_resno<=fm_gamma->maxSep()) && chain_no[i]==chain_no[j] && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) )
	table_fragment_memory(i, j);
    }
  }
	
  if (frag_mem_tb_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Table_Frag_Mem: %d\n", ntimestep);
    fprintf(dout, "Table_Frag_Mem_Energy: %f\n\n", energy[ET_FRAGMEM]);
    print_forces();
  }
	
  timerEnd(TIME_FRAGMEM);

  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      if (ssb_flag && ( i_chno!=j_chno || j_resno-i_resno>=ssb_ij_sep ) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	//			if (ssb_flag && j_resno-i_resno>=ssb_ij_sep && res_info[i]==LOCAL)
	compute_solvent_barrier(i, j);
    }
  }

  if (ssb_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "SSB: %d\n", ntimestep);
    fprintf(dout, "SSB_Energy: %f\n\n", energy[ET_SSB]);
    print_forces();
  }
  
  timerEnd(TIME_SSB);

  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      if (huckel_flag) {
        compute_DebyeHuckel_Interaction(i, j);
      }
    }
  }
	
  if (huckel_flag && ntimestep >=sStep && ntimestep <=eStep) {
    fprintf(dout, "DH: %d\n", ntimestep);
    fprintf(dout, "DH_Elect_Energy: %f\n\n", energy[ET_DH]);
    print_forces();
  }

  timerEnd(TIME_DH);
	
  for (i=0;i<nn;i++) {
    if (burial_flag && res_info[i]==LOCAL)
      compute_burial_potential(i);
  }
	
  if (burial_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Burial: %d\n", ntimestep);
    fprintf(dout, "Burial_Energy: %f\n\n", energy[ET_BURIAL]);
    print_forces();
  }
	
  timerEnd(TIME_BURIAL);

  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    if (helix_flag && i_resno<(ch_pos[i_chno]+ch_len[i_chno]-1)-helix_i_diff-1 && res_info[i]==LOCAL)
      //		if (helix_flag && i<nn-helix_i_diff-1 && i_resno==res_no[i+helix_i_diff]-helix_i_diff && res_info[i]==LOCAL)
      compute_helix_potential(i);
  }
	
  if (helix_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Helix: %d\n", ntimestep);
    fprintf(dout, "Helix_Energy: %f\n\n", energy[ET_HELIX]);
    print_forces();
  }
	
  timerEnd(TIME_HELIX);

  for (i=0;i<nn;i++) {
    if (memb_flag && res_info[i]==LOCAL)
      compute_membrane_potential(i);
  }

  if (memb_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Membrane: %d\n", ntimestep);
    fprintf(dout, "Membrane_Energy: %f\n\n", energy[ET_MEMB]);
    print_forces();
  }

  timerEnd(TIME_MEMB);

	
  for (i=0;i<nn;i++) {
    if (frag_mem_flag && res_info[i]==LOCAL)
      compute_fragment_memory_potential(i);
  }
	
  if (frag_mem_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Frag_Mem: %d\n", ntimestep);
    fprintf(dout, "Frag_Mem_Energy: %f\n\n", energy[ET_FRAGMEM]);
    print_forces();
  }
	
  timerEnd(TIME_FRAGMEM);
  
  for (i=0;i<nn;i++) {
    if (vec_frag_mem_flag && res_info[i]==LOCAL)
      compute_vector_fragment_memory_potential(i);
  }
	
  if (vec_frag_mem_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Vec_Frag_Mem: %d\n", ntimestep);
    fprintf(dout, "Vec_Frag_Mem_Energy: %f\n\n", energy[ET_VFRAGMEM]);
    print_forces();
  }
	
  timerEnd(TIME_VFRAGMEM);
  
  if (amh_go_flag)
    compute_amh_go_model();
    	
  if (amh_go_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "AMH-Go: %d\n", ntimestep);
    fprintf(dout, "AMH-Go_Energy: %f\n\n", energy[ET_AMHGO]);
    print_forces();
  }
	
  timerEnd(TIME_AMHGO);

  if (excluded_flag)
    compute_excluded_volume();

  if (p_excluded_flag)
    compute_p_degree_excluded_volume();

  if (r6_excluded_flag)
    compute_r6_excluded_volume();
	
  timerEnd(TIME_VEXCLUDED);
	
  if (ntimestep>=sStep && ntimestep<=eStep)
    fprintf(dout, "\n\n");

#else
  
  for (i=0;i<nn;i++) {
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    
    if (chain_flag && res_info[i]==LOCAL)
      compute_chain_potential(i); //i_resno is the global residue index in parallelization, i is local

    if (!isFirst(i) && !isLast(i) && chi_flag && res_info[i]==LOCAL && se[i_resno]!='G')
      compute_chi_potential(i);

    if (shake_flag && res_info[i]==LOCAL)
      compute_shake(i);

    if (!isFirst(i) && !isLast(i) && rama_flag && res_info[i]==LOCAL && se[i_resno]!='G')
      compute_rama_potential(i);

    if (memb_flag && res_info[i]==LOCAL)
      compute_membrane_potential(i);

    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      		
      //if (!isLast(i) && !isFirst(j) && ( i_chno!=j_chno || abs(j_resno-i_resno)>2 ) && dssp_hdrgn_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) && se[j_resno]!='P')
      if (i_resno-2>=0 && i_resno+2<n && j_resno-2>=0 && j_resno+2<n && !isLast(i) && !isFirst(j) && ( i_chno!=j_chno || abs(j_resno-i_resno)>2 ) && dssp_hdrgn_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) && se[j_resno]!='P')
	compute_dssp_hdrgn(i, j);
				
      // Need to change
      if (j_resno>i_resno+i_med_min && p_ap_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_P_AP_potential(i, j);

      //if (water_flag && ( i_chno!=j_chno || j_resno-i_resno>=contact_cutoff ) && res_info[i]==LOCAL)
      if (water_flag && ( (i_chno!=j_chno && j_resno > i_resno ) || ( i_chno == j_chno && j_resno-i_resno>=contact_cutoff) ) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_water_potential(i, j);
			  
      if (frag_mem_tb_flag && j_resno-i_resno>=fm_gamma->minSep() && (fm_gamma->maxSep()==-1 || j_resno-i_resno<=fm_gamma->maxSep()) && chain_no[i]==chain_no[j] && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) )
	table_fragment_memory(i, j);

      if (ssb_flag && ( i_chno!=j_chno || j_resno-i_resno>=ssb_ij_sep ) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_solvent_barrier(i, j);

      if (huckel_flag && (j > i || i_chno!=j_chno) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
        compute_DebyeHuckel_Interaction(i, j);
    }
		
    if (burial_flag && res_info[i]==LOCAL)
      compute_burial_potential(i);
    	
    //    if (helix_flag && i<nn-helix_i_diff-1 && i_resno==res_no[i+helix_i_diff]-helix_i_diff && res_info[i]==LOCAL)
    if (helix_flag && i_resno<(ch_pos[i_chno]+ch_len[i_chno]-1)-helix_i_diff-1 && res_info[i]==LOCAL)
      compute_helix_potential(i);
			
    if (frag_mem_flag && res_info[i]==LOCAL)
      compute_fragment_memory_potential(i);
      
    if (vec_frag_mem_flag && res_info[i]==LOCAL)
      compute_vector_fragment_memory_potential(i);

  }

  // if the fragment frustratometer is on and it is time to compute the fragment frustration, do so!
  if (frag_frust_flag && ntimestep % frag_frust_output_freq == 0) {
    // if running in shuffle mode...
    if (frag_frust_shuffle_flag) {
      // zero decoy energy array (because new decoy energies are generated at every calculation in "shuffle" mode)
      for (i=0; i<n; i++) {
	for (j=0; j<num_decoy_calcs; j++) {
	  decoy_energy[i][j] = 0.0;
	}
      }
      // loop over decoys to be generated
      for (int idecoy=0; idecoy<num_decoy_calcs; idecoy++) {
	// compute and store decoy energies for each residue
	for (i=0;i<n;i++) {
	  compute_decoy_memory_potential(i,idecoy);
	}
	// randomize decoy memory positions
	randomize_decoys();
      }
    }
    // if running in read mode ...
    else if (frag_frust_read_flag) {
      // zero native energies (because only the native energies are recalculated at every calculation in "read" mode)
      for (i=0;i<n;i++) {
	decoy_energy[i][0] = 0.0;
      }
      // compute and store native energies
      for (i=0;i<n;i++) {
	compute_decoy_memory_potential(i,0);
      }
      // on time step zero, compute the decoy energy distribution
      if (ntimestep == 0) {
	compute_generated_decoy_energies();
      }
    }
    else {
      // throw an error because only shuffle and read are valid modes
      error->all(FLERR,"Fragment_Frustratometer: only shuffle and read are valid modes.");
    }
    // regardless of what mode you are using, calculate and output per-residue frustration index
    compute_fragment_frustration();
  }

  // if it is time to compute the tertiary frustration, do it
  if (tert_frust_flag && ntimestep % tert_frust_output_freq == 0) {
    fprintf(tert_frust_output_file,"# timestep: %d\n", ntimestep);
    fprintf(tert_frust_vmd_script,"# timestep: %d\n", ntimestep);
    if (strcmp(tert_frust_mode, "configurational")==0 || strcmp(tert_frust_mode, "mutational")==0) {
      compute_tert_frust();
    }
    else if (strcmp(tert_frust_mode, "singleresidue")==0) {
      compute_tert_frust_singleresidue();
    }
  }

  // if it is time to compute the nmer frustration, do it
  if (nmer_frust_flag && ntimestep % nmer_frust_output_freq == 0) {
    fprintf(nmer_frust_output_file,"# timestep: %d\n", ntimestep);
    fprintf(nmer_frust_vmd_script,"# timestep: %d\n", ntimestep);
    if (strcmp(nmer_frust_mode, "pairwise")==0) {
      compute_nmer_frust();
    }
    else if (strcmp(nmer_frust_mode, "singlenmer")==0) {
      compute_singlenmer_frust();
    }
  }

  // if it is time to output the selection temperature file, do it
  if (selection_temperature_flag && ntimestep % selection_temperature_output_frequency == 0) {
    if (selection_temperature_output_interaction_energies_flag) {
      fprintf(selection_temperature_file,"# timestep: %d\n", ntimestep);
    }
    output_selection_temperature_data();
  }

  // if it is time to do the mcso, do it
  if (monte_carlo_seq_opt_flag) {
    compute_mcso();
  }

  // if it is time to output energies for contact potential optimization DO IT
  if (optimization_flag && ntimestep % optimization_output_freq == 0) {
    compute_optimization();
  }
  // if it is time to output energies for burial potential optimization DO IT
  if (burial_optimization_flag && ntimestep % burial_optimization_output_freq == 0) {
    compute_burial_optimization();
  }
  // if it is time to output energies for DebyeHuckel potential optimization DO IT
  if (debyehuckel_optimization_flag && ntimestep % debyehuckel_optimization_output_freq == 0) {
    compute_debyehuckel_optimization();
  }
  // if it is time to output energies for average sequence optimization DO IT
  if (average_sequence_optimization_flag && ntimestep % average_sequence_optimization_output_freq == 0) {
    compute_average_sequence_optimization();
  }
  // if collecting energies for optimization, shuffle the sequence.  (native sequence used on step 0)
  if ((optimization_flag || burial_optimization_flag || debyehuckel_optimization_flag) && (shuffler_flag)){
    shuffler();
  }
  // if mutating sequence to evaluate energy of mutants, call function to mutate the sequence
  if (mutate_sequence_flag && ntimestep != update->laststep) {
    mutate_sequence();
  }
  if (amh_go_flag)
    compute_amh_go_model();

  if (excluded_flag)
    compute_excluded_volume();

  if (p_excluded_flag)
    compute_p_degree_excluded_volume();

  if (r6_excluded_flag)
    compute_r6_excluded_volume();

  if (output_per_residue_contacts_flag && ntimestep % output_per_residue_contacts_frequency == 0) {
    output_per_residue_contacts();
  }

#endif
	
  for (int i=1;i<nEnergyTerms;++i) energy[ET_TOTAL] += energy[i];
	
  if (ntimestep%output->thermo_every==0) {
    if (force_flag == 0) {
      MPI_Allreduce(energy,energy_all,nEnergyTerms,MPI_DOUBLE,MPI_SUM,world);
      force_flag = 1;
    }
	
    fprintf(efile, "%d ", ntimestep);
    for (int i=1;i<nEnergyTerms;++i) fprintf(efile, "\t%8.6f", energy_all[i]);
    fprintf(efile, "\t%8.6f\n", energy_all[ET_TOTAL]);
  }
}

/* ---------------------------------------------------------------------- */

void FixBackbone::pre_force(int vflag)
{
  if (amylometer_flag) {
    compute_amylometer();
  }
  else {
    compute_backbone();
  }
}

/* ---------------------------------------------------------------------- */

void FixBackbone::pre_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBackbone::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ----------------------------------------------------------------------
   return total potential energy
   ------------------------------------------------------------------------- */

double FixBackbone::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(energy,energy_all,nEnergyTerms,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return energy_all[ET_TOTAL];
}

/* ----------------------------------------------------------------------
   return potential energies of terms computed in this fix
   ------------------------------------------------------------------------ */

double FixBackbone::compute_vector(int nv)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(energy,energy_all,nEnergyTerms,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return energy_all[nv+1];
}
