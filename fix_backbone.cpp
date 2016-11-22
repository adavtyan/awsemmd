/* ----------------------------------------------------------------------
   Copyright (2010) Aram Davtyan and Garegin Papoian

   Papoian's Group, University of Maryland at Collage Park
   http://papoian.chem.umd.edu/

   Solvent Separated Barrier Potential was contributed by Nick Schafer

   Last Update: 03/23/2011
   ------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "atom_vec_awsemmd.h"
#include "comm.h"
#include "timer.h"
#include <fstream>
#include <time.h>

using std::ifstream;

#define delta 0.00001
#define DEBUGFORCES
#define vfm_small 0.0001

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
	
  char buff[5];
  char forcefile[20]="";
  itoa(comm->me+1,buff,10);
  strcpy(forcefile,"forces\0");
  if (comm->nprocs>1) strcat(forcefile, buff);
  strcat(forcefile, ".dat");
  dout = fopen(forcefile, "w");
	
  char eheader[] = "Step   \tChain   \tShake   \tChi     \tRama    \tExcluded\tDSSP    \tP_AP    \tWater   \tBurial  \tHelix   \tAMH-Go  \tFrag_Mem\tVec_FM  \tSSB     \tVTotal\n";
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
  ssb_flag = frag_mem_tb_flag = vec_frag_mem_tb_flag = phosph_flag = 0;
  epsilon = 1.0; // general energy scale
  p = 2; // for excluded volume
	
  int i, j;
	
  for (i=0;i<12;i++) ssweight[i] = false;
	
  for (i=0;i<TIME_N;i++) ctime[i] = 0.0;

  // backbone geometry coefficients
  an = 0.4831806; bn = 0.7032820; cn = -0.1864262;
  ap = 0.4436538; bp = 0.2352006; cp = 0.3211455;
  ah = 0.8409657; bh = 0.8929599; ch = -0.7338894;
  
  // Default value for fm_sigma_exp
  fm_sigma_exp = 0.15;

  n_wells = 0;
  n_helix_wells = 0;
	
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
      if (comm->me==0) print_log("Vector_Fragment_Memory flag on\n");
      in >> k_vec_frag_mem;
      in >> vfm_sigma;
      vfm_sigma_sq = vfm_sigma*vfm_sigma;
    } else if (strcmp(varsection, "[Vector_Fragment_Memory_Table]")==0) {
      vec_frag_mem_tb_flag = 1;
      if (comm->me==0) print_log("Vector_Fragment_Memory_Table flag on\n");
      in >> k_vec_frag_mem;
      in >> vfm_sigma;
      in >> vfm_tb_size;
      vfm_sigma_sq = vfm_sigma*vfm_sigma;
      vfm_tb_vmin = 0.0;
      vfm_tb_vmax = M_PI;
      vfm_tb_dv = (vfm_tb_vmax - vfm_tb_vmin)/(double)vfm_tb_size;
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
    } else if (strcmp(varsection, "[Phosphorylation]")==0) {
      if (!water_flag) error->all(FLERR,"Cannot run phosphorylation without water potential");
      phosph_flag = 1;
      if (comm->me==0) print_log("Phosphorylation flag on\n");
      in >> k_hypercharge;
      in >> n_phosph_res;
      if (n_phosph_res > 20) error->all(FLERR,"Number of phosphorylated residues may not exceed 20");
      for (int i=0;i<n_phosph_res;++i)
	in >> phosph_res[i];
      
    } else if (strcmp(varsection, "[Epsilon]")==0)
      in >> epsilon;
    varsection[0]='\0'; // Clear buffer
  }
  in.close();
  if (comm->me==0) print_log("\n");
	
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

  // Read sequance file
  ifstream ins(arg[6]);
  if (!ins) error->all(FLERR,"Sequence file was not found");
  char *buf = new char[n+2];
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
  delete [] buf;

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
  
  // Allocate FM the table
  if (frag_mem_tb_flag) {
    if (fm_gamma->maxSep()!=-1)
      tb_nbrs = fm_gamma->maxSep()-fm_gamma->minSep()+1;
    else
      tb_nbrs = n - fm_gamma->minSep();
		
    fm_table = new TBV*[4*n*tb_nbrs];
		
    for (i=0; i<4*n*tb_nbrs; ++i) {
      fm_table[i] = NULL;
    }
	
    if (comm->me==0) print_log("Computing FM table...\n");
    compute_fragment_memory_table();
  }
  
  // Allocate VFM the table
  if (vec_frag_mem_tb_flag) {
  	if (fm_gamma->maxSep()!=-1)
      tb_nbrs = fm_gamma->maxSep()-fm_gamma->minSep()+1;
    else
      tb_nbrs = n - fm_gamma->minSep();
  	
  	vfm_table = new TBV*[n*tb_nbrs];
  	
  	for (i=0; i<n*tb_nbrs; ++i) {
      vfm_table[i] = NULL;
    }
  	
  	if (comm->me==0) print_log("Computing VFM table...\n");
    compute_vector_fragment_memory_table();
  }
  
  sStep=0, eStep=0;
  ifstream in_rs("record_steps");
  in_rs >> sStep >> eStep;
  in_rs.close();
  
  // Debug
//  tmpmax = 0.0;
//  tmpmax2 = 0.0;
}

void FixBackbone::final_log_output()
{
  double time, tmp;
  char txt_timer[][11] = {"Chain", "Shake", "Chi", "Rama", "Vexcluded", "DSSP", "PAP", "Water", "Burial", "Helix", "AHM-Go", "Frag_Mem", "Vec_FM", "SSB"};
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
  
  // Debug
/*  printf("\n\ntmpmax=%f\n", tmpmax);
  printf("iresmax=%d imax=%d jmax=%d steptmp=%d\n\n", iresmax, imax, jmax, steptmp);
  printf("\n\ntmpmax2=%f\n", tmpmax2);
  printf("iresmax2=%d jresmax2=%d steptmp2=%d\n\n", iresmax2, jresmax2, steptmp2);*/
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
    delete [] se;

    delete p_ap;
    delete R;
		
    if (amh_go_flag) {
      for (int i=0;i<3*n;i++) {
	delete [] amh_go_force[i];
      }
      
      delete [] amh_go_force;
      delete [] amh_go_force_map;
			
      delete m_amh_go;
      delete amh_go_gamma;
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
  }
  
  if (vec_frag_mem_tb_flag) {
    for (int i=0; i<n*tb_nbrs; ++i) {
      if (vfm_table[i])
	delete [] vfm_table[i];
    }
    delete [] vfm_table;
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
  res_info = new int[n];
  chain_no = new int[n];
  se = new char[n+2];

  xca = new double*[n];
  xcb = new double*[n];
  xo = new double*[n];
  xn = new double*[n];
  xcp = new double*[n];
  xh = new double*[n];
	
  water_par = WPV(water_kappa, water_kappa_sigma, treshold, n_wells, well_flag, well_r_min, well_r_max);
  helix_par = WPV(helix_kappa, helix_kappa_sigma, helix_treshold, n_helix_wells, helix_well_flag, helix_well_r_min, helix_well_r_max);

  p_ap = new cP_AP<double, FixBackbone>(n, n, &ntimestep, this);
  R = new cR<double, FixBackbone>(n, n, &ntimestep, this);
  well = new cWell<double, FixBackbone>(n, n, n_wells, water_par, &ntimestep, this);
  helix_well = new cWell<double, FixBackbone>(n, n, n_helix_wells, helix_par, &ntimestep, this);

  for (i = 0; i < n; ++i) {
    // Ca, Cb and O coordinates
    xca[i] = new double [3];
    xcb[i] = new double [3];
    xo[i] = new double [3];
		
    // Nitrogen and C prime coordinates
    xn[i] = new double [3];
    xcp[i] = new double [3];
    xh[i] = new double [3];
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
  int *res_tag = avec->residue;

	
  int i;
  for (i=0; i<n; ++i){
  	res_no_l[i] =-1;
  }

//  printf("proc: %d, n: %d\n", comm->me, n);

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
    res_no[nn] = amin;
    res_no_l[res_no[nn]-1]=nn; //local i=res_no_l[i_resno]; 
    last = amin;
    nn++;
  }

//  printf("proc: %d, nn: %d\n", comm->me, nn);

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
	} else if (oxygens[i]==-1 || beta_atoms[i]==-1) {
	  res_info[i] = OFF;
	} else res_info[i] = GHOST;
      }
			
    } else res_info[i] = OFF;
		
    if (i>0 && res_info[i-1]==LOCAL && res_info[i]==OFF) {
      error->all(FLERR,"Missing neighbor atoms in fix backbone (Code 004)");
    }
  }

/*  for (i = 0; i < nn; ++i) {
    printf("proc: %d Ca: %d Cb: %d O: %d Res#: %d Chain# %d ResI: %d\n", comm->me, alpha_carbons[i], beta_atoms[i], oxygens[i], res_no[i], chain_no[i], res_info[i]);
  }*/
	
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
  avec = (AtomVecAWSEM *) atom->style_match("awsemmd");
  if (!avec) error->all(FLERR,"Fix backbone requires atom style awsemmd");

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

void FixBackbone::setup_pre_force(int vflag)
{
  Construct_Computational_Arrays();
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
  char ln[100], *line, *str[10];
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
      mems_array[n_mems-1] = new Fragment_Memory(tpos, fpos, len, weight, str[0], (vec_frag_mem_flag || vec_frag_mem_tb_flag));
      
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
  int im1_resno;
  // N(i) - Cb(i)
  int *res_tag = avec->residue;
  if (!isFirst(i) && se[i_resno]!='G') {
    im1 = res_no_l[i_resno-1];
    im1_resno = res_no[im1]-1;
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
		/*
		if(ntimestep==1){
			fprintf(dout,"%d %d %d %d ", i_resno, res_tag[alpha_carbons[i]]-1, res_tag[oxygens[i]]-1,  res_tag[beta_atoms[i]]-1);
			for(int k=0; k<3; ++k){
				fprintf(dout,"%.9f %.9f ", f[alpha_carbons[i]][k]+bn*dx[k]*force, f[alpha_carbons[i]][k]);
			}
			fprintf(dout,"%.9f %.9f %.9f\n", f[oxygens[i]][0], f[beta_atoms[i]][0]-dx[0]*force, f[beta_atoms[i]][0]);
		}*/
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
  int i_resno = res_no[i]-1 ;
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
    if(im1==-1){
	fprintf(stderr,"im1=-1!\n");
    }
    else {
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
    if(ip1==-1){
	fprintf(stderr,"ip1=-1!\n");
    }
    else {
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
	
  if (iAng==PHI) {
    a[0] = xcp[i][0] - xca[i][0];
    a[1] = xcp[i][1] - xca[i][1];
    a[2] = xcp[i][2] - xca[i][2];
		
    b[0] = xca[i][0] - xn[i][0];
    b[1] = xca[i][1] - xn[i][1];
    b[2] = xca[i][2] - xn[i][2];
		
    c[0] = xn[i][0] - xcp[i-1][0];
    c[1] = xn[i][1] - xcp[i-1][1];
    c[2] = xn[i][2] - xcp[i-1][2];
  } else {
    a[0] = xn[i+1][0] - xcp[i][0];
    a[1] = xn[i+1][1] - xcp[i][1];
    a[2] = xn[i+1][2] - xcp[i][2];
		
    b[0] = xcp[i][0] - xca[i][0];
    b[1] = xcp[i][1] - xca[i][1];
    b[2] = xcp[i][2] - xca[i][2];
		
    c[0] = xca[i][0] - xn[i][0];
    c[1] = xca[i][1] - xn[i][1];
    c[2] = xca[i][2] - xn[i][2];
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
	
  if ( isLast(j) || se[j_resno+1]=='P' ) i_repulsive = false;
  if ( isFirst(i) || isLast(j) || se[i_resno]=='P' ) i_AP = false;
  if ( i_resno==(ch_pos[i_chno]+ch_len[i_chno]-1)-2 || isLast(j) || se[i_resno+2]=='P' ) i_P = false;
	
  /*	if ( j_resno==n-1 || se[j_resno+1]=='P' ) i_repulsive = false;
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
    if (abs(j_resno - i_resno) < 4) {
      lambda[0] = -hbscl[0][0];
      lambda[1] = -hbscl[0][1];
      lambda[2] = 0.0;
      lambda[3] = 0.0;
      hb_class = 1;
    } else if (abs(j_resno - i_resno) < 18) {
      if (aps[n_rama_par-1][i_resno]==1.0 && aps[n_rama_par-1][j_resno]==1.0) {
	lambda[0] = -hbscl[1][0];
	lambda[1] = -hbscl[1][1];
	lambda[2] = -hbscl[1][2]
	  -hbscl[1][3]*theta_seq_anti_HB[0]
	  -hbscl[1][4]*theta_seq_anti_NHB[0]
	  -hbscl[1][5]*(anti_one(se[i_resno])+anti_one(se[j_resno]));
	lambda[3] = -hbscl[1][6];
      } else {
	lambda[0] = 0.0;
	lambda[1] = -hbscl[1][1];
	lambda[2] = 0.0;
	lambda[3] = 0.0;
      }
      hb_class = 2;
    } else if (abs(j_resno - i_resno) < 45) {
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
  } else {
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
    R_NO[1]=R->rNO(i, j+1);
    R_HO[1]=R->rHO(i, j+1);
	
    xNO[1][0] = xo[i][0] - xn[j+1][0];
    xNO[1][1] = xo[i][1] - xn[j+1][1];
    xNO[1][2] = xo[i][2] - xn[j+1][2];
	
    xHO[1][0] = xo[i][0] - xh[j+1][0];
    xHO[1][1] = xo[i][1] - xh[j+1][1];
    xHO[1][2] = xo[i][2] - xh[j+1][2];
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
    R_NO[3]=R->rNO(j, i+2);
    R_HO[3]=R->rHO(j, i+2);
	
    xNO[3][0] = xo[j][0] - xn[i+2][0];
    xNO[3][1] = xo[j][1] - xn[i+2][1];
    xNO[3][2] = xo[j][2] - xn[i+2][2];
	
    xHO[3][0] = xo[j][0] - xh[i+2][0];
    xHO[3][1] = xo[j][1] - xh[i+2][1];
    xHO[3][2] = xo[j][2] - xh[i+2][2];
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

  if (i-2 > 0 && !isFirst(i-1) && !isFirst(i-2) && i+2 < nn && !isLast(i+1) && hb_class!=2) {
    dxnu[0][0] = xca[i+2][0]-xca[i-2][0];
    dxnu[0][1] = xca[i+2][1]-xca[i-2][1];
    dxnu[0][2] = xca[i+2][2]-xca[i-2][2];

    r_nu[0] = sqrt (pow(dxnu[0][0], 2) + pow(dxnu[0][1], 2) + pow(dxnu[0][2], 2) );
		
    th = tanh(pref[0]*(r_nu[0] - d_nu0));
    nu[0] = 0.5*(1+th);

    prdnu[0] = 0.5*pref[0]*(1-pow(th,2))/r_nu[0];
  } else nu[0] = 1.0;

  if (j-2 > 0 && !isFirst(j-1) && !isFirst(j-2) && j+2 < nn && !isLast(j+1) && hb_class!=2) {
    dxnu[1][0] = xca[j+2][0]-xca[j-2][0];
    dxnu[1][1] = xca[j+2][1]-xca[j-2][1];
    dxnu[1][2] = xca[j+2][2]-xca[j-2][2];

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

  if (i-2 > 0 && !isFirst(i-1) && !isFirst(i-2) && i+2 < nn && !isLast(i+1) && hb_class!=2) {
    force = k_dssp*epsilon*theta_sum*prdnu[0]*nu[1];
    f[alpha_carbons[i-2]][0] -= -force*dxnu[0][0];
    f[alpha_carbons[i-2]][1] -= -force*dxnu[0][1];
    f[alpha_carbons[i-2]][2] -= -force*dxnu[0][2];

    f[alpha_carbons[i+2]][0] -= force*dxnu[0][0];
    f[alpha_carbons[i+2]][1] -= force*dxnu[0][1];
    f[alpha_carbons[i+2]][2] -= force*dxnu[0][2];
  }

  if (j-2 > 0 && !isFirst(j-1) && !isFirst(j-2) && j+2 < nn && !isLast(j+1) && hb_class!=2) {
    force = k_dssp*epsilon*theta_sum*nu[0]*prdnu[1];
    f[alpha_carbons[j-2]][0] -= -force*dxnu[1][0];
    f[alpha_carbons[j-2]][1] -= -force*dxnu[1][1];
    f[alpha_carbons[j-2]][2] -= -force*dxnu[1][2];

    f[alpha_carbons[j+2]][0] -= force*dxnu[1][0];
    f[alpha_carbons[j+2]][1] -= force*dxnu[1][1];
    f[alpha_carbons[j+2]][2] -= force*dxnu[1][2];
  }

  f[alpha_carbons[j-1]][0] -= -VTotal*(an*prd_theta[0][0]*xNO[0][0] + ah*prd_theta[0][1]*xHO[0][0]);
  f[alpha_carbons[j-1]][1] -= -VTotal*(an*prd_theta[0][0]*xNO[0][1] + ah*prd_theta[0][1]*xHO[0][1]);
  f[alpha_carbons[j-1]][2] -= -VTotal*(an*prd_theta[0][0]*xNO[0][2] + ah*prd_theta[0][1]*xHO[0][2]);

  f[alpha_carbons[j]][0] -= -VTotal*(bn*prd_theta[0][0]*xNO[0][0] + bh*prd_theta[0][1]*xHO[0][0]);
  f[alpha_carbons[j]][1] -= -VTotal*(bn*prd_theta[0][0]*xNO[0][1] + bh*prd_theta[0][1]*xHO[0][1]);
  f[alpha_carbons[j]][2] -= -VTotal*(bn*prd_theta[0][0]*xNO[0][2] + bh*prd_theta[0][1]*xHO[0][2]);

  f[oxygens[j-1]][0] -= -VTotal*(cn*prd_theta[0][0]*xNO[0][0] + ch*prd_theta[0][1]*xHO[0][0]);
  f[oxygens[j-1]][1] -= -VTotal*(cn*prd_theta[0][0]*xNO[0][1] + ch*prd_theta[0][1]*xHO[0][1]);
  f[oxygens[j-1]][2] -= -VTotal*(cn*prd_theta[0][0]*xNO[0][2] + ch*prd_theta[0][1]*xHO[0][2]);

  f[oxygens[i]][0] -= VTotal*(prd_theta[0][0]*xNO[0][0] + prd_theta[0][1]*xHO[0][0]);
  f[oxygens[i]][1] -= VTotal*(prd_theta[0][0]*xNO[0][1] + prd_theta[0][1]*xHO[0][1]);
  f[oxygens[i]][2] -= VTotal*(prd_theta[0][0]*xNO[0][2] + prd_theta[0][1]*xHO[0][2]);

  if (i_repulsive) {
    f[alpha_carbons[j]][0] -= -V[1]*(an*prd_theta[1][0]*xNO[1][0] + ah*prd_theta[1][1]*xHO[1][0]);
    f[alpha_carbons[j]][1] -= -V[1]*(an*prd_theta[1][0]*xNO[1][1] + ah*prd_theta[1][1]*xHO[1][1]);
    f[alpha_carbons[j]][2] -= -V[1]*(an*prd_theta[1][0]*xNO[1][2] + ah*prd_theta[1][1]*xHO[1][2]);
	
    f[alpha_carbons[j+1]][0] -= -V[1]*(bn*prd_theta[1][0]*xNO[1][0] + bh*prd_theta[1][1]*xHO[1][0]);
    f[alpha_carbons[j+1]][1] -= -V[1]*(bn*prd_theta[1][0]*xNO[1][1] + bh*prd_theta[1][1]*xHO[1][1]);
    f[alpha_carbons[j+1]][2] -= -V[1]*(bn*prd_theta[1][0]*xNO[1][2] + bh*prd_theta[1][1]*xHO[1][2]);
	
    f[oxygens[j]][0] -= -V[1]*(cn*prd_theta[1][0]*xNO[1][0] + ch*prd_theta[1][1]*xHO[1][0]);
    f[oxygens[j]][1] -= -V[1]*(cn*prd_theta[1][0]*xNO[1][1] + ch*prd_theta[1][1]*xHO[1][1]);
    f[oxygens[j]][2] -= -V[1]*(cn*prd_theta[1][0]*xNO[1][2] + ch*prd_theta[1][1]*xHO[1][2]);
	
    f[oxygens[i]][0] -= V[1]*(prd_theta[1][0]*xNO[1][0] + prd_theta[1][1]*xHO[1][0]);
    f[oxygens[i]][1] -= V[1]*(prd_theta[1][0]*xNO[1][1] + prd_theta[1][1]*xHO[1][1]);
    f[oxygens[i]][2] -= V[1]*(prd_theta[1][0]*xNO[1][2] + prd_theta[1][1]*xHO[1][2]);
  }


  if (i_AP) {
    f[alpha_carbons[i-1]][0] -= -V[2]*(an*prd_theta[2][0]*xNO[2][0] + ah*prd_theta[2][1]*xHO[2][0]);
    f[alpha_carbons[i-1]][1] -= -V[2]*(an*prd_theta[2][0]*xNO[2][1] + ah*prd_theta[2][1]*xHO[2][1]);
    f[alpha_carbons[i-1]][2] -= -V[2]*(an*prd_theta[2][0]*xNO[2][2] + ah*prd_theta[2][1]*xHO[2][2]);
	
    f[alpha_carbons[i]][0] -= -V[2]*(bn*prd_theta[2][0]*xNO[2][0] + bh*prd_theta[2][1]*xHO[2][0]);
    f[alpha_carbons[i]][1] -= -V[2]*(bn*prd_theta[2][0]*xNO[2][1] + bh*prd_theta[2][1]*xHO[2][1]);
    f[alpha_carbons[i]][2] -= -V[2]*(bn*prd_theta[2][0]*xNO[2][2] + bh*prd_theta[2][1]*xHO[2][2]);
	
    f[oxygens[i-1]][0] -= -V[2]*(cn*prd_theta[2][0]*xNO[2][0] + ch*prd_theta[2][1]*xHO[2][0]);
    f[oxygens[i-1]][1] -= -V[2]*(cn*prd_theta[2][0]*xNO[2][1] + ch*prd_theta[2][1]*xHO[2][1]);
    f[oxygens[i-1]][2] -= -V[2]*(cn*prd_theta[2][0]*xNO[2][2] + ch*prd_theta[2][1]*xHO[2][2]);
	
    f[oxygens[j]][0] -= V[2]*(prd_theta[2][0]*xNO[2][0] + prd_theta[2][1]*xHO[2][0]);
    f[oxygens[j]][1] -= V[2]*(prd_theta[2][0]*xNO[2][1] + prd_theta[2][1]*xHO[2][1]);
    f[oxygens[j]][2] -= V[2]*(prd_theta[2][0]*xNO[2][2] + prd_theta[2][1]*xHO[2][2]);
  }


  if (i_P) {
    f[alpha_carbons[i+1]][0] -= -V[3]*(an*prd_theta[3][0]*xNO[3][0] + ah*prd_theta[3][1]*xHO[3][0]);
    f[alpha_carbons[i+1]][1] -= -V[3]*(an*prd_theta[3][0]*xNO[3][1] + ah*prd_theta[3][1]*xHO[3][1]);
    f[alpha_carbons[i+1]][2] -= -V[3]*(an*prd_theta[3][0]*xNO[3][2] + ah*prd_theta[3][1]*xHO[3][2]);
	
    f[alpha_carbons[i+2]][0] -= -V[3]*(bn*prd_theta[3][0]*xNO[3][0] + bh*prd_theta[3][1]*xHO[3][0]);
    f[alpha_carbons[i+2]][1] -= -V[3]*(bn*prd_theta[3][0]*xNO[3][1] + bh*prd_theta[3][1]*xHO[3][1]);
    f[alpha_carbons[i+2]][2] -= -V[3]*(bn*prd_theta[3][0]*xNO[3][2] + bh*prd_theta[3][1]*xHO[3][2]);
	
    f[oxygens[i+1]][0] -= -V[3]*(cn*prd_theta[3][0]*xNO[3][0] + ch*prd_theta[3][1]*xHO[3][0]);
    f[oxygens[i+1]][1] -= -V[3]*(cn*prd_theta[3][0]*xNO[3][1] + ch*prd_theta[3][1]*xHO[3][1]);
    f[oxygens[i+1]][2] -= -V[3]*(cn*prd_theta[3][0]*xNO[3][2] + ch*prd_theta[3][1]*xHO[3][2]);
	
    f[oxygens[j]][0] -= V[3]*(prd_theta[3][0]*xNO[3][0] + prd_theta[3][1]*xHO[3][0]);
    f[oxygens[j]][1] -= V[3]*(prd_theta[3][0]*xNO[3][1] + prd_theta[3][1]*xHO[3][1]);
    f[oxygens[j]][2] -= V[3]*(prd_theta[3][0]*xNO[3][2] + prd_theta[3][1]*xHO[3][2]);
  }
}

void FixBackbone::compute_P_AP_potential(int i, int j)
{
  double K, force[2], dx[2][3];
  //	double nu_P_AP[100][100], prd_nu_P_AP[100][100];
  bool i_AP_med, i_AP_long, i_P;

  int i_resno = res_no[i]-1;
  int j_resno = res_no[j]-1;

  // Need to change
  i_AP_med = i_resno<n-(i_med_min+2*i_diff_P_AP) && j_resno>=i_resno+(i_med_min+2*i_diff_P_AP) && j_resno<=MIN(i_resno+i_med_max+2*i_diff_P_AP,n-1);
  i_AP_long = i_resno<n-(i_med_max+2*i_diff_P_AP+1) && j_resno>=i_resno+(i_med_max+2*i_diff_P_AP+1) && j_resno<n;
  i_P = i_resno<n-(i_med_max+1+i_diff_P_AP) && j_resno>=i_resno+(i_med_max+1) && j_resno<n-i_diff_P_AP;

  //	if (ntimestep==0) fprintf(dout, "(%d %d) (%d %d) %d %d %d\n", i, j, i_resno, j_resno, i_AP_med, i_AP_long, i_P);

  if (i_AP_med || i_AP_long) {
    if (aps[n_rama_par-1][i_resno]==1.0 && aps[n_rama_par-1][j_resno]==1.0) {
      K = (i_AP_med ? k_P_AP[0] : 0.0) + (i_AP_long ? k_P_AP[1]*k_betapred_P_AP : 0.0);
    } else {
      K = (i_AP_med ? k_P_AP[0] : 0.0) + (i_AP_long ? k_P_AP[1] : 0.0);
    }

    energy[ET_PAP] += -k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->nu(i+i_diff_P_AP, j-i_diff_P_AP);

    dx[0][0] = xca[i][0] - xca[j][0];
    dx[0][1] = xca[i][1] - xca[j][1];
    dx[0][2] = xca[i][2] - xca[j][2];

    dx[1][0] = xca[i+i_diff_P_AP][0] - xca[j-i_diff_P_AP][0];
    dx[1][1] = xca[i+i_diff_P_AP][1] - xca[j-i_diff_P_AP][1];
    dx[1][2] = xca[i+i_diff_P_AP][2] - xca[j-i_diff_P_AP][2];

    force[0] = k_global_P_AP*epsilon*K*p_ap->prd_nu(i, j)*p_ap->nu(i+i_diff_P_AP, j-i_diff_P_AP);
    force[1] = k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->prd_nu(i+i_diff_P_AP, j-i_diff_P_AP);
	
    f[alpha_carbons[i]][0] -= force[0]*dx[0][0];
    f[alpha_carbons[i]][1] -= force[0]*dx[0][1];
    f[alpha_carbons[i]][2] -= force[0]*dx[0][2];
	
    f[alpha_carbons[j]][0] -= -force[0]*dx[0][0];
    f[alpha_carbons[j]][1] -= -force[0]*dx[0][1];
    f[alpha_carbons[j]][2] -= -force[0]*dx[0][2];

    f[alpha_carbons[i+i_diff_P_AP]][0] -= force[1]*dx[1][0];
    f[alpha_carbons[i+i_diff_P_AP]][1] -= force[1]*dx[1][1];
    f[alpha_carbons[i+i_diff_P_AP]][2] -= force[1]*dx[1][2];

    f[alpha_carbons[j-i_diff_P_AP]][0] -= -force[1]*dx[1][0];
    f[alpha_carbons[j-i_diff_P_AP]][1] -= -force[1]*dx[1][1];
    f[alpha_carbons[j-i_diff_P_AP]][2] -= -force[1]*dx[1][2];
  }
    
  if (i_P) {
    if (aps[n_rama_par-1][i_resno]==1.0 && aps[n_rama_par-1][j_resno]==1.0) {
      K = k_P_AP[2]*k_betapred_P_AP;
    } else {
      K = k_P_AP[2];
    }

    energy[ET_PAP] += -k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->nu(i+i_diff_P_AP, j+i_diff_P_AP);

    dx[0][0] = xca[i][0] - xca[j][0];
    dx[0][1] = xca[i][1] - xca[j][1];
    dx[0][2] = xca[i][2] - xca[j][2];

    dx[1][0] = xca[i+i_diff_P_AP][0] - xca[j+i_diff_P_AP][0];
    dx[1][1] = xca[i+i_diff_P_AP][1] - xca[j+i_diff_P_AP][1];
    dx[1][2] = xca[i+i_diff_P_AP][2] - xca[j+i_diff_P_AP][2];

    force[0] = k_global_P_AP*epsilon*K*p_ap->prd_nu(i, j)*p_ap->nu(i+i_diff_P_AP, j+i_diff_P_AP);
    force[1] = k_global_P_AP*epsilon*K*p_ap->nu(i, j)*p_ap->prd_nu(i+i_diff_P_AP, j+i_diff_P_AP);

    f[alpha_carbons[i]][0] -= force[0]*dx[0][0];
    f[alpha_carbons[i]][1] -= force[0]*dx[0][1];
    f[alpha_carbons[i]][2] -= force[0]*dx[0][2];

    f[alpha_carbons[j]][0] -= -force[0]*dx[0][0];
    f[alpha_carbons[j]][1] -= -force[0]*dx[0][1];
    f[alpha_carbons[j]][2] -= -force[0]*dx[0][2];

    f[alpha_carbons[i+i_diff_P_AP]][0] -= force[1]*dx[1][0];
    f[alpha_carbons[i+i_diff_P_AP]][1] -= force[1]*dx[1][1];
    f[alpha_carbons[i+i_diff_P_AP]][2] -= force[1]*dx[1][2];

    f[alpha_carbons[j+i_diff_P_AP]][0] -= -force[1]*dx[1][0];
    f[alpha_carbons[j+i_diff_P_AP]][1] -= -force[1]*dx[1][1];
    f[alpha_carbons[j+i_diff_P_AP]][2] -= -force[1]*dx[1][2];
  }
}

void FixBackbone::compute_water_potential(int i, int j)
{	
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
  else { xj = xcb[j]; jatom  = beta_atoms[j]; if(jatom==-1)return; }
	
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
	else { katom  = beta_atoms[k]; if(katom==-1)continue; xk = xcb[k]; }
	//else { xk = xcb[k]; katom  = beta_atoms[k]; if(katom==-1)continue;}
				
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

void FixBackbone::compute_helix_potential(int i, int j)
{
  if (R->rNO(i, j)>helix_cutoff) return;
	
  double R_NO, R_HO, xNO[3], xHO[3], dx[3];
  double pair_theta, prd_pair_theta[2], prob_sum;
  double pair_theta_gamma, sigmma_gamma, V;
  double force;
  double *xi, *xj, *xk;
  double hp1, hp2;
  int iatom, jatom, katom, k;
  int k_resno, k_chno;

  int i_resno = res_no[i]-1;
  int j_resno = res_no[j]-1;
	
  int i_chno = chain_no[i]-1;
  int j_chno = chain_no[j]-1;

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
  else { xi = xcb[i]; iatom  = beta_atoms[i]; }
  if (se[j_resno]=='G') { xj = xca[j]; jatom = alpha_carbons[j]; }
  else { xj = xcb[j]; jatom  = beta_atoms[j]; }
	
  dx[0] = xi[0] - xj[0];
  dx[1] = xi[1] - xj[1];
  dx[2] = xi[2] - xj[2];
	
  sigmma_gamma = helix_gamma_p*(1.0-helix_well->sigma(i, j)) + helix_gamma_w*helix_well->sigma(i, j);

  pair_theta_gamma = -epsilon*k_helix*(helix_gamma_w - helix_gamma_p)*pair_theta;
	
  V = -epsilon*k_helix*sigmma_gamma*pair_theta;

  energy[ET_HELIX] += V;
	
  f[alpha_carbons[j-1]][0] -= -V*(an*prd_pair_theta[0]*xNO[0] + ah*prd_pair_theta[1]*xHO[0]);
  f[alpha_carbons[j-1]][1] -= -V*(an*prd_pair_theta[0]*xNO[1] + ah*prd_pair_theta[1]*xHO[1]);
  f[alpha_carbons[j-1]][2] -= -V*(an*prd_pair_theta[0]*xNO[2] + ah*prd_pair_theta[1]*xHO[2]);

  f[alpha_carbons[j]][0] -= -V*(bn*prd_pair_theta[0]*xNO[0] + bh*prd_pair_theta[1]*xHO[0]);
  f[alpha_carbons[j]][1] -= -V*(bn*prd_pair_theta[0]*xNO[1] + bh*prd_pair_theta[1]*xHO[1]);
  f[alpha_carbons[j]][2] -= -V*(bn*prd_pair_theta[0]*xNO[2] + bh*prd_pair_theta[1]*xHO[2]);

  f[oxygens[j-1]][0] -= -V*(cn*prd_pair_theta[0]*xNO[0] + ch*prd_pair_theta[1]*xHO[0]);
  f[oxygens[j-1]][1] -= -V*(cn*prd_pair_theta[0]*xNO[1] + ch*prd_pair_theta[1]*xHO[1]);
  f[oxygens[j-1]][2] -= -V*(cn*prd_pair_theta[0]*xNO[2] + ch*prd_pair_theta[1]*xHO[2]);

  f[oxygens[i]][0] -= V*(prd_pair_theta[0]*xNO[0] + prd_pair_theta[1]*xHO[0]);
  f[oxygens[i]][1] -= V*(prd_pair_theta[0]*xNO[1] + prd_pair_theta[1]*xHO[1]);
  f[oxygens[i]][2] -= V*(prd_pair_theta[0]*xNO[2] + prd_pair_theta[1]*xHO[2]);

  for (k=0;k<nn;++k) {
    if (res_info[k]==OFF) continue;
	
    k_resno = res_no[k]-1;
    k_chno = chain_no[k]-1;
    
    if (se[res_no[k]-1]=='G') { xk = xca[k]; katom = alpha_carbons[k]; }
    else { xk = xcb[k]; katom  = beta_atoms[k]; }
		
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
    ires = avec->residue[i];
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
        jres = avec->residue[j];
        jmol = atom->molecule[j];
        jres_type = se_map[se[jres-1]-'A'];
        
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
	    
//	    if (update->ntimestep==1000 && (i_resno==36 || j_resno==36)) {
//	    	printf("compute_vector_fragment_memory_potential\n");
	//    	printf("i_resno=%d j_resno=%d gc=%f\n", i_resno, j_resno, gc);
	  //  	printf("force=%f forcei={%f %f %f} forcej={%f %f %f}\n\n", force, forcei[0], forcei[1], forcei[2], forcej[0], forcej[1], forcej[2]);
//	    }
//	    if (update->ntimestep==346 && i_resno==33) printf("%f ", -forcei[2]);
//		if (update->ntimestep==346 && j_resno==33) printf("%f ", -forcej[2]);
	    
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
	    
	    //Debug
/*	    tmpforce1[alpha_carbons[i]][0] += -forcei[0];
	    tmpforce1[alpha_carbons[i]][1] += -forcei[1];
	    tmpforce1[alpha_carbons[i]][2] += -forcei[2];
	    
	    tmpforce1[beta_atoms[i]][0] += forcei[0];
	    tmpforce1[beta_atoms[i]][1] += forcei[1];
	    tmpforce1[beta_atoms[i]][2] += forcei[2];
	    
	    tmpforce1[alpha_carbons[j]][0] += -forcej[0];
	    tmpforce1[alpha_carbons[j]][1] += -forcej[1];
	    tmpforce1[alpha_carbons[j]][2] += -forcej[2];
	    
	    tmpforce1[beta_atoms[j]][0] += forcej[0];
	    tmpforce1[beta_atoms[j]][1] += forcej[1];
	    tmpforce1[beta_atoms[j]][2] += forcej[2];*/
	  }
    }
  }
}

void FixBackbone::compute_vector_fragment_memory_table()
{
  int i, j, js, je, i_fm, itb, ig;
  int i_resno, j_resno, ires_type, jres_type;
  double gc, gf, dg;
  double V, epsilon_k_weight, force;
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
		if (je>=n) error->all(FLERR,"Missing residues in memory potential");
		
		for (j=js;j<=je;++j) {
			j_resno = j;
			jres_type = se_map[se[j_resno]-'A'];
		  
//			if (chain_no[i]!=chain_no[j]) error->all(FLERR,"Vector Fragment Memory: Interaction between residues of different chains");
	
			if (se[i_resno]!='G' && se[j_resno]!='G' && frag->getSe(i_resno)!='G' && frag->getSe(j_resno)!='G') {
				gf = frag->VMf(i_resno, j_resno);
				if (frag->error==frag->ERR_CALL || frag->error==frag->ERR_VFM_GLY)
					error->all(FLERR,"Vector_Fragment_Memory: Wrong call of VMf() function");

//				fprintf(dout, "%f ", gf);
		   
				itb = tb_nbrs*i + (j-js);
				if (!vfm_table[itb])
					vfm_table[itb] = new TBV[vfm_tb_size+1];
				
				for (ig=0;ig<=vfm_tb_size;++ig) {
					gc = vfm_tb_vmin + ig*vfm_tb_dv;

					if (1.0-gc<vfm_small) gc = 1.0 - vfm_small;
					if (gc + 1.0<vfm_small) gc = -1.0 + vfm_small;
						
					dg = gc - gf;
							
					V = -epsilon_k_weight*exp(-dg*dg/(2*vfm_sigma_sq));
							
					vfm_table[itb][ig].energy += V;
					
					vfm_table[itb][ig].force += -V*dg/(vfm_sigma_sq*fabs(sin(gc)));

//					fprintf(dout, "%f ", -V*dg/(vfm_sigma_sq*fabs(sin(gc))));
				}
//				fprintf(dout, "\n");
			}
		}
	}
  }
}


void FixBackbone::table_vector_fragment_memory(int i, int j)
{
  int tb_i, tb_j, itb, ig, ig1;
  int i_resno, j_resno, ires_type, jres_type;
  double vi[3], vj[3], vmi, vmj, vmsqi, vmsqj, vp, vpn, gc, g1, g2;
  double v1, v2, f1, f2, ff;
  double V, epsilon_k_weight, forcei[3], forcej[3];
  Fragment_Memory *frag;
  
  i_resno = res_no[i]-1;
  j_resno = res_no[j]-1;
  
  tb_i = i_resno;
  tb_j = j_resno - i_resno - fm_gamma->minSep();

  itb = tb_nbrs*tb_i + tb_j;
  if (!vfm_table[itb]) return;
  
  if ( j_resno-i_resno<fm_gamma->minSep() ) return;
  if ( fm_gamma->maxSep()!=-1 && j_resno-i_resno>fm_gamma->maxSep() ) return;
  
  
  epsilon_k_weight = epsilon*k_vec_frag_mem;
  
//  if (chain_no[i]!=chain_no[j]) error->all(FLERR,"Fragment Memory: Interaction between residues of different chains");

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
  
  if (gc<vfm_tb_vmin) gc=vfm_tb_vmin;
  if (gc>vfm_tb_vmax) gc=vfm_tb_vmax;
  
  ig = int((gc-vfm_tb_vmin)/vfm_tb_dv);
  
  if (ig<0 || ig>vfm_tb_size) error->all(FLERR,"Table Vector Fragment Memory: ig is out of range.");
  
  ig1=ig+1;
  if (ig1>vfm_tb_size) ig1=vfm_tb_size;
  
  // Energy and force values are obtained from trangle interpolation
  g1 = vfm_tb_vmin + (double)ig*vfm_tb_dv;
  g2 = vfm_tb_vmin + (double)(ig+1)*vfm_tb_dv;
  
  v1 = vfm_table[itb][ig].energy;
  v2 = vfm_table[itb][ig1].energy;

  V = ((v2-v1)*gc + v1*g2 - v2*g1)/(g2-g1);
  
  f1 = vfm_table[itb][ig].force;
  f2 = vfm_table[itb][ig1].force;
  
  ff = ((f2-f1)*gc + f1*g2 - f2*g1)/(g2-g1);
  ff /= vmi*vmj;

  energy[ET_VFRAGMEM] += V;
  
  forcei[0] = ff*(vj[0]-vi[0]*vp/vmsqi);
  forcei[1] = ff*(vj[1]-vi[1]*vp/vmsqi);
  forcei[2] = ff*(vj[2]-vi[2]*vp/vmsqi);

  forcej[0] = ff*(vi[0]-vj[0]*vp/vmsqj);
  forcej[1] = ff*(vi[1]-vj[1]*vp/vmsqj);
  forcej[2] = ff*(vi[2]-vj[2]*vp/vmsqj);

   // Debug
/*  if (fabs(forcei[0])>tmpmax2) {
	tmpmax2 = fabs(forcei[0]);
	steptmp2 = update->ntimestep;
	iresmax2 = i;
	jresmax2 = j;
  }
  if (fabs(forcei[1])>tmpmax2) {
	tmpmax2 = fabs(forcei[1]);
	steptmp2 = update->ntimestep;
	iresmax2 = i;
	jresmax2 = j;
  }
  if (fabs(forcei[2])>tmpmax2) {
	tmpmax2 = fabs(forcei[2]);
	steptmp2 = update->ntimestep;
	iresmax2 = i;
	jresmax2 = j;
  }
  if (fabs(forcej[0])>tmpmax2) {
	tmpmax2 = fabs(forcej[0]);
	steptmp2 = update->ntimestep;
	iresmax2 = i;
	jresmax2 = j;
  }
  if (fabs(forcej[1])>tmpmax2) {
	tmpmax2 = fabs(forcej[1]);
	steptmp2 = update->ntimestep;
	iresmax2 = i;
	jresmax2 = j;
  }
  if (fabs(forcej[2])>tmpmax2) {
	tmpmax2 = fabs(forcej[2]);
	steptmp2 = update->ntimestep;
	iresmax2 = i;
	jresmax2 = j;
  }*/
  
/* if (update->ntimestep==0 && (i_resno==1 || j_resno==1)) {
  		printf("table_vector_fragment_memory\n");
		printf("i_resno=%d j_resno=%d gc=%f\n", i_resno, j_resno, gc);
		printf("ig=%d g1=%f g2=%f\n", ig, g1, g2);
		printf("v1=%f v2=%f V=%f\n", v1, v2, V);
		printf("f1=%f f2=%f ff=%f\n", f1, f2, ff);
		printf("forcei={%f %f %f}\n forcej={%f %f %f}\n\n", forcei[0], forcei[1], forcei[2], forcej[0], forcej[1], forcej[2]);
	}
	if (update->ntimestep==346 && i_resno==33) printf("%f ", -forcei[2]);
	if (update->ntimestep==346 && j_resno==33) printf("%f ", -forcej[2]);*/

/*  printf("%d %d\n", i, j);
  printf("%d %d %d %d\n", alpha_carbons[i], beta_atoms[i], alpha_carbons[j], beta_atoms[j]);
  printf("%f %f %f\n", forcei[0], forcei[1], forcei[2]);
  printf("%f %f %f\n", forcej[0], forcej[1], forcej[2]);
  printf("%f %f %f\n", f[alpha_carbons[i]][0], f[alpha_carbons[i]][1], f[alpha_carbons[i]][2]);
  printf("%f %f %f\n", f[beta_atoms[i]][0], f[beta_atoms[i]][1], f[beta_atoms[i]][2]);
  printf("%f %f %f\n", f[alpha_carbons[j]][0], f[alpha_carbons[j]][1], f[alpha_carbons[j]][2]);
  printf("%f %f %f\n", f[beta_atoms[j]][0], f[beta_atoms[j]][1], f[beta_atoms[j]][2]);

  if (i>80 || j>80 || alpha_carbons[i]>242 || beta_atoms[i]>242 || alpha_carbons[j]>242 || beta_atoms[j]>242) {
	printf("\n\nERROR!!!!\n");
	exit(0);
  }*/
//  if (ff<1000.0) {
/*  f[alpha_carbons[i]][0] += -forcei[0];
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
  f[beta_atoms[j]][2] += forcej[2];*/
//  }
  
  //Debug
/*	tmpforce2[alpha_carbons[i]][0] += -forcei[0];
	tmpforce2[alpha_carbons[i]][1] += -forcei[1];
	tmpforce2[alpha_carbons[i]][2] += -forcei[2];
	
	tmpforce2[beta_atoms[i]][0] += forcei[0];
	tmpforce2[beta_atoms[i]][1] += forcei[1];
	tmpforce2[beta_atoms[i]][2] += forcei[2];
	
	tmpforce2[alpha_carbons[j]][0] += -forcej[0];
	tmpforce2[alpha_carbons[j]][1] += -forcej[1];
	tmpforce2[alpha_carbons[j]][2] += -forcej[2];
	
	tmpforce2[beta_atoms[j]][0] += forcej[0];
	tmpforce2[beta_atoms[j]][1] += forcej[1];
	tmpforce2[beta_atoms[j]][2] += forcej[2];*/
}

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
	  }
	}
      }
    }
  }
}

void FixBackbone::table_fragment_memory(int i, int j)
{
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

      v1 = fm_table[itb][ir].energy;
      v2 = fm_table[itb][ir+1].energy;
    	
      V = ((v2-v1)*r + v1*r2 - v2*r1)/(r2-r1);

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
      error->all(FLERR,"Table Fragment Memory: r is out of computed range.");
      fprintf(screen, "r=%f\n", r);
      fprintf(logfile, "r=%f\n", r);
    }	    
  }
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

void FixBackbone::print_forces(int coord)
{
  int index;

  if (coord==1) {
    fprintf(dout, "rca = \n");
    for (int i=0;i<nn;i++) {
      index = alpha_carbons[i];
      if (index!=-1) {
	fprintf(dout, "%.8f %.8f %.8f", x[index][0], x[index][1], x[index][2]);
	if (i!=nn-1) fprintf(dout, "\n");
      }
    }
    fprintf(dout, "\n\n");

    fprintf(dout, "rcb = \n");
    for (int i=0;i<nn;i++) {
      index = beta_atoms[i];
      if (index!=-1) {
	fprintf(dout, "%.8f %.8f %.8f", x[index][0], x[index][1], x[index][2]);
	if (i!=nn-1) fprintf(dout, "\n");
      }
    }
    fprintf(dout, "\n\n");

    fprintf(dout, "ro = \n");
    for (int i=0;i<nn;i++) {
      index = oxygens[i];
      if (index!=-1) {
        int i_resno=res_no[i]-1;
        fprintf(dout, "%d %d %d %.8f %.8f %.8f", comm->me, i_resno, res_info[i], x[index][0], x[index][1], x[index][2]);
	if (i!=nn-1) fprintf(dout, "\n");
      }
    }
    fprintf(dout, "\n\n");

    fprintf(dout, "rn = \n");
    for (int i=0;i<nn;i++) {
      fprintf(dout, "%.8f %.8f %.8f", xn[i][0], xn[i][1], xn[i][2]);
      if (i!=nn-1) fprintf(dout, "\n");
    }
    fprintf(dout, "\n\n");

    fprintf(dout, "rcp = \n");
    for (int i=0;i<nn;i++) {
      fprintf(dout, "%.8f %.8f %.8f", xcp[i][0], xcp[i][1], xcp[i][2]);
      if (i!=nn-1) fprintf(dout, "\n");
    }
    fprintf(dout, "\n\n");
        
    fprintf(dout, "rh = \n");
    for (int i=0;i<nn;i++) {
      fprintf(dout, "%.8f %.8f %.8f", xh[i][0], xh[i][1], xh[i][2]);
      if (i!=nn-1) fprintf(dout, "\n");
    }
    fprintf(dout, "\n\n\n");
  }
	
  fprintf(dout, "fca = \n");
  for (int i=0;i<nn;i++) {
    index = alpha_carbons[i];
    if (index!=-1) {
      int i_resno = res_no[i] - 1;
      fprintf(dout, "%.8f %.8f %.8f", f[index][0], f[index][1], f[index][2]);
      if (i!=nn-1) fprintf(dout, "\n");
    }
  }
  fprintf(dout, "\n\n");

  fprintf(dout, "fcb = \n");
  for (int i=0;i<nn;i++) {
    index = beta_atoms[i];
    if (index!=-1) {
      int i_resno = res_no[i] - 1;
      fprintf(dout, "%.8f %.8f %.8f", f[index][0], f[index][1], f[index][2]);
      if (i!=nn-1) fprintf(dout, "\n");
    }
  }
  fprintf(dout, "\n\n");

  fprintf(dout, "fo = \n");
  for (int i=0;i<nn;i++) {
    index = oxygens[i];
    if (index!=-1) {
      int i_resno = res_no[i] - 1;
      fprintf(dout, "%d %.8f %.8f %.8f", i_resno, f[index][0], f[index][1], f[index][2]);
      if (i!=nn-1) fprintf(dout, "\n");
    }
  }
  fprintf(dout, "\n\n\n\n");
}

void FixBackbone::compute_backbone()
{
  ntimestep = update->ntimestep;

  //if(atom->nlocal==0) return;
    force_flag = 0;
  if(atom->nlocal==0){
        for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;
        for (int i=1;i<nEnergyTerms;++i) energy[ET_TOTAL] += energy[i];
        if (ntimestep%output->thermo_every==0) {
           if (force_flag == 0) {
              MPI_Allreduce(energy,energy_all,nEnergyTerms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
              force_flag = 1;
           }
        }
        return;
  }

  if (comm->nprocs>1 || ntimestep==0)
    Construct_Computational_Arrays();

  x = atom->x;
  f = atom->f;
  image = atom->image;

  int i, j, xbox, ybox, zbox;
  int i_resno, j_resno;
  int i_chno, j_chno;
	
  for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;

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
    if (im1!=-1 && i_resno>0 && !isFirst(i) && (res_info[i]==LOCAL || res_info[i]==GHOST))	{
/*      if (im1==-1){
        printf("proc: %d i: %d i_resno: %d im1: %d isF: %d resI: %d\n", comm->me, i, i_resno, im1, isFirst(i), res_info[i]);
        fprintf(stderr,"Warning: In compute_backbone(), likely the bond was stretched for too long, im1=%d on processor %d, Exit!\n", im1, comm->me);
      	//error->all(FLERR,"In compute_backbone, im1==-1!");
      }	*/
      if (res_info[im1]==LOCAL || res_info[im1]==GHOST){
	      xn[i][0] = an*xca[im1][0] + bn*xca[i][0] + cn*xo[im1][0];
	      xn[i][1] = an*xca[im1][1] + bn*xca[i][1] + cn*xo[im1][1];
	      xn[i][2] = an*xca[im1][2] + bn*xca[i][2] + cn*xo[im1][2];

	      xh[i][0] = ah*xca[im1][0] + bh*xca[i][0] + ch*xo[im1][0];
	      xh[i][1] = ah*xca[im1][1] + bh*xca[i][1] + ch*xo[im1][1];
	      xh[i][2] = ah*xca[im1][2] + bh*xca[i][2] + ch*xo[im1][2];
      }
    } else {
      xn[i][0] = xn[i][1] = xn[i][2] = 0.0;
      xh[i][0] = xh[i][1] = xh[i][2] = 0.0;
    }
		
    if (im1!=-1 && i_resno>0 && !isFirst(i) && (res_info[i]==LOCAL || res_info[i]==GHOST)) {
/*      if (im1==-1){        
        printf("proc: %d i: %d i_resno: %d im1: %d isF: %d resI: %d\n", comm->me, i, i_resno, im1, isFirst(i), res_info[i]);
	fprintf(stderr,"Warning: In compute_backbone(), likely the bond was stretched for too long, im1=%d on processor %d, Exit!\n", im1, comm->me);
      	//error->all(FLERR,"In compute_backbone, im1==-1!");
      }	*/
      if ( (res_info[im1]==LOCAL || res_info[im1]==GHOST) ) {
	xcp[im1][0] = ap*xca[im1][0] + bp*xca[i][0] + cp*xo[im1][0];
	xcp[im1][1] = ap*xca[im1][1] + bp*xca[i][1] + cp*xo[im1][1];
	xcp[im1][2] = ap*xca[im1][2] + bp*xca[i][2] + cp*xo[im1][2];
      } else {
	xcp[im1][0] = xcp[im1][1] = xcp[im1][2] = 0.0;
      }
    }

  }
/*  if(nn<=0){
    	fprintf(stderr, "nn = %d rank = %d\n", nn, comm->me);    	
  	error->all(FLERR,"In compute_backbone, nn <=0 on one processor!");
  }*/
  if (nn>0) xcp[nn-1][0] = xcp[nn-1][1] = xcp[nn-1][2] = 0.0;

 // Debug  
/* for (i=0;i<atom->nlocal;i++) {
 	tmpforce1[i][0] = tmpforce1[i][1] = tmpforce1[i][2] = 0.0;
 	tmpforce2[i][0] = tmpforce2[i][1] = tmpforce2[i][2] = 0.0;
 }*/

#ifdef DEBUGFORCES

  if (ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "AtStart: %d\n", ntimestep);
    fprintf(dout, "Number of residues %d\n", n);
    fprintf(dout, "Local Number of residues %d\n\n", nn);
    print_forces(1);
  }
	
  timerBegin();
	
  for (i=0;i<nn;i++) {
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
      if (!isLast(i) && !isFirst(j) && ( i_chno!=j_chno || abs(j_resno-i_resno)>2 ) && dssp_hdrgn_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) && se[j_resno]!='P')
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
    for (j=0;j<nn;j++) {
      if (i<n-i_med_min && j>=i+i_med_min && p_ap_flag && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
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
      if (vec_frag_mem_tb_flag && j_resno-i_resno>=fm_gamma->minSep() && (fm_gamma->maxSep()==-1 || j_resno-i_resno<=fm_gamma->maxSep()) && chain_no[i]==chain_no[j] && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) )
	table_vector_fragment_memory(i, j);
    }
  }
	
  if (vec_frag_mem_tb_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Table_Vec_Frag_Mem: %d\n", ntimestep);
    fprintf(dout, "Table_Vec_Frag_Mem_Energy: %f\n\n", energy[ET_VFRAGMEM]);
    print_forces();
  }
	
  timerEnd(TIME_VFRAGMEM);
  
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
    if (helix_flag && i_resno<(ch_pos[i_chno]+ch_len[i_chno]-1)-helix_i_diff-1 && i<nn-helix_i_diff && 
	i_chno==chain_no[i+helix_i_diff]-1 && i_resno==res_no[i+helix_i_diff]-helix_i_diff-1 && res_info[i]==LOCAL)
      //		if (helix_flag && i<nn-helix_i_diff-1 && i_resno==res_no[i+helix_i_diff]-helix_i_diff && res_info[i]==LOCAL)
      compute_helix_potential(i, i+helix_i_diff);
  }
	
  if (helix_flag && ntimestep>=sStep && ntimestep<=eStep) {
    fprintf(dout, "Helix: %d\n", ntimestep);
    fprintf(dout, "Helix_Energy: %f\n\n", energy[ET_HELIX]);
    print_forces();
  }
	
  timerEnd(TIME_HELIX);
	
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
      compute_chain_potential(i);

    if (!isFirst(i) && !isLast(i) && chi_flag && res_info[i]==LOCAL && se[i_resno]!='G')
      compute_chi_potential(i);

    if (shake_flag && res_info[i]==LOCAL)
      compute_shake(i);

    if (!isFirst(i) && !isLast(i) && rama_flag && res_info[i]==LOCAL && se[i_resno]!='G')
      compute_rama_potential(i);

    for (j=0;j<nn;j++) {
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      		
      if (dssp_hdrgn_flag && !isLast(i) && !isFirst(j) && ( i_chno!=j_chno || abs(j_resno-i_resno)>2 ) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) && se[j_resno]!='P')
	compute_dssp_hdrgn(i, j);
				
      // Need to change
      if (p_ap_flag && i<n-i_med_min && j>=i+i_med_min && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_P_AP_potential(i, j);

      //if (water_flag && ( i_chno!=j_chno || j_resno-i_resno>=contact_cutoff ) && res_info[i]==LOCAL)
      if (water_flag && ( (i_chno!=j_chno && j_resno > i_resno ) || ( i_chno == j_chno && j_resno-i_resno>=contact_cutoff) ) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_water_potential(i, j);
			  
      if (frag_mem_tb_flag && j_resno-i_resno>=fm_gamma->minSep() && (fm_gamma->maxSep()==-1 || j_resno-i_resno<=fm_gamma->maxSep()) && chain_no[i]==chain_no[j] && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) )
	table_fragment_memory(i, j);

      if (vec_frag_mem_tb_flag && j_resno-i_resno>=fm_gamma->minSep() && (fm_gamma->maxSep()==-1 || j_resno-i_resno<=fm_gamma->maxSep()) && chain_no[i]==chain_no[j] && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) )
	table_vector_fragment_memory(i, j);

      if (ssb_flag && ( i_chno!=j_chno || j_resno-i_resno>=ssb_ij_sep ) && res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST))
	compute_solvent_barrier(i, j);
    }
		
    if (burial_flag && res_info[i]==LOCAL)
      compute_burial_potential(i);
    	
    //    if (helix_flag && i<nn-helix_i_diff-1 && i_resno==res_no[i+helix_i_diff]-helix_i_diff && res_info[i]==LOCAL)
    if (helix_flag && i_resno<(ch_pos[i_chno]+ch_len[i_chno]-1)-helix_i_diff-1 && i<nn-helix_i_diff && 
	i_chno==chain_no[i+helix_i_diff]-1 && i_resno==res_no[i+helix_i_diff]-helix_i_diff-1 && res_info[i]==LOCAL)
      compute_helix_potential(i, i+helix_i_diff);
			
    if (frag_mem_flag && res_info[i]==LOCAL)
      compute_fragment_memory_potential(i);
      
    if (vec_frag_mem_flag && res_info[i]==LOCAL)
      compute_vector_fragment_memory_potential(i);
  }

  if (amh_go_flag)
    compute_amh_go_model();

  if (excluded_flag)
    compute_excluded_volume();

  if (p_excluded_flag)
    compute_p_degree_excluded_volume();

  if (r6_excluded_flag)
    compute_r6_excluded_volume();

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
  
  // Debug
/*  int ii;
  double dtmp;
  for (i=0;i<nn;i++) {
  	for (j=0;j<3;j++) {
  	    ii=alpha_carbons[i];    
  		dtmp=fabs(tmpforce2[ii][j]-tmpforce1[ii][j]);
  		if (dtmp>tmpmax) {
  			tmpmax = dtmp;
  			iresmax=i;
  			imax = 1;
  			jmax = j;
  			steptmp = ntimestep;
  		}
  		
  		ii=beta_atoms[i];    
  		dtmp=fabs(tmpforce2[ii][j]-tmpforce1[ii][j]);
  		if (dtmp>tmpmax) {
  			tmpmax = dtmp;
  			iresmax=i;
  			imax = 2;
  			jmax = j;
  			steptmp = ntimestep;
  		}
  		
  		ii=oxygens[i];    
  		dtmp=fabs(tmpforce2[ii][j]-tmpforce1[ii][j]);
  		if (dtmp>tmpmax) {
  			tmpmax = dtmp;
  			iresmax=i;
  			imax = 3;
  			jmax = j;
  			steptmp = ntimestep;
  		}
  	}
  }*/
  
/*  if (ntimestep==0) {
    printf("tmpforce1 CA\n");
	for (i=0;i<nn;i++) {
		printf("%f %f %f\n", tmpforce1[alpha_carbons[i]][0], tmpforce1[alpha_carbons[i]][1], tmpforce1[alpha_carbons[i]][2]);
	}
	printf("\n");
	printf("tmpforce1 CB\n");
	for (i=0;i<nn;i++) {
		printf("%f %f %f\n", tmpforce1[beta_atoms[i]][0], tmpforce1[beta_atoms[i]][1], tmpforce1[beta_atoms[i]][2]);
	}
	printf("\n\n");
	printf("tmpforce2 CA\n");
	for (i=0;i<nn;i++) {
		printf("%f %f %f\n", tmpforce2[alpha_carbons[i]][0], tmpforce2[alpha_carbons[i]][1], tmpforce2[alpha_carbons[i]][2]);
	}
	printf("\n");
	printf("tmpforce2 CB\n");
	for (i=0;i<nn;i++) {
		printf("%f %f %f\n", tmpforce2[beta_atoms[i]][0], tmpforce2[beta_atoms[i]][1], tmpforce2[beta_atoms[i]][2]);
	}
	printf("\n");
  }*/
}

/* ---------------------------------------------------------------------- */

void FixBackbone::pre_force(int vflag)
{
  compute_backbone();
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
