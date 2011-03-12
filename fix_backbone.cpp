/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 03/04/2011
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_backbone.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "group.h"
#include "domain.h"
#include "fstream.h"

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define delta 0.00001

double max_Edssp = 0.0;
int iEStep = 0;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

// {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
// {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
int se_map[] = {0, 0, 4, 3, 6, 13, 7, 8, 9, 0, 11, 10, 12, 2, 0, 14, 5, 1, 15, 16, 0, 19, 17, 0, 18, 0};


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

FixBackbone::FixBackbone(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
	if (narg != 7) error->all("Illegal fix backbone command");

	scalar_flag = 1;
	vector_flag = 1;
	size_vector = 3;
	extscalar = 1;
	extvector = 1;
	
	abc_flag = chain_flag = shake_flag = chi_flag = rama_flag = rama_p_flag = excluded_flag = p_excluded_flag = r6_excluded_flag = 0;
	ssweight_flag = dssp_hdrgn_flag = p_ap_flag = water_flag = burial_flag = helix_flag = amh_go_flag = frag_mem_flag = 0;
	ssb_flag = 0;
	epsilon = 1.0; // general energy scale
	p = 2; // for excluded volume

	// backbone geometry coefficients
	an = 0.4831806; bn = 0.7032820; cn = -0.1864262;
	ap = 0.4436538; bp = 0.2352006; cp = 0.3211455;
	ah = 0.8409657; bh = 0.8929599; ch = -0.7338894;
	
	igroup2 = group->find(arg[3]);
	if (igroup2 == -1) 
		error->all("Could not find fix backbone beta atoms group ID"); 
	igroup3 = group->find(arg[4]);
	if (igroup3 == -1) 
		error->all("Could not find fix backbone oxygen atoms group ID"); 
	if (igroup2 == igroup || igroup3 == igroup || igroup2 == igroup3) 
		error->all("Two groups cannot be the same in fix backbone"); 
	if (group->count(igroup)!=group->count(igroup2) || group->count(igroup2)!=group->count(igroup3))
		error->all("All groups must contain the same # of atoms in fix backbone");
	group2bit = group->bitmask[igroup2];
	group3bit = group->bitmask[igroup3];
	
	char varsection[30];
	ifstream in(arg[5]);
	if (!in) error->all("Coefficient file was not found!");
	while (!in.eof()) {
		in >> varsection;
    if (strcmp(varsection, "[ABC]")==0) {
      abc_flag = 1;
      fprintf(screen, "ABC flag on\n");
      in >> an >> bn >> cn;
      in >> ap >> bp >> cp;
      in >> ah >> bh >> ch;
		} else if (strcmp(varsection, "[Chain]")==0) {
			chain_flag = 1;
			fprintf(screen, "Chain flag on\n");
			in >> k_chain[0] >> k_chain[1] >> k_chain[2]; 
			in >> r_ncb0 >> r_cpcb0 >> r_ncp0;
		} else if (strcmp(varsection, "[Shake]")==0) {
			shake_flag = 1;
			fprintf(screen, "Shake flag on\n");
			in >> k_shake >> r_sh1 >> r_sh2 >> r_sh3;
		} else if (strcmp(varsection, "[Chi]")==0) {
			chi_flag = 1;
			fprintf(screen, "Chi flag on\n");
			in >> k_chi >> chi0;
		} else if (strcmp(varsection, "[Excluded]")==0) {
			excluded_flag = 1;
			fprintf(screen, "Excluded flag on\n");
			in >> k_excluded_C >> rC_ex0;
			in >> k_excluded_O >> rO_ex0;
		} else if (strcmp(varsection, "[Excluded_P]")==0) {
			p_excluded_flag = 1;
			fprintf(screen, "Excluded_P flag on\n");
			in >> p;
			in >> k_excluded_C >> rC_ex0;
			in >> k_excluded_O >> rO_ex0;
		} else if (strcmp(varsection, "[Excluded_R6]")==0) {
			r6_excluded_flag = 1;
			fprintf(screen, "Excluded_R6 flag on\n");
			in >> k_excluded_C >> rC_ex0;
			in >> k_excluded_O >> rO_ex0;
		} else if (strcmp(varsection, "[Rama]")==0) {
			rama_flag = 1;
			fprintf(screen, "Rama flag on\n");
			in >> k_rama;
			in >> n_rama_par;
			for (int j=0;j<n_rama_par;j++) {
				in >> w[j] >> sigma[j] >> phiw[j] >> phi0[j] >> psiw[j] >> psi0[j];
			}
		} else if (strcmp(varsection, "[Rama_P]")==0) {
			rama_p_flag = 1;
			fprintf(screen, "Rama_P flag on\n");
			in >> n_rama_p_par;
			for (int j=0;j<n_rama_p_par;j++) {
				in >> w[j+i_rp] >> sigma[j+i_rp] >> phiw[j+i_rp] >> phi0[j+i_rp] >> psiw[j+i_rp] >> psi0[j+i_rp];
			}
		} else if (strcmp(varsection, "[SSWeight]")==0) {
			ssweight_flag = 1;
			fprintf(screen, "SSWeight flag on\n");
			for (int j=0;j<12;++j)
				in >> ssweight[j];
		} else if (strcmp(varsection, "[Dssp_Hdrgn]")==0) {
			dssp_hdrgn_flag = 1;
			fprintf(screen, "Dssp_Hdrgn flag on\n");
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
			fprintf(screen, "P_AP flag on\n");
			in >> k_P_AP[0] >> k_P_AP[1] >> k_P_AP[2];
			in >> P_AP_cut;
			in >> P_AP_pref;
			in >> i_med_min >> i_med_max;
			in >> i_diff_P_AP;
		} else if (strcmp(varsection, "[Water]")==0) {
			water_flag = 1;
			fprintf(screen, "Water flag on\n");
			in >> k_water;
			in >> water_kappa >> water_kappa_sigma;
			in >> treshold;
			in >> contact_cutoff;
			in >> n_wells;
			for (int j=0;j<n_wells;++j)
				in >> well_r_min[j] >> well_r_max[j] >> well_flag[j];
		} else if (strcmp(varsection, "[Burial]")==0) {
			burial_flag = 1;
			fprintf(screen, "Burial flag on\n");
			in >> k_burial;
			in >> burial_kappa;
			in >> burial_ro_min[0] >> burial_ro_max[0];
			in >> burial_ro_min[1] >> burial_ro_max[1];
			in >> burial_ro_min[2] >> burial_ro_max[2];
		} else if (strcmp(varsection, "[Helix]")==0) {
			helix_flag = 1;
			fprintf(screen, "Helix flag on\n");
			in >> k_helix;
			in >> helix_gamma_p >> helix_gamma_w;
			in >> helix_kappa >> helix_kappa_sigma;
			in >> helix_treshold;
			in >> helix_i_diff;
			in >> helix_cutoff;
			for (int j=0;j<20;++j)
				in >> h4prob[j];
			in >> helix_sigma_HO >> helix_sigma_NO;
			in >> helix_HO_zero >> helix_NO_zero;
		} else if (strcmp(varsection, "[AMH-Go]")==0) {
      amh_go_flag = 1;
      fprintf(screen, "AMH-Go flag on\n");
      in >> k_amh_go;
      in >> amh_go_p;
      in >> amh_go_rc;
    } else if (strcmp(varsection, "[Fragment_Memory]")==0) {
      frag_mem_flag = 1;
      fprintf(screen, "Fragment_Memory flag on\n");
      in >> k_frag_mem;
      in >> fmem_file;
      in >> fm_gamma_file;
	} else if (strcmp(varsection, "[Solvent_Barrier]")==0) {
	  ssb_flag = 1;
	  fprintf(screen, "Solvent separated barrier flag on\n");
	  in >> k_solventb;
	  in >> ssb_kappa;
	  in >> ssb_rmin0 >> ssb_rmax0;
	  in >> ssb_ij_sep;
	  in >> ssb_rad_cor;
	  for (int j=0;j<20;++j)
		in >> ssb_rshift[j];
		} else if (strcmp(varsection, "[Epsilon]")==0)
			in >> epsilon;
		varsection[0]='\0'; // Clear buffer
	}
	in.close();
	fprintf(screen, "\n");
		
	ifstream ins(arg[6]);
	if (!ins) error->all("Sequence file was not found");
	ins >> se;
	ins.close();
	
	force_flag = 0;
	n = (int)(group->count(igroup)+1e-12);
	foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
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
	
	int i, j;	
	
	if (dssp_hdrgn_flag) {
    ifstream in_anti_HB("anti_HB");
    ifstream in_anti_NHB("anti_NHB");
    ifstream in_para_HB("para_HB");
    ifstream in_para_one("para_one");
    ifstream in_anti_one("anti_one");
    
    if (!in_anti_HB) error->all("File anti_HB doesn't exist");
    if (!in_anti_NHB) error->all("File anti_NHB doesn't exist");
    if (!in_para_HB) error->all("File para_HB doesn't exist");
    if (!in_para_one) error->all("File para_one doesn't exist");
    if (!in_anti_one) error->all("File anti_one doesn't exist");
    
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
		if (!in_ssw) error->all("File ssweight doesn't exist");
		for (j=0;j<n;++j) {
			for (i=0;i<12;++i) {
				if (ssweight[i]) in_ssw >> aps[i][j]; else aps[i][j] = 0.0;
			}
		}
		in_ssw.close();
	}

	if (water_flag) {
		ifstream in_wg("gamma.dat");
		if (!in_wg) error->all("File gamma.dat doesn't exist");
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
	
	if (burial_flag) {
		ifstream in_brg("burial_gamma.dat");
		if (!in_brg) error->all("File burial_gamma.dat doesn't exist");
		for (i=0;i<20;++i) {
			in_brg >> burial_gamma[i][0] >> burial_gamma[i][1] >> burial_gamma[i][2];
		}
		in_brg.close();
	}
	
	if (amh_go_flag) {
	  char amhgo_gamma_file[] = "amh-go.gamma";
		amh_go_gamma = new Gamma_Array(amhgo_gamma_file);
		if (amh_go_gamma->error==amh_go_gamma->ERR_FILE) error->all("Cannot read file amh-go.gamma");
		if (amh_go_gamma->error==amh_go_gamma->ERR_CLASS_DEF) error->all("AMH_Go: Wrong definition of sequance separation classes");
		if (amh_go_gamma->error==amh_go_gamma->ERR_GAMMA) error->all("AMH_Go: Incorrect entery in gamma file");
		if (amh_go_gamma->error==amh_go_gamma->ERR_G_CLASS) error->all("AMH_Go: Wrong sequance separation class tag");
		if (amh_go_gamma->error==amh_go_gamma->ERR_ASSIGN) error->all("AMH_Go: Cannot build gamma array");
    
    char amhgo_mem_file[] = "amh-go.gro";
		m_amh_go = new Fragment_Memory(0, 0, n, amhgo_mem_file);
		if (m_amh_go->error==m_amh_go->ERR_FILE) error->all("Cannot read file amh-go.gro");
		if (m_amh_go->error==m_amh_go->ERR_ATOM_COUNT) error->all("AMH_Go: Wrong atom count in memory file");
		if (m_amh_go->error==m_amh_go->ERR_RES) error->all("AMH_Go: Unknown residue");
	}
	
	if (frag_mem_flag) {
    fm_gamma = new Gamma_Array(fm_gamma_file);
		if (fm_gamma->error==fm_gamma->ERR_FILE) error->all("Fragment_Memory: Cannot read gamma file");
		if (fm_gamma->error==fm_gamma->ERR_CLASS_DEF) error->all("Fragment_Memory: Wrong definition of sequance separation classes");
		if (fm_gamma->error==fm_gamma->ERR_GAMMA) error->all("Fragment_Memory: Incorrect entery in gamma file");
		if (fm_gamma->error==fm_gamma->ERR_G_CLASS) error->all("Fragment_Memory: Wrong sequance separation class tag");
		if (fm_gamma->error==fm_gamma->ERR_ASSIGN) error->all("Fragment_Memory: Cannot build gamma array");
	}
}

/* ---------------------------------------------------------------------- */

FixBackbone::~FixBackbone()
{
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
		delete [] res_info;
		delete [] xca;
		delete [] xcb;
		delete [] xo;
		delete [] xn;
		delete [] xcp;
		delete [] xh;

		delete p_ap;
		delete R;
		
		if (amh_go_flag) {
      int nall = atom->nlocal + atom->nghost;

      for (int i=0;i<nall;i++) {
        delete [] amh_go_force[i];
      }
      
      delete [] amh_go_force;
      delete [] amh_go_force_map;
      
     delete m_amh_go;
     delete amh_go_gamma;
    }    
	}
}

void FixBackbone::allocate()
{
	alpha_carbons = new int[n];
	beta_atoms = new int[n];
	oxygens = new int[n];
	res_no = new int[n];
	res_info = new int[n];
	
	xca = new double*[n];
	xcb = new double*[n];
	xo = new double*[n];
	xn = new double*[n];
	xcp = new double*[n];
	xh = new double*[n];
	
	water_par = WPV(water_kappa, water_kappa_sigma, treshold);
	helix_par = WPV(helix_kappa, helix_kappa_sigma, helix_treshold);

	p_ap = new cP_AP<double, FixBackbone>(n, n, &ntimestep, this);
	R = new cR<double, FixBackbone>(n, n, &ntimestep, this);
	well = new cWell<double, FixBackbone>(n, n, n_wells, water_par, &ntimestep, this);
	helix_well = new cWell<double, FixBackbone>(n, n, n_wells, helix_par, &ntimestep, this);

	for (int i = 0; i < n; ++i) {
		// Ca, Cb and O coordinates
		xca[i] = new double [3];
		xcb[i] = new double [3];
		xo[i] = new double [3];
		
		// Nitrogen and C prime coordinates
		xn[i] = new double [3];
		xcp[i] = new double [3];
		xh[i] = new double [3];
	}

	for (int i = 0; i < 12; ++i) {
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
    int nall = atom->nlocal + atom->nghost;
    
    amh_go_force = new double*[nall];
    amh_go_force_map = new int[nall];
    for (int i=0;i<nall;i++) {
      amh_go_force[i] = new double[3];
    }
	}
	
	allocated = true;
}

/* ---------------------------------------------------------------------- */
inline int MIN(int a, int b)
{
	if ((a<b && a!=-1) || b==-1) return a;
	else return b;
}

inline bool FixBackbone::isFirst(int index)
{
	if (res_no[index]==1) return true;
	if (res_no[index]!=res_no[index-1]+1) return true;
	return false;
}

inline bool FixBackbone::isLast(int index)
{
	if (res_no[index]==n) return true;
	if (res_no[index]!=res_no[index+1]-1) return true;
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

	int i;

	// Creating index arrays for Alpha_Carbons, Beta_Atoms and Oxygens
	nn = 0;
	int last = 0;
	for (i = 0; i < n; ++i) {
		int min[3] = {-1, -1, -1}, jm[3] = {-1, -1, -1}, amin = -1;
		for (int j = 0; j < nall; ++j) {
			if (i==0 && mol_tag[j]<=0)
				error->all("Molecular tag must be positive in fix backbone");
			
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
					error->all("Missing neighbor atoms in fix backbone (Code 001)");
				}
				if ( !isFirst(i) && (i==0 || res_info[i-1]==OFF) ) {
					error->all("Missing neighbor atoms in fix backbone (Code 002)");
				}
				res_info[i] = LOCAL;
			} else {
				if ( i>0 && !isFirst(i) && res_info[i-1]==LOCAL ) res_info[i] = GHOST;
				else if (i<nn-1 && !isLast(i) && alpha_carbons[i+1]<nlocal && alpha_carbons[i+1]!=-1) {
					if (oxygens[i]==-1) {
						error->all("Missing neighbor atoms in fix backbone (Code 003)");
					}
					res_info[i] = GHOST;
				} else res_info[i] = OFF;
			}
			
		} else res_info[i] = OFF;
		
		if (i>0 && res_info[i-1]==LOCAL && res_info[i]==OFF) {
			error->all("Missing neighbor atoms in fix backbone (Code 004)");
		}
	}
}

/* ---------------------------------------------------------------------- */

int FixBackbone::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	mask |= THERMO_ENERGY;
	mask |= POST_FORCE_RESPA;
	mask |= MIN_POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackbone::init()
{
	if (strcmp(update->integrate_style,"respa") == 0)
		nlevels_respa = ((Respa *) update->integrate)->nlevels;
	
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
//  neighbor->requests[irequest]->half = 1;
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
	if (strcmp(update->integrate_style,"verlet") == 0)
		post_force(vflag);
	else {
		((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
		post_force_respa(vflag,nlevels_respa-1,0);
		((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
	}
}

/* ---------------------------------------------------------------------- */

void FixBackbone::min_setup(int vflag)
{
	post_force(vflag);
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

void FixBackbone::compute_chain_potential(int i)
{	
	double dx[3], r, dr, force;

	int i_resno = res_no[i]-1;
	
	// N(i) - Cb(i)
	if (!isFirst(i) && se[i_resno]!='G') {
		dx[0] = xn[i][0] - xcb[i][0];
		dx[1] = xn[i][1] - xcb[i][1];
		dx[2] = xn[i][2] - xcb[i][2];
		r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
		dr = r - r_ncb0;
		force = 2*epsilon*k_chain[0]*dr/r;
	
		foriginal[0] += epsilon*k_chain[0]*dr*dr;
	
		f[alpha_carbons[i-1]][0] -= an*dx[0]*force;
		f[alpha_carbons[i-1]][1] -= an*dx[1]*force;
		f[alpha_carbons[i-1]][2] -= an*dx[2]*force;
		
		f[oxygens[i-1]][0] -= cn*dx[0]*force;
		f[oxygens[i-1]][1] -= cn*dx[1]*force;
		f[oxygens[i-1]][2] -= cn*dx[2]*force;
	
		f[alpha_carbons[i]][0] -= bn*dx[0]*force;
		f[alpha_carbons[i]][1] -= bn*dx[1]*force;
		f[alpha_carbons[i]][2] -= bn*dx[2]*force;
	
		f[beta_atoms[i]][0] -= -dx[0]*force;
		f[beta_atoms[i]][1] -= -dx[1]*force;
		f[beta_atoms[i]][2] -= -dx[2]*force;
	}
	
	
	// Cp(i) - Cb(i)
	if (!isLast(i) && se[i_resno]!='G') {
		dx[0] = xcp[i][0] - xcb[i][0];
		dx[1] = xcp[i][1] - xcb[i][1];
		dx[2] = xcp[i][2] - xcb[i][2];
		r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
		dr = r - r_cpcb0;
		force = 2*epsilon*k_chain[1]*dr/r;
	
		foriginal[0] += epsilon*k_chain[1]*dr*dr;
	
		f[alpha_carbons[i+1]][0] -= bp*dx[0]*force;
		f[alpha_carbons[i+1]][1] -= bp*dx[1]*force;
		f[alpha_carbons[i+1]][2] -= bp*dx[2]*force;
	
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


	// N(i) - Cp(i)
	if (!isFirst(i) && !isLast(i)) {
		dx[0] = xn[i][0] - xcp[i][0];
		dx[1] = xn[i][1] - xcp[i][1];
		dx[2] = xn[i][2] - xcp[i][2];
		r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
		dr = r - r_ncp0;
		force = 2*epsilon*k_chain[2]*dr/r;
	
		foriginal[0] += epsilon*k_chain[2]*dr*dr;
	
		f[alpha_carbons[i-1]][0] -= an*dx[0]*force;
		f[alpha_carbons[i-1]][1] -= an*dx[1]*force;
		f[alpha_carbons[i-1]][2] -= an*dx[2]*force;
			
		f[oxygens[i-1]][0] -= cn*dx[0]*force;
		f[oxygens[i-1]][1] -= cn*dx[1]*force;
		f[oxygens[i-1]][2] -= cn*dx[2]*force;
		
		f[alpha_carbons[i+1]][0] -= -bp*dx[0]*force;
		f[alpha_carbons[i+1]][1] -= -bp*dx[1]*force;
		f[alpha_carbons[i+1]][2] -= -bp*dx[2]*force;
	
		f[alpha_carbons[i]][0] -= (bn-ap)*dx[0]*force;
		f[alpha_carbons[i]][1] -= (bn-ap)*dx[1]*force;
		f[alpha_carbons[i]][2] -= (bn-ap)*dx[2]*force;
	
		f[oxygens[i]][0] -= -cp*dx[0]*force;
		f[oxygens[i]][1] -= -cp*dx[1]*force;
		f[oxygens[i]][2] -= -cp*dx[2]*force;
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
		
		foriginal[0] += epsilon*k_shake*dr*dr;
		
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
	
	foriginal[0] += epsilon*k_shake*dr*dr;
	
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
		
		foriginal[0] += epsilon*k_shake*dr*dr;
		
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
	
	foriginal[0] += epsilon*k_chi*dchi*dchi;

	if (!isFirst(i)) {
		f[alpha_carbons[i-1]][0] -= -an*bprl[0]*force;
		f[alpha_carbons[i-1]][1] -= -an*bprl[1]*force;
		f[alpha_carbons[i-1]][2] -= -an*bprl[2]*force;
		
		f[oxygens[i-1]][0] -= -cn*bprl[0]*force;
		f[oxygens[i-1]][1] -= -cn*bprl[1]*force;
		f[oxygens[i-1]][2] -= -cn*bprl[2]*force;
	}
	
	if (!isLast(i)) {
		f[alpha_carbons[i+1]][0] -= bp*aprl[0]*force;
		f[alpha_carbons[i+1]][1] -= bp*aprl[1]*force;
		f[alpha_carbons[i+1]][2] -= bp*aprl[2]*force;
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
	if (isFirst(i) || isLast(i)) return;

	double V, phi, psi;
	double force, force1[nAngles];
	int jStart, nEnd;
	int j, ia, l;
	
	int i_resno = res_no[i]-1;
	
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
			if (aps[j][i]==0.0) continue;
			V *= aps[j][i];
		}

		force = 2*V*sigma[j];
		force1[PHI] = force*phiw[j]*(cos(phi + phi0[j]) - 1)*sin(phi + phi0[j]);
		force1[PSI] = force*psiw[j]*(cos(psi + psi0[j]) - 1)*sin(psi + psi0[j]);
			
		foriginal[0] += -V;
		
		for (ia=0; ia<nAngles; ia++) {
			for (l=0; l<3; l++) {
				f[alpha_carbons[i-1]][l] += force1[ia]*(y_slope[ia][CA0][l] + x_slope[ia][CA0][l]);
				f[alpha_carbons[i]][l] += force1[ia]*(y_slope[ia][CA1][l] + x_slope[ia][CA1][l]);
				f[alpha_carbons[i+1]][l] += force1[ia]*(y_slope[ia][CA2][l] + x_slope[ia][CA2][l]);
				
				f[oxygens[i-1]][l] += force1[ia]*(y_slope[ia][O0][l] + x_slope[ia][O0][l]);
				f[oxygens[i]][l] += force1[ia]*(y_slope[ia][O1][l] + x_slope[ia][O1][l]);
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
			
				foriginal[0] += epsilon*k_excluded_C*dr*dr;
				
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
			
				foriginal[0] += epsilon*k_excluded_C*dr*dr;
				
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
			
				foriginal[0] += epsilon*k_excluded_C*dr*dr;
				
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
			
				foriginal[0] += epsilon*k_excluded_O*dr*dr;
				
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
			
				foriginal[0] += factorC*epsilon*k_excluded_C*pow(dr, p);
				
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
			
				foriginal[0] += factorC*epsilon*k_excluded_C*pow(dr, p);
				
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
			
				foriginal[0] += factorC*epsilon*k_excluded_C*pow(dr, p);
				
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
			
				foriginal[0] += factorO*epsilon*k_excluded_O*pow(dr, p);
				
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
			
				foriginal[0] += epsilon*k_excluded_C/pow(rsq, 3);
				
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
			
				foriginal[0] += epsilon*k_excluded_C/pow(rsq, 3);
				
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
			
				foriginal[0] += epsilon*k_excluded_C/pow(rsq, 3);
				
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
			
				foriginal[0] += epsilon*k_excluded_O/pow(rsq, 3);
				
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
	
	if ( j_resno==n-1 || se[j_resno+1]=='P' ) i_repulsive = false;
	if ( i_resno==0 || j_resno==n-1 || se[i_resno]=='P' ) i_AP = false;
	if ( i_resno==n-2 || j_resno==n-1 || se[i_resno+2]=='P' ) i_P = false;

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

	V[0] = epsilon*lambda[0]*theta[0]*nu[0]*nu[1];
	V[1] = epsilon*lambda[1]*theta[0]*theta[1]*nu[0]*nu[1];
	V[2] = epsilon*lambda[2]*theta[0]*theta[2]*nu[0]*nu[1];
	V[3] = epsilon*lambda[3]*theta[0]*theta[3]*nu[0]*nu[1];

	VTotal = V[0] + V[1] + V[2] + V[3];

	foriginal[0] +=  epsilon*VTotal;

	if (i-2 > 0 && !isFirst(i-1) && !isFirst(i-2) && i+2 < nn && !isLast(i+1) && hb_class!=2) {
		force = epsilon*theta_sum*prdnu[0]*nu[1];
		f[alpha_carbons[i-2]][0] -= -force*dxnu[0][0];
		f[alpha_carbons[i-2]][1] -= -force*dxnu[0][1];
		f[alpha_carbons[i-2]][2] -= -force*dxnu[0][2];

		f[alpha_carbons[i+2]][0] -= force*dxnu[0][0];
		f[alpha_carbons[i+2]][1] -= force*dxnu[0][1];
		f[alpha_carbons[i+2]][2] -= force*dxnu[0][2];
	}

	if (j-2 > 0 && !isFirst(j-1) && !isFirst(j-2) && j+2 < nn && !isLast(j+1) && hb_class!=2) {
		force = epsilon*theta_sum*nu[0]*prdnu[1];
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

	i_AP_med = i_resno<n-(i_med_min+2*i_diff_P_AP) && j_resno>=i_resno+(i_med_min+2*i_diff_P_AP) && j_resno<=MIN(i_resno+i_med_max+2*i_diff_P_AP,n-1);
	i_AP_long = i_resno<n-(i_med_max+2*i_diff_P_AP+1) && j_resno>=i_resno+(i_med_max+2*i_diff_P_AP+1) && j_resno<n;
	i_P = i_resno<n-(i_med_max+1+i_diff_P_AP) && j_resno>=i_resno+(i_med_max+1) && j_resno<n-i_diff_P_AP;

	if (i_AP_med || i_AP_long) {
		K = (i_AP_med ? k_P_AP[0] : 0.0) + (i_AP_long ? k_P_AP[1] : 0.0);
	
//		foriginal[0] += -K*nu_P_AP[i][j]*nu_P_AP[i+i_diff_P_AP][j-i_diff_P_AP];
		foriginal[0] += -epsilon*K*p_ap->nu(i, j)*p_ap->nu(i+i_diff_P_AP, j-i_diff_P_AP);
	
		dx[0][0] = xca[i][0] - xca[j][0];
		dx[0][1] = xca[i][1] - xca[j][1];
		dx[0][2] = xca[i][2] - xca[j][2];
	
		dx[1][0] = xca[i+i_diff_P_AP][0] - xca[j-i_diff_P_AP][0];
		dx[1][1] = xca[i+i_diff_P_AP][1] - xca[j-i_diff_P_AP][1];
		dx[1][2] = xca[i+i_diff_P_AP][2] - xca[j-i_diff_P_AP][2];
	
//		force[0] = K*prd_nu_P_AP[i][j]*nu_P_AP[i+i_diff_P_AP][j-i_diff_P_AP];
//		force[1] = K*nu_P_AP[i][j]*prd_nu_P_AP[i+i_diff_P_AP][j-i_diff_P_AP];
		force[0] = epsilon*K*p_ap->prd_nu(i, j)*p_ap->nu(i+i_diff_P_AP, j-i_diff_P_AP);
		force[1] = epsilon*K*p_ap->nu(i, j)*p_ap->prd_nu(i+i_diff_P_AP, j-i_diff_P_AP);
	
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
		K = k_P_AP[2];

		foriginal[0] += -epsilon*K*p_ap->nu(i, j)*p_ap->nu(i+i_diff_P_AP, j+i_diff_P_AP);
	
		dx[0][0] = xca[i][0] - xca[j][0];
		dx[0][1] = xca[i][1] - xca[j][1];
		dx[0][2] = xca[i][2] - xca[j][2];
	
		dx[1][0] = xca[i+i_diff_P_AP][0] - xca[j+i_diff_P_AP][0];
		dx[1][1] = xca[i+i_diff_P_AP][1] - xca[j+i_diff_P_AP][1];
		dx[1][2] = xca[i+i_diff_P_AP][2] - xca[j+i_diff_P_AP][2];
	
		force[0] = epsilon*K*p_ap->prd_nu(i, j)*p_ap->nu(i+i_diff_P_AP, j+i_diff_P_AP);
		force[1] = epsilon*K*p_ap->nu(i, j)*p_ap->prd_nu(i+i_diff_P_AP, j+i_diff_P_AP);
	
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
	if (abs(res_no[j]-res_no[i])<contact_cutoff) return;
	
	double dx[3], energy, sigma_gamma, theta_gamma, force;
	double *xi, *xj, *xk;
	int iatom, jatom, katom, i_well, k;
	int k_resno;

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
		
	for (i_well=0;i_well<n_wells;++i_well) {
		if (!well_flag[i_well]) continue;
		if (fabs(well->theta(i, j, i_well))<delta) continue;
		if (water_gamma[i_well][ires_type][jres_type][0] - water_gamma[i_well][ires_type][jres_type][1]<delta) continue;
			
		sigma_gamma = (1.0 - well->sigma(i, j))*water_gamma[i_well][ires_type][jres_type][0] + well->sigma(i, j)*water_gamma[i_well][ires_type][jres_type][1];
		theta_gamma = (water_gamma[i_well][ires_type][jres_type][1] - water_gamma[i_well][ires_type][jres_type][0])*well->theta(i, j, i_well);
				
		foriginal[0] += -epsilon*k_water*sigma_gamma*well->theta(i, j, i_well);
		
		force = epsilon*k_water*sigma_gamma*well->prd_theta(i, j, i_well);
		
		f[iatom][0] += force*dx[0];
		f[iatom][1] += force*dx[1];
		f[iatom][2] += force*dx[2];
		
		f[jatom][0] += -force*dx[0];
		f[jatom][1] += -force*dx[1];
		f[jatom][2] += -force*dx[2];
		
		for (k=0;k<n;++k) {
			if (se[res_no[k]-1]=='G') { xk = xca[k]; katom = alpha_carbons[k]; }
			else { xk = xcb[k]; katom  = beta_atoms[k]; }
			
			k_resno = res_no[k]-1;
			
			if (abs(k_resno-i_resno)>1) {
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
			if (abs(k_resno-j_resno)>1) {
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

void FixBackbone::compute_burial_potential(int i)
{
  double t[3][2], dx[3], force[3], force2, *xi, *xk;
  int iatom, katom, k, k_resno;

  int i_resno = res_no[i]-1;
  
  int ires_type = se_map[se[i_resno]-'A'];
  
  if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
	else { xi = xcb[i]; iatom  = beta_atoms[i]; }
  
  t[0][0] = tanh( burial_kappa*(well->ro(i) - burial_ro_min[0]) );
  t[0][1] = tanh( burial_kappa*(burial_ro_max[0] - well->ro(i)) );
  t[1][0] = tanh( burial_kappa*(well->ro(i) - burial_ro_min[1]) );
  t[1][1] = tanh( burial_kappa*(burial_ro_max[1] - well->ro(i)) );
  t[2][0] = tanh( burial_kappa*(well->ro(i) - burial_ro_min[2]) );
  t[2][1] = tanh( burial_kappa*(burial_ro_max[2] - well->ro(i)) );
  
  foriginal[0] += -0.5*epsilon*k_burial*burial_gamma[ires_type][0]*(t[0][0] + t[0][1]);
  foriginal[0] += -0.5*epsilon*k_burial*burial_gamma[ires_type][1]*(t[1][0] + t[1][1]);
  foriginal[0] += -0.5*epsilon*k_burial*burial_gamma[ires_type][2]*(t[2][0] + t[2][1]);
  
  force[0] = 0.5*epsilon*k_burial*burial_gamma[ires_type][0]*burial_kappa*( t[0][1]*t[0][1] - t[0][0]*t[0][0] );
  force[1] = 0.5*epsilon*k_burial*burial_gamma[ires_type][1]*burial_kappa*( t[1][1]*t[1][1] - t[1][0]*t[1][0] );
  force[2] = 0.5*epsilon*k_burial*burial_gamma[ires_type][2]*burial_kappa*( t[2][1]*t[2][1] - t[2][0]*t[2][0] );
  
  for (k=0;k<n;++k) {
    k_resno = res_no[k]-1;
    
    if (abs(k_resno-i_resno)>1) {
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
	int iatom, jatom, katom, k;
	int k_resno;

	int i_resno = res_no[i]-1;
	int j_resno = res_no[j]-1;

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
	
	prob_sum = h4prob[ires_type] + h4prob[jres_type];
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

	foriginal[0] += V;

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

	for (k=0;k<n;++k) {
    k_resno = res_no[k]-1;
    
		if (se[res_no[k]-1]=='G') { xk = xca[k]; katom = alpha_carbons[k]; }
		else { xk = xcb[k]; katom  = beta_atoms[k]; }
		
		if (abs(k_resno-i_resno)>1) {
			dx[0] = xi[0] - xk[0];
			dx[1] = xi[1] - xk[1];
			dx[2] = xi[2] - xk[2];
			
			force = pair_theta_gamma*helix_well->prd_H(i)*helix_well->H(j)*helix_well->prd_theta(i, k, 0);
			
			f[iatom][0] += force*dx[0];
			f[iatom][1] += force*dx[1];
			f[iatom][2] += force*dx[2];
			
			f[katom][0] += -force*dx[0];
			f[katom][1] += -force*dx[1];
			f[katom][2] += -force*dx[2];
		}
		if (abs(k_resno-j_resno)>1) {
			dx[0] = xj[0] - xk[0];
			dx[1] = xj[1] - xk[1];
			dx[2] = xj[2] - xk[2];

			force = pair_theta_gamma*helix_well->H(i)*helix_well->prd_H(j)*helix_well->prd_theta(j, k, 0);
			
			f[jatom][0] += force*dx[0];
			f[jatom][1] += force*dx[1];
			f[jatom][2] += force*dx[2];
			
			f[katom][0] += -force*dx[0];
			f[katom][1] += -force*dx[1];
			f[katom][2] += -force*dx[2];
		}
	}
}

void FixBackbone::compute_amh_go_model()
{
  int i, j, k, ii, jj, inum, jnum, imol, jmol, iatom, jatom, ires_type, jres_type;
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
    imol = atom->molecule[i];
    ires_type = se_map[se[imol-1]-'A'];
    
    // atom i is either C-Alpha or C-Bata and is LOCAL
    if ( (mask[i]&groupbit || (mask[i]&group2bit && se[imol-1]!='G') ) && i<nlocal ) {
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
        jmol = atom->molecule[j];
        jres_type = se_map[se[jmol-1]-'A'];
        
        // atom j is either C-Alpha or C-Bata
        if ( (mask[j]&groupbit || (mask[j]&group2bit && se[jmol-1]!='G') ) && abs(imol-jmol)>2 ) {
          xj[0] = x[j][0];
          xj[1] = x[j][1];
          xj[2] = x[j][2];
          
          if (domain->xperiodic) xj[0] += prd[0]*((image[j] & 1023) - 512);
          if (domain->yperiodic) xj[1] += prd[1]*((image[j] >> 10 & 1023) - 512);
          if (domain->zperiodic) xj[2] += prd[2]*((image[j] >> 20) - 512);
          
          dx[0] = xi[0] - xj[0];
          dx[1] = xi[1] - xj[1];
          dx[2] = xi[2] - xj[2];
          
          r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
          
          if (r<amh_go_rc) {            
            amhgo_sigma_sq = pow(abs(imol-jmol), 0.3);
//            amhgo_gamma = (abs(imol-jmol)<5 ? k_amh_go[0] : k_amh_go[1]);
            amhgo_gamma = amh_go_gamma->getGamma(ires_type, jres_type, imol-1, jmol-1);
            if (amh_go_gamma->error==amh_go_gamma->ERR_CALL) error->all("AMH-Go: Wrong call of getGamma() function");
            
            if (mask[i]&groupbit) iatom = m_amh_go->FM_CA; else iatom = m_amh_go->FM_CB;
            if (mask[j]&groupbit) jatom = m_amh_go->FM_CA; else jatom = m_amh_go->FM_CB;
            
            rnative = m_amh_go->Rf(imol-1, iatom, jmol-1, jatom);
            dr = r - rnative;
            drsq = dr*dr;
            
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
      
      factor = -0.5*epsilon*k_amh_go*amh_go_p*pow(Ei, amh_go_p-1);
      for (k=0;k<nforces;k++) {
        f[amh_go_force_map[k]][0] += amh_go_force[k][0];
        f[amh_go_force_map[k]][1] += amh_go_force[k][1];
        f[amh_go_force_map[k]][2] += amh_go_force[k][2];
      }
      
      E += -0.5*epsilon*k_amh_go*pow(Ei, amh_go_p);
    }
  }
}

void FixBackbone::compute_solvent_barrier(int i, int j)
{
  if (abs(res_no[j]-res_no[i])<ssb_ij_sep) return;

  double dx[3], force;
  double *xi, *xj, r, rmin, rmax, rshift;
  double t_min, t_max, theta;
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

  r=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

  rmin=ssb_rmin0;
  rmax=ssb_rmax0;
  if(ssb_rad_cor){
	rshift=ssb_rshift[ires_type]+ssb_rshift[jres_type];
	rmin+=rshift;
	rmax+=rshift;
  }

  // apply a distance cutoff criterion, cutoff = rmax + 10/kappa
  if(r>rmax+10/ssb_kappa) return;

  t_min=tanh(ssb_kappa*(r-rmin));
  t_max=tanh(ssb_kappa*(rmax-r));

  theta=0.5*(t_min+t_max);

  foriginal[0] += -epsilon*k_solventb*theta;

  force = epsilon*k_solventb*ssb_kappa*theta*(t_max-t_min)/r;

  f[iatom][0] += force*dx[0];
  f[iatom][1] += force*dx[1];
  f[iatom][2] += force*dx[2];

  f[iatom][0] += -force*dx[0];
  f[iatom][1] += -force*dx[1];
  f[iatom][2] += -force*dx[2];
}

void FixBackbone::compute_backbone()
{
	ntimestep = update->ntimestep;

	if(atom->nlocal==0) return;

	Construct_Computational_Arrays();

	x = atom->x;
	f = atom->f;
	image = atom->image;

	int i, j, xbox, ybox, zbox;
	int i_resno, j_resno;
	
	foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
	force_flag = 0;
	
	for (i=0;i<nn;++i) {
		if (res_info[i]==LOCAL) {
			foriginal[1] += f[alpha_carbons[i]][0] + f[beta_atoms[i]][0] + f[oxygens[i]][0];
			foriginal[2] += f[alpha_carbons[i]][1] + f[beta_atoms[i]][1] + f[oxygens[i]][1];
			foriginal[3] += f[alpha_carbons[i]][2] + f[beta_atoms[i]][2] + f[oxygens[i]][2];
		}
		
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
		
		if ( ( res_info[i]==LOCAL || (res_info[i]==GHOST && i>0 && res_info[i-1]==LOCAL) ) && !isFirst(i) )  {
			xn[i][0] = an*xca[i-1][0] + bn*xca[i][0] + cn*xo[i-1][0];
			xn[i][1] = an*xca[i-1][1] + bn*xca[i][1] + cn*xo[i-1][1];
			xn[i][2] = an*xca[i-1][2] + bn*xca[i][2] + cn*xo[i-1][2];

			xh[i][0] = ah*xca[i-1][0] + bh*xca[i][0] + ch*xo[i-1][0];
			xh[i][1] = ah*xca[i-1][1] + bh*xca[i][1] + ch*xo[i-1][1];
			xh[i][2] = ah*xca[i-1][2] + bh*xca[i][2] + ch*xo[i-1][2];
		} else {
			xn[i][0] = xn[i][1] = xn[i][2] = 0.0;

			xh[i][0] = xh[i][1] = xh[i][2] = 0.0;
		}

		if (i>0) {
			if ( ( res_info[i-1]==LOCAL || (res_info[i-1]==GHOST && res_info[i]==LOCAL ) ) && !isFirst(i) )  {
				xcp[i-1][0] = ap*xca[i-1][0] + bp*xca[i][0] + cp*xo[i-1][0];
				xcp[i-1][1] = ap*xca[i-1][1] + bp*xca[i][1] + cp*xo[i-1][1];
				xcp[i-1][2] = ap*xca[i-1][2] + bp*xca[i][2] + cp*xo[i-1][2];
			} else
				xcp[i-1][0] = xcp[i-1][1] = xcp[i-1][2] = 0.0;
		}

	}
	xcp[nn-1][0] = xcp[nn-1][1] = xcp[nn-1][2] = 0.0;

	for (i=0;i<nn;i++) {
    i_resno = res_no[i];
    j_resno = res_no[j];
    
		if (chain_flag && res_info[i]==LOCAL)
			compute_chain_potential(i);

		if (!isFirst(i) && !isLast(i) && chi_flag && res_info[i]==LOCAL && se[i]!='G')
			compute_chi_potential(i);

		if (shake_flag && res_info[i]==LOCAL)
			compute_shake(i);

		if (!isFirst(i) && !isLast(i) && rama_flag && res_info[i]==LOCAL && se[i]!='G')
			compute_rama_potential(i);

		for (j=0;j<nn;j++) {
			if (!isLast(i) && !isFirst(j) && abs(j_resno-i_resno)>2 && dssp_hdrgn_flag && res_info[i]==LOCAL && res_info[j]==LOCAL && se[j]!='P')
				compute_dssp_hdrgn(i, j);

			if (i<n-i_med_min && j>=i+i_med_min && p_ap_flag && res_info[i]==LOCAL && res_info[j]==LOCAL)
				compute_P_AP_potential(i, j);

			if (water_flag && abs(j_resno-i_resno)>=contact_cutoff && res_info[i]==LOCAL)
			  compute_water_potential(i, j);

			if (ssb_flag && abs(j_resno-i_resno)>=ssb_ij_sep && res_info[i]==LOCAL)
			  compute_solvent_barrier(i, j);
		}
		
 	  if (burial_flag && res_info[i]==LOCAL)
      compute_burial_potential(i);
    
    if (helix_flag && i<nn-helix_i_diff-1 && i_resno==res_no[i+helix_i_diff]-helix_i_diff && res_info[i]==LOCAL)
			compute_helix_potential(i, i+helix_i_diff);
	}

	if (amh_go_flag)
			compute_amh_go_model();

	if (excluded_flag)
		compute_excluded_volume();

	if (p_excluded_flag)
		compute_p_degree_excluded_volume();

	if (r6_excluded_flag)
		compute_r6_excluded_volume();
}

/* ---------------------------------------------------------------------- */

void FixBackbone::post_force(int vflag)
{
	compute_backbone();
}

/* ---------------------------------------------------------------------- */

void FixBackbone::post_force_respa(int vflag, int ilevel, int iloop)
{
	if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBackbone::min_post_force(int vflag)
{
	post_force(vflag);
}

/* ----------------------------------------------------------------------
	 potential energy of added force
------------------------------------------------------------------------- */

double FixBackbone::compute_scalar()
{
	// only sum across procs one time

	if (force_flag == 0) {
		MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
		force_flag = 1;
	}
	return foriginal_all[0];
}

/* ----------------------------------------------------------------------
	 return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixBackbone::compute_vector(int n)
{
	// only sum across procs one time

	if (force_flag == 0) {
		MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
		force_flag = 1;
	}
	return foriginal_all[n+1];
}
