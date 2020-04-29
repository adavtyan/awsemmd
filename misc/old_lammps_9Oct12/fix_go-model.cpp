/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Gaussian Contacts Potential was contributed by Weihua Zheng

Last Update: 03/23/2011
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_go-model.h"
#include "atom.h"
#include "timer.h"
#include "output.h"
#include "update.h"
#include "respa.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "domain.h"
#include "random_park.h"

#include <fstream>
#include <time.h>

using std::ifstream;

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

// {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
// {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
//int se_map[] = {0, 0, 4, 3, 6, 13, 7, 8, 9, 0, 11, 10, 12, 2, 0, 14, 5, 1, 15, 16, 0, 19, 17, 0, 18, 0};

inline void FixGoModel::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

FixGoModel::FixGoModel(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
	if (narg != 4) error->all(FLERR,"Illegal fix go-model command");
	
	efile = fopen("energyGO.log", "w");
	
	char eheader[] = "Step\tBond\tAngle\tDihedral\tContacts\tNative\tVTotal\n";
	fprintf(efile, "%s", eheader);
	
	restart_global = 1;

	char force_file_name[] = "forcesGO.dat";
	fout = fopen(force_file_name,"w");
	efout = fopen("epsilon.dat","w");
	
	scalar_flag = 1;
	vector_flag = 1;
	thermo_energy = 1;
	size_vector = nEnergyTerms;
	global_freq = 1;
	extscalar = 1;
	extvector = 1;

	force_flag = 0;
	for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;
	n = (int)(group->count(igroup)+1e-12);
	
	seed = time(NULL);
//	seed = 1;

	allocated = false;
	allocate();

	bonds_flag = angles_flag = dihedrals_flag = contacts_flag = contacts_dev_flag = lj_contacts_flag = gaussian_contacts_flag = 0;
	dev_type = DT_NONE;
	epsilon = epsilon2 = 1.0;
	n_basins = 1;
	rmin_cutoff = 4.0 ;
	int i, j, k;
	double r0_tmp, theta0_tmp, phi0_tmp;

	char varsection[100];
	ifstream in(arg[3]);
	if (!in) error->all(FLERR,"Coefficient file was not found!");
	while (!in.eof()) {
		in >> varsection;
		if (strcmp(varsection, "[Go-Model_LJ]")==0) {
			if(gaussian_contacts_flag || bonds_flag || angles_flag || dihedrals_flag || contacts_flag ) error->all(FLERR,"Conflict in definition of contact potential !!");
			lj_contacts_flag = 1;
			print_log("LJ Go-Model flag on\n");
			allocate_contact();
			in >> epsilon >> epsilon2 ;
		} else if (strcmp(varsection, "[Go-Model_Gaussian]")==0) {
			if(lj_contacts_flag || bonds_flag || angles_flag || dihedrals_flag || contacts_flag ) error->all(FLERR,"Conflict in definition of contact potential !!");
			gaussian_contacts_flag = 1;
			print_log("Gaussian Go-Model flag on\n");
			in >> epsilon >> epsilon2 ;
			in >> n_basins;
			allocate_contact();
			for(i=0; i<n_basins; ++i){
				in >> A[i];
			}
			in >> gaussian_width;
			g_w_sq_inv = 1.0/gaussian_width/gaussian_width ;
			in >> rmin_cutoff ;
		} else if (strcmp(varsection, "[Bonds]")==0) {
			bonds_flag = 1;
			print_log("Bonds flag on\n");
			in >> k_bonds;
			for(i=0;i<n_basins; ++i)	{
				for (j=0;j<n-1;++j) {
					if(i==0) r0[j]=0.0;
					in >> r0_tmp;
					r0[j] += r0_tmp;
					if(i==n_basins -1) r0[j] /= n_basins;
				}
			}
		} else if (strcmp(varsection, "[Angles]")==0) {
			angles_flag = 1;
			print_log("Angles flag on\n");
			in >> k_angles;
			for(i=0;i<n_basins; ++i)	{
				for (j=0;j<n-2;++j) {
					if(i==0) theta0[j]=0.0;
					in >> theta0_tmp;
					theta0[j] += theta0_tmp;
					if(i==n_basins -1) theta0[j] /= n_basins;
				}
			}
		} else if (strcmp(varsection, "[Dihedrals]")==0) {
			dihedrals_flag = 1;
			print_log("Dihedrals flag on\n");
			in >> k_dihedrals[0] >> k_dihedrals[1];
			for(i=0;i<n_basins; ++i)	{
				for (j=0;j<n-3;++j) {
					if(i==0) phi0[j]=0.0;
					in >> phi0_tmp;
					phi0[j] += phi0_tmp;
					if(i==n_basins -1) phi0[j] /= n_basins;
				}
			}
		} else if (strcmp(varsection, "[Contacts]")==0) {
			if (!lj_contacts_flag && !gaussian_contacts_flag) error->all(FLERR,"Conflict in definition of contact potential !!");
		
			//allocate_contact();
			contacts_flag = 1;
			print_log("Contacts flag on\n");
			if(lj_contacts_flag){
				for (i=0;i<n-4;++i) for (j=i;j<n-4;++j) in >> isNative[i][j];
				for (i=0;i<n-4;++i) for (j=i;j<n-4;++j) in >> sigma[i][j];
				
				for (i=0;i<n-4;++i) for (j=i;j<n-4;++j) sigma_sq[i][j] = sigma[i][j]*sigma[i][j];
			} else { //gaussian_contacts_flag
				for(k=0;k<n_basins;++k){
					for (i=0;i<n-4;++i) for (j=i;j<n-4;++j) in >> isNative_mb[k][i][j];
					for (i=0;i<n-4;++i) for (j=i;j<n-4;++j) in >> sigma_mb[k][i][j];
				}
			}
		} else if (strcmp(varsection, "[Contacts_Deviation]")==0) {
			contacts_dev_flag = 1;
			dev_type = DT_CORR;
			print_log("Contacts_Deviation flag on\n");
			in >> sdivf; // Standart deviation in epsilon fractions
			in >> tcorr; // Correlation time in femtoseconds
			in >> dev0;  // Deviation on t=0
		} else if (strcmp(varsection, "[Harmonic_Contacts_Deviation]")==0) {
			contacts_dev_flag = 1;
			dev_type = DT_SIN;
			print_log("Harmonic_Contacts_Deviation flag on\n");
			in >> sdivf; // Amplitud
			in >> tcorr; // Period
			in >> dev0;  // Phase on t=0 in half periods
		} else if (strcmp(varsection, "[Constant_Contacts_Deviation]")==0) {
			contacts_dev_flag = 1;
			dev_type = DT_CONST;
			print_log("Constant_Contacts_Deviation flag on\n");
			in >> sdivf; // Deviation in epsilon fractions
		}
		varsection[0]='\0';
	}
	in.close();
	print_log("\n");

	if (dev_type==DT_CORR) {
		xi = 1/tcorr;
		w = sqrt(xi)*epsilon*sdivf;
		dev = dev0;

		devA = w*sqrt(2*update->dt);
		devB = xi*update->dt;
		devC = (1 - xi*update->dt/2);
	} else if (dev_type==DT_SIN) {	
		devA = epsilon*sdivf/sqrt(0.5);
		devB = 2*M_PI*update->dt/tcorr;
		devC = M_PI*dev0;
		
		dev = devA*sin(devC);
	} else if (dev_type==DT_CONST) {
		dev = devA = epsilon*sdivf;
	}
	
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

	Step = 0;
	sStep=0, eStep=0;
	ifstream in_rs("record_steps");
	in_rs >> sStep >> eStep;
	in_rs.close();
	
	Construct_Computational_Arrays();
}

/* ---------------------------------------------------------------------- */

FixGoModel::~FixGoModel()
{
	if (allocated) {
		for (int i=0;i<n;i++) {
			delete [] xca[i];
		}

		delete [] r0;
		delete [] theta0;
		delete [] phi0;

		delete [] alpha_carbons;
		delete [] xca;
		delete [] res_no;
		delete [] chain_no;
		delete [] res_info;

		delete random;
	}
	
	if (contacts_allocated) {
		if(lj_contacts_flag){
			for (int i=0;i<n-4;i++) {
				delete [] sigma[i];
				delete [] sigma_sq[i];
				delete [] isNative[i];
			}
			delete [] sigma;
			delete [] sigma_sq;
			delete [] isNative;
		} else { //gaussian_contacts_flag
			delete [] G;
			delete [] A;
			for(int k=0; k<n_basins; ++k){
				for(int i=0; i<n-4; ++i){
					delete [] sigma_mb[k][i];
					delete [] isNative_mb[k][i];
				}
				delete [] sigma_mb[k];
				delete [] isNative_mb[k];
			}
			delete [] sigma_mb;
			delete [] isNative_mb;
		}
	}
}

/* ---------------------------------------------------------------------- */
inline bool FixGoModel::isFirst(int index)
{
	if (res_no[index]==1) return true;
	return false;
}

inline bool FixGoModel::isLast(int index)
{
	if (res_no[index]==n) return true;
	return false;
}

int FixGoModel::Tag(int index) {
	if (index==-1) return -1;
	return atom->tag[index];
}

inline void FixGoModel::Construct_Computational_Arrays()
{
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	int nall = atom->nlocal + atom->nghost;
	int *mol_tag = atom->molecule;
	int *res_tag = atom->residue;

	int i, j, js;

	// Creating index arrays for Alpha_Carbons
	nn = 0;
	int last = 0;
	for (i = 0; i < n; ++i) {
		int min = -1, jm = -1;
		for (int j = 0; j < nall; ++j) {
			if (i==0 && res_tag[j]<=0)
				error->all(FLERR,"Residue index must be positive in fix go-model");
			
			if ( (mask[j] & groupbit) && res_tag[j]>last ) {
				if (res_tag[j]<min || min==-1) {
					min = res_tag[j];
					jm = j;
				}
			}
		}
		
		if (min==-1) break;

		alpha_carbons[nn] = jm;
		res_no[nn] = min;
		last = min;
		nn++;
	}

	int nMinNeighbours = 3;
	int iLastLocal = -1;
	int lastResNo = -1;
	int lastResType = NONE;
	int nlastType = 0;

	// Checking sequance and marking residues
	for (i = 0; i < nn; ++i) {
		chain_no[i] = -1;
	
		if (lastResNo!=-1 && lastResNo!=res_no[i]-1) {
			if (lastResType==LOCAL && res_no[i]!=n)
				error->all(FLERR,"Missing neighbor atoms in fix go-model (code: 001)");
			if (lastResType==GHOST) {
				if (iLastLocal!=-1 && i-nMinNeighbours<=iLastLocal)
					error->all(FLERR,"Missing neighbor atoms in fix go-model (code: 002)");
			}
			
			iLastLocal = -1;
			lastResNo = -1;
			lastResType = NONE;
			nlastType = 0;
		}

		if (alpha_carbons[i]!=-1) {
			chain_no[i] = mol_tag[alpha_carbons[i]];
			
			if (alpha_carbons[i]<nlocal) {
				if ( lastResType==OFF || (lastResType==GHOST && nlastType<nMinNeighbours && nlastType!=res_no[i]-1) ) {
					error->all(FLERR,"Missing neighbor atoms in fix go-model  (code: 003)");
				}
				iLastLocal = i;
				res_info[i] = LOCAL;
			} else {
				res_info[i] = GHOST;
			}
		} else res_info[i] = OFF;

		if (lastResNo == res_no[i]) nlastType++; else nlastType = 0;

		lastResNo = res_no[i];
		lastResType = res_info[i];
	}
	if (lastResType==LOCAL && res_no[nn-1]!=n)
		error->all(FLERR,"Missing neighbor atoms in fix go-model  (code: 004)");
	if (lastResType==GHOST) {
		if (iLastLocal!=-1 && nn-nMinNeighbours<=iLastLocal)
			error->all(FLERR,"Missing neighbor atoms in fix go-model  (code: 005)");
	}
}

void FixGoModel::allocate()
{
	alpha_carbons = new int[n];
	xca = new double*[n];
	res_no = new int[n];
	chain_no = new int[n];
	res_info = new int[n];

	r0 = new double[n-1];
	theta0 = new double[n-2];
	phi0 = new double[n-3];

	for (int i = 0; i < n; ++i) {
		xca[i] = new double [3];
	}

	random = new RanPark(lmp,seed);
	
	allocated = true;
}

void FixGoModel::allocate_contact()
{
	if(lj_contacts_flag){	
		sigma = new double*[n-4];
		sigma_sq = new double*[n-4];
		isNative = new bool*[n-4];
		for (int i = 0; i < n-4; ++i) {
			sigma[i] = new double[n-4];
			sigma_sq[i] = new double[n-4];
			isNative[i] = new bool[n-4];
		}
	} else { //gaussian_contacts_flag
		G = new double[n_basins];
		A = new double[n_basins];
		sigma_mb = new double**[n_basins];
		isNative_mb = new bool**[n_basins];
		for (int k=0; k<n_basins; ++k){
			sigma_mb[k]    = new double*[n-4];
			isNative_mb[k] = new bool*[n-4];
			for(int i=0; i< n-4; ++i){
				sigma_mb[k][i]    = new double[n-4];
				isNative_mb[k][i] = new bool[n-4];
			}
		}
	}
	contacts_allocated = true;
}

/* ---------------------------------------------------------------------- */

int FixGoModel::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	mask |= THERMO_ENERGY;
	mask |= POST_FORCE_RESPA;
	mask |= MIN_POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixGoModel::init()
{
	if (strstr(update->integrate_style,"respa"))
		nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixGoModel::setup(int vflag)
{
	if (strstr(update->integrate_style,"verlet"))
		post_force(vflag);
	else {
		((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
		post_force_respa(vflag,nlevels_respa-1,0);
		((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
	}
}

/* ---------------------------------------------------------------------- */

void FixGoModel::min_setup(int vflag)
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

inline double FixGoModel::PeriodicityCorrection(double d, int i)
{
	if (!periodicity[i]) return d;
	else return ( d > half_prd[i] ? d - prd[i] : (d < -half_prd[i] ? d + prd[i] : d) );
}

void FixGoModel::compute_bond(int i) 
{
	int i_resno = res_no[i]-1;

	if (i_resno>=n-1) error->all(FLERR,"Wrong use of compute_bond() in fix go-model");

	double dx[3], r, dr, force;

	dx[0] = xca[i+1][0] - xca[i][0];
	dx[1] = xca[i+1][1] - xca[i][1];
	dx[2] = xca[i+1][2] - xca[i][2];
	r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
	dr = r - r0[i_resno];
	force = 2*epsilon*k_bonds*dr/r;
	
	energy[ET_BOND] += epsilon*k_bonds*dr*dr;

	f[alpha_carbons[i+1]][0] -= dx[0]*force;
	f[alpha_carbons[i+1]][1] -= dx[1]*force;
	f[alpha_carbons[i+1]][2] -= dx[2]*force;

	f[alpha_carbons[i]][0] -= -dx[0]*force;
	f[alpha_carbons[i]][1] -= -dx[1]*force;
	f[alpha_carbons[i]][2] -= -dx[2]*force;
}

void FixGoModel::compute_angle(int i) 
{
	int i_resno = res_no[i]-1;

	if (i_resno>=n-2) error->all(FLERR,"Wrong use of compute_angle() in fix go-model");

	double a[3], b[3], a2, b2, A, B, AB, adb, alpha;
	double theta, dtheta, force, factor[2][3];

	a[0] = xca[i][0] - xca[i+1][0];
	a[1] = xca[i][1] - xca[i+1][1];
	a[2] = xca[i][2] - xca[i+1][2];

	b[0] = xca[i+2][0] - xca[i+1][0];
	b[1] = xca[i+2][1] - xca[i+1][1];
	b[2] = xca[i+2][2] - xca[i+1][2];

	a2 = adotb(a, a);
	b2 = adotb(b, b);
	adb = adotb(a, b);
	A = sqrt(a2);
	B = sqrt(b2);
	AB = A*B;
	alpha = adb/AB;

	theta = acos(alpha);
	dtheta = theta - theta0[i_resno];

	force = 2*epsilon*k_angles*dtheta/(AB*sqrt(1-alpha*alpha));

	energy[ET_ANGLE] += epsilon*k_angles*dtheta*dtheta;

	factor[0][0] = (b[0] - a[0]*adb/a2);
	factor[0][1] = (b[1] - a[1]*adb/a2);
	factor[0][2] = (b[2] - a[2]*adb/a2);

	factor[1][0] = (a[0] - b[0]*adb/b2);
	factor[1][1] = (a[1] - b[1]*adb/b2);
	factor[1][2] = (a[2] - b[2]*adb/b2);

	f[alpha_carbons[i]][0] -= -factor[0][0]*force;
	f[alpha_carbons[i]][1] -= -factor[0][1]*force;
	f[alpha_carbons[i]][2] -= -factor[0][2]*force;

	f[alpha_carbons[i+1]][0] -= (factor[1][0] + factor[0][0])*force;
	f[alpha_carbons[i+1]][1] -= (factor[1][1] + factor[0][1])*force;
	f[alpha_carbons[i+1]][2] -= (factor[1][2] + factor[0][2])*force;

	f[alpha_carbons[i+2]][0] -= -factor[1][0]*force;
	f[alpha_carbons[i+2]][1] -= -factor[1][1]*force;
	f[alpha_carbons[i+2]][2] -= -factor[1][2]*force;
}

void FixGoModel::compute_dihedral(int i)
{
	int i_resno = res_no[i]-1;

	if (i_resno>=n-3) error->all(FLERR,"Wrong use of compute_dihedral() in fix go-model");

	double phi, dphi, y_slope[4][3], x_slope[4][3], force, V;
	double a[3], b[3], c[3];
	double bxa[3], cxa[3], cxb[3];
	double adb, bdc, adc, b2, bm, cdbxa;
	double X, Y, X2Y2;
	double dAngle_y, dAngle_x;
	double h1, h2, h3;
	int j, l;

	a[0] = xca[i+1][0] - xca[i][0];
	a[1] = xca[i+1][1] - xca[i][1];
	a[2] = xca[i+1][2] - xca[i][2];
	
	b[0] = xca[i+2][0] - xca[i+1][0];
	b[1] = xca[i+2][1] - xca[i+1][1];
	b[2] = xca[i+2][2] - xca[i+1][2];

	c[0] = xca[i+3][0] - xca[i+2][0];
	c[1] = xca[i+3][1] - xca[i+2][1];
	c[2] = xca[i+3][2] - xca[i+2][2];

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
	
	Y = -bm*cdbxa;
	X = adb*bdc - b2*adc;
	
	phi = atan2(Y, X); 

	X2Y2 = (X*X + Y*Y);
	dAngle_y = X/X2Y2;
	dAngle_x = -Y/X2Y2;

	for (l=0;l<3;++l) {
		y_slope[0][l] = dAngle_y*bm*cxb[l];
		y_slope[1][l] = dAngle_y*( b[l]*cdbxa/bm - bm*(cxb[l] + cxa[l]) );
		y_slope[2][l] = dAngle_y*( -b[l]*cdbxa/bm + bm*(bxa[l] + cxa[l]) );
		y_slope[3][l] = -dAngle_y*bm*bxa[l];
	
		h1 = b[l]*bdc - c[l]*b2;
		h2 = a[l]*bdc - 2*b[l]*adc + c[l]*adb;
		h3 = b[l]*adb - a[l]*b2;
		x_slope[0][l] = -dAngle_x*h1;
		x_slope[1][l] = dAngle_x*(h1 - h2);
		x_slope[2][l] = dAngle_x*(h2 - h3);
		x_slope[3][l] = dAngle_x*h3;
	}	
	for (j=0;j<2;j++) { 
		dphi = (2*j+1)*(phi - phi0[i_resno]);
		V = epsilon*k_dihedrals[j]*(1-cos(dphi));

		force = (2*j+1)*epsilon*k_dihedrals[j]*sin(dphi);

		energy[ET_DIHEDRAL] += V;

		for (l=0; l<3; l++) {
			f[alpha_carbons[i]][l] -= force*(y_slope[0][l] + x_slope[0][l]);
			f[alpha_carbons[i+1]][l] -= force*(y_slope[1][l] + x_slope[1][l]);
			f[alpha_carbons[i+2]][l] -= force*(y_slope[2][l] + x_slope[2][l]);
			f[alpha_carbons[i+3]][l] -= force*(y_slope[3][l] + x_slope[3][l]);
		}
	}
	
}


/*  Correlated noise generator for contact potential
    dev(t+dt) = dev(t) - 0.5*xi*(dev(t) + devp(t))*dt + w*sqrt(2*dt)*rand
    devp(t) = dev(t) - xi*dev(t)*dt + w*sqrt(2*dt)*rand

    The formula above was simplified to
    dev += (w*sqrt(2*dt)*rand - xi*dt*dev)*(1 - xi*dt/2)

    devA = w*sqrt(2*dt)
    devB = xi*dt
    devC = (1 - xi*dt/2) 

    <dev(t)dev(t+dt)> = w^2/xi * Exp[-xi*dt]
*/
void FixGoModel::compute_contact_deviation()
{
  if (dev_type==DT_CORR) {
    if (ntimestep!=0) {
		rand = random->gaussian();
		dev += (devA*rand - devB*dev)*devC;
	}
  } else if (dev_type==DT_SIN) {
    dev = devA*sin(devB*ntimestep + devC);
  } else if (dev_type==DT_CONST) {
    dev = devA;
  }
}

void FixGoModel::compute_contact(int i, int j)
{
	int i_resno = res_no[i]-1;
	int j_resno = res_no[j]-1;

	if (i_resno>=n-4 || j_resno<=i_resno+3) error->all(FLERR,"Wrong use of compute_contact() in fix go-model");

	double dx[3], rsq, sgrinvsq, sgrinv12, sgrinv10;
	double V, force, contact_epsilon;

	dx[0] = xca[j][0] - xca[i][0];
	dx[1] = xca[j][1] - xca[i][1];
	dx[2] = xca[j][2] - xca[i][2];

	rsq = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
//	r = sqrt(rsq);

	if (rsq>sigma_sq[i_resno][j_resno-4]*9) return;

	sgrinvsq = sigma_sq[i_resno][j_resno-4]/rsq;
	sgrinv12 = pow(sgrinvsq, 6);
	sgrinv10 = pow(sgrinvsq, 5);

//	sgrinv = sigma[i_resno][j_resno-4]/r;
//	sgrinv12 = pow(sgrinv, 12);
//	sgrinv10 = pow(sgrinv, 10);

	if (isNative[i_resno][j_resno-4]) {
		contact_epsilon = epsilon;
		if (contacts_dev_flag) contact_epsilon = epsilon + dev;

		V = contact_epsilon*(5*sgrinv12 - 6*sgrinv10);
		force = -60*contact_epsilon*(sgrinv12 - sgrinv10)/rsq;
		
		energy[ET_NCONTS] += V;
	} else {
		V = epsilon2*sgrinv12;
		force = -12*epsilon2*sgrinv12/rsq;
	}

	energy[ET_CONTACTS] += V;
	
	f[alpha_carbons[j]][0] -= force*dx[0];
	f[alpha_carbons[j]][1] -= force*dx[1];
	f[alpha_carbons[j]][2] -= force*dx[2];

	f[alpha_carbons[i]][0] -= -force*dx[0];
	f[alpha_carbons[i]][1] -= -force*dx[1];
	f[alpha_carbons[i]][2] -= -force*dx[2];
}

void FixGoModel::compute_contact_gaussian(int i, int j)
{
	int i_resno = res_no[i]-1;
	int j_resno = res_no[j]-1;
	int k, l;

	if (i_resno>=n-4 || j_resno<=i_resno+3) error->all(FLERR,"Wrong use of compute_contact_gaussian() in fix go-model");

	double dx[3], r, dr, rsq;
	double V, force, contact_epsilon;
	double w_sq_inv, VTotal ;
	double force_tmp, force_tmp_sum;

	dx[0] = xca[j][0] - xca[i][0];
	dx[1] = xca[j][1] - xca[i][1];
	dx[2] = xca[j][2] - xca[i][2];
	rsq = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
	r = sqrt(rsq);
	w_sq_inv = 1.0/pow( abs(i_resno - j_resno), 0.3 ) * g_w_sq_inv ; // / gaussian_width / gaussian_width ;

	//V == epsilon * {Product_k(1+V_k[k]) - 1}
	R = pow( (rmin_cutoff/r) , 12); //rmin_cutoff = 4.0 A
	
	V = 1.0;
	for(k=0; k<n_basins; ++k){
		if(isNative_mb[k][i_resno][j_resno-4]){
			dr = r - sigma_mb[k][i_resno][j_resno-4] ;
			G[k]= - A[k] * exp(-dr*dr*w_sq_inv/2.0);
			V *= 1.0 + G[k];
			//fprintf(screen, "%f %f %f\n", xca[i][0],xca[i][1], xca[i][2]);
			//fprintf(screen, "%f %f %f\n", xca[j][0],xca[j][1], xca[j][2]);
			//fprintf(screen, "i %d j %d sigma %f\n", i_resno, j_resno, sigma_mb[k][i_resno][j_resno-4]);
			//fprintf(screen, "r %f dr %f A %f G %e w_sq_inv %e\n", r, dr, A[k], G[k], w_sq_inv);
		} else {
			G[k] = 0.0;
		}
	}
	VTotal = epsilon*((1.0 + R)*V - 1.0);

	energy[ET_CONTACTS] += VTotal;

	force = 12.0 * R * V / r;
	//fprintf(screen, "R %f V %f f %f\n", R, V, force);
	force_tmp_sum = 0.0 ;
	for(k=0; k<n_basins; ++k){
		if(isNative_mb[k][i_resno][j_resno-4]){
			dr = r - sigma_mb[k][i_resno][j_resno-4] ;
			force_tmp = G[k] * dr;
			for(l=0; l<n_basins; ++l){
				if(l != k ) {
					force_tmp *= 1.0 + G[l];
				}
			}
			force_tmp_sum += force_tmp ;
		}
	}
	force_tmp_sum *= 1.0 + R ;
	force += force_tmp_sum ;
	//fprintf(screen, "f_tmp_sum %f force %f eps %f r %f\n", force_tmp_sum, force, epsilon, r);	
	force *= - epsilon / r ;
	//fprintf(screen, "f_tmp_sum %f force %f\n", force_tmp_sum, force);	
	
	f[alpha_carbons[j]][0] -= force*dx[0];
	//fprintf(screen, "f_tmp_sum %f force_final %f\n", force_tmp_sum, force*dx[0]);	
	f[alpha_carbons[j]][1] -= force*dx[1];
	f[alpha_carbons[j]][2] -= force*dx[2];

	f[alpha_carbons[i]][0] -= -force*dx[0];
	f[alpha_carbons[i]][1] -= -force*dx[1];
	f[alpha_carbons[i]][2] -= -force*dx[2];
}

void FixGoModel::out_xyz_and_force(int coord)
{
//	out.precision(12);
	
	fprintf(fout, "%d\n", Step);
	fprintf(fout, "%d%d%d%d\n", bonds_flag, angles_flag, dihedrals_flag, contacts_flag);
	fprintf(fout, "Number of atoms %d\n", n);

	int index;	

	if (coord==1) {
		fprintf(fout, "rca = {");
		for (int i=0;i<nn;i++) {
			index = alpha_carbons[i];
			if (index!=-1) {
				fprintf(fout, "{%.12f, %.12f, %.12f}", x[index][0], x[index][1], x[index][2]);
				if (i!=nn-1) fprintf(fout, ",\n");
			}
		}
		fprintf(fout, "};\n\n\n");
	}
	

	fprintf(fout, "fca = {");
	for (int i=0;i<nn;i++) {
		index = alpha_carbons[i];
		if (index!=-1) {
			fprintf(fout, "{%.12f, %.12f, %.12f}", f[index][0], f[index][1], f[index][2]);
			if (i!=nn-1) fprintf(fout, ",\n");

		}
	}
	fprintf(fout, "};\n\n\n\n");
}

void FixGoModel::compute_goModel()
{
	ntimestep = update->ntimestep;

	Step++;

	if(atom->nlocal==0) return;
	
//	int me,nprocs;
//  	MPI_Comm_rank(world,&me);
//  	MPI_Comm_size(world,&nprocs);

//	Construct_Computational_Arrays();

	x = atom->x;
        f = atom->f;
        image = atom->image;

	int i, j, xbox, ybox, zbox;
	int i_resno, j_resno;
	
	for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;
	force_flag = 0;

	for (i=0;i<nn;++i) {

		// Calculating xca Ca atoms coordinates array
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
		}
	}
	
//	seed = random->state();
	
/*	if (contacts_dev_flag)
		compute_contact_deviation();

	for (i=0;i<nn;++i) {
		if (res_info[i]!=LOCAL) continue;

		if (bonds_flag && res_no[i]<=n-1)
			compute_bond(i);

		if (angles_flag && res_no[i]<=n-2)
    	    compute_angle(i);

		if (dihedrals_flag && res_no[i]<=n-3)
	        compute_dihedral(i);

		for (j=i+1;contacts_flag && j<nn;++j) {
			if (res_info[j]!=LOCAL && res_info[j]!=GHOST) continue;
			
			if (contacts_flag) {
				if (lj_contacts_flag && res_no[i]<res_no[j]-3)
					compute_contact(i, j);
					
				if (gaussian_contacts_flag && res_no[i]<res_no[j]-3)
					compute_contact_gaussian(i, j);
			}
        	}
	}*/
	
/*	double tmp, tmp2;
	double tmp_time;
	int me,nprocs;
  	MPI_Comm_rank(world,&me);
  	MPI_Comm_size(world,&nprocs);
	if (Step>=sStep && Step<=eStep) {
		fprintf(fout, "At Start %d:\n", nn);
		out_xyz_and_force(1);
	}
	
	if (contacts_dev_flag)
		compute_contact_deviation();

	tmp = energy[ET_BOND];
	for (i=0;i<nn;i++) {
		if (bonds_flag && res_info[i]==LOCAL && res_no[i]<=n-1)
			compute_bond(i);
	}

	if (bonds_flag && Step>=sStep && Step<=eStep) {
		fprintf(fout, "Bonds %d:\n", nn);
		fprintf(fout, "Bonds_Energy: %.12f\n", energy[ET_BOND]-tmp);
		out_xyz_and_force();
	}

	tmp = energy[ET_ANGLE];
	for (i=0;i<nn;i++) {
		if (angles_flag && res_info[i]==LOCAL && res_no[i]<=n-2)
			compute_angle(i);
	}

	if (angles_flag && Step>=sStep && Step<=eStep) {
		fprintf(fout, "Angles %d:\n", nn);
		fprintf(fout, "Angles_Energy: %.12f\n", energy[ET_ANGLE]-tmp);
		out_xyz_and_force();
	}

	tmp = energy[ET_DIHEDRAL];
	for (i=0;i<nn;i++) {
		if (dihedrals_flag && res_info[i]==LOCAL && res_no[i]<=n-3)
			compute_dihedral(i);
	}

	if (dihedrals_flag && Step>=sStep && Step<=eStep) {
		fprintf(fout, "Dihedrals %d:\n", nn);
		fprintf(fout, "Dihedrals_Energy: %.12f\n", energy[ET_DIHEDRAL]-tmp);
		out_xyz_and_force();
	}
	
//	fprintf(efout, "%.15f\n", epsilon + dev);

	tmp = energy[ET_CONTACTS];
//	if (contacts_flag && Step>=sStep && Step<=eStep) fprintf(fout, "\n{");
	for (i=0;i<nn;i++) {
		for (j=0;j<nn;j++) {
			tmp2 = energy[ET_CONTACTS];
			if (res_info[i]==LOCAL && (res_info[j]==LOCAL || res_info[j]==GHOST) && res_no[i]<res_no[j]-3){
				if (lj_contacts_flag)
					compute_contact(i, j);
				
				if (gaussian_contacts_flag)
					compute_contact_gaussian(i, j);
			}
//			if (contacts_flag && Step>=sStep && Step<=eStep)
//				fprintf(fout, "%f, ", energy[ET_CONTACTS]-tmp2);
		}
	}
//	if (contacts_flag && Step>=sStep && Step<=eStep) fprintf(fout, "}\n");

	if (contacts_flag && Step>=sStep && Step<=eStep) {
		fprintf(fout, "Contacts %d:\n", nn);
		fprintf(fout, "Contacts_Energy: %.12f\n", energy[ET_CONTACTS]-tmp);
		out_xyz_and_force();
	}

	if (Step>=sStep && Step<=eStep) {
		fprintf(fout, "All:\n");
		out_xyz_and_force(1);
		fprintf(fout, "\n\n\n");
	}
	
	*/
		
	
	if (contacts_dev_flag)
		compute_contact_deviation();
		
//	fprintf(efout, "%f\n", dev);
	
	for (i=0;i<nn;++i) {
		if (res_info[i]!=LOCAL) continue;
		
		if (bonds_flag && res_no[i]<=n-1)
			compute_bond(i);

		if (angles_flag && res_no[i]<=n-2)
    	    compute_angle(i);

		if (dihedrals_flag && res_no[i]<=n-3)
	        compute_dihedral(i);

	    for (j=i+1;contacts_flag && j<nn;++j) {
			if (res_info[j]!=LOCAL && res_info[j]!=GHOST) continue;
			
			if (res_no[i]<res_no[j]-3){
				if (lj_contacts_flag)
					compute_contact(i, j);
				
				if (gaussian_contacts_flag)
					compute_contact_gaussian(i, j);
			}
		}
	}
	
	// Last term is the energy of native contacts, already included in contacts energy
	for (int i=1;i<nEnergyTerms-1;++i) energy[ET_TOTAL] += energy[i];
	
	if (ntimestep%output->thermo_every==0) {
		fprintf(efile, "%d ", ntimestep);
		for (int i=1;i<nEnergyTerms;++i) fprintf(efile, "\t%.6f", energy[i]);
			fprintf(efile, "\t%.6f\n", energy[ET_TOTAL]);
	}
}

/* ---------------------------------------------------------------------- */

void FixGoModel::post_force(int vflag)
{
	compute_goModel();
}

/* ---------------------------------------------------------------------- */

void FixGoModel::post_force_respa(int vflag, int ilevel, int iloop)
{
	if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGoModel::min_post_force(int vflag)
{
	post_force(vflag);
}

/* ----------------------------------------------------------------------
	 return total potential energy
------------------------------------------------------------------------- */

double FixGoModel::compute_scalar()
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
	 and return contact_epsilon
------------------------------------------------------------------------- */

double FixGoModel::compute_vector(int nv)
{
	// contact_epsilon output
	if (nv==nEnergyTerms-1) {
		if (contacts_dev_flag) return epsilon + dev;
		else return epsilon;
	}
	
	// Energy output
	
	// only sum across procs one time
	if (force_flag == 0) {
		MPI_Allreduce(energy,energy_all,nEnergyTerms,MPI_DOUBLE,MPI_SUM,world);
		force_flag = 1;
	}
	return energy_all[nv+1];
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */


void FixGoModel::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];

  list[n++] = dev_type;
  list[n++] = random->state();
//  list[n++] = seed;
  if (dev_type==DT_CORR) {
    list[n++] = sdivf;
    list[n++] = tcorr;
    list[n++] = dev0;
    list[n++] = dev;
  } else if (dev_type==DT_SIN) {
    list[n++] = sdivf;
    list[n++] = tcorr;
    list[n++] = devB*ntimestep + devC;
    list[n++] = dev;
  } else if (dev_type==DT_CONST) {
    list[n++] = sdivf;
    list[n++] = dev;
  }

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixGoModel::restart(char *buf)
{	
  int n = 0;
  double *list = (double *) buf;
  int r_dev_type = DT_NONE;
  double r_sdivf, r_tcorr, r_dev0, r_dev;

  r_dev_type = static_cast<int> (list[n++]);
  seed = static_cast<int> (list[n++]);
  
  if (r_dev_type==DT_CORR || r_dev_type==DT_SIN) {
    r_sdivf = static_cast<double> (list[n++]);
    r_tcorr = static_cast<double> (list[n++]);
    r_dev0 = static_cast<double> (list[n++]);
    r_dev = static_cast<double> (list[n++]);
  } else if (r_dev_type==DT_CONST) {
    r_sdivf = static_cast<double> (list[n++]);
    r_dev = static_cast<double> (list[n++]);
  }
  
  if (contacts_dev_flag && r_dev_type==dev_type) {
  	if (r_dev_type==DT_CORR) {
    	dev = r_dev;
  	} else if (r_dev_type==DT_SIN) {
  		devC = r_dev0;
    	dev = devA*sin(devC);
  	} else if (r_dev_type==DT_CONST) {
  		// Nothing need to be done
  	}
  	
  	random->reset(seed);
  }
}
