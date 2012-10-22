/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_qbias.h"
#include "atom.h"
#include "timer.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "output.h"
#include "group.h"
#include "domain.h"

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using std::ifstream;

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

inline void FixQBias::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

FixQBias::FixQBias(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
	if (narg != 4) error->all(FLERR,"Illegal fix qbias command");
	
	efile = fopen("energyQ.log", "w");
	
	char eheader[] = "Step\tQBias\tVTotal\n";
	fprintf(efile, "%s", eheader);

	char force_file_name[] = "forcesQ.dat";
	fout = fopen(force_file_name,"w");
	
	scalar_flag = 1;
	vector_flag = 1;
	thermo_energy = 1;
	size_vector = nEnergyTerms;
	global_freq = 1;
	extscalar = 1;
	extvector = 1;

	force_flag = 0;
	n = (int)(group->count(igroup)+1e-12);
	for (int i=0;i<nEnergyTerms;++i) energy[i] = 0.0;

	allocated = false;
	allocate();

	qbias_flag = 0;
	qbias_exp_flag = 0;
	qobias_flag = 0;
	qobias_exp_flag = 0;
	epsilon = 1.0;

	int i, j;
	char varsection[30];
	ifstream in(arg[3]);
	if (!in) error->all(FLERR,"Coefficient file was not found!");
	while (!in.eof()) {
		in >> varsection;
		if (strcmp(varsection, "[Epsilon]")==0) {
			in >> epsilon;
		} else if (strcmp(varsection, "[QBias]")==0) {
			qbias_flag = 1;
			print_log("QBias flag on\n");
			in >> k_qbias;
			in >> q0;
			in >> l;
			in >> sigma;
		}  else if (strcmp(varsection, "[QBias_Exp]")==0) {
			qbias_exp_flag = 1;
			print_log("QBias_Exp flag on\n");
			in >> k_qbias;
			in >> q0;
			in >> l;
			in >> sigma_exp;
		} else if (strcmp(varsection, "[QOBias]")==0) {
			qobias_flag = 1;
			print_log("QOBias flag on\n");
			in >> k_qbias;
			in >> q0;
			in >> l;
			in >> sigma;
			in >> cutoff;
		}  else if (strcmp(varsection, "[QOBias_Exp]")==0) {
			qobias_exp_flag = 1;
			print_log("QOBias_Exp flag on\n");
			in >> k_qbias;
			in >> q0;
			in >> l;
			in >> sigma_exp;
			in >> cutoff;
		}
		varsection[0]='\0';
	}
	in.close();
	print_log("\n");

	ifstream in_rnative("rnative.dat");
	if (!in_rnative) error->all(FLERR,"File rnative.dat can't be read");
	for (i=0;i<n;++i)
		for (j=0;j<n;++j)
			in_rnative >> rN[i][j];
	in_rnative.close();
	
	// Minimal sequence separation
	min_sep = 3;
	if (qobias_flag || qobias_exp_flag) min_sep=4;

	for (i=0;i<n;i++) sigma_sq[i] = Sigma(i)*Sigma(i);

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
}

/* ---------------------------------------------------------------------- */

FixQBias::~FixQBias()
{
	if (allocated) {
		for (int i=0;i<n;i++) {
			delete [] r[i];
			delete [] rN[i];
			delete [] q[i];

			delete [] xca[i];
		}

		delete [] r;
		delete [] rN;
		delete [] q;

		delete [] alpha_carbons;
		delete [] xca;
		delete [] res_no;
		delete [] res_info;
		delete [] chain_no;

		delete [] sigma_sq;
	}
}

/* ---------------------------------------------------------------------- */
inline bool FixQBias::isFirst(int index)
{
	if (res_no[index]==1) return true;
	return false;
}

inline bool FixQBias::isLast(int index)
{
	if (res_no[index]==n) return true;
	return false;
}

int FixQBias::Tag(int index) {
	if (index==-1) return -1;
	return atom->tag[index];
}

inline void FixQBias::Construct_Computational_Arrays()
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
				error->all(FLERR,"Residue index must be positive in fix qbias");
			
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
				error->all(FLERR,"Missing neighbor atoms in fix qbias (code: 001)");
			if (lastResType==GHOST) {
				if (iLastLocal!=-1 && i-nMinNeighbours<=iLastLocal)
					error->all(FLERR,"Missing neighbor atoms in fix qbias (code: 002)");
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
					error->all(FLERR,"Missing neighbor atoms in fix qbias  (code: 003)");
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
		error->all(FLERR,"Missing neighbor atoms in fix qbias  (code: 004)");
	if (lastResType==GHOST) {
		if (iLastLocal!=-1 && nn-nMinNeighbours<=iLastLocal)
			error->all(FLERR,"Missing neighbor atoms in fix qbias  (code: 005)");
	}
}

void FixQBias::allocate()
{
	r = new double*[n];
	rN = new double*[n];
	q = new double*[n];

	sigma_sq = new double[n];

	alpha_carbons = new int[n];
	xca = new double*[n];
	res_no = new int[n];
	res_info = new int[n];
	chain_no = new int[n];

	for (int i = 0; i < n; ++i) {
		r[i] = new double [n];
		rN[i] = new double [n];
		q[i] = new double [n];

		xca[i] = new double [3];
	}
	
	allocated = true;
}

/* ---------------------------------------------------------------------- */

int FixQBias::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	mask |= THERMO_ENERGY;
	mask |= POST_FORCE_RESPA;
	mask |= MIN_POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixQBias::init()
{
	if (strstr(update->integrate_style,"respa"))
		nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixQBias::setup(int vflag)
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

void FixQBias::min_setup(int vflag)
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

inline double FixQBias::PeriodicityCorrection(double d, int i)
{
	if (!periodicity[i]) return d;
	else return ( d > half_prd[i] ? d - prd[i] : (d < -half_prd[i] ? d + prd[i] : d) );
}

double FixQBias::Sigma(int sep)
{
	if (qbias_exp_flag || qobias_exp_flag)
		return pow(1+sep, sigma_exp);

	return sigma;
}

void FixQBias::compute_qbias() 
{
	double qsum, a, dx[3], dr, dql1, dql;
	double force, force1;
	int i, j;

	qsum = 0.0;
	a = 0.0;
//	a = 2.0/((n-2)*(n-3));
	for (i=0;i<n;++i) {
		for (j=i+min_sep;j<n;++j) {
			if ( (qobias_flag || qobias_exp_flag) && rN[i][j]>=cutoff ) continue;
		
			dx[0] = xca[i][0] - xca[j][0];
			dx[1] = xca[i][1] - xca[j][1];
			dx[2] = xca[i][2] - xca[j][2];
			
			r[i][j] = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
			dr = r[i][j] - rN[i][j];
			q[i][j] = exp(-dr*dr/(2*sigma_sq[j-i]));

			qsum += q[i][j];
			a +=1.0;
		}
	}
	qsum = qsum/a;

	dql1 = pow(qsum-q0, l-1);
	dql = dql1*(qsum-q0);
	
	energy[ET_QBIAS] += epsilon*k_qbias*dql;

	force = epsilon*k_qbias*dql1*l/a;

	for (i=0;i<n;++i) {
		for (j=i+min_sep;j<n;++j) {
			if ( (qobias_flag || qobias_exp_flag) && rN[i][j]>=cutoff ) continue;
		
			dx[0] = xca[i][0] - xca[j][0];
			dx[1] = xca[i][1] - xca[j][1];
			dx[2] = xca[i][2] - xca[j][2];

			dr = r[i][j] - rN[i][j];

			force1 = force*q[i][j]*dr/r[i][j]/sigma_sq[j-i];

			f[alpha_carbons[i]][0] -= -dx[0]*force1;
			f[alpha_carbons[i]][1] -= -dx[1]*force1;
			f[alpha_carbons[i]][2] -= -dx[2]*force1;

			f[alpha_carbons[j]][0] -= dx[0]*force1;
			f[alpha_carbons[j]][1] -= dx[1]*force1;
			f[alpha_carbons[j]][2] -= dx[2]*force1;
		}
	}
}

void FixQBias::out_xyz_and_force(int coord)
{
//	out.precision(12);
	
	fprintf(fout, "%d\n", Step);
	fprintf(fout, "%d%d\n", qbias_flag, qbias_exp_flag, qobias_flag, qobias_exp_flag);
	fprintf(fout, "Number of atoms %d\n", n);
	fprintf(fout, "Energy: %d\n\n", energy[ET_TOTAL]);

	int index;	

	if (coord==1) {
		fprintf(fout, "rca = {");
		for (int i=0;i<nn;i++) {
			index = alpha_carbons[i];
			if (index!=-1) {
				fprintf(fout, "{%.8f, %.8f, %.8f}", x[index][0], x[index][1], x[index][2]);
				if (i!=nn-1) fprintf(fout, ",\n");
			}
		}
		fprintf(fout, "};\n\n\n");
	}
	

	fprintf(fout, "fca = {");
	for (int i=0;i<nn;i++) {
		index = alpha_carbons[i];
		if (index!=-1) {
			fprintf(fout, "{%.8f, %.8f, %.8f}", f[index][0], f[index][1], f[index][2]);
			if (i!=nn-1) fprintf(fout, ",\n");

		}
	}
	fprintf(fout, "};\n\n\n\n");
}

void FixQBias::compute()
{
	ntimestep = update->ntimestep;

	Step++;

	if(atom->nlocal==0) return;

	Construct_Computational_Arrays();

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

	double tmp, tmp2;
	double tmp_time;
	int me,nprocs;
  	MPI_Comm_rank(world,&me);
  	MPI_Comm_size(world,&nprocs);
	if (Step>=sStep && Step<=eStep) {
		fprintf(fout, "At Start %d:\n", nn);
		out_xyz_and_force(1);
	}

	tmp = energy[ET_QBIAS];
	if (qbias_flag || qbias_exp_flag || qobias_flag || qobias_exp_flag)
		compute_qbias();
	if ((qbias_flag || qbias_exp_flag || qobias_flag || qobias_exp_flag) && Step>=sStep && Step<=eStep) {
		fprintf(fout, "Qbias %d:\n", nn);
		fprintf(fout, "Qbias_Energy: %.12f\n", energy[ET_QBIAS]-tmp);
		out_xyz_and_force();
	}

	if (Step>=sStep && Step<=eStep) {
		fprintf(fout, "All:\n");
		out_xyz_and_force(1);
		fprintf(fout, "\n\n\n");
	}
	
	for (int i=1;i<nEnergyTerms;++i) energy[ET_TOTAL] += energy[i];
	
	if (ntimestep%output->thermo_every==0) {
		fprintf(efile, "%d ", ntimestep);
		for (int i=1;i<nEnergyTerms;++i) fprintf(efile, "\t%.6f", energy[i]);
			fprintf(efile, "\t%.6f\n", energy[ET_TOTAL]);
	}
}

/* ---------------------------------------------------------------------- */

void FixQBias::post_force(int vflag)
{
	compute();
}

/* ---------------------------------------------------------------------- */

void FixQBias::post_force_respa(int vflag, int ilevel, int iloop)
{
	if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQBias::min_post_force(int vflag)
{
	post_force(vflag);
}

/* ----------------------------------------------------------------------
	 return total potential energy
------------------------------------------------------------------------- */

double FixQBias::compute_scalar()
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
------------------------------------------------------------------------- */

double FixQBias::compute_vector(int nv)
{
	// only sum across procs one time

	if (force_flag == 0) {
		MPI_Allreduce(energy,energy_all,nEnergyTerms,MPI_DOUBLE,MPI_SUM,world);
		force_flag = 1;
	}
	return energy_all[nv+1];
}
