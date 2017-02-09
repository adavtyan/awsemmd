/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pair_excluded_volume.h"
#include "atom.h"
#include "atom_vec_awsemmd.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairExcludedVolume::PairExcludedVolume(LAMMPS *lmp) : Pair(lmp)
{
  PI = 4.0*atan(1.0);
}

/* ---------------------------------------------------------------------- */

PairExcludedVolume::~PairExcludedVolume()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(lambda);
    memory->destroy(prefactor);
    memory->destroy(cut_short);
    memory->destroy(cut_long);
  }
}

/* ---------------------------------------------------------------------- */

void PairExcludedVolume::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,sign,imol,jmol,ires,jres;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rsq,dr,factor_lj,pref,rcut;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int ntypes = atom->ntypes;
  sign = (p%2==0 ? 1 : -1);
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      prefactor[i][j] = lambda[i][j] * sign;
    }
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      imol = atom->molecule[i];
      jmol = atom->molecule[j];
      
      ires = avec->residue[i];
      jres = avec->residue[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);
      jtype = type[j];
      
      if (abs(ires-jres)<5 && imol==jmol) rcut = cut_short[itype][jtype];
      else rcut = cut_long[itype][jtype];

      if (r < rcut) {
	dr = r - rcut;
	pref = factor_lj * prefactor[itype][jtype] * pow(dr, p-1);
	if (r > 0.0) fpair = - pref * p / r;
	else fpair = 0.0;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}

	if (eflag)
	  evdwl = pref * dr;

	if (evflag) ev_tally(i,j,nlocal,newton_pair,
			     evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairExcludedVolume::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1,n+1,"pair:cutsq");

   memory->create(lambda, n+1,n+1,"pair:lambda");
   memory->create(prefactor, n+1,n+1,"pair:prefactor");
   memory->create(cut_short, n+1,n+1,"pair:cut_short");
   memory->create(cut_long, n+1,n+1,"pair:cut_long");

  // init lambda to 0.0
  // since pair_hybrid can use all types even if pair_excluded_volume sub-class
  //   never sets them

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) {
      lambda[i][j] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairExcludedVolume::settings(int narg, char **arg)
{ 
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  p = atoi(arg[0]);
  cut_global[0] = atof(arg[1]);
  cut_global[1] = atof(arg[2]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) {
          cut_short[i][j] = cut_global[0];
          cut_long[i][j] = cut_global[1];
	}
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairExcludedVolume::coeff(int narg, char **arg)
{
  if (narg != 3 && narg != 5) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  double cut_one[2];

  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double lambda_one = atof(arg[2]);

  cut_one[0] = cut_global[0];
  cut_one[1] = cut_global[1];
  if (narg == 5) {
    cut_one[0] = atof(arg[3]);
    cut_one[1] = atof(arg[4]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      lambda[i][j] = lambda_one;
      cut_short[i][j] = cut_one[0];
      cut_long[i][j] = cut_one[1];
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairExcludedVolume::init_style()
{
  avec = (AtomVecAWSEM *) atom->style_match("awsemmd");
  if (!avec) error->all(FLERR,"Pair excluded_volume requires atom style awsemmd");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairExcludedVolume::init_one(int i, int j)
{
  // always mix prefactors geometrically

  if (setflag[i][j] == 0) {
    lambda[i][j] = sqrt(lambda[i][i]*lambda[j][j]);
    cut_short[i][j] = mix_distance(cut_short[i][i],cut_short[j][j]);
    cut_long[i][j] = mix_distance(cut_long[i][i],cut_long[j][j]);
  }

  lambda[j][i] = lambda[i][j];
  cut_short[j][i] = cut_short[i][j];
  cut_long[j][i] = cut_long[i][j];

  return cut_long[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairExcludedVolume::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&lambda[i][j],sizeof(double),1,fp);
	fwrite(&cut_short[i][j],sizeof(double),1,fp);
	fwrite(&cut_long[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairExcludedVolume::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&lambda[i][j],sizeof(double),1,fp);
	  fread(&cut_short[i][j],sizeof(double),1,fp);
	  fread(&cut_long[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&lambda[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_short[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_long[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairExcludedVolume::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global[0],sizeof(double),1,fp);
  fwrite(&cut_global[1],sizeof(double),1,fp);
  fwrite(&p,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairExcludedVolume::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global[0],sizeof(double),1,fp);
    fread(&cut_global[1],sizeof(double),1,fp);
    fread(&p,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global[0],1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global[1],1,MPI_DOUBLE,0,world);
  MPI_Bcast(&p,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairExcludedVolume::single(int i, int j, int itype, int jtype, double rsq,
			double factor_coul, double factor_lj,
			double &fforce)
{
  double r,dr,philj, rcut;
  int imol, jmol, ires, jres;

  imol = atom->molecule[i];
  jmol = atom->molecule[j];
  
  ires = avec->residue[i];
  jres = avec->residue[j];

  if (abs(ires-jres)<5 && imol==jmol) rcut = cut_short[itype][jtype];
  else rcut = cut_long[itype][jtype];

  r = sqrt(rsq);
  dr = r - rcut;
  fforce = - factor_lj * prefactor[itype][jtype] * p * pow(dr, p-1)/r;

  philj = prefactor[itype][jtype] * pow(dr, p);
  return factor_lj*philj;
}
