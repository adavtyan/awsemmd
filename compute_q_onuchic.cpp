/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/


Last Update: 08/18/2011
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "compute_q_onuchic.h"
#include "atom.h"
#include "update.h"
//#include "force.h"
//#include "pair.h"
#include "domain.h"
//#include "neighbor.h"
//#include "neigh_request.h"
//#include "neigh_list.h"
#include "group.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeQOnuchic::ComputeQOnuchic(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 7) error->all("Illegal compute qonuchic command");

  scalar_flag = 1;
  extscalar = 1;
  
  allocated = false;
  
  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];
  if (igroup == -1) 
		error->all("Could not find compute qonuchic group ID"); 
  
  int len = strlen(arg[3]) + 1;
  filename = new char[len];
  strcpy(filename,arg[3]);
  
  r_contact = atof(arg[5]);
  factor = atof(arg[6]);
  
  r_consq = r_contact*r_contact; 
  
  int i,j;
  FILE *fnative;
  
  fnative = fopen(arg[4], "r");
  if (!fnative) error->all("Compute qonuchic: can't read native coordinates");
  fscanf(fnative, "%d",&nAtoms);
  allocate();
  for (i=0;i<nAtoms;++i) {
  	fscanf(fnative, "%lf %lf %lf",&x_native[i][0], &x_native[i][1], &x_native[i][2]);
  }
  fclose(fnative);
  
  if (group->count(igroup)!=nAtoms)
  	error->all("Compute qonuchic: atom number mismatch");
  
  double dx[3];
  qnorm = 0.0;
  for (i=0;i<nAtoms;++i) {
  	for (j=0;j<nAtoms;++j) {
  		dx[0] = x_native[j][0] - x_native[i][0];
  		dx[1] = x_native[j][1] - x_native[i][1];
  		dx[2] = x_native[j][2] - x_native[i][2];
  		
  		rsq_native[i][j] = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
  		
  		if (rsq_native[i][j] < r_consq && j-i>3) {
  			is_native[i][j] = true;
  			qnorm += 1.0;
  		} else is_native[i][j] = false;
  		  		
  		if (is_native[i][j]) rsq_native[i][j] *= factor*factor;  		
  	}
  }
  qnorm = 1/qnorm;
  
  fout = fopen(filename, "w");
}

/* ---------------------------------------------------------------------- */

void ComputeQOnuchic::allocate()
{
	int i;
	x_native = new double*[nAtoms];
	rsq_native = new double*[nAtoms];
	is_native = new bool*[nAtoms];
	for (i=0; i<nAtoms; ++i) {
		x_native[i] = new double[3];
		rsq_native[i] = new double[nAtoms];
		is_native[i] = new bool[nAtoms];
	}
	
	allocated = true;
}

/* ---------------------------------------------------------------------- */

ComputeQOnuchic::~ComputeQOnuchic()
{
	if (allocated) {
		int i;
		for (i=0;i<nAtoms;++i) {
			delete [] x_native[i];
			delete [] rsq_native[i];
			delete [] is_native[i];
		}

		delete [] x_native;
	}
	
	fclose(fout);
}

/* ---------------------------------------------------------------------- */

void ComputeQOnuchic::init()
{
  if (atom->tag_enable == 0)
    error->all("Cannot use compute qonuchic unless atoms have IDs");
/*  if (force->pair == NULL)
    error->all("Compute qonuchic requires a pair style be defined");
  if (r_contact > force->pair->cutforce)
    error->all("Compute qonuchic cutoff is longer than pairwise cutoff");

  // need an occasional half neighbor list
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;
  neighbor->requests[irequest]->occasional = 1;*/
}

/* ---------------------------------------------------------------------- */

/*void ComputeQOnuchic::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeQOnuchic::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  
  int i, j, itype, jtype, ires, jres;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double q=0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  
  for (i=0;i<nlocal;i++) {
    if (!(mask[i] & groupbit)) continue;
    
  	itype = type[i];
    ires = tag[i]-1;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    
    for (j=i+1;j<nall;j++) {
      if (!(mask[j] & groupbit)) continue;
      
      jtype = type[j];
      jres = tag[j]-1;
      
      if (is_native[ires][jres] && abs(jres-ires)>3) {
        delx = xtmp - x[j][0];
		dely = ytmp - x[j][1];
		delz = ztmp - x[j][2];
		rsq = delx*delx + dely*dely + delz*delz;
		  
		if (rsq<rsq_native[ires][jres]) q += 1.0;
      }
    }
  }
  
  MPI_Allreduce(&q,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= qnorm;
  
  fprintf(fout, "%d %f\n", invoked_scalar, scalar);
  
  return scalar;
}

/*double ComputeQOnuchic::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  
  int i, j, ii, jj, itype, jtype, inum, jnum, ires, jres;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double q=0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  
  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list->index);
  
  fprintf(screen, "force->pair->cutforce=%f\n", force->pair->cutforce);
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  for (ii = 0; ii < inum; ii++) {  	
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    ires = tag[i]-1;
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];
      jres = tag[j]-1;
      
      delx = xtmp - x[j][0];
		dely = ytmp - x[j][1];
		delz = ztmp - x[j][2];
		rsq = delx*delx + dely*dely + delz*delz;
      
      fprintf(fout, "%d %d %d %d %f\n", i, j, ires, jres, sqrt(rsq));
      
      if (is_native[ires][jres] && abs(jres-ires)>3) {
		delx = xtmp - x[j][0];
		dely = ytmp - x[j][1];
		delz = ztmp - x[j][2];
		rsq = delx*delx + dely*dely + delz*delz;
		  
		if (rsq<rsq_native[ires][jres]) q += 1.0;
	  }
    }	
  }
  
  fprintf(fout, "%d %f\n", invoked_scalar, q);
  fprintf(screen, "%d %f ", invoked_scalar, q);

  MPI_Allreduce(&q,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= qnorm;
  
  fprintf(fout, "%d %f\n", invoked_scalar, scalar);
  fprintf(screen, "%f\n", scalar);
  
  return scalar;
}*/
