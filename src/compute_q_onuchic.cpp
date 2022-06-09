/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/


Last Update: 08/18/2011
------------------------------------------------------------------------- */

#include "mpi.h"
#include <math.h>
#include <string.h>
#include "compute_q_onuchic.h"
#include "atom.h"
#include "atom_vec_awsemmd.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

// Possible use of Compute QOnuchic
// compute 	1 alpha_carbons qonuchic cutoff r_cutoff ca_xyz_native.dat tolerance_factor
// compute 	1 alpha_carbons qonuchic shadow shadow_map_file tolerance_factor
// compute 	1 alpha_carbons qonuchic cutoff/gauss r_cutoff ca_xyz_native.dat
// compute 	1 alpha_carbons qonuchic shadow/gauss shadow_map_file

ComputeQOnuchic::ComputeQOnuchic(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 5 && narg != 6 && narg != 7) error->all(FLERR,"Illegal compute qonuchic command");

  int len;
  char *ctype;

  scalar_flag = 1;
  extscalar = 1;
  
  allocated = false;
  
  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];
  if (igroup == -1) 
		error->all(FLERR,"Could not find compute qonuchic group ID"); 
  
  len = strlen(arg[3]) + 1;
  ctype = new char[len];
  strcpy(ctype,arg[3]);
  
  if (strcmp(ctype, "cutoff")==0 || strcmp(ctype, "cutoff/gauss")==0) {
  	if (strcmp(ctype, "cutoff")==0) {
	  	cp_type = T_CUTOFF;
	  	if (narg != 7) error->all(FLERR,"Illegal compute qonuchic command");
	} else {
		cp_type = T_CUTOFF_GAUSS;
		if (narg != 6) error->all(FLERR,"Illegal compute qonuchic command");
	}

  	r_contact = atof(arg[4]);
  	r_consq = r_contact*r_contact;
 	
  	len = strlen(arg[5]) + 1;
  	datafile = new char[len];
  	strcpy(datafile,arg[5]);
  	
 	if (cp_type==T_CUTOFF) factor = atof(arg[6]);
  } else if (strcmp(ctype, "shadow")==0 || strcmp(ctype, "shadow/gauss")==0) {
  	if (strcmp(ctype, "shadow")==0) {
  		cp_type = T_SHADOW;
  		if (narg != 6) error->all(FLERR,"Illegal compute qonuchic command");
  	} else {
  		cp_type = T_SHADOW_GAUSS;
  		if (narg != 5) error->all(FLERR,"Illegal compute qonuchic command");
  	}
  	
  	len = strlen(arg[4]) + 1;
  	datafile = new char[len];
  	strcpy(datafile,arg[4]);
  	
 	if (cp_type==T_SHADOW) factor = atof(arg[5]);
  } else {
  	error->all(FLERR,"Wrong type for compute qonuchic"); 
  }
  
  sigmaexp = 0.15;
  
  createContactArrays();
  
//  fout = fopen(filename, "w");
}

/* ---------------------------------------------------------------------- */

void ComputeQOnuchic::allocate()
{
	int i;
	if (cp_type==T_CUTOFF || cp_type==T_CUTOFF_GAUSS) x_native = new double*[nAtoms];
	rsq_native = new double*[nAtoms];
	is_native = new bool*[nAtoms];
	for (i=0; i<nAtoms; ++i) {
		if (cp_type==T_CUTOFF || cp_type==T_CUTOFF_GAUSS) x_native[i] = new double[3];
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
			if (cp_type==T_CUTOFF || cp_type==T_CUTOFF_GAUSS) delete [] x_native[i];
			delete [] rsq_native[i];
			delete [] is_native[i];
		}

		if (cp_type==T_CUTOFF || cp_type==T_CUTOFF_GAUSS) delete [] x_native;
		delete [] rsq_native;
		delete [] is_native;
	}
	
//	fclose(fout);
}

/* ---------------------------------------------------------------------- */

void ComputeQOnuchic::init()
{
  avec = (AtomVecAWSEM *) atom->style_match("awsemmd");
  if (!avec) error->all(FLERR,"Compute qonuchic requires atom style awsemmd");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute qonuchic unless atoms have IDs");
}

/* ---------------------------------------------------------------------- */

void ComputeQOnuchic::createContactArrays()
{
  int i, j;

  if (cp_type==T_CUTOFF || cp_type==T_CUTOFF_GAUSS) {
	  FILE *fnative;
	  
	  fnative = fopen(datafile, "r");
	  if (!fnative) error->all(FLERR,"Compute qonuchic: can't read native coordinates");
	  fscanf(fnative, "%d",&nAtoms);
	  allocate();
	  for (i=0;i<nAtoms;++i) {
		fscanf(fnative, "%lf %lf %lf",&x_native[i][0], &x_native[i][1], &x_native[i][2]);
	  }
	  fclose(fnative);
	  
	  if (group->count(igroup)!=nAtoms)
		error->all(FLERR,"Compute qonuchic: atom number mismatch");
	  
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
			
			if (cp_type==T_CUTOFF && is_native[i][j]) rsq_native[i][j] *= factor*factor;
		}
	  }
	  qnorm = 1/qnorm;
  } else if (cp_type==T_SHADOW || cp_type==T_SHADOW_GAUSS) {
  	FILE *fshadow;
  	int ires, jres;
  	double rn;
  	
  	nAtoms = group->count(igroup);
  	allocate();
  	
  	for (i=0;i<nAtoms;++i) {
  		for (j=0;j<nAtoms;++j) {
  			is_native[i][j] = false;
  			rsq_native[i][j] = 0.0;
  		}
  	}
  	
  	qnorm = 0.0;  	
  	fshadow = fopen(datafile, "r");
  	if (!fshadow) error->all(FLERR,"Compute qonuchic: can't read contact map file");
  	while (!feof(fshadow)) {
  		fscanf(fshadow, "%d %d %lf\n",&ires, &jres, &rn);
  		ires--;
  		jres--;
  		
  		if (ires>=nAtoms || jres>=nAtoms || ires<0 || jres<0 || abs(jres-ires)<=3)
  			error->all(FLERR,"Compute qonuchic: wrong residue index in shadow contact map file");
  		
  		is_native[ires][jres] = true;
  		rsq_native[ires][jres] = rn*rn;
  		if (cp_type==T_SHADOW) rsq_native[ires][jres] *= factor*factor;
  		qnorm += 1.0;
  	}
  	qnorm = 1/qnorm;
  	fclose(fshadow);
  }
}

/* ---------------------------------------------------------------------- */

double ComputeQOnuchic::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  
  int i, j, itype, jtype, ires, jres;
  double xi[3],xj[3],delx,dely,delz,rsq;
  int xbox, ybox, zbox;
  double q=0.0, sigma_sq;

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int *tag = atom->tag;
  int *res = avec->residue;
  int *image = atom->image;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double prd[3];
  
  prd[0] = domain->xprd;
  prd[1] = domain->yprd;
  prd[2] = domain->zprd;
  
  for (i=0;i<nlocal;i++) {
    if (!(mask[i] & groupbit)) continue;
    
  	itype = type[i];
    ires = res[i]-1;

    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    
    if (domain->xperiodic) {
    	xbox = (image[i] & 1023) - 512;
		xi[0] += xbox*prd[0];
	}
	if (domain->yperiodic) {
		ybox = (image[i] >> 10 & 1023) - 512;
		xi[1] += ybox*prd[1];
	}
	if (domain->zperiodic) {
		zbox = (image[i] >> 20) - 512;
		xi[2] += zbox*prd[2];
	}
	
    
    for (j=i+1;j<nall;j++) {
      if (!(mask[j] & groupbit)) continue;
      
      jtype = type[j];
      jres = res[j]-1;
      
      if (is_native[ires][jres] && abs(jres-ires)>3) {
      	xj[0] = x[j][0];
      	xj[1] = x[j][1];
      	xj[2] = x[j][2];
      	
      	if (domain->xperiodic) {
			xbox = (image[j] & 1023) - 512;
			xj[0] += xbox*prd[0];
		}
		if (domain->yperiodic) {
			ybox = (image[j] >> 10 & 1023) - 512;
			xj[1] += ybox*prd[1];
		}
		if (domain->zperiodic) {
			zbox = (image[j] >> 20) - 512;
			xj[2] += zbox*prd[2];
		}
      
        delx = xi[0] - xj[0];
		dely = xi[1] - xj[1];
		delz = xi[2] - xj[2];
		rsq = delx*delx + dely*dely + delz*delz;
		  
	    if (cp_type==T_CUTOFF || cp_type==T_SHADOW) {
		  if (rsq<rsq_native[ires][jres]) q += 1.0;
	    } else {
		  sigma_sq = pow(1+abs(ires-jres),2.0*sigmaexp);
		  q += exp(-pow(sqrt(rsq) - sqrt(rsq_native[ires][jres]), 2)/(2.0*sigma_sq));
	    }
      }
    }
  }
  
  MPI_Allreduce(&q,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= qnorm;
  
  return scalar;
}
