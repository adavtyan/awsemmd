/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

// pair_style gocontacts pair_gomodel_coeff.data 24
// pair_coeff * * 24

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_go-contacts.h"
#include "atom.h"
#include "atom_vec_awsemmd.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"

#include <fstream>
#include <time.h>

using std::ifstream;

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGoContacts::PairGoContacts(LAMMPS *lmp) : Pair(lmp)
{
}

/* ---------------------------------------------------------------------- */

inline void PairGoContacts::print_log(char *line)
{
  if (screen) fprintf(screen, line);
  if (logfile) fprintf(logfile, line);
}

PairGoContacts::~PairGoContacts()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut);
    memory->destroy(cutsq);
    
    if (contacts_flag) {
    	memory->destroy(isNative);
    	memory->destroy(sigma);
    	memory->destroy(sigma_sq);
    }

	delete rand;
  }
}

/* ---------------------------------------------------------------------- */

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
void PairGoContacts::compute_contact_deviation()
{
  if (dev_type==DT_CORR) {
    if (update->ntimestep!=0) {
		rnd = rand->gaussian();
		dev += (devA*rnd - devB*dev)*devC;
	}
  } else if (dev_type==DT_SIN) {
    dev = devA*sin(devB*update->ntimestep + devC);
  } else if (dev_type==DT_CONST) {
    dev = devA;
  }
}

void PairGoContacts::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,sign,imol,jmol,ires,jres;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rsq,factor_lj,V;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double sgrinv, sgrinv12, sgrinv10, contact_epsilon;
  
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int ntypes = atom->ntypes;

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
  
  compute_contact_deviation();
  
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

      imol = atom->molecule[i]-1;
      jmol = atom->molecule[j]-1;
      
      ires = avec->residue[i]-1;
      jres = avec->residue[j]-1;
      
      if (abs(jres-ires)<=3) continue;

      if (j < nall) factor_lj = 1.0;
      else {
		factor_lj = special_lj[j/nall];
		j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      jtype = type[j];
      
      rsq = delx*delx + dely*dely + delz*delz;
      
      if (rsq < cutsq[itype][jtype] && rsq<=sigma_sq[ires][jres]*9) {
//      if (rsq < cutsq[itype][jtype]) {
//      	r = sqrt(rsq);
      	
      	sgrinv = sigma_sq[ires][jres]/rsq;
		sgrinv12 = pow(sgrinv, 6);
		sgrinv10 = pow(sgrinv, 5);
		
		if (isNative[ires][jres]) {
			contact_epsilon = epsilon;
			if (contacts_dev_flag) contact_epsilon = epsilon + dev;
			V = contact_epsilon*(5*sgrinv12 - 6*sgrinv10);
			fpair = 60*contact_epsilon*(sgrinv12 - sgrinv10)/rsq;
		} else {
			V = epsilon2*sgrinv12;
			fpair = 12*epsilon2*sgrinv12/rsq;
		}

		fpair *=factor_lj;
		V *= factor_lj;
		
		f[i][0] += fpair*delx;
		f[i][1] += fpair*dely;
		f[i][2] += fpair*delz;
		if (newton_pair || j < nlocal) {
			f[j][0] -= fpair*delx;
			f[j][1] -= fpair*dely;
			f[j][2] -= fpair*delz;
		}

		if (eflag)
		  evdwl = V;

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

void PairGoContacts::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  
  memory->create(cut, n+1,n+1,"pair:cut");
  memory->create(cutsq, n+1,n+1,"pair:cutsq");
  
  rand = new RanPark(lmp,seed);  
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairGoContacts::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  int i,j,len;

  len = strlen(arg[0]) + 1;
  parfile = new char[len];
  strcpy(parfile, arg[0]);
  
  cut_global = atof(arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut[i][j] = cut_global;
  }
  
  
  // Read the parameter file
  
  // Initalize
  lj_contacts_flag = contacts_flag = contacts_dev_flag = 0;
  dev_type = DT_NONE;
  seed = time(NULL);
//  seed = 1;
  
  char varsection[100];
  ifstream in(parfile);
  if (!in) error->all(FLERR,"Coefficient file was not found!");
  while (!in.eof()) {
  	in >> varsection;
  	if (strcmp(varsection, "[Go-Model_LJ]")==0) {
  		lj_contacts_flag = 1;
		print_log("LJ Go-Model flag on\n");
		in >> epsilon >> epsilon2;
  	} else if (strcmp(varsection, "[Contacts]")==0) {
  		contacts_flag = 1;
  		print_log("Contacts flag on\n");
  		in >> nres >> ncont;
  		in >> sigma0;
  			
  		memory->create(isNative, nres+1,nres+1,"pair:isNative");
  		memory->create(sigma, nres+1,nres+1,"pair:sigma");
  		memory->create(sigma_sq, nres+1,nres+1,"pair:sigma_sq");
  		
  		for (i=0;i<nres;++i) for (j=0;j<nres;++j) { isNative[i][j]=0; sigma[i][j]=sigma0; sigma_sq[i][j]=sigma0*sigma0; }
  		
  		int ires, jres;
  		double rn;
  		for (i=0; i<ncont;++i) {
  			in >> ires >> jres >> rn;
  			
  			if (ires>=nres || jres>=nres || ires<0 || jres<0 || abs(jres-ires)<=3)
  				error->all(FLERR,"Pair style go-contacts: wrong residue index in contact map file");
  			
  			isNative[ires][jres] = isNative[jres][ires] = 1;
  			sigma[ires][jres] = sigma[jres][ires] = rn;
  			sigma_sq[ires][jres] = sigma_sq[jres][ires] = rn*rn;
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
	
	fout = fopen("forcesPGO.dat","w");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGoContacts::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  double cut_one;
  
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  cut_one = cut_global;
  if (narg == 3)
    cut_one = atof(arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGoContacts::init_style()
{
  avec = (AtomVecAWSEM *) atom->style_match("awsemmd");
  if (!avec) error->all(FLERR,"Pair go-contacts requires atom style awsemmd");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGoContacts::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGoContacts::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
		fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGoContacts::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  if (!allocated) allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&cut[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGoContacts::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  
  fwrite(&dev_type,sizeof(int),1,fp);
  if (dev_type==DT_CORR) {
	int tseed=rand->state();
  	fwrite(&tseed,sizeof(int),1,fp);
  	fwrite(&sdivf,sizeof(double),1,fp);
  	fwrite(&tcorr,sizeof(double),1,fp);
  	fwrite(&dev0,sizeof(double),1,fp);
  	fwrite(&dev,sizeof(double),1,fp);
  } else if (dev_type==DT_SIN) {
  	double tdev0=devB*update->ntimestep + devC;
  	fwrite(&sdivf,sizeof(double),1,fp);
  	fwrite(&tcorr,sizeof(double),1,fp);
  	fwrite(&tdev0,sizeof(double),1,fp);
  	fwrite(&dev,sizeof(double),1,fp);
  } else if (dev_type==DT_CONST) {
  	fwrite(&sdivf,sizeof(double),1,fp);
  	fwrite(&dev,sizeof(double),1,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGoContacts::read_restart_settings(FILE *fp)
{
  if (!allocated) allocate();

  int r_dev_type, r_seed;
  double r_sdivf, r_tcorr, r_dev0, r_dev;
  if (comm->me == 0) {
	fread(&cut_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    
    fread(&r_dev_type,sizeof(int),1,fp);
    if (r_dev_type==DT_CORR) {
    	fread(&r_seed,sizeof(int),1,fp);
    	fread(&r_sdivf,sizeof(double),1,fp);
    	fread(&r_tcorr,sizeof(double),1,fp);
    	fread(&r_dev0,sizeof(double),1,fp);
    	fread(&r_dev,sizeof(double),1,fp);
    } else if (r_dev_type==DT_SIN) {
    	fread(&r_sdivf,sizeof(double),1,fp);
    	fread(&r_tcorr,sizeof(double),1,fp);
    	fread(&r_dev0,sizeof(double),1,fp);
    	fread(&r_dev,sizeof(double),1,fp);
    } else if (r_dev_type==DT_CONST) {
    	fread(&r_sdivf,sizeof(double),1,fp);
    	fread(&r_dev,sizeof(double),1,fp);
    }
    
    if (contacts_dev_flag && r_dev_type==dev_type) {
    	if (r_dev_type==DT_CORR) {
    		rand->reset(r_seed);
    		dev = r_dev;
    	} else if (r_dev_type==DT_SIN) {
    		devC = r_dev0;
    		dev = devA*sin(devC);
    	} else if (r_dev_type==DT_CONST) {
    		// Nothing need to be done
    	}
    }
  }
  
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairGoContacts::single(int i, int j, int itype, int jtype, double rsq,
			double factor_coul, double factor_lj,
			double &fforce)
{
  double r,V,contact_epsilon;
  double sgrinv, sgrinv12, sgrinv10;
  int imol, jmol, ires, jres;

  imol = atom->molecule[i]-1;
  jmol = atom->molecule[j]-1;
  
  ires = avec->residue[i]-1;
  jres = avec->residue[j]-1;
  
// || rsq>sigma_sq[ires][jres]*9
  if (abs(jres-ires)<=3 || rsq>=cutsq[itype][jtype] || rsq>sigma_sq[ires][jres]*9) return 0.0;

//  r = sqrt(rsq);
  sgrinv = sigma_sq[ires][jres]/rsq;
  sgrinv12 = pow(sgrinv, 6);
  sgrinv10 = pow(sgrinv, 5);
  
  if (isNative[ires][jres]) {
    contact_epsilon = epsilon;
	if (contacts_dev_flag) contact_epsilon = epsilon + dev;
	V = contact_epsilon*(5*sgrinv12 - 6*sgrinv10);
	fforce = 60*contact_epsilon*(sgrinv12 - sgrinv10)/rsq;
  } else {
	V = epsilon2*sgrinv12;
	fforce = 12*epsilon2*sgrinv12/rsq;
  }
  fforce *= factor_lj;
  V *= factor_lj;

  return V;
}
