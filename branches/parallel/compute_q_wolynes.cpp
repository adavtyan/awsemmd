/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Last Update: 09/15/2011
   ------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "compute_q_wolynes.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"

// for reading in the native distances
#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- 
   This routine can be used to compute Q_wolynes:
   Q_W=(2/((n-2)(n-3)))\sum_{all CA with |i-j|>sep}exp(-(rij-rijn)^2/(2*sigma^2))
   n is the number of residues
   sep is the minimum number of residues separating a pair of residues that are
   to be used in computing Q_W
   i and j are the residue indices
   sigma=|i-j|^exp where exp is the exponent specified in the compute command
   rij is the instantaneous separation of residues i and j
   rijn is the separation of residues i and j in the native state

   To output this compute to a file:
   compute 	1 alpha_carbons qwolynes rnative.dat 2 0.15
   variable	qw equal c_1
   fix		qw all print 100 "${qw}" file qw.dat screen no

   Note that 2 is the default separation and 0.15 is the default exponent.
   As written, this group is intended to be applied to the group alpha_carbons,
   as is shown above.
   ---------------------------------------------------------------------- */

ComputeQWolynes::ComputeQWolynes(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // If the incorrect number of arguments are given in the compute command, quit
  // the 6 arguments that come after "compute" in the input file are:
  // compute-ID group-ID qwolynes datafile sep exp
  // where datafile is rnative.dat, which gives the CA-CA distances
  // and sep and exp are explain above
  if (narg != 6) error->all(FLERR,"Illegal compute qwolynes command");

  int len; // used below to store lengths of strings
  
  // tell LAMMPS that we are computing a scalar
  scalar_flag = 1;
  extscalar = 1;

  allocated = false; // we haven't allocated arrays yet
  igroup = group->find(arg[1]); // the ID of the group of atoms we are using (CA atoms)
  groupbit = group->bitmask[igroup]; // not really sure what this is

  // Send an error if the ID is not properly specified
  if (igroup == -1) 
    error->all(FLERR,"Could not find compute qwolynes group ID"); 

  // find the number of residues in the group
  numres = (int)(group->count(igroup)+1e-12);

  // make datafile variable based on argument to compute function
  len = strlen(arg[3]) + 1;
  datafile = new char[len];
  strcpy(datafile,arg[3]);
  
  // make variables based on sep and exp
  sep = atoi(arg[4]);
  sigmaexp = atof(arg[5]);

  // allocate memory for arrays
  readNativeDistances();
  
}

/* ---------------------------------------------------------------------- */

void ComputeQWolynes::allocate()
{
  r_native = new double*[numres]; // an array of the native distances squared
  int i; // loop variable
  for (i=0; i<numres; ++i) {
    r_native[i] = new double[numres]; // allocate a numres x numres square matrix
  }
	
  allocated = true;
}

/* ---------------------------------------------------------------------- */

ComputeQWolynes::~ComputeQWolynes()
{
  // if allocated, deallocate the whole square matrix
  if (allocated) {
    int i;
    for (i=0;i<numres;++i) {
      delete [] r_native[i];
    }
    delete [] r_native;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeQWolynes::init()
{
  // check to make sure tags are enabled
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute qwolynes unless atoms have IDs");
}

/* ---------------------------------------------------------------------- */

void ComputeQWolynes::readNativeDistances()
{
  int i, j; // loop variables
  
  // allocate arrays
  allocate();

  ifstream in_rnative(datafile);
  if (!in_rnative) error->all(FLERR,"Native distance file can't be read");
  for (i=0;i<numres;++i)
    for (j=0;j<numres;++j)
      in_rnative >> r_native[i][j];

  in_rnative.close();

}

/* ---------------------------------------------------------------------- */

double ComputeQWolynes::compute_scalar()
{
  int i, j, ires, jres; // loop variables and residue numbers
  double rij; 			// instantaneous distance of atoms i and j
  double rijn;			// native distance of atoms i and j
  double sigmaij; 		// sigma for atoms i and j
  double q=0.0; 		// q_wolynes
  double xi[3], xj[3]; 	// temporary coordinates
  double xbox, ybox, zbox; // instantaneous periodic box indexes 
  double prd[3];		// periodic box sizes

  double **x = atom->x; // atom positions
  int *mask = atom->mask; // atom mask (?)
  int *image = atom->image; // atom image
  int *tag = atom->tag; // atom index
  int *residue = atom->residue; // atom's residue index
  int nlocal = atom->nlocal; // number of atoms on this processor
  int nall = atom->nlocal + atom->nghost; // total number of atoms
  
  // Get simulation box sizes
  prd[0] = domain->xprd;
  prd[1] = domain->yprd;
  prd[2] = domain->zprd;
  
  
  // loop over all atoms
  for (i=0;i<nlocal;i++) {
    // check to make sure the atom is in the group
    if (!(mask[i] & groupbit)) continue;
    // get residue number of atom i
    ires = residue[i]-1;
    // loop over all pairs of atoms
    for (j=i+1;j<nall;j++) {
      // check to make sure this atom is also in the group
      if (!(mask[j] & groupbit)) continue;
      // get residue number of atom j
      jres = residue[j]-1;
      // check to make sure the atoms are separated by sep residues
      if (abs(jres-ires)>sep) {
	// set sigma
	sigmaij=pow(1+abs(ires-jres),sigmaexp);
	// get the native distance
	rijn=r_native[ires][jres];
	
	// taking into account periodicity
	if (domain->xperiodic) {
    	xbox = (image[i] & 1023) - 512;
		xi[0] = x[i][0] + xbox*prd[0];
	} xi[0] = x[i][0];
	if (domain->yperiodic) {
		ybox = (image[i] >> 10 & 1023) - 512;
		xi[1] = x[i][1] + ybox*prd[1];
	} xi[1] = x[i][1];
	if (domain->zperiodic) {
		zbox = (image[i] >> 20) - 512;
		xi[2] = x[i][2] + zbox*prd[2];
	} xi[2] = x[i][2];
	
	if (domain->xperiodic) {
    	xbox = (image[j] & 1023) - 512;
		xj[0] = x[j][0] + xbox*prd[0];
	} xj[0] = x[j][0];
	if (domain->yperiodic) {
		ybox = (image[j] >> 10 & 1023) - 512;
		xj[1] = x[j][1] + ybox*prd[1];
	} xj[1] = x[j][1];
	if (domain->zperiodic) {
		zbox = (image[j] >> 20) - 512;
		xj[2] = x[j][2] + zbox*prd[2];
	} xj[2] = x[j][2];
	
	// get the instantaneous distance
	rij=sqrt(pow(xi[0]-xj[0],2)+pow(xi[1]-xj[1],2)+pow(xi[2]-xj[2],2));
	// add the contribution to q
	q += exp(-pow(rij-rijn,2)/(2.0*pow(sigmaij,2)));
      }
    }
  }
  
  // reduce q across processors, store in scalar
  MPI_Allreduce(&q,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  // normalize q
  qnorm=2.0/(((double)numres-(double)sep)*((double)numres-((double)sep+1.0)));
  scalar *= qnorm;
  // return the normalized q
  return scalar;
}
