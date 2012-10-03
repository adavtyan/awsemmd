/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Last Update: 09/23/2011
   ------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "compute_totalcontacts.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- 

   This routine can be used to compute the total number of contacts.

   To output this compute to a file:
   compute 	1 beta_atoms totalcontacts 6.5 2
   variable	tc equal c_1
   fix  		 tc all print 100 "${tc}" file tc.dat screen no

   Note that 6.5 and 2 above are examples of the distance threshold
   and sequence separation, respectively.

   ---------------------------------------------------------------------- */

ComputeTotalcontacts::ComputeTotalcontacts(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // If the incorrect number of arguments are given in the compute command, quit
  // the 5 arguments that come after "compute" in the input file are:
  // compute-ID group-ID totalcontacts cutoff sep
  // cutoff is the threshold distance for two CA atoms to be considered to
  // be in contact and sep is the minimum number of residues separating any two
  // residues that can be considered to be in contact
  if (narg != 5) error->all(FLERR,"Illegal compute totalcontacts command");

  int len; // used below to store lengths of strings
  
  // tell LAMMPS that we are computing a scalar
  scalar_flag = 1;
  extscalar = 1;

  igroup = group->find(arg[1]); // the ID of the group of atoms we are using (CA atoms)
  groupbit = group->bitmask[igroup]; // not really sure what this is

  // Send an error if the ID is not properly specified
  if (igroup == -1) 
    error->all(FLERR,"Could not find compute totalcontacts group ID"); 

  // find the number of residues in the group
  numres = (int)(group->count(igroup)+1e-12);

  // make variable based on cutoff and sep
  cutoff = atof(arg[3]);
  sep = atoi(arg[4]);
}

/* ---------------------------------------------------------------------- */

ComputeTotalcontacts::~ComputeTotalcontacts()
{
  
}

/* ---------------------------------------------------------------------- */

void ComputeTotalcontacts::init()
{
  // check to make sure tags are enabled
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute totalcontacts unless atoms have IDs");
}

/* ---------------------------------------------------------------------- */

double ComputeTotalcontacts::compute_scalar()
{
  // loop variables and residue numbers
  int i, j, ires, jres;
  // instantaneous distance of atoms i and j
  double rij;
  // total contacts
  double tc=0.0;

  double **x = atom->x; // atom positions
  int *mask = atom->mask; // atom mask (?)
  int *tag = atom->tag; // atom index
  int *residue = atom->residue; // atom's residue index
  int nlocal = atom->nlocal; // number of atoms on this processor
  int nall = atom->nlocal + atom->nghost; // total number of atoms
  
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
	// get the instantaneous distance
	rij=sqrt(pow(x[i][0]-x[j][0],2)+pow(x[i][1]-x[j][1],2)+pow(x[i][2]-x[j][2],2));
	// check to see if instantaneous distance is less than threshold
	if (rij<cutoff) {
	  // add the contribution to tc
	  tc=tc+1;
	}
      }
    }
  }

  // reduce tc across processors, store in scalar
  MPI_Allreduce(&tc,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  // return the number of totalcontacts
  return scalar;
}
