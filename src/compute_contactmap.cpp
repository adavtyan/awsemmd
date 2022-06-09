/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Created on: 2/22/12
   ------------------------------------------------------------------------- */

#include "mpi.h"
#include <math.h>
#include <string.h>
#include "compute_contactmap.h"
#include "atom.h"
#include "atom_vec_awsemmd.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- 

   This routine can be used to compute contact maps on the fly.

   To perform this compute, insert a line like this in your input file:
   compute 	contactmap alpha_carbons contactmap 8.5

   Note that 8.5 is an example of the distance threshold to be considered
   in contact.

   ---------------------------------------------------------------------- */

ComputeContactmap::ComputeContactmap(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // If the incorrect number of arguments are given in the compute command, quit
  // the 4 arguments that come after "compute" in the input file are:
  // compute-ID group-ID contactmap cutoff
  // cutoff is the threshold distance for two CA atoms to be considered to
  // be in contact 
  if (narg != 4) error->all(FLERR,"Illegal compute contactmap command");

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

  // create file to store contact map timeseries
  contactmapfile = fopen("contactmaptimeseries", "w");
}

/* ---------------------------------------------------------------------- */

ComputeContactmap::~ComputeContactmap()
{
  
}

/* ---------------------------------------------------------------------- */

void ComputeContactmap::init()
{
  avec = (AtomVecAWSEM *) atom->style_match("awsemmd");
  if (!avec) error->all(FLERR,"Compute contactmap requires atom style awsemmd");

  // check to make sure tags are enabled
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute contactmap unless atoms have IDs");
}

/* ---------------------------------------------------------------------- */

double ComputeContactmap::compute_scalar()
{
  // loop variables and residue numbers
  int i, j, ires, jres;
  // instantaneous distance of atoms i and j
  double rij;
  // contact map array
  int** contactmap = new int*[numres];
  for(int i = 0; i < numres; ++i)
    contactmap[i] = new int[numres];


  double **x = atom->x; // atom positions
  int *mask = atom->mask; // atom mask (?)
  int *tag = atom->tag; // atom index
  int *residue = avec->residue; // atom's residue index
  int nlocal = atom->nlocal; // number of atoms on this processor
  int nall = atom->nlocal + atom->nghost; // total number of atoms
  
  // loop over all atoms
  for (i=0;i<nlocal;i++) {
    // check to make sure the atom is in the group
    if (!(mask[i] & groupbit)) continue;
    // get residue number of atom i
    ires = residue[i]-1;
    // loop over all pairs of atoms
    for (j=i;j<nall;j++) {
      // check to make sure this atom is also in the group
      if (!(mask[j] & groupbit)) continue;
      // get residue number of atom j
      jres = residue[j]-1;
      // get the instantaneous distance
      rij=sqrt(pow(x[i][0]-x[j][0],2)+pow(x[i][1]-x[j][1],2)+pow(x[i][2]-x[j][2],2));
      // check to see if instantaneous distance is less than threshold
      if (rij<cutoff) {
	// add the contribution to tc
	contactmap[ires][jres]=1;
	contactmap[jres][ires]=1;
      }
      else {
	contactmap[ires][jres]=0;
	contactmap[jres][ires]=0;
      }    
    }
  }

  int ntimestep=update->ntimestep;

  fprintf(contactmapfile,"timestep %d\n",ntimestep); 
  // print contact map to file
  for(i=0; i<numres; i++) {
    for(j=0; j<numres; j++) {
      fprintf(contactmapfile,"%d ", contactmap[i][j]);
    }
    fprintf(contactmapfile,"\n");
  }
  
  // make a variable to return
  scalar=1; 

  // return the number of totalcontacts
  return scalar;
}
