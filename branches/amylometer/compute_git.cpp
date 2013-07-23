/* ----------------------------------------------------------------------
   Wolynes Group, Rice University

   Last Update: 07/23/13
   ------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "compute_git.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"
// External source for computing git vectors
#include "git.h"

// for reading in the native distances
#include <iostream>
#include <fstream>
using namespace std;
using namespace git;

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- 
   This routine can be used to compute GIT vectors:
   Example usage in input file:
   compute  gitvec alpha_carbons git
   variable gitvec equal c_gitvec
   fix	    gitprint all print 1 "${gitvec}" file gitvec.dat screen no title "# gitvec"
   For full details, see the pdb2git and pleiades packages in Phaistos:
   http://www.phaistos.org/home/Phaistos
   ---------------------------------------------------------------------- */

ComputeGit::ComputeGit(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // Describe inputs to compute function
  if (narg != 3) error->all(FLERR,"Illegal compute git command");

  // tell LAMMPS that we are computing a scalar
  scalar_flag = 1;
  extscalar = 1;

  allocated = false;
  igroup = group->find(arg[1]); // the ID of the group of atoms we are using (CA atoms)
  groupbit = group->bitmask[igroup]; // not really sure what this is

  // Send an error if the ID is not properly specified
  if (igroup == -1) 
    error->all(FLERR,"Could not find compute git group ID"); 
}

/* ---------------------------------------------------------------------- */

void ComputeGit::allocate()
{
  allocated = true;
}

/* ---------------------------------------------------------------------- */

ComputeGit::~ComputeGit()
{

}

/* ---------------------------------------------------------------------- */

void ComputeGit::init()
{
  // check to make sure tags are enabled
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute git unless atoms have IDs");
  snapshot = 1;
  gitoutputfile = fopen("git.dat", "w");
}

/* ---------------------------------------------------------------------- */

double ComputeGit::compute_scalar()
{
  int i;
  double **x = atom->x; // atom positions
  int *mask = atom->mask; // atom mask (?)
  int nall = atom->nlocal + atom->nghost; // total number of atoms
  
  bool smoothen_backbone_mode = 1;
  std::string average_gauss_integral_file = "";
  // Create git object
  Git git(average_gauss_integral_file, smoothen_backbone_mode);

  // Make a vector to store CA coordinates
  std::vector<Vec3> polyl;

  // loop over all atoms
  for (i=0;i<nall;i++) {
    // check to make sure the atom is in the group
    if ((mask[i] & groupbit)){
      // Load CA coordinates into vector
      polyl.push_back(Vec3(x[i][0],x[i][1],x[i][2]));
    }
  }
  // Compute GIT vector
  std::vector<double> git_vector = git.generate_gauss_integrals(polyl);
  // Write to file
  fprintf(gitoutputfile, "%5d %d %.4f ",snapshot, polyl.size()-1, 19.11*pow(polyl.size()-1,1.0/3.0));
  for(std::vector<int>::size_type i = 0; i != git_vector.size(); i++) {
    fprintf(gitoutputfile, "%.4f ", git_vector[i]);
  }
  fprintf(gitoutputfile, "\n");
  // increment snapshot
  snapshot += 1;
  // return dummy value
  return 1.0;
}
