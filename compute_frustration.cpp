/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer and Bobby Kim

   Wolynes Group, Rice University

   Created on: 3/9/12
   ------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "compute_frustration.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- 

   This routine can be used to compute local frustration.

   ADD NOTES WHEN COMPLETE

   ---------------------------------------------------------------------- */

ComputeFrustration::ComputeFrustration(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  // If the incorrect number of arguments are given in the compute command, quit
  // How to call this compute:
  // compute frustratometer all frustration type contactcutoff numberofdecoys
  if (narg != 6) error->all("Illegal compute frustration command");

  int len; // used below to store lengths of strings
  
  // tell LAMMPS that we are computing a scalar
  scalar_flag = 1;
  extscalar = 1;

  igroup = group->find(arg[1]); // the ID of the group of atoms we are using (CA atoms)
  groupbit = group->bitmask[igroup]; // not really sure what this is

  // Send an error if the ID is not properly specified
  if (igroup == -1) 
    error->all("Could not find compute frustration group ID"); 

  // find the number of residues in the group
  numres = (int)(group->count(igroup)+1e-12);

  // type of frustration calculation: configurational, pairmutational, singlemutational
  len = strlen(arg[3]) + 1;
  frustrationtype = new char[len];
  strcpy(frustrationtype,arg[3]);

  // distance cutoff
  cutoff = atof(arg[4]);

  // number of decoys
  numdecoys = atoi(arg[5]);

  // create file to store frustration timeseries
  frustrationfile = fopen("frustrationtimeseries", "w");
}

/* ---------------------------------------------------------------------- */

ComputeFrustration::~ComputeFrustration()
{
  
}

/* ---------------------------------------------------------------------- */

void ComputeFrustration::init()
{
  // check to make sure tags are enabled
  if (atom->tag_enable == 0)
    error->all("Cannot use compute frustration unless atoms have IDs");
}

/* ---------------------------------------------------------------------- */

double ComputeFrustration::compute_scalar()
{
  // loop variables and residue numbers
  int i, j, k, ires, jres;
  // instantaneous distance of atoms i and j
  double rij;
  // frustration map array
  double** frustrationmap = new double*[numres];
  for(i = 0; i < numres; ++i){
    frustrationmap[i] = new double[numres];
    for(j = 0; j < numres; ++j){
      frustrationmap[i][j]=0.0;
    }
  }
  // decoy energy array
  double* decoyEnergyArray = new double[numdecoys];
  for(i = 0; i < numdecoys; ++i){
      decoyEnergyArray[i] = 0.0;
  }
  double **x = atom->x; // atom positions
  int *mask = atom->mask; // atom mask (?)
  int *tag = atom->tag; // atom index
  int *residue = atom->residue; // atom's residue index
  int nlocal = atom->nlocal; // number of atoms on this processor
  int nall = atom->nlocal + atom->nghost; // total number of atoms  
  double dx[3]; // difference vector
  double *xi, *xj; // pointers to i and j positions
  int iatom, jatom; // atom of interest in residue i,j
  int i_resno,j_resno; // the residue number of residue i,j
  int i_chno,jchno; // the chain number of residue i,j
  int ires_type,jres_type; // the residue type of residue i,j
  
  double nativeEnergy = 0.0;
  double burial1,burial2;

  int decoyResNum1,decoyResNum2;
  double decoyEnergy = 0.0;
  int decoyResType1,decoyResType2;
  double decoyDist;
  double decoyBurial1,decoyBurial2;
  double totalDecoyEnergy = 0.0;
  double averageDecoyEnergy = 0.0;
  double decoyEnergyVariance = 0.0;
  double energyGap = 0.0;

  // Loop over all residue pairs
  for (i=0;i<numres;i++) {
    // get residue i information
    ires = residue[i]-1;
    i_resno = res_no[i]-1;
    i_chno = chain_no[i]-1;
    ires_type = se_map[se[i_resno]-'A'];
    if (se[i_resno]=='G') { xi = xca[i]; iatom = alpha_carbons[i]; }
    else { xi = xcb[i]; iatom  = beta_atoms[i]; }
    // loop over residue index j
    for (j=i+1;j<numres;j++) {
      // get residue j information
      jres = residue[j]-1;
      j_resno = res_no[j]-1;
      j_chno = chain_no[j]-1;
      jres_type = se_map[se[j_resno]-'A'];
      if (se[j_resno]=='G') { xj = xca[j]; jatom = alpha_carbons[j]; }
      else { xj = xcb[j]; jatom  = beta_atoms[j]; }
      // get the instantaneous distance
      dx[0] = xi[0] - xj[0];
      dx[1] = xi[1] - xj[1];
      dx[2] = xi[2] - xj[2];
      rij=sqrt(pow(dx[0],2)+pow(dx[1],2)+pow(dx[2],2));
      // check to see if instantaneous distance is less than threshold
      // If residues are within cutoff distance, compute frustration index
      if (rij<cutoff) {
	// Compute pairwise energy - the sum of direct, water and burial
	nativeEnergy = 0.0; // reset native energy
	// Compute direct contact interaction
	nativeEnergy += computeDirectEnergy(i_resno,j_resno,ires_type,jres_type,rij);
	// Compute degree of burial for each residue
	burial1 = 0.0;
	burial2 = 0.0;
	// Compute water interaction
	nativeEnergy += computeWaterEnergy(i_resno,j_resno,ires_type,jres_type,rij,burial1,burial2);
	// Compute burial interaction for both residues
	nativeEnergy += computeBurialEnergyEnergy(i_resno,ires_type,burial1);
	nativeEnergy += computeBurialEnergyEnergy(j_resno,jres_type,burial2);
	// Compute decoy energies
	for (k=0;k<numdecoys;k++) {
	  // Choose two random residue positions
	  decoyResNum1 = randResNum();
	  decoyResNum2 = randResNum();
	  // Choose two random amino acids based on protein composition
	  decoyResType1 = randResType();
	  decoyResType2 = randResType();
	  // Choose a random pairwise distance from current protein configuration
	  decoyDist = randPairDist();
	  // Choose two random degrees of burial based on current protein configuration
	  decoyBurial1 = randBurial();
	  decoyBurial2 = randBurial();
	  // Compute pairwise energy - the sum of direct, water and burial
	  decoyEnergy = 0.0; // reset decoy energy
	  // Compute direct contact interaction
	  decoyEnergy += computeDirectEnergy(decoyResNum1,decoyResNum2,decoyResType1,decoyResType2,decoyDist);
	  // Compute water interaction
	  decoyEnergy += computeWaterEnergy(decoyResNum1,decoyResNum2,decoyResType1,decoyResType2,decoyDist,decoyBurial1,decoyBurial2);
	  // Compute burial interaction for both residues
	  decoyEnergy += computeBurialEnergyEnergy(decoyResNum1,decoyResType1,decoyBurial1);
	  decoyEnergy += computeBurialEnergyEnergy(decoyResNum2,decoyResType2,decoyBurial2);
	  // Store decoy energy in decoy energy array
	  decoyEnergyArray[k] = decoyEnergy;
	} // Repeat for all decoys
	// Use equation to compute frustration index for each i,j pair in contact
	// Take average of decoy energy array
	totalDecoyEnergy = 0.0;
	for (k=0;k<numdecoys;k++) {
	  totalDecoyEnergy += decoyEnergyArray[k];
	}
	averageDecoyEnergy = totalDecoyEnergy/(double)numdecoys;
	// Take variance of decoy energy array
	decoyEnergyVariance = 0.0;
	for (k=0;k<numdecoys;k++) {
	  decoyEnergyVariance += pow(decoyEnergyArray[k]-averageDecoyEnergy,2);
	}
	decoyEnergyVariance = decoyEnergyVariance/(double)numdecoys;
	// Compute energy gap = native energy - decoy average
	energyGap = nativeEnergy-averageDecoyEnergy;
	// Divide energy gap by standard deviation of decoy energies
	frustrationmap[i][j]= energyGap/sqrt(decoyEnergyVariance);
	frustrationmap[j][i]=frustrationmap[i][j];
      }
    }
  }

  // write results to output file
  int ntimestep=update->ntimestep;

  fprintf(frustrationfile,"timestep %d\n",ntimestep); 
  // print contact map to file
  for(i=0; i<numres; i++) {
    for(j=0; j<numres; j++) {
      fprintf(frustrationfile,"%1.2f ", frustrationmap[i][j]);
    }
    fprintf(frustrationfile,"\n");
  }
  
  // make a variable to return
  scalar=1.0; 

  // return the number of totalcontacts
  return scalar;
}

// randResNum
int ComputeFrustration::randResNum()
{
  // generate random integer from 0 to number of residues
  return 0;
}

// randResType
int ComputeFrustration::randResType(int randResNum)
{
  // return sequence[randResNum]
  return 'A';
}

// randPairDist
double ComputeFrustration::randPairDist(int randResNum1, int randResNum2)
{
  // return distance between residues randResNum1 and randResNum2
  return 5.0;
}

// randBurial
double ComputeFrustration::randBurial(int randResNum)
{
  // return degree of burial of randResNum in current configuration
  return 2.0;
}

// computeDirectEnergy
double ComputeFrustration::computeDirectEnergy(int i_resno, int j_resno, int ires_type, int jres_type, double rij)
{
  double gamma = get_water_gamma(i_resno,j_resno,0,ires_type,jres_type,0);
  double theta = (1+tanh(water_kappa*(rij-well_r_min[0])))*(1+tanh(water_kappa*(well_r_max[0]-rij)));
  double directEnergy = -k_water*gamma*0.25*theta;
  return directEnergy;
}

// computeWaterEnergy
double ComputeFrustration::computeWaterEnergy(int i_resno, int j_resno, int ires_type, int jres_type, double rij, double burial1, double burial2)
{
  double theta = (1+tanh(water_kappa*(rij-well_r_min[1])))(1+tanh(water_kappa*(well_r_max[1]-rij)));
  double sigmawat = 0.25*(1-tanh(water_kappa_sigma*(burial1-treshold)))*(1-tanh(water_kappa_sigma(burial2-treshold)));
  double gammawat = get_water_gamma(i_resno,j_resno,1,ires_type,jres_type,0);
  double sigmaprot = 1.0-sigmawat;
  double gammaprot = get_water_gamma(i_resno,j_resno,1,ires_type,jres_type,1);
  double waterEnergy = -k_water*theta*(sigmawat*gammawat+sigmaprot*gammaprot);
  return waterEnergy;
}

// computeBurialEnergyEnergy
double ComputeFrustration::computeBurialEnergy(int i_resno, int ires_type, double burial)
{
  int mu = 0;
  double burialEnergy = 0.0;
  double theta = 0.0;
  for(mu=0;mu<3;mu++) {
    theta = (tanh(burial_kappa*(burial-burial_ro_min[mu])))*(tanh(burial_kappa*(burial_ro_max[mu]-burial)));
    burialEnergy += -0.5*get_burial_gamma(i_resno,ires_type,mu)*theta;
  }
  return burialEnergy;
}
