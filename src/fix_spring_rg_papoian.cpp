/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
                        Paul Crozier (SNL)
------------------------------------------------------------------------- */

/*-----------------------------------------------------------------------
    Rg biasing potential
    Modified from fix_spring.cpp
    Author: Hao Wu and Garegin Papoian
    Papoian Research Group
    University of Maryland, College Park
    http://papoian.chem.umd.edu/

    Please cite:
    "AWSEM-IDP: A Coarse-Grained Force Field for Intrinsically Disordered Proteins"
    Hao Wu, Peter G. Wolynes and Garegin A. Papoian, J. Phys. Chem. B, 2018
    Check out https://github.com/adavtyan/awsemmd for detailed usage and examples

    Last update: Aug 15 2018
    Output e_rg: print Vrg energy value
 -------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_spring_rg_papoian.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpringRGPapoian::FixSpringRGPapoian(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // Input: ID group-ID spring/rg/papoian N rg0 gamma D alpha beta
  if (narg != 9) error->all(FLERR,"Illegal fix spring/rg/papoian command");

  N = utils::numeric(FLERR,arg[3],false,lmp); // N: residue number
  rg0_flag = 0;
  if (strcmp(arg[4],"NULL") == 0) rg0_flag = 1;
  else rg0 = utils::numeric(FLERR,arg[4],false,lmp); // rg0: rg at bottom of well
  gamma = utils::numeric(FLERR,arg[5],false,lmp); // gamma: rescaling factor for rg0
  D = utils::numeric(FLERR,arg[6],false,lmp); // D: well depth
  alpha = utils::numeric(FLERR,arg[7],false,lmp); // alpha: width between two maximum
  beta = utils::numeric(FLERR,arg[8],false,lmp); // beta: well width

  dynamic_group_allow = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

int FixSpringRGPapoian::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringRGPapoian::init()
{
  masstotal = group->mass(igroup);

  // if rg0 was specified as NULL, compute current Rg
  // only occurs on 1st run

  if (rg0_flag) {
    double xcm[3];
    group->xcm(igroup,masstotal,xcm);
    rg0 = group->gyration(igroup,masstotal,xcm);
    rg0_flag = 0;
  }

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>( update->integrate))->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringRGPapoian::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>( update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    (dynamic_cast<Respa *>( update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringRGPapoian::post_force(int vflag)
{
  // compute current Rg and center-of-mass

  double xcm[3];
  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,xcm);
  double rg = group->gyration(igroup,masstotal,xcm);

  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double dx,dy,dz,term1;
  double drg; // define drg = rg - rg0

  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double massfrac;
  double unwrap[3];

  // define drg
  drg = rg-gamma*rg0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];
      // term1: dV/drg/rg (demand of integration by parts)
      // term1 = 2.0 * k * (1.0 - rg0/rg);
      term1 = (2*alpha*drg-2*alpha*beta*pow(drg,5)-4*beta*D*N*pow(drg,3))/(1+beta*pow(drg,4))/(1+beta*pow(drg,4))/rg;
      if (masstotal > 0.0) {
        if (rmass) massfrac = rmass[i]/masstotal;
        else  massfrac = mass[type[i]]/masstotal;

        f[i][0] -= term1*dx*massfrac;
        f[i][1] -= term1*dy*massfrac;
        f[i][2] -= term1*dz*massfrac;
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixSpringRGPapoian::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}
