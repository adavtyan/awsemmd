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

#ifdef FIX_CLASS

FixStyle(spring/rg/cylindrical,FixSpringRGCylindrical)

#else

#ifndef LMP_FIX_SPRING_RG_CYLINDRICAL_H
#define LMP_FIX_SPRING_RG_CYLINDRICAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringRGCylindrical : public Fix {
 public:
  FixSpringRGCylindrical(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);

 private:
  int nlevels_respa,rg0_flag;
  double rg0,k,masstotal;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
