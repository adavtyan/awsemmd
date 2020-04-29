/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*-----------------------------------------------------------------------
    Rg biasing potential
    Modified from fix_spring.h
    Author: Hao Wu and Garegin Papoian
    Papoian Research Group
    University of Maryland, College Park
    http://papoian.chem.umd.edu/

    Please cite:
    "AWSEM-IDP: A Coarse-Grained Force Field for Intrinsically Disordered Proteins"
    Hao Wu, Peter G. Wolynes and Garegin A. Papoian, J. Phys. Chem. B, 2018
    Check out https://github.com/adavtyan/awsemmd for detailed usage and examples

    Last update: Aug 15 2018
 -------------------------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(spring/rg/papoian,FixSpringRGPapoian)

#else

#ifndef LMP_FIX_SPRING_RG_PAPOIAN_H
#define LMP_FIX_SPRING_RG_PAPOIAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpringRGPapoian : public Fix {
 public:
  FixSpringRGPapoian(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);

 private:
  int ilevel_respa,rg0_flag;
  int N;
  double rg0,gamma,D,alpha,beta;
  double masstotal;
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
