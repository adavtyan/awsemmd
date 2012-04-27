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

#ifdef PAIR_CLASS

PairStyle(ex/gauss/coul/cut,PairExGaussCoulCut)

#else

#ifndef LMP_PAIR_EX_GAUSS_COUL_CUT_H
#define LMP_PAIR_EX_GAUSS_COUL_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

struct parameters{
  int ng;
  int l; // V_ex = A/r^l;
  double A, *B, *C, *R;
  bool ex_flag;
  
  parameters(): ng(0), ex_flag(false) { B=C=R=NULL; }
  parameters(int nng) {
  	ng = nng;
  	ex_flag = false;
  	B = new double[ng];
  	C = new double[ng];
  	R = new double[ng];
  }
  ~parameters() {
  	if (B) delete [] B;
  	if (C) delete [] C;
  	if (R) delete [] R;
  }
};

class PairExGaussCoulCut : public Pair {
 public:
  PairExGaussCoulCut(class LAMMPS *);
  virtual ~PairExGaussCoulCut();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  
  void read_parameters();
  bool isEmptyString(char *str);
  char *ltrim(char *s);
  char *rtrim(char *s);
  char *trim(char *s);

 protected:
  double cut_ex_global,cut_coul_global;
  double **cut_ex,**cut_exsq;
  double **cut_coul,**cut_coulsq;
//  double **epsilon,**sigma;
  double **offset;
  
  parameters **par;
  bool el_flag;
  char *parfile;

  void allocate();
  bool read_par_flag;
  
  FILE *fout;
};

}

#endif
#endif
