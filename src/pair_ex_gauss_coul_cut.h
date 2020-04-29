/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Gaussian Contacts Potential was contributed by Weihua Zheng

Last Update: 05/01/2012
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
  void *extract(const char *, int &);
  
  void read_parameters();
  bool isEmptyString(char *str);
  char *ltrim(char *s);
  char *rtrim(char *s);
  char *trim(char *s);

 protected:
  double cut_ex_global,cut_coul_global;
  double **cut_ex,**cut_exsq;
  double **cut_coul,**cut_coulsq;
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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style coul/cut requires atom attribute q

The atom style defined does not have these attributes.

*/
