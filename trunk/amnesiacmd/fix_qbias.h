/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

#ifndef FIX_QBIAS_H
#define FIX_QBIAS_H

#include "fix.h"

//#include "smart_matrix_lib.h"

namespace LAMMPS_NS {

class FixQBias : public Fix {
 public:
  FixQBias(class LAMMPS *, int, char **);
  ~FixQBias();
//  FixBackbone(class LAMMPS *, int, char **);
//  ~FixBackbone();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

private:
  double epsilon;
  double k_qbias;

  double foriginal[4], foriginal_all[4];
  int force_flag;
  int nlevels_respa;
  bool allocated;

  int ntimestep;
  int n, nn;
  int *alpha_carbons;
  bool qbias_flag, qbias_exp_flag;
  double q0, sigma, sigma_exp;
  double *sigma_sq;
  int l;
  double **rN, **r, **q;
  int *res_no, *res_info;
  double **x, **f;
  double **xca;
  int *image;
  double prd[3], half_prd[3];
  int *periodicity;

  enum ResInfo{NONE=0, LOCAL, GHOST, OFF};

private:
  void compute();
  void compute_qbias();

  double Sigma(int sep);

  void allocate();
  int Tag(int index);
  inline void Construct_Computational_Arrays();
  inline double PeriodicityCorrection(double d, int i);
  inline bool isFirst(int index);
  inline bool isLast(int index);

  int Step;
  int sStep, eStep;
  FILE *fout;
  void out_xyz_and_force(int coord=0);
};

}

#endif
