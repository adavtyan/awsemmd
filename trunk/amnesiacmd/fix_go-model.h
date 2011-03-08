/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

#ifndef FIX_GOMODEL_H
#define FIX_GOMODEL_H

#include "fix.h"
#include "random_park.h"

//#include "smart_matrix_lib.h"

namespace LAMMPS_NS {

class FixGoModel : public Fix {
 public:
  FixGoModel(class LAMMPS *, int, char **);
  ~FixGoModel();
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
  void write_restart(FILE *);
  void restart(char *);

private:
  double epsilon, epsilon2;
  double k_bonds, k_angles, k_dihedrals[2];

  double foriginal[4], foriginal_all[4];
  int force_flag;
  int nlevels_respa;
  bool allocated;

  int ntimestep;
  int n, nn;
  int *alpha_carbons;
  bool bonds_flag, angles_flag, dihedrals_flag, contacts_flag, contacts_dev_flag;
  double *r0, *theta0, *phi0, **sigma;
  double dev, devA, devB, devC;
  double sdivf, tcorr, dev0;
  double w, xi;
  bool **isNative;
  int *res_no, *res_info;
  double **x, **f;
  double **xca;
  int *image;
  double prd[3], half_prd[3];
  int *periodicity;

  RanPark *random;
  int seed;
  double rand;

  enum ResInfo{NONE=0, LOCAL, GHOST, OFF};

private:
  void compute_goModel();
  void compute_bond(int i);
  void compute_angle(int i);
  void compute_dihedral(int i);
  void compute_contact(int i, int j);
  void compute_contact_deviation();

  void allocate();
  int Tag(int index);
  inline void Construct_Computational_Arrays();
  inline double PeriodicityCorrection(double d, int i);
  inline bool isFirst(int index);
  inline bool isLast(int index);

  int Step;
  int sStep, eStep;
  FILE *fout;
  FILE *efout;
  void out_xyz_and_force(int coord=0);
};

}

#endif
