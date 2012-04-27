/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Gaussian Contacts Potential was contributed by Weihua Zheng

Last Update: 03/23/2011
------------------------------------------------------------------------- */
#ifdef FIX_CLASS
FixStyle(gomodel, FixGoModel)
#else

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
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  void write_restart(FILE *);
  void restart(char *);
  double compute_scalar();
  double compute_vector(int);

private:
  int force_flag;
  int nlevels_respa;
  int ntimestep;
  int n, nn;
  int *alpha_carbons;
  int *res_no, *res_info, *chain_no;
  int *image;
  int *periodicity;
  int seed;

  double epsilon, epsilon2;
  double k_bonds, k_angles, k_dihedrals[2];
  double *r0, *theta0, *phi0, **sigma, **sigma_sq;
  double dev, devA, devB, devC;
  double sdivf, tcorr, dev0;
  double w, xi;
  double **x, **f;
  double **xca;
  double prd[3], half_prd[3];
  double rand;  

   //Gaussian contacts, for multi-basin
  int n_basins;
  double R, *G, *A, ***sigma_mb;
  double gaussian_width;
  double g_w_sq_inv;
  double rmin_cutoff;

  bool allocated, contacts_allocated;
  bool bonds_flag, angles_flag, dihedrals_flag, contacts_flag, contacts_dev_flag, lj_contacts_flag, gaussian_contacts_flag;
  bool **isNative, ***isNative_mb; 
  int dev_type;

  RanPark *random;

  enum ResInfo{NONE=0, LOCAL, GHOST, OFF};
  
  double energy[5], energy_all[5];
  enum EnergyTerms{ET_TOTAL=0, ET_BOND, ET_ANGLE, ET_DIHEDRAL, ET_CONTACTS, ET_NCONTS, nEnergyTerms};
  
  enum ContactsDevType{DT_NONE=0, DT_CORR, DT_SIN, DT_CONST};

private:
  void compute_goModel();
  void compute_bond(int i);
  void compute_angle(int i);
  void compute_dihedral(int i);
  void compute_contact(int i, int j);
  void compute_contact_deviation();
  void compute_contact_gaussian(int i, int j);

  void allocate();
  void allocate_contact();

  int Tag(int index);
  inline void Construct_Computational_Arrays();
  inline double PeriodicityCorrection(double d, int i);
  inline bool isFirst(int index);
  inline bool isLast(int index);
  inline void print_log(char *line);

  int Step;
  int sStep, eStep;
  FILE *fout;
  FILE *efout;
  FILE *efile;
  void out_xyz_and_force(int coord=0);
};

}

#endif
#endif
