/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/


Last Update: 08/18/2011
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(qonuchic,ComputeQOnuchic)

#else

#ifndef LMP_COMPUTE_QONUCHIC_H
#define LMP_COMPUTE_QONUCHIC_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeQOnuchic : public Compute {
 public:
  ComputeQOnuchic(class LAMMPS *, int, char **);
  ~ComputeQOnuchic();
  void init();
//  void init_list(int, class NeighList *);
  void allocate();
  double compute_scalar();
  void createContactArrays();

 private:
  int nAtoms;
  int cp_type;
  int igroup,groupbit;
  bool allocated;
  bool **is_native;
  double r_contact, r_consq, factor;
  double **x_native;
  double **rsq_native;
  double qnorm;
  double sigmaexp;
  char *datafile;
  
  enum type{T_CUTOFF=0, T_SHADOW=1, T_CUTOFF_GAUSS=2, T_SHADOW_GAUSS=3};
  
  class NeighList *list;

  class AtomVecAWSEM *avec;
};

}

#endif
#endif
