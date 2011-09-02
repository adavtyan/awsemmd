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
  int type;
  int igroup,groupbit;
  bool allocated;
  bool **is_native;
  double r_contact, r_consq, factor;
  double **x_native;
  double **rsq_native;
  double qnorm;
  char *filename, *datafile;
  
  enum type{T_CUTOFF=0, T_SHADOW=1};
  
  FILE *fout;
  
  class NeighList *list;
};

}

#endif
#endif
