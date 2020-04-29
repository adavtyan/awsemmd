/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Created on: 2/22/12
   ------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(pairdistmat,ComputePairdistmat)

#else

#ifndef LMP_COMPUTE_PAIRDISTMAT_H
#define LMP_COMPUTE_PAIRDISTMAT_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputePairdistmat : public Compute {
  public:
    ComputePairdistmat(class LAMMPS *, int, char **);
    ~ComputePairdistmat();
    void init();
    //  void init_list(int, class NeighList *);
    double compute_scalar();

  private:
    double cutoff;
    int numres;
    int igroup,groupbit;
    FILE *pairdistmatfile;
  
    class NeighList *list;

    class AtomVecAWSEM *avec;
  };

}

#endif
#endif
