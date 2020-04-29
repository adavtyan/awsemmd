/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Last Update: 09/15/2011
   ------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(qwolynes,ComputeQWolynes)

#else

#ifndef LMP_COMPUTE_QWOLYNES_H
#define LMP_COMPUTE_QWOLYNES_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeQWolynes : public Compute {
  public:
    ComputeQWolynes(class LAMMPS *, int, char **);
    ~ComputeQWolynes();
    void init();
    //  void init_list(int, class NeighList *);
    void allocate();
    double compute_scalar();
    void readNativeDistances();

  private:
    int sep;
    double sigmaexp;
    int numres;
    int igroup,groupbit;
    bool allocated;
    double **x_native;
    double **r_native;
    double qnorm;
    char *filename, *datafile;
  
    class NeighList *list;

    class AtomVecAWSEM *avec;
  };

}

#endif
#endif
