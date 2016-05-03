/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Last Update: 09/23/2011
   ------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(totalcontacts,ComputeTotalcontacts)

#else

#ifndef LMP_COMPUTE_TOTALCONTACTS_H
#define LMP_COMPUTE_TOTALCONTACTS_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeTotalcontacts : public Compute {
  public:
    ComputeTotalcontacts(class LAMMPS *, int, char **);
    ~ComputeTotalcontacts();
    void init();
    //  void init_list(int, class NeighList *);
    double compute_scalar();

  private:
    double cutoff;
    int sep;
    int numres;
    int igroup,groupbit;
  
    class NeighList *list;
  };

}

#endif
#endif
