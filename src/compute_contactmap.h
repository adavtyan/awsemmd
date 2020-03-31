/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Created on: 2/22/12
   ------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(contactmap,ComputeContactmap)

#else

#ifndef LMP_COMPUTE_CONTACTMAP_H
#define LMP_COMPUTE_CONTACTMAP_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeContactmap : public Compute {
  public:
    ComputeContactmap(class LAMMPS *, int, char **);
    ~ComputeContactmap();
    void init();
    //  void init_list(int, class NeighList *);
    double compute_scalar();

  private:
    double cutoff;
    int numres;
    int igroup,groupbit;
    FILE *contactmapfile;
  
    class NeighList *list;

    class AtomVecAWSEM *avec;
  };

}

#endif
#endif
