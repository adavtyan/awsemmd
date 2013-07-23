/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer

   Wolynes Group, Rice University

   Last Update: 09/15/2011
   ------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(git,ComputeGit)

#else

#ifndef LMP_COMPUTE_GIT_H
#define LMP_COMPUTE_GIT_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeGit : public Compute {
  public:
    ComputeGit(class LAMMPS *, int, char **);
    ~ComputeGit();
    void init();
    void allocate();
    double compute_scalar();

  private:
    int igroup,groupbit;
    bool allocated;
    FILE *gitoutputfile;
    int snapshot;
  };

}

#endif
#endif
