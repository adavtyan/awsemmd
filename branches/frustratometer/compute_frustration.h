/* ----------------------------------------------------------------------
   Copyright (2011) Nick Schafer and Bobby Kim

   Wolynes Group, Rice University

   Created on: 3/9/12
   ------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(frustration,ComputeFrustration)

#else

#ifndef LMP_COMPUTE_FRUSTRATION_H
#define LMP_COMPUTE_FRUSTRATION_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeFrustration : public Compute {
  public:
    ComputeFrustration(class LAMMPS *, int, char **);
    ~ComputeFrustration();
    void init();
    //  void init_list(int, class NeighList *);
    double compute_scalar();

  private:
    double cutoff;
    int numres;
    int igroup,groupbit;
    FILE *frustrationfile;
    int numdecoys;
    char *frustrationtype;
    int randResType();
    int randResNum();
    double randPairDist();
    double randBurial();
    double computeDirectEnergy(int i_resno, int j_resno, int ires_type, int jres_type, double rij);
    double computeWaterEnergy(int i_resno, int j_resno, int ires_type, int jres_type, double rij, double burial1, double burial2);
    double computeBurialEnergy(int i_resno, int ires_type, double burial);

  };

}

#endif
#endif
