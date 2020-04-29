/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(vexcluded,PairExcludedVolume)

#else

#ifndef LMP_PAIR_EXCLUDED_VOLUME_H
#define LMP_PAIR_EXCLUDED_VOLUME_H

#include "pair.h"

namespace LAMMPS_NS {

class PairExcludedVolume : public Pair {
  friend class Pair;

 public:
  PairExcludedVolume(class LAMMPS *);
  ~PairExcludedVolume();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 private:
  double PI;
  double cut_global[2];
  int p;
  double **lambda;
  double **prefactor;
  double **cut_short;
  double **cut_long;
  double **cutsq_short;
  double **cutsq_long;

  void allocate();

  class AtomVecAWSEM *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
