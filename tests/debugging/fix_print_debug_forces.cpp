/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aram Davtyan
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_print_debug_forces.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "atom.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPrintDebugForces::FixPrintDebugForces(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // cannot have 0 atoms in group
  if (group->count(igroup) == 0) error->all(FLERR,"Fix print_debug_coords group has no atoms");

  if (narg < 5) error->all(FLERR,"Illegal fix print command");
  nstep = force->inumeric(FLERR,arg[3]);
  if (nstep < 0) error->all(FLERR,"Illegal fix print command");

  MPI_Comm_rank(world,&me);

  int n = strlen(arg[4]) + 1;
  filename = new char[n];
  strcpy(filename,arg[4]);

  x_flag = 0;
  v_flag = 0;
  q_flag = 0;

  // Read extra arguments

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"xyz") == 0) x_flag = 1;
    else if (strcmp(arg[iarg],"vel") == 0) v_flag = 1;
    else if (strcmp(arg[iarg],"charge") == 0) q_flag = 1;
    else error->all(FLERR,"Illegal fix print command");
    iarg++;
  }

  // Create the header

  int len = strlen("index type Q X Y Z VX VY VZ FX FY FZ") + 1;
  header = new char[len];
  strcpy(header, "index type");
  if (q_flag) strcat(header, " Q");
  if (x_flag) strcat(header, " X Y Z");
  if (v_flag) strcat(header, " VX VY VZ");
  strcat(header, " FX FY FZ");

  buf = new double[20];

  if (me == 0) {
    fp = fopen(filename,"w");
    if (!fp) error->all(FLERR,"Error opening file for fix_debug_coords");
  }
}

/* ---------------------------------------------------------------------- */

FixPrintDebugForces::~FixPrintDebugForces()
{
  if (fp && me == 0) fclose(fp);

  delete [] header;
  delete [] filename;
  delete [] buf;
}

/* ---------------------------------------------------------------------- */

int FixPrintDebugForces::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrintDebugForces::init()
{
}

/* ---------------------------------------------------------------------- */

void FixPrintDebugForces::end_of_step()
{
}

/* ---------------------------------------------------------------------- */

void FixPrintDebugForces::post_force(int vflag)
{
  if (update->ntimestep==nstep) {

    int j,m;
    bigint i;
    double *jquat;
    int nprocs;
    int m_rec;

    MPI_Comm_size(world,&nprocs);

    // write header
    if (me == 0) fprintf(fp,"%s\n",header);

    bigint natoms = atom->natoms;
    int nlocal = atom->nlocal;
    tagint *tag = atom->tag;
    int *type = atom->type;
    double *q = atom->q;
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    int *mask = atom->mask;

    for (i=1;i<=natoms;i++) {
      m = 0;
      for (j=0;j<nlocal;j++) {
        if (mask[j] & groupbit && i==tag[j]) {

          m = 0;
          buf[m++] = ubuf(type[j]).d;
          if (q_flag) buf[m++] = q[j]; 
          if (x_flag) {
            buf[m++] = x[j][0];
            buf[m++] = x[j][1];
            buf[m++] = x[j][2];
          }
          if (v_flag) {
            buf[m++] = v[j][0];
            buf[m++] = v[j][1];
            buf[m++] = v[j][2];
          }
          buf[m++] = f[j][0];
          buf[m++] = f[j][1];
          buf[m++] = f[j][2];

          break;
        }
      }

      // Broadcast buf if needed 

      if (me == 0) {
        for (int iproc = 1; iproc < nprocs; iproc++) {
          MPI_Recv(&m_rec,1,MPI_INT,iproc,0,world,MPI_STATUS_IGNORE);
          if (m_rec>0) {
            MPI_Recv(&buf[0],m_rec,MPI_DOUBLE,iproc,0,world,MPI_STATUS_IGNORE);
            m = m_rec;
          }
        }
      } else {
          MPI_Send(&m,1,MPI_INT,0,0,world);
          if (m>0) MPI_Send(&buf[0],m,MPI_DOUBLE,0,0,world);
      }

      // Output to the file

      if (me == 0 && m>0) {
        fprintf(fp,"%d %d",(int)i, (int)ubuf(buf[0]).i);
        for (j=1;j<m;j++) {
          fprintf(fp," %.12f", buf[j]);
        }
        fprintf(fp,"\n");
      }
    }
  }
}
