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
#include "fix_print_wzero.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPrintWithZero::FixPrintWithZero(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix print/wzero command");
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix print/wzero command");

  MPI_Comm_rank(world,&me);

  int n = strlen(arg[4]) + 1;
  string = new char[n];
  strcpy(string,arg[4]);

  copy = (char *) memory->smalloc(n*sizeof(char),"fix/print/wzero:copy");
  work = (char *) memory->smalloc(n*sizeof(char),"fix/print/wzero:work");
  maxcopy = maxwork = n;

  // parse optional args

  fp = NULL;
  screenflag = 1;
  char *title = NULL;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print/wzero command");
      if (me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix print/wzero file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print/wzero command");
      if (strcmp(arg[iarg+1],"yes") == 0) screenflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) screenflag = 0;
      else error->all(FLERR,"Illegal fix print/wzero command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"title") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print/wzero command");
      delete [] title;
      int n = strlen(arg[iarg+1]) + 1;
      title = new char[n];
      strcpy(title,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix print/wzero command");
  }

  // print file comment line

  if (fp && me == 0) {
    if (title) fprintf(fp,"%s\n",title);
    else fprintf(fp,"# Fix print/wzero output for fix %s\n",id);
  }

  delete [] title;

  // add nfirst to all computes that store invocation times
  // since don't know a priori which are invoked via variables by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  int nfirst = (update->ntimestep/nevery)*nevery + nevery;
  modify->addstep_compute_all(nfirst);
}

/* ---------------------------------------------------------------------- */

FixPrintWithZero::~FixPrintWithZero()
{
  delete [] string;
  memory->sfree(copy);
  memory->sfree(work);

  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixPrintWithZero::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrintWithZero::end_of_step()
{
  // make a copy of string to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  strcpy(copy,string);
  input->substitute(copy,work,maxcopy,maxwork,0);

  modify->addstep_compute(update->ntimestep + nevery);

  if (me == 0) {
    if (screenflag && screen) fprintf(screen,"%s\n",copy);
    if (screenflag && logfile) fprintf(logfile,"%s\n",copy);
    if (fp) {
      fprintf(fp,"%s\n",copy);
      fflush(fp);
    }
  }
}

void FixPrintWithZero::setup(int vflag)
{
  if (update->ntimestep==0) {
    end_of_step();
  }
}
