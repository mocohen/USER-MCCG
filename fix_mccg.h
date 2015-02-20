/* -*- c++ -*- ----------------------------------------------------------
Multi-Configurational Coarse-Graining Method 
Voth Group
http://vothgroup.uchicago.edu
Morris Cohen (C) 2015
mocohen@uchicago.edu

------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(mccg,FixMCCG)

#else

#ifndef LMP_FIX_MCCG_H
#define LMP_FIX_MCCG_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMCCG : public Fix {
 public:
  FixMCCG(class LAMMPS *, int, char **);
  ~FixMCCG();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
//  double compute_vector(int);
//  double memory_usage();

 private:
  //double xvalue,yvalue,zvalue;
  int varflag,iregion;
  int numMolecules;
  //char *xstr,*ystr,*zstr;
  //char *idregion;
  //int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  //double foriginal[3],foriginal_all[3];
  //int force_flag;
//  int nlevels_respa;

  int maxatom;
  //double **sforce;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix setforce does not exist

Self-explanatory.

E: Variable name for fix setforce does not exist

Self-explanatory.

E: Variable for fix setforce is invalid style

Only equal-style variables can be used.

E: Cannot use non-zero forces in an energy minimization

Fix setforce cannot be used in this manner.  Use fix addforce
instead.

*/
