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
#include "../Plumed.h"

#include "compute_pe_atom.h"



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
  void min_post_force(int);
  void post_integrate();
  void readCouplingTable(char * file);
  void compute_peratom();
  void createPlumedObject(int narg, char **arg);
  int get_CV_index(double cv1_val);
  int get_CV_index(double cv1_val, double cv2_val);
  void readRealMols(char * file);
  double getEnergy(int molid);
//  double compute_vector(int);
//  double memory_usage();

 protected:
    class ComputePEAtom *compute_pe_atom;

 private:
  //double xvalue,yvalue,zvalue;
  int varflag,iregion;
  int numMolecules;
  int numCVs;
  int table_num_points, cv1_num_points, cv2_num_points;
  double cv1_min, cv1_max, cv1_delta;
  double cv2_min, cv2_max, cv2_delta;
  
  double * cv_array;
  
  double * table_v12;     // Table of coupling
  double * table_f_cv1;   // Table of force from first CV
  double * table_f_cv2;   // Table of force from second CV
  
  int * real_mols;
  int * fake_mols;
  int * corresponding_atom_tags;
  int * num_mccg_atoms;
  
  double * v11_list;
  double * v22_list;
  double * v12_index;
  
  double * e_vector1;
  double * e_vector2;
  double * e_value;
  
  double ** f_coupling;
  
  
  int compute_pe_ID;
  

  
  double *energy;
  int pairflag,bondflag,angleflag,dihedralflag,improperflag,kspaceflag;
  int nmax;
  
  
  
  //PLUMED VARS
  PLMD::Plumed*plumed;
  int nlocal;
// array of atom indexes local to this process:
  int*gatindex;
// array of masses for local atoms:
  double*masses;
// array of charges for local atoms:
  double*charges;
// this is something to enable respa
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
