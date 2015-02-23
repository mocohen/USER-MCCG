/* ----------------------------------------------------------------------
Multi-Configurational Coarse-Graining Method 
Voth Group
http://vothgroup.uchicago.edu
Morris Cohen (C) 2015
mocohen@uchicago.edu


NOTE: YOU CURRENTLY CANNOT HAVE ATOMS OF SAME TYPE ON SAME MOLECULE
------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>

#include "string.h"
#include "stdlib.h"
#include "fix_mccg.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define MAXLINE 1024

/*
Order of Methods being called:
Before MD Begins:
  Fix
  init
  setup
  
For each MD timestep
  post_force

*/

/* ---------------------------------------------------------------------- */

FixMCCG::FixMCCG(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  printf("Hello world MCCG");
  if (narg < 7) error->all(FLERR,"Illegal fix mccg command");
  

  if (strstr(arg[4], "ncvs") == arg[4])
  {
  		int i;
  		sscanf(arg[5], "%d", &numCVs);
  		printf("%d", i);
  		//numCVs = *(int *)arg[5];
  		printf("NumCVs %d \n", numCVs);
  		std::cout << numCVs << std::endl;
  }
  else error->all(FLERR,"Illegal fix mccg command");
  
  
  
  
  readCouplingTable(arg[3]);
  readRealMols(arg[6]);
	/*
  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else if (strcmp(arg[4],"NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else if (strcmp(arg[5],"NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix setforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix setforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix setforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = atom->nmax;
  memory->create(sforce,maxatom,3,"setforce:sforce");
  */
}

/* ---------------------------------------------------------------------- */

FixMCCG::~FixMCCG()
{

  delete [] table_v12;
  delete [] table_f_cv1;
  delete [] table_f_cv2;
  delete [] real_mols;
  delete [] fake_mols;
  delete [] v11_list;
  delete [] v12_list;
  delete [] v22_list; 
  //delete [] xstr;
  //delete [] ystr;
  //delete [] zstr;
  //delete [] idregion;
  //memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixMCCG::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_INTEGRATE;
  //mask |= POST_FORCE_RESPA;
  //mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMCCG::init()
{
  printf("init MCCG!\n num molecule %d \n", atom->nmolecule);
  //First calculate number of molecules

  // check variables
/*
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix setforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix setforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // cannot use non-zero forces for a minimization since no energy is integrated
  // use fix addforce instead

  int flag = 0;
  if (update->whichflag == 2) {
    if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
    if (ystyle == EQUAL || ystyle == ATOM) flag = 1;
    if (zstyle == EQUAL || zstyle == ATOM) flag = 1;
    if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
    if (ystyle == CONSTANT && yvalue != 0.0) flag = 1;
    if (zstyle == CONSTANT && zvalue != 0.0) flag = 1;
  }
  if (flag)
    error->all(FLERR,"Cannot use non-zero forces in an energy minimization");
*/
}

/* ---------------------------------------------------------------------- */

void FixMCCG::setup(int vflag)
{
  printf("setup MCCG\n");
  // WHAT DO I NEED TO DO HERE!!!!
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
/*  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }*/
}

/* ---------------------------------------------------------------------- */

void FixMCCG::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
//THIS REALLY NEEDS TO BE REWRITTEN. THIS IS NOT GOOD CODING.
void FixMCCG::post_integrate()
{
  printf("post integrate!\n num molecule %d \n", atom->nmolecule);
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;

  for (int i = 0; i < nlocal; i++)
  {
  	  int molid = molecule[i];
  	  int atomType = type[i]; 
      for (int j = 0; j < nlocal; j++)
      {
      		if((atomType + num_mccg_atoms == type[j]) && (molid == molecule[j] - 1) )
      		{
      			x[i][0] = x[j][0];
      			x[i][1] = x[j][1];
      			x[i][2] = x[j][2];
      		}
      
      }
  
  }

}


/* ---------------------------------------------------------------------- */

void FixMCCG::post_force(int vflag)
{
  printf("mccg post force\n");
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;

  
  printf("Num local %d\n", nlocal);
  for (int i = 0; i < nlocal; i++)
  {
  ;
  	//printf("coords %f %f %f molind %d molatom %d", x[i][0], x[i][1], x[i][2], molindex[i], molatom[i]);
    printf("index %d coords %f %f %f mask %d molecule %d\n", i, x[i][0], x[i][1], x[i][2], mask[i], molecule[i]);

  }
  
  for (int i = 0; i < len(real_mols); i++)
  {
  		
  }
/*
  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        foriginal[0] += f[i][0];
        foriginal[1] += f[i][1];
        foriginal[2] += f[i][2];
        if (xstyle) f[i][0] = xvalue;
        if (ystyle) f[i][1] = yvalue;
        if (zstyle) f[i][2] = zvalue;
      }

  // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],3,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],3,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        foriginal[0] += f[i][0];
        foriginal[1] += f[i][1];
        foriginal[2] += f[i][2];
        if (xstyle == ATOM) f[i][0] = sforce[i][0];
        else if (xstyle) f[i][0] = xvalue;
        if (ystyle == ATOM) f[i][1] = sforce[i][1];
        else if (ystyle) f[i][1] = yvalue;
        if (zstyle == ATOM) f[i][2] = sforce[i][2];
        else if (zstyle) f[i][2] = zvalue;
      }
  }*/
}

int FixMCCG::get_CV_index(double cv1_val)
{
	return (cv1_val - cv1_min) / cv1_delta;
}

int FixMCCG::get_CV_index(double cv1_val, double cv2_val)
{
	int cv1_index, cv2_index;
	cv1_index = (cv1_val - cv1_min) / cv1_delta;
	cv2_index = (cv2_val - cv2_min) / cv2_delta;
	
	return cv2_num_points * cv1_index + cv2_index;

}

void FixMCCG::readRealMols(char * file)
{
	char line[MAXLINE];
	FILE * fp;
	fp = fopen(file, "r");
	fgets(line,MAXLINE,fp); // comment
	fgets(line,MAXLINE,fp);
	int numMols;
	sscanf(line,"NUM %i",&numMols);
	
	memory->create(real_mols,   numMols / 2, "mccg:real_mols");
	memory->create(fake_mols,   numMols / 2, "mccg:fake_mols");
	memory->create(v11_list,   numMols / 2, "mccg:v11_list");
	memory->create(v22_list,   numMols / 2, "mccg:v22_list");
	memory->create(v12_list,   numMols / 2, "mccg:v12_list");

    for(int i=0; i<table_num_points; i++) 
    {
     	fgets(line,MAXLINE,fp);
    	sscanf(line,"%i %i %i", &real_mols[i], &fake_mols[i], &num_mccg_atoms);
    }


}

void FixMCCG::readCouplingTable(char * file)
{
	printf("read %s\n", file);
	/*std::ifstream inf(file);
	
	if (!inf)
	{
		error->all(FLERR,"Coupling file cannot be opened for reading!");
	}
	while(inf)
	{
		std::string strInput;		
		getline(inf, strInput);
		std::cout << strInput << std::endl;
		
		
		if (strstr(strInput,"NCV") == strInput) 
		{
			std::cout << "NCV" << std::endl;
		}
	
	}*/
	
	
	char line[MAXLINE];
	FILE * fp;
	fp = fopen(file, "r");
	fgets(line,MAXLINE,fp); // comment
	fgets(line,MAXLINE,fp); // Blank line
	
	
	int num_q;
	fgets(line,MAXLINE,fp);
	sscanf(line,"N %i NCV %i",&table_num_points,&num_q);
	//printf("NumCVs %d num_q %d", numCVs, num_q);
	if(numCVs != num_q) error->one(FLERR,"Inconsistent number of CVs.");
	
	if (numCVs == 1)
	{
		int indx;
    	fgets(line,MAXLINE,fp);
    	sscanf(line,"CV %i min %lg max %lg bins %i",&indx,&cv1_min,&cv1_max,&cv1_num_points);
    	if(indx != 1) error->one(FLERR,"Expecting CVs in correct order");

	
		memory->create(table_v12,   table_num_points, "mccg:table_v12");
  		memory->create(table_f_cv1, table_num_points, "mccg:table_f_cv1");
	
		fgets(line,MAXLINE,fp); // Blank line
  
    	// Read tables: outer_loop == CV1 and inner_loop == CV2
    	double cv1, cv2;
    	for(int i=0; i<table_num_points; i++) 
    	{
     		fgets(line,MAXLINE,fp);
    		sscanf(line,"%i %lg %lg %lg",&indx, &cv1, &table_v12[i], &table_f_cv1[i]);
    	}
    	
    	cv1_delta = (cv1_max - cv1_min) / double(cv1_num_points - 1);
	}
	else if(numCVs == 2){
	    int indx;
    	fgets(line,MAXLINE,fp);
    	sscanf(line,"CV %i min %lg max %lg bins %i",&indx,&cv1_min,&cv1_max,&cv1_num_points);
    	if(indx != 1) error->one(FLERR,"Expecting CVs in correct order");

    	fgets(line,MAXLINE,fp);
    	sscanf(line,"CV %i min %lg max %lg bins %i",&indx,&cv2_min,&cv2_max,&cv2_num_points);
    	if(indx != 2) error->one(FLERR,"Expecting CVs in correct order");
	
		memory->create(table_v12,   table_num_points, "mccg:table_v12");
  		memory->create(table_f_cv1, table_num_points, "mccg:table_f_cv1");
  		memory->create(table_f_cv2, table_num_points, "mccg:table_f_cv2");
	
		fgets(line,MAXLINE,fp); // Blank line
  
    	// Read tables: outer_loop == CV1 and inner_loop == CV2
    	double cv1, cv2;
    	for(int i=0; i<table_num_points; i++) 
    	{
     		fgets(line,MAXLINE,fp);
    		sscanf(line,"%i %lg %lg %lg %lg %lg",&indx, &cv1, &cv2, &table_v12[i], &table_f_cv1[i], &table_f_cv2[i]);
    	}
    	
    	cv1_delta = (cv1_max - cv1_min) / double(cv1_num_points - 1);
		cv2_delta = (cv2_max - cv2_min) / double(cv2_num_points - 1);
	
	}
	
	//while (())

}


/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

/*double FixSetForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}*/

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

/*double FixSetForce::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}*/
