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
  else error->all(FLERR,"Illegal fix mccg command - fix name mccg v12File ncvs #numCVs mccgParamFile");
  
  
  char **computeArgs;
  memory->create(computeArgs,   4, 15, "mccg:computeArgs");
  strcpy(computeArgs[0], "compute");
  strcpy(computeArgs[1],"pe_comp_zyx");
  strcpy(computeArgs[2], "all");
  strcpy(computeArgs[3], "pe");
  
  readCouplingTable(arg[3]);
  printf("finished reading table\n");
  readRealMols(arg[6]);
  post_integrate();
  
  ComputePEAtom pe_atom(lmp, 4, computeArgs);
  

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
  delete [] v12_index;
  delete [] v22_list; 
  delete [] e_vector1;
  delete [] e_vector2;
  delete [] e_value;
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
  //printf("init MCCG!\n num molecule %d \n", atom->nmolecule);
  //First calculate number of molecules


}

/* ---------------------------------------------------------------------- */

void FixMCCG::setup(int vflag)
{
  printf("setup MCCG\n");
  // WHAT DO I NEED TO DO HERE!!!!
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);

}

/* ---------------------------------------------------------------------- */

void FixMCCG::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
void FixMCCG::post_integrate()
{
  printf("post integrate!\n num molecule %d \n", atom->nmolecule);
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  tagint *tag = atom->tag;

  for (int i = 0; i < nlocal; i++)
  {
  	  int molid = molecule[i];
  	  int atomTag = tag[i]; 
  	  //printf("atomTag %d \n", atomTag);
      for (int j = 0; j < nlocal; j++)
      {
      		//printf("corresponding tag %d \n", tag[j]);
      		if(corresponding_atom_tags[atomTag - 1] == tag[j])
      		{
      			//printf("pos %f %f %f \n", x[i][0], x[i][1], x[i][2]);
      			//offset positions slightly so that vdw don't produce nan
      			x[i][0] = x[j][0] + 0.00001;
      			x[i][1] = x[j][1] + 0.00001;
      			x[i][2] = x[j][2] + 0.00001;
      		}
      
      }
  
  }

}


/* ---------------------------------------------------------------------- */

void FixMCCG::post_force(int vflag)
{
  //need to compute per atom energies.
  printf("mccg post force\n\n\n\n\n");
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  tagint *tag = atom->tag;

  printf("Loop through mccg mols\n");
  for (int i = 0; i < sizeof(real_mols)/sizeof(int)-1; i++)
  {
  		
  		double v11, v22, v12, c1, c2, d1, d2, cv_val, discrim, norm_factor;
  		int cv_index;
  		v11 = -getEnergy(i);
  		v22 = -getEnergy(i + 1);
  		
  		if (numCVs == 1)
  		{
  			cv_index = get_CV_index(1.0); 
  		}
  		if (numCVs == 2)
  		{
			cv_index = get_CV_index(1.0, 2.0);
  		}
  		
  		
  		
  		v12 = table_v12[cv_index];
  		
  		// Calculate Eigenvector for 2x2 matrix for lower E'val
  		discrim = pow(v11, 2) + 4*pow(v12, 2) - (2* v11*v22) + pow(v22, 2);
  		d1 = (v11 -v22 -sqrt(discrim)) / (2*v12);
  		d2 = 1;
  		norm_factor = sqrt(pow(d1,2) + pow(d2,2));
  		c1 = d1 / norm_factor;
  		c2 = d2 / norm_factor;
  		
  		v11_list[i] = v11;
  		v12_index[i] = cv_index;
  		v22_list[i] = v22;
  		
  		e_vector1[i] = c1;
  		e_vector2[i] = c2;
  		
  		
  		double discr = pow(v12,2) + pow((v11 - v22)/2, 2 );
  		e_value[i] = 0.5 * (v11+v22) - sqrt(discr);
  		printf("interation %d through mccg mols\nd1 %f d2 %f c1 %f c2 %f eval %f\n", i, d1, d2, c1, c2, e_value[i]);
  }
  printf("Loop through atoms to change forces\n");
  //calculate hellman-feyman forces for 2x2
  for (int i = 0; i < nlocal; i++)
  {
  	    if ((mask[i] & groupbit) && (corresponding_atom_tags[i] != -1))
  	    {
  	   		int molind = (molecule[i] - 1) / 2;
  	   		int v12_ind;
  	   		double v11, v22, c1, c2;
  	   
  	   		printf("interation %d through atoms\n molid %d tag %d corr %d\n", i, molind, tag[i], corresponding_atom_tags[i]);
  	   
  	   		v11 = v11_list[molind];
  	   		v12_ind = v12_index[molind];
  	   		v22 = v22_list[molind];
  	   
  	   		c1 = e_vector1[molind];
  	   		c2 = e_vector2[molind];
  	   		printf("v11 %f v12 %f v22 %f c1 %f  c2 %f\n", v11, table_v12[v12_ind], v22, c1, c2);
  	   		
  			//find other atom
  			int corr_ind;
  			for (int j = 0; j < nlocal; j++){
  				//printf("tag %d, \n", tag[j]);
  				if(corresponding_atom_tags[i] == tag[j])
  				{ 
  					//printf("corr_ind %d\n", j);
  					corr_ind = j;
  				}
  			}
  			
  			
  	   		for (int j = 0; j < 3; j++)
  	   		{
  	   		
  	   		
  	   			//printf("x y z\n");
  	   			double term1, term2, term3;
  	   			double dcvdq = 1.0;
  	   			double dv12 = 1.0;
  	   			//printf("force %f force \n", f[i][j]);
  	   			term1 = pow(c1, 2)*f[i][j];
  	   			term2 = pow(c2, 2)*f[corr_ind][j];
  	   			term3 = 2*c1*c2*dcvdq*dv12; 
  	   			f[i][j] = term1 + term2 + (2*term3);
  	   			printf("term 1 %f term2 %f term3 %f hf-force %f\n",term1, term2, term3, term1 + term2 + (2*term3));
  	   		}
  	   
  	   }
  
  }

}

void FixMCCG::min_post_force(int vflag)
{
  post_force(vflag);
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

double FixMCCG::getEnergy(int molid)
{
	return rand() % 10;
}

void FixMCCG::readRealMols(char * file)
{
	printf("Begin read mols\n");
	char line[MAXLINE];
	FILE * fnp;
	fnp = fopen(file, "r");
	fgets(line,MAXLINE,fnp); // comment
	fgets(line,MAXLINE,fnp);
	int numMols;
	sscanf(line,"NUM %i",&numMols);
	printf("\nNum MOls %d\n", numMols);
	
	
	//int numMols = 1;
	memory->create(real_mols,   numMols, "mccg:real_mols");
	memory->create(fake_mols,   numMols, "mccg:fake_mols");
	memory->create(num_mccg_atoms,   numMols, "mccg:num_mccg_atoms");
	memory->create(v11_list,   	numMols, "mccg:v11_list");
	memory->create(v22_list,   	numMols, "mccg:v22_list");
	memory->create(v12_index,   numMols, "mccg:v12_list");
	memory->create(e_vector1, 	numMols, "mccg:e_vector1");
	memory->create(e_vector2, 	numMols, "mccg:e_vector2");
	memory->create(e_value, 	numMols, "mccg:e_value");
	memory->create(corresponding_atom_tags,(atom->natoms), "mccg:corresponding_atom_tags");
	
	for(int i = 0; i < atom->natoms; i++ ){
		corresponding_atom_tags[i] = -1;
	}
	
	printf("Done allocating\n");

    for(int i=0; i<numMols; i++) 
    {
    	printf("loop");
     	fgets(line,MAXLINE,fnp);
     	//fake_mols[i] = 2;
     	//real_mols[i] = 1;
     	//num_mccg_atoms = 1;
     	//int mccg_atoms;
    	sscanf(line,"MOL %d %d %d", &real_mols[i], &fake_mols[i], &num_mccg_atoms[i]);
    	for (int j=0; j< num_mccg_atoms[i]; j++)
    	{
    		fgets(line,MAXLINE,fnp);
    		int a, b;
    		sscanf(line,"%d %d", &a, &b);
    		corresponding_atom_tags[a-1] = b;
    		corresponding_atom_tags[b-1] = a;
    	}
    }
    printf("done read\n");


}

void FixMCCG::readCouplingTable(char * file)
{
	printf("read %s\n", file);

	
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
	printf("numCVs %d \n", numCVs);
	if (numCVs == 1)
	{
		printf("cv1 \n");
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
		printf("cv2 \n");
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
	printf("done reading coupling table\n");

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
