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
//#include "group.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "pair.h"
#include "comm.h"

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4



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
  printf("Hello world MCCG\n");
  if (narg < 7) error->all(FLERR,"Illegal fix mccg command");
  

  if (strstr(arg[4], "ncvs") == arg[4])
  {
  		//int i;
  		sscanf(arg[5], "%d", &numCVs);
  		//printf("%d", i);
  		//numCVs = *(int *)arg[5];
  		//printf("NumCVs %d \n", numCVs);
  		std::cout << numCVs << std::endl;
  }
  else error->all(FLERR,"Illegal fix mccg command - fix name mccg v12File ncvs #numCVs mccgParamFile cvFile");
  
  
  char **computeArgs = new char*[3];
  memory->create(computeArgs,   3, 15, "mccg:computeArgs");
  //strcpy(computeArgs[0], "compute");
  
  strcpy(computeArgs[0],"pe_comp_zyx");
  strcpy(computeArgs[1], "all");
  strcpy(computeArgs[2], "pe");
  
fflush(stdout);
  
  //ComputePEAtom pe_atom(lmp, 3, computeArgs);
  //compute_pe_atom = new ComputePEAtom(lmp, 3, computeArgs);
  
  char **newarg = new char*[3];
  newarg[0] = (char *) "pe_comp_zy";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";

  modify->add_compute(3,newarg);
fflush(stdout);  
  compute_pe_ID = modify->find_compute(newarg[0]);
  printf("compute id %d\n", compute_pe_ID);
  delete [] newarg;
  delete [] computeArgs;
  
  Compute * cmop = modify->compute[compute_pe_ID];
  compute_pe_atom =  (ComputePEAtom*)(cmop);
  compute_pe_atom->peatomflag = 1;
  printf("init compte pe atom %p %p\n", compute_pe_atom, cmop);
fflush(stdout);  
  readCouplingTable(arg[3]);
  printf("finished reading table\n");
  readRealMols(arg[6]);
  post_integrate();
    //modify->addstep_compute(update->ntimestep + 1);
  
  /*char **plumedArgs;
  memory->create(plumedArgs,   7, 20, "mccg:plumedArgs");
  strcpy(plumedArgs[1],"plumed_zyx");
  strcpy(plumedArgs[2],"plumed");
  strcpy(plumedArgs[3],"plumedfile");
  strcpy(plumedArgs[4], arg[7]);
  strcpy(plumedArgs[5],"outfile");
  strcpy(plumedArgs[6], "cv.out");  */
  
  //FROM COMPUTE PE ATOM
  peratom_flag = 1;
  size_peratom_cols = 0;
  //peatomflag = 1;
  //timeflag = 1;
  comm_reverse = 1;

  pairflag = 1;
  bondflag = angleflag = dihedralflag = improperflag = 1;
  kspaceflag = 1;

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
  //delete compute_pe_atom;
  delete [] energy;
  //delete plumed;
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
  printf("init MCCG!\n");
  //First calculate number of molecules

  //modify->addstep_compute(update->ntimestep + 1);

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
  //
  compute_pe_atom->peatomflag = 1;
  //compute_pe_atom->invoked_flag |= INVOKED_ARRAY;
  //modify->addstep_compute(update->ntimestep +1);
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
  printf("done post integrate \n");

}


/* ---------------------------------------------------------------------- */

void FixMCCG::post_force(int vflag)
{
  //need to compute per atom energies.
  printf("mccg post force\n\n\n\n\n");
  compute_pe_atom->invoked_flag |= INVOKED_ARRAY;
  modify->addstep_compute(update->ntimestep + 1);
  //compute per atom energies:
  //compute_pe_atom->compute_peratom();
  compute_peratom();
  //printf("compute per atom done");
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
  		v11 = getEnergy(real_mols[i]);
  		v22 = getEnergy(fake_mols[i]);
  		
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
  		printf("atom %d\n", i);
  	 	fflush(stdout);
  	    if ((mask[i] & groupbit) && (corresponding_atom_tags[i] != -1))
  	    {
  	    	printf("corresponding atoms \n");
  	    	fflush(stdout);
  	   		int molind = (molecule[i] - 1) / 2;
  	   	    printf("compute pe atom %f \n", energy[i]);
  	    	fflush(stdout);
  	   		int v12_ind;
  	   		double v11, v22, c1, c2;
  	   		
  	   		printf("interation %d through atoms\n molid %d tag %d corr %d energy %f\n", i, molind, tag[i], corresponding_atom_tags[i], energy[i]);
  	   
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
	double sum_energy = 0;
	tagint *molecule = atom->molecule;
	for (int i = 0; i < atom->nlocal; i++) 
	{
		if (molecule[i] == molid) sum_energy += energy[i];
	} 
	return sum_energy;
}

/* ---------------------------------------------------------------------- */

void FixMCCG::compute_peratom()
{
  printf("Compute per atom \n");
  int i;

  /*invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom energy was not tallied on needed timestep");

  // grow local energy array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(energy);
    nmax = atom->nmax;
    memory->create(energy,nmax,"pe/atom:energy");
    vector_atom = energy;
  }

  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set
  */
  int nlocal = atom->nlocal;
  int npair = nlocal;
  int nbond = nlocal;
  int ntotal = nlocal;
  int nkspace = nlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton_bond) nbond += atom->nghost;
  if (force->newton) ntotal += atom->nghost;
  if (force->kspace && force->kspace->tip4pflag) nkspace += atom->nghost;


  // clear local energy array
  printf("clear local energy array\n");
  for (i = 0; i < ntotal; i++) energy[i] = 0.0;
  printf("add in per-atom contributions from each force %d %p\n",pairflag,  force->pair);
  // add in per-atom contributions from each force
  if (force->pair){
     printf("force pair\n");
  }
  printf("pairs\n");
  
  if (pairflag && force->pair) {
    double *eatom = force->pair->eatom;
    printf("in if statement\n");
    fflush(stdout);
    printf("energy %f\n", energy[0]);
    fflush(stdout);
    printf("eatom %f\n", eatom[0]);
    fflush(stdout);
    for (i = 0; i < npair; i++) energy[i] += eatom[i];
  }
  
  printf("pair energy\n");
  if (bondflag && force->bond) {
    double *eatom = force->bond->eatom;
    printf("eatom %f", eatom[0]);
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }
  printf("pair bond\n");
  if (angleflag && force->angle) {
    double *eatom = force->angle->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }
  printf("pair angle\n");
  if (dihedralflag && force->dihedral) {
    double *eatom = force->dihedral->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }
  printf("pair dihedral\n");
  if (improperflag && force->improper) {
    double *eatom = force->improper->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }

  if (kspaceflag && force->kspace) {
    double *eatom = force->kspace->eatom;
    for (i = 0; i < nkspace; i++) energy[i] += eatom[i];
  }

  // communicate ghost energy between neighbor procs
  //printf("communicate ghost energy between neighbor procs\n");
  //if (force->newton || (force->kspace && force->kspace->tip4pflag))
  //  comm->reverse_comm_compute(this);

  // zero energy of atoms not in group
  // only do this after comm since ghost contributions must be included

  int *mask = atom->mask;

  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) energy[i] = 0.0;
  printf("Done compute");
}

/* ---------------------------------------------------------------------- */

/*void FixMCCG::createPlumedObject(int narg, char **arg)
{
// Not sure this is really necessary:
  if (!atom->tag_enable) error->all(FLERR,"fix plumed requires atom tags");
// Initialize plumed:
  plumed=new PLMD::Plumed;
  plumed->cmd("setMPIComm",&world);

// Set up units
// LAMMPS units wrt kj/mol - nm - ps
// Set up units

  if (force->boltz == 1.0){
// LAMMPS units lj
    plumed->cmd("setNaturalUnits");
  } else {
    double energyUnits=1.0;
    double lengthUnits=1.0;
    double timeUnits=1.0;
    if (force->boltz == 0.0019872067){
// LAMMPS units real :: kcal/mol; angstrom; fs
      energyUnits=4.184;
      lengthUnits=0.1;
      timeUnits=0.001;
    } else if (force->boltz == 8.617343e-5){
// LAMMPS units metal :: eV; angstrom; ps
      energyUnits=96.48530749925792;
      lengthUnits=0.1;
      timeUnits=1.0;
    } else if (force->boltz == 1.3806504e-23){
// LAMMPS units si :: Joule, m; s
      energyUnits=0.001;
      lengthUnits=1.e-9;
      timeUnits=1.e-12;
    } else if (force->boltz == 1.3806504e-16){
// LAMMPS units cgs :: erg; cms;, s
      energyUnits=6.0221418e13;
      lengthUnits=1.e-7;
      timeUnits=1.e-12;
    } else if (force->boltz == 3.16681534e-6){
// LAMMPS units electron :: Hartree, bohr, fs
      energyUnits=2625.5257;
      lengthUnits=0.052917725;
      timeUnits=0.001;
    } else error->all(FLERR,"Odd LAMMPS units, plumed cannot work with that");
    plumed->cmd("setMDEnergyUnits",&energyUnits);
    plumed->cmd("setMDLengthUnits",&lengthUnits);
    plumed->cmd("setMDTimeUnits",&timeUnits);
  }

// Read fix parameters:
  int next=0;
  for(int i=3;i<narg;++i){
    if(!strcmp(arg[i],"outfile")) next=1;
    else if(next==1){
      plumed->cmd("setLogFile",arg[i]);
      next=0;
    }
    else if(!strcmp(arg[i],"plumedfile"))next=2;
    else if(next==2){
      plumed->cmd("setPlumedDat",arg[i]);
      next=0;
    }
    else error->all(FLERR,"syntax error in fix plumed - use 'fix name plumed plumedfile plumed.dat outfile plumed.out' ");
  }
  if(next==1) error->all(FLERR,"missing argument for outfile option");
  if(next==2) error->all(FLERR,"missing argument for plumedfile option");

  plumed->cmd("setMDEngine","LAMMPS");

  int natoms=int(atom->natoms);
  plumed->cmd("setNatoms",&natoms);

  double dt=update->dt;
  plumed->cmd("setTimestep",&dt);

  virial_flag=1;

// This is the real initialization:
  plumed->cmd("init");

}*/

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
	memory->create(energy, numMols, "mccg:energy" );
	
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
    printf("done read\n print ids\n");
    for (int i=0; i<atom->natoms; i++) 
    {
    	printf("%d %d\n", i , corresponding_atom_tags[i]);
    }


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
