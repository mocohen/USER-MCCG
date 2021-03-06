/* ----------------------------------------------------------------------
Multi-Configurational Coarse-Graining Method 
Voth Group
http://vothgroup.uchicago.edu
Morris Cohen Sharp (C) 2015
mocohen@uchicago.edu


NOTE: YOU CURRENTLY CANNOT HAVE ATOMS OF SAME TYPE ON SAME MOLECULE
------------------------------------------------------------------------- */
//#include <fstream>
//#include <iostream>

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
using namespace PLMD;

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
  Fix(lmp, narg, arg),
  plumed(NULL),
  nlocal(0),
  gatindex(NULL),
  masses(NULL),
  charges(NULL),
  isUmbrellaSampling(0)
{
  setbuf(stdout, NULL);
  printf("Hello world MCCG\n");
  if (narg < 5) error->all(FLERR,"Illegal fix mccg command - fix name mccg controlFile outputFreq");
  
  sscanf(arg[4], "%d", &outputFreq);

  char **computeArgs = new char*[3];
  computeArgs[0] = (char *) "pe_comp_zyx";
  computeArgs[1]= (char *)  "all";
  computeArgs[2]= (char *)  "pe";
    
  //Add compute to force the computation of potential energy of each atom at every time step

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

  compute_pe_atom =  (ComputePEAtom*)(modify->compute[compute_pe_ID]);
  compute_pe_atom->peatomflag = 1;


  readControlFile(arg[3]);



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
  mccg_output.close();
  delete [] cv_array;
  delete [] table_v12;
  delete [] table_f_cv1;
  delete [] table_f_cv2;
  delete [] real_mols;
  delete [] fake_mols;
  delete [] other_mols;
  delete [] corresponding_atom_tags;
  delete [] num_mccg_atoms;
  delete [] v11_list;
  delete [] v22_list; 
  delete [] v12_index;
  delete [] e_vector1;
  delete [] e_vector2;
  delete [] e_value;
  delete [] deltaFs;
  delete [] f_coupling;
  delete [] energy;
  delete plumed;

}

/* ---------------------------------------------------------------------- */

int FixMCCG::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_INTEGRATE;
  mask |= END_OF_STEP;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMCCG::init()
{
  printf("init MCCG!\n");

}

/* ---------------------------------------------------------------------- */

void FixMCCG::setup(int vflag)
{

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

  compute_pe_atom->peatomflag = 1;


  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int numlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  tagint *tag = atom->tag;

  for (int i = 0; i < numlocal; i++)
  {
	  int molid = molecule[i];
	  int atomTag = tag[i]; 
	  f_coupling[i][0] = 0.0;
    f_coupling[i][1] = 0.0;
    f_coupling[i][2] = 0.0;
    if (corresponding_atom_tags[tag[i]-1] != -1){
      for (int j = 0; j < numlocal; j++)
      {

        if(corresponding_atom_tags[atomTag - 1] == tag[j])
        { 
          x[j][0] = x[i][0] + 0.00001;
          x[j][1] = x[i][1] + 0.00001;
          x[j][2] = x[i][2] + 0.00001;
        } 
        }
    }
  }
}


/* ---------------------------------------------------------------------- */

void FixMCCG::post_force(int vflag)
{
  //need to compute per atom energies.
  compute_pe_atom->invoked_flag |= INVOKED_ARRAY;
  modify->addstep_compute(update->ntimestep + 1);
  
  //compute per atom energies:
  compute_peratom();


  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  //nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  tagint *tag = atom->tag;

  //PLUMED Initializations
  /*---------------------------------------------------------------------*/
  int update_gatindex = 0 ;
  if(nlocal!=atom->nlocal){
    if(charges) {
      delete [] charges; 
    }
    if(masses) delete [] masses;
    if(gatindex) delete [] gatindex;
    nlocal=atom->nlocal;
    gatindex=new int [nlocal];
    masses=new double [nlocal];
    charges=new double [nlocal];
    update_gatindex=1;
  } else {
    for(int i=0;i<nlocal;i++){
      if(gatindex[i]!=atom->tag[i]-1){
        update_gatindex=1;
        break;
      }
    }
  }
  if(update_gatindex){
    for(int i=0;i<nlocal;i++){
      gatindex[i]=atom->tag[i]-1;
      masses[i]=atom->mass[atom->type[i]];
      if(atom->q) charges[i]=atom->q[atom->type[i]];
    }
    plumed->cmd("setAtomsNlocal",&nlocal);
    plumed->cmd("setAtomsGatindex",gatindex);
  }
  // set up local virial/box. plumed uses full 3x3 matrices
  double virial[3][3];
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) virial[i][j]=0.0;
  double box[3][3];
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) box[i][j]=0.0;
  box[0][0]=domain->h[0];
  box[1][1]=domain->h[1];
  box[2][2]=domain->h[2];
  box[2][1]=domain->h[3];
  box[2][0]=domain->h[4];
  box[1][0]=domain->h[5];
  
// local variable with timestep:
  int step=update->ntimestep;
  
  double energ = 0.0;
  double *potener = &energ;

// pass all pointers to plumed:
  plumed->cmd("setStep",&step);
  plumed->cmd("setPositions",&atom->x[0][0]);
  plumed->cmd("setBox",&box[0][0]);
  plumed->cmd("setForces",&f_coupling[0][0]);
  plumed->cmd("setMasses",&masses[0]);
  if(atom->q) plumed->cmd("setCharges",&charges[0]);
  plumed->cmd("setVirial",&virial[0][0]);
  plumed->cmd("setEngineArgs", &cv_array[0]);
// do the real calculation:
  plumed->cmd("setEnergy",&potener);  
  plumed->cmd("calc");

// retransform virial to lammps representation:
  Fix::virial[0]=-virial[0][0];
  Fix::virial[1]=-virial[1][1];
  Fix::virial[2]=-virial[2][2];
  Fix::virial[3]=-virial[0][1];
  Fix::virial[4]=-virial[0][2];
  Fix::virial[5]=-virial[1][2];
  /*---------------------------------------------------------------------*/

  nlocal = atom->nlocal;

  double potentialEnergy;
  potentialEnergy = 0.0;
  for (int i = 0; i < sizeof(real_mols)/sizeof(int)-1; i++)
  {
		
		double v11, v22, v12, c1, c2, d1, d2, cv_val, discrim, norm_factor;
		int cv_index;

    //Get PE for both states
		v11 = getEnergy(real_mols[i]);
		v22 = getEnergy(fake_mols[i]) + deltaFs[i];

		if (numCVs == 1)
		{
			cv_index = get_CV_index(cv_array[0]); 
		}
		if (numCVs == 2)
		{
		  cv_index = get_CV_index(cv_array[0], cv_array[1]);
		}
		
		
		
		v12 = table_v12[cv_index];

    // Do not bother doing calculation if coupling is ~ zero
    if(fabs(v12) < 0.00001){
      if(v11 > v22){
	      d1 = 0;
	      d2 = 1;
        e_value[i] = v22;
      }
      else{
        d1 = 1;
        d2 = 0;
        e_value[i] = v11;
      }
    }
    else{
      // Calculate E'val and E'vector
      double trace = v11+v22;
      double determ = (v11*v22) - (v12*v12);
      discrim = (pow(trace,2) / 4.0) - determ;

      e_value[i] = (trace / 2.0 ) - sqrt(discrim);
      d1 = e_value[i] - v22;
      d2 = v12;

    }
    
		
		// Normalize E'vectors and store values

		norm_factor = sqrt(pow(d1,2) + pow(d2,2));
		c1 = d1 / norm_factor;
		c2 = d2 / norm_factor;
		
		v11_list[i] = v11;
		v12_index[i] = cv_index;
		v22_list[i] = v22;
		
		e_vector1[i] = c1;
		e_vector2[i] = c2;
		
		
    potentialEnergy += e_value[i];
    if(step%outputFreq == 0)
    {
      char outString[150];
      sprintf(outString, "%8d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", step, v11, v22, v12, c1, c2, e_value[i], cv_array[0], cv_array[1]);
      mccg_output << outString;
      mccg_output.flush();
    }
  }


  //calculate hellman-feyman forces for 2x2
  for (int atom_ind = 0; atom_ind < nlocal; atom_ind++)
  {
    if ((mask[atom_ind] & groupbit) && (corresponding_atom_tags[tag[atom_ind]-1] != -1))
    {
   		int molind = (molecule[atom_ind] - 1) / 2;

   		int v12_ind;
   		double v11, v22, c1, c2;

      // READ MATRIX ELEMENTS FROM ARRAYS
   		  	   
   		v11 = v11_list[molind];
   		v12_ind = v12_index[molind];
   		v22 = v22_list[molind];
   
   		c1 = e_vector1[molind];
   		c2 = e_vector2[molind];
   		
  		//find other atom
  		int corr_ind;
  		for (int j = 0; j < nlocal; j++)
      {
  			if(corresponding_atom_tags[tag[atom_ind] - 1] == tag[j])
  			{ 
  				corr_ind = j;
  			}
		  }
		
      //GET FORCES ON COUPLING TERM FROM TABLE
  		double dv12_1, dv12_2, dv12;
  		if(numCVs == 1){
  			dv12 = table_f_cv1[v12_ind];
  		}
  		else if(numCVs == 2){
  			dv12_1 = table_f_cv1[v12_ind];
  			dv12_2 = table_f_cv2[v12_ind];
  		}
  		
   		// LOOP THROUGH X,Y,Z COORDINATES
   		for (int dim = 0; dim < 3; dim++)
   		{

   			double term1, term2, term3;
   			 
      	double dTot, dcvdq_1, dcvdq_2, dv12, dcvdq;
      	if(numCVs == 1){
      		dcvdq = f_coupling[atom_ind][dim];
      		dTot = dv12*dcvdq;
      	}
  			else if(numCVs == 2){
  		   	dcvdq_1 = f_coupling[atom_ind][dim];
  		   	dcvdq_2 = f_coupling[corr_ind][dim];
  		  	dTot = (dv12_1*dcvdq_1) + (dv12_2*dcvdq_2);
  			} 
   			term1 = pow(c1, 2)*f[atom_ind][dim];
   			term2 = pow(c2, 2)*f[corr_ind][dim];
   			term3 = 2*c1*c2*dTot;

   			double totForce = term1 + term2 + term3;

        if(isUmbrellaSampling == 1){
          if (numCVs == 1){
            error->all(FLERR,"Cannot do 1 CV");
          }
          else if(numCVs == 2){

            
            double f1 = umbrellaForce1*(cv_array[0] - umbrellaCenter1)*dcvdq_1;
            double f2 = umbrellaForce2*(cv_array[1] - umbrellaCenter2)*dcvdq_2;
            double sumForce = f1 + f2 + totForce;
            totForce += (f1 + f2);
            char outString[200];
            
            sprintf(outString, "%8d %10.3f %10.3f %10.3f %10.3f\n", step, f1, f2, totForce, sumForce);
            mccg_output << outString;
          }
        }

   			f[atom_ind][dim] = totForce;
   			f[corr_ind][dim] = totForce;
   		} 
	  }
  }
  if(step%outputFreq == 0) 
  {
    int nMols = (sizeof(other_mols)/sizeof(*other_mols));
    for(int i = 0; i < nMols; i ++){
      if(other_mols[i] == 1){
        potentialEnergy += getEnergy(i+1);
      }
    }
    printf("MCCG: PotentialEnergy %f\n",potentialEnergy);
  }

}

void FixMCCG::min_post_force(int vflag)
{
  post_force(vflag);
}

int FixMCCG::get_CV_index(double cv1_val)
{
  if (cv1_val < cv1_min){
    mccg_output << "WARNING: CV1 is less than CV1_min\n";
    cv1_val = cv1_min;
  } 
  else if(cv1_val > cv1_max) {
    mccg_output << "WARNING: CV1 is greater than CV1_min\n";
    cv1_val = cv1_max;
  }
  return (cv1_val - cv1_min) / cv1_delta;
}

int FixMCCG::get_CV_index(double cv1_val, double cv2_val)
{
	int cv1_index, cv2_index;
  if (cv1_val < cv1_min){
    mccg_output << "WARNING: CV1 is less than CV1_min\n";
    cv1_val = cv1_min;
  } 
  else if(cv1_val > cv1_max) {
    mccg_output << "WARNING: CV1 is greater than CV1_min\n";
    cv1_val = cv1_max;
  }

  if (cv2_val < cv2_min){
    mccg_output << "WARNING: CV2 is less than CV2_min\n";
    cv2_val = cv2_min;
  } 
  else if(cv2_val > cv2_max) {
    mccg_output << "WARNING: CV2 is greater than CV2_min\n";
    cv2_val = cv2_max;
  }

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
	return sum_energy / 2.0;
}

/* ---------------------------------------------------------------------- */

void FixMCCG::compute_peratom()
{
  int i;



  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set
  */
  int numlocal = atom->nlocal;
  int npair = numlocal;
  int nbond = numlocal;
  int ntotal = numlocal;
  int nkspace = numlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton_bond) nbond += atom->nghost;
  if (force->newton) ntotal += atom->nghost;
  if (force->kspace && force->kspace->tip4pflag) nkspace += atom->nghost;


  // clear local energy array
  //printf("clear local energy array\n");
  for (i = 0; i < ntotal; i++) energy[i] = 0.0;
  //printf("add in per-atom contributions from each force %d %p\n",pairflag,  force->pair);
  // add in per-atom contributions from each force
  /*if (force->pair){
     printf("force pair\n");
  }*/
  //printf("pairs\n");
  
  if (pairflag && force->pair) {
    double *eatom = force->pair->eatom;

    for (i = 0; i < npair; i++) energy[i] += eatom[i];
  }
  
  if (bondflag && force->bond) {
    double *eatom = force->bond->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }
  if (angleflag && force->angle) {
    double *eatom = force->angle->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }
  if (dihedralflag && force->dihedral) {
    double *eatom = force->dihedral->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }
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

  for (i = 0; i < numlocal; i++)
    if (!(mask[i] & groupbit)) energy[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixMCCG::createPlumedObject(int narg, char **arg)
{

  printf("create plumed object\n\n");
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
    printf("energy %f length %f time %f\n",energyUnits,lengthUnits,timeUnits );
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

}

void FixMCCG::readRealMols(char * file)
{
	printf("Begin read mols\n");
	char line[MAXLINE];
	FILE * fnp;
	fnp = fopen(file, "r");
	fgets(line,MAXLINE,fnp); // comment
	fgets(line,MAXLINE,fnp);
	int mccgMols;
  int regMols;
	sscanf(line,"NUM %i OTHER %i",&mccgMols, &regMols);
	printf("\nNum MOls %d\n", mccgMols);
	
	
	//int numMols = 1;
	memory->create(real_mols,   mccgMols, "mccg:real_mols");
	memory->create(fake_mols,   mccgMols, "mccg:fake_mols");
  memory->create(other_mols,  mccgMols*2 + regMols, "mccg:other_mols");
	memory->create(num_mccg_atoms,   mccgMols, "mccg:num_mccg_atoms");
	memory->create(v11_list,   	mccgMols, "mccg:v11_list");
	memory->create(v22_list,   	mccgMols, "mccg:v22_list");
	memory->create(v12_index,   mccgMols, "mccg:v12_list");
	memory->create(e_vector1, 	mccgMols, "mccg:e_vector1");
	memory->create(e_vector2, 	mccgMols, "mccg:e_vector2");
	memory->create(e_value, 	mccgMols, "mccg:e_value");
	memory->create(corresponding_atom_tags,(atom->natoms), "mccg:corresponding_atom_tags");
	memory->create(energy, (atom->natoms), "mccg:energy" );
	memory->create(f_coupling, (atom->natoms), 3, "mccg:f_coupling");
  memory->create(deltaFs, mccgMols, "mccg:deltaFs");
	
	
	for(int i = 0; i < atom->natoms; i++ ){
		corresponding_atom_tags[i] = -1;
	}
	
	printf("Done allocating\n");

  for(int i=0; i<mccgMols; i++) 
  {
  	printf("loop");
   	fgets(line,MAXLINE,fnp);

  	sscanf(line,"MOL %d %d %d %lg", &real_mols[i], &fake_mols[i], &num_mccg_atoms[i], &deltaFs[i]);
  	for (int j=0; j< num_mccg_atoms[i]; j++)
  	{
  		fgets(line,MAXLINE,fnp);
  		int a, b;
  		sscanf(line,"%d %d", &a, &b);
  		//only need to add corr atoms once to prevent double evaluation.
  		if (b > a) corresponding_atom_tags[a-1] = b;

  	}
  }
  printf("done read\n print ids\n");
  for (int i=0; i<atom->natoms; i++) 
  {
  	printf("%d %d\n", i , corresponding_atom_tags[i]);
  }

  int nMols = (sizeof(other_mols)/sizeof(*other_mols));
  for (int i=0; i<nMols; i++)
  {
    other_mols[i] = 1;
    for (int j=0; j<mccgMols; j++){
      if((i + 1) == real_mols[j] || (i + 1) == fake_mols[j])
      {
        other_mols[i] = 0;
      }
    }
  }
  printf("delF %f\n", deltaFs[0]);

}

void FixMCCG::readControlFile(char * file)
{

  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;

  fp = fopen(file, "r");
  if (fp == NULL)
  {
    error->all(FLERR,"Error reading mccg input file\n");
  }

  char couplingTableFile[100];
  char correspondMolsFile[100];
  char plumedFile[100];
  char outFile[100];
  numCVs = 0;
  umbrellaForce1 = umbrellaForce2 = umbrellaCenter1 = umbrellaCenter2 = 0;



  while ((read = getline(&line, &len, fp)) != -1) {
    char inputArg[100];
    char inputParam[100];
    int numArgs = sscanf(line, "%s %s", inputArg, &inputParam);


    if(strcmp(inputArg, "couplingTable") == 0) {

      sscanf(inputParam, "%s", couplingTableFile);

    }
    else if(strcmp(inputArg, "numberCvs") == 0) {
      sscanf(inputParam, "%i", &numCVs);

    }
    else if(strcmp(inputArg, "correspondingMoleculesInput") == 0) {
      sscanf(inputParam, "%s", correspondMolsFile);

    }
    else if(strcmp(inputArg, "plumedInput") == 0) {
      sscanf(inputParam, "%s", plumedFile);

    }    
    else if(strcmp(inputArg, "outputFile") == 0) {
      sscanf(inputParam, "%s", outFile);

    }
    else if(strcmp(inputArg, "umbrellaSampling") == 0) {
      if(strcmp(inputParam, "yes") == 0){
        isUmbrellaSampling = 1;
        printf("MCCG: Umbrella sampling is turned on\n Note: the 1/2 is not included in the force constant\n");
      }
      else{
        printf("MCCG: Umbrella sampling is turned off. Specify \"umbrellaSampling yes\" to use umbrella sampling \n");
      }
    }
    else if(strcmp(inputArg, "umbrellaForce1") == 0) {
      sscanf(inputParam, "%lf", &umbrellaForce1);

    }
    else if(strcmp(inputArg, "umbrellaForce2") == 0) {
      sscanf(inputParam, "%lf", &umbrellaForce2);

    }
        else if(strcmp(inputArg, "umbrellaCenter1") == 0) {
      sscanf(inputParam, "%lf", &umbrellaCenter1);

    }
    else if(strcmp(inputArg, "umbrellaCenter2") == 0) {
      sscanf(inputParam, "%lf", &umbrellaCenter2);

    }

  }
  fclose(fp);
  if (line) {
    free(line);
  }
  printf(" umb %d f1 %f f2 %f c1 %f c2 %f", isUmbrellaSampling, umbrellaForce1, umbrellaForce2, umbrellaCenter1, umbrellaCenter2);
  if (numCVs != 1 && numCVs != 2) {
    error->all(FLERR,"MCCG Error: Error reading input. numCVs should be 1 or 2\n");
  }

  if (outFile != NULL){
    mccg_output.open (outFile);
    mccg_output << "#Timestep       V11        V22        V12      Evec1      Evec2       Eval        CV1        CV2\n";

  }
  else {
    error->all(FLERR,"MCCG Error: Error reading input. Please specify *outputFile* in control file\n");
  }

  if (couplingTableFile != NULL) {
    readCouplingTable(couplingTableFile);    
  }
  else {
    error->all(FLERR,"MCCG Error: Error reading input. Please specify *couplingTable* in control file\n");
  }
  if (correspondMolsFile != NULL) {
    readRealMols(correspondMolsFile);    
  }
  else {
    error->all(FLERR,"MCCG Error: Error reading input. Please specify *correspondingMoleculesInput* in control file\n");
  }

  
  post_integrate();

  if (plumedFile != NULL) { 
    // CREATE PLUMED OBJECT
    char **plumedArgs = new char*[7];
    plumedArgs[1] = (char* )"plumed_zyx";
    plumedArgs[2] = (char*)"plumed";
    plumedArgs[3] = (char*)"plumedfile";
    plumedArgs[4] = (char*) plumedFile;
    plumedArgs[5] = (char*)"outfile";
    plumedArgs[6] = (char*) "cv.out";  /**/

    createPlumedObject(7, plumedArgs);
   
  }
  else {
    error->all(FLERR,"MCCG Error: Error reading input. Please specify *correspondingMoleculesInput* in control file\n");
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
  		memory->create(cv_array, 1, "mccg:cv_array");
	
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
  		memory->create(cv_array, 2, "mccg:cv_array");
	
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
  else{
    error->one(FLERR,"Expecting 1 or 2 CVs");
  }
	
	//while (())
	printf("done reading coupling table\n");

}

void FixMCCG::end_of_step()
{
  // double **f = atom->f;
  // int numlocal = atom->nlocal;
  // for (int atom_ind = 0; atom_ind < numlocal; atom_ind++)
  // {
  //   printf("end: atom i: %d fx: %f fy: %f fz: %f\n" , atom_ind, f[atom_ind][0], f[atom_ind][1], f[atom_ind][2]);
  // }

}


