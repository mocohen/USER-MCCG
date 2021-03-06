/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "wrapper/Plumed.h"
#include "tools/Communicator.h"
#include "tools/Random.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include "tools/Units.h"
#include "tools/PDB.h"

// when using molfile plugin
#ifdef __PLUMED_HAS_MOLFILE
#ifdef __PLUMED_INTERNAL_MOLFILE_PLUGINS
/* Use the internal ones. Alternatively:
 *    ifneq (,$(findstring __PLUMED_INTERNAL_MOLFILE_PLUGINS,$(CPPFLAGS)))
 *    CPPFLAGS+=-I../molfile  
 */
#include "molfile/libmolfile_plugin.h"
#include "molfile/molfile_plugin.h"
using namespace PLMD::molfile;
#else
#include <libmolfile_plugin.h>
#include <molfile_plugin.h>
#endif
#endif

using namespace std;

namespace PLMD {
namespace cltools{

//+PLUMEDOC TOOLS driver
/*
driver is a tool that allows one to to use plumed to post-process an existing trajectory.

The input to driver is specified using the command line arguments described below.

In addition, you can use the special \subpage READ command inside your plumed input
to read in colvar files that were generated during your MD simulation.  The values
read in can then be treated like calculated colvars. 

\par Examples

The following command tells plumed to postprocess the trajectory contained in trajectory.xyz
 by performing the actions described in the input file plumed.dat.  If an action that takes the
stride keyword is given a stride equal to \f$n\f$ then it will be performed only on every \f$n\f$th
frame in the trajectory file.
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz
\endverbatim

The following command tells plumed to postprocess the trajectory contained in trajectory.xyz.
 by performing the actions described in the input file plumed.dat. Here though
--trajectory-stride is set equal to the frequency with which frames were output during the trajectory
and the --timestep is equal to the simulation timestep.  As such the STRIDE parameters in the plumed.dat
files are no longer ignored and any files output resemble those that would have been generated
had we run the calculation we are running with driver when the MD simulation was running.
\verbatim
plumed driver --plumed plumed.dat --ixyz trajectory.xyz --trajectory-stride 100 --timestep 0.001
\endverbatim

By default you have access to a subset of the trajectory file formats
supported by VMD, e.g. xtc and dcd:

\verbatim
plumed driver --plumed plumed.dat --pdb diala.pdb --mf_xtc traj.xtc --trajectory-stride 100 --timestep 0.001
\endverbatim

where --mf_ prefixes the extension of one of the accepted molfile
plugin format.

To have support of all of VMD's plugins you need to recompile
PLUMED. You need to download the SOURCE of VMD, which contains
a plugins directory. Adapt build.sh and compile it. At
the end, you should get the molfile plugins compiled as a static
library libmolfile_plugin.a. Locate said file and libmolfile_plugin.h, 
and customize the configure command with something along
the lines of:

\verbatim
configure [...] LDFLAGS="-ltcl8.5 -L/mypathtomolfilelibrary/ -L/mypathtotcl" CPPFLAGS="-I/mypathtolibmolfile_plugin.h/"
\endverbatim

and rebuild. Check the available molfile plugins and limitations at http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/.

*/
//+ENDPLUMEDOC
//

#ifdef __PLUMED_HAS_MOLFILE
static vector<molfile_plugin_t *> plugins;
static map <string, unsigned> pluginmap;
static int register_cb(void *v, vmdplugin_t *p){
  //const char *key = p->name;
  std::pair<std::map<string,unsigned>::iterator,bool> ret; 
  ret = pluginmap.insert ( std::pair<string,unsigned>(string(p->name),plugins.size()) );
  if (ret.second==false) { 
	//cerr<<"MOLFILE: found duplicate plugin for "<<key<<" : not inserted "<<endl; 
  }else{
	//cerr<<"MOLFILE: loading plugin "<<key<<" number "<<plugins.size()-1<<endl;
  	plugins.push_back((molfile_plugin_t *)p);
  } 
  return VMDPLUGIN_SUCCESS;
}
#endif

template<typename real>
class Driver : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  Driver(const CLToolOptions& co );
  int main(FILE* in,FILE*out,Communicator& pc);
  string description()const;
};

template<typename real>
void Driver<real>::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys ); keys.isDriver();
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("compulsory","--plumed","plumed.dat","specify the name of the plumed input file");
  keys.add("compulsory","--timestep","1.0","the timestep that was used in the calculation that produced this trajectory in picoseconds");
  keys.add("compulsory","--trajectory-stride","1","the frequency with which frames were output to this trajectory during the simulation");
  keys.add("compulsory","--multi","0","set number of replicas for multi environment (needs mpi)");
  keys.addFlag("--noatoms",false,"don't read in a trajectory.  Just use colvar files as specified in plumed.dat");
  keys.add("atoms","--ixyz","the trajectory in xyz format");
  keys.add("atoms","--igro","the trajectory in gro format");
  keys.add("optional","--length-units","units for length, either as a string or a number");
  keys.add("optional","--dump-forces","dump the forces on a file");
  keys.add("optional","--dump-forces-fmt","( default=%%f ) the format to use to dump the forces");
  keys.add("optional","--pdb","provides a pdb with masses and charges");
  keys.add("optional","--box","comma-separated box dimensions (3 for orthorombic, 9 for generic)");
  keys.add("hidden","--debug-float","turns on the single precision version (to check float interface)");
  keys.add("hidden","--debug-dd","use a fake domain decomposition");
  keys.add("hidden","--debug-pd","use a fake particle decomposition");
  keys.add("hidden","--debug-grex","use a fake gromacs-like replica exchange, specify exchange stride");
  keys.add("hidden","--debug-grex-log","log file for debug=grex");
#ifdef __PLUMED_HAS_MOLFILE
  MOLFILE_INIT_ALL
  MOLFILE_REGISTER_ALL(NULL, register_cb)
  for(int i=0;i<plugins.size();i++){
	string kk="--mf_"+string(plugins[i]->name);
	string mm=" molfile: the trajectory in "+string(plugins[i]->name)+" format " ;
	//cerr<<"REGISTERING "<<kk<<mm<<endl;
  	keys.add("atoms",kk,mm);
  }
#endif
}
template<typename real>
Driver<real>::Driver(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}
template<typename real>
string Driver<real>::description()const{ return "analyze trajectories with plumed"; }

template<typename real>
int Driver<real>::main(FILE* in,FILE*out,Communicator& pc){

  Units units;
  PDB pdb;

// Parse everything
  bool printhelpdebug; parseFlag("--help-debug",printhelpdebug);
  if( printhelpdebug ){
      fprintf(out,"%s",
         "Additional options for debug (only to be used in regtest):\n"
         "  [--debug-float]         : turns on the single precision version (to check float interface)\n"
         "  [--debug-dd]            : use a fake domain decomposition\n"
         "  [--debug-pd]            : use a fake particle decomposition\n"
      );
      return 0;
  }
  // Are we reading trajectory data
  bool noatoms; parseFlag("--noatoms",noatoms);

  std::string fakein; 
  bool debugfloat=parse("--debug-float",fakein);
  if(debugfloat && sizeof(real)!=sizeof(float)){
      CLTool* cl=cltoolRegister().create(CLToolOptions("driver-float"));    //new Driver<float>(*this);
      cl->setInputData(this->getInputData());
      int ret=cl->main(in,out,pc);
      delete cl;
      return ret;
  }

  bool debug_pd=parse("--debug-pd",fakein);
  bool debug_dd=parse("--debug-dd",fakein);
  if( debug_pd || debug_dd ){
    if(noatoms) error("cannot debug without atoms");
  }

// set up for multi replica driver:
  int multi=0;
  parse("--multi",multi);
  Communicator intracomm;
  Communicator intercomm;
  if(multi){
    int ntot=pc.Get_size();
    int nintra=ntot/multi;
    if(multi*nintra!=ntot) error("invalid number of processes for multi environment");
    pc.Split(pc.Get_rank()/nintra,pc.Get_rank(),intracomm);
    pc.Split(pc.Get_rank()%nintra,pc.Get_rank(),intercomm);
  } else {
    intracomm.Set_comm(pc.Get_comm());
  }

// set up for debug replica exchange:
  bool debug_grex=parse("--debug-grex",fakein);
  int  grex_stride=0;
  FILE*grex_log=NULL;
  if(debug_grex){
    if(noatoms) error("must have atoms to debug_grex");
    if(multi<2)  error("--debug_grex needs --multi with at least two replicas");
    Tools::convert(fakein,grex_stride);
    string n; Tools::convert(intercomm.Get_rank(),n);
    string file;
    parse("--debug-grex-log",file);
    if(file.length()>0){
      file+="."+n;
      grex_log=fopen(file.c_str(),"w");
    }
  }

// Read the plumed input file name  
  string plumedFile; parse("--plumed",plumedFile);
// the timestep
  double t; parse("--timestep",t);
  real timestep=real(t);
// the stride
  unsigned stride; parse("--trajectory-stride",stride);
// are we writing forces
  string dumpforces(""), dumpforcesFmt("%f");; 
  if(!noatoms) parse("--dump-forces",dumpforces);
  if(dumpforces!="") parse("--dump-forces-fmt",dumpforcesFmt);

  string trajectory_fmt;

  bool use_molfile=false; 
#ifdef __PLUMED_HAS_MOLFILE
  molfile_plugin_t *api=NULL;      
  void *h_in=NULL;
  molfile_timestep_t ts_in; // this is the structure that has the timestep 
  ts_in.coords=NULL;
#endif

// Read in an xyz file
  string trajectoryFile(""), pdbfile("");
  bool pbc_cli_given=false; vector<double> pbc_cli_box(9,0.0);
  if(!noatoms){
     std::string traj_xyz; parse("--ixyz",traj_xyz);
     std::string traj_gro; parse("--igro",traj_gro);
#ifdef __PLUMED_HAS_MOLFILE 
     for(int i=0;i<plugins.size();i++){ 	
	string molfile_key="--mf_"+string(plugins[i]->name);
	string traj_molfile;
	parse(molfile_key,traj_molfile);
	if(traj_molfile.length()>0){
		fprintf(out,"\nDRIVER: Found molfile format trajectory %s with name %s\n",plugins[i]->name,traj_molfile.c_str());
		trajectoryFile=traj_molfile;
		trajectory_fmt=string(plugins[i]->name);
		use_molfile=true;
		api = plugins[i];
	}
     } 
#endif
     if(traj_xyz.length()>0 && traj_gro.length()>0){
       fprintf(stderr,"ERROR: cannot provide more than one trajectory file\n");
       if(grex_log)fclose(grex_log);
       return 1;
     }
     if(traj_xyz.length()>0 && trajectoryFile.length()==0){
       trajectoryFile=traj_xyz;
       trajectory_fmt="xyz";
     }
     if(traj_gro.length()>0 && trajectoryFile.length()==0){
       trajectoryFile=traj_gro;
       trajectory_fmt="gro";
     }
     if(trajectoryFile.length()==0){
       fprintf(stderr,"ERROR: missing trajectory data\n"); 
       if(grex_log)fclose(grex_log);
       return 1;
     }
     string lengthUnits(""); parse("--length-units",lengthUnits);
     if(lengthUnits.length()>0) units.setLength(lengthUnits);
  
     parse("--pdb",pdbfile);
     if(pdbfile.length()>0){
       bool check=pdb.read(pdbfile,false,1.0);
       if(!check) error("error reading pdb file");
     }

     string pbc_cli_list; parse("--box",pbc_cli_list);
     if(pbc_cli_list.length()>0) {
       pbc_cli_given=true;
       vector<string> words=Tools::getWords(pbc_cli_list,",");
       if(words.size()==3){
         for(int i=0;i<3;i++) sscanf(words[i].c_str(),"%100lf",&(pbc_cli_box[4*i]));
       } else if(words.size()==9) {
         for(int i=0;i<9;i++) sscanf(words[i].c_str(),"%100lf",&(pbc_cli_box[i]));
       } else {
         string msg="ERROR: cannot parse command-line box "+pbc_cli_list;
         fprintf(stderr,"%s\n",msg.c_str());
         return 1;
       }

     }
  }

  if( debug_dd && debug_pd ) error("cannot use debug-dd and debug-pd at the same time");
  if(debug_pd || debug_dd){
    if( !Communicator::initialized() ) error("needs mpi for debug-pd");
  }

  Plumed p;
  int rr=sizeof(real);
  p.cmd("setRealPrecision",&rr);
  int checknatoms=-1;
  int step=0;
  if(Communicator::initialized()){
    if(multi){
      if(intracomm.Get_rank()==0) p.cmd("GREX setMPIIntercomm",&intercomm.Get_comm());
      p.cmd("GREX setMPIIntracomm",&intracomm.Get_comm());
      p.cmd("GREX init");
    } 
    p.cmd("setMPIComm",&intracomm.Get_comm());
  } 
  p.cmd("setMDLengthUnits",&units.getLength());
  p.cmd("setMDEngine","driver");
  p.cmd("setTimestep",&timestep);
  p.cmd("setPlumedDat",plumedFile.c_str());
  p.cmd("setLog",out);

  if(multi){
    string n;
    Tools::convert(intercomm.Get_rank(),n);
    trajectoryFile+="."+n;
  }


  int natoms;

  FILE* fp=NULL; FILE* fp_forces=NULL;
  if(!noatoms){
     if (trajectoryFile=="-") 
       fp=in;
     else {
       if(use_molfile==true){
#ifdef __PLUMED_HAS_MOLFILE
        h_in = api->open_file_read(trajectoryFile.c_str(), trajectory_fmt.c_str(), &natoms);
        ts_in.coords = new float [3*natoms];
#endif
       }else{
         fp=fopen(trajectoryFile.c_str(),"r");
         if(!fp){
           string msg="ERROR: Error opening trajectory file "+trajectoryFile;
           fprintf(stderr,"%s\n",msg.c_str());
           return 1;
         }
       }
     }
     if(dumpforces.length()>0){
       if(Communicator::initialized() && pc.Get_size()>1){
         string n;
         Tools::convert(pc.Get_rank(),n);
         dumpforces+="."+n;
       }
       fp_forces=fopen(dumpforces.c_str(),"w");
     }
  }

  std::string line;
  std::vector<real> coordinates;
  std::vector<real> forces;
  std::vector<real> masses;
  std::vector<real> charges;
  std::vector<real> cell;
  std::vector<real> virial;

// variables to test particle decomposition
  int pd_nlocal;
  int pd_start;
// variables to test random decomposition (=domain decomposition)
  std::vector<int>  dd_gatindex;
  std::vector<int>  dd_g2l;
  std::vector<real> dd_masses;
  std::vector<real> dd_charges;
  std::vector<real> dd_forces;
  std::vector<real> dd_coordinates;
  int dd_nlocal;
// random stream to choose decompositions
  Random rnd;

  while(true){
    if(!noatoms){
       if(use_molfile==true){	
#ifdef __PLUMED_HAS_MOLFILE
          int rc; 
    	  rc = api->read_next_timestep(h_in, natoms, &ts_in);
          //if(rc==MOLFILE_SUCCESS){
	  //       printf(" read this one :success \n");
	  //}
	  if(rc==MOLFILE_EOF){
	         //printf(" read this one :eof or error \n");
	         break;
	  }
#endif
       }else{ 
         if(!Tools::getline(fp,line)) break;
       }
    }

    bool first_step=false;
    if(!noatoms){
      if(use_molfile==false){
        if(trajectory_fmt=="gro") if(!Tools::getline(fp,line)) error("premature end of trajectory file");
        sscanf(line.c_str(),"%100d",&natoms);
      }
    }
    if(checknatoms<0 && !noatoms){
      pd_nlocal=natoms;
      pd_start=0;
      first_step=true;
      masses.assign(natoms,real(1.0));
      charges.assign(natoms,real(0.0));
//case pdb: structure
      if(pdbfile.length()>0){
        for(unsigned i=0;i<pdb.size();++i){
          AtomNumber an=pdb.getAtomNumbers()[i];
          unsigned index=an.index();
          if( index>=unsigned(natoms) ) error("atom index in pdb exceeds the number of atoms in trajectory");
          masses[index]=pdb.getOccupancy()[i];
          charges[index]=pdb.getBeta()[i];
        }
      }
    } else if( checknatoms<0 && noatoms ){ 
      natoms=0; 
    }
    if( checknatoms<0 ){
      checknatoms=natoms;
      p.cmd("setNatoms",&natoms);
      p.cmd("init");
    }
    if(checknatoms!=natoms){
       std::string stepstr; Tools::convert(step,stepstr);
       error("number of atoms in frame " + stepstr + " does not match number of atoms in first frame");
    }

    coordinates.assign(3*natoms,real(0.0));
    forces.assign(3*natoms,real(0.0));
    cell.assign(9,real(0.0));
    virial.assign(9,real(0.0));

    if( first_step || rnd.U01()>0.5){
      if(debug_pd){
        int npe=intracomm.Get_size();
        vector<int> loc(npe,0);
        vector<int> start(npe,0);
        for(int i=0;i<npe-1;i++){
          int cc=(natoms*2*rnd.U01())/npe;
          if(start[i]+cc>natoms) cc=natoms-start[i];
          loc[i]=cc;
          start[i+1]=start[i]+loc[i];
        }
        loc[npe-1]=natoms-start[npe-1];
        intracomm.Bcast(loc,0);
        intracomm.Bcast(start,0);
        pd_nlocal=loc[intracomm.Get_rank()];
        pd_start=start[intracomm.Get_rank()];
        if(intracomm.Get_rank()==0){
          fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
          fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",loc[i]); printf("\n");
          fprintf(out,"DRIVER: "); for(int i=0;i<npe;i++) fprintf(out,"%d ",start[i]); printf("\n");
        }
        p.cmd("setAtomsNlocal",&pd_nlocal);
        p.cmd("setAtomsContiguous",&pd_start);
      } else if(debug_dd){
        int npe=intracomm.Get_size();
        int rank=intracomm.Get_rank();
        dd_charges.assign(natoms,0.0);
        dd_masses.assign(natoms,0.0);
        dd_gatindex.assign(natoms,-1);
        dd_g2l.assign(natoms,-1);
        dd_coordinates.assign(3*natoms,0.0);
        dd_forces.assign(3*natoms,0.0);
        dd_nlocal=0;
        for(int i=0;i<natoms;++i){
          double r=rnd.U01()*npe;
          int n; for(n=0;n<npe;n++) if(n+1>r)break;
          plumed_assert(n<npe);
          if(n==rank){
            dd_gatindex[dd_nlocal]=i;
            dd_g2l[i]=dd_nlocal;
            dd_charges[dd_nlocal]=charges[i];
            dd_masses[dd_nlocal]=masses[i];
            dd_nlocal++;
          }
        }
        if(intracomm.Get_rank()==0){
          fprintf(out,"\nDRIVER: Reassigning particle decomposition\n");
        }
        p.cmd("setAtomsNlocal",&dd_nlocal);
        p.cmd("setAtomsGatindex",&dd_gatindex[0]);
      }
    }

    int plumedStopCondition=0;
    p.cmd("setStep",&step);
    p.cmd("setStopFlag",&plumedStopCondition);
    if(!noatoms){
       if(use_molfile){
#ifdef __PLUMED_HAS_MOLFILE
    	   if(pbc_cli_given==false) {
    		   // info on the cell: convert using pbcset.tcl from pbctools in vmd distribution
    		   real cosBC=cos(ts_in.alpha*pi/180.);
    		   //double sinBC=sin(ts_in.alpha*pi/180.);
    		   real cosAC=cos(ts_in.beta*pi/180.);
    		   real cosAB=cos(ts_in.gamma*pi/180.);
    		   real sinAB=sin(ts_in.gamma*pi/180.);
    		   real Ax=ts_in.A;
    		   real Bx=ts_in.B*cosAB;
    		   real By=ts_in.B*sinAB;
                   real Cx=ts_in.C*cosAC;
                   real Cy=(ts_in.C*ts_in.B*cosBC-Cx*Bx)/By;
                   real Cz=sqrt(ts_in.C*ts_in.C-Cx*Cx-Cy*Cy);
    		   cell[0]=Ax/10.;cell[1]=0.;cell[2]=0.;
    		   cell[3]=Bx/10.;cell[4]=By/10.;cell[5]=0.;
    		   cell[6]=Cx/10.;cell[7]=Cy/10.;cell[8]=Cz/10.;
    		   //cerr<<"CELL "<<cell[0]<<" "<<cell[1]<<" "<<cell[2]<<" "<<cell[3]<<" "<<cell[4]<<" "<<cell[5]<<" "<<cell[6]<<" "<<cell[7]<<" "<<cell[8]<<endl;
    	   }else{
    		   for(unsigned i=0;i<9;i++)cell[i]=pbc_cli_box[i];
    	   }
    	   // info on coords
    	   // the order is xyzxyz...
    	   for(unsigned i=0;i<3*natoms;i++){
    		   coordinates[i]=real(ts_in.coords[i]/10.); //convert to nm
    		   //cerr<<"COOR "<<coordinates[i]<<endl;
    	   }
#endif
       }else{
       if(trajectory_fmt=="xyz"){
         if(!Tools::getline(fp,line)) error("premature end of trajectory file");

         std::vector<double> celld(9,0.0);
         if(pbc_cli_given==false) {
           std::vector<std::string> words;
           words=Tools::getWords(line);
           if(words.size()==3){
             sscanf(line.c_str(),"%100lf %100lf %100lf",&celld[0],&celld[4],&celld[8]);
           } else if(words.size()==9){
             sscanf(line.c_str(),"%100lf %100lf %100lf %100lf %100lf %100lf %100lf %100lf %100lf",
                    &celld[0], &celld[1], &celld[2],
                    &celld[3], &celld[4], &celld[5],
                    &celld[6], &celld[7], &celld[8]);
           } else error("needed box in second line of xyz file");
         } else {			// from command line
           celld=pbc_cli_box;
         }
         for(unsigned i=0;i<9;i++)cell[i]=real(celld[i]);
       }
  	   int ddist=0;
       // Read coordinates
       for(int i=0;i<natoms;i++){
         bool ok=Tools::getline(fp,line);
         if(!ok) error("premature end of trajectory file");
         double cc[3];
         if(trajectory_fmt=="xyz"){
           char dummy[1000];
           int ret=std::sscanf(line.c_str(),"%999s %100lf %100lf %100lf",dummy,&cc[0],&cc[1],&cc[2]);
           if(ret!=4) error("cannot read line"+line);
         } else if(trajectory_fmt=="gro"){
           // do the gromacs way
           if(!i){
        	   //
        	   // calculate the distance between dots (as in gromacs gmxlib/confio.c, routine get_w_conf )
        	   //
        	   const char      *p1, *p2, *p3;
        	   p1 = strchr(line.c_str(), '.');
        	   if (p1 == NULL) error("seems there are no coordinates in the gro file");
        	   p2 = strchr(&p1[1], '.');
        	   if (p2 == NULL) error("seems there is only one coordinates in the gro file");
        	   ddist = p2 - p1;
        	   p3 = strchr(&p2[1], '.');
        	   if (p3 == NULL)error("seems there are only two coordinates in the gro file");
        	   if (p3 - p2 != ddist)error("not uniform spacing in fields in the gro file");
           }
           Tools::convert(line.substr(20,ddist),cc[0]);
           Tools::convert(line.substr(20+ddist,ddist),cc[1]);
           Tools::convert(line.substr(20+ddist+ddist,ddist),cc[2]);
         } else plumed_error();
         if(!debug_pd || ( i>=pd_start && i<pd_start+pd_nlocal) ){
           coordinates[3*i]=real(cc[0]);
           coordinates[3*i+1]=real(cc[1]);
           coordinates[3*i+2]=real(cc[2]);
         }
       }
       if(trajectory_fmt=="gro"){
         if(!Tools::getline(fp,line)) error("premature end of trajectory file");
         std::vector<string> words=Tools::getWords(line);
         if(words.size()<3) error("cannot understand box format");
         Tools::convert(words[0],cell[0]);
         Tools::convert(words[1],cell[4]);
         Tools::convert(words[2],cell[8]);
         if(words.size()>3) Tools::convert(words[3],cell[1]);
         if(words.size()>4) Tools::convert(words[4],cell[2]);
         if(words.size()>5) Tools::convert(words[5],cell[3]);
         if(words.size()>6) Tools::convert(words[6],cell[5]);
         if(words.size()>7) Tools::convert(words[7],cell[6]);
         if(words.size()>8) Tools::convert(words[8],cell[7]);
       }

     }

       if(debug_dd){
         for(int i=0;i<dd_nlocal;++i){
           int kk=dd_gatindex[i];
           dd_coordinates[3*i+0]=coordinates[3*kk+0];
           dd_coordinates[3*i+1]=coordinates[3*kk+1];
           dd_coordinates[3*i+2]=coordinates[3*kk+2];
         }
         p.cmd("setForces",&dd_forces[0]);
         p.cmd("setPositions",&dd_coordinates[0]);
         p.cmd("setMasses",&dd_masses[0]);
         p.cmd("setCharges",&dd_charges[0]);
       } else {
         p.cmd("setForces",&forces[3*pd_start]);
         p.cmd("setPositions",&coordinates[3*pd_start]);
         p.cmd("setMasses",&masses[pd_start]);
         p.cmd("setCharges",&charges[pd_start]);
       }
       p.cmd("setBox",&cell[0]);
       p.cmd("setVirial",&virial[0]);
   }
   p.cmd("calc");

// this is necessary as only processor zero is adding to the virial:
   intracomm.Bcast(virial,0);
   if(debug_pd) intracomm.Sum(forces);
   if(debug_dd){
     for(int i=0;i<dd_nlocal;i++){
       forces[3*dd_gatindex[i]+0]=dd_forces[3*i+0];
       forces[3*dd_gatindex[i]+1]=dd_forces[3*i+1];
       forces[3*dd_gatindex[i]+2]=dd_forces[3*i+2];
     }
     dd_forces.assign(3*natoms,0.0);
     intracomm.Sum(forces);
   }
   if(debug_grex &&step%grex_stride==0){
     p.cmd("GREX savePositions");
     if(intracomm.Get_rank()>0){
       p.cmd("GREX prepare");
     } else {
       int r=intercomm.Get_rank();
       int n=intercomm.Get_size();
       int partner=r+(2*((r+step/grex_stride)%2))-1;
       if(partner<0)partner=0;
       if(partner>=n) partner=n-1;
       p.cmd("GREX setPartner",&partner);
       p.cmd("GREX calculate");
       p.cmd("GREX shareAllDeltaBias");
       for(int i=0;i<n;i++){
         string s; Tools::convert(i,s);
         real a; s="GREX getDeltaBias "+s; p.cmd(s.c_str(),&a);
         if(grex_log) fprintf(grex_log," %f",a);
       }
       if(grex_log) fprintf(grex_log,"\n");
     }
   }


   if(fp_forces){
     fprintf(fp_forces,"%d\n",natoms);
     string fmt=dumpforcesFmt+" "+dumpforcesFmt+" "+dumpforcesFmt+"\n";
     fprintf(fp_forces,fmt.c_str(),virial[0],virial[4],virial[8]);
     fmt="X "+fmt;
     for(int i=0;i<natoms;i++)
       fprintf(fp_forces,fmt.c_str(),forces[3*i],forces[3*i+1],forces[3*i+2]);
   }

    if(noatoms && plumedStopCondition) break;

    step+=stride;
  }
  p.cmd("runFinalJobs");

  if(fp_forces) fclose(fp_forces);
  if(fp && fp!=in)fclose(fp);
#ifdef __PLUMED_HAS_MOLFILE
  if(h_in) api->close_file_read(h_in);
  if(ts_in.coords) delete [] ts_in.coords;
#endif
  if(grex_log) fclose(grex_log);

  return 0;
}

}
}
