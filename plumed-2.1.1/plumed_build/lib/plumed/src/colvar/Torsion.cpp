/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Torsion.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR TORSION
/*
Calculate a torsional angle.

This command can be used to compute the torsion between four atoms or alternatively
to calculate the angle between two vectors projected on the plane
orthogonal to an axis. 

\par Examples

This input tells plumed to print the torsional angle between atoms 1, 2, 3 and 4
on file COLVAR.
\verbatim
t: TORSION ATOMS=1,2,3,4
# this is an alternative, equivalent, definition:
# t: TORSION VECTOR1=2,1 AXIS=2,3 VECTOR2=3,4
PRINT ARG=t FILE=COLVAR
\endverbatim

If you are working with a protein you can specify the special named torsion angles \f$\phi\f$, \f$\psi\f$, \f$\omega\f$ and \f$\chi_1\f$
by using TORSION in combination with the \ref MOLINFO command.  This can be done by using the following 
syntax.

\verbatim
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
t1: TORSION ATOMS=@phi-3
t2: TORSION ATOMS=@psi-4
PRINT ARG=t1,t2 FILE=colvar STRIDE=10
\endverbatim

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the 4th residue of the protein.
*/
//+ENDPLUMEDOC
   
class Torsion : public Colvar {
  bool pbc;
  bool do_cosine;

public:
  Torsion(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Torsion,"TORSION")

void Torsion::registerKeywords(Keywords& keys){
   Colvar::registerKeywords( keys );
   keys.add("atoms-1","ATOMS","the four atoms involved in the torsional angle");
   keys.add("atoms-2","AXIS","two atoms that define an axis.  You can use this to find the angle in the plane perpendicular to the axis between the vectors specified using the VECTOR1 and VECTOR2 keywords."); 
   keys.add("atoms-2","VECTOR1","two atoms that define a vector.  You can use this in combination with VECTOR2 and AXIS");
   keys.add("atoms-2","VECTOR2","two atoms that define a vector.  You can use this in combination with VECTOR1 and AXIS");
   keys.addFlag("COSINE",false,"calculate cosine instead of dihedral");
}

Torsion::Torsion(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
do_cosine(false)
{
  vector<AtomNumber> atoms,v1,v2,axis;
  parseAtomList("ATOMS",atoms);
  parseAtomList("VECTOR1",v1);
  parseAtomList("VECTOR2",v2);
  parseAtomList("AXIS",axis);

  parseFlag("COSINE",do_cosine);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  if(atoms.size()==4){
    if(!(v1.empty()) && (v2.empty()) && (axis.empty()))
      error("ATOMS keyword is not compatible with VECTOR1, VECTOR2 and AXIS keywords");
    log.printf("  between atoms %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
    atoms.resize(6);
    atoms[5]=atoms[3];
    atoms[4]=atoms[2];
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  }else if(atoms.empty()){
    if(!(v1.size()==2 && v2.size()==2 && axis.size()==2))
      error("VECTOR1, VECTOR2 and AXIS should specify 2 atoms each");
    log.printf("  between lines %d-%d and %d-%d, projected on the plane orthogonal to line %d-%d\n",
                v1[0].serial(),v1[1].serial(),v2[0].serial(),v2[1].serial(),axis[0].serial(),axis[1].serial());
    atoms.resize(6);
    atoms[0]=v1[1];
    atoms[1]=v1[0];
    atoms[2]=axis[0];
    atoms[3]=axis[1];
    atoms[4]=v2[0];
    atoms[5]=v2[1];
  }else error("ATOMS should specify 4 atoms");

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  if(do_cosine) log.printf("  calculating cosine instead of torsion\n");

  addValueWithDerivatives();
  if(!do_cosine) setPeriodic("-pi","pi");
  else setNotPeriodic();
  requestAtoms(atoms);
}

// calculator
void Torsion::calculate(){

  Vector d0,d1,d2;
  if(pbc){
    d0=pbcDistance(getPosition(1),getPosition(0));
    d1=pbcDistance(getPosition(3),getPosition(2));
    d2=pbcDistance(getPosition(5),getPosition(4));
  } else {
    d0=delta(getPosition(1),getPosition(0));
    d1=delta(getPosition(3),getPosition(2));
    d2=delta(getPosition(5),getPosition(4));
  }
  Vector dd0,dd1,dd2;
  PLMD::Torsion t;
  double torsion=t.compute(d0,d1,d2,dd0,dd1,dd2);
  if(do_cosine){
   dd0 *= -sin(torsion);
   dd1 *= -sin(torsion);
   dd2 *= -sin(torsion);
   torsion = cos(torsion);
  }
  setAtomsDerivatives(0,dd0);
  setAtomsDerivatives(1,-dd0);
  setAtomsDerivatives(2,dd1);
  setAtomsDerivatives(3,-dd1);
  setAtomsDerivatives(4,dd2);
  setAtomsDerivatives(5,-dd2);

  setValue           (torsion);
  setBoxDerivatives  (-(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)));
}

}
}



