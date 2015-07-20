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
//***** MCCG
//***** END MCCG 
#ifndef __PLUMED_core_PassEngineArgs_h
#define __PLUMED_core_PassEngineArgs_h

#include "ActionPilot.h"
#include "ActionAtomistic.h"
#include "ActionWithArguments.h"

#define PLUMED_EARGS_INIT(ao) Action(ao),PassEngineArgs(ao)

namespace PLMD{
namespace EArgs{



class PassEngineArgs :
  public ActionPilot,
  public ActionAtomistic,
  public ActionWithArguments
{
  //std::vector<double> outputForces;
//protected:
  //void resetOutputForces();
  //void setOutputForce(int i,double g);
public:
  static void registerKeywords(Keywords&);
  PassEngineArgs(const ActionOptions&ao);
  void apply();
  void calculate();
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ){ plumed_error(); }
  //unsigned getNumberOfDerivatives();
  //void turnOnDerivatives();
};


inline 
void PassEngineArgs::lockRequests(){
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
} 

inline
void PassEngineArgs::unlockRequests(){ 
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

}
}

#endif