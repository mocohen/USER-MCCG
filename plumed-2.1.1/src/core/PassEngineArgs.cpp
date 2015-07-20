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
#include "PassEngineArgs.h"
#include "ActionRegister.h"

using namespace std;

namespace PLMD{
namespace EArgs{

PassEngineArgs::PassEngineArgs(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionAtomistic(ao),
ActionWithArguments(ao)
{
  if(getStride()>1) error("Using Pass Engine Args with stride!=1 is not currently supported");

}

PLUMED_REGISTER_ACTION(PassEngineArgs,"PASSENGINEARGS")

void PassEngineArgs::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
}


void PassEngineArgs::apply(){
  log.printf("MCCG: Engine Args Apply\n");
  log.flush();
  log.printf("MCCG: Check if on step\n");
  log.flush();
  if(onStep()){
    log.printf("MCCG: Is on step\n");
    log.flush();
  	vector<double> &localEngineArgs(modifyEngineArgs());
    log.printf("MCCG: Before Loop\n");
    log.flush();
  	for(unsigned i=0; i<getNumberOfArguments(); i++){
      log.printf("MCCG: In Loop %d %d %f %d\n", i, getNumberOfArguments(), getArgument(i), localEngineArgs.size());
      log.flush();
  		localEngineArgs[i] = getArgument(i);
  	}
  }

}

void PassEngineArgs::calculate(){
	
}

}
}