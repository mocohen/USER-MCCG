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
#include "Action.h"
#include "ActionWithValue.h"
#include "PlumedMain.h"
#include "tools/Log.h"
#include "tools/Exception.h"
#include "Atoms.h"
#include "ActionSet.h"
#include <iostream>

namespace PLMD{

Keywords ActionOptions::emptyKeys;

ActionOptions::ActionOptions(PlumedMain&p,const std::vector<std::string>&l):
plumed(p),
line(l),
keys(emptyKeys)
{
}

ActionOptions::ActionOptions(const ActionOptions&ao,const Keywords&keys):
plumed(ao.plumed),
line(ao.line),
keys(keys)
{
}

void Action::registerKeywords( Keywords& keys ){
  plumed_assert( keys.size()==0 );
  keys.add( "hidden", "LABEL", "a label for the action so that its output can be referenced in the input to other actions.  Actions with scalar output are referenced using their label only.  Actions with vector output must have a separate label for every component.  Individual componets are then refered to using label.component" );
}

Action::Action(const ActionOptions&ao):
  name(ao.line[0]),
  line(ao.line),
  active(false),
  plumed(ao.plumed),
  log(plumed.getLog()),
  comm(plumed.comm),
  multi_sim_comm(plumed.multi_sim_comm),
  keywords(ao.keys)
{
  line.erase(line.begin());
  log.printf("Action %s\n",name.c_str());

  if ( keywords.exists("LABEL") ){ parse("LABEL",label); }

  if(label.length()==0){
    std::string s; Tools::convert(plumed.getActionSet().size(),s);
    label="@"+s;
  }
  if( plumed.getActionSet().selectWithLabel<Action*>(label) ) error("label " + label + " has been already used");
  log.printf("  with label %s\n",label.c_str());
}

Action::~Action(){
  if(files.size()!=0){
    std::cerr<<"WARNING: some files open in action "+getLabel()+" where not properly closed. This could lead to data loss!!\n";
  }
}

FILE* Action::fopen(const char *path, const char *mode){
  bool write(false);
  for(const char*p=mode;*p;p++) if(*p=='w' || *p=='a' || *p=='+') write=true;
  FILE* fp;
  if(write && comm.Get_rank()!=0) fp=plumed.fopen("/dev/null",mode);
  else      fp=plumed.fopen(path,mode);
  files.insert(fp);
  return fp;
}

int Action::fclose(FILE*fp){
  files.erase(fp);
  return plumed.fclose(fp);
}

void Action::fflush(){
  for(files_iterator p=files.begin();p!=files.end();++p){
    std::fflush((*p));
  }
}

void Action::parseFlag(const std::string&key,bool & t){
  // Check keyword has been registered
  plumed_massert(keywords.exists(key), "keyword " + key + " has not been registered");
  // Check keyword is a flag
  if(!keywords.style(key,"nohtml")){
     plumed_massert(keywords.style(key,"flag") || keywords.style(key,"hidden"), "keyword " + key + " is not a flag");
  }

  // Read in the flag otherwise get the default value from the keywords object
  if(!Tools::parseFlag(line,key,t)){
     if( keywords.style(key,"nohtml") ){ 
        t=false; 
     } else if ( !keywords.getLogicalDefault(key,t) ){
        log.printf("ERROR in action %s with label %s : flag %s has no default",name.c_str(),label.c_str(),key.c_str() );
        plumed_error();
     } 
  }
}

void Action::addDependency(Action*action){
  after.push_back(action);
}

void Action::activate(){
// preparation step is called only the first time an Action is activated.
// since it could change its dependences (e.g. in an ActionAtomistic which is
// accessing to a virtual atom), this is done just before dependencies are
// activated
  if(!active){
    this->unlockRequests();
    prepare();
    this->lockRequests();
  }
  for(Dependencies::iterator p=after.begin();p!=after.end();++p) (*p)->activate();
  active=true;
}

void Action::setOption(const std::string &s){
// This overloads the action and activate some options  
  options.insert(s);
  for(Dependencies::iterator p=after.begin();p!=after.end();++p) (*p)->setOption(s);
}

void Action::clearOptions(){
// This overloads the action and activate some options  
  options.clear();
}


void Action::clearDependencies(){
  after.clear();
}

std::string Action::getDocumentation()const{
  return std::string("UNDOCUMENTED ACTION");
}

void Action::checkRead(){
  if(!line.empty()){
    std::string msg="cannot understand the following words from the input line : ";
    for(unsigned i=0;i<line.size();i++) msg = msg + line[i] + ", ";
    error(msg);
  }
}

long int Action::getStep()const{
  return plumed.getStep();
}

double Action::getTime()const{
  return plumed.getAtoms().getTimeStep()*getStep();
}

double Action::getTimeStep()const{
  return plumed.getAtoms().getTimeStep();
}



void Action::exit(int c){
  plumed.exit(c);
}

void Action::calculateNumericalDerivatives( ActionWithValue* a ){
  plumed_merror("if you get here it means that you are trying to use numerical derivatives for a class that does not implement them");
}

void Action::prepare(){
  return;
}

void Action::error( const std::string & msg ) const {
  log.printf("ERROR in input to action %s with label %s : %s \n \n", name.c_str(), label.c_str(), msg.c_str() );
  if( !line.empty() ) keywords.print( log );
  plumed_merror("ERROR in input to action " + name + " with label " + label + " : " + msg );
}

void Action::warning( const std::string & msg ){
  log.printf("WARNING for action %s with label %s : %s \n", name.c_str(), label.c_str(), msg.c_str() ); 
}

void Action::calculateFromPDB( const PDB& pdb ){
  activate();
  for(Dependencies::iterator p=after.begin();p!=after.end();++p){
     ActionWithValue*av=dynamic_cast<ActionWithValue*>(*p);
     if(av){ av->clearInputForces(); av->clearDerivatives(); }
     (*p)->readAtomsFromPDB( pdb ); 
     (*p)->calculate();
  }
  readAtomsFromPDB( pdb );
  calculate();
}

bool Action::getExchangeStep()const{
  return plumed.getExchangeStep();
}

std::string Action::cite(const std::string&s){
  return plumed.cite(s);
}


}
