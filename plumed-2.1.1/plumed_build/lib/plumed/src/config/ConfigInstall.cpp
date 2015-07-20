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

#include "Config.h"
#include "version.h"
#include <cstring>

namespace PLMD{
namespace config{

std::string getSoExt(){
  return "so";
}

bool isInstalled(){
  return true;
}

std::string getPlumedRoot(){
  return "/project/gavoth/mocohen/source/plumed/plumed_mccg/plumed-2.1.1/plumed_build/lib/plumed/";
}

std::string getVersion(){
  return PLUMED_VERSION_SHORT;
}

std::string getVersionLong(){
  return PLUMED_VERSION_LONG;
}

std::string getVersionGit(){
  return PLUMED_VERSION_GIT;
}

std::string getMakefile(){
  static const char conf [] ={
#include "Makefile.conf.xxd"
    , 0x00 };
  return std::string(conf,conf+std::strlen(conf));
}

bool hasMatheval(){
#if __PLUMED_HAS_MATHEVAL
      return true;
#else
      return false;
#endif
}

bool hasDlopen(){
#if __PLUMED_HAS_DLOPEN
      return true;
#else
      return false;
#endif
}

bool hasAlmost(){
#if __PLUMED_HAS_ALMOST
      return true;
#else
      return false;
#endif
}


bool hasCregex(){
#if __PLUMED_HAS_CREGEX
      return true;
#else
      return false;
#endif
}

bool hasMolfile(){
#if __PLUMED_HAS_MOLFILE
      return true;
#else
      return false;
#endif
}

bool hasZlib(){
#if __PLUMED_HAS_ZLIB
      return true;
#else
      return false;
#endif
}




}
}

