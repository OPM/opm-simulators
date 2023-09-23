/*
  Copyright 2021 Equinor.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_DAMARISKEYWORDS_HEADER_INCLUDED
#define OPM_DAMARISKEYWORDS_HEADER_INCLUDED

#include <string>
#include <map>
#include <mpi.h>
/*
    Below is the std::map with the keywords that are supported by Damaris.

    Most entries in the map below are not critical ('static') and will not
    be changed. We only allow changing FileMode together with output directory
*/


namespace Opm::DamarisOutput
{
    
    /** 
*   Returns true if the file exists. 
*     Tests to see if filename string is empty 
*     or the "#" character and if so returns false.
*/
bool FileExists(const std::string filename_in, const MPI_Comm comm, const int rank) ;


/** 
*   Returns true if the file exists. Tests to see if filename string is empty or the "#" character and if so returns false.
*/
std::map<std::string, std::string>
DamarisKeywords(MPI_Comm comm, std::string OutputDir, 
                    bool enableDamarisOutputCollective, 
                    bool saveToHDF5, 
                    int  nDamarisCores,
                    int  nDamarisNodes,
                    long shmemSizeBytes,
                    std::string pythonFilename, 
                    std::string simName, 
                    std::string logLevel,
                    std::string paraviewPythonFilename ) ;

} // namespace Opm::DamarisOutput


#endif
