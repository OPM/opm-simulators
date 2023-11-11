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

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <ebos/damariswriter.hh>

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
*     Tests for file existance on ranl 0 and 
*     passes result via MPI to all other ranks.
*/
bool FileExists(const std::string& filename_in,
                const Parallel::Communication& comm);

struct DamarisSettings {
    bool enableDamarisOutputCollective = true ;
    bool saveToDamarisHDF5 = true ;
    std::string pythonFilename;
    std::string paraviewPythonFilename;

    std::string damarisSimName; // empty defaults to opm-sim-<magic_number>
    std::string damarisLogLevel = "info";
    std::string damarisDaskFile = "" ;
    int nDamarisCores   = 1 ;
    int nDamarisNodes   = 0 ;
    long shmemSizeBytes = 536870912 ;

    std::map<std::string, std::string>
    getKeywords(const Parallel::Communication& comm,
                const std::string& OutputDir);
};

/**
*   Creates the map of search strings and repacement strings that will be used to
*   modify a templated Damaris XML file which will be used to intialize Damaris.
*   This function will access all the OPM flow comand line arguments related to
*   Damaris and perform checks and logic so as to create a valid XML file.
*   N.B. The created XML file can be overridden using an environment variable
*   FLOW_DAMARIS_XML_FILE that points to a Damaris XML file.
*/
template<class TypeTag>
std::map<std::string, std::string>
DamarisKeywords(const Parallel::Communication& comm, const std::string& OutputDir)
{
    DamarisSettings settings;
    // Get all of the Damaris keywords (except for --enable-damaris, which is used in simulators/flow/Main.hpp)
    // These command line arguments are defined in ebos/damariswriter.hh and defaults are set in ebos/eclproblem_properties.hh
    settings.enableDamarisOutputCollective = EWOMS_GET_PARAM(TypeTag, bool, EnableDamarisOutputCollective) ;
    settings.saveToDamarisHDF5 = EWOMS_GET_PARAM(TypeTag, bool, DamarisSaveToHdf);
    settings.pythonFilename = EWOMS_GET_PARAM(TypeTag, std::string, DamarisPythonScript);
    settings.paraviewPythonFilename = EWOMS_GET_PARAM(TypeTag, std::string, DamarisPythonParaviewScript);
    settings.damarisSimName = EWOMS_GET_PARAM(TypeTag, std::string, DamarisSimName);
    settings.nDamarisCores = EWOMS_GET_PARAM(TypeTag, int, DamarisDedicatedCores);
    settings.nDamarisNodes = EWOMS_GET_PARAM(TypeTag, int, DamarisDedicatedNodes);
    settings.shmemSizeBytes = EWOMS_GET_PARAM(TypeTag, long, DamarisSharedMemeorySizeBytes);
    settings.damarisLogLevel = EWOMS_GET_PARAM(TypeTag, std::string, DamarisLogLevel);
    settings.damarisDaskFile = EWOMS_GET_PARAM(TypeTag, std::string, DamarisDaskFile);
    return settings.getKeywords(comm, OutputDir);
}

} // namespace Opm::DamarisOutput


#endif