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

#include <opm/models/utils/parametersystem.hh>

#include <opm/simulators/flow/DamarisParameters.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_set>

/*
    Below is the std::map with the keywords that are supported by Damaris.

    Most entries in the map below are not critical ('static') and will not
    be changed. We only allow changing FileMode together with output directory
*/

namespace Opm::DamarisOutput {
    
/** 
 * Returns true if the file exists. 
 * Tests to see if filename string is empty 
 * or the "#" character and if so returns false.
 * Tests for file existance on rank 0 and 
 * passes result via MPI to all other ranks.
 */
bool FileExists(const std::string& filename_in,
                const Parallel::Communication& comm);

struct DamarisSettings
{
    bool enableDamarisOutputCollective_ = true;
    bool saveToDamarisHDF5_ = true;
    // if saveMeshToDamarisHDF5 is true, requires enableDamarisOutputCollective to be false 
    // (until offsets are are added to mesh data for collective writing)
    bool saveMeshToHDF5_ = false;  
    std::string pythonFilename_;
    std::string paraviewPythonFilename_;

    std::string damarisSimName_; // empty and set to "opm-flow-<random-number>" if none provided on command line. Used as prefix to HDF5 filenames
    std::string shmemName_;      // empty and needs to be unique if multiple simulations are running on the same server/node. Used to name the Damaris shared memory region.
    std::string damarisLogLevel_ = "info";
    std::string damarisDaskFile_ = "";
    // std::string damarisLimitVars_ = "";
    int nDamarisCores_ = 1;  // this is the number of (Damaris server) cores per node
    int nDamarisNodes_ = 0;
    long shmemSizeBytes_ = 536870912;  // 512 MB
    
    std::string rand_value_str_ ;  // to be added to sheared memory name to make unique 

    std::map<std::string, std::string>
    getKeywords(const Parallel::Communication& comm,
                const std::string& OutputDir);
                
    void SetRandString(void);  // sets the value of rand_value_str_
};

/**
 * Creates the map of search strings and repacement strings that will be used to
 * modify a templated Damaris XML file which will be used to intialize Damaris.
 * This function will access all the OPM flow comand line arguments related to
 * Damaris and perform checks and logic so as to create a valid XML file.
 * N.B. The created XML file can be overridden using an environment variable
 * FLOW_DAMARIS_XML_FILE that points to a Damaris XML file.
 *
 * N.B. This needs to be called before damaris_init()
 */
template<class TypeTag>
std::map<std::string, std::string>
getDamarisKeywords(const Parallel::Communication& comm, const std::string& OutputDir)
{
    DamarisSettings settings;
    // Get all of the Damaris keywords (except for --enable-damaris,
    // which is used in simulators/flow/Main.hpp)
    // These command line arguments are defined in opm/simulators/flow/DamarisWriter.hpp and
    // defaults are set in opm/simulators/flow/FlowProblemProperties.hpp
    settings.enableDamarisOutputCollective_ = Parameters::get<TypeTag, Parameters::DamarisOutputHdfCollective>();
    settings.saveMeshToHDF5_ = Parameters::get<TypeTag, Parameters::DamarisSaveMeshToHdf>();
    settings.saveToDamarisHDF5_ = Parameters::get<TypeTag, Parameters::DamarisSaveToHdf>();
    settings.pythonFilename_ = Parameters::get<TypeTag, Parameters::DamarisPythonScript>();
    settings.paraviewPythonFilename_ = Parameters::get<TypeTag, Parameters::DamarisPythonParaviewScript>();
    settings.damarisSimName_ = Parameters::get<TypeTag, Parameters::DamarisSimName>();
    settings.nDamarisCores_ = Parameters::get<TypeTag, Parameters::DamarisDedicatedCores>();
    settings.nDamarisNodes_ = Parameters::get<TypeTag, Parameters::DamarisDedicatedNodes>();
    settings.shmemSizeBytes_ = Parameters::get<TypeTag, Parameters::DamarisSharedMemorySizeBytes>();
    settings.shmemName_ = Parameters::get<TypeTag, Parameters::DamarisSharedMemoryName>();
    settings.damarisLogLevel_ = Parameters::get<TypeTag, Parameters::DamarisLogLevel>();
    settings.damarisDaskFile_ = Parameters::get<TypeTag, Parameters::DamarisDaskFile>();

    return settings.getKeywords(comm, OutputDir);
}

template<class TypeTag>
std::unordered_set<std::string>
getSetOfIncludedVariables(void) 
{
    std::unordered_set<std::string> resuset ;
    std::string tstr;
    // The --damaris-limit-variables command line option (defaults to empty string)
    std::string damarisLimitVars = Parameters::get<TypeTag, Parameters::DamarisLimitVariables>();
    std::stringstream ss(damarisLimitVars); 

    // Use while loop to check the getline() function condition.  
    while (std::getline(ss, tstr, ',')) {
        //remove whitespace
        std::string::iterator end_pos = std::remove(tstr.begin(), tstr.end(), ' ');
        tstr.erase(end_pos, tstr.end());
        // place in set (no duplicates possible in set and no empty string)
        if (tstr != "") {
            resuset.insert(tstr) ;
        }
    }
    return resuset;
}

} // namespace Opm::DamarisOutput


#endif
