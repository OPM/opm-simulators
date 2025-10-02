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

#include <opm/models/utils/parametersystem.hpp>

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
std::map<std::string, std::string>
getDamarisKeywords(const Parallel::Communication& comm, const std::string& OutputDir);

std::unordered_set<std::string> getSetOfIncludedVariables();

} // namespace Opm::DamarisOutput


#endif
