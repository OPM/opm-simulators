/*
  Copyright 2021 Equinor.
  Copyright 2023 Inria.

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

#include <config.h>
#include <opm/simulators/utils/DamarisKeywords.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <Damaris.h>
#include <damaris/env/Environment.hpp>

#include <fmt/format.h>

#include <fstream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>

/**
    Below are the Damaris Keywords supported by Damaris to be filled
    in the built-in XML file.

    The entries in the map below will be filled by the corresponding
    Damaris Keywords.

    The command line arguments are defined in opm/simulators/flow/DamarisWriter.hpp
    and defaults are set in opm/simulators/flow/FlowProblemProperties.hpp
*/

namespace Opm::DamarisOutput {

bool FileExists(const std::string& filename_in,
                const Parallel::Communication& comm)
{
    // From c++17  : std::filesystem::exists(filename_in);

    int retint = 1;
    std::ifstream filestr;
    bool file_exists = false;

    if ((filename_in.length() == 0) || (filename_in == "#") ) {
        return file_exists;
    }

    if (comm.rank() == 0) {
        filestr.open(filename_in);
        if(filestr.fail()) {
            retint = 0;
        } else {
            retint = 1;
            filestr.close();
        }
    }

    comm.broadcast(&retint, 1, 0);
    if (retint == 1) {
        file_exists = true;
    } else {
        file_exists = false;
    }

    return (file_exists);
}

void DamarisSettings::SetRandString(void)
{
    // rand_value_str_ =  damaris::Environment::GetMagicNumber(comm);  // requires Damaris >= v1.9.2
    // We will create a random value.
    // Seed with a real random value, if available
    std::random_device r;
    // Choose a random number between 0 and MAX_INT
    std::default_random_engine e1(r());
    std::uniform_int_distribution<int> uniform_dist(0, std::numeric_limits<int>::max());
    int rand_int = uniform_dist(e1);

    rand_value_str_ = std::to_string(rand_int) ;
}

std::tuple<DamarisSettings::StartEnd,
           DamarisSettings::StartEnd,
           std::string>
DamarisSettings::checkFiles(const Parallel::Communication& comm)
{
    StartEnd python{"!--", "--"};
    StartEnd paraview{"!--", "--"};
    std::string publishToPython_str("#");

    // Test if input Python file exists and set the name of the script for <variable ...  script="" > )XML elements
    if (!pythonFilename_.empty()) {
        if (FileExists(pythonFilename_, comm)) {
             publishToPython_str = "PythonScript"; // the name of the PyScript XML element
             python.start.clear();
             python.end.clear();
        }
        else {
            pythonFilename_.clear(); // set to empty if it does not exist
        }
    }

     // Test if input Paraview Python file exists
    if (!paraviewPythonFilename_.empty()) {
        if (FileExists(paraviewPythonFilename_, comm)) {
            paraview.start.clear();
            paraview.end.clear();
        } else  {
            paraviewPythonFilename_.clear(); // set to empty if it does not exist
        }
    }

    // Flag error if both scripts are enabled
    if (!pythonFilename_.empty() && !paraviewPythonFilename_.empty()) {
        // A work around of this issue is to remove the Paraview mpi4py library (use print(inspect.getfile(mpi4py)))
        // and then possibly not use mpi4py in the Paraview script code. OR try to install paraview mpi4py with headers.
        // The Python script will issue: ERROR: Initialiazation via BOOST_PYTHON_MODULE(server) failed at import_mpi4py()
        OPM_THROW(std::runtime_error, "ERROR: Both the Python (--damaris-python-script command line argument) and Paraview Python "
                                      "(--damaris-python-paraview-script command line argument) scripts are valid, however only one "
                                      "type of analysis is supported in a single simulation (due to Paraview installing mpi4py library "
                                      "locally and without header files). "
                                      "Please choose one or the other method of analysis for now. Exiting." );
    }

    return {python, paraview, publishToPython_str};
}

std::map<std::string, std::string>
DamarisSettings::getKeywords(const Parallel::Communication& comm,
                             const std::string& OutputDir)
{
    SetRandString() ;  // sets rand_value_str_ used for naming things that might need a unique name

    const std::string saveToHDF5_str = saveToDamarisHDF5_ ? "MyStore" : "#";

    // These strings are used to comment out an XML element if it is not reqired
    const auto& [python, paraview, publishToPython_str] = checkFiles(comm);

    const std::string saveMeshToHDF5_str = saveMeshToHDF5_ ? "MyStore" : "#";
    if (saveMeshToHDF5_) {
        enableDamarisOutputCollective_ = false;
    }

    const std::string damarisOutputCollective_str =
        enableDamarisOutputCollective_ ? "Collective"
                                       : "FilePerCore";

    OpmLog::info(fmt::format(
        "Opm::DamarisOutput::DamarisKeywords() : <option key=\"FileMode\"> {} </option> ",
        damarisOutputCollective_str
    ));

    const std::string simName_str =
        damarisSimName_.empty() ? "opm-flow-" + rand_value_str_
                                : damarisSimName_;

    // Check if simulation name was given on command line
    // The simulation name is used as a prefix to name HDF5 files
    OpmLog::info(fmt::format(
        "Opm::DamarisOutput::DamarisKeywords() : <simulation name={} ",
        simName_str
    ));

    // A different shared memory buffer name is important if multiple simulations
    // are running on the same node, as one simulation will remove the buffer when it exits,
    // which will remove the buffer for other simulations.
    const std::string shmemName_str =
        shmemName_.empty() ? "opm-damaris-" + rand_value_str_
                           : shmemName_;

    const std::string shmemSizeBytes_str =
        shmemSizeBytes_ != 0 ? std::to_string(shmemSizeBytes_)
                             : "536870912"; // 512 MB

    OpmLog::info(fmt::format(
        "Opm::DamarisOutput::DamarisKeywords() : <buffer name={} size={} ",
        shmemName_str,
        shmemSizeBytes_str
    ));

    if (nDamarisCores_ > 0 && nDamarisNodes_ > 0) {
        nDamarisNodes_ = 0; // Default is to use Damaris Cores
    }

    const std::string nDamarisCores_str = std::to_string(nDamarisCores_);
    const std::string nDamarisNodes_str = std::to_string(nDamarisNodes_);

    OpmLog::info(fmt::format(
        "Opm::DamarisOutput::DamarisKeywords() : <dedicated cores={} nodes={} ",
        nDamarisCores_str,
        nDamarisNodes_str
    ));

    const std::string logFlush_str =
        (damarisLogLevel_ == "debug" ||
         damarisLogLevel_ == "trace") ? "true" : "false";

    OpmLog::info(fmt::format(
        "Opm::DamarisOutput::DamarisKeywords() : <log FileName={}/damaris_log/{} Flush={}  LogLevel={} ",
        OutputDir,
        simName_str,
        logFlush_str,
        damarisLogLevel_
    ));

    std::map<std::string, std::string> damaris_keywords = {
        {"_SHMEM_BUFFER_BYTES_REGEX_",    shmemSizeBytes_str},
        {"_DC_REGEX_",                    nDamarisCores_str},
        {"_DN_REGEX_",                    nDamarisNodes_str},
        {"_File_Mode",                    damarisOutputCollective_str},
        {"_MORE_VARIABLES_REGEX_",        ""},
        {"_PATH_REGEX_",                  OutputDir},   /* Do Not change the string "_PATH_REGEX_" as it is used to search for the output path */
        {"_MYSTORE_OR_EMPTY_REGEX_",      saveToHDF5_str},
        {"_MYSTORE_MESH_OR_EMPTY_REGEX_", saveMeshToHDF5_str},
        {"_PARAVIEW_PYTHON_SCRIPT_",      paraviewPythonFilename_},  /* this has to be before _PYTHON_SCRIPT_ entry */
        {"_PYTHON_SCRIPT_",               pythonFilename_}, /* if a Python script is specified then assume that we want to publish the data to Python */
        {"_PRESSURE_UNIT_",               "Pa"},
        {"_MAKE_AVAILABLE_IN_PYTHON_",    publishToPython_str},  /* must match  <pyscript name="PythonScript" */
        {"_SIM_NAME_",                    simName_str},
        {"_SHMEM_NAME_",                  shmemName_str},
        {"_LOG_LEVEL_",                   damarisLogLevel_},
        {"_LOG_FLUSH_",                   logFlush_str},
        {"_DISABLEPYTHONSTART_",          python.start},
        {"_DISABLEPYTHONFIN_",            python.end},
        {"_DISABLEPARAVIEWSTART_",        paraview.start},
        {"_DISABLEPARAVIEWFIN_",          paraview.end},
        {"_DASK_SCHEDULER_FILE_",         damarisDaskFile_},
    };

    return damaris_keywords;
}

std::map<std::string, std::string>
getDamarisKeywords(const Parallel::Communication& comm, const std::string& OutputDir)
{
    DamarisSettings settings;
    // Get all of the Damaris keywords (except for --enable-damaris,
    // which is used in simulators/flow/Main.hpp)
    // These command line arguments are defined in opm/simulators/flow/DamarisWriter.hpp and
    // defaults are set in opm/simulators/flow/FlowProblemProperties.hpp
    settings.enableDamarisOutputCollective_ = Parameters::Get<Parameters::DamarisOutputHdfCollective>();
    settings.saveMeshToHDF5_ = Parameters::Get<Parameters::DamarisSaveMeshToHdf>();
    settings.saveToDamarisHDF5_ = Parameters::Get<Parameters::DamarisSaveToHdf>();
    settings.pythonFilename_ = Parameters::Get<Parameters::DamarisPythonScript>();
    settings.paraviewPythonFilename_ = Parameters::Get<Parameters::DamarisPythonParaviewScript>();
    settings.damarisSimName_ = Parameters::Get<Parameters::DamarisSimName>();
    settings.nDamarisCores_ = Parameters::Get<Parameters::DamarisDedicatedCores>();
    settings.nDamarisNodes_ = Parameters::Get<Parameters::DamarisDedicatedNodes>();
    settings.shmemSizeBytes_ = Parameters::Get<Parameters::DamarisSharedMemorySizeBytes>();
    settings.shmemName_ = Parameters::Get<Parameters::DamarisSharedMemoryName>();
    settings.damarisLogLevel_ = Parameters::Get<Parameters::DamarisLogLevel>();
    settings.damarisDaskFile_ = Parameters::Get<Parameters::DamarisDaskFile>();

    return settings.getKeywords(comm, OutputDir);
}

std::unordered_set<std::string>
getSetOfIncludedVariables()
{
    std::unordered_set<std::string> resuset ;
    std::string tstr;
    // The --damaris-limit-variables command line option (defaults to empty string)
    std::string damarisLimitVars = Parameters::Get<Parameters::DamarisLimitVariables>();
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
