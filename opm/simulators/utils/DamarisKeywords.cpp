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
    
    The command line arguments are defined in ebos/damariswriter.hh 
    and defaults are set in ebos/eclproblem_properties.hh
*/

namespace Opm::DamarisOutput
{

bool FileExists(const std::string& filename_in,
                const Parallel::Communication& comm)
{
    // From c++17  : std::filesystem::exists(filename_in);
    
    int retint = 0;
    std::ifstream filestr;
    bool file_exists = false;
    
    if ((filename_in.length() == 0) || (filename_in == "#") ) {
        return file_exists;
    }

    if (comm.rank() == 0) {
        filestr.open(filename_in);
        file_exists = true;
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

std::map<std::string, std::string>
DamarisSettings::getKeywords([[maybe_unused]] const Parallel::Communication& comm,
                             const std::string& OutputDir)
{
    SetRandString() ;  // sets rand_value_str_ used for naming things that might need a unique name

    std::string saveToHDF5_str("MyStore");
    if (! saveToDamarisHDF5_ ){
        saveToHDF5_str = "#";
    }

    // These strings are used to comment out an XML element if it is not reqired
    std::string disablePythonXMLstart("!--");
    std::string disablePythonXMLfin("--");
    std::string disableParaviewXMLstart("!--");
    std::string disableParaviewXMLfin("--");

    std::string publishToPython_str("#"); // to be changed to the name of the PyScript XML element
#ifdef HAVE_PYTHON_ENABLED
    // Test if input Python file exists and set the name of the script for <variable ...  script="" > )XML elements
    if (pythonFilename_ != ""){
        if (FileExists(pythonFilename_, comm)) {
             publishToPython_str="PythonScript"; // the name of the PyScript XML element
             disablePythonXMLstart.clear();
             disablePythonXMLfin.clear();
        } else {
            pythonFilename_.clear(); // set to empty if it does not exist
            disablePythonXMLstart = std::string("!--");
            disablePythonXMLfin = std::string("--");
        }
    }
#else
     OpmLog::info(fmt::format("INFO: Opm::DamarisOutput::DamarisKeywords() : Python is not enabled in the Damaris library. "
                              "The commandline --damaris-python-script={} will be set to empty string", pythonFilename_));
     pythonFilename_.clear();
#endif

#ifdef HAVE_PARAVIEW_ENABLED
     // Test if input Paraview Python file exists 
    if (paraviewPythonFilename_ != ""){
        if (FileExists(paraviewPythonFilename_, comm)) {
            disableParaviewXMLstart.clear();
            disableParaviewXMLfin.clear();
        } else  {
            paraviewPythonFilename_.clear(); // set to empty if it does not exist
            disableParaviewXMLstart = std::string("!--");
            disableParaviewXMLfin = std::string("--");
        }
    }
#else
     OpmLog::info(fmt::format("Opm::DamarisOutput::DamarisKeywords() : Paraview is not enabled in the Damaris library. "
                              "The commandline --damaris-python-paraview-script={} will be set to empty string", 
                              paraviewPythonFilename_));
     paraviewPythonFilename_.clear();
#endif

    // Flag error if both scripts are enabled 
    if ((pythonFilename_.size() > 0) && (paraviewPythonFilename_.size() > 0) )
    {
        // A work around of this issue is to remove the Paraview mpi4py library (use print(inspect.getfile(mpi4py)))
        // and then possibly not use mpi4py in the Paraview script code. OR try to install paraview mpi4py with headers.
        OPM_THROW(std::runtime_error, "ERROR: Both the Python (--damaris-python-script command line argument) and Paraview Python "
                                      "(--damaris-python-paraview-script command line argument) scripts are valid, however only one "
                                      "type of analysis is supported in a single simulation (due to Paraview installing mpi4py library "
                                      "locally and without header files). "
                                      "Please choose one or the other method of analysis for now. Exiting." );
    }
    
    std::string saveMeshToHDF5_str("#");
    if (saveMeshToHDF5_ == true) {
        enableDamarisOutputCollective_ = false ;
        saveMeshToHDF5_str = "MyStore" ;
    }
    std::string damarisOutputCollective_str;
    if (enableDamarisOutputCollective_) {
        damarisOutputCollective_str = "Collective";
    } else {
        damarisOutputCollective_str = "FilePerCore";
    }
    OpmLog::info(fmt::format("Opm::DamarisOutput::DamarisKeywords() : <option key=\"FileMode\"> {} </option> ", 
                              damarisOutputCollective_str)); 

    std::string  simName_str;
    // Check if simulation name was given on command line
    // The simulation name is used as a prefix to name HDF5 files
    if (damarisSimName_.empty()) {
        simName_str = "opm-flow-" + rand_value_str_;
    } else {
        simName_str = damarisSimName_;
    }
    OpmLog::info(fmt::format("Opm::DamarisOutput::DamarisKeywords() : <simulation name={} ", 
                              simName_str));
    
    // A different shared memory buffer name is important if multiple simulations 
    // are running on the same node, as one simulation will remove the buffer when it exits,
    // which will remove the buffer for other simulations.
    std::string  shmemName_str;
    if ( shmemName_.empty()) {
        shmemName_str = "opm-damaris-" + rand_value_str_;
    } else {
        shmemName_str = shmemName_ ;
    }
    
    std::string shmemSizeBytes_str;
    if (shmemSizeBytes_ != 0) {
        shmemSizeBytes_str = std::to_string(shmemSizeBytes_);
    } else {
        shmemSizeBytes_str = "536870912";  // 512 MB
    }

    OpmLog::info(fmt::format("Opm::DamarisOutput::DamarisKeywords() : <buffer name={} size={} ", 
                              shmemName_str, shmemSizeBytes_str));
    
    if ((nDamarisCores_ > 0) && (nDamarisNodes_ > 0))
    {
        nDamarisNodes_ = 0; // Default is to use Damaris Cores
    }
    std::string nDamarisCores_str;
    if ( nDamarisCores_ != 0 ) {
        nDamarisCores_str = std::to_string(nDamarisCores_);
    } else {
        nDamarisCores_str = "0";
    }

    std::string nDamarisNodes_str;
    if ( nDamarisNodes_ != 0 ) {
        nDamarisNodes_str = std::to_string(nDamarisNodes_);
    } else {
        nDamarisNodes_str = "0";
    }
    
    OpmLog::info(fmt::format("Opm::DamarisOutput::DamarisKeywords() : <dedicated cores={} nodes={} ", 
                              nDamarisCores_str, nDamarisNodes_str));

   
                              
    std::string logLevel_str(damarisLogLevel_);
    std::string logFlush_str("false");
    if ((logLevel_str == "debug") || (logLevel_str == "trace") ) {
        logFlush_str = "true";
    }
    OpmLog::info(fmt::format("Opm::DamarisOutput::DamarisKeywords() : <log FileName={}/damaris_log/{} Flush={}  LogLevel={} ", 
                              OutputDir, simName_str, logFlush_str, logLevel_str));
    
    std::map<std::string, std::string> damaris_keywords = {
        {"_SHMEM_BUFFER_BYTES_REGEX_", shmemSizeBytes_str},
        {"_DC_REGEX_", nDamarisCores_str},
        {"_DN_REGEX_", nDamarisNodes_str},
        {"_File_Mode", damarisOutputCollective_str},
        {"_MORE_VARIABLES_REGEX_", ""},
        {"_PATH_REGEX_", OutputDir},   /* Do Not change the string "_PATH_REGEX_" as it is used to search for the output path */
        {"_MYSTORE_OR_EMPTY_REGEX_", saveToHDF5_str},
        {"_MYSTORE_MESH_OR_EMPTY_REGEX_", saveMeshToHDF5_str},
        {"_PARAVIEW_PYTHON_SCRIPT_",paraviewPythonFilename_},  /* this has to be before _PYTHON_SCRIPT_ entry */
        {"_PYTHON_SCRIPT_",pythonFilename_}, /* if a Python script is specified then assume that we want to publish the data to Python */
        {"_PRESSURE_UNIT_","Pa"},
        {"_MAKE_AVAILABLE_IN_PYTHON_",publishToPython_str},  /* must match  <pyscript name="PythonScript" */
        {"_SIM_NAME_",simName_str},
        {"_SHMEM_NAME_",shmemName_str},
        {"_LOG_LEVEL_",logLevel_str},
        {"_LOG_FLUSH_",logFlush_str},
        {"_DISABLEPYTHONSTART_",disablePythonXMLstart},
        {"_DISABLEPYTHONFIN_",disablePythonXMLfin},
        {"_DISABLEPARAVIEWSTART_",disableParaviewXMLstart},
        {"_DISABLEPARAVIEWFIN_",disableParaviewXMLfin},
        {"_DASK_SCHEDULER_FILE_",damarisDaskFile_},
    };
                            
    return damaris_keywords;
}

} // namespace Opm::DamarisOutput
