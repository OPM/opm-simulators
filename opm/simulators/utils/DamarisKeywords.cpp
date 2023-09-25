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

#include <config.h>
#include <opm/simulators/flow/Main.hpp>



#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/simulators/utils/DamarisKeywords.hpp>
#include <damaris/env/Environment.hpp>
#include <string>
#include <map>
#include <random>
#include <fstream>


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

bool FileExists(const std::string filename_in, const MPI_Comm comm, const int rank) {
    // From c++17  : std::filesystem::exists(filename_in);
    
    int retint = 0 ;
    std::ifstream filestr ;
    bool file_exists = false ;
    
    if ((filename_in.length() == 0) || (filename_in == "#") ) {
        return file_exists ;
    }

    if (rank == 0) {
        
        filestr.open(filename_in);
        file_exists = true ;
        if(filestr.fail()) {
            retint = 0 ;
        } else {
            retint = 1 ;
            filestr.close() ;
        }
        MPI_Bcast(&retint,1,MPI_INT,0,comm);
    } else {
        MPI_Bcast(&retint,1,MPI_INT,0,comm);  // recieve the value from rank 0
    }

    if (retint == 1) {
        file_exists = true ;
    }

    return (file_exists) ;
}


std::map<std::string, std::string>
DamarisKeywords(MPI_Comm comm, std::string OutputDir)
{
    typedef Properties::TTag::FlowEarlyBird PreTypeTag;
    bool enableDamarisOutputCollective = true ;
    bool saveToDamarisHDF5 = true ;
    std::string pythonFilename = "" ;
    std::string paraviewPythonFilename = "" ;

    std::string damarisSimName  = ""       ; // empty defaults to opm-sim-<magic_number>
    std::string damarisLogLevel = "info"   ;
    int nDamarisCores   = 1 ;
    int nDamarisNodes   = 0 ;
    long shmemSizeBytes = 536870912 ;

    // Get all of the Damaris keywords (except for --enable-damaris, which is used in simulators/flow/Main.hpp)
    // These command line arguments are defined in ebos/damariswriter.hh and defaults are set in ebos/eclproblem_properties.hh
    enableDamarisOutputCollective = EWOMS_GET_PARAM(PreTypeTag, bool, EnableDamarisOutputCollective) ;
    saveToDamarisHDF5 = EWOMS_GET_PARAM(PreTypeTag, bool, DamarisSaveToHdf);
    pythonFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, DamarisPythonScript);
    paraviewPythonFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, DamarisPythonParaviewScript);
    damarisSimName = EWOMS_GET_PARAM(PreTypeTag, std::string, DamarisSimName);
    nDamarisCores = EWOMS_GET_PARAM(PreTypeTag, int, DamarisDedicatedCores);
    nDamarisNodes = EWOMS_GET_PARAM(PreTypeTag, int, DamarisDedicatedNodes);
    shmemSizeBytes = EWOMS_GET_PARAM(PreTypeTag, long, DamarisSharedMemeorySizeBytes);
    damarisLogLevel = EWOMS_GET_PARAM(PreTypeTag, std::string, DamarisLogLevel);
    
    int rank ;
    MPI_Comm_rank(comm, &rank) ;
    std::string saveToHDF5_str("MyStore") ;
    if (! saveToDamarisHDF5 )  saveToHDF5_str = "#" ;

    // These strings are used to comment out an XML element if it is not reqired
    std::string disablePythonXMLstart("!--") ;
    std::string disablePythonXMLfin("--") ;
    std::string disableParaviewXMLstart("!--") ;
    std::string disableParaviewXMLfin("--") ;

    // Test if input Python file exists and set the name of the script for <variable ...  script="" > )XML elements
    std::string publishToPython_str("#") ;
    if (pythonFilename != ""){
        if (FileExists(pythonFilename, comm, rank)) {
             publishToPython_str="PythonScript" ; // the name of the PyScript XML element
             disablePythonXMLstart = std::string("")  ;
             disablePythonXMLfin = std::string("")  ;
        } else {
            pythonFilename = "" ; // set to empty if it does not exist
            std::string disablePythonXMLstart("!--") ;
            std::string disablePythonXMLfin("--") ;
        }
    }

     // Test if input Paraview Python file exists 
    if (paraviewPythonFilename != ""){
        if (FileExists(paraviewPythonFilename, comm, rank)) {
            disableParaviewXMLstart = std::string("")  ;
            disableParaviewXMLfin = std::string("")  ;
        } else  {
            paraviewPythonFilename = "" ; // set to empty if it does not exist
            disableParaviewXMLstart = std::string("!--")  ;
            disableParaviewXMLfin = std::string("--")  ;
        }
    }

    // Flag error if both scripts are enabled 
    // It would be good to know if Damaris has either one compiled into the library - this is currently not possible
    if ((pythonFilename.size() > 0) && (paraviewPythonFilename.size() > 0) )
    {
        std::cerr << "ERROR: Both the Python (--damaris-python-script command line argument) and Paraview Python " <<
            "(--damaris-python-paraview-script command line argument) scripts are valid, however only one type "
            "of analysis is supported in a single simulation. Please choose one. Exiting." << std::endl ;
        std::exit(-1) ;
    }

    std::string damarisOutputCollective_str("") ;
    if (enableDamarisOutputCollective) {
        damarisOutputCollective_str="Collective" ;
    } else {
        damarisOutputCollective_str="FilePerCore" ;
    }

    std::string  simName_str("")  ;
    if (damarisSimName == "") {
        // Having a different simulation name is important if multiple simulations 
        // are running on the same node, as it is used to name the simulations shmem area
        // and when one sim finishes it removes its shmem file.
        // simName_str =  damaris::Environment::GetMagicNumber(comm) ;
        if (simName_str == "") {
            // We will add a random value as GetMagicNumber(comm) requires Damaris v1.9.3
            // Seed with a real random value, if available
            std::random_device r;
            // Choose a random number between 0 and MAX_INT
            std::default_random_engine e1(r());
            std::uniform_int_distribution<int> uniform_dist(0, std::numeric_limits<int>::max());
            int rand_int = uniform_dist(e1);
            simName_str = "opm-flow-" + std::to_string(rand_int) ;
        } else {
            simName_str = "opm-flow-" + simName_str ;
        }
    } else {
        simName_str = damarisSimName ;
    }

    if ((nDamarisCores > 0) && (nDamarisNodes > 0))
    {
        nDamarisNodes = 0 ; // Default is to use Damaris Cores
    }
    std::string nDamarisCores_str  ;
    if ( nDamarisCores != 0 ) {
        nDamarisCores_str = std::to_string(nDamarisCores);
    } else {
        nDamarisCores_str = "0" ;
    }

    std::string nDamarisNodes_str  ;
    if ( nDamarisNodes != 0 ) {
        nDamarisNodes_str = std::to_string(nDamarisNodes);
    } else {
        nDamarisNodes_str = "0" ;
    }

    std::string shmemSizeBytes_str  ;
    if ( shmemSizeBytes != 0 ) {
        shmemSizeBytes_str = std::to_string(shmemSizeBytes);
    } else {
        shmemSizeBytes_str = "536870912" ;
    }

    std::string logLevel_str(damarisLogLevel) ;
    std::string logFlush_str("false") ;
    if ((logLevel_str == "debug") || (logLevel_str == "trace") ) {
        logFlush_str = "true" ;
    }

    std::map<std::string, std::string> damaris_keywords = {
        {"_SHMEM_BUFFER_BYTES_REGEX_", shmemSizeBytes_str},
        {"_DC_REGEX_", nDamarisCores_str},
        {"_DN_REGEX_", nDamarisNodes_str},
        {"_File_Mode", damarisOutputCollective_str},
        {"_MORE_VARIABLES_REGEX_", ""},
        {"_PATH_REGEX_", OutputDir},   /* Do Not change the string "_PATH_REGEX_" as it is used to search for the output path */
        {"_MYSTORE_OR_EMPTY_REGEX_", saveToHDF5_str},
        {"_PARAVIEW_PYTHON_SCRIPT_",paraviewPythonFilename},  /* this has to be before _PYTHON_SCRIPT_ entry */
        {"_PYTHON_SCRIPT_",pythonFilename}, /* if a Python script is specified then assume that we want to publish the data to Python */
        {"_PRESSURE_UNIT_","Pa"},
        {"_MAKE_AVAILABLE_IN_PYTHON_",publishToPython_str},  /* must match  <pyscript name="PythonScript" */
        {"_SIM_NAME_",simName_str},
        {"_LOG_LEVEL_",logLevel_str},
        {"_LOG_FLUSH_",logFlush_str},
        {"_DISABLEPYTHONSTART_",disablePythonXMLstart},
        {"_DISABLEPYTHONFIN_",disablePythonXMLfin},
        {"_DISABLEPARAVIEWSTART_",disableParaviewXMLstart},
        {"_DISABLEPARAVIEWFIN_",disableParaviewXMLfin},
    };
    return damaris_keywords;
    
}

} // namespace Opm::DamarisOutput
