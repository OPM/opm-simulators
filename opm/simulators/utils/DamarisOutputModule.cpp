/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2023 INRIA

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

#define XSD_CXX11_TEMPLATE_ALIAS 1

#include <damaris/model/ModifyModel.hpp>
#include <opm/simulators/utils/DamarisKeywords.hpp>
#include <opm/simulators/utils/DamarisOutputModule.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

namespace Opm::DamarisOutput
{

std::string initDamarisXmlFile(); // Defined in initDamarisXMLFile.cpp, to avoid messing up this file.


void
initializeDamaris(MPI_Comm comm, int mpiRank, std::string outputDir, bool enableDamarisOutputCollective)
{
    if (outputDir.empty()) {
        outputDir = ".";
    }
    // Prepare the XML file
    std::string damaris_config_xml = initDamarisXmlFile();
    damaris::model::ModifyModel myMod = damaris::model::ModifyModel(damaris_config_xml);
    // The map will make it precise the output directory and FileMode (either FilePerCore or Collective storage)
    // The map file find all occurences of the string in position 1 and repalce it/them with string in position 2
    std::map<std::string, std::string> find_replace_map = DamarisKeywords(outputDir, enableDamarisOutputCollective);
    myMod.RepalceWithRegEx(find_replace_map);
    std::string damaris_xml_filename_str = outputDir + "/damaris_config.xml";

    if (mpiRank == 0) {
        myMod.SaveXMLStringToFile(damaris_xml_filename_str);
    }

    int damaris_err;

    /* Get the name of the Damaris input file from an environment variable if available */
    const char* cs_damaris_xml_file = getenv("FLOW_DAMARIS_XML_FILE");
    if (cs_damaris_xml_file != NULL) {
        std::cout << "INFO: initializing Damaris from environment variable FLOW_DAMARIS_XML_FILE: "
                  << cs_damaris_xml_file << std::endl;
        damaris_err = damaris_initialize(cs_damaris_xml_file, MPI_COMM_WORLD);
        if (damaris_err != DAMARIS_OK) {
            std::cerr << "ERROR: damaris_initialize() error via FLOW_DAMARIS_XML_FILE=" << cs_damaris_xml_file
                      << std::endl;
        }
    } else {
        std::cout << "INFO: initializing Damaris using internally built file:" << damaris_xml_filename_str << std::endl;
        damaris_err = damaris_initialize(damaris_xml_filename_str.c_str(), comm);
        if (damaris_err != DAMARIS_OK) {
            std::cerr << "ERROR: damaris_initialize() error via built file:" << std::endl << myMod.GetConfigString();
        }
    }
}

} // namespace Opm::DamarisOutput
