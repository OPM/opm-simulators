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
#include <opm/simulators/utils/DamarisOutputModule.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/utils/DamarisKeywords.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#define XSD_CXX11_TEMPLATE_ALIAS 1

#include <Damaris.h>
#include <damaris/model/ModifyModel.hpp>
#include <fmt/format.h>

namespace Opm::DamarisOutput
{

std::string initDamarisXmlFile(); // Defined in initDamarisXMLFile.cpp, to avoid messing up this file.

/**
 * Initialize Damaris by either reading a file specified by the environment variable FLOW_DAMARIS_XML_FILE or
 * by  filling in the XML file and storing it in the chosen directory
 */
void
initializeDamaris(const Parallel::Communication comm, const int mpiRank, const std::map<std::string, std::string>& find_replace_map )
{
    int dam_err;

    /* Get the name of the Damaris input file from an environment variable if available */
    const char* cs_damaris_xml_file = getenv("FLOW_DAMARIS_XML_FILE");
    if (cs_damaris_xml_file != NULL)
    {
        OpmLog::info(std::string("Initializing Damaris from environment variable FLOW_DAMARIS_XML_FILE: ") + cs_damaris_xml_file);
        dam_err = damaris_initialize(cs_damaris_xml_file, comm);
        if (dam_err != DAMARIS_OK) {
           OpmLog::error(fmt::format("ERORR: damariswriter::initializeDamaris()       : ( rank:{}) "
                                               "damaris_initialize({}, comm), Damaris Error: {}  ",
                                               mpiRank, cs_damaris_xml_file, damaris_error_string(dam_err) ));
        }
    } else {
        // Prepare the inbuilt XML file
        std::string damaris_config_xml = initDamarisXmlFile();  // This is the template for a Damaris XML file
        damaris::model::ModifyModel myMod = damaris::model::ModifyModel(damaris_config_xml);
        // The map file find all occurences of the string in position 1 and replace it/them with string in position 2
        // std::map<std::string, std::string> find_replace_map = DamarisKeywords(outputDir, enableDamarisOutputCollective);
        myMod.RepalceWithRegEx(find_replace_map);

        std::string outputDir = find_replace_map.at("_PATH_REGEX_");
        std::string damaris_xml_filename_str = outputDir + "/damaris_config.xml";

        if (mpiRank == 0) {
            myMod.SaveXMLStringToFile(damaris_xml_filename_str);
        }

        OpmLog::info("Initializing Damaris using internally built file: " + damaris_xml_filename_str +
                     " (N.B. use environment variable FLOW_DAMARIS_XML_FILE to override)");

        dam_err = damaris_initialize(damaris_xml_filename_str.c_str(), comm);
        if (dam_err != DAMARIS_OK) {
            OpmLog::error(fmt::format("damariswriter::initializeDamaris()       : ( rank:{}) "
                                      "damaris_initialize({}, comm), Damaris Error: {}.  Error via OPM internally built file:",
                                      mpiRank, cs_damaris_xml_file, damaris_error_string(dam_err) ));
        }
    }
}

} // namespace Opm::DamarisOutput
