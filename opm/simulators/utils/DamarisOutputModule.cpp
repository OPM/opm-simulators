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

 // Initialize Damaris by filling in th XML file and storing it in the chosen directory
void
initializeDamaris(MPI_Comm comm, int mpiRank, std::map<std::string, std::string>& find_replace_map )
{
    int dam_err_;

    /* Get the name of the Damaris input file from an environment variable if available */
    const char* cs_damaris_xml_file = getenv("FLOW_DAMARIS_XML_FILE");
    if (cs_damaris_xml_file != NULL) 
    {
        std::cout << "INFO: Initializing Damaris from environment variable FLOW_DAMARIS_XML_FILE: "
                  << cs_damaris_xml_file << std::endl;
        dam_err_ = damaris_initialize(cs_damaris_xml_file, MPI_COMM_WORLD);
        if (dam_err_ != DAMARIS_OK) {
           OpmLog::error(fmt::format("ERORR: damariswriter::initializeDamaris()       : ( rank:{}) "
                                               "damaris_initialize({}, MPI_COMM_WORLD), Damaris Error: {}  ",
                                               mpiRank, cs_damaris_xml_file, damaris_error_string(dam_err_) ));
        }
    } else {
        // Prepare the XML file
        std::string damaris_config_xml = initDamarisXmlFile();  // This is the template for a Damaris XML file
        damaris::model::ModifyModel myMod = damaris::model::ModifyModel(damaris_config_xml);
        // The map will make it precise the output directory and FileMode (either FilePerCore or Collective storage)
        // The map file find all occurences of the string in position 1 and replace it/them with string in position 2
        // std::map<std::string, std::string> find_replace_map = DamarisKeywords(outputDir, enableDamarisOutputCollective);
        myMod.RepalceWithRegEx(find_replace_map);

        std::string outputDir = find_replace_map["_PATH_REGEX_"] ;
        std::string damaris_xml_filename_str = outputDir + "/damaris_config.xml";

        if (mpiRank == 0) {
            myMod.SaveXMLStringToFile(damaris_xml_filename_str);
        }
        std::cout << "INFO: Initializing Damaris using internally built file:" << damaris_xml_filename_str << " (N.B. use FLOW_DAMARIS_XML_FILE to override)" << std::endl;
        dam_err_ = damaris_initialize(damaris_xml_filename_str.c_str(), comm);
        if (dam_err_ != DAMARIS_OK) {
            OpmLog::error(fmt::format("ERORR: damariswriter::initializeDamaris()       : ( rank:{}) "
                                               "damaris_initialize({}, MPI_COMM_WORLD), Damaris Error: {}.  Error via OPM internally built file:",
                                               mpiRank, cs_damaris_xml_file, damaris_error_string(dam_err_) ));
        }
    }
}

} // namespace Opm::DamarisOutput
