/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.

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

void
setupDamarisWritingPars(Parallel::Communication comm, const int n_elements_local_grid)
{
    int damaris_err = DAMARIS_OK;

    const int nranks = comm.size();
    const int rank = comm.rank();

    std::vector<unsigned long long> elements_rank_sizes(nranks); // one for each rank -- to be gathered from each client rank
    std::vector<unsigned long long> elements_rank_offsets(nranks); // one for each rank, first one 0 -- to be computed - Probably could use MPI_Scan()!

    // n_elements_local_grid should be the full model size
    const unsigned long long n_elements_local = n_elements_local_grid;

    // Left in for debugging, but commented out to avoid spamming the terminal from non-output ranks.
    // std::cout << "INFO (" << rank << "): n_elements_local_grid   = " << n_elements_local_grid << std::endl;

    // This gets the n_elements_local from all ranks and copies them to a std::vector of all the values on all ranks
    // (elements_rank_sizes[]).
    comm.allgather(&n_elements_local, 1, elements_rank_sizes.data());
    elements_rank_offsets[0] = 0ULL;
    // This scan makes the offsets to the start of each ranks grid section if each local grid data was concatenated (in
    // rank order)
    for (int t1 = 1; t1 < nranks; t1++) {
        elements_rank_offsets[t1] = elements_rank_offsets[t1 - 1] + elements_rank_sizes[t1 - 1];
    }

    // find the global/total size
    unsigned long long n_elements_global_max = elements_rank_offsets[nranks - 1];
    n_elements_global_max += elements_rank_sizes[nranks - 1]; // add the last ranks size to the already accumulated offset values

    if (rank == 0) {
        OpmLog::debug(fmt::format("In setupDamarisWritingPars(): n_elements_global_max = {}", n_elements_global_max));
    }

    // Set the paramater so that the Damaris servers can allocate the correct amount of memory for the variabe
    // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
    // ToDo: Do we need to check that local ranks are 0 based ?
    int temp_int = static_cast<int>(elements_rank_sizes[rank]);
    damaris_err = damaris_parameter_set("n_elements_local", &temp_int, sizeof(int));
    if (damaris_err != DAMARIS_OK && rank == 0) {
        OpmLog::error("Damaris library produced an error result for "
                      "damaris_parameter_set(\"n_elements_local\", &temp_int, sizeof(int));");
    }
    // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
    // ToDo: Do we need to check that n_elements_global_max will fit in a C int type (INT_MAX)
    temp_int = static_cast<int>(n_elements_global_max);
    damaris_err = damaris_parameter_set("n_elements_total", &temp_int, sizeof(int));
    if (damaris_err != DAMARIS_OK && rank == 0) {
        OpmLog::error("Damaris library produced an error result for "
                      "damaris_parameter_set(\"n_elements_total\", &temp_int, sizeof(int));");
    }

    // Use damaris_set_position to set the offset in the global size of the array.
    // This is used so that output functionality (e.g. HDF5Store) knows global offsets of the data of the ranks
    int64_t temp_int64_t[1];
    temp_int64_t[0] = static_cast<int64_t>(elements_rank_offsets[rank]);
    damaris_err = damaris_set_position("PRESSURE", temp_int64_t);
    if (damaris_err != DAMARIS_OK && rank == 0) {
        OpmLog::error("Damaris library produced an error result for "
                      "damaris_set_position(\"PRESSURE\", temp_int64_t);");
    }
    damaris_err = damaris_set_position("GLOBAL_CELL_INDEX", temp_int64_t);
    if (damaris_err != DAMARIS_OK && rank == 0) {
        OpmLog::error("Damaris library produced an error result for "
                      "damaris_set_position(\"GLOBAL_CELL_INDEX\", temp_int64_t);");
    }
}
} // namespace Opm::DamarisOutput
