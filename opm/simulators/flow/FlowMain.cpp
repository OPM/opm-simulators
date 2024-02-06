/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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
#include <opm/simulators/flow/FlowMain.hpp>

#include <opm/common/utility/String.hpp>

#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>

#include <opm/simulators/utils/ParallelFileMerger.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>

namespace Opm {
namespace detail {

void mergeParallelLogFiles(std::string_view output_dir,
                           std::string_view deckFilename,
                           bool enableLoggingFalloutWarning)
{
    namespace fs = ::std::filesystem;
    fs::path output_path(output_dir);
    fs::path deck_filename(deckFilename);
    std::string basename;
    // Strip extension "." and ".DATA"
    std::string extension = uppercase(deck_filename.extension().string());
    if (extension == ".DATA" || extension == ".") {
        basename = uppercase(deck_filename.stem().string());
    } else {
        basename = uppercase(deck_filename.filename().string());
    }
    std::for_each(fs::directory_iterator(output_path),
                  fs::directory_iterator(),
                  detail::ParallelFileMerger(output_path, basename,
                                             enableLoggingFalloutWarning));
}

void handleExtraConvergenceOutput(SimulatorReport& report,
                                  std::string_view option,
                                  std::string_view optionName,
                                  std::string_view output_dir,
                                  std::string_view base_name)
{
    const auto extraConvOutput = ConvergenceOutputConfiguration {
        option, optionName
    };

    if (extraConvOutput.want(ConvergenceOutputConfiguration::Option::Steps)) {
      namespace fs = ::std::filesystem;

      const auto infostep = fs::path{output_dir} / fs::path{base_name}.concat(".INFOSTEP");

      std::ofstream os(infostep);
      report.fullReports(os);
    }
}
void checkAllMPIProcesses()
{
#if HAVE_MPI
    const auto& comm = EclGenericVanguard::comm();
    if (comm.size() > 1)
    {
        // we try to prevent the abort here.
        // For that we need a signal that each process is here.
        // Each process sends  a message to rank 0.
        const int tag = 357912;
        if (comm.rank() == 0)
        {
            // wait for a message from all processes.
            std::vector<MPI_Request> requests(comm.size() - 1, MPI_REQUEST_NULL);
            std::vector<int> data(comm.size()-1);

            for(decltype(comm.size()) i = 1; i < comm.size(); ++i)
            {
                if (auto error = MPI_Irecv(data.data() + (i - 1), 1, MPI_INT, i, tag, comm, requests.data() + (i - 1));
                    error != MPI_SUCCESS) {
                    OpmLog::error(fmt::format("Error: Could not set up MPI receive (error code : {})", error));
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
            }
            std::size_t msgs = comm.size() - 1;
            for(std::size_t tries = 0; msgs >0 && tries < 3; ++tries)
            {
                sleep(3);
                int flag, idx;
                for(auto left_msgs = msgs; left_msgs > 0; --left_msgs)
                {
                    if( auto error = MPI_Testany(comm.size()-1, requests.data(), &idx, &flag, MPI_STATUS_IGNORE);
                        error != MPI_SUCCESS) {
                        OpmLog::error(fmt::format("Error: Could not test for MPI message (error code : {})", error));
                        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                    }
                    if (flag)
                    {
                        --msgs;
                    }
                }
            }
            if (msgs) {
                // seems like some processes are stuck. Abort just to be save
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
        else
        {
            int data= 3;
            MPI_Request request = MPI_REQUEST_NULL;
            if (auto error = MPI_Isend(&data, 1, MPI_INT, 0, tag, comm, &request);
                error != MPI_SUCCESS) {
                OpmLog::error(fmt::format("Error: Could send MPI message (error code : {})", error));
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            bool completed = false;
            for(std::size_t tries = 0; !completed && tries < 3; tries++)
            {
                sleep(3);
                int flag;
                if( auto error = MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
                    error != MPI_SUCCESS) {
                    OpmLog::error(fmt::format("Error: Could not test for MPI message (error code : {})", error));
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
                if (flag)
                {
                    completed = true;
                }
            }
            if (!completed) {
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
    }
#endif
}

} // namespace detail
} // namespace Opm
