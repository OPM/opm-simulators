/*
  Copyright 2025, SINTEF Digital

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

#include <opm/simulators/flow/NlddReporting.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <fmt/format.h>

#include <sstream>
#include <vector>

namespace Opm {

/**
 * Reports NLDD statistics after simulation.
 *
 * @param domain_reports The accumulated reports per domain
 * @param local_report The accumulated reports per rank
 * @param output_cout Whether to output to cout
 * @param comm The communication object for parallel runs
 */
void reportNlddStatistics(const std::vector<SimulatorReport>& domain_reports,
                          const SimulatorReport& local_report,
                          const bool output_cout,
                          const Parallel::Communication& comm)
{
    const auto& mpi_rank = comm.rank();
    {
        // Create a deferred logger and add the report with rank as tag
        DeferredLogger local_log;
        std::ostringstream ss;

        // Accumulate reports per domain
        const auto dr_size = domain_reports.size();
        for (auto i = 0*dr_size; i < dr_size; ++i) {
            const auto& dr = domain_reports[i];
            ss << "======  Accumulated local solve data for domain " << i << " on rank " << mpi_rank << " ======\n";
            dr.reportNLDD(ss);
            // Use combined rank and domain index as tag to ensure global uniqueness and correct ordering
            local_log.debug(fmt::format("R{:05d}D{:05d}", mpi_rank, i), ss.str());
            ss.str(""); // Clear the stringstream
            ss.clear(); // Clear any error flags
        }
        // Gather all logs and output them in sorted order
        auto global_log = gatherDeferredLogger(local_log, comm);
        if (output_cout) {
            global_log.logMessages();
        }
    }

    {
        // Create a deferred logger and add the report with rank as tag
        DeferredLogger local_log;
        std::ostringstream ss;
        ss << "======  Accumulated local solve data for rank " << mpi_rank << " ======\n";
        local_report.reportNLDD(ss);
        // Use rank number as tag to ensure correct ordering
        local_log.debug(fmt::format("{:05d}", mpi_rank), ss.str());

        // Gather all logs and output them in sorted order
        auto global_log = gatherDeferredLogger(local_log, comm);
        if (output_cout) {
            global_log.logMessages();
        }
    }
}

} // namespace Opm
