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
#ifndef OPM_READDECK_HEADER_INCLUDED
#define OPM_READDECK_HEADER_INCLUDED

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <string>

namespace Opm {
    class EclipseState;
    class ErrorGuard;
    class ParseContext;
    class Python;
    class Schedule;
    class SummaryConfig;
    class UDQState;
    class WellTestState;
} // end namespace Opm

namespace Opm {

namespace Action {
class State;
}

enum class FileOutputMode {
    //! \brief No file output.
    OUTPUT_NONE = 0,

    //! \brief Output only to log files, no ECLIPSE output.
    OUTPUT_LOG_ONLY = 1,

    //! \brief Output to all files.
    OUTPUT_ALL = 3,
};

// Ensure that a directory exists, creating it if it does not.
void
ensureOutputDirExists(const std::string& cmdline_output_dir);

// Prepare the result ouptut directory by removing files from previous simulation runs.
void
prepareResultOutputDirectory(const std::string& baseName, const std::filesystem::path& outputDir);

std::unique_ptr<ParseContext> setupParseContext(const bool exitOnAllErrors);

// Setup the OpmLog backends
FileOutputMode
setupLogging(Parallel::Communication& comm,
             const std::string&       deck_filename,
             const std::string&       cmdline_output_dir,
             const std::string&       cmdline_output,
             bool                     output_cout_,
             const std::string&       stdout_log_id,
             const bool               allRanksDbgLog);

/// \brief Reads the deck and creates all necessary objects if needed
///
/// If pointers already contains objects then they are used otherwise they
/// are created and can be used outside later.
void readDeck(Parallel::Communication         comm,
              const std::string&              deckFilename,
              std::shared_ptr<EclipseState>&  eclipseState,
              std::shared_ptr<Schedule>&      schedule,
              std::unique_ptr<UDQState>&      udqState,
              std::unique_ptr<Action::State>& actionState,
              std::unique_ptr<WellTestState>& wtestState,
              std::shared_ptr<SummaryConfig>& summaryConfig,
              std::shared_ptr<Python>         python,
              const std::string&              parsingStrictness,
              bool                            initFromRestart,
              bool                            checkDeck,
              const std::optional<int>&       outputInterval);

void verifyValidCellGeometry(Parallel::Communication comm,
                             const EclipseState&     eclipseState);
} // end namespace Opm

#endif // OPM_READDECK_HEADER_INCLUDED
