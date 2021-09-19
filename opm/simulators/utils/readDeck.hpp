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

#include <memory>
#include <optional>
#include <string>

namespace Opm 
{

class Deck;
class EclipseState;
class ErrorGuard;
class ParseContext;
class Python;
class Schedule;
class SummaryConfig;
class UDQState;

namespace Action {
class State;
}

enum class FileOutputMode {
    //! \brief No output to files.
    OUTPUT_NONE = 0,
    //! \brief Output only to log files, no eclipse output.
    OUTPUT_LOG_ONLY = 1,
    //! \brief Output to all files.
    OUTPUT_ALL = 3
};

// Setup the OpmLog backends
FileOutputMode setupLogging(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id);

/// \brief Reads the deck and creates all necessary objects if needed
///
/// If pointers already contains objects then they are used otherwise they are created and can be used outside later.
void readDeck(int rank, std::string& deckFilename, std::unique_ptr<Deck>& deck, std::unique_ptr<EclipseState>& eclipseState,
              std::unique_ptr<Schedule>& schedule, std::unique_ptr<UDQState>& udqState, std::unique_ptr<Action::State>& actionState, std::unique_ptr<SummaryConfig>& summaryConfig,
              std::unique_ptr<ErrorGuard> errorGuard, std::shared_ptr<Python>& python, std::unique_ptr<ParseContext> parseContext,
              bool initFromRestart, bool checkDeck, const std::optional<int>& outputInterval);
} // end namespace Opm

#endif // OPM_READDECK_HEADER_INCLUDED
