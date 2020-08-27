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

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Python/Python.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MessageLimits.hpp>

#include <memory>
#include <string>

namespace Opm 
{
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

void setupMessageLimiter(const Opm::MessageLimits msgLimits,  const std::string& stdout_log_id);


/// \brief Reads the deck and creates all necessary objects if needed
///
/// If pointers already contains objects then they are used otherwise they are created and can be used outside later.
void readDeck(int rank, std::string& deckFilename, std::unique_ptr<Opm::Deck>& deck, std::unique_ptr<Opm::EclipseState>& eclipseState,
              std::unique_ptr<Opm::Schedule>& schedule, std::unique_ptr<Opm::SummaryConfig>& summaryConfig,
              std::unique_ptr<ErrorGuard> errorGuard, std::shared_ptr<Opm::Python>& python, std::unique_ptr<ParseContext> parseContext,
              bool initFromRestart, bool checkDeck);
} // end namespace Opm

#endif // OPM_READDECK_HEADER_INCLUDED
