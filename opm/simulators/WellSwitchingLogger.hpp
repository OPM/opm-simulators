/*
  Copyright 2016 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 Statoil AS

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

#ifndef OPM_WELLSWITCHINGLOGGER_HEADER_INCLUDED
#define OPM_WELLSWITCHINGLOGGER_HEADER_INCLUDED

#include <array>
#include <map>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/core/well_controls.h>

namespace Opm
{
namespace wellhelpers
{

/// \brief Utility class to handle the log messages about well switching.
///
/// In parallel all the messages will be send to a root processor
/// and logged there.
class WellSwitchingLogger
{
    typedef std::map<std::string, std::array<char,2> > SwitchMap;

public:
    /// \brief The type of the collective communication used.
    typedef Dune::CollectiveCommunication<typename Dune::MPIHelper::MPICommunicator>
    Communication;

    /// \brief Constructor.
    ///
    /// \param cc The collective communication to use.
    explicit WellSwitchingLogger(const Communication& cc =
                                 Dune::MPIHelper::getCollectiveCommunication())
        : cc_(cc)
    {}

    /// \brief Log that a well switched.
    /// \param name The name of the well.
    /// \param from The control of the well before the switch.
    /// \param to The control of the well after the switch.
    void wellSwitched(std::string name,
                      WellControlType from,
                      WellControlType to)
    {
        if( cc_.size() > 1 )
        {
            using Pair = typename SwitchMap::value_type;
            switchMap_.insert(Pair(name, {{char(from), char(to)}}));
        }
        else
        {
            std::ostringstream ss;
            ss << "    Switching control mode for well " << name
               << " from " << modestring[from]
               << " to " <<  modestring[to];
            OpmLog::info(ss.str());
        }
    }

    /// \brief Destructor send does the actual logging.
    ~WellSwitchingLogger();

private:

#if HAVE_MPI
    void unpackDataAndLog(std::vector<char>& recv_buffer,
                          const std::vector<int>& displ);

    void packData(std::vector<int>& well_name_length,
                  std::vector<char>& buffer);

    int calculateMessageSize(std::vector<int>& well_name_length);

    void logSwitch(const char* name, std::array<char,2> fromto,
                   int rank);

#endif // HAVE_MPI

    void gatherDataAndLog();
    
    /// \brief A map containing the local switches
    SwitchMap switchMap_;
    /// \brief Collective communication object.
    Communication cc_;
    /// \brief The strings for printing.
    const std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };
};
} // end namespace wellhelpers
} // end namespace Opm
#endif
