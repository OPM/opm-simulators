/*
  Copyright 2020 Equinor AS.

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
#ifndef PARALLEL_SERIALIZATION_HPP
#define PARALLEL_SERIALIZATION_HPP

#include <dune/common/parallel/mpihelper.hh>

namespace Opm {

class EclipseState;
class Schedule;
class SummaryConfig;

using CollCommType = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>;
/*! \brief Broadcasts an eclipse state from root node in parallel runs.
 *! \param eclState EclipseState to broadcast
 *! \param schedule Schedule to broadcast
 *! \param summaryConfig SummaryConfig to broadcast
*/
void eclStateBroadcast(CollCommType comm, EclipseState& eclState, Schedule& schedule,
                       SummaryConfig& summaryConfig);

/// \brief Broadcasts an schedule from root node in parallel runs.
void eclScheduleBroadcast(CollCommType comm, Schedule& schedule);

} // end namespace Opm

#endif // PARALLEL_SERIALIZATION_HPP
