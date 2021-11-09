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

#include <opm/simulators/utils/ParallelCommunication.hpp>

namespace Opm {

class EclipseState;
class Schedule;
class SummaryConfig;
class UDQState;
class WellTestState;
class TransMult;

namespace Action {
class State;
}


/*! \brief Broadcasts an eclipse state from root node in parallel runs.
 *! \param eclState EclipseState to broadcast
 *! \param schedule Schedule to broadcast
 *! \param summaryConfig SummaryConfig to broadcast
*/
void eclStateBroadcast(Parallel::Communication  comm, EclipseState& eclState, Schedule& schedule,
                       SummaryConfig& summaryConfig,
                       UDQState& udqState,
                       Action::State& actionState,
                       WellTestState& wtestState);


template <class T>
void eclBroadcast(Parallel::Communication, T& )
#if HAVE_MPI
;
#else
{}
#endif


} // end namespace Opm

#endif // PARALLEL_SERIALIZATION_HPP
