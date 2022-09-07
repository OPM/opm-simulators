/*
  Copyright 2019 Equinor AS.

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
#ifndef PARALLEL_RESTART_HPP
#define PARALLEL_RESTART_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <vector>

namespace Opm
{

class EclipseIO;
class SummaryState;
class RestartKey;
class RestartValue;

namespace Action
{
class State;
}

RestartValue loadParallelRestart(const EclipseIO* eclIO,
                                 Action::State& actionState,
                                 SummaryState& summaryState,
                                 const std::vector<RestartKey>& solutionKeys,
                                 const std::vector<RestartKey>& extraKeys,
                                 Parallel::Communication comm);

} // end namespace Opm

#endif // PARALLEL_RESTART_HPP
