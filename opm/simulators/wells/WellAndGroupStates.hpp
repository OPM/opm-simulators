/*
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_WELLANDGROUPSTATES_HEADER_INCLUDED
#define OPM_WELLANDGROUPSTATES_HEADER_INCLUDED

#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/SingleGroupState.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm
{

/// Encapsulate all data needed to represent the state of all wells
/// and groups for persistence purposes, i.e. from one timestep to the
/// next or for restarting, and for accumulating rates to get correct
/// group rate limits.
/// In a parallel context, the well states will be those for local
/// wells only, while the group information will be global (and
/// communicated as needed).
template <int NumActiveComponents, int NumActivePhases>
class WellAndGroupStates
{
public:
private:
    using WellState = SingleWellState<NumActiveComponents, NumActivePhases>;
    using GroupState = SingleGroupState<NumActiveComponents, NumActivePhases>;
    std::vector<WellState> well_states_;
    std::map<std::string, GroupState> group_states_;
};

}

#endif // OPM_WELLANDGROUPSTATES_HEADER_INCLUDED
