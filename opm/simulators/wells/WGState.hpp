/*
  Copyright 2021 Equinor ASA

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

#ifndef OPM_WGSTATE_HEADER_INCLUDED
#define OPM_WGSTATE_HEADER_INCLUDED

#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>

namespace Opm {

/*
  Microscopic class to handle well , group and well test state.
*/

struct PhaseUsage;

struct WGState {
    WGState(const PhaseUsage& pu);

    WellState well_state;
    GroupState group_state;
    WellTestState well_test_state;

};

}

#endif
