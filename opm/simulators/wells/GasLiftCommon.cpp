/*
  Copyright 2022 Equinor ASA.

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
#include <opm/simulators/wells/GasLiftCommon.hpp>

namespace Opm {

GasLiftCommon::
GasLiftCommon(
    WellState &well_state,
    DeferredLogger &deferred_logger,
    bool glift_debug
) :
    well_state_{well_state},
    deferred_logger_{deferred_logger},
    debug{glift_debug}
{

}

/****************************************
 * Protected methods in alphabetical order
 ****************************************/

int
GasLiftCommon::
debugUpdateGlobalCounter_() const
{
    auto count = this->well_state_.gliftUpdateDebugCounter();
    const std::string msg = fmt::format("global counter = {}", count);
    displayDebugMessage_(msg);
    return count;
}


/****************************************
 * Private methods in alphabetical order
 ****************************************/


} // namespace Opm
