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

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

namespace Opm {

GasLiftCommon::
GasLiftCommon(
    WellState &well_state,
    const GroupState &group_state,
    DeferredLogger &deferred_logger,
    const Parallel::Communication& comm,
    bool glift_debug
) :
    well_state_{well_state},
    group_state_{group_state},
    deferred_logger_{deferred_logger},
    comm_{comm},
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

void
GasLiftCommon::
displayDebugMessageOnRank0_(const std::string &msg) const
{
    // This output should be identical for all ranks.

    if (   (!this->debug_output_only_on_rank0)
        || (this->debug_output_only_on_rank0 && this->comm_.rank() == 0) ) {
        displayDebugMessage_(msg);
    }
}

void
GasLiftCommon::
logMessage_(
    const std::string& prefix, const std::string& msg, MessageType msg_type) const
{
    std::string rank = "";
    if (this->comm_.size() > 1) {
        rank = fmt::format(" Rank #{} :", this->comm_.rank());
    }
    std::string type_str;
    switch (msg_type) {
    case MessageType::INFO:
        type_str = "DEBUG";
        break;
    case MessageType::WARNING:
        type_str = "WARNING";
        break;
    default:
        throw std::runtime_error("This should not happen");
    }
    const std::string message = fmt::format(
        "  {} ({}) :{} {}", prefix, type_str, rank, msg);
    switch (msg_type) {
    case MessageType::INFO:
        this->deferred_logger_.info(message);
        break;
    case MessageType::WARNING:
        this->deferred_logger_.info(message);
        break;
    default:
        throw std::runtime_error("This should not happen");
    }
}

/****************************************
 * Private methods in alphabetical order
 ****************************************/


} // namespace Opm
