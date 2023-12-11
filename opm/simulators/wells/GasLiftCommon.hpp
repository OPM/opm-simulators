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

#ifndef OPM_GASLIFT_COMMON_HEADER_INCLUDED
#define OPM_GASLIFT_COMMON_HEADER_INCLUDED

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <string>

namespace Opm
{

class DeferredLogger;
class GroupState;
class WellState;

class GasLiftCommon
{
public:
    virtual ~GasLiftCommon() = default;

protected:
    GasLiftCommon(
        WellState &well_state,
        const GroupState &group_state,
        DeferredLogger &deferred_logger,
        const Parallel::Communication& comm,
        bool debug
    );
    enum class MessageType { INFO, WARNING };

    int debugUpdateGlobalCounter_() const;
    virtual void displayDebugMessage_(const std::string& msg) const = 0;
    void displayDebugMessageOnRank0_(const std::string &msg) const;
    void logMessage_(
        const std::string& prefix,
        const std::string& msg,
        MessageType msg_type = MessageType::INFO) const;

    WellState &well_state_;
    const GroupState& group_state_;
    DeferredLogger &deferred_logger_;
    const Parallel::Communication& comm_;
    bool debug;
    // By setting this variable to true we restrict some debug output
    // to only be printed for rank 0. By setting this variable to false we keep
    // the output on all ranks. This can in some cases be helpful as a debugging
    // aid to check that the output is in fact identical over all ranks
    bool debug_output_only_on_rank0 = false;
};

} // namespace Opm

#endif // OPM_GASLIFT_COMMON_INCLUDED
