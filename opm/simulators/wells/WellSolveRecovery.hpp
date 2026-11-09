/*
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

#ifndef OPM_WELL_SOLVE_RECOVERY_HPP
#define OPM_WELL_SOLVE_RECOVERY_HPP

#include <opm/common/Exceptions.hpp>

#include <string>
#include <utility>

namespace Opm {

//! \brief Run a solve of a well and recover from a numerical failure on all ranks of the well.
//!
//! A well can be distributed over several ranks, and the solve and the recovery contain collective
//! operations on the communicator of the well. If only the ranks where the solve threw a
//! NumericalProblem would run the recovery, the collective operations of the ranks would no longer
//! match. The failure is therefore agreed on across the communicator: if the solve threw on any
//! rank, the recovery is run on every rank.
//!
//! All ranks of the communicator must call this function, and the solve and the recovery must
//! handle all exceptions other than NumericalProblem themselves, as these are not agreed on.
//!
//! \param comm The communicator of the well
//! \param solve The solve, called on every rank
//! \param recover Called on every rank if the solve threw on any rank, with arguments
//!        (bool failed_on_this_rank, const std::string& reason). The reason is the message of the
//!        exception if it was thrown on this rank, otherwise it is empty.
//! \return Whether the solve failed on any rank, and the recovery was run
template<class Communication, class Solve, class Recover>
bool solveWithSharedRecovery(const Communication& comm, Solve&& solve, Recover&& recover)
{
    bool failed_on_this_rank = false;
    std::string reason;
    try {
        std::forward<Solve>(solve)();
    } catch (const NumericalProblem& e) {
        failed_on_this_rank = true;
        reason = e.what();
    }

    const bool failed = (comm.size() > 1)
        ? (comm.max(static_cast<int>(failed_on_this_rank)) != 0)
        : failed_on_this_rank;
    if (failed) {
        std::forward<Recover>(recover)(failed_on_this_rank, std::as_const(reason));
    }
    return failed;
}

} // namespace Opm

#endif // OPM_WELL_SOLVE_RECOVERY_HPP
