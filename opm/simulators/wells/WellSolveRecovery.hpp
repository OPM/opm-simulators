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

#include <exception>
#include <stdexcept>
#include <string>
#include <utility>

namespace Opm {

//! \brief Whether a condition holds on any rank of the communicator.
//!
//! Used to take decisions that are followed by collective operations on all ranks of a well.
//! All ranks of the communicator must call this function.
template<class Communication>
bool anyRank(const Communication& comm, const bool condition)
{
    return (comm.size() > 1) ? (comm.max(static_cast<int>(condition)) != 0) : condition;
}

//! \brief Run a solve of a well and recover from a numerical failure on all ranks of the well.
//!
//! A well can be distributed over several ranks, and the solve and the recovery contain collective
//! operations on the communicator of the well. If only the ranks where the solve threw would run the
//! recovery, or leave the function, the collective operations of the ranks would no longer match.
//! The outcome is therefore agreed on across the communicator:
//! - if the solve threw a NumericalProblem on any rank, the recovery is run on every rank.
//! - if the solve threw any other exception on any rank, the exception is thrown on that rank and
//!   a std::runtime_error on the other ranks. The recovery is not run, as the state of the solve is
//!   not known. This is a failure that is handled by the caller, for instance by failing the time step.
//!
//! All ranks of the communicator must call this function. This only covers failures of the solve itself,
//! a rank that throws while the other ranks are inside a collective operation of the solve is not helped.
//!
//! \param comm The communicator of the well
//! \param solve The solve, called on every rank
//! \param recover Called on every rank if the solve threw a NumericalProblem on any rank, with
//!        arguments (bool failed_on_this_rank, const std::string& reason). The reason is the message
//!        of the exception if it was thrown on this rank, otherwise it is empty.
//! \return Whether the solve failed with a NumericalProblem on any rank, and the recovery was run
template<class Communication, class Solve, class Recover>
bool solveWithSharedRecovery(const Communication& comm, Solve&& solve, Recover&& recover)
{
    constexpr int ok = 0;
    constexpr int numerical_problem = 1;
    constexpr int fatal = 2;

    int outcome = ok;
    std::string reason;
    std::exception_ptr fatal_exception;
    try {
        std::forward<Solve>(solve)();
    } catch (const NumericalProblem& e) {
        outcome = numerical_problem;
        reason = e.what();
    } catch (...) {
        outcome = fatal;
        fatal_exception = std::current_exception();
    }

    // the most severe outcome on any rank
    const int agreed = (comm.size() > 1) ? comm.max(outcome) : outcome;
    if (agreed == fatal) {
        if (fatal_exception) {
            std::rethrow_exception(fatal_exception);
        }
        throw std::runtime_error("The solve of a well failed on another rank of the well.");
    }
    if (agreed == numerical_problem) {
        std::forward<Recover>(recover)(outcome == numerical_problem, std::as_const(reason));
        return true;
    }
    return false;
}

//! \brief Solve a list of wells, and recover from numerical failures, with the failures of the wells deferred.
//!
//! The ranks can have different lists of wells, for instance a well on one rank only and a well
//! distributed over two ranks. The collective operations of a well must then be reached by all ranks
//! of the well, also if an earlier well failed on this rank. An exception from the solve of a well is
//! therefore kept until all wells have been solved, and then thrown. The first exception is thrown.
//! The ranks that did not fail reach the end of the loop normally, so the code that handles the
//! exception collectively, such as OPM_END_PARALLEL_TRY_CATCH, can follow the call.
//!
//! \param wells The wells of this rank, in the same order for the ranks of a well
//! \param get_communication Called with a well, gives the communicator of the well
//! \param solve Called with a well, solves it
//! \param recover Called with a well, a bool and a string, as in solveWithSharedRecovery()
//! \return Whether any well failed with a NumericalProblem on any of its ranks
template<class Wells, class GetCommunication, class Solve, class Recover>
bool solveWellsWithSharedRecovery(Wells&& wells,
                                  GetCommunication&& get_communication,
                                  Solve&& solve,
                                  Recover&& recover)
{
    bool any_recovered = false;
    std::exception_ptr deferred;
    for (auto& well : wells) {
        try {
            const bool recovered = solveWithSharedRecovery(
                get_communication(well),
                [&]() { solve(well); },
                [&](const bool failed_on_this_rank, const std::string& reason) {
                    recover(well, failed_on_this_rank, reason);
                });
            any_recovered = any_recovered || recovered;
        } catch (...) {
            if (!deferred) {
                deferred = std::current_exception();
            }
        }
    }
    if (deferred) {
        std::rethrow_exception(deferred);
    }
    return any_recovered;
}

} // namespace Opm

#endif // OPM_WELL_SOLVE_RECOVERY_HPP
