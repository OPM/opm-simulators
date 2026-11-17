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
#ifndef OPM_PARALLEL_EXCEPTION_SYNC_HPP
#define OPM_PARALLEL_EXCEPTION_SYNC_HPP

// A safety net for exceptions that only some MPI ranks throw.
//
// Every collective MPI call must be made by all ranks, in the same order. A
// rank that throws leaves that order. If the other ranks do not throw too,
// they wait in one collective call and the throwing rank waits in another,
// forever. The functions below let a rank that holds an exception wait for the
// other ranks for a bounded time only, and abort the run with the exception's
// message when they do not show up.
//
// The waiting is done with nonblocking collective calls on a communicator of
// its own (a duplicate of the simulation's communicator). On that communicator
// all ranks make the same calls in the same order, also when one rank has left
// the order on the simulation's communicator. MPI does not allow a nonblocking
// collective call to be matched by a blocking one, so the two kinds cannot
// share a communicator here.

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <string>

namespace Opm::detail {

/// Turn the safety net on for \c comm.
///
/// Collective: all ranks of \c comm must call it, at a point where they are
/// in step. Until it is called, and after disableParallelExceptionSync(), the
/// functions below behave as plain blocking calls without a time limit.
///
/// \param comm The communicator of the simulation.
/// \param timeout_seconds How long a rank that holds an exception waits for
///        the other ranks. Zero or negative leaves the safety net off.
void enableParallelExceptionSync(const Parallel::Communication& comm,
                                 double timeout_seconds);

/// Turn the safety net off and free its communicator.
void disableParallelExceptionSync();

/// Whether the safety net is on for \c comm.
bool parallelExceptionSyncEnabled(const Parallel::Communication& comm);

/// Called by a rank that has caught an exception which all ranks are expected
/// to have caught. Returns when every rank has called it. Aborts the run, with
/// \c message, when the other ranks have not done so within the time limit.
///
/// \param message The text of the exception.
/// \param where Where in the program this is, for the abort message.
void abortIfNotAllRanksFailed(const Parallel::Communication& comm,
                              const std::string& message,
                              const std::string& where);

/// The maximum of \c value over all ranks; all ranks must call it.
///
/// A rank that holds an exception (\c has_exception) waits for the other
/// ranks for the time limit only, and then aborts the run with \c message.
/// A rank without an exception waits without a limit.
int maxOverRanksOrAbort(const Parallel::Communication& comm,
                        int value,
                        bool has_exception,
                        const std::string& message,
                        const std::string& where);

} // namespace Opm::detail

#endif // OPM_PARALLEL_EXCEPTION_SYNC_HPP
