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
// The waiting is done with nonblocking collective calls on communicators of
// their own (duplicates of the simulation's communicator), one for each kind
// of call: the barrier of abortIfNotAllRanksFailed(), the barrier of
// abortIfNotAllRanksEndedRun(), and the maximum of maxOverRanksOrAbort(). MPI
// requires all ranks to make the collective calls of a communicator in the
// same order. A rank that has thrown leaves that order, so:
// - the waiting cannot share the simulation's communicator, where the other
//   ranks may be in blocking calls (MPI does not allow a nonblocking
//   collective call to be matched by a blocking one);
// - a barrier and the maximum cannot share one communicator: a rank that
//   failed a time step alone waits in the barrier, while the other ranks may
//   reach the end of a parallel try-catch region and take the maximum;
// - the two barriers cannot share one communicator either. After the barrier
//   in the failure handler of a time step, the rank retries the step; after
//   the barrier in checkAllMPIProcesses(), it ends the run. If some ranks
//   failed a time step while the others ended the run, the two barriers would
//   complete together. The ranks that retry would then wait, without a time
//   limit, for ranks that have stopped.

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <string>

namespace Opm::detail {

/// Turn the safety net on for \c comm.
///
/// Collective: all ranks of \c comm must call it, at a point where they are
/// in step. Until it is called, and after disableParallelExceptionSync(),
/// maxOverRanksOrAbort() is a plain blocking maximum without a time limit,
/// and abortIfNotAllRanksFailed() and abortIfNotAllRanksEndedRun() do
/// nothing.
///
/// \param comm The communicator of the simulation.
/// \param timeout_seconds How long a rank that holds an exception waits for
///        the other ranks. Zero or negative leaves the safety net off.
void enableParallelExceptionSync(const Parallel::Communication& comm,
                                 double timeout_seconds);

/// Turn the safety net off and free its communicators.
void disableParallelExceptionSync();

/// Whether the safety net is on for \c comm.
bool parallelExceptionSyncEnabled(const Parallel::Communication& comm);

/// Called by a rank that has caught an exception which all ranks are expected
/// to have caught: the failure handler of a time step (handleFailure in
/// AdaptiveTimeStepping::SubStepIteration::runSubStep_()), before the step is
/// retried with a smaller time step. Returns when every rank has called it.
/// Aborts the run, with \c message, when the other ranks have not done so
/// within the time limit.
///
/// \param message The text of the exception.
/// \param where Where in the program this is, for the abort message.
void abortIfNotAllRanksFailed(const Parallel::Communication& comm,
                              const std::string& message,
                              const std::string& where);

/// Called by a rank that ends the run because of an exception:
/// checkAllMPIProcesses() in FlowUtils.cpp, which FlowMain calls before it
/// returns with an error. Like abortIfNotAllRanksFailed(), but with a barrier
/// on a communicator of its own, so that it never completes together with
/// the barrier of a rank that is about to retry a time step.
///
/// \param message The text of the exception.
void abortIfNotAllRanksEndedRun(const Parallel::Communication& comm,
                                const std::string& message);

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
