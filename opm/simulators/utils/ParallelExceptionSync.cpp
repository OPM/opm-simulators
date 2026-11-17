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
#include <config.h>

#include <opm/simulators/utils/ParallelExceptionSync.hpp>

#if HAVE_MPI
#include <mpi.h>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <thread>
#endif

namespace {

#if HAVE_MPI
// The state of the safety net in this process.
//
// The waiting is done on two duplicates of the simulation's communicator (see
// the header for why), and they must be made once, at a point where all ranks
// are in step. So they have to be kept between the calls: this struct holds
// them, together with the time limit.
//
// The fields are set by enableParallelExceptionSync(), which the simulator
// calls once at startup (FlowMain::setupParallelism()), and reset by
// disableParallelExceptionSync(), called when the simulator cleans up. Between
// the two they do not change.
//
// One instance per process is enough: in the simulator each MPI process takes
// part in one run of one model, which uses one communicator. (In a reservoir
// coupling run the master and each slave are separate groups of processes,
// each with its own communicator, so this still holds.)
struct SyncState
{
    // The simulation's communicator, as given to enableParallelExceptionSync().
    MPI_Comm base_comm = MPI_COMM_NULL;
    // Our duplicates of it: one for the barrier of abortIfNotAllRanksFailed(),
    // one for the maximum of maxOverRanksOrAbort(). Both are set, or both are
    // MPI_COMM_NULL, which means that the safety net is off.
    MPI_Comm barrier_comm = MPI_COMM_NULL;
    MPI_Comm max_comm = MPI_COMM_NULL;
    // How long a rank that holds an exception waits for the other ranks.
    double timeout_seconds = 0.0;
};

SyncState& state()
{
    static SyncState s;
    return s;
}

// Whether the safety net is on, and on for this communicator.
//
// Asked at the start of every waiting function below. The callers pass
// whatever communicator they have at hand. The safety net can only help on
// the one communicator it was enabled for: only there does it have its
// duplicates, on which all ranks make the same calls. For anything else (the
// safety net was never enabled, as in a unit test or another program that
// uses this library without calling enableParallelExceptionSync(), or the
// caller uses another communicator) the answer is no, and the caller falls
// back to a plain blocking call.
bool enabledFor(const MPI_Comm comm)
{
    const auto& s = state();
    if (s.max_comm == MPI_COMM_NULL) {
        return false;
    }
    int result = MPI_UNEQUAL;
    MPI_Comm_compare(comm, s.base_comm, &result);
    return result == MPI_IDENT;
}

[[noreturn]] void abortOnDivergentException(const int rank,
                                            const std::string& message,
                                            const std::string& where)
{
    // Written to std::cerr, not to OpmLog. This rank may be any rank, and on
    // ranks other than 0 OpmLog has no output by default. And MPI_Abort()
    // ends the process without flushing the log files.
    std::cerr << "\nError: rank " << rank << " failed, and the other ranks of this "
              << "process group did not report a failure within "
              << state().timeout_seconds << " s (" << where << ").\n"
              << "The ranks can no longer communicate in step, so the run is aborted.\n"
              << "The failure on rank " << rank << " was: " << message << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    std::abort(); // not reached
}

// Wait for a nonblocking request. With \c bounded set, give up after the time
// limit and abort the run.
void waitOrAbort(MPI_Request& request,
                 const bool bounded,
                 const int rank,
                 const std::string& message,
                 const std::string& where)
{
    if (!bounded) {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        return;
    }
    const auto start = std::chrono::steady_clock::now();
    const std::chrono::duration<double> limit{state().timeout_seconds};
    while (true) {
        int done = 0;
        MPI_Test(&request, &done, MPI_STATUS_IGNORE);
        if (done) {
            return;
        }
        if (std::chrono::steady_clock::now() - start > limit) {
            abortOnDivergentException(rank, message, where);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds{1});
    }
}
#endif // HAVE_MPI

} // anonymous namespace

namespace Opm::detail {

void abortIfNotAllRanksFailed([[maybe_unused]] const Parallel::Communication& comm,
                              [[maybe_unused]] const std::string& message,
                              [[maybe_unused]] const std::string& where)
{
#if HAVE_MPI
    if (comm.size() > 1 && enabledFor(comm)) {
        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Ibarrier(state().barrier_comm, &request);
        waitOrAbort(request, /*bounded=*/true, comm.rank(), message, where);
    }
#endif
}

void disableParallelExceptionSync()
{
#if HAVE_MPI
    auto& s = state();
    if (s.barrier_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&s.barrier_comm);
    }
    if (s.max_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&s.max_comm);
    }
    s = SyncState{};
#endif
}

void enableParallelExceptionSync([[maybe_unused]] const Parallel::Communication& comm,
                                 [[maybe_unused]] const double timeout_seconds)
{
#if HAVE_MPI
    disableParallelExceptionSync();
    if (comm.size() > 1 && timeout_seconds > 0.0) {
        auto& s = state();
        s.base_comm = comm;
        s.timeout_seconds = timeout_seconds;
        MPI_Comm_dup(comm, &s.barrier_comm);
        MPI_Comm_dup(comm, &s.max_comm);
    }
#endif
}

int maxOverRanksOrAbort(const Parallel::Communication& comm,
                        const int value,
                        [[maybe_unused]] const bool has_exception,
                        [[maybe_unused]] const std::string& message,
                        [[maybe_unused]] const std::string& where)
{
#if HAVE_MPI
    if (comm.size() > 1 && enabledFor(comm)) {
        int local = value;
        int global = value;
        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Iallreduce(&local, &global, 1, MPI_INT, MPI_MAX, state().max_comm, &request);
        waitOrAbort(request, has_exception, comm.rank(), message, where);
        return global;
    }
#endif
    return comm.max(value);
}

bool parallelExceptionSyncEnabled([[maybe_unused]] const Parallel::Communication& comm)
{
#if HAVE_MPI
    return enabledFor(comm);
#else
    return false;
#endif
}

} // namespace Opm::detail
