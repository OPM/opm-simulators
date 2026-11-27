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

#define BOOST_TEST_MODULE TestParallelExceptionSync
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/FlowUtils.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/utils/ParallelExceptionSync.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <chrono>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

Opm::Parallel::Communication comm()
{
    return Opm::Parallel::Communication{Dune::MPIHelper::getCommunicator()};
}

// Turns the safety net on for the lifetime of the object.
struct Enabled
{
    explicit Enabled(const double timeout_seconds)
    {
        Opm::detail::enableParallelExceptionSync(comm(), timeout_seconds);
    }

    ~Enabled()
    {
        Opm::detail::disableParallelExceptionSync();
    }
};

bool init_unit_test_func()
{
    return true;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(OffByDefault)
{
    const auto cc = comm();
    BOOST_CHECK(!Opm::detail::parallelExceptionSyncEnabled(cc));

    // Plain blocking behavior: the maximum over the ranks, and no waiting.
    const int result = Opm::detail::maxOverRanksOrAbort(cc,
                                                        /*value=*/cc.rank(),
                                                        /*has_exception=*/cc.rank() == 0,
                                                        /*message=*/"msg",
                                                        /*where=*/"test");
    BOOST_CHECK_EQUAL(result, cc.size() - 1);
    Opm::detail::abortIfNotAllRanksFailed(cc, /*message=*/"msg", /*where=*/"test");
}

BOOST_AUTO_TEST_CASE(ZeroTimeoutLeavesItOff)
{
    const Enabled enabled{/*timeout_seconds=*/0.0};
    BOOST_CHECK(!Opm::detail::parallelExceptionSyncEnabled(comm()));
}

BOOST_AUTO_TEST_CASE(MaxOverRanks)
{
    const auto cc = comm();
    const Enabled enabled{/*timeout_seconds=*/5.0};
    BOOST_CHECK_EQUAL(Opm::detail::parallelExceptionSyncEnabled(cc), cc.size() > 1);

    // No rank holds an exception.
    BOOST_CHECK_EQUAL(Opm::detail::maxOverRanksOrAbort(cc,
                                                       /*value=*/0,
                                                       /*has_exception=*/false,
                                                       /*message=*/"msg",
                                                       /*where=*/"test"),
                      0);

    // Rank 0 alone holds one. All ranks take part, so nothing is aborted and
    // every rank learns about it.
    const bool failed = cc.rank() == 0;
    BOOST_CHECK_EQUAL(Opm::detail::maxOverRanksOrAbort(cc,
                                                       /*value=*/failed ? 3 : 0,
                                                       /*has_exception=*/failed,
                                                       /*message=*/"msg",
                                                       /*where=*/"test"),
                      3);
}

BOOST_AUTO_TEST_CASE(AllRanksFailed)
{
    const auto cc = comm();
    const Enabled enabled{/*timeout_seconds=*/5.0};

    // The ranks arrive at different times, all well within the time limit.
    std::this_thread::sleep_for(std::chrono::milliseconds{50 * cc.rank()});
    Opm::detail::abortIfNotAllRanksFailed(cc, /*message=*/"msg", /*where=*/"test");
    // In step again afterwards.
    BOOST_CHECK_EQUAL(cc.sum(1), cc.size());
}

// A "parallel try-catch region" is the code between the macros
// OPM_BEGIN_PARALLEL_TRY_CATCH and OPM_END_PARALLEL_TRY_CATCH: an exception
// thrown in it on any rank is thrown again on every rank at its end.
BOOST_AUTO_TEST_CASE(ParallelTryCatchRegion)
{
    const auto cc = comm();
    const Enabled enabled{/*timeout_seconds=*/5.0};

    // An exception on the last rank only is rethrown on every rank.
    const auto region = [&cc]() {
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        if (cc.rank() == cc.size() - 1) {
            throw std::runtime_error("failure on the last rank");
        }
        OPM_END_PARALLEL_TRY_CATCH(/*prefix=*/"region: ", cc);
    };
    BOOST_CHECK_THROW(region(), std::runtime_error);

    // The same with the variant that also gathers log messages.
    const auto logged_region = [&cc]() {
        Opm::DeferredLogger logger;
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        if (cc.rank() == cc.size() - 1) {
            throw std::runtime_error("failure on the last rank");
        }
        OPM_END_PARALLEL_TRY_CATCH_LOG(logger, /*prefix=*/"region: ", /*terminal_output=*/false, cc);
    };
    BOOST_CHECK_THROW(logged_region(), std::runtime_error);

    // In step again afterwards.
    BOOST_CHECK_EQUAL(cc.sum(1), cc.size());
}

// Rank 0 fails alone: the other ranks never report a failure, so rank 0 must
// abort the run after the time limit. The run ends in MPI_Abort(), which ends
// the whole test program: the test cases after it would not run, and the
// program would report a failure. So this case is disabled in the normal run
// of the test, and has a ctest entry of its own, which checks the abort
// message. ("* disabled()" is how Boost.Test
// attaches a decorator to a test case: the case is skipped unless it is asked
// for by name with --run_test.)
BOOST_AUTO_TEST_CASE(LoneFailureAborts, * boost::unit_test::disabled())
{
    const auto cc = comm();
    BOOST_REQUIRE(cc.size() > 1);
    const Enabled enabled{/*timeout_seconds=*/1.0};

    if (cc.rank() == 0) {
        Opm::detail::abortIfNotAllRanksFailed(cc,
                                              /*message=*/"the message of the exception",
                                              /*where=*/"in the unit test");
    }
    else {
        // Stay away from the agreement for much longer than rank 0's time
        // limit of 1 s. The abort by rank 0 ends this sleep early.
        std::this_thread::sleep_for(std::chrono::seconds{20});
    }
    BOOST_ERROR("The run should have been aborted");
}

// Rank 0 fails a time step alone and waits in the barrier of
// abortIfNotAllRanksFailed(), while the other ranks go on and reach the end
// of a parallel try-catch region, where maxOverRanksOrAbort() takes the
// maximum. The two calls must not be matched with each other: rank 0 must
// abort after its time limit, while the other ranks (which have no exception,
// and so no time limit) wait until the abort ends them. Disabled in the normal
// run for the same reason as LoneFailureAborts, with a ctest entry of its own.
BOOST_AUTO_TEST_CASE(BarrierMeetsRegionEnd, * boost::unit_test::disabled())
{
    const auto cc = comm();
    BOOST_REQUIRE(cc.size() > 1);
    const Enabled enabled{/*timeout_seconds=*/1.0};

    if (cc.rank() == 0) {
        Opm::detail::abortIfNotAllRanksFailed(cc,
                                              /*message=*/"the message of the exception",
                                              /*where=*/"when a time step failed");
    }
    else {
        Opm::detail::maxOverRanksOrAbort(cc,
                                         /*value=*/0,
                                         /*has_exception=*/false,
                                         /*message=*/"",
                                         /*where=*/"at the end of a region");
    }
    BOOST_ERROR("The run should have been aborted");
}

// Rank 0 throws inside a parallel try-catch region, while the other ranks
// wait in a blocking collective call further down in the same region, on the
// simulation's communicator. Rank 0 never makes that call, so the other ranks
// never reach the end of the region. Rank 0 must abort at the end of the
// region after its time limit, with the message of its exception. Disabled in
// the normal run for the same reason as LoneFailureAborts, with a ctest entry
// of its own.
BOOST_AUTO_TEST_CASE(RegionEndFailureAborts, * boost::unit_test::disabled())
{
    const auto cc = comm();
    BOOST_REQUIRE(cc.size() > 1);
    const Enabled enabled{/*timeout_seconds=*/1.0};

    OPM_BEGIN_PARALLEL_TRY_CATCH();
    if (cc.rank() == 0) {
        throw std::runtime_error("failure on rank 0 only");
    }
    // Without rank 0, this call never returns.
    [[maybe_unused]] const int total = cc.sum(1);
    OPM_END_PARALLEL_TRY_CATCH(/*prefix=*/"region: ", cc);

    BOOST_ERROR("The run should have been aborted");
}

// The last rank alone ends the run with an exception, and checks the other
// ranks in checkAllMPIProcesses(), as FlowMain does. The other ranks never
// get there. The last rank must abort after its time limit, with the message
// of its exception. On a rank other than 0 this needs the safety net: the
// check without it sends one int to rank 0, and MPI reports a send that small
// as complete as soon as it has taken the data, even if rank 0 never posts the
// receive. Disabled in the normal run for the same reason as
// LoneFailureAborts, with a ctest entry of its own.
BOOST_AUTO_TEST_CASE(LoneFailureOnLastRankAborts, * boost::unit_test::disabled())
{
    const auto cc = comm();
    BOOST_REQUIRE(cc.size() > 1);
    Opm::FlowGenericVanguard::setCommunication(
        std::make_unique<Opm::Parallel::Communication>(cc));
    const Enabled enabled{/*timeout_seconds=*/1.0};

    if (cc.rank() == cc.size() - 1) {
        Opm::detail::checkAllMPIProcesses("failure on the last rank only");
    }
    else {
        // Much longer than the time limit of 1 s, and than the 9 s of the
        // check without the safety net. The abort ends this sleep early.
        std::this_thread::sleep_for(std::chrono::seconds{20});
    }
    BOOST_ERROR("The run should have been aborted");
}

// Rank 0 ends the run with an exception and checks the other ranks in
// checkAllMPIProcesses(), while the other ranks failed a time step and wait
// in the barrier of abortIfNotAllRanksFailed() before they retry it. The two
// waits must not be matched with each other: after them, rank 0 would stop
// and the other ranks would go on, out of step. Both sides must abort after
// the time limit instead. Disabled in the normal run for the same reason as
// LoneFailureAborts, with a ctest entry of its own.
BOOST_AUTO_TEST_CASE(BarrierMeetsRunEnd, * boost::unit_test::disabled())
{
    const auto cc = comm();
    BOOST_REQUIRE(cc.size() > 1);
    Opm::FlowGenericVanguard::setCommunication(
        std::make_unique<Opm::Parallel::Communication>(cc));
    const Enabled enabled{/*timeout_seconds=*/1.0};

    if (cc.rank() == 0) {
        Opm::detail::checkAllMPIProcesses("the run ended on rank 0");
    }
    else {
        Opm::detail::abortIfNotAllRanksFailed(cc,
                                              /*message=*/"the time step failed",
                                              /*where=*/"when a time step failed");
    }
    BOOST_ERROR("The run should have been aborted");
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
