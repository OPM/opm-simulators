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

#define BOOST_TEST_MODULE TestWellSolveRecovery
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/simulators/wells/WellSolveRecovery.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <string>
#include <utility>
#include <vector>

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e) {}
    std::string errorstring;
    int errorcode;
};

void MPI_err_handler(MPI_Comm*, int* err_code, ...)
{
    std::vector<char> err_string(MPI_MAX_ERROR_STRING);
    int err_length;
    MPI_Error_string(*err_code, err_string.data(), &err_length);
    std::string s(err_string.data(), err_length);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *err_code);
}
#endif

bool
init_unit_test_func()
{
    return true;
}

namespace {

// What the recovery observed on this rank.
struct Outcome
{
    int recoveries = 0;
    bool failed_on_this_rank = false;
    std::string reason;
    // result of a collective operation done in the recovery
    int ranks_in_recovery = 0;
};

// The recovery of a well contains collective operations, so all ranks must run it.
template<class Comm, class Solve>
std::pair<bool, Outcome> run(const Comm& comm, Solve&& solve)
{
    Outcome outcome;
    const bool failed = Opm::solveWithSharedRecovery(
        comm,
        std::forward<Solve>(solve),
        [&](const bool failed_on_this_rank, const std::string& reason) {
            ++outcome.recoveries;
            outcome.failed_on_this_rank = failed_on_this_rank;
            outcome.reason = reason;
            outcome.ranks_in_recovery = comm.sum(1);
        });
    return {failed, outcome};
}

} // namespace

BOOST_AUTO_TEST_CASE(NoFailureOnAnyRank)
{
    const auto& comm = Dune::MPIHelper::getCommunication();
    BOOST_REQUIRE_GE(comm.size(), 2);

    int solves = 0;
    const auto [failed, outcome] = run(comm, [&]() { ++solves; });

    BOOST_CHECK(!failed);
    BOOST_CHECK_EQUAL(solves, 1);
    BOOST_CHECK_EQUAL(outcome.recoveries, 0);
}

BOOST_AUTO_TEST_CASE(FailureOnOneRankRecoversOnAllRanks)
{
    const auto& comm = Dune::MPIHelper::getCommunication();
    BOOST_REQUIRE_GE(comm.size(), 2);

    // The solve throws on the last rank only.
    const bool throws_here = comm.rank() == comm.size() - 1;
    const auto [failed, outcome] = run(comm, [&]() {
        if (throws_here) {
            OPM_THROW(Opm::NumericalProblem, "singular matrix");
        }
    });

    BOOST_CHECK(failed);
    BOOST_CHECK_EQUAL(outcome.recoveries, 1);
    BOOST_CHECK_EQUAL(outcome.failed_on_this_rank, throws_here);
    if (throws_here) {
        BOOST_CHECK(outcome.reason.find("singular matrix") != std::string::npos);
    } else {
        BOOST_CHECK(outcome.reason.empty());
    }
    // all ranks took part in the collective operation of the recovery
    BOOST_CHECK_EQUAL(outcome.ranks_in_recovery, comm.size());
}

BOOST_AUTO_TEST_CASE(FailureOnFirstRankRecoversOnAllRanks)
{
    const auto& comm = Dune::MPIHelper::getCommunication();
    BOOST_REQUIRE_GE(comm.size(), 2);

    const bool throws_here = comm.rank() == 0;
    const auto [failed, outcome] = run(comm, [&]() {
        if (throws_here) {
            OPM_THROW(Opm::NumericalProblem, "nan in residual");
        }
    });

    BOOST_CHECK(failed);
    BOOST_CHECK_EQUAL(outcome.recoveries, 1);
    BOOST_CHECK_EQUAL(outcome.failed_on_this_rank, throws_here);
    BOOST_CHECK_EQUAL(outcome.ranks_in_recovery, comm.size());
}

BOOST_AUTO_TEST_CASE(FailureOnAllRanks)
{
    const auto& comm = Dune::MPIHelper::getCommunication();
    BOOST_REQUIRE_GE(comm.size(), 2);

    const auto [failed, outcome] = run(comm, []() {
        OPM_THROW(Opm::NumericalProblem, "failed everywhere");
    });

    BOOST_CHECK(failed);
    BOOST_CHECK_EQUAL(outcome.recoveries, 1);
    BOOST_CHECK(outcome.failed_on_this_rank);
    BOOST_CHECK_EQUAL(outcome.ranks_in_recovery, comm.size());
}

BOOST_AUTO_TEST_CASE(SolveIsRepeatableAfterRecovery)
{
    // A well is solved again in later iterations, so the helper keeps no state between calls.
    const auto& comm = Dune::MPIHelper::getCommunication();
    BOOST_REQUIRE_GE(comm.size(), 2);

    const bool throws_here = comm.rank() == 0;
    const auto first = run(comm, [&]() {
        if (throws_here) {
            OPM_THROW(Opm::NumericalProblem, "first solve");
        }
    });
    const auto second = run(comm, []() {});

    BOOST_CHECK(first.first);
    BOOST_CHECK(!second.first);
    BOOST_CHECK_EQUAL(second.second.recoveries, 0);
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    // Throw an exception instead of aborting on MPI errors
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
