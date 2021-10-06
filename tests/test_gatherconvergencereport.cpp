/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Equinor.

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

#define BOOST_TEST_MODULE TestGatherConvergenceReport
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <dune/common/version.hh>
#include <opm/simulators/timestepping/gatherConvergenceReport.hpp>
#include <dune/common/parallel/mpihelper.hh>

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
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

bool operator==(const Opm::ConvergenceReport::WellFailure& wf1,
                const Opm::ConvergenceReport::WellFailure& wf2)
{
    return wf1.type() == wf2.type()
        && wf1.severity() == wf2.severity()
        && wf1.phase() == wf2.phase()
        && wf1.wellName() == wf2.wellName();
}


BOOST_AUTO_TEST_CASE(AllHaveFailure)
{
    auto cc = Dune::MPIHelper::getCollectiveCommunication();
    std::ostringstream name;
    name << "WellRank" << cc.rank() << std::flush;
    using CR = Opm::ConvergenceReport;
    CR cr;
    cr.setWellFailed({CR::WellFailure::Type::ControlBHP, CR::Severity::Normal, -1, name.str()});
    CR global_cr = gatherConvergenceReport(cr, cc);
    BOOST_CHECK(global_cr.wellFailures().size() == std::size_t(cc.size()));
    BOOST_CHECK(global_cr.wellFailures()[cc.rank()] == cr.wellFailures()[0]);
    // Extra output for debugging.
    if (cc.rank() == 0) {
        for (const auto& wf : global_cr.wellFailures()) {
            std::cout << "Well name of failure: " << wf.wellName() << std::endl;
        }
    }
}

BOOST_AUTO_TEST_CASE(EvenHaveFailure)
{
    auto cc = Dune::MPIHelper::getCollectiveCommunication();
    using CR = Opm::ConvergenceReport;
    CR cr;
    if (cc.rank() % 2 == 0) {
        std::ostringstream name;
        name << "WellRank" << cc.rank() << std::flush;
        cr.setWellFailed({CR::WellFailure::Type::ControlBHP, CR::Severity::Normal, -1, name.str()});
    }
    CR global_cr = gatherConvergenceReport(cr, cc);
    BOOST_CHECK(global_cr.wellFailures().size() == std::size_t((cc.size())+1) / 2);
    if (cc.rank() % 2 == 0) {
        BOOST_CHECK(global_cr.wellFailures()[cc.rank()/2] == cr.wellFailures()[0]);
    }
    // Extra output for debugging.
    if (cc.rank() == 0) {
        for (const auto& wf : global_cr.wellFailures()) {
            std::cout << "Well name of failure, should be only even: " << wf.wellName() << std::endl;
        }
    }
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
