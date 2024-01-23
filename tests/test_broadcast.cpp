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

#define BOOST_TEST_MODULE TestBroadCast
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/utils/MPISerializer.hpp>
#include <dune/common/parallel/mpihelper.hh>

#include <numeric>

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

BOOST_AUTO_TEST_CASE(BroadCast)
{
    const auto& cc = Dune::MPIHelper::getCommunication();

    std::vector<double> d(3);
    if (cc.rank() == 1)
        std::iota(d.begin(), d.end(), 1.0);

    std::vector<int> i(3);
    if (cc.rank() == 1)
        std::iota(i.begin(), i.end(), 4);

    double d1 = cc.rank() == 1 ? 7.0 : 0.0;
    size_t i1 = cc.rank() == 1 ? 8 : 0;

    Opm::Parallel::MpiSerializer ser(cc);
    ser.broadcast(1, d, i, d1, i1);

    for (size_t c = 0; c < 3; ++c) {
        BOOST_CHECK_EQUAL(d[c], 1.0+c);
        BOOST_CHECK_EQUAL(i[c], 4+c);
    }
    BOOST_CHECK_EQUAL(d1, 7.0);
    BOOST_CHECK_EQUAL(i1, 8);
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
