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

#define BOOST_TEST_MODULE TestRstConv
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 > 66
#include <boost/test/data/test_case.hpp>
#endif

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/input/eclipse/Schedule/RSTConfig.hpp>
#include <opm/simulators/flow/RSTConv.hpp>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <random>

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

struct TestCase {
    std::size_t N;
    std::array<int,3> phase;
};

std::ostream& operator<<(std::ostream& os, const TestCase& t)
{
    os << t.N;
    for (int i : t.phase)
        os << " " << i;
    os << std::endl;

    return os;
}

static const std::vector<TestCase> tests = {
    {1, {0, 1, 2}},
    {3, {0, 1, 2}},
    {1, {0, -1, 1}},
    {5, {0, -1, 1}},
    {1, {-1, -1, 0}},
    {4, {-1, -1, 0}},
};

#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 > 66
BOOST_DATA_TEST_CASE(RstConvTest, tests)
#else
BOOST_AUTO_TEST_CASE(RstConvTest)
#endif
{
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 67
    for (const auto& sample : tests) {
#endif
    const auto& cc = Dune::MPIHelper::getCommunication();

    Opm::RSTConfig rst;
    rst.keywords["CONV"] = sample.N;

    std::vector<int> cellMapping(10);
    std::iota(cellMapping.begin(), cellMapping.end(), cc.rank()*10);

    Dune::BlockVector<Dune::FieldVector<double,3>> residual(10);

    // generate data
    std::vector<std::vector<int>> max(3);
    if (cc.rank() == 0) {
        std::random_device rng_device;
        std::mt19937 mersenne_engine{rng_device()};
        std::uniform_int_distribution<int> dist{0, 10*cc.size()-1};

        for (int c = 0; c < 3; ++c) {
            if (sample.phase[c] == -1) {
                continue;
            }
            std::vector<int>& v = max[c];
            while (v.size() < sample.N) {
                int m = dist(mersenne_engine);
                if (std::find(v.begin(), v.end(), m) == v.end()) {
                    v.push_back(m);
                }
            }
        }
    }

    for (int c = 0; c < 3; ++c) {
        std::size_t size = max[c].size();
        cc.broadcast(&size, 1, 0);
        if (cc.rank() != 0) {
            max[c].resize(size);
        }
        cc.broadcast(max[c].data(), max[c].size(), 0);
    }

    for (int i = 0; i < 10; ++i) {
        for (int c = 0; c < 3; ++c) {
            if (sample.phase[c] != -1) {
                bool inMax = std::find(max[c].begin(),
                                       max[c].end(),
                                       cellMapping[i]) != max[c].end();
                residual[i][sample.phase[c]] = inMax ? 1.0 : 1.0 / (i+2);
            }
        }
    }

    Opm::RSTConv cnv(cellMapping, cc);
    cnv.init(10*cc.size(), rst, sample.phase);

    cnv.update(residual);
    cnv.update(residual);

    if (cc.rank() == 0) {
        BOOST_CHECK_EQUAL(cnv.getData().size(), 3);
        BOOST_CHECK_EQUAL(cnv.getData()[0].size(),
                          sample.phase[0] == -1 ? 0 : cc.size() * 10);
        BOOST_CHECK_EQUAL(cnv.getData()[1].size(),
                          sample.phase[1] == -1 ? 0 : cc.size() * 10);
        BOOST_CHECK_EQUAL(cnv.getData()[2].size(),
                          sample.phase[2] == -1 ? 0 : cc.size() * 10);

        for (int i = 0; i < cc.size() * 10; ++i) {
            for (int c = 0; c < 3; ++c) {
                if (sample.phase[c] != -1) {
                    bool inMax = std::find(max[c].begin(),
                                           max[c].end(), i) != max[c].end();
                    BOOST_CHECK_EQUAL(cnv.getData()[c][i], inMax ? 2 : 0);
                }
            }
        }
    }
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 67
    }
#endif
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
