/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2021 Equinor

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

#define BOOST_TEST_MODULE OPM_test_openclSolver
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6) && \
    BOOST_VERSION / 100 % 1000 > 48

#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

class PlatformInitException : public std::logic_error
{
public:
    PlatformInitException(std::string msg) : logic_error(msg){};
};

template <int bz>
Dune::BlockVector<Dune::FieldVector<double, bz>>
testOpenclSolver(const boost::property_tree::ptree& prm, const std::string& matrix_filename, const std::string& rhs_filename)
{
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;
    Matrix matrix;
    {
        std::ifstream mfile(matrix_filename);
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        readMatrixMarket(matrix, mfile);
    }
    Vector rhs;
    {
        std::ifstream rhsfile(rhs_filename);
        if (!rhsfile) {
            throw std::runtime_error("Could not read rhs file");
        }
        readMatrixMarket(rhs, rhsfile);
    }

    const int linear_solver_verbosity = prm.get<int>("verbosity");
    const int maxit = prm.get<int>("maxiter");
    const double tolerance = prm.get<double>("tol");
    const std::string opencl_ilu_reorder("none");
    const int platformID = 0;
    const int deviceID = 0;
    const std::string gpu_mode("opencl");
    const std::string fpga_bitstream("empty");    // unused
    Dune::InverseOperatorResult result;

    Vector x(rhs.size());
    auto wellContribs = Opm::WellContributions::create("opencl", false);
    std::unique_ptr<Opm::BdaBridge<Matrix, Vector, bz> > bridge;
    try {
        bridge = std::make_unique<Opm::BdaBridge<Matrix, Vector, bz> >(gpu_mode, fpga_bitstream, linear_solver_verbosity, maxit, tolerance, platformID, deviceID, opencl_ilu_reorder);
    } catch (const std::logic_error& error) {
        BOOST_WARN_MESSAGE(true, error.what());
        throw PlatformInitException(error.what());
    }
    bridge->solve_system(&matrix, rhs, *wellContribs, result);
    bridge->get_result(x);

    return x;
}

namespace pt = boost::property_tree;

void test3(const pt::ptree& prm)
{
    const int bz = 3;
    auto sol = testOpenclSolver<bz>(prm, "matr33.txt", "rhs3.txt");
    Dune::BlockVector<Dune::FieldVector<double, bz>> expected {{-1.30307e-2, -3.58263e-6, 1.13836e-9},
            {-1.25425e-3, -1.4167e-4, -3.2213e-3},
                {-4.5436e-4, 1.28682e-5, 4.7644e-6}};
    BOOST_REQUIRE_EQUAL(sol.size(), expected.size());
    for (size_t i = 0; i < sol.size(); ++i) {
        for (int row = 0; row < bz; ++row) {
            BOOST_CHECK_CLOSE(sol[i][row], expected[i][row], 1e-3);
        }
    }
}


BOOST_AUTO_TEST_CASE(TestDefaultPreconditionerFactory)
{
    pt::ptree prm;

    // Read parameters.
    {
        std::ifstream file("options_flexiblesolver.json");
        pt::read_json(file, prm);
    }

    try {
        // Test with 3x3 block solvers.
        test3(prm);
    } catch(const PlatformInitException& ) {
        BOOST_WARN_MESSAGE(true, "Problem with initializing Platform. skipping test");
    }
}


#else

// Do nothing if we do not have at least Dune 2.6.
BOOST_AUTO_TEST_CASE(DummyTest)
{
    BOOST_REQUIRE(true);
}

#endif
