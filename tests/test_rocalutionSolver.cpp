/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2022 Equinor

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

#define BOOST_TEST_MODULE OPM_test_rocalutionSolver
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#include <rocalution.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

template <int bz>
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
template <int bz>
using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

template <int bz>
void readLinearSystem(const std::string& matrix_filename, const std::string& rhs_filename, Matrix<bz>& matrix, Vector<bz>& rhs)
{
    {
        std::ifstream mfile(matrix_filename);
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        readMatrixMarket(matrix, mfile);
    }
    {
        std::ifstream rhsfile(rhs_filename);
        if (!rhsfile) {
            throw std::runtime_error("Could not read rhs file");
        }
        readMatrixMarket(rhs, rhsfile);
    }
}

template <int bz>
Dune::BlockVector<Dune::FieldVector<double, bz>>
getDuneSolution(Matrix<bz>& matrix, Vector<bz>& rhs)
{
    Dune::InverseOperatorResult result;

    Vector<bz> x(rhs.size());

    typedef Dune::MatrixAdapter<Matrix<bz>,Vector<bz>,Vector<bz> > Operator;
    Operator fop(matrix);
    double relaxation = 0.9;
    Dune::SeqILU<Matrix<bz>,Vector<bz>,Vector<bz> > prec(matrix, relaxation);
    double reduction = 1e-2;
    int maxit = 10;
    int verbosity = 0;
    Dune::BiCGSTABSolver<Vector<bz> > solver(fop, prec, reduction, maxit, verbosity);
    solver.apply(x, rhs, result);
    return x;
}

template <int bz>
Dune::BlockVector<Dune::FieldVector<double, bz>>
testRocalutionSolver(const boost::property_tree::ptree& prm, Matrix<bz>& matrix, Vector<bz>& rhs)
{
    const int linear_solver_verbosity = prm.get<int>("verbosity");
    const int maxit = prm.get<int>("maxiter");
    const double tolerance = prm.get<double>("tol");
    const bool opencl_ilu_parallel(true);
    const int platformID = 0;
    const int deviceID = 0;
    const std::string accelerator_mode("rocalution");
    const std::string linsolver("ilu0");
    Dune::InverseOperatorResult result;

    Vector<bz> x(rhs.size());
    auto wellContribs = Opm::WellContributions::create(accelerator_mode, true);
    std::unique_ptr<Opm::BdaBridge<Matrix<bz>, Vector<bz>, bz> > bridge;
    try {
        bridge = std::make_unique<Opm::BdaBridge<Matrix<bz>, Vector<bz>, bz> >(accelerator_mode,
                                                                               linear_solver_verbosity,
                                                                               maxit,
                                                                               tolerance,
                                                                               platformID,
                                                                               deviceID,
                                                                               opencl_ilu_parallel,
                                                                               linsolver);
    } catch (const std::logic_error& error) {
        BOOST_WARN_MESSAGE(true, error.what());
    }
    auto mat2 = matrix; // deep copy to make sure nnz values are in contiguous memory
                        // matrix created by readMatrixMarket() did not have contiguous memory
    bridge->solve_system(&mat2, &mat2, /*numJacobiBlocks=*/0, rhs, *wellContribs, result);
    bridge->get_result(x);

    return x;
}

namespace pt = boost::property_tree;

void test3(const pt::ptree& prm)
{
    const int bz = 3;
    Matrix<bz> matrix;
    Vector<bz> rhs;
    readLinearSystem("matr33.txt", "rhs3.txt", matrix, rhs);
    Vector<bz> rhs2 = rhs; // deep copy, getDuneSolution() changes values in rhs vector
    auto duneSolution = getDuneSolution<bz>(matrix, rhs);
    auto sol = testRocalutionSolver<bz>(prm, matrix, rhs2);

    BOOST_REQUIRE_EQUAL(sol.size(), duneSolution.size());
    for (size_t i = 0; i < sol.size(); ++i) {
        for (int row = 0; row < bz; ++row) {
            BOOST_CHECK_CLOSE(sol[i][row], duneSolution[i][row], 1e-3);
        }
    }
}


BOOST_AUTO_TEST_CASE(TestRocalutionSolver)
{
    pt::ptree prm;

    // Read parameters.
    {
        std::ifstream file("options_flexiblesolver.json");
        pt::read_json(file, prm);
    }

    rocalution::init_rocalution();
    auto rocalution_backend_descriptor = rocalution::_get_backend_descriptor();

    if (rocalution_backend_descriptor->accelerator) {
        // test rocalution with 3x3 blocks
        test3(prm);
    } else {
        BOOST_ERROR("Problem with initializing a device.");
    }
}
