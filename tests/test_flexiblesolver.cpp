/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE OPM_test_FlexibleSolver
#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>

#include <fstream>
#include <iostream>


template <int bz>
Dune::BlockVector<Dune::FieldVector<double, bz>>
testSolver(const Opm::PropertyTree& prm, const std::string& matrix_filename, const std::string& rhs_filename)
{
    using Matrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;
    Matrix matrix;
    {
        std::ifstream mfile(matrix_filename);
        if (!mfile) {
            throw std::runtime_error("Could not read matrix file");
        }
        using M = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
        readMatrixMarket(reinterpret_cast<M&>(matrix), mfile); // Hack to avoid hassle
    }
    Vector rhs;
    {
        std::ifstream rhsfile(rhs_filename);
        if (!rhsfile) {
            throw std::runtime_error("Could not read rhs file");
        }
        readMatrixMarket(rhs, rhsfile);
    }
    bool transpose = false;

    if(prm.get<std::string>("preconditioner.type") == "cprt"){
        transpose = true;
    }
    auto wc = [&matrix, transpose]()
    {
        return Opm::Amg::getQuasiImpesWeights<Matrix, Vector>(matrix, 1, transpose);
    };

    using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
    SeqOperatorType op(matrix);
    Dune::FlexibleSolver<SeqOperatorType> solver(op, prm, wc, 1);
    Vector x(rhs.size());
    Dune::InverseOperatorResult res;
    solver.apply(x, rhs, res);
    return x;
}

BOOST_AUTO_TEST_CASE(TestFlexibleSolver)
{
    // Read parameters.
    Opm::PropertyTree prm("options_flexiblesolver.json");

    // Test with 1x1 block solvers.
    {
        const int bz = 1;
        auto sol = testSolver<bz>(prm, "matr33.txt", "rhs3.txt");
        Dune::BlockVector<Dune::FieldVector<double, bz>> expected {-1.62493,
                                                                   -1.76435e-06,
                                                                   1.86991e-10,
                                                                   -458.542,
                                                                   2.28308e-06,
                                                                   -2.45341e-07,
                                                                   -1.48005,
                                                                   -5.02264e-07,
                                                                   -1.049e-05};
        BOOST_REQUIRE_EQUAL(sol.size(), expected.size());
        for (size_t i = 0; i < sol.size(); ++i) {
            for (int row = 0; row < bz; ++row) {
                BOOST_CHECK_CLOSE(sol[i][row], expected[i][row], 1e-3);
            }
        }
    }

    // Test with 3x3 block solvers.
    {
        const int bz = 3;
        auto sol = testSolver<bz>(prm, "matr33.txt", "rhs3.txt");
        Dune::BlockVector<Dune::FieldVector<double, bz>> expected {{-1.62493, -1.76435e-06, 1.86991e-10},
                                                                   {-458.542, 2.28308e-06, -2.45341e-07},
                                                                   {-1.48005, -5.02264e-07, -1.049e-05}};
        BOOST_REQUIRE_EQUAL(sol.size(), expected.size());
        for (size_t i = 0; i < sol.size(); ++i) {
            for (int row = 0; row < bz; ++row) {
                BOOST_CHECK_CLOSE(sol[i][row], expected[i][row], 1e-3);
            }
        }
    }
}
