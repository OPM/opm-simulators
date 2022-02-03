/*
  Copyright 2021 Equinor ASA

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

#define BOOST_TEST_MODULE SolveTransposed3x3
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>

#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <dune/istl/bcrsmatrix.hh>

#include <opm/simulators/linalg/bda/opencl/CPR.hpp>

BOOST_AUTO_TEST_CASE(testsolvetransposed3x3)
{
    const unsigned numTests = 3;
    const unsigned blockSize = 3;

    std::vector<std::vector<double> > A = {{4, 2, 1,
                                            3, 4, 2,
                                            2, 4, 3},
                                           {1, 2, 4,
                                            1, 3, 2,
                                            2, 4, 2},
                                           {1, 2, 2,
                                            1, 3, 5,
                                            3, 2, 4}};

    std::vector<std::vector<double> > b = {{0, 1, 0},
                                           {1, 3, 5},
                                           {2, 4, 5}};

    std::vector<std::vector<double> > x_expected = {{-0.5,   1, -0.5},
                                                    {   1,   1, -0.5},
                                                    { 1.3, 0.4,  0.1}};

    for (unsigned testId = 0; testId < numTests; ++testId) {
        std::vector<double> x(blockSize);
        Opm::Accelerator::solve_transposed_3x3(A[testId].data(), b[testId].data(), x.data());

        for (unsigned i = 0; i < blockSize; ++i) {
            BOOST_CHECK_CLOSE(x[i], x_expected[testId][i], 1e-6);
        }
    }

    // retest cases using Dune methods
    using Mat3 = Dune::FieldMatrix<double, 3, 3>;
    using Vec3 = Dune::FieldVector<double, 3>;
    for (unsigned testId = 0; testId < numTests; ++testId) {
        Mat3 a3 = {{ A[testId][0], A[testId][1], A[testId][2] },
                   { A[testId][3], A[testId][4], A[testId][5] },
                   { A[testId][6], A[testId][7], A[testId][8] } };
        Vec3 y({b[testId][0], b[testId][1], b[testId][2]});
        Mat3 b3 = a3;

        // b3 = inv(a3^T)
        Dune::FMatrixHelp::invertMatrix_retTransposed(a3, b3);

        // x = b3 * y
        Vec3 x = Dune::FMatrixHelp::mult(b3, y);

        for (unsigned i = 0; i < blockSize; ++i) {
            BOOST_CHECK_CLOSE(x[i], x_expected[testId][i], 1e-6);
        }
    }
}
