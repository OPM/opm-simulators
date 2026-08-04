/*
  Copyright 2017 IRIS AS

  This file is part of the Open Porous Media Project (OPM).

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

#define BOOST_TEST_MODULE InvertSpecializationTest
#include <boost/test/unit_test.hpp>
#include <opm/simulators/linalg/matrixblock.hh>


void checkIdentity(Dune::FieldMatrix<double, 4, 4> M) {
    double diag = 0.0;
    double offDiag = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j)
                diag += M[i][j];
            else
                offDiag += M[i][j];
        }
    }
    BOOST_CHECK_CLOSE(4, diag, 1e-14);
    BOOST_CHECK_SMALL(offDiag, 1e-14);
}

BOOST_AUTO_TEST_CASE(Invert4x4)
{
    using BaseType = Dune::FieldMatrix<double, 4, 4>;
    BaseType matrix;
    BaseType inverse;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
                matrix[i][j] = i + 4*j + 1;
        }
    }
    BaseType matrix_sing (matrix);
    // make matrix non-singular
    matrix[3][0] = 5;
    matrix[0][3] = 14;

    double det = Opm::detail::invertMatrix4<Opm::detail::FMat4>(matrix, inverse);
    BOOST_CHECK_CLOSE(4, det, 1e-14);

    // check matrix * inverse close to identiy
    checkIdentity(matrix.rightmultiply(inverse));

    // check singular matrix
    BOOST_CHECK_THROW(Opm::detail::invertMatrix4<Opm::detail::FMat4>(matrix_sing, inverse),
                      Dune::MatrixBlockError);
}

BOOST_AUTO_TEST_CASE(Invert4x4WithUnderflowingDeterminant)
{
    using BaseType = Dune::FieldMatrix<double, 4, 4>;

    // A well conditioned matrix (condition number about 9.5) scaled such that its
    // determinant, which is homogeneous of degree one in every row, falls below the
    // 1e-40 threshold used to decide whether the cofactor determinant can be divided
    // by. Scaling leaves the matrix just as invertible as it was, so this must be
    // inverted through the pivoted fallback rather than reported as singular.
    BaseType matrix;
    const double base[4][4] = {{2, 1, 0, 0},
                               {1, 2, 1, 0},
                               {0, 1, 2, 1},
                               {0, 0, 1, 2}};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            matrix[i][j] = base[i][j] * 1e-11;
        }
    }

    BaseType inverse;
    double det = 0.0;
    BOOST_CHECK_NO_THROW(det = Opm::detail::invertMatrix4<Opm::detail::FMat4>(matrix, inverse));
    BOOST_CHECK_LT(std::abs(det), 1e-40);

    // the inverse is still accurate even though the determinant underflowed
    const BaseType product = matrix.rightmultiply(inverse);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            BOOST_CHECK_SMALL(product[i][j] - (i == j ? 1.0 : 0.0), 1e-12);
        }
    }
}
