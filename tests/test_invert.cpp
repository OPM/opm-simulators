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
                      Opm::NumericalProblem);
}
