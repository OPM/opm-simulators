/*
  Copyright Equinor ASA 2026

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
#define BOOST_TEST_MODULE OPM_test_MswCprWellDiagonal
#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/MSWellHelpers.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <cstddef>

namespace {

using Scalar = double;
constexpr int numWellEq = 4; // three conservation equations plus segment pressure
constexpr int pressureColumn = numWellEq - 1;

using Block = Dune::FieldMatrix<Scalar, numWellEq, numWellEq>;
using DiagMatWell = Dune::BCRSMatrix<Block>;
using WellWeight = Dune::FieldVector<Scalar, 3>;

// Two segments, each with an entry on the other: the tridiagonal pattern an
// inlet/outlet pair produces. Only the pressure column is filled, since that is
// all the contraction reads.
DiagMatWell makeD(const Scalar d00, const Scalar d01,
                  const Scalar d10, const Scalar d11)
{
    DiagMatWell D(2, 2, 4, DiagMatWell::row_wise);
    for (auto row = D.createbegin(); row != D.createend(); ++row) {
        row.insert(0);
        row.insert(1);
    }
    D = 0.0;
    // Put the same value on every conservation row, so the contraction against
    // a weight vector is just (sum of weights) times that value.
    for (int i = 0; i < 3; ++i) {
        D[0][0][i][pressureColumn] = d00;
        D[0][1][i][pressureColumn] = d01;
        D[1][0][i][pressureColumn] = d10;
        D[1][1][i][pressureColumn] = d11;
    }
    return D;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(ContractionSumsEveryBlockOfD)
{
    const auto D = makeD(2.0, -0.5, -0.25, 3.0);
    const WellWeight lambda {0.5, 0.25, 0.125};
    const Scalar weightSum = 0.5 + 0.25 + 0.125;

    const auto diag = Opm::mswellhelpers::contractCprWellDiagonal(D, lambda, pressureColumn);

    // (2.0 - 0.5 - 0.25 + 3.0) * sum(lambda)
    BOOST_CHECK_CLOSE(diag, 4.25 * weightSum, 1e-12);
}

BOOST_AUTO_TEST_CASE(OffDiagonalSegmentCouplingIsKept)
{
    // Same diagonal, different segment-to-segment coupling. The row-sum
    // convention never reads D at all and cannot tell these two apart; the
    // contraction must.
    const WellWeight lambda {1.0, 1.0, 1.0};
    const auto weak = Opm::mswellhelpers::contractCprWellDiagonal(makeD(2.0, 0.0, 0.0, 3.0),
                                                                  lambda, pressureColumn);
    const auto strong = Opm::mswellhelpers::contractCprWellDiagonal(makeD(2.0, -0.5, -0.25, 3.0),
                                                                    lambda, pressureColumn);

    BOOST_CHECK_CLOSE(weak, 5.0 * 3.0, 1e-12);
    BOOST_CHECK_CLOSE(strong, 4.25 * 3.0, 1e-12);
    BOOST_CHECK_LT(strong, weak);
}

BOOST_AUTO_TEST_CASE(ExactCancellationIsRegularised)
{
    // A zero here would make the coarse pressure system singular.
    const auto D = makeD(2.0, -1.0, -1.0, 0.0);
    const WellWeight lambda {1.0, 1.0, 1.0};

    const auto diag = Opm::mswellhelpers::contractCprWellDiagonal(D, lambda, pressureColumn);

    BOOST_CHECK_EQUAL(diag, 1.0);
}
