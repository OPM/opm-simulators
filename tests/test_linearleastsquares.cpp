// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#include "config.h"

#define BOOST_TEST_MODULE TestLinearLeastSquares

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <opm/models/utils/linearleastsquares.hpp>

using namespace Opm;

using ScalarTypes = boost::mpl::list<double, float>;

BOOST_AUTO_TEST_CASE_TEMPLATE(SimpleOverdeterminedSystem, Scalar, ScalarTypes)
{
    // Solve an overdetermined system Ax=b using linear least squares (i.e. normal equations)
    // Example from: https://en.wikipedia.org/wiki/Linear_least_squares
    Dune::DynamicMatrix<Scalar> A(3, 2);
    A[0][0] = 1.0;
    A[0][1] = 0.0;
    A[1][0] = 0.0;
    A[1][1] = 1.0;
    A[2][0] = 1.0;
    A[2][1] = 1.0;

    Dune::DynamicVector<Scalar> b(3);
    b[0] = 1.0;
    b[1] = 1.0;
    b[2] = 0.0;

    // Instantiate LinearLeastSquares and solve normal equations
    LinearLeastSquares<Scalar> lsq(A, b);
    lsq.solve();

    // Check results
    Dune::DynamicVector<Scalar> xExpected(2);
    xExpected[0] = 1.0 / 3.0;
    xExpected[1] = 1.0 / 3.0;

    const auto& x = lsq.x();
    BOOST_CHECK_EQUAL(x.N(), xExpected.N());
    for (std::size_t i = 0; i < x.N(); ++i) {
        BOOST_CHECK_CLOSE(x[i], xExpected[i], 1e-4);
    }

    // Check sum of squared residuals
    BOOST_CHECK_CLOSE(lsq.residualSumOfSquares(), 4.0 / 3.0, 1e-4);

    // Check total sum of squares
    BOOST_CHECK_CLOSE(lsq.totalSumOfSquares(), 2.0 / 3.0, 1e-4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(LinearRegression, Scalar, ScalarTypes)
{
    // Linear regression for line: y = beta_0 + beta_1*x, with some y-data
    // Example from: https://en.wikipedia.org/wiki/Linear_least_squares#Example
    Dune::DynamicMatrix<Scalar> X(4, 2);
    X[0][0] = 1.0;
    X[0][1] = 1.0;
    X[1][0] = 1.0;
    X[1][1] = 2.0;
    X[2][0] = 1.0;
    X[2][1] = 3.0;
    X[3][0] = 1.0;
    X[3][1] = 4.0;

    Dune::DynamicVector<Scalar> y(4);
    y[0] = 6.0;
    y[1] = 5.0;
    y[2] = 7.0;
    y[3] = 10.0;

    // Solve linear least squares
    LinearLeastSquares<Scalar> lsq(X, y);
    lsq.solve();

    // Check against solution
    Dune::DynamicVector<Scalar> betaExpected(2);
    betaExpected[0] = 3.5;
    betaExpected[1] = 1.4;

    const auto& beta = lsq.x();
    BOOST_CHECK_EQUAL(beta.N(), betaExpected.N());
    for (std::size_t i = 0; i < beta.N(); ++i) {
        BOOST_CHECK_CLOSE(beta[i], betaExpected[i], 1e-4);
    }

    // Check sum of squared residuals
    BOOST_CHECK_CLOSE(lsq.residualSumOfSquares(), 4.2, 1e-4);

    // Check total sum of squares
    BOOST_CHECK_CLOSE(lsq.totalSumOfSquares(), 14.0, 1e-4);

    // Check R-squared, goodness-of-fit
    BOOST_CHECK_CLOSE(lsq.RSquared(), 0.7, 1e-4);

    // Check interpolation, i.e. y = beta^T * x for an arbitrary x
    Dune::DynamicVector<Scalar> xPred(2);
    xPred[0] = 1.0;
    xPred[1] = 5.0;
    BOOST_CHECK_CLOSE(lsq.evaluate(xPred), 10.5, 1e-4);
}
