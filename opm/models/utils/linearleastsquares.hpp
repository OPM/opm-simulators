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
#ifndef LINEAR_LEAST_SQUARES_HPP
#define LINEAR_LEAST_SQUARES_HPP

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/SmallDenseMatrixUtils.hpp>

#include <stdexcept>

namespace Opm
{

/*!
 * @brief Linear least squares calculations and properties
 *
 * @warning Do not use with large matrices since matrix inversion is LU decomposition
 */
template <class Scalar>
class LinearLeastSquares
{
    using Matrix = Dune::DynamicMatrix<Scalar>;
    using Vector = Dune::DynamicVector<Scalar>;

public:
    /*!
     * @brief Constructor
     *
     * @param A Coefficient matrix
     * @param b Right-hand side (data) vector
     */
    LinearLeastSquares(const Matrix& A, const Vector& b)
        : A_(A), b_(b), x_(A.M())
    {
    }

    /*!
     * @brief Solve linear least squares system
     *
     * Done by solving the normal equations: (A^T*A)*x = A^T*b
     */
    void solve()
    {
        solveNormalEquations_();
    }

    /*!
     * @brief Read-only vector of calculated coefficient vector
     *
     * @return Coeff. vector
     */
    const Vector& x() const
    {
        return x_;
    }

    /*!
     * @brief Evaluate regression model at input x vector
     *
     * @param x Input vector
     * @return Regression model value
     */
    Scalar evaluate(const Vector& x) const
    {
        return x_ * x;
    }

    /*!
     * @brief Measure of the discrepancy between data and regression model
     *
     * @return Sum of the square of residual
     */
    Scalar residualSumOfSquares() const
    {
        Vector r(b_);
        A_.mmv(x_, r);
        return r.two_norm2();
    }

    /*!
     * @brief Value for how well regression model represents the model
     *
     * @return Model sum of squares
     */
    Scalar explainedSumOfSquares() const
    {
        Scalar ymean = 0.0;
        const auto ndata = b_.N();
        for (size_t i = 0; i < ndata; ++i) {
            ymean += b_[i];
        }
        ymean /= ndata;

        Vector r(ndata, ymean);
        A_.mmv(x_, r);

        return r.two_norm2();
    }

    /*!
     * @brief Sum of all squared differences
     *
     * Total sum of squares = explained sum of squares + residual sum of squares
     *
     * @return Total sum of squares
     */
    Scalar totalSumOfSquares() const
    {
        Scalar ymean = 0.0;
        const auto ndata = b_.N();
        for (size_t i = 0; i < ndata; ++i) {
            ymean += b_[i];
        }
        ymean /= ndata;

        Vector r(ndata, ymean);
        r -= b_;

        return r.two_norm2();
    }

    /*!
     * @brief Coefficient of determination
     *
     * Typical value for how well a regression model fits the data
     *
     * @return R^2 value
     */
    Scalar RSquared() const
    {
        const auto tss = totalSumOfSquares();
        if (tss < 1e-16) {
            OPM_THROW(std::runtime_error,
                      "Total sum of squares is close to zero, hence R^2 is undefined!");
        }
        return explainedSumOfSquares() / tss;
    }

private:
    /*!
     * @brief Solve the linear least squares system
     */
    void solveNormalEquations_()
    {
        // Right-hand side, bhat = A^T*b
        Vector bhat(A_.M());
        A_.mtv(b_, bhat);

        // Normal matrix Ahat = A^T*A
        Matrix Ahat(A_.M(), A_.M());
        detail::multMatrixTransposed(A_, A_, Ahat);

        // Solve normal equations x = Ahat^{-1}*bhat
        Ahat.solve(x_, bhat, /*doPivoting=*/true);
    }

    Matrix A_{}; ///< Input matrix
    Vector b_{}; ///< Data vector
    Vector x_{}; ///< Coefficient vector
};

} // Opm

#endif // LINEAR_LEAST_SQUARES_HPP
