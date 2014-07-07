/*
  Copyright (C) 2013 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \copydoc Opm::UniformTabulated2DFunction
 */
#ifndef OPM_UNIFORM_TABULATED_2D_FUNCTION_HPP
#define OPM_UNIFORM_TABULATED_2D_FUNCTION_HPP

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <vector>

#include <assert.h>

namespace Opm {

/*!
 * \brief Implements a scalar function that depends on two variables and which is sampled
 *        on an uniform X-Y grid.
 *
 * This class can be used when the sampling points are calculated at
 * run time.
 */
template <class Scalar>
class UniformTabulated2DFunction
{
public:
    UniformTabulated2DFunction()
    { }

     /*!
      * \brief Constructor where the tabulation parameters are already
      *        provided.
      */
    UniformTabulated2DFunction(Scalar xMin, Scalar xMax, int m,
                               Scalar yMin, Scalar yMax, int n)
    {
        resize(xMin, xMax, m, yMin, yMax, n);
    }

    /*!
     * \brief Resize the tabulation to a new range.
     */
    void resize(Scalar xMin, Scalar xMax, int m,
                Scalar yMin, Scalar yMax, int n)
    {
        samples_.resize(m*n);

        m_ = m;
        n_ = n;

        xMin_ = xMin;
        xMax_ = xMax;

        yMin_ = yMin;
        yMax_ = yMax;
    }

    /*!
     * \brief Returns the minimum of the X coordinate of the sampling points.
     */
    Scalar xMin() const
    { return xMin_; }

    /*!
     * \brief Returns the maximum of the X coordinate of the sampling points.
     */
    Scalar xMax() const
    { return xMax_; }

    /*!
     * \brief Returns the minimum of the Y coordinate of the sampling points.
     */
    Scalar yMin() const
    { return yMin_; }

    /*!
     * \brief Returns the maximum of the Y coordinate of the sampling points.
     */
    Scalar yMax() const
    { return yMax_; }

    /*!
     * \brief Returns the number of sampling points in X direction.
     */
    int numX() const
    { return m_; }

    /*!
     * \brief Returns the number of sampling points in Y direction.
     */
    int numY() const
    { return n_; }

    /*!
     * \brief Return the position on the x-axis of the i-th interval.
     */
    Scalar iToX(int i) const
    {
        assert(0 <= i && i < numX());

        return xMin() + i*(xMax() - xMin())/(numX() - 1);
    }

    /*!
     * \brief Return the position on the y-axis of the j-th interval.
      */
    Scalar jToY(int j) const
    {
        assert(0 <= j && j < numY());

        return yMin() + j*(yMax() - yMin())/(numY() - 1);
    }

    /*!
     * \brief Return the interval index of a given position on the x-axis.
     *
     * This method returns a *floating point* number. The integer part
     * should be interpreted as interval, the decimal places are the
     * position of the x value between the i-th and the (i+1)-th
     * sample point.
      */
    Scalar xToI(Scalar x) const
    { return (x - xMin())/(xMax() - xMin())*(numX() - 1); }

    /*!
     * \brief Return the interval index of a given position on the y-axis.
     *
     * This method returns a *floating point* number. The integer part
     * should be interpreted as interval, the decimal places are the
     * position of the y value between the j-th and the (j+1)-th
     * sample point.
     */
    Scalar yToJ(Scalar y) const
    { return (y - yMin())/(yMax() - yMin())*(numY() - 1); }

    /*!
     * \brief Returns true iff a coordinate lies in the tabulated range
     */
    bool applies(Scalar x, Scalar y) const
    { return xMin() <= x && x <= xMax() && yMin() <= y && y <= yMax(); }

    /*!
     * \brief Evaluate the function at a given (x,y) position.
     *
     * If this method is called for a value outside of the tabulated
     * range, a \c Opm::NumericalProblem exception is thrown.
     */
    Scalar eval(Scalar x, Scalar y) const
    {
#ifndef NDEBUG
        if (!applies(x,y))
        {
            OPM_THROW(NumericalProblem,
                       "Attempt to get tabulated value for ("
                       << x << ", " << y
                       << ") on a table of extend "
                       << xMin() << " to " << xMax() << " times "
                       << yMin() << " to " << yMax());
        };
#endif

        Scalar alpha = xToI(x);
        Scalar beta = yToJ(y);

        int i = std::max(0, std::min(numX() - 2, static_cast<int>(alpha)));
        int j = std::max(0, std::min(numY() - 2, static_cast<int>(beta)));

        alpha -= i;
        beta -= j;

        // bi-linear interpolation
        Scalar s1 = getSamplePoint(i, j)*(1.0 - alpha) + getSamplePoint(i + 1, j)*alpha;
        Scalar s2 = getSamplePoint(i, j + 1)*(1.0 - alpha) + getSamplePoint(i + 1, j + 1)*alpha;
        return s1*(1.0 - beta) + s2*beta;
    }

    /*!
     * \brief Get the value of the sample point which is at the
     *         intersection of the \f$i\f$-th interval of the x-Axis
     *         and the \f$j\f$-th of the y-Axis.
     */
    Scalar getSamplePoint(int i, int j) const
    {
        assert(0 <= i && i < m_);
        assert(0 <= j && j < n_);

        return samples_[j*m_ + i];
    }

    /*!
     * \brief Set the value of the sample point which is at the
     *        intersection of the \f$i\f$-th interval of the x-Axis
     *        and the \f$j\f$-th of the y-Axis.
     */
    void setSamplePoint(int i, int j, Scalar value)
    {
        assert(0 <= i && i < m_);
        assert(0 <= j && j < n_);

        samples_[j*m_ + i] = value;
    }

private:
    // the vector which contains the values of the sample points
    // f(x_i, y_j). don't use this directly, use getSamplePoint(i,j)
    // instead!
    std::vector<Scalar> samples_;

    // the number of sample points in x direction
    int m_;

    // the number of sample points in y direction
    int n_;

    // the range of the tabulation on the x axis
    Scalar xMin_;
    Scalar xMax_;

    // the range of the tabulation on the y axis
    Scalar yMin_;
    Scalar yMax_;
};
} // namespace Opm

#endif
