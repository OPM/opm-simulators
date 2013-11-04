// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \copydoc Opm::DynamicTabulated2DFunction
 */
#ifndef OPM_DYNAMIC_TABULATED_2D_FUNCTION_HH
#define OPM_DYNAMIC_TABULATED_2D_FUNCTION_HH

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/material/Tabulated2dFunction.hpp>

#include <vector>

#include <assert.h>

namespace Opm {

/*!
 * \copydoc Opm::Tabulated2DFunction
 *
 * This class can be used when the sampling points are calculated at
 * run time.
 */
template <class Scalar>
class DynamicTabulated2DFunction
    : public Tabulated2DFunction<Scalar, DynamicTabulated2DFunction<Scalar> >
{
public:
    DynamicTabulated2DFunction()
    { }

     /*!
      * \brief Constructor where the tabulation parameters are already
      *        provided.
      */
    DynamicTabulated2DFunction(Scalar xMin, Scalar xMax, int m,
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
