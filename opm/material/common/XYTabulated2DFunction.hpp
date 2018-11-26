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
/*!
 * \file
 *
 * \copydoc Opm::XYTabulated2DFunction.hpp
 */
#ifndef OPM_XY_TABULATED_2D_FUNCTION_HPP
#define OPM_XY_TABULATED_2D_FUNCTION_HPP

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <iostream>
#include <vector>
#include <limits>
#include <tuple>
#include <sstream>
#include <cassert>

namespace Opm {

/*!
 * \brief Implements a scalar function that depends on two variables and which is sampled
 *        in X and Y direction.
 */

// TODO: for now, we should not allow extraploation on the lower end, because we do not deal with
// TODO: negative values
// TODO: remove all the functions that will not be called from outside, or put them as private members

template <class Scalar>
class XYTabulated2DFunction
{

public:
    XYTabulated2DFunction()
    { }

    XYTabulated2DFunction(const std::vector<Scalar>& x_pos,
                          const std::vector<Scalar>& y_pos,
                          const std::vector<std::vector<Scalar> >& data,
                          const bool x_extrapolate,
                          const bool y_extrapolate)
       : xPos_(x_pos)
       , yPos_(y_pos)
       , samples_(data)
       , xExtrapolate_(x_extrapolate)
       , yExtrapolate_(y_extrapolate)
    {
        // make sure the size is correct
        if (numX() != samples_.size()) {
            throw std::runtime_error("numX() is not equal with the number of rows of the sample points");
        }

        for (unsigned ix = 0; ix < numX(); ++ix) {
            if (samples_[ix].size() != numY()) {
                std::ostringstream oss;
                oss << "the " << ix << "th row of the smaples_ has different number of data from numY() ";
                throw std::runtime_error(oss.str());
            }
        }
    }

    /*!
     * \brief Evaluate the function at a given (x,y) position.
     *
     * If this method is called for a value outside of the tabulated
     * range, a \c Opm::NumericalIssue exception is thrown.
     */
    template <typename DataType>
    void eval(const DataType& x, const DataType& y, DataType& result) const
    {
        if ( (!xExtrapolate_ && !appliesX(x)) ||
             (!yExtrapolate_ && !appliesY(y)) ) {
            std::ostringstream oss;
            oss << "Attempt to get undefined table value (" << x << ", " << y << ")";
            throw NumericalIssue(oss.str());
        };

        // bi-linear interpolation: first, calculate the x and y indices in the lookup
        // table ...
        const unsigned i = xSegmentIndex(x);

        // std::cout << " i is " << i << std::endl;

        const unsigned j = ySegmentIndex(y);
        // std::cout << " j is " << j << std::endl;

        if ((i == numX() - 1) && (j == numY() - 1)) {
            result = valueAt(i, j);
            return;
        }

        if (i == numX() - 1) {
            const DataType beta = yToBeta(y, j);
            result = valueAt(i, j) * (1.0 - beta) + valueAt(i, j + 1) * beta;
            Valgrind::CheckDefined(result);
            return;
        }

        if (j == numY() - 1) {
            const DataType alpha = xToAlpha(x, i);
            result = valueAt(i, j) * (1.0 - alpha) + valueAt(i + 1, j) * alpha;
            Valgrind::CheckDefined(result);
            return;
        }

        // normal bi-linear interpolation / extrapolation
        const DataType alpha = xToAlpha(x, i);
        const DataType beta = yToBeta(y, j);

        const DataType s1 = valueAt(i, j) * (1.0 - beta) + valueAt(i, j + 1) * beta;
        const DataType s2 = valueAt(i + 1, j) * (1.0 - beta) + valueAt(i + 1, j + 1) * beta;

        Valgrind::CheckDefined(s1);
        Valgrind::CheckDefined(s2);

        // ... and combine them using the x position
        result = s1 * (1.0 - alpha) + s2 * alpha;
        Valgrind::CheckDefined(result);
    }

private:
    // the sampling points in the x-drection
    std::vector<Scalar> xPos_;
    // the sampling points in the y-drection
    std::vector<Scalar> yPos_;

    std::vector<std::vector<Scalar> > samples_;

    bool xExtrapolate_ = false;
    bool yExtrapolate_ = false;

    /*!
     * \brief Returns the minimum of the X coordinate of the sampling points.
     */
    Scalar xMin() const
    { return xPos_.front(); }

    /*!
     * \brief Returns the maximum of the X coordinate of the sampling points.
     */
    Scalar xMax() const
    { return xPos_.back(); }

    /*!
     * \brief Returns the minimum of the Y coordinate of the sampling points.
     */
    Scalar yMin() const
    { return yPos_.front(); }

    /*!
     * \brief Returns the maximum of the Y coordinate of the sampling points.
     */
    Scalar yMax() const
    { return yPos_.back(); }

    /*!
     * \brief Returns the value of a sampling point.
     */
    Scalar valueAt(size_t i, size_t j) const
    { return samples_[i][j]; }

    /*!
     * \brief Returns true if a coordinate lies in the tabulated range on the x direction
     */
    template <class Evaluation>
    bool appliesX(const Evaluation& x) const
    {
        if (x < xMin() || xMax() < x) {
            return false;
        } else {
            return true;
        }
    }

    /*!
     * \brief Returns true if a coordinate lies in the tabulated range on the y direction
     */
    template <class Evaluation>
    bool appliesY(const Evaluation& y) const
    {
        if (y < yMin() || yMax() < y) {
            return false;
        } else {
            return true;
        }
    }

    /*!
     * \brief Returns the number of sampling points in X direction.
     */
    size_t numX() const
    { return xPos_.size(); }

    /*!
     * \brief Returns the number of sampling points in Y direction.
     */
    size_t numY() const
    { return yPos_.size(); }

    /*!
     * \brief Return the interval index of a given position on the x-axis.
     */
    template <class Evaluation>
    unsigned xSegmentIndex(const Evaluation& x) const
    {
        assert(xExtrapolate_ || (xMin() <= x && x <= xMax()));

        // TODO: we should have better way to enforce no extrapolation in the low end
        // assert(x >= xMin());

        // we need at least two sampling points!
        // TODO: not necessary, right?
        assert(numX() >= 2);

        if (x <= xMin())
            return 0;
        else if (x >= xMax())
            if (xExtrapolate_)
                return numX() - 1;
            else
                return numX() - 2;
        else {
            assert(numX() >= 3);

            // bisection
            unsigned lowerIdx = 0;
            unsigned upperIdx = numX() - 1;
            while (lowerIdx + 1 < upperIdx) {
                unsigned pivotIdx = (lowerIdx + upperIdx) / 2;
                if (x < xPos_[pivotIdx])
                    upperIdx = pivotIdx;
                else
                    lowerIdx = pivotIdx;
            }

            return lowerIdx;
        }
    }

    /*!
     * \brief Return the interval index of a given position on the y-axis.
     */
    template <class Evaluation>
    unsigned ySegmentIndex(const Evaluation& y) const
    {
        assert(yExtrapolate_ || (yMin() <= y && y <= yMax()));

        // we need at least two sampling points!
        assert(numY() >= 2);

        if (y <= yMin())
            return 0;
        else if (y >= yMax())
            if (yExtrapolate_)
                return numY() - 2;
            else
                return numY() - 1;
        else {
            assert(numY() >= 3);

            // bisection
            unsigned lowerIdx = 0;
            unsigned upperIdx = numY() - 1;
            while (lowerIdx + 1 < upperIdx) {
                const unsigned pivotIdx = (lowerIdx + upperIdx) / 2;
                if (y < yPos_[pivotIdx])
                    upperIdx = pivotIdx;
                else
                    lowerIdx = pivotIdx;
            }

            return lowerIdx;
        }
    }

    /*!
     * \brief Return the relative position of an x value in an intervall
     *
     * The returned value can be larger than 1 or smaller than zero if it is outside of
     * the range of the segment. In particular this happens for the extrapolation case.
     */
    template <class Evaluation>
    Evaluation xToAlpha(const Evaluation& x, unsigned xSegmentIdx) const
    {
        if ( xSegmentIdx == numX() -1 )
            return 0.;

        Scalar x1 = xPos_[xSegmentIdx];
        Scalar x2 = xPos_[xSegmentIdx + 1];
        return (x - x1)/(x2 - x1);
    }

    /*!
     * \brief Return the relative position of an y value in an interval
     *
     * The returned value can be larger than 1 or smaller than zero if it is outside of
     * the range of the segment. In particular this happens for the extrapolation case.
     */
    template <class Evaluation>
    Evaluation yToBeta(const Evaluation& y, unsigned ySegmentIdx) const
    {
        if ( ySegmentIdx == numY() - 1)
            return 0.;

        Scalar y1 = yPos_[ySegmentIdx];
        Scalar y2 = yPos_[ySegmentIdx + 1];
        return (y - y1)/(y2 - y1);
    }

};
} // namespace Opm

#endif
