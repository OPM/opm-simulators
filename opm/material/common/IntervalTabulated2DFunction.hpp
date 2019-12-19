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
 * \copydoc Opm::IntervalTabulated2DFunction
 */
#ifndef OPM_INTERVAL_TABULATED_2D_FUNCTION_HPP
#define OPM_INTERVAL_TABULATED_2D_FUNCTION_HPP

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <vector>
#include <limits>
#include <sstream>
#include <cassert>
#include <algorithm>

namespace Opm {

/*!
 * \brief Implements a function that depends on two variables.
 *
 * The function is sampled in regular intervals in both directions, i.e., the
 * interpolation cells are rectangles. The table can be extrapolated in either direction.
 */
template <class Scalar>
class IntervalTabulated2DFunction
{
public:
    IntervalTabulated2DFunction()
    { }

    template <class DataContainer>
    IntervalTabulated2DFunction(const std::vector<Scalar>& xPos,
                                const std::vector<Scalar>& yPos,
                                const DataContainer& data,
                                const bool xExtrapolate = false,
                                const bool yExtrapolate = false)
        : xPos_(xPos)
        , yPos_(yPos)
        , samples_(data)
        , xExtrapolate_(xExtrapolate)
        , yExtrapolate_(yExtrapolate)
    {
#ifndef NDEBUG
        // in debug mode, ensure that the x and y positions arrays are strictly
        // mononically increasing.
        for (unsigned i = 0; i < xPos.size() - 1; ++ i) {
            if (xPos[i + 1] <= xPos[i])
                throw std::runtime_error("The array for the x-positions is not strictly increasing!");
        }

        for (unsigned i = 0; i < yPos.size() - 1; ++ i) {
            if (yPos[i + 1] <= yPos[i])
                throw std::runtime_error("The array for the y-positions is not strictly increasing!");
        }
#endif

        // make sure the size is correct
        if (numX() != samples_.size())
            throw std::runtime_error("numX() is not equal to the number of rows of the sampling points");

        for (unsigned xIdx = 0; xIdx < numX(); ++xIdx) {
            if (samples_[xIdx].size() != numY()) {
                std::ostringstream oss;
                oss << "The " << xIdx << "-th row of the sampling points has different size than numY() ";
                throw std::runtime_error(oss.str());
            }
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

    const std::vector<Scalar>& xPos() const
    { return xPos_; }

    const std::vector<Scalar>& yPos() const
    { return yPos_; }

    const std::vector<std::vector<Scalar>>& samples() const
    { return samples_; }

    bool xExtrapolate() const
    { return xExtrapolate_; }

    bool yExtrapolate() const
    { return yExtrapolate_; }

    bool operator==(const IntervalTabulated2DFunction<Scalar>& data) const {
        return this->xPos() == data.xPos() &&
               this->yPos() == data.yPos() &&
               this->samples() == data.samples() &&
               this->xExtrapolate() == data.xExtrapolate() &&
               this->yExtrapolate() == data.yExtrapolate();
    }

    /*!
     * \brief Returns the value of a sampling point.
     */
    Scalar valueAt(size_t i, size_t j) const
    { return samples_[i][j]; }

    /*!
     * \brief Returns true if a coordinate lies in the tabulated range
     */
    template <class Evaluation>
    bool applies(const Evaluation& x, const Evaluation& y) const
    { return appliesX(x) && appliesY(y); }

    /*!
     * \brief Returns true if a coordinate lies in the tabulated range on the x direction
     */
    template <class Evaluation>
    bool appliesX(const Evaluation& x) const
    { return xMin() <= x && x <= xMax(); }

    /*!
     * \brief Returns true if a coordinate lies in the tabulated range on the y direction
     */
    template <class Evaluation>
    bool appliesY(const Evaluation& y) const
    { return yMin() <= y && y <= yMax(); }


    /*!
     * \brief Evaluate the function at a given (x,y) position.
     *
     * If this method is called for a value outside of the tabulated
     * range, and extrapolation is not allowed in the corresponding direction,
     * a \c Opm::NumericalIssue exception is thrown.
     */
    template <typename Evaluation>
    Evaluation eval(const Evaluation& x, const Evaluation& y) const
    {
        if ((!xExtrapolate_ && !appliesX(x)) || (!yExtrapolate_ && !appliesY(y))) {
            std::ostringstream oss;
            oss << "Attempt to get undefined table value (" << x << ", " << y << ")";
            throw NumericalIssue(oss.str());
        };

        // bi-linear interpolation: first, calculate the x and y indices in the lookup
        // table ...
        const unsigned i = xSegmentIndex_(x);
        const unsigned j = ySegmentIndex_(y);

        // bi-linear interpolation / extrapolation
        const Evaluation alpha = xToAlpha(x, i);
        const Evaluation beta = yToBeta(y, j);

        const Evaluation s1 = valueAt(i, j) * (1.0 - beta) + valueAt(i, j + 1) * beta;
        const Evaluation s2 = valueAt(i + 1, j) * (1.0 - beta) + valueAt(i + 1, j + 1) * beta;

        Valgrind::CheckDefined(s1);
        Valgrind::CheckDefined(s2);

        // ... and combine them using the x position
        return s1*(1.0 - alpha) + s2*alpha;
    }

private:
    // the sampling points in the x-drection
    std::vector<Scalar> xPos_;
    // the sampling points in the y-drection
    std::vector<Scalar> yPos_;
    // data at the sampling points
    std::vector<std::vector<Scalar> > samples_;

    bool xExtrapolate_ = false;
    bool yExtrapolate_ = false;

    /*!
     * \brief Return the interval index of a given position on the x-axis.
     */
    template <class Evaluation>
    unsigned xSegmentIndex_(const Evaluation& x) const
    {
        assert(xExtrapolate_ || appliesX(x) );

        return segmentIndex_(x, xPos_);
    }

    /*!
     * \brief Return the interval index of a given position on the y-axis.
     */
    template <class Evaluation>
    unsigned ySegmentIndex_(const Evaluation& y) const
    {
        assert(yExtrapolate_ || appliesY(y) );

        return segmentIndex_(y, yPos_);
    }


    template <class Evaluation>
    static unsigned segmentIndex_(const Evaluation& v, const std::vector<Scalar>& vPos)
    {
        const unsigned n = vPos.size();
        assert(n >= 2);

        if (v <= vPos.front() || n == 2)
            return 0;
        else if (v >= vPos.back())
            return n - 2;

        assert(n > 2 && v > vPos.front() && v < vPos.back());

        // bisection. this assumes that the vPos array is strictly mononically
        // increasing.
        size_t lowerIdx = 0;
        size_t upperIdx = vPos.size() - 1;
        while (lowerIdx + 1 < upperIdx) {
            size_t pivotIdx = (lowerIdx + upperIdx) / 2;
            if (v < vPos[pivotIdx])
                upperIdx = pivotIdx;
            else
                lowerIdx = pivotIdx;
        }

        assert(vPos[lowerIdx] <= v);
        assert(v <= vPos[lowerIdx + 1]);
        return lowerIdx;
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
        Scalar y1 = yPos_[ySegmentIdx];
        Scalar y2 = yPos_[ySegmentIdx + 1];
        return (y - y1)/(y2 - y1);
    }

};
} // namespace Opm

#endif
