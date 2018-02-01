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
 * \copydoc Opm::UniformXTabulated2DFunction
 */
#ifndef OPM_UNIFORM_X_TABULATED_2D_FUNCTION_HPP
#define OPM_UNIFORM_X_TABULATED_2D_FUNCTION_HPP

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
 *        uniformly in the X direction, but non-uniformly on the Y axis-
 *
 * "Uniform on the X-axis" means that all Y sampling points must be located along a line
 * for this value. This class can be used when the sampling points are calculated at run
 * time.
 */
template <class Scalar>
class UniformXTabulated2DFunction
{
    typedef std::tuple</*x=*/Scalar, /*y=*/Scalar, /*value=*/Scalar> SamplePoint;

public:
    UniformXTabulated2DFunction()
    { }

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
     * \brief Returns the value of the X coordinate of the sampling points.
     */
    Scalar xAt(size_t i) const
    { return xPos_[i]; }

    /*!
     * \brief Returns the value of the Y coordinate of a sampling point.
     */
    Scalar yAt(size_t i, size_t j) const
    { return std::get<1>(samples_[i][j]); }

    /*!
     * \brief Returns the value of a sampling point.
     */
    Scalar valueAt(size_t i, size_t j) const
    { return std::get<2>(samples_[i][j]); }

    /*!
     * \brief Returns the number of sampling points in X direction.
     */
    size_t numX() const
    { return xPos_.size(); }

    /*!
     * \brief Returns the minimum of the Y coordinate of the sampling points for a given column.
     */
    Scalar yMin(unsigned i) const
    { return std::get<1>(samples_.at(i).front()); }

    /*!
     * \brief Returns the maximum of the Y coordinate of the sampling points for a given column.
     */
    Scalar yMax(unsigned i) const
    { return std::get<1>(samples_.at(i).back()); }

    /*!
     * \brief Returns the number of sampling points in Y direction a given column.
     */
    size_t numY(unsigned i) const
    { return samples_.at(i).size(); }

    /*!
     * \brief Return the position on the x-axis of the i-th interval.
     */
    Scalar iToX(unsigned i) const
    {
        assert(0 <= i && i < numX());

        return xPos_.at(i);
    }

    /*!
     * \brief Return the position on the y-axis of the j-th interval.
      */
    Scalar jToY(unsigned i, unsigned j) const
    {
        assert(0 <= i && i < numX());
        assert(0 <= j && size_t(j) < samples_[i].size());

        return std::get<1>(samples_.at(i).at(j));
    }

    /*!
     * \brief Return the interval index of a given position on the x-axis.
     */
    template <class Evaluation>
    unsigned xSegmentIndex(const Evaluation& x, bool extrapolate OPM_OPTIM_UNUSED = false) const
    {
        assert(extrapolate || (xMin() <= x && x <= xMax()));

        // we need at least two sampling points!
        assert(xPos_.size() >= 2);

        if (x <= xPos_[1])
            return 0;
        else if (x >= xPos_[xPos_.size() - 2])
            return xPos_.size() - 2;
        else {
            assert(xPos_.size() >= 3);

            // bisection
            unsigned lowerIdx = 1;
            unsigned upperIdx = xPos_.size() - 2;
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
     * \brief Return the relative position of an x value in an intervall
     *
     * The returned value can be larger than 1 or smaller than zero if it is outside of
     * the range of the segment. In particular this happens for the extrapolation case.
     */
    template <class Evaluation>
    Evaluation xToAlpha(const Evaluation& x, unsigned segmentIdx) const
    {
        Scalar x1 = xPos_[segmentIdx];
        Scalar x2 = xPos_[segmentIdx + 1];
        return (x - x1)/(x2 - x1);
    }

    /*!
     * \brief Return the interval index of a given position on the y-axis.
     */
    template <class Evaluation>
    unsigned ySegmentIndex(const Evaluation& y, unsigned xSampleIdx, bool extrapolate OPM_OPTIM_UNUSED = false) const
    {
        assert(0 <= xSampleIdx && xSampleIdx < numX());
        const auto& colSamplePoints = samples_.at(xSampleIdx);

        assert(colSamplePoints.size() >= 2);
        assert(extrapolate || (yMin(xSampleIdx) <= y && y <= yMax(xSampleIdx)));

        if (y <= std::get<1>(colSamplePoints[1]))
            return 0;
        else if (y >= std::get<1>(colSamplePoints[colSamplePoints.size() - 2]))
            return colSamplePoints.size() - 2;
        else {
            assert(colSamplePoints.size() >= 3);

            // bisection
            unsigned lowerIdx = 1;
            unsigned upperIdx = colSamplePoints.size() - 2;
            while (lowerIdx + 1 < upperIdx) {
                unsigned pivotIdx = (lowerIdx + upperIdx) / 2;
                if (y < std::get<1>(colSamplePoints[pivotIdx]))
                    upperIdx = pivotIdx;
                else
                    lowerIdx = pivotIdx;
            }

            return lowerIdx;
        }
    }

    /*!
     * \brief Return the relative position of an y value in an interval
     *
     * The returned value can be larger than 1 or smaller than zero if it is outside of
     * the range of the segment. In particular this happens for the extrapolation case.
     */
    template <class Evaluation>
    Evaluation yToBeta(const Evaluation& y, unsigned xSampleIdx, unsigned ySegmentIdx) const
    {
        assert(0 <= xSampleIdx && xSampleIdx < numX());
        assert(0 <= ySegmentIdx && ySegmentIdx < numY(xSampleIdx) - 1);

        const auto& colSamplePoints = samples_.at(xSampleIdx);

        Scalar y1 = std::get<1>(colSamplePoints[ySegmentIdx]);
        Scalar y2 = std::get<1>(colSamplePoints[ySegmentIdx + 1]);

        return (y - y1)/(y2 - y1);
    }

    /*!
     * \brief Returns true iff a coordinate lies in the tabulated range
     */
    template <class Evaluation>
    bool applies(const Evaluation& x, const Evaluation& y) const
    {
        if (x < xMin() || xMax() < x)
            return false;

        unsigned i = xSegmentIndex(x, /*extrapolate=*/false);
        Scalar alpha = xToAlpha(Opm::decay<Scalar>(x), i);

        const auto& col1SamplePoints = samples_.at(i);
        const auto& col2SamplePoints = samples_.at(i + 1);

        Scalar minY =
                alpha*std::get<1>(col1SamplePoints.front()) +
                (1 - alpha)*std::get<1>(col2SamplePoints.front());

        Scalar maxY =
                alpha*std::get<1>(col1SamplePoints.back()) +
                (1 - alpha)*std::get<1>(col2SamplePoints.back());

        return minY <= y && y <= maxY;
    }
    /*!
     * \brief Evaluate the function at a given (x,y) position.
     *
     * If this method is called for a value outside of the tabulated
     * range, a \c Opm::NumericalIssue exception is thrown.
     */
    template <class Evaluation>
    Evaluation eval(const Evaluation& x, const Evaluation& y, bool extrapolate=false) const
    {
#ifndef NDEBUG
        if (!extrapolate && !applies(x, y)) {
            std::ostringstream oss;
            oss << "Attempt to get undefined table value (" << x << ", " << y << ")";
            throw NumericalIssue(oss.str());
        };
#endif

        // bi-linear interpolation: first, calculate the x and y indices in the lookup
        // table ...
        unsigned i = xSegmentIndex(x, extrapolate);
        const Evaluation& alpha = xToAlpha(x, i);

        unsigned j1 = ySegmentIndex(y, i, extrapolate);
        unsigned j2 = ySegmentIndex(y, i + 1, extrapolate);
        const Evaluation& beta1 = yToBeta(y, i, j1);
        const Evaluation& beta2 = yToBeta(y, i + 1, j2);

        // evaluate the two function values for the same y value ...
        const Evaluation& s1 = valueAt(i, j1)*(1.0 - beta1) + valueAt(i, j1 + 1)*beta1;
        const Evaluation& s2 = valueAt(i + 1, j2)*(1.0 - beta2) + valueAt(i + 1, j2 + 1)*beta2;

        Valgrind::CheckDefined(s1);
        Valgrind::CheckDefined(s2);

        // ... and combine them using the x position
        const Evaluation& result = s1*(1.0 - alpha) + s2*alpha;
        Valgrind::CheckDefined(result);

        return result;
    }

    /*!
     * \brief Set the x-position of a vertical line.
     *
     * Returns the i index of that line.
     */
    size_t appendXPos(Scalar nextX)
    {
        if (xPos_.empty() || xPos_.back() < nextX) {
            xPos_.push_back(nextX);
            samples_.resize(xPos_.size());
            return xPos_.size() - 1;
        }
        else if (xPos_.front() > nextX) {
            // this is slow, but so what?
            xPos_.insert(xPos_.begin(), nextX);
            samples_.insert(samples_.begin(), std::vector<SamplePoint>());
            return 0;
        }
        throw std::invalid_argument("Sampling points should be specified either monotonically "
                                    "ascending or descending.");
    }

    /*!
     * \brief Append a sample point.
     *
     * Returns the i index of that line.
     */
    size_t appendSamplePoint(size_t i, Scalar y, Scalar value)
    {
        assert(0 <= i && i < numX());

        Scalar x = iToX(i);
        if (samples_[i].empty() || std::get<1>(samples_[i].back()) < y) {
            samples_[i].push_back(SamplePoint(x, y, value));
            return samples_[i].size() - 1;
        }
        else if (std::get<1>(samples_[i].front()) > y) {
            // slow, but we still don't care...
            samples_[i].insert(samples_[i].begin(), SamplePoint(x, y, value));
            return 0;
        }

        throw std::invalid_argument("Sampling points must be specified in either monotonically "
                                    "ascending or descending order.");
    }

    /*!
     * \brief Print the table for debugging purposes.
     *
     * It will produce the data in CSV format on stdout, so that it can be visualized
     * using e.g. gnuplot.
     */
    void print(std::ostream& os = std::cout) const
    {
        Scalar x0 = xMin();
        Scalar x1 = xMax();
        int m = numX();

        Scalar y0 = 1e30;
        Scalar y1 = -1e30;
        int n = 0;
        for (int i = 0; i < m; ++ i) {
            y0 = std::min(y0, yMin(i));
            y1 = std::max(y1, yMax(i));
            n = std::max(n, numY(i));
        }

        m *= 3;
        n *= 3;
        for (int i = 0; i <= m; ++i) {
            Scalar x = x0 + (x1 - x0)*i/m;
            for (int j = 0; j <= n; ++j) {
                Scalar y = y0 + (y1 - y0)*j/n;
                os << x << " " << y << " " << eval(x, y) << "\n";
            }
            os << "\n";
        }
    }

private:
    // the vector which contains the values of the sample points
    // f(x_i, y_j). don't use this directly, use getSamplePoint(i,j)
    // instead!
    std::vector<std::vector<SamplePoint> > samples_;

    // the position of each vertical line on the x-axis
    std::vector<Scalar> xPos_;
};
} // namespace Opm

#endif
