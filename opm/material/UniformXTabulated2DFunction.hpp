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
 * \copydoc Opm::UniformXTabulated2DFunction
 */
#ifndef OPM_UNIFORM_X_TABULATED_2D_FUNCTION_HPP
#define OPM_UNIFORM_X_TABULATED_2D_FUNCTION_HPP

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <iostream>
#include <vector>
#include <limits>
#include <tuple>

#include <assert.h>

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
     * \brief Returns the number of sampling points in X direction.
     */
    int numX() const
    { return xPos_.size(); }

    /*!
     * \brief Returns the minimum of the Y coordinate of the sampling points for a given column.
     */
    Scalar yMin(int i) const
    { return std::get<1>(samples_.at(i).front()); }

    /*!
     * \brief Returns the maximum of the Y coordinate of the sampling points for a given column.
     */
    Scalar yMax(int i) const
    { return std::get<1>(samples_.at(i).back()); }

    /*!
     * \brief Returns the number of sampling points in Y direction a given column.
     */
    int numY(int i) const
    { return samples_.at(i).size(); }

    /*!
     * \brief Return the position on the x-axis of the i-th interval.
     */
    Scalar iToX(int i) const
    {
        assert(0 <= i && i < numX());

        return xPos_.at(i);
    }

    /*!
     * \brief Return the position on the y-axis of the j-th interval.
      */
    Scalar jToY(int i, int j) const
    {
        assert(0 <= i && i < numX());
        assert(0 <= j && size_t(j) < samples_[i].size());

        return std::get<1>(samples_.at(i).at(j));
    }

    /*!
     * \brief Return the interval index of a given position on the x-axis.
     *
     * This method returns a *floating point* number. The integer part
     * should be interpreted as interval, the decimal places are the
     * position of the x value between the i-th and the (i+1)-th
     * sample point.
      */
    Scalar xToI(Scalar x, bool extrapolate = false) const
    {
        assert(extrapolate || (xMin() <= x && x <= xMax()));

        // interval halving
        int lowerIdx = 0;
        int upperIdx = xPos_.size() - 2;
        int pivotIdx = (lowerIdx + upperIdx) / 2;
        while (lowerIdx + 1 < upperIdx) {
            if (x < xPos_[pivotIdx])
                upperIdx = pivotIdx;
            else
                lowerIdx = pivotIdx;

            pivotIdx = (lowerIdx + upperIdx) / 2;
        }

        Scalar x1 = xPos_[lowerIdx];
        Scalar x2 = xPos_[lowerIdx + 1];
        return lowerIdx + (x - x1)/(x2 - x1);
    }

    /*!
     * \brief Return the interval index of a given position on the y-axis.
     *
     * This method returns a *floating point* number. The integer part
     * should be interpreted as interval, the decimal places are the
     * position of the y value between the j-th and the (j+1)-th
     * sample point.
     */
    Scalar yToJ(int i, Scalar y, bool extrapolate = false) const
    {
        assert(0 <= i && i < numX());
        const auto &colSamplePoints = samples_.at(i);

        assert(extrapolate || (yMin(i) <= y && y <= yMax(i)));

        // interval halving
        int lowerIdx = 0;
        int upperIdx = colSamplePoints.size() - 2;
        int pivotIdx = (lowerIdx + upperIdx) / 2;
        while (lowerIdx + 1 < upperIdx) {
            if (y < std::get<1>(colSamplePoints[pivotIdx]))
                upperIdx = pivotIdx;
            else
                lowerIdx = pivotIdx;

            pivotIdx = (lowerIdx + upperIdx) / 2;
        }

        Scalar y1 = std::get<1>(colSamplePoints[lowerIdx]);
        Scalar y2 = std::get<1>(colSamplePoints[lowerIdx + 1]);
        return lowerIdx + (y - y1)/(y2 - y1);
    }

    /*!
     * \brief Returns true iff a coordinate lies in the tabulated range
     */
    bool applies(Scalar x, Scalar y) const
    {
        if (x < xMin() || xMax() < x)
            return false;

        Scalar i = xToI(x, /*extrapolate=*/false);
        const auto &col1SamplePoints = samples_.at(int(i));
        const auto &col2SamplePoints = samples_.at(int(i));
        Scalar alpha = i - int(i);

        Scalar yMin =
                alpha*std::get<1>(col1SamplePoints.front()) +
                (1 - alpha)*std::get<1>(col2SamplePoints.front());

        Scalar yMax =
                alpha*std::get<1>(col1SamplePoints.back()) +
                (1 - alpha)*std::get<1>(col2SamplePoints.back());

        return yMin <= y && y <= yMax;
    }

    /*!
     * \brief Evaluate the function at a given (x,y) position.
     *
     * If this method is called for a value outside of the tabulated
     * range, a \c Opm::NumericalProblem exception is thrown.
     */
    Scalar eval(Scalar x, Scalar y, bool extrapolate = true) const
    {
#ifndef NDEBUG
        if (!extrapolate && !applies(x,y))
        {
            OPM_THROW(NumericalProblem,
                       "Attempt to get tabulated value for ("
                       << x << ", " << y
                       << ") on table");
        };
#endif

        Scalar alpha = xToI(x, extrapolate);
        int i = std::max(0, std::min(numX() - 2, static_cast<int>(alpha)));
        alpha -= i;

        Scalar beta1 = yToJ(i, y, extrapolate);
        Scalar beta2 = yToJ(i + 1, y, extrapolate);

        int j1 = std::max(0, std::min(numY(i) - 2, static_cast<int>(beta1)));
        int j2 = std::max(0, std::min(numY(i + 1) - 2, static_cast<int>(beta2)));

        beta1 -= j1;
        beta2 -= j2;

        // bi-linear interpolation
        Scalar s1 = getSamplePoint(i, j1)*(1.0 - beta1) + getSamplePoint(i, j1 + 1)*beta1;
        Scalar s2 = getSamplePoint(i + 1, j2)*(1.0 - beta2) + getSamplePoint(i + 1, j2 + 1)*beta2;
        return s1*(1.0 - alpha) + s2*alpha;
    }

    /*!
     * \brief Get the value of the sample point which is at the
     *         intersection of the \f$i\f$-th interval of the x-Axis
     *         and the \f$j\f$-th of the y-Axis.
     */
    Scalar getSamplePoint(int i, int j) const
    {
        assert(0 <= i && i < numX());

        const auto &colSamples = samples_[i];

        assert(0 <= j && size_t(j) < colSamples.size());

        return std::get<2>(colSamples[j]);
    }

    /*!
     * \brief Set the x-position of a vertical line.
     */
    void appendXPos(Scalar nextX)
    {
#ifndef NDEBUG
        if (xPos_.size())
            assert(xPos_.back() < nextX);
#endif

        xPos_.push_back(nextX);
        samples_.resize(xPos_.size());
    }

    /*!
     * \brief Append a sample point.
     */
    void appendSamplePoint(int i, Scalar y, Scalar value)
    {
        assert(0 <= i && i < numX());

        Scalar x = iToX(i);
        samples_[i].push_back(SamplePoint(x, y, value));
    }

    /*!
     * \brief Print the table for debugging purposes.
     *
     * It will produce the data in CSV format on stdout, so that it can be visualized
     * using e.g. gnuplot.
     */
    void print(std::ostream &os = std::cout) const
    {
        Scalar x0 = xMin();
        Scalar x1 = xMax();
        int m = numX();

        Scalar y0 = 1e100;
        Scalar y1 = -1e100;
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
