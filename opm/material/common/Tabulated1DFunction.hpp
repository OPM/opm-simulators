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
 * \copydoc Opm::Tabulated1DFunction
 */
#ifndef OPM_TABULATED_1D_FUNCTION_HPP
#define OPM_TABULATED_1D_FUNCTION_HPP

#include <opm/material/densead/Math.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <tuple>
#include <vector>

namespace Opm {
/*!
 * \brief Implements a linearly interpolated scalar function that depends on one
 *        variable.
 */
template <class Scalar>
class Tabulated1DFunction
{
public:
    /*!
     * \brief Default constructor for a piecewise linear function.
     *
     * To specfiy the acutal curve, use one of the set() methods.
     */
    Tabulated1DFunction()
    {}

    /*!
     * \brief Convenience constructor for a piecewise linear function.
     *
     * \param nSamples The number of sampling points (must be >= 2)
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     */
    template <class ScalarArrayX, class ScalarArrayY>
    Tabulated1DFunction(size_t nSamples,
                        const ScalarArrayX& x,
                        const ScalarArrayY& y,
                        bool sortInputs = true)
    { this->setXYArrays(nSamples, x, y, sortInputs); }

    /*!
     * \brief Convenience constructor for a piecewise linear function.
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     *          (must have a size() method)
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     *          (must have a size() method)
     */
    template <class ScalarContainer>
    Tabulated1DFunction(const ScalarContainer& x,
                        const ScalarContainer& y,
                        bool sortInputs = true)
    { this->setXYContainers(x, y, sortInputs); }

    /*!
     * \brief Convenience constructor for a piecewise linear function.
     *
     * \param points A container of \f$(x,y)\f$ tuples of the spline's sampling points (must
     *               have a size() method)
     */
    template <class PointContainer>
    Tabulated1DFunction(const PointContainer& points,
                        bool sortInputs = true)
    { this->setContainerOfTuples(points, sortInputs); }

    /*!
     * \brief Set the sampling points for the piecewise linear function
     *
     * This method takes C-style arrays (pointers) as arguments.
     */
    template <class ScalarArrayX, class ScalarArrayY>
    void setXYArrays(size_t nSamples,
                     const ScalarArrayX& x,
                     const ScalarArrayY& y,
                     bool sortInputs = true)
    {
        assert(nSamples > 1);

        resizeArrays_(nSamples);
        for (size_t i = 0; i < nSamples; ++i) {
            xValues_[i] = x[i];
            yValues_[i] = y[i];
        }

        if (sortInputs)
            sortInput_();
        else if (xValues_[0] > xValues_[numSamples() - 1])
            reverseSamplingPoints_();
    }

    /*!
     * \brief Set the sampling points for the piecewise linear function
     *
     * This method takes STL compatible containers as arguments.
     */
    template <class ScalarContainerX, class ScalarContainerY>
    void setXYContainers(const ScalarContainerX& x,
                         const ScalarContainerY& y,
                         bool sortInputs = true)
    {
        assert(x.size() == y.size());
        assert(x.size() > 1);

        resizeArrays_(x.size());
        std::copy(x.begin(), x.end(), xValues_.begin());
        std::copy(y.begin(), y.end(), yValues_.begin());

        if (sortInputs)
            sortInput_();
        else if (xValues_[0] > xValues_[numSamples() - 1])
            reverseSamplingPoints_();
    }

    /*!
     * \brief Set the sampling points for the piecewise linear function
     */
    template <class PointArray>
    void setArrayOfPoints(size_t nSamples,
                          const PointArray& points,
                          bool sortInputs = true)
    {
        // a linear function with less than two sampling points? what an incredible
        // bad idea!
        assert(nSamples > 1);

        resizeArrays_(nSamples);
        for (size_t i = 0; i < nSamples; ++i) {
            xValues_[i] = points[i][0];
            yValues_[i] = points[i][1];
        }

        if (sortInputs)
            sortInput_();
        else if (xValues_[0] > xValues_[numSamples() - 1])
            reverseSamplingPoints_();
    }

    /*!
     * \brief Set the sampling points of the piecewise linear function using a
     *        STL-compatible container of tuple-like objects.
     *
     * This method uses a single STL-compatible container of sampling
     * points, which are assumed to be tuple-like objects storing the
     * X and Y coordinates.  "tuple-like" means that the objects
     * provide access to the x values via std::get<0>(obj) and to the
     * y value via std::get<1>(obj) (e.g. std::tuple or
     * std::pair). "STL-compatible" means that the container provides
     * access to iterators using the begin(), end() methods and also
     * provides a size() method. Also, the number of entries in the X
     * and the Y containers must be equal and larger than 1.
     */
    template <class XYContainer>
    void setContainerOfTuples(const XYContainer& points,
                              bool sortInputs = true)
    {
        // a linear function with less than two sampling points? what an incredible
        // bad idea!
        assert(points.size() > 1);

        resizeArrays_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (unsigned i = 0; it != endIt; ++i, ++it) {
            xValues_[i] = std::get<0>(*it);
            yValues_[i] = std::get<1>(*it);
        }

        if (sortInputs)
            sortInput_();
        else if (xValues_[0] > xValues_[numSamples() - 1])
            reverseSamplingPoints_();
    }

    /*!
     * \brief Returns the number of sampling points.
     */
    size_t numSamples() const
    { return xValues_.size(); }

    /*!
     * \brief Return the x value of the leftmost sampling point.
     */
    Scalar xMin() const
    { return xValues_[0]; }

    /*!
     * \brief Return the x value of the rightmost sampling point.
     */
    Scalar xMax() const
    { return xValues_[numSamples() - 1]; }

    /*!
     * \brief Return the x value of the a sample point with a given index.
     */
    Scalar xAt(size_t i) const
    { return xValues_[i]; }

    const std::vector<Scalar>& xValues() const
    { return xValues_; }

    const std::vector<Scalar>& yValues() const
    { return yValues_; }

    /*!
     * \brief Return the value of the a sample point with a given index.
     */
    Scalar valueAt(size_t i) const
    { return yValues_[i]; }

    /*!
     * \brief Return true iff the given x is in range [x1, xn].
     */
    template <class Evaluation>
    bool applies(const Evaluation& x) const
    { return xValues_[0] <= x && x <= xValues_[numSamples() - 1]; }

    /*!
     * \brief Evaluate the spline at a given position.
     *
     * \param x The value on the abscissa where the function ought to be evaluated
     * \param extrapolate If this parameter is set to true, the function will be extended
     *                    beyond its range by straight lines, if false calling
     *                    extrapolate for \f$ x \not [x_{min}, x_{max}]\f$ will cause a
     *                    failed assertation.
     */
    template <class Evaluation>
    Evaluation eval(const Evaluation& x, bool extrapolate = false) const
    {
        size_t segIdx = findSegmentIndex_(x, extrapolate);

        Scalar x0 = xValues_[segIdx];
        Scalar x1 = xValues_[segIdx + 1];

        Scalar y0 = yValues_[segIdx];
        Scalar y1 = yValues_[segIdx + 1];

        return y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    }

    /*!
     * \brief Evaluate the spline's derivative at a given position.
     *
     * \param x The value on the abscissa where the function's
     *          derivative ought to be evaluated
     * \param extrapolate If this parameter is set to true, the spline
     *                    will be extended beyond its range by
     *                    straight lines, if false calling extrapolate
     *                    for \f$ x \not [x_{min}, x_{max}]\f$ will
     *                    cause a failed assertation.
     */
    template <class Evaluation>
    Evaluation evalDerivative(const Evaluation& x, bool extrapolate = false) const
    {
        unsigned segIdx = findSegmentIndex_(x, extrapolate);
        return evalDerivative_(x, segIdx);
    }

    /*!
     * \brief Evaluate the function's second derivative at a given position.
     *
     * Since this class represents a piecewise linear function, this method will always
     * return 0.
     *
     * \param x The value on the abscissa where the function's
     *          derivative ought to be evaluated
     * \param extrapolate If this parameter is set to true, the function
     *                    will be extended beyond its range by
     *                    straight lines, if false calling extrapolate
     *                    for \f$ x \not [x_{min}, x_{max}]\f$ will
     *                    cause a failed assertation.
     */
    template <class Evaluation>
    Evaluation evalSecondDerivative(const Evaluation& x OPM_UNUSED, bool extrapolate OPM_UNUSED = false) const
    { return 0.0; }

    /*!
     * \brief Evaluate the function's third derivative at a given position.
     *
     * Since this class represents a piecewise linear function, this method will always
     * return 0.
     *
     * \param x The value on the abscissa where the function's
     *          derivative ought to be evaluated
     * \param extrapolate If this parameter is set to true, the function
     *                    will be extended beyond its range by
     *                    straight lines, if false calling extrapolate
     *                    for \f$ x \not [x_{min}, x_{max}]\f$ will
     *                    cause a failed assertation.
     */
    template <class Evaluation>
    Evaluation evalThirdDerivative(const Evaluation& x OPM_UNUSED, bool extrapolate OPM_UNUSED = false) const
    { return 0.0; }

    /*!
     * \brief Returns 1 if the function is monotonically increasing, -1
     *        if the function is mononously decreasing and 0 if the
     *        function is not monotonous in the interval (x0, x1).
     *
     * In the corner case that the function is constant within the given
     * interval, this method returns 3.
     */
    int monotonic(Scalar x0, Scalar x1, bool extrapolate OPM_OPTIM_UNUSED = false) const
    {
        assert(x0 != x1);

        // make sure that x0 is smaller than x1
        if (x0 > x1)
            std::swap(x0, x1);

        assert(x0 < x1);

        int r = 3;
        if (x0 < xMin()) {
            assert(extrapolate);

            x0 = xMin();
        };

        size_t i = findSegmentIndex_(x0, extrapolate);
        if (xValues_[i + 1] >= x1) {
            // interval is fully contained within a single function
            // segment
            updateMonotonicity_(i, r);
            return r;
        }

        // the first segment overlaps with the specified interval
        // partially
        updateMonotonicity_(i, r);
        ++ i;

        // make sure that the segments which are completly in the
        // interval [x0, x1] all exhibit the same monotonicity.
        size_t iEnd = findSegmentIndex_(x1, extrapolate);
        for (; i < iEnd - 1; ++i) {
            updateMonotonicity_(i, r);
            if (!r)
                return 0;
        }

        // if the user asked for a part of the function which is
        // extrapolated, we need to check the slope at the function's
        // endpoint
        if (x1 > xMax()) {
            assert(extrapolate);

            Scalar m = evalDerivative_(xMax(), /*segmentIdx=*/numSamples() - 2);
            if (m < 0)
                return (r < 0 || r==3)?-1:0;
            else if (m > 0)
                return (r > 0 || r==3)?1:0;

            return r;
        }

        // check for the last segment
        updateMonotonicity_(iEnd, r);

        return r;
    }

    /*!
     * \brief Same as monotonic(x0, x1), but with the entire range of the
     *        function as interval.
     */
    int monotonic() const
    { return monotonic(xMin(), xMax()); }

    /*!
     * \brief Prints k tuples of the format (x, y, dx/dy, isMonotonic)
     *        to stdout.
     *
     * If the function does not apply for parts of [x0, x1] it is
     * extrapolated using a straight line. The result can be inspected
     * using the following commands:
     *
     ----------- snip -----------
     ./yourProgramm > function.csv
     gnuplot

     gnuplot> plot "function.csv" using 1:2 w l ti "Curve", \
     "function.csv" using 1:3 w l ti "Derivative"
     ----------- snap -----------
    */
    void printCSV(Scalar xi0, Scalar xi1, unsigned k, std::ostream& os = std::cout) const
    {
        Scalar x0 = std::min(xi0, xi1);
        Scalar x1 = std::max(xi0, xi1);
        const int n = numSamples() - 1;
        for (int i = 0; i <= k; ++i) {
            double x = i*(x1 - x0)/k + x0;
            double y;
            double dy_dx;
            if (!applies(x)) {
                if (x < xValues_[0]) {
                    dy_dx = evalDerivative(xValues_[0]);
                    y = (x - xValues_[0])*dy_dx + yValues_[0];
                }
                else if (x > xValues_[n]) {
                    dy_dx = evalDerivative(xValues_[n]);
                    y = (x - xValues_[n])*dy_dx + yValues_[n];
                }
                else {
                    throw std::runtime_error("The sampling points given to a function must be sorted by their x value!");
                }
            }
            else {
                y = eval(x);
                dy_dx = evalDerivative(x);
            }

            os << x << " " << y << " " << dy_dx << "\n";
        }
    }

    bool operator==(const Tabulated1DFunction<Scalar>& data) const {
        return xValues_ == data.xValues_ &&
               yValues_ == data.yValues_;
    }

private:
    template <class Evaluation>
    size_t findSegmentIndex_(const Evaluation& x, bool extrapolate = false) const
    {
        if (!extrapolate && !applies(x))
            throw Opm::NumericalIssue("Tried to evaluate a tabulated function outside of its range");

        // we need at least two sampling points!
        assert(xValues_.size() >= 2);

        if (x <= xValues_[1])
            return 0;
        else if (x >= xValues_[xValues_.size() - 2])
            return xValues_.size() - 2;
        else {
            // bisection
            size_t lowerIdx = 1;
            size_t upperIdx = xValues_.size() - 2;
            while (lowerIdx + 1 < upperIdx) {
                size_t pivotIdx = (lowerIdx + upperIdx) / 2;
                if (x < xValues_[pivotIdx])
                    upperIdx = pivotIdx;
                else
                    lowerIdx = pivotIdx;
            }

            assert(xValues_[lowerIdx] <= x);
            assert(x <= xValues_[lowerIdx + 1]);
            return lowerIdx;
        }
    }

    template <class Evaluation>
    Evaluation evalDerivative_(const Evaluation& x, size_t segIdx) const
    {
        Scalar x0 = xValues_[segIdx];
        Scalar x1 = xValues_[segIdx + 1];

        Scalar y0 = yValues_[segIdx];
        Scalar y1 = yValues_[segIdx + 1];

        Evaluation ret = blank(x);
        ret = (y1 - y0)/(x1 - x0);
        return ret;
    }

    // returns the monotonicity of a segment
    //
    // The return value have the following meaning:
    //
    // 3: function is constant within interval [x0, x1]
    // 1: function is monotonously increasing in the specified interval
    // 0: function is not monotonic in the specified interval
    // -1: function is monotonously decreasing in the specified interval
    int updateMonotonicity_(size_t i, int& r) const
    {
        if (yValues_[i] < yValues_[i + 1]) {
            // monotonically increasing?
            if (r == 3 || r == 1)
                r = 1;
            else
                r = 0;
            return 1;
        }
        else if (yValues_[i] > yValues_[i + 1]) {
            // monotonically decreasing?
            if (r == 3 || r == -1)
                r = -1;
            else
                r = 0;
            return -1;
        }

        return 3;
    }

    /*!
     * \brief Helper class needed to sort the input sampling points.
     */
    struct ComparatorX_
    {
        ComparatorX_(const std::vector<Scalar>& x)
            : x_(x)
        {}

        bool operator ()(size_t idxA, size_t idxB) const
        { return x_.at(idxA) < x_.at(idxB); }

        const std::vector<Scalar>& x_;
    };

    /*!
     * \brief Sort the sample points in ascending order of their x value.
     */
    void sortInput_()
    {
        size_t n = numSamples();

        // create a vector containing 0...n-1
        std::vector<unsigned> idxVector(n);
        for (unsigned i = 0; i < n; ++i)
            idxVector[i] = i;

        // sort the indices according to the x values of the sample
        // points
        ComparatorX_ cmp(xValues_);
        std::sort(idxVector.begin(), idxVector.end(), cmp);

        // reorder the sample points
        std::vector<Scalar> tmpX(n), tmpY(n);
        for (size_t i = 0; i < idxVector.size(); ++ i) {
            tmpX[i] = xValues_[idxVector[i]];
            tmpY[i] = yValues_[idxVector[i]];
        }
        xValues_ = tmpX;
        yValues_ = tmpY;
    }

    /*!
     * \brief Reverse order of the elements in the arrays which
     *        contain the sampling points.
     */
    void reverseSamplingPoints_()
    {
        // reverse the arrays
        size_t n = numSamples();
        for (size_t i = 0; i <= (n - 1)/2; ++i) {
            std::swap(xValues_[i], xValues_[n - i - 1]);
            std::swap(yValues_[i], yValues_[n - i - 1]);
        }
    }

    /*!
     * \brief Resizes the internal vectors to store the sample points.
     */
    void resizeArrays_(size_t nSamples)
    {
        xValues_.resize(nSamples);
        yValues_.resize(nSamples);
    }

    std::vector<Scalar> xValues_;
    std::vector<Scalar> yValues_;
};
} // namespace Opm

#endif
