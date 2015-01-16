/*
  Copyright (C) 2015 by Andreas Lauser

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
 * \copydoc Opm::Tabulated1DFunction
 */
#ifndef OPM_TABULATED_1D_FUNCTION_HPP
#define OPM_TABULATED_1D_FUNCTION_HPP

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
    Tabulated1DFunction(int nSamples,
                        const ScalarArrayX &x,
                        const ScalarArrayY &y,
                        bool sortInputs = false)
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
    Tabulated1DFunction(const ScalarContainer &x,
                        const ScalarContainer &y,
                        bool sortInputs = false)
    { this->setXYContainers(x, y, sortInputs); }

    /*!
     * \brief Convenience constructor for a piecewise linear function.
     *
     * \param points A container of \f$(x,y)\f$ tuples of the spline's sampling points (must
     *               have a size() method)
     */
    template <class PointContainer>
    Tabulated1DFunction(const PointContainer &points,
                        bool sortInputs = false)
    { this->setContainerOfPoints(points, sortInputs); }

    /*!
     * \brief Set the sampling points for the piecewise linear function
     *
     * This method takes C-style arrays (pointers) as arguments.
     */
    template <class ScalarArrayX, class ScalarArrayY>
    void setXYArrays(int nSamples,
                     const ScalarArrayX &x,
                     const ScalarArrayY &y,
                     bool sortInputs = false)
    {
        assert(nSamples > 1);

        resizeArrays_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
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
    void setXYContainers(const ScalarContainerX &x,
                         const ScalarContainerY &y,
                         bool sortInputs = false)
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
    void setArrayOfPoints(int nSamples,
                          const PointArray &points,
                          bool sortInputs = false)
    {
        // a linear function with less than two sampling points? what an incredible
        // bad idea!
        assert(nSamples > 1);

        resizeArrays_(nSamples);
        for (int i = 0; i < nSamples; ++i) {
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
    void setContainerOfTuples(const XYContainer &points,
                              bool sortInputs = false)
    {
        // a linear function with less than two sampling points? what an incredible
        // bad idea!
        assert(points.size() > 1);

        resizeArrays_(points.size());
        typename XYContainer::const_iterator it = points.begin();
        typename XYContainer::const_iterator endIt = points.end();
        for (int i = 0; it != endIt; ++i, ++it) {
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
    int numSamples() const
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
    Scalar xAt(int i) const
    { return xValues_[i]; }

    /*!
     * \brief Return the value of the a sample point with a given index.
     */
    Scalar valueAt(int i) const
    { return yValues_[i]; }

    /*!
     * \brief Return true iff the given x is in range [x1, xn].
     */
    bool applies(Scalar x) const
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
    Scalar eval(Scalar x, bool extrapolate=false) const
    {
        int segIdx;
        if (extrapolate && x < xValues_.front())
            segIdx = 0;
        else if (extrapolate && x > xValues_.back())
            segIdx = numSamples() - 2;
        else {
            assert(xValues_.front() <= x && x <= xValues_.back());
            segIdx = findSegmentIndex_(x);
        }

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
    Scalar evalDerivative(Scalar x, bool extrapolate=false) const
    {
        int segIdx = findSegmentIndex_(x);

        return evalDerivative(x, segIdx);
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
    Scalar evalSecondDerivative(Scalar x, bool extrapolate=false) const
    {
        assert(extrapolate || applies(x));
        return 0.0;
    }

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
    Scalar evalThirdDerivative(Scalar x, bool extrapolate=false) const
    {
        assert(extrapolate || applies(x));
        return 0.0;
    }

    /*!
     * \brief Returns 1 if the function is monotonically increasing, -1
     *        if the function is mononously decreasing and 0 if the
     *        function is not monotonous in the interval (x0, x1).
     *
     * In the corner case that the function is constant within the given
     * interval, this method returns 3.
     */
    int monotonic(Scalar x0, Scalar x1, bool extrapolate=false) const
    {
        assert(x0 != x1);

        // make sure that x0 is smaller than x1
        if (x0 > x1)
            std::swap(x0, x1);

        assert(x0 < x1);

        int r = 3;
        if (x0 < xMin()) {
            static_cast<void>(extrapolate);
            assert(extrapolate);

            x0 = xMin();
        };

        int i = findSegmentIndex_(x0);
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
        int iEnd = findSegmentIndex_(x1);
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
    void printCSV(Scalar xi0, Scalar xi1, int k, std::ostream &os = std::cout) const
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
                    OPM_THROW(std::runtime_error,
                              "The sampling points given to a function must be sorted by their x value!");
                }
            }
            else {
                y = eval(x);
                dy_dx = evalDerivative(x);
            }

            os << x << " " << y << " " << dy_dx << "\n";
        }
    }

private:
    int findSegmentIndex_(Scalar x) const
    {
        int n = xValues_.size() - 1;
        assert(n >= 1); // we need at least two sampling points!
        if (xValues_[n] < x)
            return n - 1;
        else if (xValues_[0] > x)
            return 0;

        // bisection
        int lowIdx = 0, highIdx = n;
        while (lowIdx + 1 < highIdx) {
            int curIdx = (lowIdx + highIdx)/2;
            if (xValues_[curIdx] < x)
                lowIdx = curIdx;
            else
                highIdx = curIdx;
        }

        return lowIdx;
    }

    Scalar evalDerivative_(Scalar x, int segIdx) const
    {
        Scalar x0 = xValues_[segIdx];
        Scalar x1 = xValues_[segIdx + 1];

        Scalar y0 = yValues_[segIdx];
        Scalar y1 = yValues_[segIdx + 1];

        return (y1 - y0)/(x1 - x0);
    }

    // returns the monotonicity of a segment
    //
    // The return value have the following meaning:
    //
    // 3: function is constant within interval [x0, x1]
    // 1: function is monotonously increasing in the specified interval
    // 9: function is not monotonic in the specified interval
    // -1: function is monotonously decreasing in the specified interval
    int updateMonotonicity_(int i, int &r) const
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
        ComparatorX_(const std::vector<Scalar> &x)
            : x_(x)
        {};

        bool operator ()(int idxA, int idxB) const
        { return x_.at(idxA) < x_.at(idxB); }

        const std::vector<Scalar> &x_;
    };

    /*!
     * \brief Sort the sample points in ascending order of their x value.
     */
    void sortInput_()
    {
        int n = numSamples();

        // create a vector containing 0...n-1
        std::vector<int> idxVector(n);
        for (int i = 0; i < n; ++i)
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
        int n = numSamples();
        for (int i = 0; i <= (n - 1)/2; ++i) {
            std::swap(xValues_[i], xValues_[n - i - 1]);
            std::swap(yValues_[i], yValues_[n - i - 1]);
        }
    }

    /*!
     * \brief Resizes the internal vectors to store the sample points.
     */
    void resizeArrays_(int nSamples)
    {
        xValues_.resize(nSamples);
        yValues_.resize(nSamples);
    }

    std::vector<Scalar> xValues_;
    std::vector<Scalar> yValues_;
};
} // namespace Opm

#endif
