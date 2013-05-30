// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
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
 * \copydoc Opm::SplineCommon_
 */
#ifndef OPM_SPLINE_COMMON__HH
#define OPM_SPLINE_COMMON__HH

#include "valgrind.hh"
#include "math.hh"

#include <dune/common/exceptions.hh>

#include <tuple>
#include <algorithm>
#include <iostream>
#include <cassert>

namespace Opm
{
/*!
 * \brief The common code for all 3rd order polynomial splines.
 */
template<class ScalarT, class ImplementationT>
class SplineCommon_
{
    typedef ScalarT Scalar;
    typedef ImplementationT Implementation;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

public:
    /*!
     * \brief Return true iff the given x is in range [x1, xn].
     */
    bool applies(Scalar x) const
    {
        return x_(0) <= x && x <= x_(numSamples_() - 1);
    }

    /*!
     * \brief Return the x value of the leftmost sampling point.
     */
    Scalar xMin() const
    { return x_(0); }

    /*!
     * \brief Return the x value of the rightmost sampling point.
     */
    Scalar xMax() const
    { return x_(numSamples_() - 1); }

    /*!
     * \brief Prints k tuples of the format (x, y, dx/dy, isMonotonic)
     *        to stdout.
     *
     * If the spline does not apply for parts of [x0, x1] it is
     * extrapolated using a straight line. The result can be inspected
     * using the following commands:
     *
     ----------- snip -----------
     ./yourProgramm > spline.csv
     gnuplot

     gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
     "spline.csv" using 1:3 w l ti "Derivative", \
     "spline.csv" using 1:4 w p ti "Monotonic"
     ----------- snap -----------
    */
    void printCSV(Scalar xi0, Scalar xi1, int k) const
    {
        Scalar x0 = std::min(xi0, xi1);
        Scalar x1 = std::max(xi0, xi1);
        const int n = numSamples_() - 1;
        for (int i = 0; i <= k; ++i) {
            double x = i*(x1 - x0)/k + x0;
            double x_p1 = x + (x1 - x0)/k;
            double y;
            double dy_dx;
            double mono = 1;
            if (!applies(x)) {
                if (x < 0) {
                    dy_dx = evalDerivative(x_(0));
                    y = (x - x_(0))*dy_dx + y_(0);
                    mono = (dy_dx>0)?1:-1;
                }
                else if (x > x_(n)) {
                    dy_dx = evalDerivative(x_(n));
                    y = (x - x_(n))*dy_dx + y_(n);
                    mono = (dy_dx>0)?1:-1;
                }
                else {
                    DUNE_THROW(Dune::InvalidStateException,
                               "The sampling points given to a spline must be sorted by their x value!");
                }
            }
            else {
                y = eval(x);
                dy_dx = evalDerivative(x);
                mono = monotonic(std::max<Scalar>(x_(0), x), std::min<Scalar>(x_(n), x_p1));
            }

            std::cout << x << " " << y << " " << dy_dx << " " << mono << "\n";
        }
    }

    /*!
     * \brief Evaluate the spline at a given position.
     *
     * \param x The value on the abscissa where the spline ought to be
     *          evaluated
     * \param extrapolate If this parameter is set to true, the spline
     *                    will be extended beyond its range by
     *                    straight lines, if false calling extrapolate
     *                    for \f$ x \not [x_{min}, x_{max}]\f$ will
     *                    cause a failed assertation.
     */
    Scalar eval(Scalar x, bool extrapolate=false) const
    {
        assert(extrapolate || applies(x));

        // handle extrapolation
        if (extrapolate) {
            if (x < xMin()) {
                Scalar m = evalDerivative(xMin(), /*segmentIdx=*/0);
                Scalar y0 = y_(0);
                return y0 + m*(x - xMin());
            }
            else if (x > xMax()) {
                Scalar m = evalDerivative(xMax(), /*segmentIdx=*/numSamples_()-2);
                Scalar y0 = y_(numSamples_() - 1);
                return y0 + m*(x - xMax());
            }
        }

        return eval_(x, segmentIdx_(x));
    }

    /*!
     * \brief Evaluate the spline's derivative at a given position.
     *
     * \param x The value on the abscissa where the spline's
     *          derivative ought to be evaluated
     *
     * \param extrapolate If this parameter is set to true, the spline
     *                    will be extended beyond its range by
     *                    straight lines, if false calling extrapolate
     *                    for \f$ x \not [x_{min}, x_{max}]\f$ will
     *                    cause a failed assertation.

     */
     Scalar evalDerivative(Scalar x, bool extrapolate=false) const
    {
        assert(extrapolate || applies(x));
        if (extrapolate) {
            if (x < xMin())
                evalDerivative_(xMin(), 0);
            else if (x > xMax())
                evalDerivative_(xMax(), numSamples_() - 2);
        }

        return evalDerivative_(x, segmentIdx_(x));
    }

    /*!
     * \brief Find the intersections of the spline with a cubic
     *        polynomial in the whole intervall, throws
     *        Dune::MathError exception if there is more or less than
     *        one solution.
     */
    Scalar intersect(Scalar a, Scalar b, Scalar c, Scalar d) const
    {
        return intersectIntervall(xMin(), xMax(), a, b, c, d);
    }

    /*!
     * \brief Find the intersections of the spline with a cubic
     *        polynomial in a sub-intervall of the spline, throws
     *        Dune::MathError exception if there is more or less than
     *        one solution.
     */
    Scalar intersectInterval(Scalar x0, Scalar x1,
                             Scalar a, Scalar b, Scalar c, Scalar d) const
    {
        assert(applies(x0) && applies(x1));

        Scalar tmpSol[3];
        int nSol = 0;
        int iFirst = segmentIdx_(x0);
        int iLast = segmentIdx_(x1);
        for (int i = iFirst; i <= iLast; ++i)
        {
            nSol += intersectSegment_(tmpSol, i, a, b, c, d, x0, x1);

            if (nSol > 1) {
                DUNE_THROW(Dune::MathError,
                           "Spline has more than one intersection"); //<<a<<"x^3 + "<<b<"x^2 + "<<c<"x + "<<d);
            }
        }

        if (nSol != 1)
            DUNE_THROW(Dune::MathError,
                       "Spline has no intersection"); //<<a<"x^3 + " <<b<"x^2 + "<<c<"x + "<<d<<"!");

        return tmpSol[0];
    }

    /*!
     * \brief Returns 1 if the spline is monotonically increasing, -1
     *        if the spline is mononously decreasing and 0 if the
     *        spline is not monotonous in the interval (x0, x1).
     *
     * In the corner case where the whole spline is flat, it returns
     * 2.
     */
    int monotonic(Scalar x0, Scalar x1) const
    {
        assert(applies(x0));
        assert(applies(x1));
        assert(x0 != x1);

        // make sure that x0 is smaller than x1
        if (x0 > x1)
            std::swap(x0, x1);

        assert(x0 < x1);

        // corner case where the whole spline is a constant
        if (moment_(0) == 0 &&
            moment_(1) == 0 &&
            y_(0) == y_(1))
        {
            // actually the is monotonically increasing as well as
            // monotonously decreasing
            return 2;
        }


        int i = segmentIdx_(x0);
        if (x_(i) <= x0 && x1 <= x_(i+1)) {
            // the interval to check is completely included in a
            // single segment
            return monotonic_(i, x0, x1);
        }

        // make sure that the segments which are completly in the
        // interval [x0, x1] all exhibit the same monotonicity.
        int iEnd = segmentIdx_(x1);
        int r = monotonic_(i, x0, x_(1));
        for (; i < iEnd - 1; ++i)
            if (r != monotonic_(i, x_(i), x_(i + 1)))
                return 0;

        // check for the last segment
        if (x_(iEnd) < x1 && r != monotonic_(iEnd, x_(iEnd), x1))
        { return 0; }

        return r;
    }

    /*!
     * \brief Same as monotonic(x0, x1), but with the entire range of the
     *        spline as interval.
     */
    int monotonic() const
    { return monotonic(xMin(), xMax()); }

protected:
    // this is an internal class, so everything is protected!
    SplineCommon_()
    { Valgrind::SetUndefined(asImp_()); }

    /*!
     * \brief Set the sampling point vectors.
     *
     * This takes care that the order of the x-values is ascending,
     * although the input must be ordered already!
     */
    template <class DestVector, class SourceVector>
    void assignSamplingPoints_(DestVector &destX,
                               DestVector &destY,
                               const SourceVector &srcX,
                               const SourceVector &srcY,
                               int numSamples)
    {
        assert(numSamples >= 2);

        // copy sample points, make sure that the first x value is
        // smaller than the last one
        for (int i = 0; i < numSamples; ++i) {
            int idx = i;
            if (srcX[0] > srcX[numSamples - 1])
                idx = numSamples - i - 1;
            destX[i] = srcX[idx];
            destY[i] = srcY[idx];
        }
    }

    template <class DestVector, class ListIterator>
    void assignFromArrayList_(DestVector &destX,
                              DestVector &destY,
                              const ListIterator &srcBegin,
                              const ListIterator &srcEnd,
                              int numSamples)
    {
        assert(numSamples >= 2);

        // find out wether the x values are in reverse order
        ListIterator it = srcBegin;
        ++it;
        bool reverse = false;
        if ((*srcBegin)[0] > (*it)[0])
            reverse = true;
        --it;

        // loop over all sampling points
        for (int i = 0; it != srcEnd; ++i, ++it) {
            int idx = i;
            if (reverse)
                idx = numSamples - i - 1;
            destX[i] = (*it)[0];
            destY[i] = (*it)[1];
        }
    }

    /*!
     * \brief Set the sampling points.
     *
     * Here we assume that the elements of the source vector have an
     * [] operator where v[0] is the x value and v[1] is the y value
     * if the sampling point.
     */
    template <class DestVector, class ListIterator>
    void assignFromTupleList_(DestVector &destX,
                              DestVector &destY,
                              ListIterator srcBegin,
                              ListIterator srcEnd,
                              int numSamples)
    {
        assert(numSamples >= 2);

        // copy sample points, make sure that the first x value is
        // smaller than the last one

        // find out wether the x values are in reverse order
        ListIterator it = srcBegin;
        ++it;
        bool reverse = false;
        if (std::get<0>(*srcBegin) > std::get<0>(*it))
            reverse = true;
        --it;

        // loop over all sampling points
        for (int i = 0; it != srcEnd; ++i, ++it) {
            int idx = i;
            if (reverse)
                idx = numSamples - i - 1;
            destX[i] = std::get<0>(*it);
            destY[i] = std::get<1>(*it);
        }
    }


    /*!
     * \brief Make the linear system of equations Mx = d which results
     *        in the moments of the periodic spline.
     *
     * When solving Mx = d, it should be noted that x[0] is trash and
     * needs to be set to x[n-1]
     */
    template <class Vector, class Matrix>
    void makePeriodicSystem_(Matrix &M, Vector &d)
    {
        // a periodic spline only makes sense if the first and the
        // last sampling point have the same y value
        assert(y_(0) == y_(numSamples_() - 1));

        makeNaturalSystem(M, d);

        int n = numSamples_() - 1;

        Scalar lambda_n = h_(1)/(h_(n) + h_(1));
        Scalar lambda_1 = M[1][2];
        Scalar mu_1 = 1 - lambda_1;
        Scalar mu_n = 1 - lambda_n;

        Scalar d_n =
            6 / (h_(n) + h_(1))
            *
            ( (y_(1) - y_(n))/h_(1) - (y_(n) - y_(n-1))/h_(n));

        // insert lambda_n and mu_1
        M[1][n] = mu_1;
        M[n][1] = lambda_n;
        M[n][n-1] = mu_n;
        d[n] = d_n;
    }

    /*!
     * \brief Make the linear system of equations Mx = d which results
     *        in the moments of the full spline.
     */
    template <class Vector, class Matrix>
    void makeFullSystem_(Matrix &M, Vector &d, Scalar m0, Scalar m1)
    {
        makeNaturalSystem_(M, d);

        int n = numSamples_() - 1;
        // first row
        M[0][1] = 1;
        d[0] = 6/h_(1) * ( (y_(1) - y_(0))/h_(1) - m0);

        // last row
        M[n][n - 1] = 1;

        // right hand side
        d[n] =
            6/h_(n)
            *
            (m1 - (y_(n) - y_(n - 1))/h_(n));
    }

    /*!
     * \brief Make the linear system of equations Mx = d which results
     *        in the moments of the natural spline.
     */
    template <class Vector, class Matrix>
    void makeNaturalSystem_(Matrix &M, Vector &d)
    {
        M = 0;
        d = 0;

        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 111
        const int n = numSamples_() - 1;

        // second to next to last rows
        for (int i = 1; i < n; ++i) {
            const Scalar lambda_i = h_(i + 1) / (h_(i) + h_(i+1));
            const Scalar mu_i = 1 - lambda_i;
            const Scalar d_i =
                6 / (h_(i) + h_(i+1))
                *
                ( (y_(i+1) - y_(i))/h_(i+1) - (y_(i) - y_(i-1))/h_(i));

            M[i][i-1] = mu_i;
            M[i][i] = 2;
            M[i][i + 1] = lambda_i;
            d[i] = d_i;
        };

        // first row
        M[0][0] = 2;

        // last row
        M[n][n] = 2;
        // right hand side
        d[0] = 0.0;
        d[n] = 0.0;
    }

    // evaluate the spline at a given the position and given the
    // segment index
    Scalar eval_(Scalar x, int i) const
    {
        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 109
        Scalar h_i1 = h_(i + 1);
        Scalar x_i = x - x_(i);
        Scalar x_i1 = x_(i+1) - x;

        Scalar A_i =
            (y_(i+1) - y_(i))/h_i1
            -
            h_i1/6*(moment_(i+1) - moment_(i));
        Scalar B_i = y_(i) - moment_(i)* (h_i1*h_i1) / 6;

        return
            moment_(i)* x_i1*x_i1*x_i1 / (6 * h_i1)
            +
            moment_(i + 1)* x_i*x_i*x_i / (6 * h_i1)
            +
            A_i*x_i
            +
            B_i;
    }

    // evaluate the derivative of a spline given the actual position
    // and the segment index
    Scalar evalDerivative_(Scalar x, int i) const
    {
        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 109
        Scalar h_i1 = h_(i + 1);
        Scalar x_i = x - x_(i);
        Scalar x_i1 = x_(i+1) - x;

        Scalar A_i =
            (y_(i+1) - y_(i))/h_i1
            -
            h_i1/6*(moment_(i+1) - moment_(i));

        return
            -moment_(i) * x_i1*x_i1 / (2 * h_i1)
            +
            moment_(i + 1) * x_i*x_i / (2 * h_i1)
            +
            A_i;
    }


    // returns the monotonicality of an interval of a spline segment
    //
    // The return value have the following meaning:
    //
    // 1: spline is monotonously increasing in the specified interval
    // 0: spline is not monotonic in the specified interval
    // -1: spline is monotonously decreasing in the specified interval
    int monotonic_(int i, Scalar x0, Scalar x1) const
    {
        // shift the interval so that it is consistent with the
        // definitions by Stoer
        x0 = x0 - x_(i);
        x1 = x1 - x_(i);

        Scalar a = a_(i);
        Scalar b = b_(i);
        Scalar c = c_(i);

        Scalar disc = 4*b*b - 12*a*c;
        if (disc < 0) {
            // discriminant is smaller than 0, i.e. the segment does
            // not exhibit any extrema.
            return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
        }
        disc = std::sqrt(disc);
        Scalar xE1 = (-2*b + disc)/(6*a);
        Scalar xE2 = (-2*b - disc)/(6*a);

        if (disc == 0) {
            // saddle point -> no extrema
            if (xE1 == x0)
                // make sure that we're not picking the saddle point
                // to determine whether we're monotonically increasing
                // or decreasing
                x0 = x1;
            return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
        };
        if ((x0 < xE1 && xE1 < x1) ||
            (x0 < xE2 && xE2 < x1))
        {
            // there is an extremum in the range (x0, x1)
            return 0;
        }
        // no extremum in range (x0, x1)
        x0 = (x0 + x1)/2; // pick point in the middle of the interval
                          // to avoid extrema on the boundaries
        return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
    }

    /*!
     * \brief Find all the intersections of a segment of the spline
     *        with a cubic polynomial within a specified interval.
     */
    int intersectSegment_(Scalar *sol,
                          int segIdx,
                          Scalar a, Scalar b, Scalar c, Scalar d,
                          Scalar x0 = -1e100, Scalar x1 = 1e100) const
    {
        int n = Opm::invertCubicPolynomial(sol,
                                             a_(segIdx) - a,
                                             b_(segIdx) - b,
                                             c_(segIdx) - c,
                                             d_(segIdx) - d);
        x0 = std::max(x_(segIdx), x0);
        x1 = std::max(x_(segIdx+1), x1);

        // filter the intersections outside of the specified intervall
        int k = 0;
        for (int j = 0; j < n; ++j) {
            sol[j] += x_(segIdx); // add the offset of the intervall. For details see Stoer
            if (x0 <= sol[j] && sol[j] <= x1) {
                sol[k] = sol[j];
                ++k;
            }
        }
        return k;
    }

    // find the segment index for a given x coordinate
    int segmentIdx_(Scalar x) const
    {
        // bisection
        int iLow = 0;
        int iHigh = numSamples_() - 1;

        while (iLow + 1 < iHigh) {
            int i = (iLow + iHigh) / 2;
            if (x_(i) > x)
                iHigh = i;
            else
                iLow = i;
        };
        return iLow;
    }

    /*!
     * \brief Returns x[i] - x[i - 1]
     */
    Scalar h_(int i) const
    {
        assert(x_(i) > x_(i-1)); // the sampling points must be given
                                 // in ascending order
        return x_(i) - x_(i - 1);
    }

    /*!
     * \brief Returns the y coordinate of the i-th sampling point.
     */
    Scalar x_(int i) const
    { return asImp_().x_(i); }

    /*!
     * \brief Returns the y coordinate of the i-th sampling point.
     */
    Scalar y_(int i) const
    { return asImp_().y_(i); }

    /*!
     * \brief Returns the moment (i.e. second derivative) of the
     *        spline at the i-th sampling point.
     */
    Scalar moment_(int i) const
    { return asImp_().moment_(i); }

    // returns the coefficient in front of the x^0 term. In Stoer this
    // is delta.
    Scalar a_(int i) const
    { return (moment_(i+1) - moment_(i))/(6*h_(i+1)); }

    // returns the coefficient in front of the x^2 term In Stoer this
    // is gamma.
    Scalar b_(int i) const
    { return moment_(i)/2; }

    // returns the coefficient in front of the x^1 term. In Stoer this
    // is beta.
    Scalar c_(int i) const
    { return (y_(i+1) - y_(i))/h_(i + 1) - h_(i+1)/6*(2*moment_(i) + moment_(i+1)); }

    // returns the coefficient in front of the x^0 term. In Stoer this
    // is alpha.
    Scalar d_(int i) const
    { return y_(i); }

    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples_() const
    { return asImp_().numSamples(); }
};

}

#endif
