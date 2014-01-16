/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
 * \copydoc Opm::Tabulated2DFunction
 */
#ifndef OPM_TABULATED_2D_FUNCTION_HPP
#define OPM_TABULATED_2D_FUNCTION_HPP

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <assert.h>

namespace Opm {
/*!
 * \brief A generic class that represents tabulated 2 dimensional functions
 *
 * This class can be used to tabulate a two dimensional function
 * \f$f(x, y)\f$ over the range \f$[x_{min}, x_{max}] \times [y_{min},
 * y_{max}]\f$. For this, the ranges of the \f$x\f$ and \f$y\f$ axes are
 * divided into \f$m\f$ and \f$n\f$ sub-intervals and the values of
 * \f$f(x_i, y_j)\f$ need to be provided. Here, \f$x_i\f$ and
 * \f$y_j\f$ are the largest positions of the \f$i\f$-th and
 * \f$j\f$-th intervall. Between these sampling points this tabulation
 * class uses linear interpolation.
 *
 * If the class is queried for a value outside of the tabulated range,
 * a \c Opm::NumericalProblem exception is thrown.
 */
template <class Scalar, class Implementation>
class Tabulated2DFunction
{
public:
    Tabulated2DFunction()
    { }

    /*!
     * \brief Return the position on the x-axis of the i-th interval.
     */
    Scalar iToX(int i) const
    {
        assert(0 <= i && i < asImp_().numX());

        return asImp_().xMin() + i*(asImp_().xMax() - asImp_().xMin())/(asImp_().numX() - 1);
    }

    /*!
     * \brief Return the position on the y-axis of the j-th interval.
      */
    Scalar jToY(int j) const
    {
        assert(0 <= j && j < asImp_().numY());

        return asImp_().yMin() + j*(asImp_().yMax() - asImp_().yMin())/(asImp_().numY() - 1);
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
    { return (x - asImp_().xMin())/(asImp_().xMax() - asImp_().xMin())*(asImp_().numX() - 1); }

    /*!
     * \brief Return the interval index of a given position on the y-axis.
     *
     * This method returns a *floating point* number. The integer part
     * should be interpreted as interval, the decimal places are the
     * position of the y value between the j-th and the (j+1)-th
     * sample point.
     */
    Scalar yToJ(Scalar y) const
    { return (y - asImp_().yMin())/(asImp_().yMax() - asImp_().yMin())*(asImp_().numY() - 1); }

    /*!
     * \brief Returns true iff a coordinate lies in the tabulated range
     */
    bool applies(Scalar x, Scalar y) const
    { return asImp_().xMin() <= x && x <= asImp_().xMax() && asImp_().yMin() <= y && y <= asImp_().yMax(); }

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
                       << asImp_().xMin() << " to " << asImp_().xMax() << " times "
                       << asImp_().yMin() << " to " << asImp_().yMax());
        };
#endif

        Scalar alpha = xToI(x);
        Scalar beta = yToJ(y);

        int i = std::max(0, std::min(asImp_().numX(), static_cast<int>(alpha)));
        int j = std::max(0, std::min(asImp_().numY(), static_cast<int>(beta)));

        alpha -= i;
        beta -= j;

        // bi-linear interpolation
        Scalar s1 = asImp_().getSamplePoint(i, j)*(1.0 - alpha) + asImp_().getSamplePoint(i + 1, j)*alpha;
        Scalar s2 = asImp_().getSamplePoint(i, j + 1)*(1.0 - alpha) + asImp_().getSamplePoint(i + 1, j + 1)*alpha;
        return s1*(1.0 - beta) + s2*beta;
    }

private:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
} // namespace Opm

#endif
