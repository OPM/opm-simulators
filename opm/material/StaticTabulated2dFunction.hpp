// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \copydoc Opm::StaticTabulated2DFunction
 */
#ifndef OPM_STATIC_TABULATED_2D_FUNCTION_HPP
#define OPM_STATIC_TABULATED_2D_FUNCTION_HPP

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/material/Tabulated2dFunction.hpp>

#include <assert.h>

namespace Opm {
/*!
 * \copydoc Opm::Tabulated2DFunction
 *
 * This class can be used when the sampling points are calculated at
 * compile time.
 */
template <class Traits>
class StaticTabulated2DFunction
    : public Tabulated2DFunction<typename Traits::Scalar, StaticTabulated2DFunction<Traits> >
{
    typedef typename Traits::Scalar Scalar;

public:
    StaticTabulated2DFunction()
    { }

    /*!
     * \brief Returns the minimum of the X coordinate of the sampling points.
     */
    Scalar xMin() const
    { return Traits::xMin; }

    /*!
     * \brief Returns the maximum of the X coordinate of the sampling points.
     */
    Scalar xMax() const
    { return Traits::xMax; }

    /*!
     * \brief Returns the minimum of the Y coordinate of the sampling points.
     */
    Scalar yMin() const
    { return Traits::yMin; }

    /*!
     * \brief Returns the maximum of the Y coordinate of the sampling points.
     */
    Scalar yMax() const
    { return Traits::yMax; }

    /*!
     * \brief Returns the number of sampling points in X direction.
     */
    int numX() const
    { return Traits::numX; }

    /*!
     * \brief Returns the number of sampling points in Y direction.
     */
    int numY() const
    { return Traits::numY; }

    /*!
     * \brief Get the value of the sample point which is at the
     *         intersection of the \f$i\f$-th interval of the x-Axis
     *         and the \f$j\f$-th of the y-Axis.
     */
    Scalar getSamplePoint(int i, int j) const
    {
#if !defined NDEBUG
        if (i < 0 || i >= Traits::numX ||
            j < 0 || j >= Traits::numY) {
            OPM_THROW(NumericalProblem,
                       "Attempt to access element ("
                       << i << ", " << j
                       << ") on a " << Traits::name << " table of size ("
                       << Traits::numX << ", " << Traits::numY
                       << ")\n");
        };
#endif
        return Traits::vals[i][j];
    }
};
} // namespace Opm

#endif
