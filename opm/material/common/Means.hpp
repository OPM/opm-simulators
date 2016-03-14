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
 * \brief Implements some common averages.
 *
 * i.e., arithmetic, geometric and harmonic averages.
 */
#ifndef OPM_MEANS_HH
#define OPM_MEANS_HH

#include <cmath>

namespace Opm {
/*!
 * \brief Computes the arithmetic average of two values.
 *
 * This uses the usual definition of the arithmethic mean:
 * \f[
 <a(x,y)> = (x+y)/2
\f]
 */
template <class Scalar>
inline Scalar arithmeticMean(Scalar x, Scalar y)
{ return (x+y)/2; }

/*!
 * \brief Computes the geometric average of two values.
 *
 * This uses the usual definition of the geometric mean:
 * \f[
 <a(x,y)> = \sqrt{x^2 + y^2}
\f]
 */
template <class Scalar>
inline Scalar geometricMean(Scalar x, Scalar y)
{
    if (x*y <= 0.0)
        return 0.0;

    return std::sqrt(x*y);
}

/*!
 * \brief Computes the harmonic average of two values.
 *
 * This uses the usual definition of the harmonic mean:
 * \f[
 <a(x,y)> = \frac{2}{1/x + 1/y}
\f]
 */
template <class Scalar>
inline Scalar harmonicMean(Scalar x, Scalar y)
{
    if (x*y <= 0)
        return 0.0;

    return (2*x*y)/(y + x);
}

} // namespace Ewoms

#endif // EWOMS_AVERAGE_HH
