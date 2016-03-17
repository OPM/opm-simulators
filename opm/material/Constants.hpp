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
 * \copydoc Opm::Constants
 */
#ifndef OPM_CONSTANTS_HPP
#define OPM_CONSTANTS_HPP

#include <cmath>

namespace Opm
{

/*!
 * \brief A central place for various physical constants occuring in
 *        some equations.
 */
template<class Scalar>
class Constants
{ public:
    /*!
     * \brief The ideal gas constant [J/(mol K)]
     */
    static const Scalar R;

    /*!
     * \brief The Avogadro constant [1/mol]
     */
    static const Scalar Na;

    /*!
     * \brief The Boltzmann constant [J/K]
     */
    static const Scalar kb;

    /*!
     * \brief Speed of light in vacuum [m/s]
     */
    static const Scalar c;

    /*!
     * \brief Newtonian constant of gravitation [m^3/(kg s^2)]
     */
    static const Scalar G;

    /*!
     * \brief Planck constant [J s]
     */
    static const Scalar h;

    /*!
     * \brief Reduced Planck constant [J s]
     */
    static const Scalar hRed;
};

template<class Scalar>
const Scalar Constants<Scalar>::R = 8.314472;
template <class Scalar>
const Scalar Constants<Scalar>::Na = 6.02214179e23;
template <class Scalar>
const Scalar Constants<Scalar>::kb = R/Na;
template <class Scalar>
const Scalar Constants<Scalar>::c = 299792458.0;
template <class Scalar>
const Scalar Constants<Scalar>::G = 6.67428e-11;
template <class Scalar>
const Scalar Constants<Scalar>::h = 6.62606896e-34;
template <class Scalar>
const Scalar Constants<Scalar>::hRed = h / (2 * M_PI);
} // namespace Opm

#endif
