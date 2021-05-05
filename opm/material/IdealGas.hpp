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
 * \copydoc Opm::IdealGas
 */
#ifndef OPM_IDEAL_GAS_HPP
#define OPM_IDEAL_GAS_HPP

#include <opm/material/Constants.hpp>

namespace Opm {
/*!
 * \brief Relations valid for an ideal gas.
 */
template <class Scalar>
class IdealGas
{
public:
    //! The ideal gas constant \f$\mathrm{[J/mol/K]}\f$
    static const Scalar R;

    /*!
     * \brief The density of the gas in \f$\mathrm{[kg/m^3]}\f$, depending on
     *        pressure, temperature and average molar mass of the gas.
     */
    template <class Evaluation>
    static Evaluation density(const Evaluation& avgMolarMass,
                              const Evaluation& temperature,
                              const Evaluation& pressure)
    { return pressure*avgMolarMass/(R*temperature); }

    /*!
     * \brief The pressure of the gas in \f$\mathrm{[N/m^2]}\f$, depending on
     *        the molar density and temperature.
     */
    template <class Evaluation>
    static Evaluation pressure(const Evaluation& temperature,
                               const Evaluation& rhoMolar)
    { return R*temperature*rhoMolar; }

    /*!
     * \brief The molar density of the gas \f$\mathrm{[mol/m^3]}\f$,
     *        depending on pressure and temperature.
     */
    template <class Evaluation>
    static Evaluation molarDensity(const Evaluation& temperature,
                                   const Evaluation& pressure)
    { return pressure/(R*temperature); }
};

template <class Scalar>
const Scalar IdealGas<Scalar>::R = Constants<Scalar>::R;

} // namespace Opm

#endif
