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
 * \copydoc Opm::DNAPL
 */
#ifndef OPM_DNAPL_HPP
#define OPM_DNAPL_HPP

#include "Component.hpp"

#include <opm/material/IdealGas.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <opm/material/common/Unused.hpp>

namespace Opm {
/*!
 * \ingroup Components
 *
 * \brief A simple implementation of a dense non-aqueous phase liquid (DNAPL).
 *
 * The parameters are chosen to roughly correspond to those of
 * trichloroethylene (TCE) at standard conditions.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class DNAPL : public Component<Scalar, DNAPL<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the TCE.
     */
    static const char* name()
    { return "DNAPL"; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of TCE.
     */
    static Scalar molarMass()
    {
        return 131.39e-3; // [kg/mol]
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure TCE
     *        at a given temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& /*T*/)
    {
        return 3900; // [Pa] (at 20C)
    }

    /*!
     * \brief The density of steam at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        return IdealGas<Scalar>::density(Evaluation(molarMass()),
                                         temperature,
                                         pressure);
    }

    /*!
     * \brief The density of pure TCE at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 1460.0; // [kg/m^3]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure TCE.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 5.7e-4; // [Pa s]
    }

    /*!
     * \brief The enthalpy of pure TCE at a given pressure and temperature \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature, const Evaluation& /*pressure*/)
    {
        return 120.0/molarMass() * temperature; // [J/kg]
    }

    /*!
     * \brief Specific isobaric heat capacity \f$[J/(kg K)]\f$ of pure
     *        liquid TCE.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature OPM_UNUSED,
                                         const Evaluation& pressure OPM_UNUSED)
    {
        return 120.0/molarMass();
    }

    /*!
     * \brief Specific heat conductivity of liquid TCE \f$\mathrm{[W/(m K)]}\f$.
     *
     * \todo The value returned here is a guess which does not necessarily correspond to reality in any way!
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 0.3;
    }
};

} // namespace Opm

#endif
