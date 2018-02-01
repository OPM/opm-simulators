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
 * \copydoc Opm::SimpleH2O
 */
#ifndef OPM_SIMPLE_H2O_HPP
#define OPM_SIMPLE_H2O_HPP

#include "Component.hpp"

#include <opm/material/IdealGas.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <opm/material/common/Unused.hpp>

#include <cmath>

namespace Opm {

/*!
 * \ingroup Components
 *
 * \brief A simple version of  pure water.
 *
 * Compared to the water formulation of IAPWS'97, this class provides
 * a much simpler component that represents the thermodynamic
 * properties of of pure water. This implies that the likelyhood for
 * bugs in this class is reduced and the numerical performance is
 * increased. (At the cost of accuracy for the representation of the
 * physical quantities, of course.)
 *
 * \tparam Scalar The type used for representing scalar values
 */
template <class Scalar>
class SimpleH2O : public Component<Scalar, SimpleH2O<Scalar> >
{
    typedef Opm::IdealGas<Scalar> IdealGas;

    static const Scalar R;  // specific gas constant of water

public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char* name()
    { return "H2O"; }

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
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static Scalar molarMass()
    { return 18e-3; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water.
     */
    static Scalar criticalTemperature()
    { return 647.096; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static Scalar criticalPressure()
    { return 22.064e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.16; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static Scalar triplePressure()
    { return 611.657; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& T)
    {
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // water is solid: We don't take sublimation into account

        static const Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        Evaluation sigma = T + n[8]/(T - n[9]);

        Evaluation A = (sigma + n[0])*sigma + n[1];
        Evaluation B = (n[2]*sigma + n[3])*sigma + n[4];
        Evaluation C = (n[5]*sigma + n[6])*sigma + n[7];

        Evaluation tmp = 2.0*C/(Opm::sqrt(B*B - 4.0*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return 1e6*tmp;
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& /*pressure*/)
    { return 1.976e3*temperature + 40.65e3/molarMass(); }


    /*!
     * \copydoc Component::gasHeatCapacity
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature OPM_UNUSED,
                                      const Evaluation& pressure OPM_UNUSED)
    { return 1.976e3; }

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature,
                                     const Evaluation& /*pressure*/)
    { return 4180*temperature; }

    /*!
     * \copydoc Component::liquidHeatCapacity
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature OPM_UNUSED,
                                         const Evaluation& pressure OPM_UNUSED)
    { return 4.184e3; }

    /*!
     * \brief Specific internal energy of steam \f$\mathrm{[J/kg]}\f$.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]
            IdealGas::R*temperature; // = pressure *spec. volume for an ideal gas
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidInternalEnergy(const Evaluation& temperature,
                                           const Evaluation& pressure)
    {
        return
            liquidEnthalpy(temperature, pressure) -
            pressure/liquidDensity(temperature, pressure);
    }

    /*!
     * \brief Specific heat conductivity of liquid water \f$\mathrm{[W/(m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& /*temperature*/,
                                                const Evaluation& /*pressure*/)
    {
        return 0.578078; // conductivity of liquid water [W / (m K ) ] IAPWS evaluated at p=.1 MPa, T=8°C
    }

    /*!
     * \brief Specific heat conductivity of steam \f$\mathrm{[W/(m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/)
    {
        return 0.028224; // conductivity of steam [W / (m K ) ] IAPWS evaluated at p=.1 MPa, T=8°C
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of steam at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        // Assume an ideal gas
        return molarMass()*IdealGas::molarDensity(temperature, pressure);
    }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation gasPressure(const Evaluation& temperature, const Evaluation& density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief The density of pure water at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 1000;
    }

    /*!
     * \brief The pressure of water in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidPressure(const Evaluation& /*temperature*/, const Evaluation& /*density*/)
    {
        throw std::logic_error("The liquid pressure is undefined for incompressible fluids");
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& /*temperature*/,
                                   const Evaluation& /*pressure*/)
    {
        return 1e-05;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 1e-03;
    }
};

template <class Scalar>
const Scalar SimpleH2O<Scalar>::R = Opm::Constants<Scalar>::R / 18e-3;

} // namespace Opm

#endif
