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
 *
 * \ingroup Components
 *
 * \copydoc Ewoms::SimpleCO2
 */
#ifndef EWOMS_SIMPLE_CO2_HH
#define EWOMS_SIMPLE_CO2_HH

#include <ewoms/material/idealgas.hh>
#include <ewoms/material/components/component.hh>

#include <cmath>

namespace Ewoms {

/*!
 * \ingroup Components
 *
 * \brief A simplistic class representing the \f$CO_2\f$ fluid properties
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SimpleCO2 : public Component<Scalar, SimpleCO2<Scalar> >
{
    typedef Ewoms::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \copydoc Component::name
     */
    static const char *name()
    { return "CO2"; }

    /*!
     * \copydoc Component::molarMass
     */
    static Scalar molarMass()
    { return 44e-3; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of \f$CO_2\f$.
     */
    static Scalar criticalTemperature()
    { return 273.15 + 30.95; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of \f$CO_2\f$.
     */
    static Scalar criticalPressure()
    { return 73.8e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at the triple point of \f$CO_2\f$.
     */
    static Scalar tripleTemperature()
    { return 273.15 - 56.35; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at the triple point of \f$CO_2\f$.
     */
    static Scalar triplePressure()
    { return 5.11e5; /* [N/m^2] */ }

    /*!
     * \copydoc Component::gasIsCompressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \copydoc Component::gasIsIdeal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \copydoc Component::gasEnthalpy
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    { return 571.3e3 + (temperature - 298.15)*0.85e3; }

    /*!
     * \copydoc Component::liquidEnthalpy
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    { return (temperature - 298.15)*5e3; }

    /*!
     * \copydoc Component::gasInternalEnergy
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]
            IdealGas::R*temperature; // = pressure * spec. volume for an ideal gas
    }

    /*!
     * \copydoc Component::gasDensity
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \copydoc Component::gasViscosity
     *
     * See:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 396-397, 667
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 93.9; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.239; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]
        const Scalar dipole = 0.0; // dipole moment [debye]

        Scalar mu_r4 = 131.3 * dipole / std::sqrt(Vc * Tc);
        mu_r4 *= mu_r4;
        mu_r4 *= mu_r4;

        Scalar Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
        Scalar Tstar = 1.2593 * temperature/Tc;
        Scalar Omega_v =
            1.16145*std::pow(Tstar, -0.14874) +
            0.52487*std::exp(- 0.77320*Tstar) +
            2.16178*std::exp(- 2.43787*Tstar);
        Scalar mu = 40.785*Fc*std::sqrt(M*temperature)/(std::pow(Vc, 2./3)*Omega_v);

        // convertion from micro poise to Pa s
        return mu/1e6 / 10;
    }
};

} // end namepace

#endif
