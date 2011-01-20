/*****************************************************************************
 *   Copyright (C) 2009-2010 by Melanie Darcis                               *
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief A simple class for the \f$CO_2\f$ fluid properties
 */
#ifndef DUMUX_SIMPLE_CO2_HH
#define DUMUX_SIMPLE_CO2_HH

#include <dumux/material/idealgas.hh>

#include "component.hh"

#include <cmath>

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief A class for the \f$CO_2\f$ fluid properties
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SimpleCO2 : public Component<Scalar, SimpleCO2<Scalar> >
{
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for the \f$CO_2\f$.
     */
    static const char *name()
    { return "CO2"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of \f$CO_2\f$.
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
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure \f$CO_2\f$
     *        at a given temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    { DUNE_THROW(Dune::NotImplemented, "vaporPressure of simple CO2"); }

    /*!
     * \brief Specific enthalpy of gaseous \f$CO_2\f$ \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    { return 571.3e3 + (temperature - 298.15)*0.85e3; }

    /*!
     * \brief Specific enthalpy of liquid \f$CO_2\f$ \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    { return (temperature - 298.15)*5e3; }

    /*!
     * \brief Specific internal energy of \f$CO_2\f$ \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            IdealGas::R*temperature; // = pressure * spec. volume for an ideal gas
    }

    /*!
     * \brief Specific internal energy of liquid \f$CO_2\f$ \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidInternalEnergy of simple CO2"); }

    /*!
     * \brief The density of \f$CO_2\f$ at a given pressure and temperature [kg/m^3].
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \brief The pressure of gaseous \f$CO_2\f$ at a given density and temperature \f$\mathrm{[Pa]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief The density of pure \f$CO_2\f$ at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidDensity of simple CO2"); }

    /*!
     * \brief The pressure of liquid \f$CO_2\f$ in \f$\mathrm{[Pa]}\f$ at a given density and
     *        temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    { DUNE_THROW(Dune::NotImplemented, "liquidPressure for simple CO2"); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of \f$CO_2\f$ at a given pressure and temperature.
     *
     *\param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
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

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure \f$CO_2\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidViscosity of simple CO2"); }
};

} // end namepace

#endif
