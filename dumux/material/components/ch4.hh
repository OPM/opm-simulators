/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser
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
 * \brief Properties of methane \f$CH_4\f$.
 */
#ifndef DUMUX_CH4_HH
#define DUMUX_CH4_HH

#include <dumux/material/idealgas.hh>

#include "component.hh"

#include <cmath>

namespace Dumux
{

/*!
 * \ingroup Components
 *
 * \brief Properties of pure molecular methane \f$CH_4\f$.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class CH4 : public Component<Scalar, CH4<Scalar> >
{
    typedef Component<Scalar, CH4<Scalar> >  ParentType;
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for methane.
     */
    static const char *name()
    { return "CH4"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of molecular methane.
     */
    static Scalar molarMass()
    { return 16.043e-3;}

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of molecular methane
     */
    static Scalar criticalTemperature()
    { return 190.4; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of molecular methane
     */
    static Scalar criticalPressure()
    { return 46e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at molecular methane's triple point.
     */
    static Scalar tripleTemperature()
    { return 90.7; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at molecular methane's triple point.
     */
    static Scalar triplePressure()
    { return 0; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure molecular methane
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    { DUNE_THROW(Dune::NotImplemented, "vaporPressure for CH4"); }


    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of \f$CH_4\f$ gas at a given pressure and temperature.
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
     * \brief The pressure of gaseous \f$CH_4\f$ in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
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
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of \f$CH_4\f$ gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidDensity for CH4"); }

    /*!
     * \brief The pressure of liquid methane in \f$\mathrm{[Pa]}\f$ at a given density and
     *        temperature.
    *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    { DUNE_THROW(Dune::NotImplemented, "liquidPressure for CH4"); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure methane gas.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
     * edition, McGraw-Hill, 1987, pp 154, 657, 671
     */
    static const Scalar gasEnthalpy(Scalar T,
                                    Scalar pressure)
    {
        // method of Joback
        const Scalar cpVapA = 19.25;
        const Scalar cpVapB = 0.05213;
        const Scalar cpVapC = 1.197e-5;
        const Scalar cpVapD = -1.132e-8;

        //Scalar cp =
        //    cpVapA + T*(cpVapB + T*(cpVapC + T*cpVapD));

        // calculate: \int_0^T c_p dT
        return
            1/molarMass()* // conversion from [J/mol] to [J/kg]

            T*(cpVapA + T*
               (cpVapB/2 + T*
                (cpVapC/3 + T*
                 (cpVapD/4))));
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure liquid \f$CH_4\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for CH4"); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure methane gas.
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
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure liquid \f$CH_4\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidInternalEnergy of CH4"); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of \f$CH_4\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids,
     * 4th edition, McGraw-Hill, 1987, pp 396-397, 670
     * 5th edition, McGraw-Hill, 2001  pp 9.7-9.8 (omega and V_c taken from p. A.5)
     *
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 98.6; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.011; // accentric factor
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
     * \brief The dynamic liquid viscosity \f$\mathrm{[Pa*s]}\f$ of pure \f$CH_4\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "liquidViscosity for CH4"); }
};

} // end namepace

#endif
