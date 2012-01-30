/*****************************************************************************
 *   Copyright (C) 2010 leopold stadler lstadler@wahyd.tu-berlin.de          *
 *   Copyright (C) 2011 by Benjamin Faigle, Holger Class                     *   
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
 * \brief A simple class for the \f$Air\f$ fluid properties
 */
#ifndef DUMUX_AIR_HH
#define DUMUX_AIR_HH

#include <dune/common/exceptions.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/idealgas.hh>

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief A class for the \f$AIR\f$ fluid properties
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Air : public Component<Scalar, Air<Scalar> >
{
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for the \f$Air\f$.
     */
    static const char *name()
    { return "Air"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of \f$AIR\f$.
     *
     * Taken from constrelair.hh.
     */
    static Scalar molarMass()
    { return 0.02896; /* [kg/mol] */ }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of \f$AIR\f$.
     */
    static Scalar criticalTemperature()
    { return 132.531 ; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of \f$AIR\f$.
     */
    static Scalar criticalPressure()
    { return 37.86e5; /* [Pa] */ }
    
    /*!
     * \brief The density of \f$AIR\f$ at a given pressure and temperature [kg/m^3].
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of gaseous \f$AIR\f$ at a given density and temperature \f$\mathrm{[Pa]}\f$.
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
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of \f$AIR\f$ at a given pressure and temperature.
     *
     *\param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids,
     * 4th edition, McGraw-Hill, 1987, pp 396-397, 667
     * 5th edition, McGraw-Hill, 2001, pp 9.7-9.8
     *
     * accentric factor taken from:
     * Journal of Energy Resources Technology, March 2005, Vol 127
     * Formulation for the Thermodynamic Properties
     * Georeg A. Abediyi
     * University, Mississippi State
     *
     * V_c = (R*T_c)/p_c
     *
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {

        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 84.525138; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.078; // accentric factor
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

    // simpler method, from old constrelAir.hh
    static Scalar simpleGasViscosity(Scalar temperature, Scalar pressure)
    {
        Scalar r;
        if(temperature < 273.15 || temperature > 660.)
        {
            DUNE_THROW(Dune::NotImplemented,
                "ConstrelAir: Temperature out of range at ");
        }
        r = 1.496*1.E-6*pow(temperature,1.5)/(temperature+120.);
        return (r);

    };

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$
     *        with 273.15 K as basis.
     * See:
     * W. Kays, M. Crawford, B. Weigand
     * Convective heat and mass transfer, 4th edition (2005)
     * p. 431ff
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return 1005*(temperature-273.15);
    }

    /*!
     * \brief Specific internal energy of \f$AIR\f$ \f$\mathrm{[J/kg]}\f$.
     *
     * Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     * Rearranging for internal energy yields: \f$u = h - pv\f$.
     * Exploiting the Ideal Gas assumption
     * (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure) 
            -
            IdealGas::R * temperature // = pressure * molar volume for an ideal gas
            / molarMass(); // conversion from [J/(mol K)] to [J/(kg K)]
    }

};

} // end namepace

#endif
