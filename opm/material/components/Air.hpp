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
 * \copydoc Opm::Air
 */
#ifndef OPM_AIR_HPP
#define OPM_AIR_HPP

#include <opm/material/components/Component.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/IdealGas.hpp>

#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

namespace Opm {

/*!
 * \ingroup Components
 *
 * \brief A simple class implementing the fluid properties of air.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Air : public Component<Scalar, Air<Scalar> >
{
    typedef Opm::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { throw std::runtime_error("Not implemented: Component::liquidIsCompressible()"); }

    /*!
     * \brief A human readable name for the \f$Air\f$.
     */
    static const char* name()
    { return "Air"; }

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
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    { return IdealGas::density(Evaluation(molarMass()), temperature, pressure); }

    /*!
     * \brief The pressure of gaseous \f$AIR\f$ at a given density and temperature \f$\mathrm{[Pa]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation gasPressure(const Evaluation& temperature, Scalar density)
    { return IdealGas::pressure(temperature, density/molarMass()); }

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
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& temperature, const Evaluation& /*pressure*/)
    {
        Scalar Tc = criticalTemperature();
        Scalar Vc = 84.525138; // critical specific volume [cm^3/mol]
        Scalar omega = 0.078; // accentric factor
        Scalar M = molarMass() * 1e3; // molar mas [g/mol]
        Scalar dipole = 0.0; // dipole moment [debye]

        Scalar mu_r4 = 131.3 * dipole / std::sqrt(Vc * Tc);
        mu_r4 *= mu_r4;
        mu_r4 *= mu_r4;

        Scalar Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
        Evaluation Tstar = 1.2593 * temperature/Tc;
        Evaluation Omega_v =
            1.16145*Opm::pow(Tstar, -0.14874) +
            0.52487*Opm::exp(- 0.77320*Tstar) +
            2.16178*Opm::exp(- 2.43787*Tstar);
        return 40.7851e-7*Fc*Opm::sqrt(M*temperature)/(std::pow(Vc, 2./3)*Omega_v);
    }

    // simpler method, from old constrelAir.hh
    template <class Evaluation>
    static Evaluation simpleGasViscosity(const Evaluation& temperature, const Evaluation& /*pressure*/)
    {
        if(temperature < 273.15 || temperature > 660.) {
            throw NumericalIssue("Air: Temperature "+std::to_string(Opm::scalarValue(temperature))+"K out of range");
        }
        return 1.496e-6*Opm::pow(temperature, 1.5)/(temperature + 120);
    }

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
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature, const Evaluation& /*pressure*/)
    {
        return 1005.0*temperature;
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
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {
        return
            gasEnthalpy(temperature, pressure)
            - (IdealGas::R*temperature/molarMass()); // <- pressure times specific volume of an ideal gas
    }

    /*!
     * \brief Specific heat conductivity of steam \f$\mathrm{[W/(m K)]}\f$.
     *
     * Isobaric Properties for Nitrogen in: NIST Standard Reference
     * Database Number 69, Eds. P.J. Linstrom and W.G. Mallard
     * evaluated at p=.1 MPa, T=8°C, does not change dramatically with
     * p,T
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/)
    {
        // Isobaric Properties for Nitrogen in: NIST Standard
        // see http://webbook.nist.gov/chemistry/fluid/
        // evaluated at p=.1 MPa, T=20°C
        // Nitrogen: 0.025398
        // Oxygen: 0.026105
        // lambda_air is approximately 0.78*lambda_N2+0.22*lambda_O2
        return 0.0255535;
    }

    /*!
     * \brief Specific isobaric heat capacity \f$[J/(kg K)]\f$ of pure
     *        air.
     *
     * This methods uses the formula for "zero-pressure" heat capacity
     * that is only dependent on temperature, because the pressure
     * dependence is rather small.  This one should be accurate for a
     * pressure of 1 atm.  Values taken from NASA Contractor Report
     * 4755, Real-Gas Flow Properties for NASA Langley Research Center
     * Aerothermodynamic Facilities Complex Wind Tunnels using data
     * from Hilsenrath et al 1955, "Tables of Thermal Properties of
     * Gases"
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature OPM_UNUSED,
                                      const Evaluation& pressure OPM_UNUSED)
    {
        return 1005.0;
    }
};

} // namespace Opm

#endif
