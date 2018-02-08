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
 * \copydoc Opm::Xylene
 */
#ifndef OPM_XYLENE_HPP
#define OPM_XYLENE_HPP

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/Component.hpp>
#include <opm/material/Constants.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <cmath>

namespace Opm {
/*!
 * \ingroup Components
 * \brief Component for Xylene
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Xylene : public Component<Scalar, Xylene<Scalar> >
{
    typedef Constants<Scalar> Consts;
    typedef Opm::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for the xylene
     */
    static const char* name()
    { return "xylene"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of xylene
     */
    static Scalar molarMass()
    { return 0.106; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of xylene
     */
    static Scalar criticalTemperature()
    { return 617.1; }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of xylene
     */
    static Scalar criticalPressure()
    { return 35.4e5; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at xylene's triple point.
     */
    static Scalar tripleTemperature()
    { throw std::runtime_error("Not implemented: tripleTemperature for xylene"); }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at xylene's triple point.
     */
    static Scalar triplePressure()
    { throw std::runtime_error("Not implemented: triplePressure for xylene"); }

    /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of pure xylene
     *        at a given temperature according to Antoine after Betz 1997 ->  Gmehling et al 1980
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& temperature)
    {
        const Scalar A = 7.00909;;
        const Scalar B = 1462.266;;
        const Scalar C = 215.110;;

        return 100*1.334*Opm::pow(10.0, (A - (B/(temperature - 273.15 + C))));
    }

    /*!
     * \brief Specific heat cap of liquid xylene \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Missenard group contrib. method (chap 5-7, Table 5-11, s. example 5-8)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation spHeatCapLiquidPhase(const Evaluation& temperature, const Evaluation& /*pressure*/)
    {
        Evaluation CH3,C6H5,H;
        // after Reid et al. : Missenard group contrib. method (s. example 5-8)
        // Xylene: C9H12  : 3* CH3 ; 1* C6H5 (phenyl-ring) ; -2* H (this was too much!)
        // linear interpolation between table values [J/(mol K)]

        if(temperature < 298.0){                              // take care: extrapolation for Temperature<273
            H = 13.4 + 1.2*(temperature - 273.0)/25.0;        // 13.4 + 1.2 = 14.6 = H(T=298K) i.e. interpolation of table values 273<T<298
            CH3 = 40.0 + 1.6*(temperature - 273.0)/25.0;    // 40 + 1.6 = 41.6 = CH3(T=298K)
            C6H5 = 113.0 + 4.2*(temperature - 273.0)/25.0; // 113 + 4.2 = 117.2 = C6H5(T=298K)
        }
        else if(temperature < 323.0){
            H = 14.6 + 0.9*(temperature - 298.0)/25.0;        // i.e. interpolation of table values 298<T<323
            CH3 = 41.6 + 1.9*(temperature - 298.0)/25.0;
            C6H5 = 117.2 + 6.2*(temperature - 298.0)/25.0;
        }
        else if(temperature < 348.0){
            H = 15.5 + 1.2*(temperature - 323.0)/25.0;        // i.e. interpolation of table values 323<T<348
            CH3 = 43.5 + 2.3*(temperature - 323.0)/25.0;
            C6H5 = 123.4 + 6.3*(temperature - 323.0)/25.0;
        }
        else {
            H = 16.7 + 2.1*(temperature - 348.0)/25.0;         // i.e. interpolation of table values 348<T<373
            CH3 = 45.8 + 2.5*(temperature - 348.0)/25.0;        // take care: extrapolation for Temperature>373
            C6H5 = 129.7 + 6.3*(temperature - 348.0)/25.0;        // most likely leads to underestimation
        }

        return (C6H5 + 2*CH3 - H)/molarMass();// J/(mol K) -> J/(kg K)
    }


    /*!
     * \copydoc Component::liquidEnthalpy
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature, const Evaluation& pressure)
    {
        // Gauss quadrature rule:
        // Interval: [0K; temperature (K)]
        // Gauss-Legendre-Integration with variable transformation:
        // \int_a^b f(T) dT  \approx (b-a)/2 \sum_i=1^n \alpha_i f( (b-a)/2 x_i + (a+b)/2 )
        // with: n=2, legendre -> x_i = +/- \sqrt(1/3), \apha_i=1
        // here: a=0, b=actual temperature in Kelvin
        // \leadsto h(T) = \int_0^T c_p(T) dT
        //                 \approx 0.5 T * (cp( (0.5-0.5*\sqrt(1/3)) T) + cp((0.5+0.5*\sqrt(1/3)) T))
        //                = 0.5 T * (cp(0.2113 T) + cp(0.7887 T) )

        // enthalpy may have arbitrary reference state, but the empirical/fitted heatCapacity function needs Kelvin as input
        return 0.5*temperature*(spHeatCapLiquidPhase(Evaluation(0.2113*temperature),pressure)
                                + spHeatCapLiquidPhase(Evaluation(0.7887*temperature),pressure));
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at xylene's boiling point (1 atm).
     */
    static Scalar boilingTemperature()
    { return 412.3; }

    /*!
     * \brief Latent heat of vaporization for xylene \f$\mathrm{[J/kg]}\f$.
     *
     * source : Reid et al. (fourth edition): Chen method (chap. 7-11, Delta H_v = Delta H_v (T) according to chap. 7-12)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation heatVap(Evaluation temperature,
                              const Evaluation& /*pressure*/)
    {
        temperature = Opm::min(temperature, criticalTemperature()); // regularization
        temperature = Opm::max(temperature, 0.0); // regularization

        const Scalar T_crit = criticalTemperature();
        const Scalar Tr1 = boilingTemperature()/criticalTemperature();
        const Scalar p_crit = criticalPressure();

        //        Chen method, eq. 7-11.4 (at boiling)
        const Scalar DH_v_boil = Consts::R * T_crit * Tr1
            * (3.978 * Tr1 - 3.958 + 1.555*std::log(p_crit * 1e-5 /*Pa->bar*/ ) )
            / (1.07 - Tr1); /* [J/mol] */

        /* Variation with temperature according to Watson relation eq 7-12.1*/
        const Evaluation& Tr2 = temperature/criticalTemperature();
        const Scalar n = 0.375;
        const Evaluation& DH_vap = DH_v_boil * Opm::pow(((1.0 - Tr2)/(1.0 - Tr1)), n);

        return (DH_vap/molarMass());          // we need [J/kg]
    }

    /*!
     * \copydoc Component::gasEnthalpy
     *
     * The relation used here is true on the vapor pressure curve, i.e. as long
     * as there is a liquid phase present.
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature, const Evaluation& pressure)
    {
        return liquidEnthalpy(temperature, pressure) + heatVap(temperature, pressure);
    }

    /*!
     * \copydoc Component::gasDensity
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        return IdealGas::density(Evaluation(molarMass()), temperature, pressure);
    }

    /*!
     * \brief The density \f$\mathrm{[mol/m^3]}\f$ of xylene gas at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation molarGasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        return gasDensity(temperature, pressure) / molarMass();
    }

    /*!
     * \brief The molar density of pure xylene at a given pressure and temperature
     * \f$\mathrm{[mol/m^3]}\f$.
     *
     * source : Reid et al. (fourth edition): Modified Racket technique (chap. 3-11, eq. 3-11.9)
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation molarLiquidDensity(Evaluation temperature, const Evaluation& /*pressure*/)
    {
        // saturated molar volume according to Lide, CRC Handbook of
        // Thermophysical and Thermochemical Data, CRC Press, 1994
        // valid for 245 < Temperature < 600
        temperature = Opm::min(temperature, 500.0); // regularization
        temperature = Opm::max(temperature, 250.0); // regularization

        const Scalar A1 = 0.25919;           // from table
        const Scalar A2 = 0.0014569;         // from table
        const Evaluation& expo = 1.0 + Opm::pow((1.0 - temperature/criticalTemperature()), (2.0/7.0));
        return 1.0/(A2*Opm::pow(A1, expo));
    }

    /*!
     * \copydoc Component::liquidDensity
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& temperature, const Evaluation& pressure)
    { return molarLiquidDensity(temperature, pressure)*molarMass(); }

    /*!
     * \copydoc Component::gasIsCompressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \copydoc Component::gasIsIdeal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \copydoc Component::liquidIsCompressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \copydoc Component::gasViscosity
     */
    template <class Evaluation>
    static Evaluation gasViscosity(Evaluation temperature, const Evaluation& /*pressure*/)
    {
        temperature = Opm::min(temperature, 500.0); // regularization
        temperature = Opm::max(temperature, 250.0); // regularization

        const Evaluation& Tr = Opm::max(temperature/criticalTemperature(), 1e-10);
        const Scalar Fp0 = 1.0;
        const Scalar xi = 0.004623;
        const Evaluation& eta_xi = Fp0*(0.807*Opm::pow(Tr, 0.618)
                                   - 0.357*Opm::exp(-0.449*Tr)
                                   + 0.34*Opm::exp(-4.058*Tr)
                                   + 0.018);
        return eta_xi/xi / 1e7; // [Pa s]
    }

    /*!
     * \copydoc Component::liquidViscosity
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(Evaluation temperature, const Evaluation& /*pressure*/)
    {
        temperature = Opm::min(temperature, 500.0); // regularization
        temperature = Opm::max(temperature, 250.0); // regularization

        const Scalar A = -3.82;
        const Scalar B = 1027.0;
        const Scalar C = -6.38e-4;
        const Scalar D = 4.52e-7;

        return 1e-3*Opm::exp(A
                             + B/temperature
                             + C*temperature
                             + D*temperature*temperature); // in [cP]
    }
};

} // namespace Opm

#endif
