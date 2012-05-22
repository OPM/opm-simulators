// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \ingroup Components
 * \brief Component for Xylene
 */
#ifndef DUMUX_XYLENE_HH
#define DUMUX_XYLENE_HH

#include <cmath>
#include <dumux/material/idealgas.hh>
#include <dumux/material/components/component.hh>


namespace Dumux
{
/*!
 * \ingroup Components
 * \brief Component for Xylene
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Xylene : public Component<Scalar, Xylene<Scalar> >
{

public:
    /*!
     * \brief A human readable name for the xylene
     */
    static const char *name()
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
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for xylene");
    }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at xylene's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for xylene");
    }

    /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of pure xylene
     *        at a given temperature according to Antoine after Betz 1997 ->  Gmehling et al 1980
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */

    static Scalar vaporPressure(Scalar temperature)
    {
        const Scalar A = 7.00909;;
        const Scalar B = 1462.266;;
        const Scalar C = 215.110;;

        Scalar T = temperature - 273.15;
        Scalar psat = 1.334*std::pow(10.0, (A - (B/(T + C))));  // in [mbar]
        psat *= 100.0;  // in [Pa] (0.001*1.E5)

        return psat;
    }

    /*!
     * \brief Specific heat cap of liquid xylene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar spHeatCapLiquidPhase(Scalar temp, Scalar pressure)
    {
        Scalar CH3,C6H5,H;
        // after Reid et al. : Missenard group contrib. method (s. example 5-8)
        // Xylene: C9H12  : 3* CH3 ; 1* C6H5 (phenyl-ring) ; -2* H (this was too much!)
        // linear interpolation between table values [J/(mol K)]

        if(temp < 298.0){                          // take care: extrapolation for Temp<273
            H = 13.4 + 1.2*(temp - 273.0)/25.0;
            CH3 = 40.0 + 1.6*(temp - 273.0)/25.0;
            C6H5 = 113.0 + 4.2*(temp - 273.0)/25.0;
        }
        else if(temp < 323.0){
            H = 14.6 + 0.9*(temp - 298.0)/25.0;
            CH3 = 41.6 + 1.9*(temp - 298.0)/25.0;
            C6H5 = 117.2 + 6.2*(temp - 298.0)/25.0;
        }
        else if(temp < 348.0){
            H = 15.5 + 1.2*(temp - 323.0)/25.0;
            CH3 = 43.5 + 2.3*(temp - 323.0)/25.0;
            C6H5 = 123.4 + 6.3*(temp - 323.0)/25.0;
        }
        else {                        // take care: extrapolation for Temp>373
            H = 16.7 + 2.1*(temp - 348.0)/25.0;          // most likely leads to underestimation
            CH3 = 45.8 + 2.5*(temp - 348.0)/25.0;
            C6H5 = 129.7 + 6.3*(temp - 348.0)/25.0;
        }

        return (C6H5 + 2*CH3 - H)/molarMass();
    }


    /*!
     * \brief Specific enthalpy of liquid xylene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(Scalar temp,
                                 Scalar pressure)
    {
        const Scalar T_ref = 273.0;  // typically T_ref=0 C, but not necessarily
        const Scalar DTemp = temp - T_ref;

        // numerical integration using 2-point Gauss */
        return 0.5*DTemp*(spHeatCapLiquidPhase(0.2113*DTemp,pressure)
                          + spHeatCapLiquidPhase(0.7887*DTemp,pressure));
    }

    /*!
     * \brief latent heat of vaporization for xylene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatVap(Scalar temp, Scalar pressure)
    {
        // Variation with Temp according to Watson relation
        temp = std::min(temp, criticalTemperature()); // regularization
        temp = std::max(temp, 0.0); // regularization

        const Scalar DH_v_b = 36195.0;       // [J/mol] Chen method (at boiling point 437.9 K
        const Scalar Tr1 = 0.668;            // Tb/T_crit
        const Scalar Tr2 = temp/criticalTemperature();
        const Scalar n = 0.375;
        const Scalar DH_vap = DH_v_b * std::pow(((1.0 - Tr2)/(1.0 - Tr1)), n);

        return (DH_vap/molarMass());          // we need [J/kg]
    }

    /*!
     * \brief Specific enthalpy of xylene vapor \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return liquidEnthalpy(temperature, pressure) + heatVap(temperature, pressure);
    }

    /*!
     * \brief
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return IdealGas<Scalar>::density(molarMass(),
                                         temperature,
                                         pressure);
    }

    static Scalar molarGasDensity(Scalar temperature, Scalar pressure)
    {
        return (gasDensity(temperature, pressure)*molarMass());
    }

    /*!
     * \brief The molar density of pure xylene at a given pressure and temperature
     * \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar molarLiquidDensity(Scalar temp, Scalar pressure)
    {
        // saturated molar volume according to Lide, CRC Handbook of
        // Thermophysical and Thermochemical Data, CRC Press, 1994
        // valid for 245 < Temp < 600
        temp = std::min(temp, 500.0); // regularization
        temp = std::max(temp, 250.0); // regularization

        const Scalar A1 = 0.25919;           // from table
        const Scalar A2 = 0.0014569;         // from table
        const Scalar expo = 1.0 + std::pow((1.0 - temp/criticalTemperature()), (2.0/7.0));
        const Scalar V = A2*std::pow(A1, expo);    // liquid molar volume [m^3/mol]

        return 1.0/V;             // molar density [mol/m^3]
    }

    /*!
     * \brief The density of pure xylene at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return molarLiquidDensity(temperature, pressure)*molarMass(); // [kg/m^3]
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
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of xylene vapor
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    static Scalar gasViscosity(Scalar temp, Scalar pressure, bool regularize=true)
    {
        temp = std::min(temp, 500.0); // regularization
        temp = std::max(temp, 250.0); // regularization

        const Scalar Tr = std::max(temp/criticalTemperature(), 1e-10);
        const Scalar Fp0 = 1.0;
        const Scalar xi = 0.004623;
        const Scalar eta_xi = Fp0*(0.807*std::pow(Tr, 0.618)
                                   - 0.357*std::exp(-0.449*Tr)
                                   + 0.34*std::exp(-4.058*Tr)
                                   + 0.018);
        Scalar r = eta_xi/xi; // [1e-6 P]
        r /= 1.0e7; // [Pa s]

        return r;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure xylene.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temp, Scalar pressure)
    {
        temp = std::min(temp, 500.0); // regularization
        temp = std::max(temp, 250.0); // regularization

        const Scalar A = -3.82;
        const Scalar B = 1027.0;
        const Scalar C = -6.38e-4;
        const Scalar D = 4.52e-7;

        Scalar r = std::exp(A + B/temp + C*temp + D*temp*temp); // in [cP]
        r *= 1.0e-3; // in [Pa s]

        return r; // [Pa s]
    }
};

} // end namespace

#endif
