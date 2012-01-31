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
#ifndef DUMUX_MESITYLENE_HH
#define DUMUX_MESITYLENE_HH

#include <dumux/material/components/component.hh>
#include <dumux/material/constants.hh>

namespace Dumux
{
/*!
 * \ingroup Components
 * \brief mesitylene
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Mesitylene : public Component<Scalar, Mesitylene<Scalar> >
{
    typedef Dumux::Constants<Scalar> Constants;

public:
    /*!
     * \brief A human readable name for the mesitylene
     */
    static const char *name()
    { return "mesitylene"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of mesitylene
     */
    static Scalar molarMass()
    { return 0.120; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of mesitylene
     */
    static Scalar criticalTemperature()
    { return 637.3; }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of mesitylene
     */
    static Scalar criticalPressure()
    { return 31.3e5; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at mesitylene's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for mesitylene");
    };

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at mesitylene's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for mesitylene");
    };

    /*!
     * \brief The saturation vapor pressure in \f$\mathrm{[Pa]}\f$ of
     *        pure mesitylene at a given temperature according to
     *        Antoine after Betz 1997, see Gmehling et al 1980
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        const Scalar A = 7.07638;
        const Scalar B = 1571.005;
        const Scalar C = 209.728;

        const Scalar T = temperature - 273.15;

        return 100 * 1.334 * std::pow<Scalar>(10.0, (A - (B / (T + C))));
    }


    /*!
     * \brief Specific enthalpy of liquid mesitylene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar DTemp = temperature-273.0; // K -> degC
        return 0.5*DTemp*(spHeatCapLiquidPhase_(0.2113*DTemp) + spHeatCapLiquidPhase_(0.7887*DTemp));
    }

    /*!
     * \brief Specific enthalpy of mesitylene vapor \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return liquidEnthalpy(temperature,pressure) + vaporizationHeat_(temperature);
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

    /*!
     * \brief The density of pure mesitylene at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return molarLiquidDensity_(temperature)*molarMass(); // [kg/m^3]

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
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of mesitylene vapor
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure, bool regularize=true)
    {
        temperature = std::min(temperature, 500.0); // regularization
        temperature = std::max(temperature, 250.0);

        // reduced temperature
        Scalar Tr = temperature/criticalTemperature();

        Scalar Fp0 = 1.0;
        Scalar xi = 0.00474;
        Scalar eta_xi = 
            Fp0*(0.807*std::pow(Tr,0.618)
                 - 0.357*std::exp(-0.449*Tr)
                 + 0.34*std::exp(-4.058*Tr)
                 + 0.018);

        return eta_xi/xi/1e7; // [Pa s]
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure mesitylene.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        temperature = std::min(temperature, 500.0); // regularization
        temperature = std::max(temperature, 250.0);

        const Scalar A = -6.749;
        const Scalar B = 2010.0;

        return std::exp(A + B/temperature)*1e-3; // [Pa s]
    }

protected:
    /*!
     * \brief The molar density of pure mesitylene at a given pressure and temperature
     * \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar molarLiquidDensity_(Scalar temperature)
    {
        temperature = std::min(temperature, 500.0); // regularization
        temperature = std::max(temperature, 250.0);

        const Scalar Z_RA = 0.2556; // from equation
        const Scalar expo = 1.0 + std::pow(1.0 - temperature/criticalTemperature(), 2.0/7.0);
        Scalar V = Constants::R*criticalTemperature()/criticalPressure()*std::pow(Z_RA, expo); // liquid molar volume [cm^3/mol]

        return 1.0/V; // molar density [mol/m^3]
    }

    /*!
     * \brief latent heat of vaporization for mesitylene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporizationHeat_(Scalar temperature)
    {
        // regularization
        temperature = std::min(temperature, criticalTemperature());
        temperature = std::max(temperature, 250.0);

        const Scalar DH_v_b = 39086.0; // [J/mol] Chen method (at boiling point 437.9 K */
        // Variation with Temp according to Watson relation
        const Scalar Tr1 = 0.687;
        const Scalar Tr2 = temperature/criticalTemperature();
        const Scalar n = 0.375;
        const Scalar DH_vap = DH_v_b * std::pow(((1.0 - Tr2)/(1.0 - Tr1)), n);

        return (DH_vap/0.120); // we need [J/kg]
    }

    /*!
     * \brief Specific heat cap of liquid mesitylene \f$\mathrm{[J/kg]}\f$.
     *
     * \param temp temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar spHeatCapLiquidPhase_(Scalar temperature)
    {
        /* according Reid et al. : Missenard group contrib. method (s. example 5-8) */
        /* Mesitylen: C9H12  : 3* CH3 ; 1* C6H5 (phenyl-ring) ; -2* H (this was to much!) */
        /* linear interpolation between table values [J/(mol K)]*/
        Scalar H, CH3, C6H5;
        if(temperature<298.) {
            // extrapolation for Temperature<273 */
            H = 13.4+1.2*(temperature-273.0)/25.;
            CH3 = 40.0+1.6*(temperature-273.0)/25.;
            C6H5 = 113.0+4.2*(temperature-273.0)/25.;
        }
        else if((temperature>=298.0)&&(temperature<323.)){
            H = 14.6+0.9*(temperature-298.0)/25.;
            CH3 = 41.6+1.9*(temperature-298.0)/25.;
            C6H5 = 117.2+6.2*(temperature-298.0)/25.;
        }
        else if((temperature>=323.0)&&(temperature<348.)){
            H = 15.5+1.2*(temperature-323.0)/25.;
            CH3 = 43.5+2.3*(temperature-323.0)/25.;
            C6H5 = 123.4+6.3*(temperature-323.0)/25.;
        }
        else {
            assert(temperature>=348.0);
        
            /* take care: extrapolation for Temperature>373 */
            H = 16.7+2.1*(temperature-348.0)/25.;          /* leads probably to underestimates    */
            CH3 = 45.8+2.5*(temperature-348.0)/25.;
            C6H5 = 129.7+6.3*(temperature-348.0)/25.;
        }

        return (C6H5 + 3*CH3 - 2*H)/molarMass();
    }
};

} // end namespace

#endif
