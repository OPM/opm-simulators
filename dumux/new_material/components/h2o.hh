/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file 
 *
 * \brief Properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_H2O_HH
#define DUMUX_H2O_HH

#include <dumux/new_material/idealgas.hh>
#include <dune/common/exceptions.hh>
#include <dumux/exceptions.hh>

#include "component.hh"

#include "iapws/common.hh"
#include "iapws/region1.hh"
#include "iapws/region2.hh"
#include "iapws/region4.hh"

#include <cmath>
#include <iostream>

namespace Dune
{
/*!
 * \brief Properties of pure water \f$H_2O\f$.
 *
 * See: 
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
template <class Scalar>
class H2O : public Component<Scalar, H2O<Scalar> >
{
    typedef Component<Scalar, H2O<Scalar> > ParentType;
    
    typedef IAPWS::Common<Scalar> Common;
    typedef IAPWS::Region1<Scalar> Region1;
    typedef IAPWS::Region2<Scalar> Region2;
    typedef IAPWS::Region4<Scalar> Region4;

    static const Scalar R = Common::R;  // specific gas constant of water
public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "H2O"; } 

    /*!
     * \brief The mass in [kg] of one mole of water.
     */
    static Scalar molarMass()
    { return Common::molarMass; } 

    /*!
     * \brief Returns the critical temperature [K] of water
     */
    static Scalar criticalTemperature()
    { return Common::criticalTemperature; } 

    /*!
     * \brief Returns the critical pressure [Pa] of water
     */
    static Scalar criticalPressure()
    { return Common::criticalPressure; } 

    /*!
     * \brief Returns the temperature [K]at water's triple point.
     */
    static Scalar tripleTemperature()
    { return Common::tripleTemperature; } 

    /*!
     * \brief Returns the pressure [Pa] at water's triple point.
     */
    static Scalar triplePressure()
    { return Common::triplePressure; } 

    /*!
     * \brief The vapor pressure in [N/m^2] of pure water
     *        at a given temperature.
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar vaporPressure(Scalar T)
    { 
        if (T > criticalTemperature())
            T = criticalTemperature();
        if (T < tripleTemperature())
            T = tripleTemperature();

        return Region4::saturationPressure(T);
    }

    /*!
     * \brief Specific enthalpy of water steam [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar gasEnthalpy(Scalar temperature, 
                                    Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Enthalpy of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }
        
        return 
            Region2::tau(temperature) *
            Region2::dgamma_dtau(temperature, pressure) *
            R*temperature; 
    }

    /*!
     * \brief Specific enthalpy of liquid water [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Enthalpy of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        return
            Region1::tau(temperature) *
            Region1::dgamma_dtau(temperature, pressure) *
            R*temperature; 
    }

    /*!
     * \brief Specific internal energy of liquid water [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Internal Energy of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        return
            R * temperature *
            ( Region1::tau(temperature)*Region1::dgamma_dtau(temperature, pressure) - 
              Region1::pi(pressure)*Region1::dgamma_dpi(temperature, pressure));
    }

    /*!
     * \brief Specific internal energy of steam and water vapor [J/kg].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
    */
    static Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Internal Energy of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        return
            R * temperature *
            ( Region2::tau(temperature)*Region2::dgamma_dtau(temperature, pressure) - 
              Region2::pi(pressure)*Region2::dgamma_dpi(temperature, pressure));
    }

    /*!
     * \brief The density of steam at a given pressure and temperature [kg/m^3].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Density of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        Scalar specificVolume = 
            Region2::pi(pressure)*
            Region2::dgamma_dpi(temperature, pressure) *
            R * temperature / pressure;
        return 1/specificVolume;
    }

    /*!
     * \brief The pressure of steam at a given density and temperature [Pa].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        Valgrind::CheckDefined(temperature);
        Valgrind::CheckDefined(density);

        // We use the newton method for this. For the initial value we
        // assume steam to be an ideal gas
        Scalar pressure = IdealGas<Scalar>::pressure(temperature, density/molarMass());
        Scalar eps = pressure*1e-7;
        
        Scalar deltaP = pressure*2;
        Valgrind::CheckDefined(pressure);
        Valgrind::CheckDefined(deltaP);
        for (int i = 0; i < 5 && std::abs(pressure*1e-9) < std::abs(deltaP); ++i) {
            Scalar f = gasDensity(temperature, pressure) - density;
            
            Scalar df_dp;
            df_dp  = gasDensity(temperature, pressure + eps);
            df_dp -= gasDensity(temperature, pressure - eps);
            df_dp /= 2*eps;
            
            deltaP = - f/df_dp;
            
            pressure += deltaP;
            Valgrind::CheckDefined(pressure);
            Valgrind::CheckDefined(deltaP);
        }
        
        return pressure;
    }

    /*!
     * \brief The density of pure water at a given pressure and temperature [kg/m^3].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Density of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }
        
        Scalar specificVolume = 
            Region1::pi(pressure)*
            Region1::dgamma_dpi(temperature, pressure) *
            R * temperature / pressure;
        return 1/specificVolume;
    }

    /*!
     * \brief The pressure of liquid water at a given density and
     *        temperature [Pa].
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        // We use the newton method for this. For the initial value we
        // assume the pressure to be 10% higher than the vapor
        // pressure
        Scalar pressure = 1.1*vaporPressure(temperature);
        Scalar eps = pressure*1e-7;
        
        Scalar deltaP = pressure*2;
        for (int i = 0; i < 5 && std::abs(pressure*1e-9) < std::abs(deltaP); ++i) {
            Scalar f = liquidDensity(temperature, pressure) - density;
            
            Scalar df_dp;
            df_dp  = liquidDensity(temperature, pressure + eps);
            df_dp -= liquidDensity(temperature, pressure - eps);
            df_dp /= 2*eps;
            
            deltaP = - f/df_dp;
            
            pressure += deltaP;
        }
        
        return pressure;
    }

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of steam.
     *
     * This method is only valid if pressure is below or at the vapour
     * pressure of water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        if (!Region2::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Viscosity of steam is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        }

        Scalar rho = gasDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    };

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of pure water.
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        if (!Region1::isValid(temperature, pressure))
        {
            DUNE_THROW(NumericalProblem,
                       "Viscosity of water is only implemented for temperatures below 623.15K and "
                       "pressures below 100MPa. (T = " << temperature << ", p=" << pressure);
        };

        Scalar rho = liquidDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    };

private:
};

} // end namepace

#endif
