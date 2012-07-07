// $Id: simpleh2o.hh 3783 2010-06-24 11:33:53Z bernd $
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
 * \brief A much simpler (and thus potentially less buggy) version of
 *        pure water.
 */
#ifndef DUMUX_SIMPLE_H2O_HH
#define DUMUX_SIMPLE_H2O_HH

#include <dumux/material/idealgas.hh>
#include <dune/common/exceptions.hh>

#include "component.hh"

#include <cmath>
#include <iostream>

namespace Dumux
{
/*!
 * \brief A much simpler (and thus potentially less buggy) version of
 *        pure water.
 */
template <class Scalar>
class SimpleH2O : public Component<Scalar, SimpleH2O<Scalar> >
{
    typedef Component<Scalar, SimpleH2O<Scalar> > ParentType;
    typedef Dumux::IdealGas<Scalar> IdealGas;

    static const double R = 461.526;  // specific gas constant of water

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
    { return 18e-3; }

    /*!
     * \brief Returns the critical temperature [K] of water
     */
    static Scalar criticalTemperature()
    { return 647.096; /* [K] */ }

    /*!
     * \brief Returns the critical pressure [Pa] of water
     */
    static Scalar criticalPressure()
    { return 22.064e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature [K]at water's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.16; /* [K] */ }

    /*!
     * \brief Returns the pressure [Pa] at water's triple point.
     */
    static Scalar triplePressure()
    { return 611.657; /* [N/m^2] */ }

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
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // water is solid: We don't take sublimation into account

        static const Scalar n[10] = {
            0.11670521452767e4,  -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5,  -0.32325550322333e7,  0.14915108613530e2,
            -0.48232657361591e4,  0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        Scalar sigma = T + n[8]/(T - n[9]);

        Scalar A = (sigma + n[0])*sigma + n[1];
        Scalar B = (n[2]*sigma + n[3])*sigma + n[4];
        Scalar C = (n[5]*sigma + n[6])*sigma + n[7];

        Scalar tmp = Scalar(2.0)*C/(std::sqrt(B*B - Scalar(4.0)*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return Scalar(1e6)*tmp;
    }

    /*!
     * \brief Specific enthalpy of water steam [J/kg].
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    { return 1976*(temperature - 293.15) + 2.45e6; }

    /*!
     * \brief Specific enthalpy of liquid water [J/kg].
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        return 4180*(temperature - 293.15);
    }

    /*!
     * \brief Specific internal energy of steam [J/kg].
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            IdealGas::R*temperature; // = pressure *spec. volume for an ideal gas
    }

    /*!
     * \brief Specific internal energy of liquid water [J/kg].
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    { return
            liquidEnthalpy(temperature, pressure) -
            pressure/liquidDensity(temperature, pressure); }

    /*!
     * \brief The density [kg/m^3] of steam at a given pressure and temperature.
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return molarMass()*IdealGas::concentration(temperature, pressure);
    }

    /*
     * \brief The pressure of steam at a given density and temperature [Pa].
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief The density of pure water at a given pressure and temperature [kg/m^3].
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1000;
    }

    /*
     * \brief The pressure of water at a given density and temperature [Pa].
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The liquid pressure is undefined for incompressible fluids");
    }

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of steam.
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure, bool regularize=true)
    {
        return 1e-05;
    };

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of pure water.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1e-03;
    };
};

} // end namepace

#endif
