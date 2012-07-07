// $Id: simplednapl.hh 3777 2010-06-24 06:46:46Z bernd $
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
#ifndef DUMUX_SIMPLE_DNAPL_HH
#define DUMUX_SIMPLE_DNAPL_HH

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
class SimpleDNAPL : public Component<Scalar, SimpleDNAPL<Scalar> >
{
    typedef Component<Scalar, SimpleDNAPL<Scalar> > ParentType;

public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "DNAPL"; }

    /*!
     * \brief The mass in [kg] of one mole of DNAPL.
     */
    static Scalar molarMass()
    {
        DUNE_THROW(Dune::NotImplemented, "molarMass for DNAPL");
        return 0;
    };

    /*!
     * \brief Returns the critical temperature [K] of DNAPL
     */
    static Scalar criticalTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalTemperature for DNAPL");
        return 0;
    };

    /*!
     * \brief Returns the critical pressure [Pa] of DNAPL
     */
    static Scalar criticalPressure()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalPressure for DNAPL");
        return 0;
    };

    /*!
     * \brief Returns the temperature [K]at DNAPL's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for DNAPL");
        return 0;
    };

    /*!
     * \brief Returns the pressure [Pa] at DNAPL's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for DNAPL");
        return 0;
    };

    /*!
     * \brief The vapor pressure in [N/m^2] of pure DNAPL
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
        DUNE_THROW(Dune::NotImplemented, "vaporPressure for DNAPL");
        return 0;
    };
    /*!
     * \brief Specific enthalpy of DNAPL steam [J/kg].
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "gasEnthalpy for DNAPL");
        return 0;
    };

    /*!
     * \brief Specific enthalpy of liquid DNAPL [J/kg].
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for DNAPL");
        return 0;
    };

    /*!
     * \brief The density of steam at a given pressure and temperature [kg/m^3].
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "gasDensity for DNAPL");
        return 0;
    };

    /*!
     * \brief The density of pure DNAPL at a given pressure and temperature [kg/m^3].
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1460.0; // [kg/m^3]
    }

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of steam.
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure, bool regularize=true)
    {
        DUNE_THROW(Dune::NotImplemented, "gasViscosity for DNAPL");
        return 0;
    };

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of pure DNAPL.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 5.7e-4;//[kg/(ms)]
    };
};

} // end namepace

#endif
