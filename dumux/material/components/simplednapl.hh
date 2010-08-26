// $Id$
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


#include "component.hh"


namespace Dumux
{
/*!
 * \brief A much simple component for an exemplary dense NAPL (TCE).
 */
template <class Scalar>
class SimpleDNAPL : public Component<Scalar, SimpleDNAPL<Scalar> >
{
    typedef Component<Scalar, SimpleDNAPL<Scalar> > ParentType;

public:
    /*!
     * \brief A human readable name for the TCE.
     */
    static const char *name()
    { return "TCE"; }

    /*!
     * \brief The mass in [kg] of one mole of TCE.
     */
    static Scalar molarMass()
    {
        DUNE_THROW(Dune::NotImplemented, "molarMass for TCE");
    };

    /*!
     * \brief Returns the critical temperature [K] of TCE
     */
    static Scalar criticalTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalTemperature for TCE");
    };

    /*!
     * \brief Returns the critical pressure [Pa] of TCE
     */
    static Scalar criticalPressure()
    {
        DUNE_THROW(Dune::NotImplemented, "criticalPressure for TCE");
    };

    /*!
     * \brief Returns the temperature [K]at TCE's triple point.
     */
    static Scalar tripleTemperature()
    {
        DUNE_THROW(Dune::NotImplemented, "tripleTemperature for TCE");
    };

    /*!
     * \brief Returns the pressure [Pa] at TCE's triple point.
     */
    static Scalar triplePressure()
    {
        DUNE_THROW(Dune::NotImplemented, "triplePressure for TCE");
    };

    /*!
     * \brief The vapor pressure in [N/m^2] of pure TCE
     *        at a given temperature.
     */
    static Scalar vaporPressure(Scalar T)
    {
        DUNE_THROW(Dune::NotImplemented, "vaporPressure for TCE");
    };
    /*!
     * \brief Specific enthalpy of TCE steam [J/kg].
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "gasEnthalpy for TCE");
    };

    /*!
     * \brief Specific enthalpy of liquid TCE [J/kg].
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy for TCE");
    };

    /*!
     * \brief The density of steam at a given pressure and temperature [kg/m^3].
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "gasDensity for TCE");
    };

    /*!
     * \brief The density of pure TCE at a given pressure and temperature [kg/m^3].
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
        DUNE_THROW(Dune::NotImplemented, "gasViscosity for TCE");
    };

    /*!
     * \brief The dynamic viscosity [N/m^3*s] of pure TCE.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 5.7e-4;//[Pa s]
    };
};

} // end namepace

#endif
