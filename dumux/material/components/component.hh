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
 * \brief Abstract base class of a pure chemical species.
 */
#ifndef DUMUX_COMPONENT_HH
#define DUMUX_COMPONENT_HH

namespace Dumux
{

/*!
 * \brief Abstract base class of a pure chemical species.
 */
template <class Scalar, class Implementation>
class Component
{
public:
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
            Scalar pressMin, Scalar pressMax, unsigned nPress)
    {   Dune::dwarn << "No init routine defined - make shure that this is not necessary!" << std::endl; }

    /*!
     * \brief A human readable name for the compoent.
     */
    static const char *name()
    { DUNE_THROW(Dune::NotImplemented, "Component::name()"); }

    /*!
     * \brief The mass in [kg] of one mole of the component.
     */
    static Scalar molarMass()
    { DUNE_THROW(Dune::NotImplemented, "Component::molarMass()"); }

    /*!
     * \brief Returns the critical temperature of the component
     */
    static Scalar criticalTemperature()
    { DUNE_THROW(Dune::NotImplemented, "Component::criticalTemperature()"); }

    /*!
     * \brief Returns the critical pressure of the component
     */
    static Scalar criticalPressure()
    { DUNE_THROW(Dune::NotImplemented, "Component::criticalPressure()"); }

    /*!
     * \brief Returns the temperature at the component's triple point.
     */
    static Scalar tripleTemperature()
    { DUNE_THROW(Dune::NotImplemented, "Component::tripleTemperature()"); }

    /*!
     * \brief Returns the pressure at the component's triple point.
     */
    static Scalar triplePressure()
    { DUNE_THROW(Dune::NotImplemented, "Component::triplePressure()"); }

    /*!
     * \brief The vapor pressure in [N/m^2] of the component at a given
     *        temperature.
     */
    static Scalar vaporPressure(Scalar T)
    { DUNE_THROW(Dune::NotImplemented, "Component::vaporPressure()"); }

    /*!
     * \brief The density [kg/m^3] of the component at a given pressure and temperature.
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::density()"); }

    /*!
     * \brief The density [kg/m^3] of the liquid component at a given pressure and temperature.
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::density()"); }

    /*!
     * \brief Specific enthalpy [J/kg] of pure the pure component in gas.
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasEnthalpy()"); }

    /*!
     * \brief Specific enthalpy [J/kg] of pure the pure component in liquid.
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidEnthalpy()"); }

    /*!
     * \brief Specific internal energy [J/kg] of pure the pure component in gas.
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasInternalEnergy()"); }

    /*!
     * \brief Specific internal energy [J/kg] of pure the pure component in liquid.
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidInternalEnergy()"); }

    /*!
     * \brief The dynamic viscosity [Pa s] of the pure component at a given pressure and temperature.
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasViscosity()"); }

    /*!
     * \brief The dynamic liquid viscosity [N/m^3*s] of the pure component.
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidViscosity()"); }

};

} // end namepace

#endif
