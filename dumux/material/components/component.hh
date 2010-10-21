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
 * \defgroup Material Material Properties and Laws
 */
/*!
 * \ingroup Material
 * \defgroup Components Components
 */
/*!
 * \file
 *
 * \ingroup Components
 * \brief Abstract base class of a pure chemical species.
 */
#ifndef DUMUX_COMPONENT_HH
#define DUMUX_COMPONENT_HH

#include <dune/common/stdstreams.hh>

namespace Dumux
{

/*!
 * \ingroup Components
 * \brief Abstract base class of a pure chemical species.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam Implementation Necessary for static polymorphism
 */
template <class Scalar, class Implementation>
class Component
{
public:
    /*!
     * \brief A default routine for initialization, not needed for components and must not be called.
     *
     * \param tempMin The minimum of the temperature range in [K]
     * \param tempMax The maximum of the temperature range in [K]
     * \param nTemp The number of entries/steps within the temperature range
     * \param pressMin The minimum of the pressure range in [Pa]
     * \param pressMax The maximum of the pressure range in [Pa]
     * \param nPress The number of entries/steps within the pressure range
     *
     * This function throws a warning when called: "No init routine defined - make sure that this is not necessary!"
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
            Scalar pressMin, Scalar pressMax, unsigned nPress)
    {   Dune::dwarn << "No init routine defined - make sure that this is not necessary!" << std::endl; }

    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    { DUNE_THROW(Dune::NotImplemented, "Component::name()"); }

    /*!
     * \brief The molar mass in [kg] of the component.
     */
    static Scalar molarMass()
    { DUNE_THROW(Dune::NotImplemented, "Component::molarMass()"); }

    /*!
     * \brief Returns the critical temperature in [K] of the component.
     */
    static Scalar criticalTemperature()
    { DUNE_THROW(Dune::NotImplemented, "Component::criticalTemperature()"); }

    /*!
     * \brief Returns the critical pressure in [Pa] of the component.
     */
    static Scalar criticalPressure()
    { DUNE_THROW(Dune::NotImplemented, "Component::criticalPressure()"); }

    /*!
     * \brief Returns the temperature in [K] at the component's triple point.
     */
    static Scalar tripleTemperature()
    { DUNE_THROW(Dune::NotImplemented, "Component::tripleTemperature()"); }

    /*!
     * \brief Returns the pressure in [Pa] at the component's triple point.
     */
    static Scalar triplePressure()
    { DUNE_THROW(Dune::NotImplemented, "Component::triplePressure()"); }

    /*!
     * \brief The vapor pressure in [Pa] of the component at a given
     *        temperature in [K].
     *
     * \param T temperature of the component in [K]
     */
    static Scalar vaporPressure(Scalar T)
    { DUNE_THROW(Dune::NotImplemented, "Component::vaporPressure()"); }

    /*!
     * \brief The density in [kg/m^3] of the component at a given pressure in [Pa] and temperature in [K].
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::density()"); }

    /*!
     * \brief The density [kg/m^3] of the liquid component at a given pressure in [Pa] and temperature in [K].
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::density()"); }

    /*!
     * \brief Specific enthalpy [J/kg] of the pure component in gas.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasEnthalpy()"); }

    /*!
     * \brief Specific enthalpy [J/kg] of the pure component in liquid.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidEnthalpy()"); }

    /*!
     * \brief Specific internal energy [J/kg] of the pure component in gas.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasInternalEnergy()"); }

    /*!
     * \brief Specific internal energy [J/kg] of pure the pure component in liquid.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidInternalEnergy()"); }

    /*!
     * \brief The dynamic viscosity [Pa*s] of the pure component at a given pressure in [Pa] and temperature in [K].
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasViscosity()"); }

    /*!
     * \brief The dynamic liquid viscosity [Pa*s] of the pure component.
     *
     * \param temperature temperature of component in [K]
     * \param pressure pressure of component in [Pa]
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidViscosity()"); }

};

} // end namepace

#endif
