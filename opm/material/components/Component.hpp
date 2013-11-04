// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Felix Bode                                        *
 *   Copyright (C) 2010 by Katherina Baber                                   *
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
 * \copydoc Opm::Component
 */
#ifndef OPM_COMPONENT_HH
#define OPM_COMPONENT_HH

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>

namespace Opm {

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
    static const bool isTabulated = false;

    /*!
     * \brief A default routine for initialization, not needed for components and must not be called.
     *
     * \param tempMin The minimum of the temperature range in \f$\mathrm{[K]}\f$
     * \param tempMax The maximum of the temperature range in \f$\mathrm{[K]}\f$
     * \param nTemp The number of entries/steps within the temperature range
     * \param pressMin The minimum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param nPress The number of entries/steps within the pressure range
     *
     * This function throws a warning when called: "No init routine defined - make sure that this is not necessary!"
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    { }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasIsCompressible()"); }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasIsIdeal()"); }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::liquidIsCompressible()"); }

    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::name()"); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg]}\f$ of the component.
     */
    static Scalar molarMass()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::molarMass()"); }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component.
     */
    static Scalar criticalTemperature()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::criticalTemperature()"); }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component.
     */
    static Scalar criticalPressure()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::criticalPressure()"); }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::tripleTemperature()"); }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure()
    { OPM_THROW(std::runtime_error, "Not implemented: Component::triplePressure()"); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature in \f$\mathrm{[K]}\f$.
     *
     * \param T temperature of the component in \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::vaporPressure()"); }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasDensity()"); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the liquid component at a given pressure in \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::liquidDensity()"); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasEnthalpy()"); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::liquidEnthalpy()"); }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasInternalEnergy()"); }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of pure the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::liquidInternalEnergy()"); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component at a given pressure in \f$\mathrm{[Pa]}\f$ and
     * temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasViscosity()"); }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::liquidViscosity()"); }

    /*!
     * \brief Thermal conductivity of the component [W/(m^2 K/m)] as a gas.
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasThermalConductivity()"); }

    /*!
     * \brief Thermal conductivity of the component [W/(m^2 K/m)] as a liquid.
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::liquidThermalConductivity()"); }

    /*!
     * \brief Specific isobaric heat capacity of the component [J/kg] as a gas.
     */
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::gasHeatCapacity()"); }

    /*!
     * \brief Specific isobaric heat capacity of the component [J/kg] as a liquid.
     */
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    { OPM_THROW(std::runtime_error, "Not implemented: Component::liquidHeatCapacity()"); }
};

} // namespace Opm

#endif
