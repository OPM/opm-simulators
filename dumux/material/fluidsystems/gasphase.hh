// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Bernd Flemisch                               *
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 *
 * \brief Represents the gas phase of a single (pseudo-) component
 */
#ifndef DUMUX_GAS_PHASE_HH
#define DUMUX_GAS_PHASE_HH

namespace Dumux {

/*!
 * \ingroup Fluidsystems
 * \brief Represents the gas phase of a single (pseudo-) component
 */
template <class Scalar, class ComponentT>
class GasPhase
{
public:
    /*!
     * \brief The type of the phase's underlying (pseudo-) component
     */
    typedef ComponentT Component;

    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    { return Component::name(); }

    /*!
     * \brief Returs whether the fluid is a liquid
     */
    static constexpr bool isLiquid()
    { return false; }

    /*!
     * \brief Returns true iff the fluid is assumed to be compressible
     */
    static constexpr bool isCompressible()
    { return Component::gasIsCompressible(); }

    /*!
     * \brief Returns true iff the fluid is assumed to be an ideal gas
     */
    static constexpr bool isIdealGas()
    { return Component::gasIsIdeal(); }

    /*!
     * \brief The mass in [kg] of one mole of the component.
     */
    static constexpr Scalar molarMass()
    {  return Component::molarMass(); }

    /*!
     * \brief Returns the critical temperature of the component
     */
    static constexpr Scalar criticalTemperature()
    {  return Component::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure of the component
     */
    static constexpr Scalar criticalPressure()
    {  return Component::criticalPressure(); }

    /*!
     * \brief Returns the temperature at the component's triple point.
     */
    static constexpr Scalar tripleTemperature()
    {  return Component::tripleTemperature(); }

    /*!
     * \brief Returns the pressure at the component's triple point.
     */
    static constexpr Scalar triplePressure()
    { return Component::triplePressure(); }

    /*!
     * \brief The vapor pressure in [N/m^2] of the component at a given
     *        temperature.
     *
     * \copydetails Doxygen::temperatureParam
     */
    static Scalar vaporPressure(Scalar temperature)
    { return Component::vaporPressure(temperature); }

    /*!
     * \brief The density [kg/m^3] of the component at a given pressure and temperature.
     *
     * \copydetails Doxygen::TpParams
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::gasDensity(temperature, pressure); }

    /*!
     * \brief The pressure [Pa] of the component at a given density and temperature.
     *
     * \param temperature The temperature of interest [K]
     * \param density The density of interest [kg / m^3]
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {  return Component::gasPressure(temperature, density); }

    /*!
     * \brief Specific enthalpy [J/kg] the pure component as a gas.
     *
     * \copydetails Doxygen::TpParams
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {  return Component::gasEnthalpy(temperature, pressure); }

    /*!
     * \brief Specific internal energy [J/kg] the pure component as a gas.
     *
     * \copydetails Doxygen::TpParams
     */
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    { return Component::gasInternalEnergy(temperature, pressure); }

    /*!
     * \brief The dynamic viscosity [Pa s] of the pure component at a given pressure and temperature.
     *
     * \copydetails Doxygen::TpParams
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::gasViscosity(temperature, pressure); }

    /*!
     * \brief Thermal conductivity of the fluid [W/(m^2 K/m)].
     *
     * \copydetails Doxygen::TpParams
     */
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    { return Component::gasThermalConductivity(temperature, pressure); }

    /*!
     * \brief Specific isobaric heat capacity of the fluid [J/kg].
     *
     * \copydetails Doxygen::TpParams
     */
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    { return Component::gasHeatCapacity(temperature, pressure); }
};
} // namespace

#endif
