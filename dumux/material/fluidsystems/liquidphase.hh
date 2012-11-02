// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2011 by Bernd Flemisch                               *
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
 * \copydoc Dumux::LiquidPhase
 */
#ifndef DUMUX_LIQUID_PHASE_HH
#define DUMUX_LIQUID_PHASE_HH

namespace Dumux {

/*!
 * \ingroup Fluidsystems
 * \brief Represents the liquid phase of a single (pseudo-) component.
 */
template <class Scalar, class ComponentT>
class LiquidPhase
{
public:
    //! \copydoc GasPhase::Component
    typedef ComponentT Component;

    //! \copydoc GasPhase::name
    static const char *name()
    { return Component::name(); }

    //! \copydoc GasPhase::isLiquid
    static constexpr bool isLiquid()
    { return true; }

    //! \copydoc GasPhase::isCompressible
    static constexpr bool isCompressible()
    { return Component::liquidIsCompressible(); }

    //! \copydoc GasPhase::isIdealGas
    static constexpr bool isIdealGas()
    { return false; /* we're a liquid! */ }

    //! \copydoc GasPhase::molarMass
    static constexpr Scalar molarMass()
    {  return Component::molarMass(); }

    //! \copydoc GasPhase::criticalTemperature
    static constexpr Scalar criticalTemperature()
    {  return Component::criticalTemperature(); }

    //! \copydoc GasPhase::criticalPressure
    static constexpr Scalar criticalPressure()
    {  return Component::criticalPressure(); }

    //! \copydoc GasPhase::tripleTemperature
    static constexpr Scalar tripleTemperature()
    {  return Component::tripleTemperature(); }

    //! \copydoc GasPhase::triplePressure
    static constexpr Scalar triplePressure()
    { return Component::triplePressure(); }

    //! \copydoc GasPhase::vaporPressure
    static Scalar vaporPressure(Scalar temperature)
    {  return Component::vaporPressure(temperature); }

    //! \copydoc GasPhase::density
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::liquidDensity(temperature, pressure); }

    //! \copydoc GasPhase::pressure
    static Scalar pressure(Scalar temperature, Scalar density)
    {  return Component::liquidPressure(temperature, density); }

    //! \copydoc GasPhase::enthalpy
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {  return Component::liquidEnthalpy(temperature, pressure); }

    //! \copydoc GasPhase::internalEnergy
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    { return Component::liquidInternalEnergy(temperature, pressure); }

    //! \copydoc GasPhase::viscosity
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::liquidViscosity(temperature, pressure); }

    //! \copydoc GasPhase::thermalConductivity
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    { return Component::liquidThermalConductivity(temperature, pressure); }

    //! \copydoc GasPhase::heatCapacity
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    { return Component::liquidHeatCapacity(temperature, pressure); }
};
} // namespace

#endif
