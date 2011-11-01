/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
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
/*!
 * \file
 *
 * \brief A gaseous phase consisting of a single component.
 */
#ifndef DUMUX_GASPHASE_HH
#define DUMUX_GASPHASE_HH

namespace Dumux
{

/*!
 * \brief gaseous phase consisting of a single component
 */
template <class Scalar, class ComponentT>
class GasPhase
{
public:
    typedef ComponentT Component;
    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    { return Component::name(); }

    /*!
     * \brief The mass in [kg] of one mole of the component.
     */
    static Scalar molarMass()
    {  return Component::molarMass(); }

    /*!
     * \brief Returns the critical temperature of the component
     */
    static Scalar criticalTemperature()
    {  return Component::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure of the component
     */
    static Scalar criticalPressure()
    {  return Component::criticalPressure(); }

    /*!
     * \brief Returns the temperature at the component's triple point.
     */
    static Scalar tripleTemperature()
    {  return Component::tripleTemperature(); }

    /*!
     * \brief Returns the pressure at the component's triple point.
     */
    static Scalar triplePressure()
    { return Component::triplePressure(); }

    /*!
     * \brief The vapor pressure in [N/m^2] of the component at a given
     *        temperature.
     */
    static Scalar vaporPressure(Scalar T)
    { return Component::vaporPressure(T); }

    /*!
     * \brief The density [kg/m^3] of the component at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::gasDensity(temperature, pressure); }

    /*!
     * \brief The pressure [Pa] of the component at a given density and temperature.
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {  return Component::gasPressure(temperature, density); }

    /*!
     * \brief Specific enthalpy [J/kg] the pure component in gas.
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    { return Component::gasEnthalpy(temperature, pressure); }

    /*!
     * \brief Specific internal energy [J/kg] the pure component in gas.
     */
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    { return Component::gasInternalEnergy(temperature, pressure); }

    /*!
     * \brief The dynamic viscosity [Pa s] of the pure component at a given pressure and temperature.
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::gasViscosity(temperature, pressure); }
};
} // namespace

#endif /* DUMUX_AIRPHASE_HH */
