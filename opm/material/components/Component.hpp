// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Component
 */
#ifndef OPM_COMPONENT_HPP
#define OPM_COMPONENT_HPP

#include <opm/material/common/Exceptions.hpp>

namespace Opm {

/*!
 * \ingroup Components
 * \brief Abstract base class of a pure chemical species.
 *
 * \tparam ScalarT The type used for scalar values
 * \tparam Implementation Necessary for static polymorphism
 */
template <class ScalarT, class Implementation>
class Component
{
public:
    typedef ScalarT Scalar;

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
    static void init(Scalar /* tempMin */, Scalar /* tempMax */, unsigned /* nTemp */,
                     Scalar /* pressMin */, Scalar /* pressMax */, unsigned /* nPress */)
    { }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { throw std::runtime_error("Not implemented: Component::gasIsCompressible()"); }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { throw std::runtime_error("Not implemented: Component::gasIsIdeal()"); }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { throw std::runtime_error("Not implemented: Component::liquidIsCompressible()"); }

    /*!
     * \brief A human readable name for the component.
     */
    static const char* name()
    { throw std::runtime_error("Not implemented: Component::name()"); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg]}\f$ of the component.
     */
    static Scalar molarMass()
    { throw std::runtime_error("Not implemented: Component::molarMass()"); }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component.
     */
    static Scalar criticalTemperature()
    { throw std::runtime_error("Not implemented: Component::criticalTemperature()"); }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component.
     */
    static Scalar criticalPressure()
    { throw std::runtime_error("Not implemented: Component::criticalPressure()"); }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature()
    { throw std::runtime_error("Not implemented: Component::tripleTemperature()"); }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure()
    { throw std::runtime_error("Not implemented: Component::triplePressure()"); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of the component in \f$\mathrm{[K]}\f$
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& /* temperature */)
    { throw std::runtime_error("Not implemented: Component::vaporPressure()"); }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::gasDensity()"); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the liquid component at a given pressure in \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::liquidDensity()"); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::gasEnthalpy()"); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::liquidEnthalpy()"); }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::gasInternalEnergy()"); }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of pure the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidInternalEnergy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::liquidInternalEnergy()"); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component at a given pressure in \f$\mathrm{[Pa]}\f$ and
     * temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::gasViscosity()"); }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::liquidViscosity()"); }

    /*!
     * \brief Thermal conductivity of the component [W/(m^2 K/m)] as a gas.
     */
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::gasThermalConductivity()"); }

    /*!
     * \brief Thermal conductivity of the component [W/(m^2 K/m)] as a liquid.
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::liquidThermalConductivity()"); }

    /*!
     * \brief Specific isobaric heat capacity of the component [J/kg] as a gas.
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::gasHeatCapacity()"); }

    /*!
     * \brief Specific isobaric heat capacity of the component [J/kg] as a liquid.
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { throw std::runtime_error("Not implemented: Component::liquidHeatCapacity()"); }
};

} // namespace Opm

#endif
