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
 * \copydoc Opm::Unit
 */
#ifndef OPM_UNIT_HPP
#define OPM_UNIT_HPP

#include "Component.hpp"

namespace Opm
{
/*!
 * \ingroup Components
 *
 * \brief A component where all quantities are fixed at 1.0
 *
 * This component is meant as a debugging tool. Do not use it in
 * real-life applications!
 *
 * \tparam Scalar  The type used for scalar values
 */
template <class Scalar>
class Unit : public Component<Scalar, Unit<Scalar> >
{

public:
    /*!
     * \copydoc Component::name
     */
    static const char* name()
    { return "Unit"; }

    /*!
     * \copydoc Component::molarMass
     */
    static Scalar molarMass()
    { return 1.0; }

    /*!
     * \copydoc Component::criticalTemperature
     */
    static Scalar criticalTemperature()
    { return 1.0; }

    /*!
     * \copydoc Component::criticalPressure
     */
    static Scalar criticalPressure()
    { return 1.0; }

    /*!
     * \copydoc Component::tripleTemperature
     */
    static Scalar tripleTemperature()
    { return 1.0; }

    /*!
     * \copydoc Component::triplePressure
     */
    static Scalar triplePressure()
    { return 1.0; }

    /*!
     * \copydoc Component::vaporPressure
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& /* temperature */)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidIsCompressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \copydoc Component::gasIsCompressible
     */
    static bool gasIsCompressible()
    { return false; }

    /*!
     * \copydoc Component::gasIsIdeal
     */
    static bool gasIsIdeal()
    { return false; }

    /*!
     * \copydoc Component::liquidDensity
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidViscosity
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::gasDensity
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::gasViscosity
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }


    /*!
     * \copydoc Component::gasEnthalpy
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidEnthalpy
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::gasInternalEnergy
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidInternalEnergy
     */
    template <class Evaluation>
    static Evaluation liquidInternalEnergy(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::gasThermalConductivity
     */
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidThermalConductivity
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::gasHeatCapacity
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidHeatCapacity
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& /* temperature */, const Evaluation& /* pressure */)
    { return 1.0; }
};

} // namespace Opm

#endif
