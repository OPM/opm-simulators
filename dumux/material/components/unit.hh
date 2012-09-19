// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
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
 * \brief A component where all quantities are fixed at 1.0
 *
 * This component is meant as a debugging tool. Do not use it in
 * real-life applications!
 */
#ifndef DUMUX_UNIT_HH
#define DUMUX_UNIT_HH

#include "component.hh"

namespace Dumux
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
    static const char *name()
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
     * \copydoc Component::liquidIsCompressible
     */
    static Scalar vaporPressure(Scalar T)
    { return 1.0; }
    /*!
     * \copydoc Component::liquidIsCompressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \copydoc Component::gasIsCompressible
     */
    static constexpr bool gasIsCompressible()
    { return false; }

    /*!
     * \copydoc Component::gasIsIdeal
     */
    static bool gasIsIdeal()
    { return false; }

    /*!
     * \copydoc Component::liquidDensity
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidViscosity
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::gasDensity
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::gasViscosity
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { return 1.0; }


    /*!
     * \copydoc Component::gasEnthalpy
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidEnthalpy
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::gasInternalEnergy
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidInternalEnergy
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::gasThermalConductivity
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidThermalConductivity
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::gasHeatCapacity
     */
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    { return 1.0; }

    /*!
     * \copydoc Component::liquidHeatCapacity
     */
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    { return 1.0; }
};

} // end namepace

#endif
