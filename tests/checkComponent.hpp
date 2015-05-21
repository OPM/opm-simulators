/*
  Copyright (C) 2015 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc checkComponent
 */
#ifndef OPM_CHECK_COMPONENT_HPP
#define OPM_CHECK_COMPONENT_HPP

#include <opm/material/common/Unused.hpp>
#include <opm/material/common/ClassName.hpp>

#include <iostream>
#include <string>

/*!
 * \brief Ensures that a class which represents a chemical components adheres to the
 *        components API.
 *
 * Note that this does *not* imply that the methods are implemented or even make sense...
 */
template <class Component, class Evaluation>
void checkComponent()
{
    std::cout << "Testing component '" << Opm::className<Component>() << "'\n";

    // make sure the necessary typedefs exist
    typedef typename Component::Scalar Scalar;

    // make sure the necessary constants are exported
    OPM_UNUSED bool isTabulated = Component::isTabulated;

    // test for the gas-phase functions
    Evaluation T=0, p=0;
    while (0) {
        { OPM_UNUSED bool b = Component::gasIsCompressible(); }
        { OPM_UNUSED bool b = Component::gasIsIdeal(); }
        { OPM_UNUSED bool b = Component::liquidIsCompressible(); }
        { OPM_UNUSED std::string s = Component::name(); }
        { OPM_UNUSED Scalar M = Component::molarMass(); }
        { OPM_UNUSED Scalar Tc = Component::criticalTemperature(); }
        { OPM_UNUSED Scalar pc = Component::criticalPressure(); }
        { OPM_UNUSED Scalar Tt = Component::tripleTemperature(); }
        { OPM_UNUSED Scalar pt = Component::triplePressure(); }
        { OPM_UNUSED Evaluation pv = Component::vaporPressure(T); }
        { OPM_UNUSED Evaluation rho = Component::gasDensity(T, p); }
        { OPM_UNUSED Evaluation rho = Component::liquidDensity(T, p); }
    }
    std::cout << "----------------------------------\n";
}

#endif
