// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 *
 * \copydoc Ewoms::Co2InjectionFlash
 */
#ifndef EWOMS_CO2_INJECTION_FLASH_HH
#define EWOMS_CO2_INJECTION_FLASH_HH

#include <opm/material/constraintsolvers/NcpFlash.hpp>

namespace Ewoms {
/*!
 * \brief Flash solver used by the CO2 injection problem.
 *
 * This class is just the NCP flash solver with the guessInitial()
 * method that is adapted to the pressure regime of the CO2 injection
 * problem.
 */
template <class Scalar, class FluidSystem>
class Co2InjectionFlash : public Opm::NcpFlash<Scalar, FluidSystem>
{
    typedef Opm::NcpFlash<Scalar, FluidSystem> ParentType;

    typedef typename FluidSystem::ParameterCache ParameterCache;

    enum { numPhases = FluidSystem::numPhases };

public:
    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState, class ComponentVector>
    static void guessInitial(FluidState &fluidState, ParameterCache &paramCache,
                             const ComponentVector &globalMolarities)
    {
        ParentType::guessInitial(fluidState, paramCache, globalMolarities);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // pressure. something close to the reservoid pressure as initial
            // guess
            fluidState.setPressure(phaseIdx, 1.0135e6);
        }

        paramCache.updateAllPressures(fluidState);
    }
};

} // namespace Ewoms

#endif // EWOMS_CO2_INJECTION_FLASH_HH
