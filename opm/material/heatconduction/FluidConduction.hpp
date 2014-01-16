/*
  Copyright (C) 2008-2013 by Andreas Lauser

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
 * \copydoc Opm::FluidHeatConduction
 */
#ifndef OPM_FLUID_HEAT_CONDUCTION_HPP
#define OPM_FLUID_HEAT_CONDUCTION_HPP

#include "FluidConductionParams.hpp"

#include <opm/core/utility/Spline.hpp>

#include <algorithm>

namespace Opm
{
/*!
 * \ingroup material
 *
 * \brief Implements a heat conduction law which just takes the conductivity of a given fluid phase.
 */
template <class FluidSystem,
          class ScalarT,
          int phaseIdx,
          class ParamsT = FluidHeatConductionParams<ScalarT> >
class FluidHeatConduction
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the effective heat conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     */
    template <class FluidState>
    static Scalar heatConductivity(const Params &params,
                                   const FluidState &fluidState)
    {
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);
        return FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);
    }
};
} // namespace Opm

#endif
