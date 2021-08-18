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
 * \copydoc Opm::FluidThermalConduction
 */
#ifndef OPM_FLUID_THERMAL_CONDUCTION_LAW_HPP
#define OPM_FLUID_THERMAL_CONDUCTION_LAW_HPP

#include "FluidThermalConductionLawParams.hpp"

#include <opm/material/common/Spline.hpp>

#include <algorithm>

namespace Opm {
/*!
 * \ingroup material
 *
 * \brief Implements a thermal conduction law which just takes the conductivity of a given fluid phase.
 */
template <class FluidSystem,
          class ScalarT,
          int phaseIdx,
          class ParamsT = FluidThermalConductionLawParams<ScalarT> >
class FluidThermalConductionLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the effective thermal conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation thermalConductivity(const Params&,
                                          const FluidState& fluidState)
    {
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);
        return FluidSystem::template thermalConductivity<FluidState, Evaluation>(fluidState,
                                                                                 paramCache,
                                                                                 phaseIdx);
    }
};
} // namespace Opm

#endif
