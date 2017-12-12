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
 * \copydoc Opm::EclHeatcrLaw
 */
#ifndef OPM_ECL_HEATCR_LAW_HPP
#define OPM_ECL_HEATCR_LAW_HPP

#include "EclHeatcrLawParams.hpp"

#include <opm/material/densead/Math.hpp>

namespace Opm
{
/*!
 * \ingroup material
 *
 * \brief Implements the volumetric interior energy relations of rock used by ECL.
 *
 * This class uses the approach defined via keywords HEATCR, HEATCRT and STCOND.
 */
template <class ScalarT,
          class FluidSystem,
          class ParamsT = EclHeatcrLawParams<ScalarT> >
class EclHeatcrLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, compute the volumetric internal energy of the rock [W/m^3].
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation solidInternalEnergy(const Params& params, const FluidState& fluidState)
    {
        const Evaluation& T = fluidState.temperature(/*phaseIdx=*/0);
        const Evaluation& deltaT = T - params.referenceTemperature();

        Scalar C0 = params.referenceRockHeatCapacity();
        Scalar C1 = params.dRockHeatCapacity_dT();

        return deltaT*(C0 + deltaT*C1 / 2.0);
    }
};
} // namespace Opm

#endif
