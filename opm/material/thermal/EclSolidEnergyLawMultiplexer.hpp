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
 * \copydoc Opm::EclSolidEnergyLawMultiplexer
 */
#ifndef OPM_ECL_SOLID_ENERGY_LAW_MULTIPLEXER_HPP
#define OPM_ECL_SOLID_ENERGY_LAW_MULTIPLEXER_HPP

#include "EclSolidEnergyLawMultiplexerParams.hpp"

#include "EclHeatcrLaw.hpp"
#include "EclSpecrockLaw.hpp"
#include "NullSolidEnergyLaw.hpp"

#include <opm/material/densead/Math.hpp>

namespace Opm
{
/*!
 * \ingroup material
 *
 * \brief Provides the energy storage relation of rock
 */
template <class ScalarT,
          class FluidSystem,
          class ParamsT = EclSolidEnergyLawMultiplexerParams<ScalarT>>
class EclSolidEnergyLawMultiplexer
{
    enum { numPhases = FluidSystem::numPhases };

    typedef Opm::EclHeatcrLaw<ScalarT, FluidSystem, typename ParamsT::HeatcrLawParams> HeatcrLaw;
    typedef Opm::EclSpecrockLaw<ScalarT, typename ParamsT::SpecrockLawParams> SpecrockLaw;
    typedef Opm::NullSolidEnergyLaw<ScalarT> NullLaw;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, compute the volumetric internal energy of the rock [W/m^3].
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation solidInternalEnergy(const Params& params, const FluidState& fluidState)
    {
        switch (params.solidEnergyApproach()) {
        case Params::heatcrApproach:
            // relevant ECL keywords: HEATCR, HEATCRT and STCOND
            return HeatcrLaw::solidInternalEnergy(params.template getRealParams<Params::heatcrApproach>(), fluidState);

        case Params::specrockApproach:
            // relevant ECL keyword: SPECROCK
            return SpecrockLaw::solidInternalEnergy(params.template getRealParams<Params::specrockApproach>(), fluidState);

        case Params::nullApproach:
            // (no relevant ECL keyword)
            return NullLaw::solidInternalEnergy(0, fluidState);

        default:
            throw std::logic_error("Invalid solid energy approach: "+std::to_string(int(params.solidEnergyApproach())));
        }
    }
};
} // namespace Opm

#endif
