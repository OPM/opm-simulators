/*
  Copyright 2024, SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "SingleCompWellState.hpp"

namespace Opm
{

// template <typename FluidSystem, typename Indices>
// void
// CompWellPrimaryVariables<FluidSystem, Indices>::
// init()
// {
//     // the following code looks like a resize function
//     value_.resize(numWellEq, 0.);
//     evaluation_.resize(numWellEq, 0.);
// }

template <typename FluidSystem, typename Indices>
void
CompWellPrimaryVariables<FluidSystem, Indices>::
update(const SingleCompWellState<Scalar>& well_state)
{
    value_[QTotal] = well_state.get_total_surface_rate();
    for (int i = 0; i < numWellConservationEq; ++i) { // should be the same with numComponents
        value_[i + 1] = well_state.total_molar_fractions[i];
    }
    value_[Bhp] = well_state.bhp;
}

template <typename FluidSystem, typename Indices>
void
CompWellPrimaryVariables<FluidSystem, Indices>::
updateEvaluation()
{
    for (std::size_t idx = 0; idx < numWellEq; ++idx) {
        evaluation_[idx] = 0.;
        evaluation_[idx].setValue(value_[idx]);
        evaluation_[idx].setDerivative(idx + numResEq, 1.);
    }
}

template <typename FluidSystem, typename Indices>
typename CompWellPrimaryVariables<FluidSystem, Indices>::FluidStateScalar
CompWellPrimaryVariables<FluidSystem, Indices>::
toFluidStateScalar() const
{
    FluidStateScalar fluid_state;
    // will be different if more connections are involved
    const auto& pressure = value_[Bhp];
    std::array<Scalar, FluidSystem::numComponents> total_molar_fractions;
    for (int i = 0; i < FluidSystem::numComponents - 1; ++i) {
        total_molar_fractions[i] = value_[i + 1];
    }
    total_molar_fractions[FluidSystem::numComponents - 1] = 1.0 - std::accumulate(total_molar_fractions.begin(),
                                                                                 total_molar_fractions.end() - 1,
                                                                                 0.0);
    for (int i = 0; i < FluidSystem::numComponents; ++i) {
        fluid_state.setMoleFraction(i, total_molar_fractions[i]);
    }
    fluid_state.setPressure(FluidSystem::oilPhaseIdx, pressure);
    fluid_state.setPressure(FluidSystem::gasPhaseIdx, pressure);

    return fluid_state;
}


}