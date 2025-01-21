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

#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include "opm/models/immiscible/immisciblemodel.hh"
#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

namespace Opm
{

template <typename TypeTag>
CompWell<TypeTag>::
CompWell(const Well& well,
         int index_of_well,
         const std::vector<CompConnectionData<Scalar>>& well_connection_data)
  : CompWellInterface<TypeTag>(well, index_of_well, well_connection_data)
{
}

template <typename TypeTag>
void
CompWell<TypeTag>::
init() {
    Base::init();
    // primary_variables_.init();
    well_equations_.init(this->number_of_connection_);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
calculateExplicitQuantities(const Simulator& simulator,
                            const SingleCompWellState<Scalar>& well_state)
{
    updatePrimaryVariables(simulator, well_state);
    updatePrimaryVariableEvaluation();

    // flash calculation in the wellbore
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;
    FluidState fluid_state = this->primary_variables_.toFluidStateScalar();
    PTFlash<Scalar, FluidSystem>::flash_solve_scalar_(fluid_state, "ssi", 1.e-6, CompositionalConfig::EOSType::PR, 3);
    // calculating the mass within the wellbore
    constexpr Scalar R = Constants<Scalar>::R;
    typename FluidSystem::template ParameterCache<Scalar> param_cache {CompositionalConfig::EOSType::PR};
    param_cache.updatePhase(fluid_state, FluidSystem::oilPhaseIdx);
    const Scalar Z_L = (param_cache.molarVolume(FluidSystem::oilPhaseIdx) * fluid_state.pressure(FluidSystem::oilPhaseIdx) )/
                       (R * fluid_state.temperature(FluidSystem::oilPhaseIdx));
    param_cache.updatePhase(fluid_state, FluidSystem::gasPhaseIdx);
    const Scalar Z_V = (param_cache.molarVolume(FluidSystem::gasPhaseIdx) * fluid_state.pressure(FluidSystem::gasPhaseIdx) )/
                       (R * fluid_state.temperature(FluidSystem::gasPhaseIdx));

    Scalar L = fluid_state.L();
    Scalar So = Opm::max((L * Z_L / ( L * Z_L + (1 - L) * Z_V)), 0.0);
    Scalar Sg = Opm::max(1 - So, 0.0);
    Scalar sumS = So + Sg;
    So /= sumS;
    Sg /= sumS;

    fluid_state.setSaturation(FluidSystem::oilPhaseIdx, So);
    fluid_state.setSaturation(FluidSystem::gasPhaseIdx, Sg);

    fluid_state.setCompressFactor(FluidSystem::oilPhaseIdx, Z_L);
    fluid_state.setCompressFactor(FluidSystem::gasPhaseIdx, Z_V);

    fluid_state.setDensity(FluidSystem::oilPhaseIdx, FluidSystem::density(fluid_state, param_cache, FluidSystem::oilPhaseIdx));
    fluid_state.setDensity(FluidSystem::gasPhaseIdx, FluidSystem::density(fluid_state, param_cache, FluidSystem::gasPhaseIdx));
    const auto density_oil = fluid_state.density(FluidSystem::oilPhaseIdx);
    const auto density_gas = fluid_state.density(FluidSystem::gasPhaseIdx);

    std::array<Scalar, FluidSystem::numComponents> oil_mass_fractions;
    std::array<Scalar, FluidSystem::numComponents> gas_mass_fractions;
    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        oil_mass_fractions[compidx] = fluid_state.massFraction(FluidSystem::oilPhaseIdx, compidx);
        gas_mass_fractions[compidx] = fluid_state.massFraction(FluidSystem::gasPhaseIdx, compidx);
    }

    // TODO: this will be a member variable of the class
    std::array<Scalar, FluidSystem::numComponents> component_masses_;
    constexpr Scalar wellbore_volume = 1.e-5;
    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        component_masses_[compidx] = (oil_mass_fractions[compidx] * density_oil * So +
                                      gas_mass_fractions[compidx] * density_gas * Sg) * wellbore_volume;
    }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updatePrimaryVariables(const Simulator& /* simulator */,
                       const SingleCompWellState<Scalar>& well_state)
{
    this->primary_variables_.update(well_state);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updatePrimaryVariableEvaluation()
{
    this->primary_variables_.updateEvaluation();
}



} // end of namespace Opm
