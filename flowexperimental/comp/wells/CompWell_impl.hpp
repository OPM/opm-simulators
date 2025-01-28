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

    {
        // flash calculation in the wellbore
        using FluidState = CompositionalFluidState<Scalar, FluidSystem>;
        FluidState fluid_state_scalar = this->primary_variables_.toFluidStateScalar();
        PTFlash<Scalar, FluidSystem>::flash_solve_scalar_(fluid_state_scalar, "ssi", 1.e-6, CompositionalConfig::EOSType::PR, 3);
        // calculating the mass within the wellbore
        constexpr Scalar R = Constants<Scalar>::R;
        typename FluidSystem::template ParameterCache<Scalar> param_cache {CompositionalConfig::EOSType::PR};
        param_cache.updatePhase(fluid_state_scalar, FluidSystem::oilPhaseIdx);
        const Scalar Z_L = (param_cache.molarVolume(FluidSystem::oilPhaseIdx) * fluid_state_scalar.pressure(FluidSystem::oilPhaseIdx) )/
                           (R * fluid_state_scalar.temperature(FluidSystem::oilPhaseIdx));
        param_cache.updatePhase(fluid_state_scalar, FluidSystem::gasPhaseIdx);
        const Scalar Z_V = (param_cache.molarVolume(FluidSystem::gasPhaseIdx) * fluid_state_scalar.pressure(FluidSystem::gasPhaseIdx) )/
                           (R * fluid_state_scalar.temperature(FluidSystem::gasPhaseIdx));

        Scalar L = fluid_state_scalar.L();
        Scalar So = Opm::max((L * Z_L / ( L * Z_L + (1 - L) * Z_V)), 0.0);
        Scalar Sg = Opm::max(1 - So, 0.0);
        Scalar sumS = So + Sg;
        So /= sumS;
        Sg /= sumS;

        fluid_state_scalar.setSaturation(FluidSystem::oilPhaseIdx, So);
        fluid_state_scalar.setSaturation(FluidSystem::gasPhaseIdx, Sg);

        fluid_state_scalar.setCompressFactor(FluidSystem::oilPhaseIdx, Z_L);
        fluid_state_scalar.setCompressFactor(FluidSystem::gasPhaseIdx, Z_V);

        fluid_state_scalar.setDensity(FluidSystem::oilPhaseIdx, FluidSystem::density(fluid_state_scalar, param_cache, FluidSystem::oilPhaseIdx));
        fluid_state_scalar.setDensity(FluidSystem::gasPhaseIdx, FluidSystem::density(fluid_state_scalar, param_cache, FluidSystem::gasPhaseIdx));
        const auto density_oil = fluid_state_scalar.density(FluidSystem::oilPhaseIdx);
        const auto density_gas = fluid_state_scalar.density(FluidSystem::gasPhaseIdx);

        std::array<Scalar, FluidSystem::numComponents> oil_mass_fractions;
        std::array<Scalar, FluidSystem::numComponents> gas_mass_fractions;
        for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
            oil_mass_fractions[compidx] = fluid_state_scalar.massFraction(FluidSystem::oilPhaseIdx, compidx);
            gas_mass_fractions[compidx] = fluid_state_scalar.massFraction(FluidSystem::gasPhaseIdx, compidx);
        }

        const Scalar wellbore_volume = this->wellbore_volume_;
        for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
            this->component_masses_[compidx] = (oil_mass_fractions[compidx] * density_oil * So +
                                          gas_mass_fractions[compidx] * density_gas * Sg) * wellbore_volume;
        }
    }
    assembleWellEq(simulator, 1.0, well_state);
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
updateSecondaryQuantities(const Simulator& simulator)
{
    updateTotalMass();
    updateSurfaceQuantities(simulator);
}


template <typename TypeTag>
void
CompWell<TypeTag>::
updateTotalMass()
{
    // flash calculation in the wellbore
    using FluidState = CompositionalFluidState<EvalWell, FluidSystem>;
    FluidState fluid_state = this->primary_variables_.toFluidState();
    PTFlash<Scalar, FluidSystem>::solve(fluid_state, "ssi", 1.e-6, CompositionalConfig::EOSType::PR, 3);
    // calculating the mass within the wellbore
    constexpr Scalar R = Constants<Scalar>::R;
    typename FluidSystem::template ParameterCache<EvalWell> param_cache {CompositionalConfig::EOSType::PR};
    param_cache.updatePhase(fluid_state, FluidSystem::oilPhaseIdx);
    const EvalWell Z_L = (param_cache.molarVolume(FluidSystem::oilPhaseIdx) * fluid_state.pressure(FluidSystem::oilPhaseIdx) )/
                       (R * fluid_state.temperature(FluidSystem::oilPhaseIdx));
    param_cache.updatePhase(fluid_state, FluidSystem::gasPhaseIdx);
    const EvalWell Z_V = (param_cache.molarVolume(FluidSystem::gasPhaseIdx) * fluid_state.pressure(FluidSystem::gasPhaseIdx) )/
                       (R * fluid_state.temperature(FluidSystem::gasPhaseIdx));

    EvalWell L = fluid_state.L();
    EvalWell So = Opm::max((L * Z_L / ( L * Z_L + (1 - L) * Z_V)), 0.0);
    EvalWell Sg = Opm::max(1 - So, 0.0);
    EvalWell sumS = So + Sg;
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

    std::array<EvalWell, FluidSystem::numComponents> oil_mass_fractions;
    std::array<EvalWell, FluidSystem::numComponents> gas_mass_fractions;
    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        oil_mass_fractions[compidx] = fluid_state.massFraction(FluidSystem::oilPhaseIdx, compidx);
        gas_mass_fractions[compidx] = fluid_state.massFraction(FluidSystem::gasPhaseIdx, compidx);
    }

    EvalWell total_mass = 0.;
    const Scalar wellbore_volume = this->wellbore_volume_;
    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        this->new_component_masses_[compidx] = (oil_mass_fractions[compidx] * density_oil * So +
                                      gas_mass_fractions[compidx] * density_gas * Sg) * wellbore_volume;
        total_mass += component_masses_[compidx];
    }
    // TODO: checking all the calculation's here
    // TODO: some properties should go to the fluid_state?
    fluid_density_ = density_oil * So + density_gas * Sg;

    // TODO: the derivative of the mass fradtions does not look correct
    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        mass_fractions_[compidx] = component_masses_[compidx] / total_mass;
    }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updateSurfaceQuantities(const Simulator& simulator)
{
    const auto& surface_cond = simulator.vanguard().eclState().getTableManager().stCond();
    std::cout << " well surface condition temperature " << surface_cond.temperature << " pressure " << surface_cond.pressure << std::endl;
    if (this->well_ecl_.isInjector()) { // we look for well stream for injection composition
        const auto& inj_composition = this->well_ecl_.getInjectionProperties().gasInjComposition();
        using FluidStateScalar = CompositionalFluidState<Scalar, FluidSystem>;
        FluidStateScalar fluid_state;
        fluid_state.setTemperature(surface_cond.temperature);
        // we can have a function to set the pressure for all the phases
        fluid_state.setPressure(FluidSystem::oilPhaseIdx, surface_cond.pressure);
        fluid_state.setPressure(FluidSystem::gasPhaseIdx, surface_cond.pressure);

        for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; ++comp_idx) {
            fluid_state.setMoleFraction(comp_idx, std::max(inj_composition[comp_idx], 1.e-10));
        }

        for (int i = 0; i < FluidSystem::numComponents; ++i) {
            fluid_state.setKvalue(i, fluid_state.wilsonK_(i));
        }

        fluid_state.setLvalue(-1.);

        PTFlash<Scalar, FluidSystem>::flash_solve_scalar_(fluid_state, "ssi", 1.e-6, CompositionalConfig::EOSType::PR, 3);

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
        this->surface_conditions_.surface_densities_[FluidSystem::oilPhaseIdx] = density_oil;
        this->surface_conditions_.surface_densities_[FluidSystem::gasPhaseIdx] = density_gas;
        this->surface_conditions_.volume_fractions_[FluidSystem::oilPhaseIdx] = So;
        this->surface_conditions_.volume_fractions_[FluidSystem::gasPhaseIdx] = Sg;
        std::cout << " oil surface density " << density_oil << " gas surface density " << density_gas
                  << " oil volume fraction " << So << " gas volume fraction " << Sg << std::endl;
        // TODO: it shows it is liquid, which is not correct
    } else { // the composition will be from the wellbore
        // here, it will use the composition from the wellbore and the pressure and temperature from the surface condition
        std::cout << " well is a producer " << std::endl;
    }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updatePrimaryVariableEvaluation()
{
    this->primary_variables_.updateEvaluation();
}

template <typename TypeTag>
void
CompWell<TypeTag>::
calculateSingleConnectionRate(const Simulator& simulator,
                              std::vector<EvalWell>& cq_s) const
{
    constexpr int con_idx = 0; // TODO: to be a function argument for multiple connection wells
    constexpr int np = 2; // TODO: this will be the number of phases
    const EvalWell& bhp = this->primary_variables_.getBhp();
    const unsigned cell_idx = this->well_cells_[0];
    const auto& int_quantities = simulator.problem().model().cachedIntensiveQuantities(cell_idx, 0);
    std::vector<EvalWell> mob(np, 0.);
    getMoblity(simulator, con_idx, mob);

    const Scalar tw = this->well_index_[0]; // only one connection

    const auto& fluid_state = int_quantities->fluidState();

    const EvalWell cell_pressure = PrimaryVariables::extendEval(fluid_state.pressure(FluidSystem::oilPhaseIdx));
    const EvalWell drawdown = cell_pressure - bhp;
    //
    // the following will go to a funciton getMobility
    if (drawdown > 0.) { // producing connection
        // for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; comp_idx++) {
        //     cq_s[comp_idx] = 0.;
        // }
        std::vector<EvalWell> cq_v(np);
        for (unsigned phase_idx = 0; phase_idx < np; ++phase_idx) {
            cq_v[phase_idx] = - mob[phase_idx] * tw * drawdown;
            for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; comp_idx++) {
                const EvalWell density = PrimaryVariables::extendEval(fluid_state.density(phase_idx));
                const EvalWell mass_fraction = PrimaryVariables::extendEval(fluid_state.massFraction(phase_idx, comp_idx));
                cq_s[comp_idx] += cq_v[phase_idx] * density * mass_fraction;
            }
        }
    } else { // injecting connection
        EvalWell total_mobility = 0.;
        for (unsigned phase_idx = 0; phase_idx < np; ++phase_idx) {
            total_mobility += mob[phase_idx];
        }
        EvalWell cq_v = - total_mobility * tw * drawdown;
        for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; comp_idx++) {
            const EvalWell density = PrimaryVariables::extendEval(fluid_state.density(FluidSystem::oilPhaseIdx));
            cq_s[comp_idx] = cq_v * fluid_density_ * mass_fractions_[comp_idx];
        }
    }
}

template <typename TypeTag>
void CompWell<TypeTag>::
getMoblity(const Simulator& simulator,
           const int connectin_idx,
           std::vector<EvalWell>& mob) const
{
    const unsigned cell_idx = this->well_cells_[connectin_idx];
    const auto& int_quants = simulator.problem().model().cachedIntensiveQuantities(cell_idx, 0);
    const auto& material_law_manager = simulator.problem().materialLawManager();

    // either use mobility of the perforation cell or calculate its own
    // based on passing the saturation table index
    const int satid = this->saturation_table_number_[connectin_idx] - 1;
    const int satid_elem = material_law_manager->satnumRegionIdx(cell_idx);

    if (satid == satid_elem) {
        for (unsigned phase_idx = 0; phase_idx < FluidSystem::numPhases; ++phase_idx) {
            mob[phase_idx] = PrimaryVariables::extendEval(int_quants->mobility(phase_idx));
        }
    }

}

template <typename TypeTag>
void
CompWell<TypeTag>::
assembleWellEq(const Simulator& simulator,
               const double dt,
               const SingleCompWellState<Scalar>& well_state)
{
    this->well_equations_.clear();

    this->updateSecondaryQuantities(simulator);

    std::vector<EvalWell> connection_rates(FluidSystem::numComponents, 0.);
    calculateSingleConnectionRate(simulator, connection_rates);

    // there will be num_comp mass balance equations for each component and one for the well control equations
    // for the mass balance equations, it will be the sum of the connection rates for each component,
    // add minus the production rate for each component, will equal to the mass change for each component

}

template <typename TypeTag>
void
CompWell<TypeTag>::
assembleSourceTerm()
{
    for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; ++comp_idx) {
        EvalWell residual = this->new_component_masses_[comp_idx] - this->component_masses_[comp_idx];

//        this->well_equations_.residual()[comp_idx][PrimaryVariables::Bhp] = 0.;
//        for (unsigned con_idx = 0; con_idx < this->number_of_connection_; ++con_idx) {
//            this->well_equations_.residual()[comp_idx][PrimaryVariables::Bhp] += this->well_equations_.residual()[comp_idx][con_idx];
//        }
//        this->well_equations_.residual()[comp_idx][PrimaryVariables::Bhp] -= this->component_masses_[comp_idx];
    }

}

} // end of namespace Opm
