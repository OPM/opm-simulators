/*
  Copyright 2024, 2026, SINTEF Digital

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

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/input/eclipse/EclipseState/Tables/StandardCond.hpp>

#include <stdexcept>

namespace Opm {

template <typename TypeTag>
CompWell<TypeTag>::
CompWell(const Well& well,
         int index_of_well,
         const std::vector<CompConnectionData>& well_connection_data)
  : CompWellInterface<TypeTag>(well, index_of_well, well_connection_data)
{
}

template <typename TypeTag>
void
CompWell<TypeTag>::
init()
{
    Base::init();
    well_equations_.init(this->number_of_connection_, this->well_cells_);
}

template <typename TypeTag>
template <typename T>
T
CompWell<TypeTag>::
waterDensity_(const T& pressure, const Scalar temperature)
{
    // the water PVT only reads pressure and temperature; a fluid state and
    // parameter cache are built to satisfy the fluid system's interface
    FluidState<T> fluid_state;
    fluid_state.setPressure(FluidSystem::waterPhaseIdx, pressure);
    fluid_state.setTemperature(temperature);
    typename FluidSystem::template ParameterCache<T> param_cache
        {CompositionalConfig::EOSType::PR};
    return FluidSystem::density(fluid_state, param_cache, FluidSystem::waterPhaseIdx);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
calculateExplicitQuantities(const Simulator& simulator,
                            const SingleWellState& well_state)
{
    updatePrimaryVariables(simulator, well_state);
    {
        // flash calculation in the wellbore to obtain the explicit
        // component masses
        auto fluid_state_scalar = this->primary_variables_.template toFluidState<Scalar>();

        flashFluidState_(fluid_state_scalar);

        this->component_masses_ = wellboreComponentMasses(fluid_state_scalar, this->wellbore_volume_);

        if constexpr (FluidSystem::waterEnabled) {
            // water occupies its volume fraction of the wellbore; the
            // hydrocarbon masses scale with the remainder
            const Scalar wfrac = std::clamp(well_state.wellbore_water_volume_fraction,
                                            Scalar{0.}, Scalar{1.});
            for (auto& mass : this->component_masses_) {
                mass *= (1. - wfrac);
            }
            const Scalar rho_w = waterDensity_(fluid_state_scalar.pressure(FluidSystem::waterPhaseIdx),
                                               fluid_state_scalar.temperature(0));
            this->water_mass_ = wfrac * rho_w * this->wellbore_volume_;
        }
    }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updatePrimaryVariables(const Simulator& /* simulator */,
                       const SingleWellState& well_state)
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
    auto fluid_state = this->primary_variables_.template toFluidState<EvalWell>();

    flashFluidState_(fluid_state);

    this->new_component_masses_ = wellboreComponentMasses(fluid_state, this->wellbore_volume_);

    const auto& so = fluid_state.saturation(FluidSystem::oilPhaseIdx);
    const auto& sg = fluid_state.saturation(FluidSystem::gasPhaseIdx);
    const auto& density_oil = fluid_state.density(FluidSystem::oilPhaseIdx);
    const auto& density_gas = fluid_state.density(FluidSystem::gasPhaseIdx);
    // TODO: some properties should go to the fluid_state?
    fluid_density_ = density_oil * so + density_gas * sg;

    if constexpr (FluidSystem::waterEnabled) {
        // water occupies its volume fraction of the wellbore, the flashed
        // hydrocarbon system the remainder
        const EvalWell wfrac = this->primary_variables_.getWaterVolumeFraction();
        for (auto& mass : this->new_component_masses_) {
            mass *= (1. - wfrac);
        }
        const EvalWell rho_w = waterDensity_(fluid_state.pressure(FluidSystem::waterPhaseIdx),
                                             getValue(fluid_state.temperature(0)));
        this->new_water_mass_ = wfrac * rho_w * this->wellbore_volume_;
        fluid_density_ = (1. - wfrac) * fluid_density_ + wfrac * rho_w;
    }

    EvalWell total_mass = this->new_water_mass_;
    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        total_mass += this->new_component_masses_[compidx];
    }

    // The AD derivatives of the component masses and mass fractions with respect
    // to the wellbore primary variables (pressure and composition), including the
    // dependence that flows through the flash, are checked against finite
    // differences in tests/test_compwell_jacobian.cpp.
    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        mass_fractions_[compidx] = this->new_component_masses_[compidx] / total_mass;
    }
    if constexpr (FluidSystem::waterEnabled) {
        this->water_mass_fraction_ = this->new_water_mass_ / total_mass;
    }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updateSurfaceQuantities(const Simulator& simulator)
{
    const auto& surface_cond = simulator.vanguard().eclState().getTableManager().stCond();
    if (this->well_ecl_.isInjector()) { // we look for well stream for injection composition
        const auto& inj_composition = this->well_ecl_.getInjectionProperties().gasInjComposition();
        FluidState<Scalar> fluid_state;
        for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; ++comp_idx) {
            fluid_state.setMoleFraction(comp_idx, std::max(inj_composition[comp_idx], 1.e-10));
        }
        // the injection stream carries no water
        updateSurfaceCondition_(surface_cond, fluid_state, Scalar{0.});
    } else { // the composition will be from the wellbore
        // here, it will use the composition from the wellbore and the pressure and temperature from the surface condition
        auto fluid_state = this->primary_variables_.template toFluidState<EvalWell>();
        updateSurfaceCondition_(surface_cond, fluid_state, this->water_mass_fraction_);
     }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
calculateSingleConnectionRate(const Simulator& simulator,
                              std::vector<EvalWell>& con_rates) const
{
    constexpr int con_idx = 0; // TODO: to be a function argument for multiple connection wells
    // The wellbore carries the two EOS phases only; a water phase in the
    // reservoir stays in the reservoir (the model can neither inject nor
    // produce it), so the rate loops below run over the miscible phases while
    // the mobility vector is sized for every phase getMobility() fills.
    constexpr int np = FluidSystem::numMisciblePhases;
    const EvalWell& bhp = this->primary_variables_.getBhp();
    const unsigned cell_idx = this->well_cells_[0];
    const auto& int_quantities = simulator.problem().model().cachedIntensiveQuantities(cell_idx, 0);
    assert(int_quantities);
    std::vector<EvalWell> mob(FluidSystem::numPhases, 0.);
    getMobility(simulator, con_idx, mob);

    const Scalar tw = this->well_index_[0]; // only one connection

    const auto& fluid_state = int_quantities->fluidState();

    const EvalWell cell_pressure = PrimaryVariables::extendEval(fluid_state.pressure(FluidSystem::oilPhaseIdx));
    const EvalWell drawdown = cell_pressure - bhp;

    if (drawdown > 0.) { // producing connection
        std::vector<EvalWell> cq_v(np);
        for (unsigned phase_idx = 0; phase_idx < np; ++phase_idx) {
            cq_v[phase_idx] = - mob[phase_idx] * tw * drawdown;
            for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; comp_idx++) {
                const EvalWell density = PrimaryVariables::extendEval(fluid_state.density(phase_idx));
                const EvalWell mass_fraction = PrimaryVariables::extendEval(fluid_state.massFraction(phase_idx, comp_idx));
                con_rates[comp_idx] += cq_v[phase_idx] * density * mass_fraction;
            }
        }
        if constexpr (FluidSystem::waterEnabled) {
            // the water phase is pure water; its mass rate fills the extra slot
            const EvalWell cq_w = - mob[FluidSystem::waterPhaseIdx] * tw * drawdown;
            const EvalWell density = PrimaryVariables::extendEval(fluid_state.density(FluidSystem::waterPhaseIdx));
            con_rates[FluidSystem::numComponents] += cq_w * density;
        }
    } else { // injecting connection
        // the injected fluid displaces every mobile phase in the cell
        EvalWell total_mobility = 0.;
        for (unsigned phase_idx = 0; phase_idx < FluidSystem::numPhases; ++phase_idx) {
            total_mobility += mob[phase_idx];
        }
        EvalWell cq_v = - total_mobility * tw * drawdown;
        for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; comp_idx++) {
            con_rates[comp_idx] = cq_v * fluid_density_ * mass_fractions_[comp_idx];
        }
        if constexpr (FluidSystem::waterEnabled) {
            con_rates[FluidSystem::numComponents] = cq_v * fluid_density_ * water_mass_fraction_;
        }
    }
}

template <typename TypeTag>
void CompWell<TypeTag>::
getMobility(const Simulator& simulator,
            const int connectin_idx,
            std::vector<EvalWell>& mob) const
{
    const unsigned cell_idx = this->well_cells_[connectin_idx];
    const auto& int_quants = simulator.problem().model().cachedIntensiveQuantities(cell_idx, 0);
    assert(int_quants);
    const auto& material_law_manager = simulator.problem().materialLawManager();

    // either use mobility of the perforation cell or calculate its own
    // based on passing the saturation table index
    const int satid = this->saturation_table_number_[connectin_idx] - 1;
    const int satid_elem = material_law_manager->satnumRegionIdx(cell_idx);

    if (satid == satid_elem) {
        for (unsigned phase_idx = 0; phase_idx < FluidSystem::numPhases; ++phase_idx) {
            mob[phase_idx] = PrimaryVariables::extendEval(int_quants->mobility(phase_idx));
        }
    } else {
        // TODO: not sure how to handle this at the moment, throw for now
        OPM_THROW(std::logic_error,
                  "CompWell::getMobility: a connection saturation table differing from the "
                  "cell saturation region is not supported yet");
    }

}

template <typename TypeTag>
void
CompWell<TypeTag>::
assembleWellEq(const Simulator& simulator,
               const SingleWellState& well_state,
               const double dt)
{
    this->well_equations_.clear();

    this->updateSecondaryQuantities(simulator);

    assembleSourceTerm(dt);

    // one equation per hydrocarbon component plus one for water when enabled;
    // the water slot in the reservoir rate vector has the same position
    // (conti0EqIdx + numComponents) as the water conservation row here
    std::vector<EvalWell> connection_rates(PrimaryVariables::numWellConservationEq, 0.);
    calculateSingleConnectionRate(simulator, connection_rates);
    // only one perforation for now
    auto& con_rates = this->connectionRates_[0];
    for (unsigned comp_idx = 0; comp_idx < PrimaryVariables::numWellConservationEq; ++comp_idx) {
        con_rates[comp_idx] = PrimaryVariables::restrictEval(connection_rates[comp_idx]);
    }

    // here we use perf index, need to check how the things are done in the StandardWellAssemble
    // assemble the well equations related to the production/injection mass rates for each component
    for (unsigned comp_idx = 0; comp_idx < PrimaryVariables::numWellConservationEq; ++comp_idx) {
        // the signs need to be checked
        this->well_equations_.residual()[0][comp_idx] += connection_rates[comp_idx].value();
        for (unsigned pvIdx = 0; pvIdx < PrimaryVariables::numWellEq; ++pvIdx) {
            // C, needs the cell_idx
            this->well_equations_.C()[0][0][pvIdx][comp_idx] -= connection_rates[comp_idx].derivative(pvIdx + PrimaryVariables::numResEq);
            this->well_equations_.D()[0][0][comp_idx][pvIdx] += connection_rates[comp_idx].derivative(pvIdx + PrimaryVariables::numResEq);
        }

        for (unsigned pvIdx = 0; pvIdx < PrimaryVariables::numResEq; ++pvIdx) {
            this->well_equations_.B()[0][0][comp_idx][pvIdx] += connection_rates[comp_idx].derivative(pvIdx);
        }
    }

    const auto& summary_state = simulator.vanguard().summaryState();
    assembleControlEq(well_state, summary_state);

    this->well_equations_.invert();
    // there will be num_comp mass balance equations for each component and one for the well control equations
    // for the mass balance equations, it will be the sum of the connection rates for each component,
    // add minus the production rate for each component, will equal to the mass change for each component

}

template <typename TypeTag>
void
CompWell<TypeTag>::
assembleControlEq(const SingleWellState& well_state,
                  const SummaryState& summary_state)
{
    EvalWell control_eq;
    if (this->well_ecl_.isProducer()) {
        const auto prod_controls = this->well_ecl_.productionControls(summary_state);
        assembleControlEqProd(well_state, prod_controls, control_eq);
    } else {
        const auto inj_controls = this->well_ecl_.injectionControls(summary_state);
        assembleControlEqInj(well_state, inj_controls, control_eq);
    }

    this->well_equations_.residual()[0][PrimaryVariables::Bhp] = control_eq.value();
    for (unsigned pvIdx = 0; pvIdx < PrimaryVariables::numWellEq; ++pvIdx) {
        this->well_equations_.D()[0][0][PrimaryVariables::Bhp][pvIdx] = control_eq.derivative(pvIdx + PrimaryVariables::numResEq);
    }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
assembleControlEqProd(const SingleWellState& well_state,
                      const Well::ProductionControls& prod_controls,
                      EvalWell& control_eq) const
{
    // TODO: we only need to pass in the current control?
    const auto current = well_state.production_cmode;

    const auto& surface_cond = this->surface_conditions_;

    switch (current) {
    case WellProducerCMode::BHP : {
        const Scalar bhp_limit = prod_controls.bhp_limit;
        control_eq = this->primary_variables_.getBhp() - bhp_limit;
        break;
    }
    case WellProducerCMode::ORAT : {
        const Scalar rate_target = prod_controls.oil_rate;
        const EvalWell& total_rate = this->primary_variables_.getTotalRate();
        const EvalWell oil_rate = total_rate * surface_cond.volume_fractions_[FluidSystem::oilPhaseIdx];
        control_eq = oil_rate + rate_target;
        break;
    }
    case WellProducerCMode::GRAT : {
        const Scalar rate_target = prod_controls.gas_rate;
        const EvalWell& total_rate = this->primary_variables_.getTotalRate();
        const EvalWell gas_rate = total_rate * surface_cond.volume_fractions_[FluidSystem::gasPhaseIdx];
        control_eq = gas_rate + rate_target;
        break;
    }
    case WellProducerCMode::WRAT : {
        if constexpr (FluidSystem::waterEnabled) {
            const Scalar rate_target = prod_controls.water_rate;
            const EvalWell& total_rate = this->primary_variables_.getTotalRate();
            const EvalWell water_rate = total_rate * surface_cond.volume_fractions_[FluidSystem::waterPhaseIdx];
            control_eq = water_rate + rate_target;
            break;
        } else {
            OPM_THROW(std::logic_error, "WRAT control requires an active water phase");
        }
    }
    default:
        OPM_THROW(std::logic_error, "only handles BHP, ORAT, GRAT and WRAT control for producers for now");
    }
}

template <typename TypeTag>
void
CompWell<TypeTag>::
assembleControlEqInj(const SingleWellState& well_state,
                      const Well::InjectionControls& inj_controls,
                      EvalWell& control_eq) const
{
    // TODO: we only need to pass in the current control?
    const auto current = well_state.injection_cmode;

    switch (current) {
    case WellInjectorCMode::BHP : {
        const Scalar bhp_limit = inj_controls.bhp_limit;
        control_eq = this->primary_variables_.getBhp() - bhp_limit;
        break;
    }
    case WellInjectorCMode::RATE : {
        const Scalar rate_target = inj_controls.surface_rate;
        const EvalWell& injection_rate = this->primary_variables_.getTotalRate();
        control_eq = injection_rate - rate_target;
        break;
    }
    default:
        OPM_THROW(std::logic_error, "only handles BHP and RATE control for injectors for now");
    }
}


template <typename TypeTag>
void
CompWell<TypeTag>::
assembleSourceTerm(const Scalar dt)
{
    // calculating the injection mass rate for each component
    const EvalWell total_surface_rate = this->primary_variables_.getTotalRate();
    const EvalWell density = this->surface_conditions_.density();
    const EvalWell total_mass_rate = total_surface_rate * density;
    std::array<EvalWell, FluidSystem::numComponents> component_mass_rates;
    for (unsigned  comp_idx = 0; comp_idx < FluidSystem::numComponents; ++comp_idx) {
        component_mass_rates[comp_idx] = total_mass_rate * this->surface_conditions_.massFraction(comp_idx);
    }

    for (unsigned comp_idx = 0; comp_idx < FluidSystem::numComponents; ++comp_idx) {
        const EvalWell residual = (this->new_component_masses_[comp_idx] - this->component_masses_[comp_idx]) / dt - component_mass_rates[comp_idx];
        // let us put it in the well equation
        for (int pvIdx = 0; pvIdx < PrimaryVariables::numWellEq; ++pvIdx) {
            this->well_equations_.D()[0][0][comp_idx][pvIdx] += residual.derivative(pvIdx + PrimaryVariables::numResEq);
        }
        this->well_equations_.residual()[0][comp_idx] += residual.value();
    }

    if constexpr (FluidSystem::waterEnabled) {
        // water mass balance of the wellbore, in the row after the components
        const EvalWell water_mass_rate = total_mass_rate * this->surface_conditions_.waterMassFraction();
        const EvalWell residual = (this->new_water_mass_ - this->water_mass_) / dt - water_mass_rate;
        constexpr int water_row = FluidSystem::numComponents;
        for (int pvIdx = 0; pvIdx < PrimaryVariables::numWellEq; ++pvIdx) {
            this->well_equations_.D()[0][0][water_row][pvIdx] += residual.derivative(pvIdx + PrimaryVariables::numResEq);
        }
        this->well_equations_.residual()[0][water_row] += residual.value();
    }
}

template <typename TypeTag>
bool
CompWell<TypeTag>::
iterateWellEq(const Simulator& simulator,
              const Scalar dt,
              SingleWellState& well_state)
{
    constexpr int max_iter = 200;

    int it = 0;
    bool converged = false;

    do {
        updateWellControl(simulator.vanguard().summaryState(), well_state);

        assembleWellEq(simulator, well_state, dt);

        // get convergence
        converged = this->getConvergence();

        if (converged) {
            break;
        }

        ++it;

        solveEqAndUpdateWellState(well_state);
    } while (it < max_iter);
    return converged;
}

template <typename TypeTag>
void
CompWell<TypeTag>::
solveEqAndUpdateWellState(SingleWellState& well_state)
{
   BVectorWell dx_well(1);

   this->well_equations_.solve(dx_well);

    this->updateWellState(dx_well, well_state);
}

template<typename TypeTag>
void
CompWell<TypeTag>::
apply(BVector& r) const
{
    this->well_equations_.apply(r);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
recoverWellSolutionAndUpdateWellState(const BVector& x,
                                      SingleWellState& well_state)
{
    BVectorWell xw(1);

    this->well_equations_.recoverSolutionWell(x, xw);

    updateWellState(xw, well_state);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updatePrimaryVariablesNewton(const BVectorWell& dwells)
{
    this->primary_variables_.updateNewton(dwells);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updateWellState(const CompWell::BVectorWell& xw,
                SingleWellState& well_state)
{
    updatePrimaryVariablesNewton(xw);
    updateWellStateFromPrimaryVariables(well_state);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updateWellStateFromPrimaryVariables(SingleWellState& well_state) const
{
    well_state.bhp = this->primary_variables_.getBhp().value();

    auto& total_molar_fractions = well_state.total_molar_fractions;
    const auto fluid_state = this->primary_variables_.template toFluidState<Scalar>();
    for (int comp_idx = 0; comp_idx < FluidSystem::numComponents - 1; ++comp_idx) {
        total_molar_fractions[comp_idx] = fluid_state.moleFraction(comp_idx);
    }
    if constexpr (FluidSystem::waterEnabled) {
        well_state.wellbore_water_volume_fraction =
            getValue(this->primary_variables_.getWaterVolumeFraction());
    }

    const Scalar total_rate = this->primary_variables_.getTotalRate().value();
    auto& surface_phase_rates = well_state.surface_phase_rates;
    if (well_state.producer) { // producer
        const auto& surface_cond = this->surface_conditions_;
        for (int p = 0; p < SurfaceConditons::num_phases; ++p) {
            surface_phase_rates[p] = total_rate * getValue(surface_cond.volume_fractions_[p]);
        }
    } else { // injector
        // only gas injection yet
        surface_phase_rates[FluidSystem::gasPhaseIdx] = total_rate;
    }
}

template <typename TypeTag>
bool
CompWell<TypeTag>::
getConvergence() const
{
    bool converged = true;
    for (const auto& val : this->well_equations_.residual()[0]) {
        converged = converged && (std::abs(val) < 1.e-8);
    }
    return converged;
}

template <typename TypeTag>
void
CompWell<TypeTag>::
addWellContributions(SparseMatrixAdapter&) const
{
    assert(false);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updateWellControl(const SummaryState& summary_state,
                  SingleWellState& well_state) const
{
    std::string from;
    if (this->well_ecl_.isInjector()) {
        from = WellInjectorCMode2String(well_state.injection_cmode);
    } else {
        from = WellProducerCMode2String(well_state.production_cmode);
    }
    bool changed = false;
    if (this->well_ecl_.isProducer()) {
        const auto production_controls = this->well_ecl_.productionControls(summary_state);
        const auto current_control = well_state.production_cmode;

        if (production_controls.hasControl(Well::ProducerCMode::BHP) && current_control != WellProducerCMode::BHP) {
            const Scalar bhp_limit = production_controls.bhp_limit;
            const Scalar current_bhp = well_state.bhp;
            if (current_bhp < bhp_limit) {
                well_state.bhp = bhp_limit;
                well_state.production_cmode = WellProducerCMode::BHP;
                changed = true;
            }
        }

        if (!changed && production_controls.hasControl(Well::ProducerCMode::ORAT) && current_control != WellProducerCMode::ORAT) {
            const Scalar current_rate = -well_state.surface_phase_rates[FluidSystem::oilPhaseIdx];
            if (current_rate > production_controls.oil_rate) {
                well_state.production_cmode = WellProducerCMode::ORAT;
                changed = true;
            }
        }

        if (!changed && production_controls.hasControl(Well::ProducerCMode::WRAT) && current_control != WellProducerCMode::WRAT) {
            const Scalar current_rate = -well_state.surface_phase_rates[FluidSystem::waterPhaseIdx];
            if (current_rate > production_controls.water_rate) {
                well_state.production_cmode = WellProducerCMode::WRAT;
                changed = true;
            }
        }

        if (!changed && production_controls.hasControl(Well::ProducerCMode::GRAT) && current_control != WellProducerCMode::GRAT) {
            const Scalar current_rate = -well_state.surface_phase_rates[FluidSystem::gasPhaseIdx];
            if (current_rate > production_controls.gas_rate) {
                well_state.production_cmode = WellProducerCMode::GRAT;
                changed = true;
            }
        }
    } else {
        const auto injection_controls = this->well_ecl_.injectionControls(summary_state);
        const auto current_control = well_state.injection_cmode;
        if (injection_controls.hasControl(Well::InjectorCMode::BHP) && current_control != WellInjectorCMode::BHP) {
            const Scalar bhp_limit = injection_controls.bhp_limit;
            const Scalar current_bhp = well_state.bhp;
            OpmLog::debug(fmt::format("Well {} BHP control check: current_bhp={:.6e}, bhp_limit={:.6e}, exceeds_limit={}",
                                      this->well_ecl_.name(), current_bhp, bhp_limit, (current_bhp > bhp_limit)));
            if (current_bhp > bhp_limit) {
                well_state.bhp = bhp_limit;
                well_state.injection_cmode = WellInjectorCMode::BHP;
                changed = true;
            }
        }
        if (!changed && injection_controls.hasControl(Well::InjectorCMode::RATE) && current_control != WellInjectorCMode::RATE) {
            // InjectorType injector_type = injection_controls.injector_type;
            const Scalar rate_limit = injection_controls.surface_rate;
            // TODO: hack to get the injection rate
            const Scalar current_rate = std::accumulate(well_state.surface_phase_rates.begin(),
                                                        well_state.surface_phase_rates.end(), 0.0);
            OpmLog::debug(fmt::format("Well {} RATE control check: current_rate={:.6e}, rate_limit={:.6e}, bhp={:.6e}, phase_rates=[{}]",
                                      this->well_ecl_.name(), current_rate, rate_limit, well_state.bhp,
                                      fmt::join(well_state.surface_phase_rates, ", ")));
            if (current_rate > rate_limit) {
                OpmLog::debug(fmt::format("Well {} RATE control TRIGGERED: current_rate={:.6e} > rate_limit={:.6e}",
                                          this->well_ecl_.name(), current_rate, rate_limit));
                well_state.injection_cmode = WellInjectorCMode::RATE;
                changed = true;
            }
        }
    }

    if (changed) {
        std::string to;
        if (this->well_ecl_.isInjector()) {
            to = WellInjectorCMode2String(well_state.injection_cmode);
        } else {
            to = WellProducerCMode2String(well_state.production_cmode);
        }
        OpmLog::info(fmt::format("Well {} changed control from {} to {} \n", this->well_ecl_.name(), from, to));
    }
}

template <typename TypeTag>
template <typename T>
void
CompWell<TypeTag>::
updateSurfaceCondition_(const StandardCond& surface_cond,
                        FluidState<T>& fluid_state,
                        const T& water_mass_fraction)
{
    static_assert(std::is_same_v<T, Scalar> || std::is_same_v<T, EvalWell>, "Unsupported type in CompWell::updateSurfaceCondition_");

    fluid_state.setTemperature(surface_cond.temperature);
    fluid_state.setPressure(FluidSystem::oilPhaseIdx, surface_cond.pressure);
    fluid_state.setPressure(FluidSystem::gasPhaseIdx, surface_cond.pressure);

    for (int i = 0; i < FluidSystem::numComponents; ++i) {
        fluid_state.setKvalue(i, fluid_state.wilsonK_(i));
    }

    flashFluidState_(fluid_state);

    for (unsigned compidx = 0; compidx < FluidSystem::numComponents; ++compidx) {
        this->surface_conditions_.mass_fractions_[FluidSystem::oilPhaseIdx][compidx] =
                fluid_state.massFraction(FluidSystem::oilPhaseIdx, compidx);
        this->surface_conditions_.mass_fractions_[FluidSystem::gasPhaseIdx][compidx] =
                fluid_state.massFraction(FluidSystem::gasPhaseIdx, compidx);
    }
    const auto& density_oil = fluid_state.density(FluidSystem::oilPhaseIdx);
    const auto& density_gas = fluid_state.density(FluidSystem::gasPhaseIdx);
    this->surface_conditions_.surface_densities_[FluidSystem::oilPhaseIdx] = density_oil;
    this->surface_conditions_.surface_densities_[FluidSystem::gasPhaseIdx] = density_gas;

    // the hydrocarbon flash splits the hydrocarbon surface volume; with water
    // in the stream, both shares shrink by the water volume fraction
    const auto& so = fluid_state.saturation(FluidSystem::oilPhaseIdx);
    const auto& sg = fluid_state.saturation(FluidSystem::gasPhaseIdx);
    if constexpr (FluidSystem::waterEnabled) {
        fluid_state.setPressure(FluidSystem::waterPhaseIdx, surface_cond.pressure);
        const T rho_w = waterDensity_(fluid_state.pressure(FluidSystem::waterPhaseIdx),
                                      surface_cond.temperature);
        // per unit mass of stream: the volumes of the water and hydrocarbon parts
        const T hc_density = so * density_oil + sg * density_gas;
        const T water_volume = water_mass_fraction / rho_w;
        const T hc_volume = (1. - water_mass_fraction) / hc_density;
        const T water_volume_fraction = water_volume / (water_volume + hc_volume);
        this->surface_conditions_.surface_densities_[FluidSystem::waterPhaseIdx] = rho_w;
        this->surface_conditions_.volume_fractions_[FluidSystem::waterPhaseIdx] = water_volume_fraction;
        this->surface_conditions_.volume_fractions_[FluidSystem::oilPhaseIdx] = (1. - water_volume_fraction) * so;
        this->surface_conditions_.volume_fractions_[FluidSystem::gasPhaseIdx] = (1. - water_volume_fraction) * sg;
    } else {
        static_cast<void>(water_mass_fraction);
        this->surface_conditions_.volume_fractions_[FluidSystem::oilPhaseIdx] = so;
        this->surface_conditions_.volume_fractions_[FluidSystem::gasPhaseIdx] = sg;
    }
}

template <typename TypeTag>
template <typename T>
void
CompWell<TypeTag>::
flashFluidState_(FluidState<T>& fluid_state)
{
    static_assert(std::is_same_v<T, Scalar> || std::is_same_v<T, EvalWell>, "Unsupported type in CompWell::flashFluidState_");

    // The wellbore flash is a free function so it can be unit tested in
    // isolation (see tests/test_compwell_jacobian.cpp).
    flashWellboreFluidState(fluid_state);
}

} // end of namespace Opm
