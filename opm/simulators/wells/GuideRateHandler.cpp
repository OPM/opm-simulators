/*
  Copyright 2025 Equinor ASA

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
#include <config.h>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>
#include <opm/common/TimingMacros.hpp>

namespace Opm {

// ---------------------------------------------
// Constructor for the GuideRateHandler class
// ---------------------------------------------
template <class Scalar>
GuideRateHandler<Scalar>::GuideRateHandler(
    const Schedule& schedule,
    const PhaseUsage& phase_usage,
    const SummaryState& summary_state,
    const Parallel::Communication& comm,
    GuideRate& guide_rate
) :
          schedule_{schedule}
        , phase_usage_{phase_usage}
        , summary_state_{summary_state}
        , comm_{comm}
        , guide_rate_{guide_rate}
{ }


// ----------------------------------------------------------
// Public methods for GuideRateHandler sorted alphabetically
// -----------------------------------------------------------

template <class Scalar>
DeferredLogger&
GuideRateHandler<Scalar>::
deferredLogger()
{
    assert(this->deferred_logger_ != nullptr);
    return *this->deferred_logger_;
}

#ifdef RESERVOIR_COUPLING_ENABLED
template <class Scalar>
void
GuideRateHandler<Scalar>::
receiveMasterGroupPotentialsFromSlaves()
{
    assert(this->isReservoirCouplingMaster());
    auto& rescoup_master = this->reservoirCouplingMaster();
    const auto& comm = rescoup_master.getComm();
    if (comm.rank() == 0) {
        rescoup_master.receivePotentialsFromSlaves();
    }
}
#endif

#ifdef RESERVOIR_COUPLING_ENABLED
template<class Scalar>
void
GuideRateHandler<Scalar>::
sendSlaveGroupPotentialsToMaster()
{
    assert(this->isReservoirCouplingSlave());
    auto& rescoup_slave = this->reservoirCouplingSlave();
    const auto& slave_master_group_map = rescoup_slave.getSlaveToMasterGroupNameMap();
    std::vector<Potentials> potentials;
    if (this->comm_.rank() == 0) {
        // NOTE: Traversing a std::map is guaranteed to iterate the keys in lexical order for
        //   std::string keys, so the master can from this order determine which potentials
        //   correspond to which slave group.
        for (const auto& item : slave_master_group_map) {
            const auto& slave_group_name = item.first;
            Potentials pot;
            // For injection groups, we do not have potentials. In that case,
            //  we will send dummy values (0.0) for the potentials.
            if (this->guide_rate_.hasPotentials(slave_group_name)) {
                const auto& gr_pot = this->guide_rate_.getPotentials(slave_group_name);
                pot.oil_rate = gr_pot.oil_rat;
                pot.gas_rate = gr_pot.gas_rat;
                pot.water_rate = gr_pot.wat_rat;
            }
            potentials.push_back(pot);
        }
        rescoup_slave.sendPotentialsToMaster(potentials);
    }
}
#endif

template <class Scalar>
void
GuideRateHandler<Scalar>::setLogger(DeferredLogger *deferred_logger)
{
    deferred_logger_ = deferred_logger;
#ifdef RESERVOIR_COUPLING_ENABLED
    if (reservoir_coupling_master_) {
        reservoir_coupling_master_->setDeferredLogger(deferred_logger);
    }
    if (reservoir_coupling_slave_) {
        reservoir_coupling_slave_->setDeferredLogger(deferred_logger);
    }
#endif
}

template <class Scalar>
void
GuideRateHandler<Scalar>::
updateGuideRates(
    const int report_step_idx,
    const double sim_time,
    WellState<Scalar> &well_state,
    GroupState<Scalar> &group_state
)
{
    OPM_TIMEFUNCTION();
    auto num_phases = this->phase_usage_.num_phases;
    UpdateGuideRates updater {
        *this, report_step_idx, sim_time, well_state, group_state, num_phases
    };
    updater.update();

}

// ------------------------------------------
// Inner class UpdateGuideRates constructor
// ------------------------------------------

template <class Scalar>
GuideRateHandler<Scalar>::UpdateGuideRates::
UpdateGuideRates(
    GuideRateHandler<Scalar>& parent,
    const int report_step_idx,
    const double sim_time,
    const WellState<Scalar>& well_state,
    const GroupState<Scalar>& group_state,
    const int num_phases
) : parent_{parent}
  , report_step_idx_{report_step_idx}
  , sim_time_{sim_time}
  , well_state_{well_state}
  , group_state_{group_state}
  , num_phases_{num_phases}
  , unit_system_{parent_.schedule_.getUnits()}
{
}

// ------------------------------------------------------------------
// Inner class UpdateGuideRates public methods sorted alphabetically
// ------------------------------------------------------------------


template <class Scalar>
void
GuideRateHandler<Scalar>::UpdateGuideRates::
update()
{
    this->guideRate().updateGuideRateExpiration(this->sim_time_, this->report_step_idx_);
    const Group& group = this->schedule().getGroup("FIELD", this->report_step_idx_);
    std::vector<Scalar> pot(this->num_phases_, 0.0);
    this->updateGuideRatesForProductionGroups_(group, pot);
    this->updateGuideRatesForInjectionGroups_(group);
    this->updateGuideRatesForWells_();
}

// --------------------------------------------------------------------
// Inner class UpdateGuideRates private methods sorted alphabetically
// --------------------------------------------------------------------

#ifdef RESERVOIR_COUPLING_ENABLED
template<class Scalar>
bool
GuideRateHandler<Scalar>::UpdateGuideRates::
isMasterGroup_(const Group& group)
{
    if (this->isReservoirCouplingMaster()) {
        return this->reservoirCouplingMaster().isMasterGroup(group.name());
    }
    else {
        return false;
    }
}
#endif

template<class Scalar>
void
GuideRateHandler<Scalar>::UpdateGuideRates::
updateGuideRatesForInjectionGroups_(const Group& group)
{
    OPM_TIMEFUNCTION();
    if (group.isProductionGroup() && !group.isInjectionGroup()) {
        // This is a pure production group which will be handled by
        // updateGuideRatesForProductionGroups_()
        return;
    }
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule().getGroup(group_name, this->report_step_idx_);
        // compute guide rates for sub groups first recursively
        this->updateGuideRatesForInjectionGroups_(group_tmp);
    }
    const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
    const auto &pu = this->phaseUsage();
    for (Phase phase : all) {
        if(!group.hasInjectionControl(phase))
            continue;

        std::optional<Scalar> guide_rate_value;
        const auto& controls = group.injectionControls(phase, this->summaryState());
        switch (controls.guide_rate_def){
        case Group::GuideRateInjTarget::RATE:
            break;
        case Group::GuideRateInjTarget::VOID:
        {
            guide_rate_value = std::max(
                Scalar(0.0), this->group_state_.injection_vrep_rate(group.name())
            );
            break;
        }
        case Group::GuideRateInjTarget::NETV:
        {
            guide_rate_value = this->group_state_.injection_vrep_rate(group.name());
            const std::vector<Scalar>& injRES
                                = this->group_state_.injection_reservoir_rates(group.name());
            if (phase != Phase::OIL && pu.phase_used[BlackoilPhases::Liquid])
                guide_rate_value = *guide_rate_value - injRES[pu.phase_pos[BlackoilPhases::Liquid]];
            if (phase != Phase::GAS && pu.phase_used[BlackoilPhases::Vapour])
                guide_rate_value = *guide_rate_value - injRES[pu.phase_pos[BlackoilPhases::Vapour]];
            if (phase != Phase::WATER && pu.phase_used[BlackoilPhases::Aqua])
                guide_rate_value = *guide_rate_value - injRES[pu.phase_pos[BlackoilPhases::Aqua]];

            guide_rate_value = std::max(Scalar(0.0), *guide_rate_value);
            break;
        }
        case Group::GuideRateInjTarget::RESV:
            OPM_DEFLOG_THROW(
                std::runtime_error,
                "GUIDE PHASE RESV not implemented. Group " + group.name(),
                this->deferredLogger()
            );
        case Group::GuideRateInjTarget::NO_GUIDE_RATE:
            break;
        default:
            OPM_DEFLOG_THROW(
                std::logic_error,
                "Invalid GuideRateInjTarget in updateGuideRatesForInjectionGroups",
                this->deferredLogger()
            );
        }

        if (guide_rate_value) {
            guide_rate_value = this->unit_system_.from_si(
                UnitSystem::measure::rate, *guide_rate_value
            );
        }
        this->guideRate().compute(
            group.name(), phase, this->report_step_idx_, guide_rate_value
        );
    }
}

template<class Scalar>
void
GuideRateHandler<Scalar>::UpdateGuideRates::
updateGuideRatesForProductionGroups_(const Group& group, std::vector<Scalar>& pot)
{
    OPM_TIMEFUNCTION();
    if (group.isInjectionGroup() && !group.isProductionGroup()) {
        // This is a pure injection group which will be handled by
        // updateGuideRatesForInjectionGroups_()
        return;
    }
#ifdef RESERVOIR_COUPLING_ENABLED
    if (this->isMasterGroup_(group)) {
        // This is a master group, so we do not need to compute potentials
        // for sub groups. We just set the guide rates for the master group
        // as submitted from its corresponding slave group.
        this->updateProductionGroupPotentialFromSlaveGroup_(group, pot);
    }
    else {
        this->updateProductionGroupPotentialFromSubGroups(group, pot);
    }
#else
    this->updateProductionGroupPotentialFromSubGroups(group, pot);
#endif
    std::array<Scalar,3> potentials{};
    const auto &pu = this->phaseUsage();
    auto& [oil_pot, gas_pot, water_pot] = potentials;
    if (pu.phase_used[BlackoilPhases::Liquid])
        oil_pot = pot[pu.phase_pos[BlackoilPhases::Liquid]];

    if (pu.phase_used[BlackoilPhases::Vapour])
        gas_pot = pot[pu.phase_pos[BlackoilPhases::Vapour]];

    if (pu.phase_used[BlackoilPhases::Aqua])
        water_pot = pot[pu.phase_pos[BlackoilPhases::Aqua]];

    // Synchronize potentials across all ranks
    this->comm().sum(potentials.data(), potentials.size());
    oil_pot = this->unit_system_.from_si(UnitSystem::measure::liquid_surface_rate, oil_pot);
    water_pot = this->unit_system_.from_si(UnitSystem::measure::liquid_surface_rate, water_pot);
    gas_pot = this->unit_system_.from_si(UnitSystem::measure::gas_surface_rate, gas_pot);
    this->guideRate().compute(
        group.name(), this->report_step_idx_, this->sim_time_, oil_pot, gas_pot, water_pot
    );
}

template<class Scalar>
void
GuideRateHandler<Scalar>::UpdateGuideRates::
updateGuideRatesForWells_()
{
    OPM_TIMEFUNCTION();
    const auto &pu = this->phaseUsage();
    for (const auto& well : this->schedule().getWells(this->report_step_idx_)) {
        std::array<Scalar,3> potentials{};
        auto& [oilpot, gaspot, waterpot] = potentials;

        const auto& well_index = this->well_state_.index(well.name());
        if (well_index.has_value() && this->well_state_.wellIsOwned(well_index.value(), well.name()))
        {
            // the well is found and owned
            const auto& ws = this->well_state_.well(well_index.value());
            const auto& wpot = ws.well_potentials;
            if (pu.phase_used[BlackoilPhases::Liquid] > 0)
                oilpot = wpot[pu.phase_pos[BlackoilPhases::Liquid]];

            if (pu.phase_used[BlackoilPhases::Vapour] > 0)
                gaspot = wpot[pu.phase_pos[BlackoilPhases::Vapour]];

            if (pu.phase_used[BlackoilPhases::Aqua] > 0)
                waterpot = wpot[pu.phase_pos[BlackoilPhases::Aqua]];
        }
        this->comm().sum(potentials.data(), potentials.size());
        oilpot = this->unit_system_.from_si(UnitSystem::measure::liquid_surface_rate, oilpot);
        waterpot = this->unit_system_.from_si(UnitSystem::measure::liquid_surface_rate, waterpot);
        gaspot = this->unit_system_.from_si(UnitSystem::measure::gas_surface_rate, gaspot);
        this->guideRate().compute(
            well.name(), this->report_step_idx_, this->sim_time_, oilpot, gaspot, waterpot
        );
    }
}

#ifdef RESERVOIR_COUPLING_ENABLED
template<class Scalar>
void
GuideRateHandler<Scalar>::UpdateGuideRates::
updateProductionGroupPotentialFromSlaveGroup_(const Group& group, std::vector<Scalar>& pot)
{
    assert(this->isReservoirCouplingMaster());
    auto& rescoup_master = this->reservoirCouplingMaster();
    const auto& slave_pot = rescoup_master.getSlaveGroupPotentials(group.name());
    auto& pu = this->phaseUsage();
    // TODO: Here we should check that the master uses the same phases as the
    //   slave.
    pot[pu.phase_pos[BlackoilPhases::Liquid]] = slave_pot.oil_rate;
    pot[pu.phase_pos[BlackoilPhases::Vapour]] = slave_pot.gas_rate;
    pot[pu.phase_pos[BlackoilPhases::Aqua]] = slave_pot.water_rate;
}
#endif

template<class Scalar>
void
GuideRateHandler<Scalar>::UpdateGuideRates::
updateProductionGroupPotentialFromSubGroups(const Group& group, std::vector<Scalar>& pot)
{
    for (const std::string& group_name : group.groups()) {
        std::vector<Scalar> this_pot(this->num_phases_, 0.0);
        const Group& group_tmp = this->schedule().getGroup(group_name, this->report_step_idx_);

        // compute group potentials and guide rates for sub groups recursively
        this->updateGuideRatesForProductionGroups_(group_tmp, this_pot);

        // If group_tmp is *not* available for group control at a higher level,
        // it should not contribute to the group potentials and guide rates of the parent group
        const auto current_group_control = this->group_state_.production_control(group_name);
        if (current_group_control != Group::ProductionCMode::FLD
                && current_group_control != Group::ProductionCMode::NONE) {
            continue;
        }

        // Apply potential for group_tmp to the parent's pot
        auto gefac = group_tmp.getGroupEfficiencyFactor();

        for (int phase = 0; phase < this->num_phases_; phase++) {
            pot[phase] += gefac*this_pot[phase];
        }

    }

    // If this is group on the lowest level in the group tree, add contribution from its wells
    for (const std::string& well_name : group.wells()) {
        const auto& well_tmp = this->schedule().getWell(well_name, this->report_step_idx_);
        const auto wefac = well_tmp.getEfficiencyFactor();

        if (well_tmp.isInjector())
            continue;

        if (well_tmp.getStatus() == Well::Status::SHUT)
            continue;
        const auto& well_index = this->well_state_.index(well_name);
        if (!well_index.has_value()) // the well is not found
            continue;

        if (!this->well_state_.wellIsOwned(well_index.value(), well_name) ) // Only sum once
        {
            continue;
        }

        const auto& ws = this->well_state_.well(well_index.value());
        for (int phase = 0; phase < this->num_phases_; phase++) {
            pot[phase] += wefac * ws.well_potentials[phase];
        }
    }
}


template class GuideRateHandler<double>;

#if FLOW_INSTANTIATE_FLOAT
template class GuideRateHandler<float>;
#endif

} // namespace Opm

