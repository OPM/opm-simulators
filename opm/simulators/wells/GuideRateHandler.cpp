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

#include <opm/simulators/wells/GuideRateHandler.hpp>

#include <opm/common/TimingMacros.hpp>

#include <opm/output/data/GuideRateValue.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace Opm {

// ---------------------------------------------
// Constructor for the GuideRateHandler class
// ---------------------------------------------
template<typename Scalar, typename IndexTraits>
GuideRateHandler<Scalar, IndexTraits>::GuideRateHandler(
    BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model,
    const Schedule& schedule,
    const SummaryState& summary_state,
    const Parallel::Communication& comm
) :
          well_model_{well_model}
        , schedule_{schedule}
        , summary_state_{summary_state}
        , comm_{comm}
        , guide_rate_{well_model_.guideRate()}
{
}


// ----------------------------------------------------------
// Public methods for GuideRateHandler sorted alphabetically
// -----------------------------------------------------------

template<typename Scalar, typename IndexTraits>
DeferredLogger&
GuideRateHandler<Scalar, IndexTraits>::
deferredLogger()
{
    assert(this->deferred_logger_ != nullptr);
    return *this->deferred_logger_;
}

// This is an alternative to plotting the summary keywords:
// - GOPGR (group oil production guide rate),
// - WOPGR (well oil production guide rate),
// - GWPGR (group water production guide rate),
// - WWPGR (well water production guide rate),
// - GGPGR (group gas production guide rate),
// - WGPGR (well gas production guide rate),
// - GGIGR (group gas injection guide rate),
// - WGIGR (well gas injection guide rate),
// - GWIGR (group water injection guide rate),
// - WWIGR (well water injection guide rate),
// - GVPGR (group reservoir volume production guide rate),
// - WVPGR (well reservoir volume production guide rate),
template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::
debugDumpGuideRates(const int report_step_idx, const double sim_time)
{
    if (this->comm_.rank() == 0) {
        GuideRateDumper dumper{*this, report_step_idx, sim_time};
        dumper.dumpGuideRates();
    }
}


#ifdef RESERVOIR_COUPLING_ENABLED
template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::
receiveMasterGroupPotentialsFromSlaves()
{
    assert(this->isReservoirCouplingMaster());
    auto& rescoup_master = this->reservoirCouplingMaster();
    rescoup_master.receivePotentialsFromSlaves();
}
#endif

#ifdef RESERVOIR_COUPLING_ENABLED
template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::
sendSlaveGroupPotentialsToMaster(const GroupState<Scalar>& group_state)
{
    assert(this->isReservoirCouplingSlave());
    if (this->comm_.rank() == 0) {
        // NOTE: Traversing a std::map is guaranteed to iterate the keys in lexical order for
        //   std::string keys, so the master can from this order determine which potentials
        //   correspond to which slave group.
        auto& rescoup_slave = this->reservoirCouplingSlave();
        const auto& slave_master_group_map = rescoup_slave.getSlaveToMasterGroupNameMap();
        std::vector<Potentials> potentials;
        for (const auto& item : slave_master_group_map) {
            const auto& slave_group_name = item.first;
            Potentials pot;
            // For injection groups, we do not have potentials. In that case,
            //  we will send dummy values (0.0) for the potentials.
            if (this->guide_rate_.hasPotentials(slave_group_name)) {
                const auto& gr_pot = group_state.get_production_group_potential(slave_group_name);
                pot[Potentials::Phase::Oil] = gr_pot.oil_rate;
                pot[Potentials::Phase::Gas] = gr_pot.gas_rate;
                pot[Potentials::Phase::Water] = gr_pot.water_rate;
            }
            potentials.push_back(pot);
        }
        rescoup_slave.sendPotentialsToMaster(potentials);
    }
}
#endif

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::setLogger(DeferredLogger *deferred_logger)
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

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::
updateGuideRates(const int report_step_idx,
                 const double sim_time,
                 const WellState<Scalar, IndexTraits>& well_state,
                 GroupState<Scalar>& group_state)
{
    OPM_TIMEFUNCTION();
    int num_phases = this->well_model_.phaseUsage().numActivePhases();
    UpdateGuideRates updater {
        *this, report_step_idx, sim_time, well_state, group_state, num_phases
    };
    updater.update();

}

// ------------------------------------------
// Inner class GuideRateDumper constructor
// ------------------------------------------

// NOTE: See debugDumpGuideRates() above for more information on the
//       purpose of this class. It is used to dump the guide rates
//       to the terminal in a human-readable format.
template<typename Scalar, typename IndexTraits>
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
GuideRateDumper(
    GuideRateHandler<Scalar, IndexTraits> &parent, const int report_step_idx, const double sim_time
) : parent_{parent}
  , report_step_idx_{report_step_idx}
  , sim_time_{sim_time}
  , well_model_{parent.wellModel()}
  , schedule_{parent.schedule()}
  , comm_{parent.getComm()}
{
}

// ---------------------------------------------------------------------
// Public methods for inner class GuideRateDumper sorted alphabetically
// ---------------------------------------------------------------------

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
dumpGuideRates()
{
    if (this->comm_.rank() == 0) {
        auto calcGR = BlackoilWellModelGuideRates{this->well_model_};
        this->well_guide_rates_ = calcGR.calculateWellGuideRates(this->report_step_idx_);
        this->group_guide_rates_ = calcGR.calculateAllGroupGuideRates(this->report_step_idx_);
        const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);
        printHeader_();
        int level = 0;  // 0 is the root level
        this->dumpGuideRatesRecursive_(group, level);
        printFooter_();
    }
}

// ----------------------------------------------------------
// Private methods for GuideRateHandler sorted alphabetically
// -----------------------------------------------------------


template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
dumpGuideRatesRecursive_(const Group& group, int level)
{
    if (group.name() != "FIELD") {
        this->printGroupGuideRates_(group, level);
    }
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule_.getGroup(group_name, this->report_step_idx_);
        this->dumpGuideRatesRecursive_(group_tmp, level+1);
    }
    for (const std::string& well_name : group.wells()) {
        const auto& well_tmp = this->schedule_.getWell(well_name, this->report_step_idx_);
        printWellGuideRates_(well_tmp, level+1);
    }
}

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
getGroupGuideRatesInjection_(
    const Group& group,
    const data::GroupGuideRates& group_guide_rate,
    std::vector<std::string>& msg_items
) const
{
    // NOTE: See comment in getGroupGuideRatesProduction_() for a discussion of the
    // different types of group guiderate enums.
    const auto& wm_guide_rate = this->well_model_.guideRate();
    const auto& name = group.name();
    using Value = data::GuideRateValue::Item;
    static const std::array<std::tuple<Phase, Value, std::string_view>, 3> items = {{
        {Phase::OIL, Value::Oil, "oil"},
        {Phase::GAS, Value::Gas, "gas"},
        {Phase::WATER, Value::Water, "water"},
    }};
    const auto& guide_rate_value = group_guide_rate.production;
    for (const auto& [phase, value, phase_str] : items) {
        if (wm_guide_rate.has(name, phase)) {
            if (guide_rate_value.has(value)) {
                msg_items.push_back(
                    fmt::format(
                        "{}(Group Inj, {}): {}",
                        name,
                        phase_str,
                        guide_rate_value.get(value)
                    )
                );
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
getGroupGuideRatesProduction_(
    const Group& group,
    const data::GroupGuideRates& group_guide_rate,
    std::vector<std::string>& msg_items
) const
{
    // NOTE: There are many different types of group guide rate phase target enums that are
    //       used in different contexts:
    //  - Group.hpp defines Group::GuideRateProdTarget and Group::GuideRateInjTarget
    //    - Production: OIL, WAT, GAS, LIQ, RES, COMB, WGA, CVAL, INJV, POTN, FORM
    //      - These are defined by item 9 of GCONPROD and used if item 10 of GCONPROD is not
    //        set to FORM.
    //    - Injection: RATE, VOID, NETV, RESV
    //  - GuideRateModel.hpp defines GuideRateModel::Target
    //    - OIL, LIQ, GAS, WAT, RES, COMB, NONE
    //    - These are defined by item 2 of GUIDERAT.
    //    - These are used if item 10 of GCONPROD is set to FORM, then the guide rate set in
    //      item 9 of GCONPROD is ignored.
    //  - Phase.hpp defines the Phase enum which is used by BlackoilWellModelGuideRates.cpp and
    //    GuideRate.cpp
    //    - OIL, GAS, WATER, SOLVENT, POLYMER, ENERGY, POLYMW, FOAM, BRINE
    //    - Only OIL, GAS, WATER are used in the context of guide rates and only for
    //      injection groups. And only GAS and WATER are checked for in
    //      getGuideRateInjectionGroupValues() in BlackoilWellModelGuideRates.cpp.
    //  - GuideRateValue.hpp defines GuideRateValue::Item which is used by
    //    BlackoilWellModelGuideRates.cpp to assign guide rates. The values are:
    //    - Oil, Gas, Water, ResV
    //    - However, ResV is not implemented, see getGuideRateValues() in
    //      BlackoilWellModelGuideRates.cpp. This, despite the summary keywords GVPGR and WVPGR
    //      are using it, see Summary.cpp. TODO: Check if this is a bug.
    const auto& wm_guide_rate = this->well_model_.guideRate();
    const auto& name = group.name();
    if (wm_guide_rate.has(name)) {  // Check if group has production guiderates
        using Value = data::GuideRateValue::Item;
        static const std::array<std::tuple<Value, std::string_view>, 4> value_types = {{
            {Value::Oil, "oil"},
            {Value::Gas, "gas"},
            {Value::Water, "water"},
            {Value::ResV, "resv"}
        }};
        const auto& guide_rate_value = group_guide_rate.production;
        for (const auto& [value_type, phase_str] : value_types) {
            if (guide_rate_value.has(value_type)) {
                msg_items.push_back(
                    fmt::format(
                        "{}(Group Prod, {}): {}",
                        name,
                        phase_str,
                        guide_rate_value.get(value_type)
                    )
                );
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
printGroupGuideRates_(const Group& group, int level)
{
    const auto& name = group.name();
    assert( level >= 1 );
    std::string indent(level - 1, ' ');
    auto gr_itr = this->group_guide_rates_.find(name);
    if (gr_itr == this->group_guide_rates_.end()) {
        this->deferredLogger().debug(
            fmt::format("{}{}: <NO GUIDERATE>", indent, name)
        );
        return;
    }
    const auto& group_guide_rate = gr_itr->second;
    std::vector<std::string> msg_items;
    getGroupGuideRatesProduction_(group, group_guide_rate, msg_items);
    getGroupGuideRatesInjection_(group, group_guide_rate, msg_items);
    if (msg_items.empty()) {
        this->deferredLogger().debug(
            fmt::format("{}{}(Group): <NO GUIDERATES FOUND>", indent, name)
        );
        return;
    }
    this->deferredLogger().debug(
        fmt::format("{}{}", indent, fmt::join(msg_items, ", "))
    );
}

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
printHeader_()
{
    this->deferredLogger().debug(
        fmt::format(
            "\n---- GUIDE RATE REPORT: step {}, simtime {} ----",
            this->report_step_idx_,
            this->sim_time_
        )
    );
}

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
printFooter_()
{
    this->deferredLogger().debug(
        fmt::format(
            "---- END GUIDE RATE REPORT: step {}, simtime {} ----\n",
            this->report_step_idx_,
            this->sim_time_
        )
    );
}

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::GuideRateDumper::
printWellGuideRates_(const Well& well, int level)
{
    const auto& name = well.name();
    assert( level >= 1 );
    std::string indent(level - 1, ' ');
    auto gr_itr = this->well_guide_rates_.find(name);
    if (gr_itr == this->well_guide_rates_.end()) {
        this->deferredLogger().debug(
            fmt::format("{}{}: <NO GUIDERATE>", indent, name)
        );
        return;
    }
    const auto& guide_rate_value = gr_itr->second;
    std::vector<std::string> msg_items;
    using Value = data::GuideRateValue::Item;
    static const std::array<std::tuple<Value, std::string_view>, 3> value_types = {{
        {Value::Oil, "oil"},
        {Value::Gas, "gas"},
        {Value::Water, "water"}
    }};
    const std::string well_type = well.isInjector() ? "Inj" : "Prod";
    for (const auto& [value_type, phase_str] : value_types) {
        if (guide_rate_value.has(value_type)) {
            msg_items.push_back(
                fmt::format(
                    "{}(Well {}, {}): {}",
                    name,
                    well_type,
                    phase_str,
                    guide_rate_value.get(value_type)
                )
            );
        }
    }
    if (msg_items.empty()) {
        this->deferredLogger().debug(
            fmt::format("{}{}(Well {}): <NO GUIDERATES FOUND>", indent, name, well_type)
        );
        return;
    }
    this->deferredLogger().debug(
        fmt::format("{}{}", indent, fmt::join(msg_items, ", "))
    );
}

// ------------------------------------------
// Inner class UpdateGuideRates constructor
// ------------------------------------------

template<typename Scalar, typename IndexTraits>
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
UpdateGuideRates(
    GuideRateHandler<Scalar, IndexTraits>& parent,
    const int report_step_idx,
    const double sim_time,
    const WellState<Scalar, IndexTraits>& well_state,
    GroupState<Scalar>& group_state,
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


template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
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
template<typename Scalar, typename IndexTraits>
bool
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
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

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
updateGuideRatesForInjectionGroups_(const Group& group)
{
    OPM_TIMEFUNCTION();
    // NOTE: Even if this is a pure production group, we still need to compute the
    // group potentials since they may be used by an injection group at a higher
    // level in the group hierarchy.
    for (const std::string& group_name : group.groups()) {
        const Group& group_tmp = this->schedule().getGroup(group_name, this->report_step_idx_);
        // compute guide rates for sub groups first recursively
        this->updateGuideRatesForInjectionGroups_(group_tmp);
    }
    const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
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
            const auto& pu = this->phaseUsage();

            guide_rate_value = this->group_state_.injection_vrep_rate(group.name());
            const std::vector<Scalar>& injRES
                                = this->group_state_.injection_reservoir_rates(group.name());
            if (phase != Phase::OIL && pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
                const int phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
                guide_rate_value = *guide_rate_value - injRES[phase_pos];
            }
            if (phase != Phase::GAS && pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
                const int phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
                guide_rate_value = *guide_rate_value - injRES[phase_pos];
            }
            if (phase != Phase::WATER && pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
                const int phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
                guide_rate_value = *guide_rate_value - injRES[phase_pos];
            }

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

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
updateGuideRatesForProductionGroups_(const Group& group, std::vector<Scalar>& pot)
{
    OPM_TIMEFUNCTION();
    // NOTE: Even if this is a pure injection group, we still need to compute the
    // group potentials since they may be used by a production group at a higher
    // level in the group hierarchy. See MOD4_UDQ_ACTIONX.DATA for an example.
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
    const auto& pu = this->phaseUsage();

    std::array<Scalar,3> potentials{};
    auto& [oil_pot, gas_pot, water_pot] = potentials;
    if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        oil_pot = pot[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)];
    }

    if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        gas_pot = pot[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)];
    }

    if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
        water_pot = pot[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)];
    }

    // Synchronize potentials across all ranks
    this->comm().sum(potentials.data(), potentials.size());
    oil_pot = this->unit_system_.from_si(UnitSystem::measure::liquid_surface_rate, oil_pot);
    water_pot = this->unit_system_.from_si(UnitSystem::measure::liquid_surface_rate, water_pot);
    gas_pot = this->unit_system_.from_si(UnitSystem::measure::gas_surface_rate, gas_pot);
    this->guideRate().compute(
        group.name(), this->report_step_idx_, this->sim_time_, oil_pot, gas_pot, water_pot
    );
    // NOTE: For reservoir coupling: We will later need to send the slave group potentials
    //   to the master. The above call to GuideRate::compute() will as a side-effect store
    //   them in the GuideRate object. But those potentials are meant to be internal to the
    //   GuideRate object, so we choose to store the potentials in the GroupState object
    //   below for this reason.
    this->group_state_.update_group_production_potential(
        group.name(), oil_pot, gas_pot, water_pot
    );
}

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
updateGuideRatesForWells_()
{
    OPM_TIMEFUNCTION();

    const auto& pu = this->phaseUsage();

    const auto o_pos = pu.phaseIsActive(IndexTraits::oilPhaseIdx)
        ? pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx) : -1;

    const auto g_pos = pu.phaseIsActive(IndexTraits::gasPhaseIdx)
        ? pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx) : -1;

    const auto w_pos = pu.phaseIsActive(IndexTraits::waterPhaseIdx)
        ? pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx) : -1;

    constexpr auto o_ix = std::size_t{0};
    constexpr auto g_ix = o_ix + 1;
    constexpr auto w_ix = g_ix + 1;
    constexpr auto npot = w_ix + 1;

    const auto& wnames = this->schedule()[this->report_step_idx_].well_order();

    auto all_well_pot = std::vector<double>(npot * wnames.size());
    auto well_pot = all_well_pot.begin();

    for (const auto& wname : wnames) {
        const auto well_index = this->well_state_.index(wname);

        if (well_index.has_value() && this->well_state_.wellIsOwned(*well_index, wname)) {
            // the well is found and owned
            const auto& wpot = this->well_state_.well(*well_index).well_potentials;

            if (o_pos >= 0) { well_pot[o_ix] = static_cast<double>(wpot[o_pos]); }
            if (g_pos >= 0) { well_pot[g_ix] = static_cast<double>(wpot[g_pos]); }
            if (w_pos >= 0) { well_pot[w_ix] = static_cast<double>(wpot[w_pos]); }
        }

        well_pot += npot;
    }

    this->comm().sum(all_well_pot.data(), all_well_pot.size());

    using M = UnitSystem::measure;

    well_pot = all_well_pot.begin();
    for (const auto& wname : wnames) {
        const auto o_pot = this->unit_system_.from_si(M::liquid_surface_rate, well_pot[o_ix]);
        const auto g_pot = this->unit_system_.from_si(M::gas_surface_rate   , well_pot[g_ix]);
        const auto w_pot = this->unit_system_.from_si(M::liquid_surface_rate, well_pot[w_ix]);

        this->guideRate().compute
            (wname, this->report_step_idx_, this->sim_time_, o_pot, g_pot, w_pot);

        well_pot += npot;
    }
}

#ifdef RESERVOIR_COUPLING_ENABLED
template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
updateProductionGroupPotentialFromSlaveGroup_(const Group& group, std::vector<Scalar>& pot)
{
    assert(this->isReservoirCouplingMaster());
    auto& rescoup_master = this->reservoirCouplingMaster();
    const auto& slave_pot = rescoup_master.getSlaveGroupPotentials(group.name());
    const auto& pu = this->phaseUsage();
    // TODO: Here we should check that the master uses the same phases as the
    //   slave.
    pot[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] = slave_pot[Potentials::Phase::Oil];
    pot[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] = slave_pot[Potentials::Phase::Gas];
    pot[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] = slave_pot[Potentials::Phase::Water];
}
#endif

template<typename Scalar, typename IndexTraits>
void
GuideRateHandler<Scalar, IndexTraits>::UpdateGuideRates::
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

        // Only include producers in group potentials for production groups
        if (well_tmp.isInjector())
            continue;

        const auto& well_index = this->well_state_.index(well_name);
        if (!well_index.has_value()) // the well is not found
            continue;

        if (!this->well_state_.wellIsOwned(well_index.value(), well_name) ) // Only sum once
        {
            continue;
        }

        const auto& ws = this->well_state_.well(well_index.value());
        if (ws.status == Well::Status::SHUT)
            continue;

        for (int phase = 0; phase < this->num_phases_; phase++) {
            pot[phase] += wefac * ws.well_potentials[phase];
        }
    }
}


template class GuideRateHandler<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class GuideRateHandler<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
