/*
  Copyright 2021 Equinor ASA.

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

#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>

#include <cmath>
#include <optional>
#include <string>

#include <fmt/format.h>

namespace Opm {

template<typename TypeTag>
GasLiftStage2<TypeTag>::
GasLiftStage2(
    const BlackoilWellModel &well_model,
    const Simulator &ebos_simulator,
    DeferredLogger &deferred_logger,
    WellState &well_state,
    GLiftProdWells &prod_wells,
    GLiftOptWells &glift_wells,
    GLiftWellStateMap &state_map
) :
    deferred_logger_{deferred_logger},
    ebos_simulator_{ebos_simulator},
    well_model_{well_model},
    well_state_{well_state},
    prod_wells_{prod_wells},
    stage1_wells_{glift_wells},
    well_state_map_{state_map},
    report_step_idx_{ebos_simulator_.episodeIndex()},
    summary_state_{ebos_simulator_.vanguard().summaryState()},
    schedule_{ebos_simulator.vanguard().schedule()},
    phase_usage_{well_model_.phaseUsage()},
    glo_{schedule_.glo(report_step_idx_)},
    debug_{false}
{
//    this->time_step_idx_
//        = this->ebos_simulator_.model().newtonMethod().currentTimeStep();
    this->nonlinear_iteration_idx_
        = this->ebos_simulator_.model().newtonMethod().numIterations();

}

/********************************************
 * Public methods in alphabetical order
 ********************************************/

// runOptimize():
//
// If a group has any production rate constraints, and/or a limit on
// its total rate of lift gas supply, allocates lift gas
// preferentially to the wells that gain the most benefit from
// it. Lift gas increments are allocated in turn to the well that
// currently has the largest weighted incremental gradient. The
// procedure takes account of any limits on the group production rate
// or lift gas supply applied to any level of group, including the FIELD level group.
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
runOptimize()
{
    const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);

    optimizeGroupsRecursive_(group);

}


/********************************************
 * Private methods in alphabetical order
 ********************************************/

// Update GasLiftWellState and WellState for "well_name" to the
//   new ALQ value and related data (the data has already been computed and
//   saved in "grad_map")
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
addOrRemoveALQincrement_(GradMap &grad_map, const std::string well_name, bool add)
{
    // only applies to wells in the well_state_map (i.e. wells on this rank)
    auto it = this->well_state_map_.find(well_name);
    if (it == this->well_state_map_.end())
        return;

    GLiftWellState &state = *(it->second.get());
    const GradInfo &gi = grad_map.at(well_name);
    if (this->debug_) {
        auto new_alq = gi.alq;
        auto old_alq = state.alq();
        const std::string msg = fmt::format("well {} : {} ALQ increment, "
            "old alq: {}, new alq: {}",
            well_name, (add ? "adding" : "subtracting"), old_alq, new_alq);
        this->displayDebugMessage_(msg);
    }
    state.update(gi.new_oil_rate, gi.oil_is_limited,
        gi.new_gas_rate, gi.gas_is_limited,
        gi.alq, gi.alq_is_limited, add);
    this->well_state_.setALQ(well_name, gi.alq);
}

template<typename TypeTag>
std::optional<typename GasLiftStage2<TypeTag>::GradInfo>
GasLiftStage2<TypeTag>::
calcIncOrDecGrad_(
    const std::string well_name, const GasLiftSingleWell &gs_well, bool increase)
{

    // only applies to wells in the well_state_map (i.e. wells on this rank)
    if(this->well_state_map_.count(well_name) == 0)
        return std::nullopt;

    GLiftWellState &state = *(this->well_state_map_.at(well_name).get());
    if (checkRateAlreadyLimited_(state, increase)) {
        /*
        const std::string msg = fmt::format(
            "well {} : not able to obtain {} gradient",
            well_name,
            (increase ? "incremental" : "decremental")
        );
        displayDebugMessage_(msg);
        */
        return std::nullopt;
    }
    else {
        auto [oil_rate, gas_rate] = state.getRates();
        auto alq = state.alq();
        auto grad = gs_well.calcIncOrDecGradient(oil_rate, gas_rate, alq, increase);
        if (grad) {
            const std::string msg = fmt::format(
              "well {} : adding {} gradient = {}",
              well_name,
              (increase ? "incremental" : "decremental"),
              grad->grad
            );
            displayDebugMessage_(msg);
        }
        return grad;
    }
}

template<typename TypeTag>
bool
GasLiftStage2<TypeTag>::
checkRateAlreadyLimited_(GLiftWellState &state, bool increase)
{
    auto current_increase = state.increase();
    bool do_check = false;
    if (current_increase) {
        if (*current_increase == increase) do_check = true;
    }
    else {
        // If current_increase is not defined, it means that stage1
        //   was unable to either increase nor decrease the ALQ. If the
        //   initial rates stored in "state" is limited, and if
        //   "increase" is true, it is not likely that adding ALQ will
        //   cause the new rates not to be limited. However, if
        //   "increase" is false, subtracting ALQ can make the new rates
        //   not limited.
        if (increase) do_check = true;
    }
    if (do_check) {
        if (state.gasIsLimited() || state.oilIsLimited() || state.alqIsLimited()) {
            const std::string msg = fmt::format(
                "{} gradient : skipping since {} was limited in previous step",
                (increase ? "incremental" : "decremental"),
                (state.oilIsLimited() ? "oil" :
                    (state.gasIsLimited() ? "gas" : "alq")));
            displayDebugMessage_(msg);
            return true;
        }
    }
    return false;
}


template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
deleteGrad_(const std::string &name, bool increase)
{
    GradMap &map = increase ? this->inc_grads_ : this->dec_grads_;
    auto value = map.at(name);
    map.erase(name);
    return value;
}

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
deleteDecGradItem_(const std::string &name)
{
    return deleteGrad_(name, /*increase=*/false);
}

template<typename TypeTag>
typename GasLiftStage2<TypeTag>::GradInfo
GasLiftStage2<TypeTag>::
deleteIncGradItem_(const std::string &name)
{
    return deleteGrad_(name, /*increase=*/true);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayWarning_(const std::string &msg, const std::string &group_name)
{
    const std::string message = fmt::format(
        "GAS LIFT OPTIMIZATION (STAGE2), GROUP: {} : {}", group_name, msg);
    this->deferred_logger_.warning("WARNING", message);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayWarning_(const std::string &msg)
{
    const std::string message = fmt::format(
        "GAS LIFT OPTIMIZATION (STAGE2) : {}", msg);
    this->deferred_logger_.warning("WARNING", message);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayDebugMessage_(const std::string &msg)
{
    if (this->debug_) {
        const std::string message = fmt::format(
            "  GLIFT2 (DEBUG) : {}", msg);
        this->deferred_logger_.info(message);
    }
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayDebugMessage2B_(const std::string &msg)
{
    if (this->debug_) {
        const std::string message = fmt::format(
            "Stage 2B : {}", msg);
        displayDebugMessage_(message);
    }
}
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
displayDebugMessage_(const std::string &msg, const std::string &group_name)
{
    if (this->debug_) {
        const std::string message = fmt::format(
            "Group {} : {}", group_name, msg);
        displayDebugMessage_(message);
    }
}

template<typename TypeTag>
std::tuple<double, double, double>
GasLiftStage2<TypeTag>::
getCurrentGroupRates_(const Group &group)
{
    auto rates = getCurrentGroupRatesRecursive_(group);
    const auto& comm = ebos_simulator_.vanguard().grid().comm();
    comm.sum(rates.data(), rates.size());
    auto [oil_rate, gas_rate, alq] = rates;
    if (this->debug_) {
        const std::string msg = fmt::format(
            "Current group rates for {} : oil: {}, gas: {}, alq: {}",
            group.name(), oil_rate, gas_rate, alq);
        displayDebugMessage2B_(msg);
    }

    return {oil_rate, gas_rate, alq};
}


template<typename TypeTag>
std::array <double, 3> GasLiftStage2<TypeTag>::
getCurrentGroupRatesRecursive_(const Group &group)
{
    double oil_rate = 0.0;
    double gas_rate = 0.0;
    double alq = 0.0;
    // NOTE: A group can either contain wells or groups, but not both
    if (group.wellgroup()) {
        for (const std::string& well_name : group.wells()) {
            auto [sw_oil_rate, sw_gas_rate, sw_alq] =
                getCurrentWellRates_(well_name, group.name());
            oil_rate += sw_oil_rate;
            gas_rate += sw_gas_rate;
            alq += sw_alq;
        }

    }
    else {
        for (const std::string& group_name : group.groups()) {
            if(this->schedule_.back().groups.has(group_name)) {
                const Group& sub_group =
                    this->schedule_.getGroup(group_name, this->report_step_idx_);
                // If groups have efficiency factors to model
                // synchronized downtime of their subordinate wells
                // (see keyword GEFAC), their lift gas injection rates
                // are multiplied by their efficiency factors when
                // they are added to the lift gas supply rate of the
                // parent group.
                const auto gefac = sub_group.getGroupEfficiencyFactor();
                auto rates = getCurrentGroupRatesRecursive_(sub_group);
                auto [sg_oil_rate, sg_gas_rate, sg_alq] = rates;
                oil_rate += (gefac * sg_oil_rate);
                gas_rate += (gefac * sg_gas_rate);
                alq += (gefac * sg_alq);
            }
        }
    }
    return {oil_rate, gas_rate, alq};
}

template<typename TypeTag>
std::tuple<double, double, double>
GasLiftStage2<TypeTag>::
getCurrentWellRates_(const std::string &well_name, const std::string &group_name)
{
    double oil_rate, gas_rate, alq;
    bool success = false;
    const WellInterface<TypeTag> *well_ptr = nullptr;
    std::string debug_info;
    if (this->stage1_wells_.count(well_name) == 1) {
        GasLiftSingleWell &gs_well = *(this->stage1_wells_.at(well_name).get());
        const WellInterface<TypeTag> &well = gs_well.getStdWell();
        well_ptr = &well;
        GLiftWellState &state = *(this->well_state_map_.at(well_name).get());
        std::tie(oil_rate, gas_rate) = state.getRates();
        success = true;
        if ( this->debug_) debug_info = "(A)";
    }
    else if (this->prod_wells_.count(well_name) == 1) {
        well_ptr = this->prod_wells_.at(well_name);
        std::tie(oil_rate, gas_rate) = getStdWellRates_(*well_ptr);
        success = true;
        if ( this->debug_) debug_info = "(B)";
    }

    if (well_ptr) {
        // we only want rates from wells owned by the rank
        if (!well_state_.wellIsOwned(well_ptr->indexOfWell(), well_name)) {
            success = false;
        }
    }

    if (success) {
        assert(well_ptr);
        assert(well_ptr->isProducer());
        alq = this->well_state_.getALQ(well_name);
        if (this->debug_) {
            const std::string msg = fmt::format(
                "Rates {} for well {} : oil: {}, gas: {}, alq: {}",
                debug_info, well_name, oil_rate, gas_rate, alq);
            displayDebugMessage_(msg, group_name);
        }
        // If wells have efficiency factors to take account of regular
        // downtime (see keyword WEFAC), their lift gas injection
        // rates are multiplied by their efficiency factors when they
        // are added to the group lift gas supply rate. This is
        // consistent with the summation of flow rates for wells with
        // downtime, and preserves the ratio of production rate to
        // lift gas injection rate.
        const auto &well_ecl = well_ptr->wellEcl();
        double factor = well_ecl.getEfficiencyFactor();
        oil_rate *= factor;
        gas_rate *= factor;
        alq *= factor;
        if (this->debug_ && (factor != 1)) {
            const std::string msg = fmt::format(
                "Well {} : efficiency factor {}. New rates : oil: {}, gas: {}, alq: {}",
                well_name, factor, oil_rate, gas_rate, alq);
            displayDebugMessage_(msg, group_name);
        }
    }
    else {
        // NOTE: This happens for wells that are not producers, or not active.
        if (this->debug_) {
            const std::string msg = fmt::format("Could not determine current rates for "
                "well {}: (not active or injector)", well_name);
            displayDebugMessage_(msg, group_name);
        }
        oil_rate = 0.0; gas_rate = 0.0; alq = 0.0;
    }
    return std::make_tuple(oil_rate, gas_rate, alq);
}

template<typename TypeTag>
std::pair<double, double>
GasLiftStage2<TypeTag>::
getStdWellRates_(const WellInterface<TypeTag> &well)
{
    const int well_index = well.indexOfWell();
    const int np = this->well_state_.numPhases();
    const auto& pu = well.phaseUsage();
    auto oil_rate =
        -this->well_state_.wellRates()[np * well_index + pu.phase_pos[Oil]];
    auto gas_rate =
        -this->well_state_.wellRates()[np * well_index + pu.phase_pos[Gas]];
    return {oil_rate, gas_rate};
}

// Find all subordinate wells of a given group.
//
// NOTE: A group can either contain wells or groups, not both.
//   If it contains groups, we have to traverse those recursively to find the wells.
//
// NOTE: This means that wells are located at the leaf nodes of the tree, and
//       groups are located at the other nodes (not leaf nodes) of the tree
//
template<typename TypeTag>
std::vector<GasLiftSingleWell<TypeTag> *>
GasLiftStage2<TypeTag>::
getGroupGliftWells_(const Group &group)
{
    std::vector<GasLiftSingleWell *> wells;
    getGroupGliftWellsRecursive_(group, wells);
    return wells;
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
getGroupGliftWellsRecursive_(const Group &group,
         std::vector<GasLiftSingleWell *> &wells)
{
    for (const std::string& group_name : group.groups()) {
    if(this->schedule_.back().groups.has(group_name)) {
            const Group& sub_group =
                this->schedule_.getGroup(group_name, this->report_step_idx_);
            getGroupGliftWellsRecursive_(sub_group, wells);
        }
    }
    for (const std::string& well_name : group.wells()) {
        if (this->stage1_wells_.count(well_name) == 1) {
            GasLiftSingleWell *well_ptr = this->stage1_wells_.at(well_name).get();
            wells.push_back(well_ptr);
        }
    }
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
optimizeGroup_(const Group &group)
{
    const auto &gl_group = this->glo_.group(group.name());
    const auto &max_glift = gl_group.max_lift_gas();
    const auto &max_total_gas = gl_group.max_total_gas();
    if (group.has_control(Group::ProductionCMode::ORAT)
                       || max_glift || max_total_gas)
    {
        displayDebugMessage_("optimizing", group.name());
        auto wells = getGroupGliftWells_(group);
        std::vector<GradPair> inc_grads;
        std::vector<GradPair> dec_grads;
        redistributeALQ_(wells, group, inc_grads, dec_grads);
        removeSurplusALQ_(group, inc_grads, dec_grads);
    }
    else {
        displayDebugMessage_("skipping", group.name());
    }
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
optimizeGroupsRecursive_(const Group &group)
{
    for (const std::string& group_name : group.groups()) {
        if(!this->schedule_.back().groups.has(group_name))
            continue;
        const Group& sub_group = this->schedule_.getGroup(
            group_name, this->report_step_idx_);
        optimizeGroupsRecursive_(sub_group);
    }
    // TODO: should we also optimize groups that do not have GLIFTOPT defined?
    //   (i.e. glo_.has_group(name) returns false)
    //   IF GLIFTOPT is not specified for the group or if item 2 of GLIFTOPT
    //   is defaulted, there is no maximum lift gas supply for the group.
    //   But even if there is no limit on the liftgas supply it can still
    //   be desireable to use as little ALQ as possible to achieve a
    //   group oil rate limit or gas rate limit.
    if (this->glo_.has_group(group.name())) // only optimize if GLIFTOPT is given
        optimizeGroup_(group);

}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
recalculateGradientAndUpdateData_(
    GradPairItr &grad_itr, bool increase,
    std::vector<GradPair> &grads, std::vector<GradPair> &other_grads)
{
    // NOTE: We make a copy of the name string instead of taking a reference
    //   since we may have to erase grad_itr (in the "else" condition below)
    const std::string name = grad_itr->first;
    std::optional<GradInfo> old_grad = std::nullopt;

    // only applies to wells in the well_state_map (i.e. wells on this rank)
    // the grads and other grads are synchronized later
    if(this->stage1_wells_.count(name) > 0) {
        GasLiftSingleWell &gs_well = *(this->stage1_wells_.at(name).get());
        auto grad = calcIncOrDecGrad_(name, gs_well, increase);
        if (grad) {
            grad_itr->second = grad->grad;
            old_grad = updateGrad_(name, *grad, increase);
        }
        else {
            grads.erase(grad_itr); // NOTE: this invalidates grad_itr
            old_grad = deleteGrad_(name, increase);
        }
    }

    if (old_grad) {
        // NOTE: Either creates a new item or reassigns
        // The old incremental gradient becomes the new decremental gradient
        //   or the old decremental gradient becomes the new incremental gradient
        updateGrad_(name, *old_grad, !increase);
        updateGradVector_(name, other_grads, old_grad->grad);
    }
}

// redistributeALQ_() :
//
// DESCRIPTION: The currently allocated increments are redistributed by
// comparing the largest weighted incremental gradient and the
// smallest weighted decremental gradient of the subordinate
// wells. Consider, for example, that the largest weighted incremental
// gradient belongs to well W1, and exceeds the smallest weighted
// decremental gradient, which belongs to well W2. It would be more
// profitable to subtract an increment of lift gas from well W2 and
// allocate it to well W1. After the exchange, the incremental and
// decremental gradients of wells W1 and W2 are recalculated. The
// exchange of increments continues until no weighted incremental
// gradient exceeds a weighted decremental gradient.
//
// TODO: maybe the redistribution can be simplified by only doing it for the
//   FIELD group: For example:
//
//
//
//                                       FIELD
//                                         |
//                                       PLAT-A
//                          ---------------+---------------------
//                         |                                    |
//                        M5S                                  M5N
//                ---------+----------                     -----+-------
//               |                   |                    |            |
//              B1                  G1                   C1           F1
//           ----+------          ---+---              ---+---       ---+---
//          |    |     |         |      |             |      |      |      |
//        B-1H  B-2H  B-3H     G-3H    G-4H         C-1H   C-2H    F-1H   F-2H
//
//
//  it is probably unecessary to first redistribute ALQ for the wells B-1H, B-2H,
//  and B-3H in group B1, and then in a next separate step, redistribute ALQ for the
//  wells G-3H, and G-4H, and similarly, for the wells in groups C1, and F1,
//  and then, for the wells under M5S, and then M5N, and finally repeat the procedure
//  for all the wells under group PLAT-A, i.e. the wells:
//    B-1H, B-2H, B-3H, G-3H, G-4H, C-1H, C-2H, F-1H, and F-2H.
//  It seems like it would be more efficient to
//  just do it once for the topmost group "PLAT-A" and then skip redistribution for
//  all sub groups of "PLAT-A"
//
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
redistributeALQ_(std::vector<GasLiftSingleWell *> &wells,  const Group &group,
    std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    OptimizeState state {*this, group};
    std::vector<GradPair> inc_grads_local;
    std::vector<GradPair> dec_grads_local;
    inc_grads_local.reserve(wells.size());
    dec_grads_local.reserve(wells.size());
    state.calculateEcoGradients(wells, inc_grads_local, dec_grads_local);

    // the gradients needs to be communicated to all ranks
    dec_grads = localToGlobalGradVector_(dec_grads_local);
    inc_grads = localToGlobalGradVector_(inc_grads_local);

    if (!state.checkAtLeastTwoWells(wells)) {
        // NOTE: Even though we here in redistributeALQ_() do not use the
        //   economic gradient if there is only a single well, we still
        //   need to calculate it since inc_grads and dec_grads are returned
        //   and will be used by removeSurplusALQ_() later.
        return;
    }
    bool stop_iteration = false;
    while (!stop_iteration && (state.it++ <= this->max_iterations_)) {
        state.debugShowIterationInfo();
        auto [min_dec_grad, max_inc_grad]
            = state.getEcoGradients(inc_grads, dec_grads);
        if (min_dec_grad) {
            assert( max_inc_grad );
            // Redistribute if the largest incremental gradient exceeds the
            //   smallest decremental gradient
            if ((*max_inc_grad)->second > (*min_dec_grad)->second) {
                state.redistributeALQ(*min_dec_grad, *max_inc_grad);
                state.recalculateGradients(
                    inc_grads, dec_grads, *min_dec_grad, *max_inc_grad);
                continue;
            }
        }
        stop_iteration = true;
    }
    if (state.it > this->max_iterations_) {
        const std::string msg = fmt::format("Max iterations {} exceeded.",
            this->max_iterations_);
        displayWarning_(msg, group.name());
    }
}

// The group has surplus lift gas if it exceeds any production rate limits
//   or a lift gas supply limit, or contains any wells that have a weighted
//   decremental gradient less than the minimum economic gradient.
// Lift gas increments are removed in turn from the well that currently has
//   the smallest weighted decremental gradient, until there is no surplus
//   lift gas in the group.
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
removeSurplusALQ_(const Group &group,
    std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    if (dec_grads.size() == 0) {
        displayDebugMessage2B_("no wells to remove ALQ from. Skipping");
        return;
    }
    assert(dec_grads.size() > 0);
    const auto &gl_group = this->glo_.group(group.name());
    const auto &max_glift = gl_group.max_lift_gas();
    const auto controls = group.productionControls(this->summary_state_);
    //const auto &max_total_gas = gl_group.max_total_gas();
    auto [oil_rate, gas_rate, alq] = getCurrentGroupRates_(group);
    auto min_eco_grad = this->glo_.min_eco_gradient();
    bool stop_iteration = false;
    if (this->debug_) {
        std::string max_glift_str = "unlimited";
        if (max_glift) max_glift_str = fmt::format("{}", *max_glift);
        const std::string msg = fmt::format("Starting iteration for group: {}. "
            "oil_rate = {}, oil_target = {}, gas_rate = {}, gas_target = {}, "
            "alq = {}, max_alq = {}", group.name(), oil_rate, controls.oil_target,
            gas_rate, controls.gas_target, alq, max_glift_str);
        displayDebugMessage2B_(msg);
    }
    SurplusState state {*this, group, oil_rate, gas_rate, alq,
            min_eco_grad, controls.oil_target, controls.gas_target, max_glift };

    while (!stop_iteration) {
        if (dec_grads.size() >= 2) {
            sortGradients_(dec_grads);
        }
        auto dec_grad_itr = dec_grads.begin();
        const auto well_name = dec_grad_itr->first;
        auto eco_grad = dec_grad_itr->second;
        bool remove = false;
        if (state.checkOilTarget() || state.checkGasTarget() || state.checkALQlimit()) {
            remove = true;
        }
        else {
            // NOTE: It is enough to check the economic gradient of the first well
            //   in dec_grads since they are sorted according to the eco. grad.
            //   If the first well's eco. grad. is greater than the minimum
            //   eco. grad. then all the other wells' eco. grad. will also be
            //   greater.
            if (state.checkEcoGradient(well_name, eco_grad)) remove = true;
        }
        if (remove) {
            state.updateRates(well_name);
            state.addOrRemoveALQincrement( this->dec_grads_, well_name, /*add=*/false);
            recalculateGradientAndUpdateData_(
                        dec_grad_itr, /*increase=*/false, dec_grads, inc_grads);

            // The dec_grads and inc_grads needs to be syncronized across ranks
            dec_grads = updateGlobalGradVector_(dec_grads);
            inc_grads = updateGlobalGradVector_(inc_grads);
            // NOTE: recalculateGradientAndUpdateData_() will remove the current gradient
            //   from dec_grads if it cannot calculate a new decremental gradient.
            //   This will invalidate dec_grad_itr and well_name
            if (dec_grads.size() == 0) stop_iteration = true;
            ++state.it;
        }
        else {
            stop_iteration = true;
        }
    }
    if (state.it >= 1) {
        if (this->debug_) {
            auto [oil_rate2, gas_rate2, alq2] = getCurrentGroupRates_(group);
            const std::string msg = fmt::format(
                 "Finished after {} iterations for group: {}."
                 " oil_rate = {}, gas_rate = {}, alq = {}", state.it,
                 group.name(), oil_rate2, gas_rate2, alq2);
            displayDebugMessage2B_(msg);
        }
    }
    else {
        displayDebugMessage2B_("Finished after 0 iterations");
    }
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
saveGrad_(GradMap &map, const std::string &name, GradInfo &grad)
{
    if (auto it = map.find(name); it == map.end()) {
        [[maybe_unused]] auto result = map.emplace(name, grad);
        assert(result.second); // the insert was successful
    }
    else {
        it->second = grad;
    }
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
saveDecGrad_(const std::string &name, GradInfo &grad)
{
    saveGrad_(this->dec_grads_, name, grad);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
saveIncGrad_(const std::string &name, GradInfo &grad)
{
    saveGrad_(this->inc_grads_, name, grad);
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
sortGradients_(std::vector<GradPair> &grads)
{
    auto cmp = [](GradPair a, GradPair b) {
         return a.second <  b.second;
    };
    std::sort(grads.begin(), grads.end(), cmp);
}

template<typename TypeTag>
std::optional<typename GasLiftStage2<TypeTag>::GradInfo>
GasLiftStage2<TypeTag>::
updateGrad_(const std::string &name, GradInfo &grad, bool increase)
{
    GradMap &map = increase ? this->inc_grads_ : this->dec_grads_;
    std::optional<GradInfo> old_value = std::nullopt;
    if (map.count(name) == 1) {
        old_value = map.at(name);
    }
    map[name] = grad;
    return old_value;
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::
updateGradVector_(const std::string &name, std::vector<GradPair> &grads, double grad)
{
    for (auto itr = grads.begin(); itr != grads.end(); itr++) {
        if (itr->first == name) {
            itr->second = grad;
            return;
        }
    }
    grads.push_back({name, grad});
    // NOTE: the gradient vector is no longer sorted, but sorting will be done
    //   later in getEcoGradients()
}

template<typename TypeTag>
std::vector<typename GasLiftStage2<TypeTag>::GradPair>
GasLiftStage2<TypeTag>::
localToGlobalGradVector_(const std::vector<GradPair> &grads_local) const
{
    const auto& comm = ebos_simulator_.vanguard().grid().comm();
    if (comm.size() == 1)
        return grads_local;

    using Pair = std::pair<int, double>;
    std::vector<Pair> grads_local_tmp;
    grads_local_tmp.reserve(grads_local.size());
    for (size_t i = 0; i < grads_local.size(); ++i) {
        if(!well_state_.wellIsOwned(grads_local[i].first))
            continue;

        grads_local_tmp.push_back(std::make_pair(well_state_.wellNameToGlobalIdx(grads_local[i].first), grads_local[i].second));
    }

    std::vector<int> sizes_(comm.size());
    std::vector<int> displ_(comm.size() + 1, 0);
    int mySize = grads_local_tmp.size();
    comm.allgather(&mySize, 1, sizes_.data());
    std::partial_sum(sizes_.begin(), sizes_.end(), displ_.begin()+1);
    std::vector<Pair> grads_global(displ_.back());

    comm.allgatherv(grads_local_tmp.data(), grads_local_tmp.size(), grads_global.data(), sizes_.data(), displ_.data());
    std::vector<GradPair> grads(grads_global.size());
    for (size_t i = 0; i < grads_global.size(); ++i) {
        grads[i] = std::make_pair(well_state_.globalIdxToWellName(grads_global[i].first), grads_global[i].second);
    }
    return grads;
}

template<typename TypeTag>
std::vector<typename GasLiftStage2<TypeTag>::GradPair>
GasLiftStage2<TypeTag>::
updateGlobalGradVector_(const std::vector<GradPair> &grads_global) const
{
    const auto& comm = ebos_simulator_.vanguard().grid().comm();
    if (comm.size() == 1)
        return grads_global;

    std::vector<GradPair> grads_local;
    for (auto itr = grads_global.begin(); itr != grads_global.end(); itr++) {
        if (well_state_map_.count(itr->first) > 0) {
            grads_local.push_back(*itr);
        }
    }
    return localToGlobalGradVector_(grads_local);
}

/***********************************************
 * Public methods declared in OptimizeState
 ***********************************************/

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
calculateEcoGradients(std::vector<GasLiftSingleWell *> &wells,
           std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    for (auto well_ptr : wells) {
        const auto &gs_well = *well_ptr;  // gs = GasLiftSingleWell
        const auto &name = gs_well.name();
        auto inc_grad = this->parent.calcIncOrDecGrad_(name, gs_well, /*increase=*/true);
        if (inc_grad) {
            inc_grads.emplace_back(std::make_pair(name, inc_grad->grad));
            this->parent.saveIncGrad_(name, *inc_grad);
        }
        auto dec_grad = this->parent.calcIncOrDecGrad_(name, gs_well, /*increase=*/false);
        if (dec_grad) {
            dec_grads.emplace_back(std::make_pair(name, dec_grad->grad));
            this->parent.saveDecGrad_(name, *dec_grad);
        }
    }
}


template<typename TypeTag>
bool
GasLiftStage2<TypeTag>::OptimizeState::
checkAtLeastTwoWells(std::vector<GasLiftSingleWell *> &wells)
{
    int numberOfwells = 0;
    for (auto well : wells){
        int index_of_wells = well->getStdWell().indexOfWell();
        if (!this->parent.well_state_.wellIsOwned(index_of_wells, well->name()))
            continue;
        numberOfwells++;
    }
    const auto& comm = this->parent.ebos_simulator_.vanguard().grid().comm();
    numberOfwells = comm.sum(numberOfwells);
    if (numberOfwells < 2) {
        const std::string msg = fmt::format(
            "skipping: too few wells ({}) to optimize.", numberOfwells);
        displayDebugMessage_(msg);
        return false;
    }
    return true;
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
debugShowIterationInfo()
{
    const std::string msg = fmt::format("iteration {}", this->it);
    displayDebugMessage_(msg);
}

template<typename TypeTag>
std::pair<std::optional<typename GasLiftStage2<TypeTag>::GradPairItr>,
          std::optional<typename GasLiftStage2<TypeTag>::GradPairItr>>
GasLiftStage2<TypeTag>::OptimizeState::
getEcoGradients(std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    if (inc_grads.size() > 0 && dec_grads.size() > 0) {
        this->parent.sortGradients_(inc_grads);
        this->parent.sortGradients_(dec_grads);
        // The largest incremental gradient is the last element
        auto inc_grad = std::prev(inc_grads.end());
        std::optional<GradPairItr> inc_grad_opt;
        std::optional<GradPairItr> dec_grad_opt;
        // The smallest decremental gradient is at the beginning
        for (auto itr = dec_grads.begin(); itr != dec_grads.end(); itr++) {
            if (itr->first == inc_grad->first) {
                // Don't consider decremental gradients with the same well name
                continue;
            }
            dec_grad_opt = itr;
            break;
        }
        if (dec_grad_opt) {
            inc_grad_opt = inc_grad;
            return { dec_grad_opt, inc_grad_opt };
        }
    }
    return {std::nullopt, std::nullopt};
}

// Recalculate gradients (and related information, see struct GradInfo in
//   GasLiftSingleWell.hpp) after an ALQ increment
//   has been given from the well with minumum decremental gradient (represented
//   by the input argument min_dec_grad_itr) to the well with the largest
//   incremental gradient (represented by input argument max_inc_grad_itr).
//
// For the well with the largest incremental gradient, we compute a new
//   incremental gradient given the new ALQ. The new decremental gradient for this
//   well is set equal to the current incremental gradient (before the ALQ is added)
// Similiarly, for the well with the smalles decremental gradient, we compute
//   a new decremental gradient given the new ALQ. The new incremental gradient
//   for this well is set equal to the current decremental gradient
//   (before the ALQ is subtracted)
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
recalculateGradients(
         std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads,
         GradPairItr &min_dec_grad_itr, GradPairItr &max_inc_grad_itr)
{
    this->parent.recalculateGradientAndUpdateData_(
        max_inc_grad_itr, /*increase=*/true, inc_grads, dec_grads);
    this->parent.recalculateGradientAndUpdateData_(
        min_dec_grad_itr, /*increase=*/false, dec_grads, inc_grads);

    // The dec_grads and inc_grads needs to be syncronized across ranks
    dec_grads = this->parent.updateGlobalGradVector_(dec_grads);
    inc_grads = this->parent.updateGlobalGradVector_(inc_grads);
}

// Take one ALQ increment from well1, and give it to well2
template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
    redistributeALQ( GradPairItr &min_dec_grad, GradPairItr &max_inc_grad)
{
    const std::string msg = fmt::format(
        "redistributing ALQ from well {} (dec gradient: {}) "
        "to well {} (inc gradient {})",
        min_dec_grad->first, min_dec_grad->second,
        max_inc_grad->first, max_inc_grad->second);
    displayDebugMessage_(msg);
    this->parent.addOrRemoveALQincrement_(
        this->parent.dec_grads_, /*well_name=*/min_dec_grad->first, /*add=*/false);
    this->parent.addOrRemoveALQincrement_(
        this->parent.inc_grads_, /*well_name=*/max_inc_grad->first, /*add=*/true);
}

/**********************************************
 * Private methods declared in OptimizeState
 **********************************************/

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
displayDebugMessage_(const std::string &msg)
{
    this->parent.displayDebugMessage_(msg, this->group.name());
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::OptimizeState::
displayWarning_(const std::string &msg)
{
    this->parent.displayWarning_(msg, this->group.name());
}

/**********************************************
 * Public methods declared in SurplusState
 **********************************************/

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::SurplusState::
addOrRemoveALQincrement(GradMap &grad_map, const std::string well_name, bool add)
{
    if (this->parent.debug_) {
        const std::string msg = fmt::format("group: {} : well {} : {} ALQ increment",
            this->group.name(), well_name, (add ? "adding" : "subtracting"));
        this->parent.displayDebugMessage2B_(msg);
    }
    this->parent.addOrRemoveALQincrement_(grad_map, well_name, add);
}

template<typename TypeTag>
bool
GasLiftStage2<TypeTag>::SurplusState::
checkALQlimit()
{
    if (this->max_glift) {
        double max_alq = *(this->max_glift);
        double increment = this->parent.glo_.gaslift_increment();
        double epsilon = 1e-6 * increment;
        if ((max_alq+epsilon) < this->alq  ) {
            if (this->parent.debug_) {
                const std::string msg = fmt::format("group: {} : "
                    "ALQ rate {} is greater than ALQ limit {}", this->group.name(),
                    this->alq, max_alq);
                this->parent.displayDebugMessage2B_(msg);
            }
            return true;
        }
    }
    return false;
}

template<typename TypeTag>
bool
GasLiftStage2<TypeTag>::SurplusState::
checkEcoGradient(const std::string &well_name, double eco_grad)
{
    if (eco_grad < this->min_eco_grad) {
        if (this->parent.debug_) {
            const std::string msg = fmt::format("group: {}, well: {} : "
                "economic gradient {} less than minimum ({})", this->group.name(),
                well_name, eco_grad, this->min_eco_grad);
            this->parent.displayDebugMessage2B_(msg);
        }
        return true;
    }
    else {
        return false;
    }
}

template<typename TypeTag>
bool
GasLiftStage2<TypeTag>::SurplusState::
checkGasTarget()
{
    if (this->group.has_control(Group::ProductionCMode::GRAT)) {
        if (this->gas_target < this->gas_rate  ) {
            if (this->parent.debug_) {
                const std::string msg = fmt::format("group: {} : "
                    "gas rate {} is greater than gas target {}", this->group.name(),
                    this->gas_rate, this->gas_target);
                this->parent.displayDebugMessage2B_(msg);
            }
            return true;
        }
    }
    return false;
}

template<typename TypeTag>
bool
GasLiftStage2<TypeTag>::SurplusState::
checkOilTarget()
{
    if (this->group.has_control(Group::ProductionCMode::ORAT)) {
        if (this->oil_target < this->oil_rate  ) {
            if (this->parent.debug_) {
                const std::string msg = fmt::format("group: {} : "
                    "oil rate {} is greater than oil target {}", this->group.name(),
                    this->oil_rate, this->oil_target);
                this->parent.displayDebugMessage2B_(msg);
            }
            return true;
        }
    }
    return false;
}

template<typename TypeTag>
void
GasLiftStage2<TypeTag>::SurplusState::
updateRates(const std::string &well_name)
{
    std::array<double, 3> delta = {0.0,0.0,0.0};
    // compute the delta on wells on own rank
    if (this->parent.well_state_map_.count(well_name) > 0) {
        const GradInfo &gi = this->parent.dec_grads_.at(well_name);
        GLiftWellState &state = *(this->parent.well_state_map_.at(well_name).get());
        GasLiftSingleWell &gs_well = *(this->parent.stage1_wells_.at(well_name).get());
        const WellInterface<TypeTag> &well = gs_well.getStdWell();
        // only get deltas for wells owned by this rank
        if (this->parent.well_state_.wellIsOwned(well.indexOfWell(), well_name)) {
            const auto &well_ecl = well.wellEcl();
            double factor = well_ecl.getEfficiencyFactor();
            auto& [delta_oil, delta_gas, delta_alq] = delta;
                    delta_oil = factor * (gi.new_oil_rate - state.oilRate());
                    delta_gas = factor * (gi.new_gas_rate - state.gasRate());
                    delta_alq = factor * (gi.alq - state.alq());
        }
    }

    // and communicate the results
    const auto& comm = this->parent.ebos_simulator_.vanguard().grid().comm();
    comm.sum(delta.data(), delta.size());

    // and update
    const auto& [delta_oil, delta_gas, delta_alq] = delta;
    this->oil_rate += delta_oil;
    this->gas_rate += delta_gas;
    this->alq += delta_alq;
}

} // namespace Opm
