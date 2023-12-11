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

#include <config.h>
#include <opm/simulators/wells/GasLiftStage2.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>

#include <fmt/format.h>

#include <cstddef>
#include <optional>
#include <string>

namespace Opm {

GasLiftStage2::GasLiftStage2(
    const int report_step_idx,
    const Parallel::Communication& comm,
    const Schedule& schedule,
    const SummaryState& summary_state,
    DeferredLogger &deferred_logger,
    WellState &well_state,
    const GroupState &group_state,
    GLiftProdWells &prod_wells,
    GLiftOptWells &glift_wells,
    GasLiftGroupInfo& group_info,
    GLiftWellStateMap &state_map,
    bool glift_debug
) :
    GasLiftCommon(well_state, group_state, deferred_logger, comm, glift_debug)
    , prod_wells_{prod_wells}
    , stage1_wells_{glift_wells}
    , group_info_{group_info}
    , well_state_map_{state_map}
    , report_step_idx_{report_step_idx}
    , summary_state_{summary_state}
    , schedule_{schedule}
    , glo_{schedule_.glo(report_step_idx_)}
{
//    this->time_step_idx_
//        = this->ebos_simulator_.model().newtonMethod().currentTimeStep();
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
void
GasLiftStage2::
runOptimize()
{
    const auto& group = this->schedule_.getGroup("FIELD", this->report_step_idx_);

    optimizeGroupsRecursive_(group);

}


/********************************************
 * Protected methods in alphabetical order
 ********************************************/

// Update GasLiftWellState and WellState for "well_name" to the
//   new ALQ value and related data (the data has already been computed and
//   saved in "grad_map")
// INPUT: grad_map : map of incremental (if "add" is true) or decremental
//        (if "add" is false) GradInfo structs for each well name.
void
GasLiftStage2::
addOrRemoveALQincrement_(GradMap &grad_map, const std::string& well_name, bool add)
{
    // only applies to wells in the well_state_map (i.e. wells on this rank)
    auto it = this->well_state_map_.find(well_name);
    if (it == this->well_state_map_.end())
        return;

    GasLiftWellState &state = *(it->second.get());
    const GradInfo &gi = grad_map.at(well_name);
    if (this->debug) {
        auto new_alq = gi.alq;
        auto old_alq = state.alq();
        const std::string msg = fmt::format("well {} : {} ALQ increment, "
            "old alq: {}, new alq: {}",
            well_name, (add ? "adding" : "subtracting"), old_alq, new_alq);
        this->displayDebugMessage_(msg);
    }
    state.update(gi.new_oil_rate, gi.oil_is_limited,
        gi.new_gas_rate, gi.gas_is_limited,
        gi.alq, gi.alq_is_limited, gi.new_water_rate, gi.water_is_limited, add);

    this->well_state_.setALQ(well_name, gi.alq);
    const auto& pu = this->well_state_.phaseUsage();
    std::vector<double> well_pot(pu.num_phases, 0.0);
    if(pu.phase_used[BlackoilPhases::PhaseIndex::Liquid])
        well_pot[pu.phase_pos[BlackoilPhases::PhaseIndex::Liquid]] = gi.new_oil_rate;
    if(pu.phase_used[BlackoilPhases::PhaseIndex::Aqua])
        well_pot[pu.phase_pos[BlackoilPhases::PhaseIndex::Aqua]] = gi.new_water_rate;
    if(pu.phase_used[BlackoilPhases::PhaseIndex::Vapour])
        well_pot[pu.phase_pos[BlackoilPhases::PhaseIndex::Vapour]] = gi.new_gas_rate;

    this->well_state_[well_name].well_potentials = well_pot;
}

std::optional<GasLiftStage2::GradInfo>
GasLiftStage2::
calcIncOrDecGrad_(
    const std::string well_name, const GasLiftSingleWell &gs_well, const std::string& gr_name_dont_limit, bool increase)
{

    // only applies to wells in the well_state_map (i.e. wells on this rank)
    if(this->well_state_map_.count(well_name) == 0)
        return std::nullopt;
    GasLiftWellState &state = *(this->well_state_map_.at(well_name).get());
    if (checkRateAlreadyLimited_(well_name, state, increase)) {
        return std::nullopt;
    }
    else {
        auto [oil_rate, gas_rate] = state.getRates();
        auto alq = state.alq();
        auto grad = gs_well.calcIncOrDecGradient(
            oil_rate, gas_rate, state.waterRate(), alq, gr_name_dont_limit, increase, /*debug_output=*/false);
        if (this->debug) {
            if (grad) {
                const std::string msg = fmt::format(
                    "well {} : alq = {} : adding {} gradient = {}",
                    well_name,
                    alq,
                    (increase ? "incremental" : "decremental"),
                    grad->grad
                );
                displayDebugMessage_(msg);
            }
            else {
                const std::string msg = fmt::format(
                    "well {} : alq = {} : failed to compute {} gradient",
                    well_name,
                    alq,
                    (increase ? "incremental" : "decremental")
                );
                displayDebugMessage_(msg);
            }
        }
        return grad;
    }
}

bool
GasLiftStage2::
checkRateAlreadyLimited_(const std::string& well_name, GasLiftWellState &state, bool increase)
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
        if (state.gasIsLimited() || state.oilIsLimited() || state.alqIsLimited() || state.waterIsLimited()) {
            const std::string msg = fmt::format(
                "Well {} : alq = {} : skipping {} gradient since {} was limited in previous step",
                well_name,
                state.alq(),
                (increase ? "incremental" : "decremental"),
                (state.oilIsLimited() ? "oil" : (state.gasIsLimited() ? "gas" : "alq"))
            );
            displayDebugMessage_(msg);
            return true;
        }
    }
    return false;
}

GasLiftStage2::GradInfo
GasLiftStage2::
deleteGrad_(const std::string &name, bool increase)
{
    GradMap &map = increase ? this->inc_grads_ : this->dec_grads_;
    auto value = map.at(name);
    map.erase(name);
    return value;
}

GasLiftStage2::GradInfo
GasLiftStage2::
deleteDecGradItem_(const std::string &name)
{
    return deleteGrad_(name, /*increase=*/false);
}

GasLiftStage2::GradInfo
GasLiftStage2::
deleteIncGradItem_(const std::string &name)
{
    return deleteGrad_(name, /*increase=*/true);
}

void
GasLiftStage2::
displayWarning_(const std::string &msg, const std::string &group_name)
{
    const std::string message = fmt::format("GROUP: {} : {}", group_name, msg);
    displayWarning_(message);
}

void
GasLiftStage2::
displayWarning_(const std::string &msg)
{
    logMessage_(/*prefix=*/"GLIFT2", msg, MessageType::WARNING);
}

void
GasLiftStage2::
displayDebugMessage_(const std::string &msg) const
{
    if (this->debug) {
        logMessage_(/*prefix=*/"GLIFT2", msg);
    }
}

void
GasLiftStage2::
displayDebugMessage2B_(const std::string &msg)
{
    if (this->debug) {
        logMessage_(/*prefix=*/"GLIFT2B", msg);
    }
}

void
GasLiftStage2::
displayDebugMessage_(const std::string &msg, const std::string &group_name)
{
    if (this->debug) {
        const std::string message = fmt::format(
            "Group {} : {}", group_name, msg);
        displayDebugMessage_(message);
    }
}

std::tuple<double, double, double, double>
GasLiftStage2::
getCurrentGroupRates_(const Group &group)
{
    return {this->group_info_.oilRate(group.name()),
            this->group_info_.gasRate(group.name()),
            this->group_info_.waterRate(group.name()),
            this->group_info_.alqRate(group.name())};
}

std::optional<double>
GasLiftStage2::getGroupMaxALQ_(const Group &group)
{
    if (this->glo_.has_group(group.name())) {
        const auto &gl_group = this->glo_.group(group.name());
        return gl_group.max_lift_gas();
    }
    return std::nullopt; // If GLIFTOPT is missing from schedule, assume unlimited alq
}

std::optional<double>
GasLiftStage2::getGroupMaxTotalGas_(const Group &group)
{
    if (this->glo_.has_group(group.name())) {
        const auto &gl_group = this->glo_.group(group.name());
        return gl_group.max_total_gas();
    }
    return std::nullopt; // If GLIFTOPT is missing from schedule, assume unlimited alq
}

// Find all subordinate wells of a given group.
//
// NOTE: A group can either contain wells or groups, not both.
//   If it contains groups, we have to traverse those recursively to find the wells.
//
// NOTE: This means that wells are located at the leaf nodes of the tree, and
//       groups are located at the other nodes (not leaf nodes) of the tree
//
std::vector<GasLiftSingleWellGeneric*>
GasLiftStage2::
getGroupGliftWells_(const Group &group)
{
    std::vector<GasLiftSingleWell *> wells;
    getGroupGliftWellsRecursive_(group, wells);
    return wells;
}

void
GasLiftStage2::
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

void
GasLiftStage2::
mpiSyncGlobalGradVector_(std::vector<GradPair> &grads_global) const
{
    if (this->comm_.size() == 1)
        return;

    std::vector<GradPair> grads_local;
    for (auto itr = grads_global.begin(); itr != grads_global.end(); itr++) {
        if (this->well_state_map_.count(itr->first) > 0) {
            grads_local.push_back(*itr);
        }
    }
    mpiSyncLocalToGlobalGradVector_(grads_local, grads_global);
}

void
GasLiftStage2::
mpiSyncLocalToGlobalGradVector_(
    const std::vector<GradPair> &grads_local, std::vector<GradPair> &grads_global) const
{
    assert(this->comm_.size() > 1);  // The parent should check if comm. size is > 1
    using Pair = std::pair<int, double>;
    std::vector<Pair> grads_local_tmp;
    grads_local_tmp.reserve(grads_local.size());
    for (std::size_t i = 0; i < grads_local.size(); ++i) {
        if (!this->well_state_.wellIsOwned(grads_local[i].first))
            continue;
        grads_local_tmp.push_back(
           std::make_pair(
              this->well_state_.wellNameToGlobalIdx(grads_local[i].first),
              grads_local[i].second));
    }

    std::vector<int> sizes_(this->comm_.size());
    std::vector<int> displ_(this->comm_.size() + 1, 0);
    int mySize = grads_local_tmp.size();
    this->comm_.allgather(&mySize, 1, sizes_.data());
    std::partial_sum(sizes_.begin(), sizes_.end(), displ_.begin()+1);
    std::vector<Pair> grads_global_tmp(displ_.back());

    this->comm_.allgatherv(grads_local_tmp.data(), grads_local_tmp.size(),
        grads_global_tmp.data(), sizes_.data(), displ_.data());

    // NOTE: This leaves the capacity of 'grads_global' unchanged, so
    //   memory is not reallocated here
    grads_global.clear();

    for (std::size_t i = 0; i < grads_global_tmp.size(); ++i) {
        grads_global.emplace_back(
            std::make_pair(
                well_state_.globalIdxToWellName(grads_global_tmp[i].first),
                grads_global_tmp[i].second));
    }
}

void
GasLiftStage2::
optimizeGroup_(const Group &group)
{
    const auto& group_name = group.name();
    const auto prod_control = this->group_state_.production_control(group_name);
    if ((prod_control != Group::ProductionCMode::NONE) && (prod_control != Group::ProductionCMode::FLD))
    {
        if (this->debug) {
            const std::string msg = fmt::format("optimizing (control = {})", Group::ProductionCMode2String(prod_control));
            displayDebugMessage_(msg, group_name);
        }
        auto wells = getGroupGliftWells_(group);
        std::vector<GradPair> inc_grads;
        std::vector<GradPair> dec_grads;
        redistributeALQ_(wells, group, inc_grads, dec_grads);
        removeSurplusALQ_(group, inc_grads, dec_grads);
    }
    else {
        if (this->debug) {
            const std::string msg = fmt::format("skipping (control = {})", Group::ProductionCMode2String(prod_control));
            displayDebugMessage_(msg, group_name);
        }
    }
}

void
GasLiftStage2::
optimizeGroupsRecursive_(const Group &group)
{
    for (const std::string& group_name : group.groups()) {
        if(!this->schedule_.back().groups.has(group_name))
            continue;
        const Group& sub_group = this->schedule_.getGroup(
            group_name, this->report_step_idx_);
        optimizeGroupsRecursive_(sub_group);
    }
    optimizeGroup_(group);
}

void
GasLiftStage2::
recalculateGradientAndUpdateData_(
        GradPairItr &grad_itr, const std::string& gr_name_dont_limit, bool increase,

        //incremental and decremental gradients, if 'grads' are incremental, then
        // 'other_grads' are decremental, or conversely, if 'grads' are decremental, then
        // 'other_grads' are incremental
        std::vector<GradPair> &grads,  std::vector<GradPair> &other_grads)
{
    // NOTE: We make a copy of the name string instead of taking a reference
    //   since we may have to erase grad_itr (in the "else" condition below)
    const std::string name = grad_itr->first;
    std::optional<GradInfo> old_grad = std::nullopt;

    // only applies to wells in the well_state_map (i.e. wells on this rank)
    // the grads and other grads are synchronized later
    if(this->stage1_wells_.count(name) > 0) {
        GasLiftSingleWell &gs_well = *(this->stage1_wells_.at(name).get());
        {
            auto grad = calcIncOrDecGrad_(name, gs_well, gr_name_dont_limit, increase);
            if (grad) {
                grad_itr->second = grad->grad;
                old_grad = updateGrad_(name, *grad, increase);
            }
            else {
                grads.erase(grad_itr); // NOTE: this invalidates grad_itr
                old_grad = deleteGrad_(name, increase);
            }
        }
        // If the old incremental/decremental gradient was defined, it becomes the new
        //   decremental/incremental gradient
        if (old_grad) {
            // NOTE: Either creates a new item or reassigns
            // The old incremental gradient becomes the new decremental gradient
            //   or the old decremental gradient becomes the new incremental gradient
            // NOTE: The gradient itself (old_grad.grad) will be equal to the new decremental
            //   gradient, but the other fields in old_grad cannot be used as they refer to
            //   the incremental gradient. E.g. old_grad.alq is the alq after an increase in alq,
            //   not a decrease, so we need to recalculate the decremental values here..
            auto grad = calcIncOrDecGrad_(name, gs_well, gr_name_dont_limit, !increase);
            if (grad) {
                updateGrad_(name, *grad, !increase);
                // NOTE: This may invalidate any iterator into 'other_grads' since
                //   updateGradVector_() will do a push_back() if 'name' is not found..
                updateGradVector_(name, other_grads, grad->grad);
            }
            else {
                for (auto it = other_grads.begin(); it != other_grads.end(); it++) {
                    if (it->first == name) {
                        other_grads.erase(it);
                        deleteGrad_(name, !increase);
                    }
                }
            }
        }
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
void
GasLiftStage2::
redistributeALQ_(std::vector<GasLiftSingleWell *> &wells,  const Group &group,
    std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    OptimizeState state {*this, group};
    if (this->comm_.size() == 1) {
        // NOTE: 'inc_grads' and 'dec_grads' can never grow larger than wells.size()
        //   By reserving space here, we can ensure that any push_back() on these
        //   will never reallocate memory and invalidate any iterators.
        inc_grads.reserve(wells.size());
        dec_grads.reserve(wells.size());
        state.calculateEcoGradients(wells, inc_grads, dec_grads);
    }
    else {
        auto max_size = this->comm_.sum(wells.size());
        inc_grads.reserve(max_size);
        dec_grads.reserve(max_size);
        std::vector<GradPair> inc_grads_local;
        std::vector<GradPair> dec_grads_local;
        inc_grads_local.reserve(wells.size());
        dec_grads_local.reserve(wells.size());
        state.calculateEcoGradients(wells, inc_grads_local, dec_grads_local);
        // the gradients needs to be communicated to all ranks
        mpiSyncLocalToGlobalGradVector_(dec_grads_local, dec_grads);
        mpiSyncLocalToGlobalGradVector_(inc_grads_local, inc_grads);
    }

    if (!state.checkAtLeastTwoWells(wells)) {
        // NOTE: Even though we here in redistributeALQ_() do not use the
        //   economic gradient if there is only a single well, we still
        //   need to calculate it (see above) since inc_grads and dec_grads are returned
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
void
GasLiftStage2::
removeSurplusALQ_(const Group &group,
    std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    if (dec_grads.empty()) {
        displayDebugMessage_("no wells to remove ALQ from. Skipping");
        return;
    }
    assert(!dec_grads.empty());
    const auto max_glift = getGroupMaxALQ_(group);
    const auto max_totalgas = getGroupMaxTotalGas_(group);
    const auto controls = group.productionControls(this->summary_state_);
    //const auto &max_total_gas = gl_group.max_total_gas();
    auto [oil_rate, gas_rate, water_rate, alq] = getCurrentGroupRates_(group);
    auto min_eco_grad = this->glo_.min_eco_gradient();
    bool stop_iteration = false;
    if (this->debug) {
        std::string max_glift_str = "unlimited";
        if (max_glift) max_glift_str = fmt::format("{}", *max_glift);
        const std::string msg = fmt::format("Starting remove surplus iteration for group: {}. "
            "oil_rate = {}, oil_target = {}, gas_rate = {}, gas_target = {}, "
            "water_rate = {}, liquid_target = {}, alq = {}, max_alq = {}",
            group.name(), oil_rate, controls.oil_target,
            gas_rate, controls.gas_target, water_rate, controls.liquid_target,
            alq, max_glift_str);
        displayDebugMessage_(msg);
    }
    SurplusState state {*this, group, oil_rate, gas_rate, water_rate, alq,
            min_eco_grad, controls.oil_target, controls.gas_target, controls.water_target,
            controls.liquid_target, max_glift, max_totalgas };

    while (!stop_iteration) {
        if (dec_grads.size() >= 2) {
            sortGradients_(dec_grads);
        }
        auto dec_grad_itr = dec_grads.begin();
        const auto well_name = dec_grad_itr->first;
        auto eco_grad = dec_grad_itr->second;
        bool remove = false;
        const auto delta = state.computeDelta(well_name);
        const auto& [delta_oil, delta_gas, delta_water, delta_alq] = delta;
        if (state.checkOilTarget(delta_oil) || state.checkGasTarget(delta_gas)
              || state.checkLiquidTarget(delta_oil + delta_water) || state.checkWaterTarget(delta_water)
              || state.checkALQlimit(delta_alq, delta_gas))
        {
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
            state.updateRates(delta);
            state.addOrRemoveALQincrement( this->dec_grads_, well_name, /*add=*/false);
            // We pass the current group rate in order to avoid limiting the rates
            // and gaslift based on the current group limits. In other words we want to reduce
            // the gasslift as much as possible as long as we are able to produce the group
            // targets
            recalculateGradientAndUpdateData_(
                        dec_grad_itr, group.name(), /*increase=*/false, dec_grads, inc_grads);

            // The dec_grads and inc_grads needs to be syncronized across ranks
            mpiSyncGlobalGradVector_(dec_grads);
            mpiSyncGlobalGradVector_(inc_grads);
            // NOTE: recalculateGradientAndUpdateData_() will remove the current gradient
            //   from dec_grads if it cannot calculate a new decremental gradient.
            //   This will invalidate dec_grad_itr and well_name
            if (dec_grads.empty()) stop_iteration = true;
            ++state.it;
        }
        else {
            stop_iteration = true;
        }
    }
    if (state.it >= 1) {
        if (this->debug) {
            auto [oil_rate2, gas_rate2, water_rate2, alq2] = getCurrentGroupRates_(group);
            const std::string msg = fmt::format(
                 "Finished after {} iterations for group: {}."
                 " oil_rate = {}, gas_rate = {}, water_rate = {}, alq = {}", state.it,
                 group.name(), oil_rate2, gas_rate2, water_rate2, alq2);
            displayDebugMessage_(msg);
        }
    }
    else {
        displayDebugMessage_("Finished after 0 iterations");
    }
}

void
GasLiftStage2::
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

void
GasLiftStage2::
saveDecGrad_(const std::string &name, GradInfo &grad)
{
    saveGrad_(this->dec_grads_, name, grad);
}

void
GasLiftStage2::
saveIncGrad_(const std::string &name, GradInfo &grad)
{
    saveGrad_(this->inc_grads_, name, grad);
}

void
GasLiftStage2::
sortGradients_(std::vector<GradPair> &grads)
{
    auto cmp = [](GradPair a, GradPair b) {
         return a.second <  b.second;
    };
    std::sort(grads.begin(), grads.end(), cmp);
}

std::optional<GasLiftStage2::GradInfo>
GasLiftStage2::
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

void
GasLiftStage2::
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

/***********************************************
 * Public methods declared in OptimizeState
 ***********************************************/

void
GasLiftStage2::OptimizeState::
calculateEcoGradients(std::vector<GasLiftSingleWell *> &wells,
           std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    for (auto well_ptr : wells) {
        const auto &gs_well = *well_ptr;  // gs = GasLiftSingleWell
        const auto &name = gs_well.name();
        auto inc_grad = this->parent.calcIncOrDecGrad_(name, gs_well, group.name(), /*increase=*/true);
        if (inc_grad) {
            inc_grads.emplace_back(std::make_pair(name, inc_grad->grad));
            this->parent.saveIncGrad_(name, *inc_grad);
        }
        auto dec_grad = this->parent.calcIncOrDecGrad_(name, gs_well, group.name(), /*increase=*/false);
        if (dec_grad) {
            dec_grads.emplace_back(std::make_pair(name, dec_grad->grad));
            this->parent.saveDecGrad_(name, *dec_grad);
        }
    }
}


bool
GasLiftStage2::OptimizeState::
checkAtLeastTwoWells(std::vector<GasLiftSingleWell *> &wells)
{
    int numberOfwells = 0;
    for (auto well : wells){
        int index_of_wells = well->getWell().indexOfWell();
        if (!this->parent.well_state_.wellIsOwned(index_of_wells, well->name()))
            continue;
        numberOfwells++;
    }
    numberOfwells = this->parent.comm_.sum(numberOfwells);
    if (numberOfwells < 2) {
        const std::string msg = fmt::format(
            "skipping: too few wells ({}) to optimize.", numberOfwells);
        displayDebugMessage_(msg);
        return false;
    }
    return true;
}

void
GasLiftStage2::OptimizeState::
debugShowIterationInfo()
{
    const std::string msg = fmt::format("redistribute ALQ iteration {}", this->it);
    displayDebugMessage_(msg);
}

std::pair<std::optional<GasLiftStage2::GradPairItr>,
          std::optional<GasLiftStage2::GradPairItr>>
GasLiftStage2::OptimizeState::
getEcoGradients(std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads)
{
    if (!inc_grads.empty() && !dec_grads.empty()) {
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
//   has been given from the well with minimum decremental gradient (represented
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
void
GasLiftStage2::OptimizeState::
recalculateGradients(
         std::vector<GradPair> &inc_grads, std::vector<GradPair> &dec_grads,
         GradPairItr &min_dec_grad_itr, GradPairItr &max_inc_grad_itr)
{
    this->parent.recalculateGradientAndUpdateData_(
        max_inc_grad_itr, this->group.name(), /*increase=*/true, inc_grads, dec_grads);
    this->parent.recalculateGradientAndUpdateData_(
        min_dec_grad_itr, this->group.name(), /*increase=*/false, dec_grads, inc_grads);

    // The dec_grads and inc_grads needs to be syncronized across ranks
    this->parent.mpiSyncGlobalGradVector_(dec_grads);
    this->parent.mpiSyncGlobalGradVector_(inc_grads);
}

// Take one ALQ increment from well1, and give it to well2
void
GasLiftStage2::OptimizeState::
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

void
GasLiftStage2::OptimizeState::
displayDebugMessage_(const std::string &msg)
{
    this->parent.displayDebugMessage_(msg, this->group.name());
}

void
GasLiftStage2::OptimizeState::
displayWarning_(const std::string &msg)
{
    this->parent.displayWarning_(msg, this->group.name());
}

/**********************************************
 * Public methods declared in SurplusState
 **********************************************/

void
GasLiftStage2::SurplusState::
addOrRemoveALQincrement(GradMap &grad_map, const std::string& well_name, bool add)
{
    if (this->parent.debug) {
        const std::string msg = fmt::format("group: {} : well {} : {} ALQ increment",
            this->group.name(), well_name, (add ? "adding" : "subtracting"));
        this->parent.displayDebugMessage_(msg);
    }
    this->parent.addOrRemoveALQincrement_(grad_map, well_name, add);
}

bool
GasLiftStage2::SurplusState::
checkALQlimit(double delta_alq, double delta_gas)
{
    if (this->max_glift) {
        double max_alq = *(this->max_glift);
        if ((max_alq) < (this->alq + delta_alq)  ) {
            if (this->parent.debug) {
                const std::string msg = fmt::format("group: {} : "
                    "ALQ rate {} is greater than ALQ limit {}", this->group.name(),
                    this->alq, max_alq);
                this->parent.displayDebugMessage_(msg);
            }
            return true;
        }
    }
    if (this->max_total_gas) {
        double max_total = *(this->max_total_gas);
        double total_gas_rate = (this->alq + delta_alq + this->gas_rate + delta_gas);
        if ((max_total) < total_gas_rate ) {
            if (this->parent.debug) {
                const std::string msg = fmt::format("group: {} : "
                    "Total gas rate {} is greater than Total gas limit {}", this->group.name(),
                    total_gas_rate, max_total);
                this->parent.displayDebugMessage_(msg);
            }
            return true;
        }
    }
    return false;
}

bool
GasLiftStage2::SurplusState::
checkEcoGradient(const std::string &well_name, double eco_grad)
{
    if (eco_grad < this->min_eco_grad) {
        if (this->parent.debug) {
            const std::string msg = fmt::format("group: {}, well: {} : "
                "economic gradient {} less than minimum ({})", this->group.name(),
                well_name, eco_grad, this->min_eco_grad);
            this->parent.displayDebugMessage_(msg);
        }
        return true;
    }
    else {
        return false;
    }
}

bool
GasLiftStage2::SurplusState::
checkGasTarget(double delta_gas)
{
    if (this->group.has_control(Group::ProductionCMode::GRAT)) {
        if (this->gas_target < (this->gas_rate + delta_gas) ) {
            if (this->parent.debug) {
                const std::string msg = fmt::format("group: {} : "
                    "gas rate {} is greater than gas target {}", this->group.name(),
                    this->gas_rate, this->gas_target);
                this->parent.displayDebugMessage_(msg);
            }
            return true;
        }
    }
    return false;
}
bool
GasLiftStage2::SurplusState::
checkLiquidTarget(double delta_liquid)
{
    if (this->group.has_control(Group::ProductionCMode::LRAT)) {
        // the change in liquid rate from the is subtracted from the rate to make sure the
        // group still can produce its target
        auto liquid_rate = this->oil_rate + this->water_rate + delta_liquid;
        if (this->liquid_target < (liquid_rate) ) {
            if (this->parent.debug) {
                const std::string msg = fmt::format("group: {} : "
                    "liquid rate {} is greater than liquid target {}", this->group.name(),
                    liquid_rate, this->liquid_target);
                this->parent.displayDebugMessage_(msg);
            }
            return true;
        }
    }
    return false;
}

bool
GasLiftStage2::SurplusState::
checkOilTarget(double delta_oil)
{
    if (this->group.has_control(Group::ProductionCMode::ORAT)) {
        // the change in oil rate from the eco_grad is subtracted from the rate to make sure the
        // group still can produce its target
        if (this->oil_target < (this->oil_rate + delta_oil) ) {
            if (this->parent.debug) {
                const std::string msg = fmt::format("group: {} : "
                    "oil rate {} is greater than oil target {}", this->group.name(),
                    this->oil_rate - delta_oil, this->oil_target);
                this->parent.displayDebugMessage_(msg);
            }
            return true;
        }
    }
    return false;
}

bool
GasLiftStage2::SurplusState::
checkWaterTarget(double delta_water)
{
    if (this->group.has_control(Group::ProductionCMode::WRAT)) {
        if (this->water_target < (this->water_rate + delta_water) ) {
            if (this->parent.debug) {
                const std::string msg = fmt::format("group: {} : "
                    "water rate {} is greater than oil target {}", this->group.name(),
                    this->water_rate, this->water_target);
                this->parent.displayDebugMessage_(msg);
            }
            return true;
        }
    }
    return false;
}

std::array<double, 4>
GasLiftStage2::SurplusState::
computeDelta(const std::string &well_name)
{
    std::array<double, 4> delta = {0.0, 0.0, 0.0, 0.0};
    // compute the delta on wells on own rank
    if (this->parent.well_state_map_.count(well_name) > 0) {
        const GradInfo &gi = this->parent.dec_grads_.at(well_name);
        GasLiftWellState &state = *(this->parent.well_state_map_.at(well_name).get());
        GasLiftSingleWell &gs_well = *(this->parent.stage1_wells_.at(well_name).get());
        const WellInterfaceGeneric &well = gs_well.getWell();
        // only get deltas for wells owned by this rank
        if (this->parent.well_state_.wellIsOwned(well.indexOfWell(), well_name)) {
            const auto &well_ecl = well.wellEcl();
            double factor = well_ecl.getEfficiencyFactor();
            auto& [delta_oil, delta_gas, delta_water, delta_alq] = delta;
            delta_oil = factor * (gi.new_oil_rate - state.oilRate());
            delta_gas = factor * (gi.new_gas_rate - state.gasRate());
            delta_water = factor * (gi.new_water_rate - state.waterRate());
            delta_alq = factor * (gi.alq - state.alq());
        }
    }

    // and communicate the results
    this->parent.comm_.sum(delta.data(), delta.size());

    return delta;
}

void
GasLiftStage2::SurplusState::
updateRates(const std::array<double, 4>& delta)
{
    const auto& [delta_oil, delta_gas, delta_water, delta_alq] = delta;
    this->oil_rate += delta_oil;
    this->gas_rate += delta_gas;
    this->water_rate += delta_water;
    this->alq += delta_alq;
}


} // namespace Opm
