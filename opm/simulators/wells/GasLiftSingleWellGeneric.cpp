/*
  Copyright 2020 Equinor ASA.

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
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>


#include <fmt/format.h>

#include <cassert>
#include <sstream>

namespace Opm
{

GasLiftSingleWellGeneric::GasLiftSingleWellGeneric(
    DeferredLogger& deferred_logger,
    WellState& well_state,
    const GroupState& group_state,
    const Well& ecl_well,
    const SummaryState& summary_state,
    GasLiftGroupInfo &group_info,
    const Schedule& schedule,
    const int report_step_idx,
    GLiftSyncGroups &sync_groups
) :
    deferred_logger_{deferred_logger}
    , well_state_{well_state}
    , group_state_{group_state}
    , ecl_well_{ecl_well}
    , summary_state_{summary_state}
    , group_info_{group_info}
    , sync_groups_{sync_groups}
    , controls_{ecl_well_.productionControls(summary_state_)}
    , num_phases_{well_state_.numPhases()}
    , debug_{false}  // extra debugging output
    , debug_limit_increase_decrease_{false}
{
    this->well_name_ = ecl_well_.name();
    const GasLiftOpt& glo = schedule.glo(report_step_idx);
    // NOTE: According to LIFTOPT, item 1:
    //   "Increment size for lift gas injection rate. Lift gas is
    //   allocated to individual wells in whole numbers of the increment
    //   size. If gas lift optimization is no longer required, it can be
    //   turned off by entering a zero or negative number."
    // NOTE: This condition was checked in doGasLiftOptimize() in StandardWell
    //   so it can be assumed that increment_ > 0
    this->increment_ = glo.gaslift_increment();
    assert( this->increment_ > 0);
    // NOTE: The manual (see LIFTOPT, item 2) does not mention
    //  any default value or restrictions on the economic gradient.
    // TODO: The value of the gradient would most likely be a positive
    //  number. Should we warn or fail on a negative value?
    //  A negative value for the economic gradient would mean that
    //  the oil production is decreasing with increased liftgas
    //  injection (which seems strange)
    this->eco_grad_ = glo.min_eco_gradient();
    gl_well_ = &glo.well(this->well_name_);

}

/****************************************
 * Public methods in alphabetical order
 ****************************************/
// NOTE: Used from GasLiftStage2
std::optional<GasLiftSingleWellGeneric::GradInfo>
GasLiftSingleWellGeneric::
calcIncOrDecGradient(double oil_rate, double gas_rate, double alq, bool increase) const
{
    auto [new_alq_opt, alq_is_limited] = addOrSubtractAlqIncrement_(alq, increase);
    // TODO: What to do if ALQ is limited and new_alq != alq?
    if (!new_alq_opt)
        return std::nullopt;
    double new_alq = *new_alq_opt;
    if (auto bhp = computeBhpAtThpLimit_(new_alq)) {
        auto new_bhp = getBhpWithLimit_(*bhp);
        // TODO: What to do if BHP is limited?
        std::vector<double> potentials(this->num_phases_, 0.0);
        computeWellRates_(new_bhp.first, potentials);
        auto [new_oil_rate, oil_is_limited] = getOilRateWithLimit_(potentials);
        auto [new_gas_rate, gas_is_limited] = getGasRateWithLimit_(potentials);
        if (!increase && new_oil_rate < 0 ) {
            return std::nullopt;
        }
        auto grad = calcEcoGradient_(
            oil_rate, new_oil_rate, gas_rate, new_gas_rate, increase);
        return GradInfo(grad, new_oil_rate, oil_is_limited,
            new_gas_rate, gas_is_limited, new_alq, alq_is_limited);
    }
    else {
        return std::nullopt;
    }
}

std::unique_ptr<GasLiftWellState>
GasLiftSingleWellGeneric::
runOptimize(const int iteration_idx)
{
    std::unique_ptr<GasLiftWellState> state;
    if (this->optimize_) {
        if (this->debug_limit_increase_decrease_) {
            state = runOptimize1_();
        }
        else {
            state = runOptimize2_();
        }
        if (state) {
            // NOTE: that state->increase() returns a std::optional<bool>, if
            //  this is std::nullopt it means that we was not able to change ALQ
            //  (either increase or decrease)
            if (state->increase()) {  // ALQ changed..
                double alq = state->alq();
                if (this->debug_)
                    logSuccess_(alq, iteration_idx);
                this->well_state_.setALQ(this->well_name_, alq);
            }
        }
    }
    return state;
}

/****************************************
 * Protected methods in alphabetical order
 ****************************************/

std::pair<std::optional<double>, bool>
GasLiftSingleWellGeneric::
addOrSubtractAlqIncrement_(double alq, bool increase) const
{
    bool limited = false;
    double orig_alq = alq;
    if (increase) {
        alq += this->increment_;
        // NOTE: if max_alq_ was defaulted in WLIFTOPT, item 3, it has
        //   already been set to the largest value in the VFP table in
        //   the contructor of GasLiftSingleWell
        if (alq > this->max_alq_) {
            alq = this->max_alq_;
            limited = true;
        }
    }
    else { // we are decreasing ALQ
        alq -= this->increment_;
        if (this->min_alq_ > 0) {
            // According to WLIFTOPT item 5: If a positive value is
            // specified (min_alq_), the well is allocated at least that amount
            // of lift gas, unless the well is unable to flow with
            // that rate of lift gas injection, or unless the well can
            // already meet one of its own rate limits before
            // receiving its minimum lift gas rate.
            if (alq < this->min_alq_) {
                alq = this->min_alq_;
                limited = true;
            }
        }
        else {
            if (alq < 0) {
                alq = 0.0;
                limited = true;
            }
        }
    }
    std::optional<double> alq_opt {alq};
    // If we were not able to change ALQ (to within rounding error), we
    //   return std::nullopt
    if (limited && checkALQequal_(orig_alq,alq))
        alq_opt = std::nullopt;

    return {alq_opt, limited};
}

double
GasLiftSingleWellGeneric::
calcEcoGradient_(double oil_rate, double new_oil_rate, double gas_rate,
                 double new_gas_rate, bool increase) const
{
    auto dqo = new_oil_rate - oil_rate;
    auto dqg = new_gas_rate - gas_rate;
    // TODO: Should the gas rate term in the denominator be subject to the constraint:
    //
    //    alpha_g_ * dqg >= 0.0
    //
    //   ?
    auto gradient = (this->alpha_w_ * dqo) / (this->increment_ + this->alpha_g_*dqg);
    // TODO: Should we do any error checks on the calculation of the
    //   gradient?

    if (!increase) gradient = -gradient;
    return gradient;
}

bool
GasLiftSingleWellGeneric::
checkALQequal_(double alq1, double alq2) const
{
    return std::fabs(alq1-alq2) < (this->increment_*ALQ_EPSILON);
}

bool
GasLiftSingleWellGeneric::
checkInitialALQmodified_(double alq, double initial_alq) const
{
    if (checkALQequal_(alq,initial_alq)) {
        return false;
    }
    else {
        const std::string msg = fmt::format("initial ALQ changed from {} "
            "to {} before iteration starts..", initial_alq, alq);
        displayDebugMessage_(msg);
        return true;
    }
}

bool
GasLiftSingleWellGeneric::
checkWellRatesViolated_(
    std::vector<double>& potentials,
    const std::function<bool(double, double, const std::string &)>& callback,
    bool increase)
{
    if (!increase) {
        auto oil_rate = -potentials[this->oil_pos_];
        if (oil_rate < 0) {
            // The well is not flowing, and it will(?) not help to reduce lift
            // gas further. Note that this assumes that the oil rates drops with
            // decreasing lift gas.
            displayDebugMessage_("Negative oil rate detected while descreasing "
                "lift gas. Stopping iteration.");
            return true;
        }
    }
    // TODO: the below checks could probably be skipped if we are decreasing
    //   lift gas (provided we can assume that rates declines monotonically with
    //   decreasing lift gas).
    if (this->controls_.hasControl(Well::ProducerCMode::ORAT)) {
        auto oil_rate = -potentials[this->oil_pos_];
        if (callback(oil_rate, this->controls_.oil_rate, "oil"))
            return true;
    }
    if (this->controls_.hasControl(Well::ProducerCMode::WRAT)) {
        auto water_rate = -potentials[this->water_pos_];
        if (callback(water_rate, this->controls_.water_rate, "water"))
            return true;
    }
    if (this->controls_.hasControl(Well::ProducerCMode::GRAT)) {
        auto gas_rate = -potentials[this->gas_pos_];
        if (callback(gas_rate, this->controls_.gas_rate, "gas"))
            return true;
    }
    if (this->controls_.hasControl(Well::ProducerCMode::LRAT)) {
        auto oil_rate = -potentials[this->oil_pos_];
        auto water_rate = -potentials[this->water_pos_];
        auto liq_rate = oil_rate + water_rate;
        if (callback(liq_rate, this->controls_.liquid_rate, "liquid"))
            return true;
    }
    // TODO: Also check RESV, see checkIndividualContraints() in
    //   WellInterface_impl.hpp
    // TODO: Check group contraints?

    return false;
}

bool
GasLiftSingleWellGeneric::
computeInitialWellRates_(std::vector<double>& potentials)
{

    if (auto bhp = computeBhpAtThpLimit_(this->orig_alq_); bhp) {
        {
            const std::string msg = fmt::format(
                "computed initial bhp {} given thp limit and given alq {}",
                *bhp, this->orig_alq_);
            displayDebugMessage_(msg);
        }
        computeWellRates_(*bhp, potentials);
        {
            const std::string msg = fmt::format(
                "computed initial well potentials given bhp, "
                "oil: {}, gas: {}, water: {}",
                -potentials[this->oil_pos_],
                -potentials[this->gas_pos_],
                -potentials[this->water_pos_]);
            displayDebugMessage_(msg);
        }
        return true;
    }
    else {
        displayDebugMessage_("Aborting optimization.");
        return false;
    }
}

void
GasLiftSingleWellGeneric::
debugCheckNegativeGradient_(double grad, double alq, double new_alq,
                            double oil_rate, double new_oil_rate,
                            double gas_rate, double new_gas_rate, bool increase) const
{
    {
        const std::string msg = fmt::format("calculating gradient: "
            "new_oil_rate = {}, oil_rate = {}, grad = {}", new_oil_rate, oil_rate, grad);
        displayDebugMessage_(msg);
    }
    if (grad < 0 ) {
        const std::string msg = fmt::format("negative {} gradient detected ({}) : "
            "alq: {}, new_alq: {}, "
            "oil_rate: {}, new_oil_rate: {}, gas_rate: {}, new_gas_rate: {}",
            (increase ? "incremental" : "decremental"),
            grad, alq, new_alq, oil_rate, new_oil_rate, gas_rate, new_gas_rate);
        displayDebugMessage_(msg);
    }
}

void
GasLiftSingleWellGeneric::
debugShowAlqIncreaseDecreaseCounts_()
{
    auto inc_count = this->well_state_.gliftGetAlqIncreaseCount(this->well_name_);
    auto dec_count = this->well_state_.gliftGetAlqDecreaseCount(this->well_name_);
    const std::string msg =
        fmt::format("ALQ increase/decrease count : {}/{}", inc_count, dec_count);
    displayDebugMessage_(msg);
}

void
GasLiftSingleWellGeneric::
debugShowBhpAlqTable_()
{
    double alq = 0.0;
    const std::string fmt_fmt1 {"{:^12s} {:^12s} {:^12s} {:^12s}"};
    const std::string fmt_fmt2 {"{:>12.5g} {:>12.5g} {:>12.5g} {:>12.5g}"};
    const std::string header = fmt::format(fmt_fmt1, "ALQ", "BHP", "oil", "gas");
    displayDebugMessage_(header);
    while (alq <= (this->max_alq_+this->increment_)) {
        auto bhp_at_thp_limit = computeBhpAtThpLimit_(alq);
        if (!bhp_at_thp_limit) {
            const std::string msg = fmt::format("Failed to get converged potentials "
                "for ALQ = {}. Skipping.", alq );
            displayDebugMessage_(msg);
        }
        else {
            std::vector<double> potentials(this->num_phases_, 0.0);
            computeWellRates_(*bhp_at_thp_limit, potentials, /*debug_out=*/false);
            auto oil_rate = -potentials[this->oil_pos_];
            auto gas_rate = -potentials[this->gas_pos_];
            const std::string msg = fmt::format(
                fmt_fmt2, alq, *bhp_at_thp_limit, oil_rate, gas_rate);
            displayDebugMessage_(msg);
        }
        alq += this->increment_;
    }
}

void
GasLiftSingleWellGeneric::
debugShowStartIteration_(double alq, bool increase, double oil_rate)
{
    const std::string msg =
        fmt::format("starting {} iteration, ALQ = {}, oilrate = {}",
            (increase ? "increase" : "decrease"),
            alq, oil_rate);
    displayDebugMessage_(msg);
}

void
GasLiftSingleWellGeneric::
debugShowTargets_()
{
    if (this->controls_.hasControl(Well::ProducerCMode::ORAT)) {
        auto target = this->controls_.oil_rate;
        const std::string msg = fmt::format("has ORAT control with target {}", target);
        displayDebugMessage_(msg);
    }
    if (this->controls_.hasControl(Well::ProducerCMode::GRAT)) {
        auto target = this->controls_.gas_rate;
        const std::string msg = fmt::format("has GRAT control with target {}", target);
        displayDebugMessage_(msg);
    }
    if (this->controls_.hasControl(Well::ProducerCMode::LRAT)) {
        auto target = this->controls_.liquid_rate;
        const std::string msg = fmt::format("has LRAT control with target {}", target);
        displayDebugMessage_(msg);
    }
}

void
GasLiftSingleWellGeneric::
displayDebugMessage_(const std::string& msg) const
{

    if (this->debug_) {
        const std::string message = fmt::format(
            "  GLIFT (DEBUG) : Well {} : {}", this->well_name_, msg);
        this->deferred_logger_.info(message);
    }
}

void
GasLiftSingleWellGeneric::
displayWarning_(const std::string& msg)
{
    const std::string message = fmt::format(
        "GAS LIFT OPTIMIZATION, WELL {} : {}", this->well_name_, msg);
    this->deferred_logger_.warning("WARNING", message);
}

std::pair<double, bool>
GasLiftSingleWellGeneric::
getBhpWithLimit_(double bhp) const
{
    bool limited = false;
    if (this->controls_.hasControl(Well::ProducerCMode::BHP)) {
        auto limit = this->controls_.bhp_limit;
        if (bhp < limit) {
            bhp = limit;
            limited = true;
        }
    }
    return {bhp, limited};
}

// TODO: what if the gas_rate_target_ has been defaulted
//   (i.e. value == 0, meaning: "No limit") but the
//   oil_rate_target_ has not been defaulted ?
//   If the new_oil_rate exceeds the oil_rate_target_ it is cut back,
//   but the same cut-back will not happen for the new_gas_rate
//   Seems like an inconsistency, since alq should in this
//   case also be adjusted (to the smaller value that would
//   give oil target rate) but then the gas rate would also be smaller?
//   The effect of not reducing the gas rate (if it should be
//   reduced?) is that a too large value is used in the
//   computation of the economic gradient making the gradient
//   smaller than it should be since the term appears in the denominator.
std::pair<double, bool>
GasLiftSingleWellGeneric::
getGasRateWithLimit_(const std::vector<double>& potentials) const
{
    double new_rate = -potentials[this->gas_pos_];
    bool limit = false;
    if (this->controls_.hasControl(Well::ProducerCMode::GRAT)) {
        auto target = this->controls_.gas_rate;
        if (new_rate > target) {
            new_rate = target;
            limit = true;
        }
    }
    return { new_rate, limit};
}

std::pair<double, bool>
GasLiftSingleWellGeneric::
getWaterRateWithLimit_(const std::vector<double>& potentials) const
{
    double new_rate = -potentials[this->water_pos_];
    bool limit = false;
    if (this->controls_.hasControl(Well::ProducerCMode::WRAT)) {
        auto target = this->controls_.water_rate;
        if (new_rate > target) {
            new_rate = target;
            limit = true;
        }
    }
    return { new_rate, limit};
}

// NOTE: If the computed oil rate is larger than the target
//   rate of the well, we reduce it to the target rate. This
//   will make the economic gradient smaller than it would be
//   if we did not reduce the rate, and it is less
//   likely that the current gas lift increment will be
//   accepted.
// TODO: If it still is accepted, we should ideally reduce the alq
//  also since we also reduced the rate. This might involve
//   some sort of iteration though..
std::pair<double, bool>
GasLiftSingleWellGeneric::
getOilRateWithLimit_(const std::vector<double>& potentials) const
{
    double oil_rate = -potentials[this->oil_pos_];
    double new_rate = oil_rate;
    bool limited = false;
    if (this->controls_.hasControl(Well::ProducerCMode::ORAT)) {
        auto target = this->controls_.oil_rate;
        if (oil_rate > target) {
            const std::string msg = fmt::format("limiting oil rate to target: "
                "computed rate: {}, target: {}", new_rate, target);
            displayDebugMessage_(msg);
            new_rate = target;
            limited = true;
        }
    }
    if (this->controls_.hasControl(Well::ProducerCMode::LRAT)) {
        auto target = this->controls_.liquid_rate;
        double water_rate = -potentials[this->water_pos_];
        double liq_rate = oil_rate + water_rate;
        if (liq_rate > target) {
            double oil_fraction = oil_rate / liq_rate;
            new_rate = std::min(new_rate, oil_fraction * target);
            limited = true;
            const std::string msg = fmt::format(
                "limiting oil rate due to LRAT target: "
                "computed rate: {}, target: {}", oil_rate, new_rate);
            displayDebugMessage_(msg);
       }
    }
    return { new_rate, limited};
}


std::pair<double, bool>
GasLiftSingleWellGeneric::
getOilRateWithGroupLimit_(const double new_oil_rate, const double oil_rate) const
{
    const double delta_oil = new_oil_rate - oil_rate;
    const auto &pairs =
        this->group_info_.getWellGroups(this->well_name_);
    for (const auto &[group_name, efficiency] : pairs) {
        auto gr_oil_target_opt = this->group_info_.oilTarget(group_name);
        if (gr_oil_target_opt) {
            double gr_oil_rate =
                this->group_info_.oilRate(group_name);
            double new_gr_oil_rate = gr_oil_rate + efficiency * delta_oil;
            if (new_gr_oil_rate > *gr_oil_target_opt) {
                const std::string msg = fmt::format("limiting oil rate to group target: "
                    "computed group rate: {}, target: {}", new_gr_oil_rate, *gr_oil_target_opt);
                displayDebugMessage_(msg);
                double new_rate = oil_rate + (*gr_oil_target_opt - gr_oil_rate) / efficiency;
                return { std::min(new_rate, new_oil_rate), /*limit=*/true};
            }
        }
    }
    return { new_oil_rate, /*limit=*/false};
}

std::pair<double, bool>
GasLiftSingleWellGeneric::
getGasRateWithGroupLimit_(const double new_gas_rate, const double gas_rate) const
{
    const double delta_gas = new_gas_rate - gas_rate;
    const auto &pairs =
        this->group_info_.getWellGroups(this->well_name_);
    for (const auto &[group_name, efficiency] : pairs) {
        auto gr_gas_target_opt = this->group_info_.gasTarget(group_name);
        if (gr_gas_target_opt) {
            double gr_gas_rate =
                this->group_info_.gasRate(group_name);
            double new_gr_gas_rate = gr_gas_rate + efficiency * delta_gas;
            if (new_gr_gas_rate > *gr_gas_target_opt) {
                const std::string msg = fmt::format("limiting gas rate to group target: "
                    "computed group rate: {}, target: {}", new_gr_gas_rate, *gr_gas_target_opt);
                displayDebugMessage_(msg);
                double new_rate = gas_rate + (*gr_gas_target_opt - gr_gas_rate) / efficiency;
                return { std::min(new_rate, new_gas_rate), /*limit=*/true};
            }
        }
    }
    return { new_gas_rate, /*limit=*/false};
}

std::pair<double, bool>
GasLiftSingleWellGeneric::
getWaterRateWithGroupLimit_(const double new_water_rate, const double water_rate) const
{
    const double delta_water = new_water_rate - water_rate;
    const auto &pairs =
        this->group_info_.getWellGroups(this->well_name_);
    for (const auto &[group_name, efficiency] : pairs) {
        auto gr_water_target_opt = this->group_info_.waterTarget(group_name);
        if (gr_water_target_opt) {
            double gr_water_rate =
                this->group_info_.waterRate(group_name);
            double new_gr_water_rate = gr_water_rate + efficiency * delta_water;
            if (new_gr_water_rate > *gr_water_target_opt) {
                const std::string msg = fmt::format("limiting water rate to group target: "
                    "computed group rate: {}, target: {}", new_gr_water_rate, *gr_water_target_opt);
                displayDebugMessage_(msg);
                double new_rate = water_rate + (*gr_water_target_opt - gr_water_rate) / efficiency;
                return { std::min(new_rate, new_water_rate), /*limit=*/true};
            }
        }
    }
    return { new_water_rate, /*limit=*/false};
}

std::tuple<double, double, bool, bool>
GasLiftSingleWellGeneric::
getLiquidRateWithGroupLimit_(const double new_oil_rate, const double oil_rate,
                             const double new_water_rate, const double water_rate) const
{
    const double delta_water = new_water_rate - water_rate;
    const double delta_oil = new_oil_rate - oil_rate;
    const auto &pairs =
        this->group_info_.getWellGroups(this->well_name_);
    for (const auto &[group_name, efficiency] : pairs) {
        auto gr_liquid_target_opt = this->group_info_.liquidTarget(group_name);
        if (gr_liquid_target_opt) {
            double gr_water_rate =
                this->group_info_.waterRate(group_name);
            double gr_oil_rate =
                this->group_info_.oilRate(group_name);
            double new_gr_water_rate = gr_water_rate + efficiency * delta_water;
            double new_gr_oil_rate = gr_oil_rate + efficiency * delta_oil;
            double new_gr_liquid_rate = new_gr_water_rate + new_gr_oil_rate;
            if (new_gr_liquid_rate > *gr_liquid_target_opt) {
                const std::string msg = fmt::format("limiting liquid rate to group target: "
                    "computed group rate: {}, target: {}", new_gr_liquid_rate, *gr_liquid_target_opt);
                displayDebugMessage_(msg);
                double oil_fraction = new_gr_oil_rate / new_gr_liquid_rate;
                double water_rate_limited = water_rate + (1.0 - oil_fraction) * (new_gr_liquid_rate - *gr_liquid_target_opt) / efficiency;
                double oil_rate_limited = oil_rate + oil_fraction * (new_gr_liquid_rate - *gr_liquid_target_opt) / efficiency;
                return { std::min(oil_rate_limited, new_oil_rate), std::min(water_rate_limited, new_water_rate), /*limit=*/true, /*limit=*/true};
            }
        }
    }
    return { new_oil_rate, new_water_rate, /*limit=*/false, /*limit=*/false};
}

std::tuple<double,double,double, bool, bool,bool>
GasLiftSingleWellGeneric::
getInitialRatesWithLimit_(const std::vector<double>& potentials)
{
    auto [oil_rate, oil_is_limited] = getOilRateWithLimit_(potentials);
    auto [gas_rate, gas_is_limited] = getGasRateWithLimit_(potentials);
    auto [water_rate, water_is_limited] = getWaterRateWithLimit_(potentials);
    if (oil_is_limited) {
        const std::string msg = fmt::format(
            "initial oil rate was limited to: {}", oil_rate);
        displayDebugMessage_(msg);
    }
    if (gas_is_limited) {
        const std::string msg = fmt::format(
            "initial gas rate was limited to: {}", gas_rate);
        displayDebugMessage_(msg);
    }
    if (water_is_limited) {
        const std::string msg = fmt::format(
            "initial water rate was limited to: {}", water_rate);
        displayDebugMessage_(msg);
    }
    return std::make_tuple(oil_rate, gas_rate, water_rate, oil_is_limited, gas_is_limited, water_is_limited);
}

std::tuple<double,double,bool,bool,double>
GasLiftSingleWellGeneric::
increaseALQtoPositiveOilRate_(double alq,
                              double oil_rate,
                              double gas_rate,
                              bool oil_is_limited,
                              bool gas_is_limited,
                              std::vector<double>& potentials)
{
    bool stop_iteration = false;
    double temp_alq = alq;
    while(!stop_iteration) {
        temp_alq += this->increment_;
        if (temp_alq > this->max_alq_) break;
        auto bhp_opt = computeBhpAtThpLimit_(temp_alq);
        if (!bhp_opt) break;
        alq = temp_alq;
        auto bhp_this = getBhpWithLimit_(*bhp_opt);
        computeWellRates_(bhp_this.first, potentials);
        oil_rate = -potentials[this->oil_pos_];
        if (oil_rate > 0) break;
    }
    std::tie(oil_rate, oil_is_limited) = getOilRateWithLimit_(potentials);
    std::tie(gas_rate, gas_is_limited) = getGasRateWithLimit_(potentials);
    return std::make_tuple(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq);
}

std::tuple<double,double,bool,bool,double>
GasLiftSingleWellGeneric::
increaseALQtoMinALQ_(double alq,
            double oil_rate,
            double gas_rate,
            bool oil_is_limited,
            bool gas_is_limited,
            std::vector<double>& potentials)
{
    auto min_alq = this->min_alq_;
    assert(min_alq >= 0);
    assert(alq < min_alq);
    assert(min_alq < this->max_alq_);
    bool stop_iteration = false;
    double temp_alq = alq;
    while(!stop_iteration) {
        temp_alq += this->increment_;
        if (temp_alq >= min_alq) break;
        auto bhp_opt = computeBhpAtThpLimit_(temp_alq);
        if (!bhp_opt) break;
        alq = temp_alq;
        auto bhp_this = getBhpWithLimit_(*bhp_opt);
        computeWellRates_(bhp_this.first, potentials);
        std::tie(oil_rate, oil_is_limited) = getOilRateWithLimit_(potentials);
        std::tie(gas_rate, gas_is_limited) = getGasRateWithLimit_(potentials);
        if (oil_is_limited || gas_is_limited) break;
    }
    return std::make_tuple(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq);
}

void
GasLiftSingleWellGeneric::
logSuccess_(double alq, const int iteration_idx)
{
    const std::string message = fmt::format(
         "GLIFT, IT={}, WELL {} : {} ALQ from {} to {}",
         iteration_idx,
         this->well_name_,
         ((alq > this->orig_alq_) ? "increased" : "decreased"),
         this->orig_alq_, alq);
    this->deferred_logger_.info(message);
}

std::tuple<double,double,double,double,bool,bool,bool>
GasLiftSingleWellGeneric::
maybeAdjustALQbeforeOptimizeLoop_(
    bool increase, double alq, double oil_rate, double gas_rate, double water_rate,
    bool oil_is_limited, bool gas_is_limited,bool water_is_limited,
    std::vector<double> &potentials)
{
    double orig_alq = alq;
    if (this->debug_) {
        const std::string msg = fmt::format("initial ALQ: {}", alq);
        displayDebugMessage_(msg);
    }
    if (!increase && oil_is_limited) {
        // NOTE: Try to decrease ALQ down to a value where the oil target is
        //   not exceeded.
        // NOTE: This may reduce ALQ below the minimum value set in WLIFTOPT
        //   item 5. However, this is OK since the rate target is met and there
        //   is no point in using a higher ALQ value then.
        std::tie(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq) =
            reduceALQtoOilTarget_(
                alq, oil_rate, gas_rate, oil_is_limited, gas_is_limited, potentials);
    }
    else {
        if (increase && oil_rate < 0) {
            // NOTE: Try to increase ALQ up to a value where oil_rate is positive
            std::tie(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq) =
               increaseALQtoPositiveOilRate_(alq, oil_rate,
                   gas_rate, oil_is_limited, gas_is_limited, potentials);
        }
        if (increase && (this->min_alq_> 0) && (alq < this->min_alq_)) {
            // NOTE: Try to increase ALQ up to the minimum limit without checking
            //   the economic gradient..
            std::tie(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq) =
                increaseALQtoMinALQ_(alq, oil_rate, gas_rate,
                    oil_is_limited, gas_is_limited, potentials);
        }
    }
    if (this->debug_ && (orig_alq != alq)) {
        const std::string msg = fmt::format("adjusted ALQ to: {}", alq);
        displayDebugMessage_(msg);
    }
    return std::make_tuple(oil_rate, gas_rate, water_rate, alq, oil_is_limited, gas_is_limited, water_is_limited);
}

std::tuple<double,double,bool,bool,double>
GasLiftSingleWellGeneric::
reduceALQtoOilTarget_(double alq,
                      double oil_rate,
                      double gas_rate,
                      bool oil_is_limited,
                      bool gas_is_limited,
                      std::vector<double>& potentials)
{
    displayDebugMessage_("Reducing ALQ to meet oil target before iteration starts..");
    double orig_oil_rate = oil_rate;
    double orig_alq = alq;
    // NOTE: This method should only be called if oil_is_limited, and hence
    //   we know that it has oil rate control
    assert(this->controls_.hasControl(Well::ProducerCMode::ORAT) || this->controls_.hasControl(Well::ProducerCMode::LRAT));
    auto target = this->controls_.oil_rate;
    bool stop_iteration = false;
    double temp_alq = alq;
    while(!stop_iteration) {
        temp_alq -= this->increment_;
        if (temp_alq <= 0) break;
        auto bhp_opt = computeBhpAtThpLimit_(temp_alq);
        if (!bhp_opt) break;
        auto bhp_this = getBhpWithLimit_(*bhp_opt);
        computeWellRates_(bhp_this.first, potentials);
        oil_rate = -potentials[this->oil_pos_];
        if (oil_rate < target) {
            break;
        }
        alq = temp_alq;
    }
    std::tie(oil_rate, oil_is_limited) = getOilRateWithLimit_(potentials);
    std::tie(gas_rate, gas_is_limited) = getGasRateWithLimit_(potentials);
    if (this->debug_) {
        assert( alq <= orig_alq );
        if (alq < orig_alq) {
            // NOTE: ALQ may drop below zero before we are able to meet the target
            const std::string msg = fmt::format(
                "Reduced (oil_rate, alq) from ({}, {}) to ({}, {}) to meet target "
                "at {}. ", orig_oil_rate, orig_alq, oil_rate, alq, target);
            displayDebugMessage_(msg);
        }
        else if (alq == orig_alq) {
            // We might not be able to reduce ALQ, for example if ALQ starts out at zero.
            const std::string msg = fmt::format("Not able to reduce ALQ {} further. "
                "Oil rate is {} and oil target is {}", orig_alq, oil_rate, target);
            displayDebugMessage_(msg);
        }
    }
    return std::make_tuple(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq);
}

// INPUT:
//  - increase (boolean) :
//   - true  : try increase the lift gas supply,
//   - false : try decrease lift gas supply.
//
// OUTPUT:
//
//  - return value: a new GasLiftWellState or nullptr
//
std::unique_ptr<GasLiftWellState>
GasLiftSingleWellGeneric::
runOptimizeLoop_(bool increase)
{
    std::vector<double> potentials(this->num_phases_, 0.0);
    std::unique_ptr<GasLiftWellState> ret_value; // nullptr initially
    if (!computeInitialWellRates_(potentials)) return ret_value;
    bool alq_is_limited = false;
    bool oil_is_limited = false;
    bool gas_is_limited = false;
    bool water_is_limited = false;
    double oil_rate, gas_rate, water_rate;
    std::tie(oil_rate, gas_rate, water_rate, oil_is_limited, gas_is_limited, water_is_limited) =
        getInitialRatesWithLimit_(potentials);
    //if (this->debug_) debugShowBhpAlqTable_();
    if (this->debug_) debugShowAlqIncreaseDecreaseCounts_();
    if (this->debug_) debugShowTargets_();
    bool success = false;  // did we succeed to increase alq?
    auto cur_alq = this->orig_alq_;
    double new_oil_rate, new_gas_rate, new_water_rate, new_alq;
    bool new_oil_is_limited, new_gas_is_limited, new_water_is_limited;
    std::tie(new_oil_rate, new_gas_rate, new_water_rate, new_alq,
             new_oil_is_limited, new_gas_is_limited, new_water_is_limited)
        = maybeAdjustALQbeforeOptimizeLoop_(
            increase, cur_alq, oil_rate, gas_rate, water_rate,
            oil_is_limited, gas_is_limited, water_is_limited, potentials);
    double delta_oil = 0.0;
    double delta_gas = 0.0;
    double delta_alq = 0.0;
    double delta_water = 0.0;
    OptimizeState state {*this, increase};

    // potentially reduce alq if group control is violated
    std::tie(new_oil_rate, new_gas_rate, new_water_rate, new_alq) =
        state.reduceALQtoGroupTarget(new_alq, new_oil_rate, new_gas_rate, new_water_rate, potentials);

    if (checkInitialALQmodified_(new_alq, cur_alq)) {
        delta_oil = new_oil_rate - oil_rate;
        delta_gas = new_gas_rate - gas_rate;
        delta_water = new_water_rate - water_rate;
        delta_alq = new_alq - cur_alq;
        if (!(state.checkGroupTargetsViolated(delta_oil, delta_gas, delta_water)) &&
            !(state.checkGroupALQrateExceeded(delta_alq)))
        {
            oil_rate = new_oil_rate;
            gas_rate = new_gas_rate;
            water_rate = new_water_rate;
            oil_is_limited = new_oil_is_limited;
            gas_is_limited = new_gas_is_limited;
            water_is_limited = new_water_is_limited;
            cur_alq = new_alq;
            success = true;
        }
        else {
            state.stop_iteration = true;
        }
    }
    // we only iterate if well is under thp control
    if (!state.checkThpControl()) {
        state.stop_iteration = true;
    }

    auto temp_alq = cur_alq;
    if (this->debug_) debugShowStartIteration_(temp_alq, increase, oil_rate);
    while (!state.stop_iteration && (++state.it <= this->max_iterations_)) {
        if (!increase && state.checkNegativeOilRate(oil_rate)) break;
        if (state.checkWellRatesViolated(potentials)) break;
        if (state.checkGroupTargetsViolated(delta_oil, delta_gas, delta_water)) break;
        if (state.checkAlqOutsideLimits(temp_alq, oil_rate)) break;
        std::optional<double> alq_opt;
        std::tie(alq_opt, alq_is_limited)
            = state.addOrSubtractAlqIncrement(temp_alq);
        if (!alq_opt) break;
        delta_alq = *alq_opt - temp_alq;
        if (state.checkGroupALQrateExceeded(delta_alq)) break;
        temp_alq = *alq_opt;
        if (this->debug_) state.debugShowIterationInfo(temp_alq);
        if (!state.computeBhpAtThpLimit(temp_alq)) break;
        // NOTE: if BHP is below limit, we set state.stop_iteration = true
        auto bhp = state.getBhpWithLimit();
        computeWellRates_(bhp, potentials);
        std::tie(new_oil_rate, new_oil_is_limited) = getOilRateWithLimit_(potentials);
        std::tie(new_oil_rate, new_oil_is_limited) = getOilRateWithGroupLimit_(new_oil_rate, oil_rate);

/*        if (this->debug_abort_if_decrease_and_oil_is_limited_) {
            if (oil_is_limited && !increase) {
                // if oil is limited we do not want to decrease
                displayDebugMessage_(
                    "decreasing ALQ and oil is limited -> aborting iteration");
                break;
            }
            }
*/
        std::tie(new_gas_rate, new_gas_is_limited) = getGasRateWithLimit_(potentials);
        std::tie(new_gas_rate, new_gas_is_limited) = getGasRateWithGroupLimit_(new_gas_rate, gas_rate);

        std::tie(new_water_rate, new_water_is_limited) = getWaterRateWithLimit_(potentials);
        std::tie(new_water_rate, new_water_is_limited) = getWaterRateWithGroupLimit_(new_water_rate, water_rate);

        std::tie(new_oil_rate, new_water_rate, new_oil_is_limited, new_water_is_limited)
            = getLiquidRateWithGroupLimit_(new_oil_rate, oil_rate, new_water_rate, water_rate);


/*        if (this->debug_abort_if_increase_and_gas_is_limited_) {
            if (gas_is_limited && increase) {
                // if gas is limited we do not want to increase
                displayDebugMessage_(
                    "increasing ALQ and gas is limited -> aborting iteration");
                break;
            }
        }
*/
        auto gradient = state.calcEcoGradient(
            oil_rate, new_oil_rate, gas_rate, new_gas_rate);
        if (this->debug_)
            debugCheckNegativeGradient_(
                gradient, cur_alq, temp_alq, oil_rate, new_oil_rate,
                gas_rate, new_gas_rate, increase);
        if (state.checkEcoGradient(gradient)) break;
        cur_alq = temp_alq;
        success = true;
        delta_oil = new_oil_rate - oil_rate;
        delta_gas = new_gas_rate - gas_rate;
        delta_water = new_water_rate - water_rate;
        oil_rate = new_oil_rate;
        gas_rate = new_gas_rate;
        water_rate = new_water_rate;
        oil_is_limited = new_oil_is_limited;
        gas_is_limited = new_gas_is_limited;
        water_is_limited = new_water_is_limited;
        state.updateGroupRates(delta_oil, delta_gas, delta_water, delta_alq);
    }
    if (state.it > this->max_iterations_) {
        warnMaxIterationsExceeded_();
    }
    std::optional<bool> increase_opt;
    if (success) {
        this->well_state_.gliftUpdateAlqIncreaseCount(this->well_name_, increase);
        increase_opt = increase;
    }
    else {
        increase_opt = std::nullopt;
    }
    ret_value = std::make_unique<GasLiftWellState>(oil_rate, oil_is_limited,
        gas_rate, gas_is_limited, cur_alq, alq_is_limited, increase_opt);
    return ret_value;
}

std::unique_ptr<GasLiftWellState>
GasLiftSingleWellGeneric::
runOptimize1_()
{
    std::unique_ptr<GasLiftWellState> state;
    int inc_count = this->well_state_.gliftGetAlqIncreaseCount(this->well_name_);
    int dec_count = this->well_state_.gliftGetAlqDecreaseCount(this->well_name_);
    if (dec_count == 0 && inc_count == 0) {
        state = tryIncreaseLiftGas_();
        if (!state || !(state->alqChanged())) {
            state = tryDecreaseLiftGas_();
        }
    }
    else if (dec_count == 0) {
        assert(inc_count > 0);
        state = tryIncreaseLiftGas_();
    }
    else if (inc_count == 0) {
        assert(dec_count > 0);
        state = tryDecreaseLiftGas_();
    }
    return state;
}

std::unique_ptr<GasLiftWellState>
GasLiftSingleWellGeneric::
runOptimize2_()
{
    std::unique_ptr<GasLiftWellState> state;
    state = tryIncreaseLiftGas_();
    if (!state || !(state->alqChanged())) {
        state = tryDecreaseLiftGas_();
    }
    return state;
}

std::unique_ptr<GasLiftWellState>
GasLiftSingleWellGeneric::
tryDecreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/ false);
}

std::unique_ptr<GasLiftWellState>
GasLiftSingleWellGeneric::
tryIncreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/ true);
}

void
GasLiftSingleWellGeneric::
setAlqMinRate_(const GasLiftOpt::Well& well)
{
    // NOTE:  According to WLIFTOPT item 5 :
    //   if min_rate() is negative, it means: allocate at least enough lift gas
    //   to enable the well to flow
    // NOTE: "to enable the well to flow" : How to interpret this?
    //   We choose to interpret it to mean a positive oil rate as returned from
    //
    //    computeWellRates_(bhp, cur_potentials);
    //
    //   So even if the well is producing gas, if the oil rate is zero
    //   we say that the "well is not flowing".
    //
    //   Note that if WECON item 2 is set, the well can be shut off
    //   before the flow rate reaches zero. Also,
    //   if bhp drops below the bhp lower limit, the well might switch to bhp
    //   control before the oil rate becomes zero.

    this->min_alq_ = well.min_rate();
    if (this->min_alq_ > 0) {
        if (this->min_alq_ >= this->max_alq_) {
            // NOTE: We reset the value to a negative value.
            //   negative value means: Allocate at least enough lift gas
            //   to allow the well to flow.
            // TODO: Consider other options for resetting the value..
            this->min_alq_ = -1;
            displayWarning_("Minimum ALQ value is larger than maximum ALQ value!"
                " Resetting value.");
        }
    }
}

// Called when we should use a fixed ALQ value
void
GasLiftSingleWellGeneric::
updateWellStateAlqFixedValue_(const GasLiftOpt::Well& well)
{
    auto& max_alq_optional = well.max_rate();
    if (max_alq_optional) {
        // According to WLIFTOPT, item 3:
        // If item 2 is NO, then item 3 is regarded as the fixed
        // lift gas injection rate for the well.
        auto new_alq = *max_alq_optional;
        this->well_state_.setALQ(this->well_name_, new_alq);
    }
    // else {
    //    // If item 3 is defaulted, the lift gas rate remains
    //    // unchanged at its current value.
    //}

}

// Determine if we should use a fixed ALQ value.
//
// From the manual for WLIFTOPT, item 2:
//   Is the well's lift gas injection rate to be calculated by the
//   optimization facility?
// - YES : The well's lift gas injection rate is calculated by the
//   optimization facility.
// - NO  : The well's lift gas injection rate remains fixed at a
//   value that can be set either in Item 3 of this keyword, or in
//   Item 12 of keyword WCONPROD, or with keyword WELTARG.
bool
GasLiftSingleWellGeneric::
useFixedAlq_(const GasLiftOpt::Well& well)
{
    auto wliftopt_item2 = well.use_glo();
    if (wliftopt_item2) {
        return false;
    }
    else {
        //  auto& max_alq_optional = well.max_rate();
        //  if (max_alq_optional) {
               // According to WLIFTOPT, item 3:
               // If item 2 is NO, then item 3 is regarded as the fixed
               // lift gas injection rate for the well.
        //  }
        //  else {
              // If item 3 is defaulted, the lift gas rate remains
              // unchanged at its current value.
        //  }
        return true;
    }
}
void
GasLiftSingleWellGeneric::
warnMaxIterationsExceeded_()
{
    const std::string msg = fmt::format(
        "Max iterations ({}) exceeded", this->max_iterations_);
    displayWarning_(msg);
}

/****************************************
 * Methods declared in OptimizeState
 ****************************************/

std::pair<std::optional<double>, bool>
GasLiftSingleWellGeneric::OptimizeState::
addOrSubtractAlqIncrement(double alq)
{
    auto [alq_opt, limited]
        = this->parent.addOrSubtractAlqIncrement_(alq, this->increase);
    if (!alq_opt) {
        const std::string msg = fmt::format(
            "iteration {}, alq = {} : not able to {} ALQ increment",
            this->it, alq, (this->increase ? "add" : "subtract"));
    }
    return {alq_opt, limited};
}

double
GasLiftSingleWellGeneric::OptimizeState::
calcEcoGradient(double oil_rate, double new_oil_rate,
                double gas_rate, double new_gas_rate)
{
    return this->parent.calcEcoGradient_(oil_rate, new_oil_rate,
                                         gas_rate, new_gas_rate, this->increase);
}

// NOTE:  According to WLIFTOPT item 5 :
//   if min_rate() is negative, it means: allocate at least enough lift gas
//   to enable the well to flow
//  We will interpret this as (see discussion above GasLiftSingleWell()
//   in this file): Allocate at least the amount of lift gas needed to
//   get a positive oil production rate.
bool
GasLiftSingleWellGeneric::OptimizeState::
checkAlqOutsideLimits(double alq, [[maybe_unused]] double oil_rate)
{
    std::ostringstream ss;
    bool result = false;

    if (this->increase) {
        if (alq >= this->parent.max_alq_) {
            ss << "ALQ >= " << this->parent.max_alq_ << " (max limit), "
               << "stopping iteration";
            result = true;
        }
        else { // checking the minimum limit...
            // NOTE: A negative min_alq_ means: allocate at least enough lift gas
            //  to enable the well to flow, see WLIFTOPT item 5.
            if (this->parent.min_alq_ < 0) {
                // - if oil rate is negative (i.e. the well is not flowing), continue to
                //    increase ALQ (according WLIFTOPT item 5) and try make the well
                //    flow.
                // - else if oil rate is already positive, there is no minimum
                //    limit for ALQ in this case
                result = false;
            }
            else {
                // NOTE: checking for a lower limit is not necessary
                //   when increasing alq. If ALQ was smaller than the minimum when
                //   we entered the runOptimizeLoop_() method,
                //   increaseALQtoMinALQ_() will ensure that ALQ >= min_alq
                assert(alq >= this->parent.min_alq_ );
                result = false;
            }
        }
    }
    else { // we are decreasing lift gas
        if ( alq == 0 ) {
            ss << "ALQ is zero, cannot decrease further. Stopping iteration.";
            return true;
        }
        else if ( alq < 0 ) {
            ss << "Negative ALQ: " << alq << ". Stopping iteration.";
            return true;
        }
        // NOTE: A negative min_alq_ means: allocate at least enough lift gas
        //  to enable the well to flow, see WLIFTOPT item 5.
        if (this->parent.min_alq_ < 0) {
            // We know that the well is flowing (oil_rate > 0) since that was
            //  already checked in runOptimizeLoop_() by calling checkNegativeOilRate()
            assert(oil_rate >= 0);
            result = false;
        }
        else {
            if (alq <= this->parent.min_alq_ ) {
                // According to WLIFTOPT item 5:
                //   "If a positive value is specified, the well is
                //   allocated at least that amount of lift gas,
                //   unless the well is unable to flow with that rate
                //   of lift gas injection, or unless the well can
                //   already meet one of its own rate limits before
                //   receiving its minimum lift gas rate."
                //
                // - We already know that the well is flowing, (oil_rate > 0),
                //    since that was already checked in runOptimizeLoop_() by calling
                //    checkNegativeOilRate().
                // - We also know that the rate limit was not exceeded since that was
                //    checked by checkWellRatesViolated()
                assert( oil_rate >= 0);
                ss << "ALQ <= " << this->parent.min_alq_ << " (min limit), "
                       << "stopping iteration";
                result = true;
            }
            else {
                // NOTE: checking for an upper limit should not be necessary
                // when decreasing alq.. so this is just to catch an
                // illegal state at an early point.
                if (alq >= this->parent.max_alq_) {
                    warn_( "unexpected: alq above upper limit when trying to "
                        "decrease lift gas. aborting iteration.");
                    result = true;
                }
                else {
                    result = false;
                }
            }
        }
    }
    if (this->parent.debug_) {
        const std::string msg = ss.str();
        if (!msg.empty())
            this->parent.displayDebugMessage_(msg);
    }
    return result;
}

bool
GasLiftSingleWellGeneric::OptimizeState::
checkGroupALQrateExceeded(double delta_alq)
{
    const auto &pairs =
        this->parent.group_info_.getWellGroups(this->parent.well_name_);
    for (const auto &[group_name, efficiency] : pairs) {
        auto max_alq_opt = this->parent.group_info_.maxAlq(group_name);
        if (max_alq_opt) {
            double alq =
                this->parent.group_info_.alqRate(group_name) + efficiency * delta_alq;
            if (alq > *max_alq_opt) {
                if (this->parent.debug_) {
                    const std::string msg = fmt::format(
                        "Group {} : alq {} exceeds max_alq {}. Stopping iteration",
                        group_name, alq, *max_alq_opt);
                    this->parent.displayDebugMessage_(msg);
                }
                return true;
            }
        }
    }
    return false;
}

bool
GasLiftSingleWellGeneric::OptimizeState::
checkGroupTargetsViolated(double delta_oil, double delta_gas, double delta_water)
{
    const auto &pairs =
        this->parent.group_info_.getWellGroups(this->parent.well_name_);
    for (const auto &[group_name, efficiency] : pairs) {
        auto oil_target_opt = this->parent.group_info_.oilTarget(group_name);
        if (oil_target_opt) {
            double oil_rate =
                this->parent.group_info_.oilRate(group_name) + efficiency * delta_oil;
            if (oil_rate > *oil_target_opt) {
                if (this->parent.debug_) {
                    const std::string msg = fmt::format(
                       "Group {} : oil rate {} exceeds oil target {}. Stopping iteration",
                       group_name, oil_rate, *oil_target_opt);
                    this->parent.displayDebugMessage_(msg);
                }
                return true;
            }
        }
        auto gas_target_opt = this->parent.group_info_.gasTarget(group_name);
        if (gas_target_opt) {
            double gas_rate =
                this->parent.group_info_.gasRate(group_name) + efficiency * delta_gas;
            if (gas_rate > *gas_target_opt) {
                if (this->parent.debug_) {
                    const std::string msg = fmt::format(
                       "Group {} : gas rate {} exceeds gas target {}. Stopping iteration",
                       group_name, gas_rate, *gas_target_opt);
                    this->parent.displayDebugMessage_(msg);
                }
                return true;
            }
        }
        auto liquid_target_opt = this->parent.group_info_.liquidTarget(group_name);
        if (liquid_target_opt) {
            double oil_rate =
                this->parent.group_info_.oilRate(group_name) + efficiency * delta_oil;
            double water_rate =
                this->parent.group_info_.waterRate(group_name) + efficiency * delta_water;
            double liquid_rate = oil_rate + water_rate;
            if (liquid_rate > *liquid_target_opt) {
                if (this->parent.debug_) {
                    const std::string msg = fmt::format(
                       "Group {} : liquid rate {} exceeds liquid target {}. Stopping iteration",
                       group_name, liquid_rate, *liquid_target_opt);
                    this->parent.displayDebugMessage_(msg);
                }
                return true;
            }
        }
        auto water_target_opt = this->parent.group_info_.waterTarget(group_name);
        if (water_target_opt) {
            double water_rate =
                this->parent.group_info_.waterRate(group_name) + efficiency * delta_water;
            if (water_rate > *water_target_opt) {
                if (this->parent.debug_) {
                    const std::string msg = fmt::format(
                       "Group {} : water rate {} exceeds water target {}. Stopping iteration",
                       group_name, water_rate, *water_target_opt);
                    this->parent.displayDebugMessage_(msg);
                }
                return true;
            }
        }
    }
    return false;
}

std::tuple<double,double,double,double>
GasLiftSingleWellGeneric::OptimizeState::
reduceALQtoGroupTarget(double alq,
                       double oil_rate,
                       double gas_rate,
                       double water_rate,
                       std::vector<double>& potentials)
{
    bool stop_this_iteration = true;
    const auto &pairs =
        this->parent.group_info_.getWellGroups(this->parent.well_name_);
    for (const auto &groups : pairs) {
        if (!this->parent.group_state_.has_production_control(groups.first))
            continue;
        const auto& current_control = this->parent.group_state_.production_control(groups.first);
        if(current_control == Group::ProductionCMode::ORAT
                || current_control == Group::ProductionCMode::LRAT
                || current_control == Group::ProductionCMode::WRAT
                || current_control == Group::ProductionCMode::GRAT){
            stop_this_iteration = false;
            this->parent.displayDebugMessage_("Reducing ALQ to meet groups target before iteration starts.");
            break;
        }
    }
    double temp_alq = alq;
    double oil_rate_orig = oil_rate;
    double gas_rate_orig = gas_rate;
    double water_rate_orig = water_rate;
    while(!stop_this_iteration) {
        temp_alq -= this->parent.increment_;
        if (temp_alq <= 0) break;
        auto bhp_opt = this->parent.computeBhpAtThpLimit_(temp_alq);
        if (!bhp_opt) break;
        auto bhp_this = this->parent.getBhpWithLimit_(*bhp_opt);
        this->parent.computeWellRates_(bhp_this.first, potentials);
        oil_rate = -potentials[this->parent.oil_pos_];
        gas_rate = -potentials[this->parent.gas_pos_];
        water_rate = -potentials[this->parent.gas_pos_];
        double delta_oil = oil_rate_orig - oil_rate;
        double delta_gas = gas_rate_orig - gas_rate;
        double delta_water = water_rate_orig - water_rate;

        if (!this->checkGroupTargetsViolated(delta_oil, delta_gas, delta_water)) {
            break;
        }
        alq = temp_alq;
    }
    return std::make_tuple(oil_rate, gas_rate, water_rate, alq);
}
bool
GasLiftSingleWellGeneric::OptimizeState::
checkNegativeOilRate(double oil_rate)
{
    if (oil_rate < 0) {
        const std::string msg = fmt::format(
            "Negative oil rate {}. Stopping iteration", oil_rate);
        this->parent.displayDebugMessage_(msg);
        return true;
    }
    return false;
}

bool
GasLiftSingleWellGeneric::OptimizeState::
checkThpControl()
{
    const int well_index = this->parent.well_state_.index(this->parent.well_name_).value();
    const Well::ProducerCMode& control_mode = this->parent.well_state_.well(well_index).production_cmode;
    return control_mode == Well::ProducerCMode::THP;
}

//
// bool checkEcoGradient(double gradient)
//
//  - Determine if the gradient has reached the limit of the economic gradient.
//
//  - If we are increasing lift gas, returns true if the gradient is smaller
//    than or equal to the economic gradient,
//
//  - If we are decreasing lift gas, returns true if the gradient is greater
//    than or equal to the economic gradient. (I.e., we assume too much lift gas
//    is being used and the gradient has become too small. We try to decrease
//    lift gas until the gradient increases and reaches the economic gradient..)
//
bool
GasLiftSingleWellGeneric::OptimizeState::
checkEcoGradient(double gradient)
{
    std::ostringstream ss;
    bool result = false;

    if (this->parent.debug_) {
        ss << "checking gradient: " << gradient;
    }
    if (this->increase) {
        if (this->parent.debug_) ss << " <= " << this->parent.eco_grad_ << " --> ";
        if (gradient <= this->parent.eco_grad_) {
            if (this->parent.debug_) ss << "yes, stopping";
            result = true;
        }
        else {
            if (this->parent.debug_) ss << "no, continue";
        }
    }
    else {  // decreasing lift gas
        if (this->parent.debug_) ss << " >= " << this->parent.eco_grad_ << " --> ";
        if (gradient >= this->parent.eco_grad_) {
            if (this->parent.debug_) ss << "yes, stopping";
            result = true;
        }
        else {
            if (this->parent.debug_) ss << "no, continue";
        }
    }
    if (this->parent.debug_) this->parent.displayDebugMessage_(ss.str());
    return result;
}

bool
GasLiftSingleWellGeneric::OptimizeState::
checkRate(double rate, double limit, const std::string& rate_str) const
{
    if (limit < rate) {
        if (this->parent.debug_) {
            const std::string msg = fmt::format(
                "iteration {} : {} rate {} exceeds target {}. Stopping iteration",
                this->it, rate_str, rate, limit);
            this->parent.displayDebugMessage_(msg);
        }
        return true;
    }
    return false;
}

bool
GasLiftSingleWellGeneric::OptimizeState::
checkWellRatesViolated(std::vector<double>& potentials)
{
    auto callback = [*this](double rate, double limit, const std::string& rate_str)
                    -> bool
                    { return this->checkRate(rate, limit, rate_str); };
    return this->parent.checkWellRatesViolated_(potentials, callback, this->increase);
}

bool
GasLiftSingleWellGeneric::OptimizeState::
computeBhpAtThpLimit(double alq)
{
    auto bhp_opt = this->parent.computeBhpAtThpLimit_(alq);
    if (bhp_opt) {
        this->bhp = *bhp_opt;
        return true;
    }
    else {
        return false;
    }
}

void
GasLiftSingleWellGeneric::OptimizeState::
debugShowIterationInfo(double alq)
{
    const std::string msg = fmt::format("iteration {}, ALQ = {}", this->it, alq);
    this->parent.displayDebugMessage_(msg);
}


//  NOTE: When calculating the gradient, determine what the well would produce if
//  the lift gas injection rate were increased by one increment. The
//  production rates are adjusted if necessary to obey
//  any rate or BHP limits that the well may be subject to. From this
//  information, calculate the well's "weighted incremental
//  gradient"
//
// TODO: What does it mean to "adjust the production rates" given a
//   BHP limit?
//
double
GasLiftSingleWellGeneric::OptimizeState::
getBhpWithLimit()
{
    auto [new_bhp, limited] = this->parent.getBhpWithLimit_(this->bhp);
    if (limited) {
        // TODO: is it possible that bhp falls below the limit when
        // adding lift gas? I.e. if this->increase == true..
        // TODO: we keep the current alq, but it should probably
        // be adjusted since we changed computed bhp. But how?

        // Stop iteration, but first check the economic gradient
        //   with the bhp_update. If the gradient looks OK (see the
        //   main optimize loop) we keep the current ALQ value.
        this->stop_iteration = true;
    }
    return new_bhp;
}

void
GasLiftSingleWellGeneric::OptimizeState::
updateGroupRates(double delta_oil, double delta_gas, double delta_water, double delta_alq)
{
    const auto &pairs =
        this->parent.group_info_.getWellGroups(this->parent.well_name_);
    for (const auto &[group_name, efficiency] : pairs) {
        int idx = this->parent.group_info_.getGroupIdx(group_name);
        this->parent.sync_groups_.insert(idx);
        this->parent.group_info_.update(group_name,
            efficiency * delta_oil, efficiency * delta_gas, efficiency * delta_water, efficiency * delta_alq);
    }
}


} // namespace Opm
