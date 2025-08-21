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
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <cassert>
#include <sstream>

namespace Opm {

template<typename Scalar, typename IndexTraits>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
GasLiftSingleWellGeneric(DeferredLogger& deferred_logger,
                         WellState<Scalar, IndexTraits>& well_state,
                         const GroupState<Scalar>& group_state,
                         const Well& ecl_well,
                         const SummaryState& summary_state,
                         GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                         const Schedule& schedule,
                         const int report_step_idx,
                         GLiftSyncGroups& sync_groups,
                         const Parallel::Communication& comm,
                         bool glift_debug)
    : GasLiftCommon<Scalar, IndexTraits>(well_state, group_state, deferred_logger, comm, glift_debug)
    , ecl_well_ {ecl_well}
    , summary_state_ {summary_state}
    , group_info_ {group_info}
    , sync_groups_ {sync_groups}
    , controls_ {ecl_well_.productionControls(summary_state_)}
    , debug_limit_increase_decrease_ {false}
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
    assert(this->increment_ > 0);
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
template<typename Scalar, typename IndexTraits>
std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::GradInfo>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
calcIncOrDecGradient(Scalar oil_rate,
                     Scalar gas_rate,
                     Scalar water_rate,
                     Scalar alq,
                     const std::string& gr_name_dont_limit,
                     bool increase,
                     bool debug_output) const
{
    auto [new_alq_opt, alq_is_limited] = addOrSubtractAlqIncrement_(alq, increase);
    // TODO: What to do if ALQ is limited and new_alq != alq?
    if (!new_alq_opt)
        return std::nullopt;

    Scalar new_alq = *new_alq_opt;

    auto delta_alq = new_alq - alq;
    if (checkGroupALQrateExceeded(delta_alq, gr_name_dont_limit))
        return std::nullopt;

    if (auto bhp = computeBhpAtThpLimit_(new_alq, debug_output)) {
        auto [new_bhp, bhp_is_limited] = getBhpWithLimit_(*bhp);
        // TODO: What to do if BHP is limited?
        auto rates = computeWellRates_(new_bhp, bhp_is_limited, debug_output);
        const auto ratesLimited = getLimitedRatesFromRates_(rates);
        BasicRates oldrates = {oil_rate, gas_rate, water_rate, false};
        const auto new_rates = updateRatesToGroupLimits_(oldrates, ratesLimited, gr_name_dont_limit);

        if (!increase && new_rates.oil < 0) {
            return std::nullopt;
        }
        auto grad = calcEcoGradient_(oil_rate, new_rates.oil, gas_rate, new_rates.gas, increase);
        return GradInfo(grad,
                        new_rates.oil,
                        new_rates.oil_is_limited,
                        new_rates.gas,
                        new_rates.gas_is_limited,
                        new_rates.water,
                        new_rates.water_is_limited,
                        new_alq,
                        alq_is_limited);
    } else {
        return std::nullopt;
    }
}

template<typename Scalar, typename IndexTraits>
std::unique_ptr<GasLiftWellState<Scalar>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
runOptimize(const int iteration_idx)
{
    OPM_TIMEFUNCTION();
    std::unique_ptr<GasLiftWellState<Scalar>> state;
    if (this->optimize_) {
        if (this->debug_limit_increase_decrease_) {
            state = runOptimize1_();
        } else {
            state = runOptimize2_();
        }
        if (state) {
            // NOTE: that state->increase() returns a std::optional<bool>, if
            //  this is std::nullopt it means that we was not able to change ALQ
            //  (either increase or decrease)
            if (state->increase()) { // ALQ changed..
                Scalar alq = state->alq();
                if (this->debug)
                    logSuccess_(alq, iteration_idx);
                this->well_state_.well(this->well_name_).alq_state.set(alq);
                const auto& pu = this->well_state_.phaseUsageInfo();

                std::vector<Scalar> well_pot(pu.numActivePhases(), 0.0);
                if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
                    const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
                    well_pot[oil_pos] = state->oilRate();
                }
                if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
                    const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
                    well_pot[water_pos] = state->waterRate();
                }
                if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
                    const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
                    well_pot[gas_pos] = state->gasRate();
                }

                this->well_state_[this->well_name_].well_potentials = well_pot;
            }
        }
    }
    return state;
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
wellTestALQ()
{
    // If WLIFTOPT item 2 is NO we don't optimize
    if (!this->optimize_) {
        return {0.0, false};
    }

    Scalar temp_alq = std::max(this->min_alq_, Scalar(0.0));
    auto cur_alq = temp_alq;
    auto init_rates = computeLimitedWellRatesWithALQ_(temp_alq);
    LimitedRates new_rates = *init_rates;
    std::optional<LimitedRates> rates;
    Scalar old_gradient = 0.0;
    bool alq_is_limited = false;
    bool increase = true;
    OptimizeState state {*this, increase};
    bool success = false;

    while (!state.stop_iteration && (++state.it <= this->max_iterations_)) {
        if (state.checkAlqOutsideLimits(temp_alq, new_rates.oil))
            break;
        std::optional<Scalar> alq_opt;
        std::tie(alq_opt, alq_is_limited) = state.addOrSubtractAlqIncrement(temp_alq);
        if (!alq_opt)
            break;

        temp_alq = *alq_opt;
        if (this->debug)
            state.debugShowIterationInfo(temp_alq);
        rates = new_rates;
        auto temp_rates = computeLimitedWellRatesWithALQ_(temp_alq);
        if (!temp_rates)
            temp_rates->oil = 0.0;
        if (temp_rates->bhp_is_limited)
            state.stop_iteration = true;

        auto gradient = state.calcEcoGradient(rates->oil, temp_rates->oil, rates->gas, temp_rates->gas);
        if (this->debug)
            debugCheckNegativeGradient_(gradient,
                                        cur_alq,
                                        temp_alq,
                                        rates->oil,
                                        temp_rates->oil,
                                        rates->gas,
                                        temp_rates->gas,
                                        increase);
        if (!success && gradient != old_gradient) {
            success = true;
        }
        if (success && gradient <= old_gradient) {
            break;
        }
        cur_alq = temp_alq;
        new_rates = *temp_rates;
        old_gradient = gradient;
    }
    if (state.it > this->max_iterations_) {
        warnMaxIterationsExceeded_();
    }
    if (success) {
        this->well_state_.well(this->well_name_).alq_state.update_count(increase);
    }
    return {cur_alq, success};
}

/****************************************
 * Protected methods in alphabetical order
 ****************************************/

template<typename Scalar, typename IndexTraits>
std::pair<std::optional<Scalar>, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
addOrSubtractAlqIncrement_(Scalar alq, bool increase) const
{
    bool limited = false;
    Scalar orig_alq = alq;
    if (increase) {
        alq += this->increment_;
        // NOTE: if max_alq_ was defaulted in WLIFTOPT, item 3, it has
        //   already been set to the largest value in the VFP table in
        //   the contructor of GasLiftSingleWell
        if (alq > this->max_alq_) {
            alq = this->max_alq_;
            limited = true;
        }
    } else { // we are decreasing ALQ
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
        } else {
            if (alq < 1e-12) {
                alq = 0.0;
                limited = true;
            }
        }
    }
    std::optional<Scalar> alq_opt {alq};
    // If we were not able to change ALQ (to within rounding error), we
    //   return std::nullopt
    if (limited && checkALQequal_(orig_alq, alq))
        alq_opt = std::nullopt;

    return {alq_opt, limited};
}

template<typename Scalar, typename IndexTraits>
Scalar
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
calcEcoGradient_(Scalar oil_rate,
                 Scalar new_oil_rate,
                 Scalar gas_rate,
                 Scalar new_gas_rate,
                 bool increase) const
{
    auto dqo = new_oil_rate - oil_rate;
    auto dqg = new_gas_rate - gas_rate;
    // TODO: Should the gas rate term in the denominator be subject to the constraint:
    //
    //    alpha_g_ * dqg >= 0.0
    //
    //   ?
    auto gradient = (this->alpha_w_ * dqo) / (this->increment_ + this->alpha_g_ * dqg);
    // TODO: Should we do any error checks on the calculation of the
    //   gradient?

    if (!increase)
        gradient = -gradient;
    return gradient;
}

template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::
checkALQequal_(Scalar alq1, Scalar alq2) const
{
    return std::fabs(alq1 - alq2) < (this->increment_ * ALQ_EPSILON);
}

template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::
checkGroupTargetsViolated(const BasicRates& rates,
                          const BasicRates& new_rates) const
{
    const auto& pairs = this->group_info_.getWellGroups(this->well_name_);
    for (const auto& [group_name, efficiency] : pairs) {
        for (const auto rate_type : {Rate::oil, Rate::gas, Rate::water, Rate::liquid}) {
            auto target_opt = this->group_info_.getTarget(rate_type, group_name);
            if (target_opt) {
                auto delta_rate = new_rates[rate_type] - rates[rate_type];
                auto new_group_rate = this->group_info_.getRate(rate_type, group_name) + efficiency * delta_rate;
                if (new_group_rate > *target_opt) {
                    if (this->debug) {
                        const std::string msg
                            = fmt::format("Group {} : {} rate {} exceeds target {}. Stopping iteration",
                                          group_name,
                                          GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(rate_type),
                                          new_group_rate,
                                          *target_opt);
                        displayDebugMessage_(msg);
                    }
                    return true;
                }
            }
        }
    }
    return false;
}

template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::
checkInitialALQmodified_(Scalar alq, Scalar initial_alq) const
{
    if (checkALQequal_(alq, initial_alq)) {
        return false;
    } else {
        const std::string msg = fmt::format("initial ALQ changed from {} "
                                            "to {} before iteration starts..",
                                            initial_alq,
                                            alq);
        displayDebugMessage_(msg);
        return true;
    }
}

template<typename Scalar, typename IndexTraits>
std::pair<std::optional<Scalar>,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
computeConvergedBhpAtThpLimitByMaybeIncreasingALQ_() const
{
    auto alq = this->orig_alq_;
    Scalar new_alq = alq;
    std::optional<Scalar> bhp;
    while ((alq < this->max_alq_) || checkALQequal_(alq, this->max_alq_)) {
        if (bhp = computeBhpAtThpLimit_(alq); bhp) {
            new_alq = alq;
            break;
        }
        alq += this->increment_;
    }
    return {bhp, new_alq};
}

template<typename Scalar, typename IndexTraits>
std::pair<std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::BasicRates>,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::computeInitialWellRates_() const
{
    std::optional<BasicRates> rates;
    Scalar initial_alq = this->orig_alq_;
    // auto alq = initial_alq;
    // if (auto bhp = computeBhpAtThpLimit_(this->orig_alq_); bhp) {
    if (auto [bhp, alq] = computeConvergedBhpAtThpLimitByMaybeIncreasingALQ_(); bhp) {
        {
            const std::string msg = fmt::format("computed initial bhp {} given thp limit and given alq {}", *bhp, alq);
            displayDebugMessage_(msg);
        }
        initial_alq = alq;
        auto [new_bhp, bhp_is_limited] = getBhpWithLimit_(*bhp);
        rates = computeWellRates_(new_bhp, bhp_is_limited);
        if (rates) {
            const std::string msg = fmt::format("computed initial well potentials given bhp, "
                                                "oil: {}, gas: {}, water: {}",
                                                rates->oil,
                                                rates->gas,
                                                rates->water);
            displayDebugMessage_(msg);
        }
    }
    else {
        displayDebugMessage_("Aborting optimization.");
    }
    return {rates, initial_alq};
}

template<typename Scalar, typename IndexTraits>
std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
computeLimitedWellRatesWithALQ_(Scalar alq) const
{
    std::optional<LimitedRates> limited_rates;
    if (auto rates = computeWellRatesWithALQ_(alq); rates) {
        limited_rates = getLimitedRatesFromRates_(*rates);
    }
    return limited_rates;
}

template<typename Scalar, typename IndexTraits>
std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::BasicRates>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
computeWellRatesWithALQ_(Scalar alq) const
{
    OPM_TIMEFUNCTION();
    std::optional<BasicRates> rates;
    auto bhp_opt = computeBhpAtThpLimit_(alq);
    if (bhp_opt) {
        auto [bhp, bhp_is_limited] = getBhpWithLimit_(*bhp_opt);
        rates = computeWellRates_(bhp, bhp_is_limited);
    }
    return rates;
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
debugCheckNegativeGradient_(Scalar grad,
                            Scalar alq,
                            Scalar new_alq,
                            Scalar oil_rate,
                            Scalar new_oil_rate,
                            Scalar gas_rate,
                            Scalar new_gas_rate,
                            bool increase) const
{
    {
        const std::string msg = fmt::format("calculating gradient: "
                                            "new_oil_rate = {}, oil_rate = {}, grad = {}",
                                            new_oil_rate,
                                            oil_rate,
                                            grad);
        displayDebugMessage_(msg);
    }
    if (grad < 0) {
        const std::string msg = fmt::format("negative {} gradient detected ({}) : "
                                            "alq: {}, new_alq: {}, "
                                            "oil_rate: {}, new_oil_rate: {}, gas_rate: {}, new_gas_rate: {}",
                                            (increase ? "incremental" : "decremental"),
                                            grad,
                                            alq,
                                            new_alq,
                                            oil_rate,
                                            new_oil_rate,
                                            gas_rate,
                                            new_gas_rate);
        displayDebugMessage_(msg);
    }
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
debugShowAlqIncreaseDecreaseCounts_()
{
    auto inc_count = this->well_state_.well(this->well_name_).alq_state.get_increment_count();
    auto dec_count = this->well_state_.well(this->well_name_).alq_state.get_decrement_count();
    const std::string msg = fmt::format("ALQ increase/decrease count : {}/{}", inc_count, dec_count);
    displayDebugMessage_(msg);
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
debugShowBhpAlqTable_()
{
    Scalar alq = 0.0;
    constexpr std::string_view fmt_fmt1 {"{:^12s} {:^12s} {:^12s} {:^12s}"};
    constexpr std::string_view fmt_fmt2 {"{:>12.5g} {:>12.5g} {:>12.5g} {:>12.5g}"};
    const std::string header = fmt::format(fmt_fmt1, "ALQ", "BHP", "oil", "gas");
    displayDebugMessage_(header);
    auto max_it = 50;
    auto it = 1;
    while (alq <= (this->max_alq_ + this->increment_)) {
        auto bhp_at_thp_limit = computeBhpAtThpLimit_(alq);
        if (!bhp_at_thp_limit) {
            const std::string msg = fmt::format("Failed to get converged potentials "
                                                "for ALQ = {}. Skipping.",
                                                alq);
            displayDebugMessage_(msg);
        } else {
            auto [bhp, bhp_is_limited] = getBhpWithLimit_(*bhp_at_thp_limit);
            auto rates = computeWellRates_(bhp, bhp_is_limited, /*debug_out=*/false);
            const std::string msg = fmt::format(fmt_fmt2, alq, bhp, rates.oil, rates.gas);
            displayDebugMessage_(msg);
        }
        alq += this->increment_;
        if (it > max_it) {
            const std::string msg = fmt::format("ALQ table : max iterations {} reached. Stopping iteration.", max_it);
            displayDebugMessage_(msg);
            break;
        }
        it++;
    }
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
debugShowLimitingTargets_(const LimitedRates& rates) const
{
    if (rates.limited()) {
        if (rates.oil_is_limited) {
            const std::string msg = fmt::format("oil rate {} is limited by {} target",
                                                rates.oil,
                                                GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(*(rates.oil_limiting_target)));
            displayDebugMessage_(msg);
        }
        if (rates.gas_is_limited) {
            const std::string msg = fmt::format("gas rate {} is limited by GRAT target", rates.gas);
            displayDebugMessage_(msg);
        }
        if (rates.water_is_limited) {
            const std::string msg = fmt::format("water rate {} is limited by {} target",
                                                rates.water,
                                                GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(*(rates.water_limiting_target)));
            displayDebugMessage_(msg);
        }
    } else {
        displayDebugMessage_("no rates are currently limited by a target");
    }
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
debugShowProducerControlMode() const
{
    const int well_index = this->well_state_.index(this->well_name_).value();
    const Well::ProducerCMode& control_mode = this->well_state_.well(well_index).production_cmode;
    const std::string msg = fmt::format("Current control mode is: {}", WellProducerCMode2String(control_mode));
    displayDebugMessage_(msg);
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
debugShowStartIteration_(Scalar alq, bool increase, Scalar oil_rate)
{
    const std::string msg = fmt::format(
        "starting {} iteration, ALQ = {}, oilrate = {}", (increase ? "increase" : "decrease"), alq, oil_rate);
    displayDebugMessage_(msg);
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
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

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
displayDebugMessage_(const std::string& msg) const
{
    if (this->debug) {
        const std::string message = fmt::format("Well {} : {}", this->well_name_, msg);
        this->logMessage_(/*prefix=*/"GLIFT", message);
    }
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
displayWarning_(const std::string& msg)
{
    const std::string message = fmt::format("WELL {} : {}", this->well_name_, msg);
    this->logMessage_(/*prefix=*/"GLIFT", msg, MessageType::WARNING);
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getBhpWithLimit_(Scalar bhp) const
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
template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getGasRateWithLimit_(const BasicRates& rates) const
{
    auto [rate, target_type] = getRateWithLimit_(Rate::gas, rates);
    bool limited = target_type.has_value();
    return {rate, limited};
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
template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getOilRateWithLimit_(const BasicRates& rates) const
{
    auto [rate, target_type] = getRateWithLimit_(Rate::oil, rates);
    bool limited = target_type.has_value();
    return {rate, limited};
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::Rate>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getOilRateWithLimit2_(const BasicRates& rates) const
{
    return getRateWithLimit_(Rate::oil, rates);
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getWaterRateWithLimit_(const BasicRates& rates) const
{
    auto [rate, target_type] = getRateWithLimit_(Rate::water, rates);
    bool limited = target_type.has_value();
    return {rate, limited};
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar,
          std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::Rate>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getWaterRateWithLimit2_(const BasicRates& rates) const
{
    return getRateWithLimit_(Rate::water, rates);
}

template<typename Scalar, typename IndexTraits>
Scalar
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getRate_(Rate rate, const BasicRates& rates) const
{
    switch (rate) {
    case Rate::oil:
        return rates.oil;
    case Rate::gas:
        return rates.gas;
    case Rate::water:
        return rates.water;
    case Rate::liquid:
        return rates.oil + rates.water;
    default:
        // Need this to avoid compiler warning : control reaches end of non-void function
        throw std::runtime_error("This should not happen");
    }
}

template<typename Scalar, typename IndexTraits>
Scalar
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getProductionTarget_(Rate rate) const
{
    switch (rate) {
    case Rate::oil:
        return this->controls_.oil_rate;
    case Rate::gas:
        return this->controls_.gas_rate;
    case Rate::water:
        return this->controls_.water_rate;
    case Rate::liquid:
        return this->controls_.liquid_rate;
    default:
        // Need this to avoid compiler warning : control reaches end of non-void function
        throw std::runtime_error("This should not happen");
    }
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar,
          std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::Rate>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getRateWithLimit_(Rate rate_type, const BasicRates& rates) const
{
    Scalar new_rate = getRate_(rate_type, rates);
    // If "target_type" is empty at the end of this method, it means the rate
    //   was not limited. Otherwise, target_type gives the reason (the type of target)
    //   for why the rate was limited.
    std::optional<Rate> target_type;

    // we also need to limit the other rate (currently only for water and oil and not gas)
    Scalar rate2 = 0.0;
    if (rate_type == Rate::oil) {
        rate2 = getRate_(Rate::water, rates);
    } else if (rate_type == Rate::water) {
        rate2 = getRate_(Rate::oil, rates);
    }

    if (hasProductionControl_(rate_type)) {
        auto target = getProductionTarget_(rate_type);
        if (new_rate > target) {
            const std::string msg = fmt::format("limiting {} rate to target: "
                                                "computed rate: {}, target: {}",
                                                GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(rate_type),
                                                new_rate,
                                                target);
            displayDebugMessage_(msg);
            rate2 *= target/new_rate;
            new_rate = target;
            target_type = rate_type; 
        }
    }
    if (((rate_type == Rate::oil) || (rate_type == Rate::water))) {
        if (rate_type == Rate::oil) {
            if(hasProductionControl_(Rate::water)) {
                auto water_target = getProductionTarget_(Rate::water);
                if (rate2 > water_target) {
                    new_rate *= (water_target / rate2);
                    target_type = Rate::water;
                    rate2 = water_target;
                    const std::string msg = fmt::format("limiting {} rate to {} due to WRAT target: "
                                                        "computed WRAT: {}, target WRAT: {}",
                                                        GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(rate_type),
                                                        new_rate,
                                                        rate2,
                                                        water_target);
                    displayDebugMessage_(msg);
                }
            }
        } else {
            if(hasProductionControl_(Rate::oil)) {
                auto oil_target = getProductionTarget_(Rate::oil);
                if (rate2 > oil_target) {
                    new_rate *= (oil_target / rate2);
                    target_type = Rate::oil;
                    rate2 = oil_target;
                    const std::string msg = fmt::format("limiting {} rate to {} due to ORAT target: "
                                                        "computed ORAT: {}, target ORAT: {}",
                                                        GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(rate_type),
                                                        new_rate,
                                                        rate2,
                                                        oil_target);
                    displayDebugMessage_(msg);
                }
            }
        }

        if(hasProductionControl_(Rate::liquid)) {
            // Note: Since "new_rate" was first updated for ORAT or WRAT, see first "if"
            //   statement in the method, the rate is limited due to LRAT only if
            //   it becomes less than the rate limited by a WRAT or ORAT target..
            Scalar liq_rate = new_rate + rate2;

            auto liq_target = getProductionTarget_(Rate::liquid);
            if (liq_rate > liq_target) {
                Scalar fraction = new_rate / liq_rate;
                // NOTE: since
                //      fraction * liq_rate = new_rate,
                //  we must have
                //      fraction * liq_target < new_rate
                //  since
                //      liq_target < liq_rate
                //  therefore new_rate will become less than it original was and
                //  limited = true.
                new_rate = fraction * liq_target;
                target_type = Rate::liquid;
                const std::string msg = fmt::format("limiting {} rate to {} due to LRAT target: "
                                                    "computed LRAT: {}, target LRAT: {}",
                                                    GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(rate_type),
                                                    new_rate,
                                                    liq_rate,
                                                    liq_target);
                displayDebugMessage_(msg);
            }
        }
    }
    // TODO: Also check RESV target?
    return {new_rate, target_type};
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getOilRateWithGroupLimit_(Scalar new_oil_rate,
                          Scalar oil_rate,
                          const std::string& gr_name_dont_limit) const
{
    [[maybe_unused]] auto [rate, gr_name, efficiency] = getRateWithGroupLimit_(Rate::oil, new_oil_rate, oil_rate, gr_name_dont_limit);
    bool limited = gr_name != nullptr;
    return {rate, limited};
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getGasRateWithGroupLimit_(Scalar new_gas_rate,
                          Scalar gas_rate,
                          const std::string& gr_name_dont_limit) const
{
    [[maybe_unused]] auto [rate, gr_name, efficiency] = getRateWithGroupLimit_(Rate::gas, new_gas_rate, gas_rate, gr_name_dont_limit);
    bool limited = gr_name != nullptr;
    return {rate, limited};
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getWaterRateWithGroupLimit_(Scalar new_water_rate,
                            Scalar water_rate,
                                const std::string& gr_name_dont_limit) const
{
    [[maybe_unused]] auto [rate, gr_name, efficiency] = getRateWithGroupLimit_(Rate::water, new_water_rate, water_rate, gr_name_dont_limit);
    bool limited = gr_name != nullptr;
    return {rate, limited};
}

template<typename Scalar, typename IndexTraits>
std::tuple<Scalar,
           Scalar, bool, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getLiquidRateWithGroupLimit_(const Scalar new_oil_rate,
                             const Scalar oil_rate,
                             const Scalar new_water_rate,
                             const Scalar water_rate,
                             const std::string& gr_name_dont_limit) const
{
    auto liquid_rate = oil_rate + water_rate;
    auto new_liquid_rate = new_oil_rate + new_water_rate;
    auto [liquid_rate_limited, group_name, efficiency]
        = getRateWithGroupLimit_(Rate::liquid, new_liquid_rate, liquid_rate, gr_name_dont_limit);
    bool limited = group_name != nullptr;
    if (limited) {
        Scalar oil_fraction = oil_rate / liquid_rate;
        Scalar delta_liquid = liquid_rate_limited - liquid_rate;
        auto limited_oil_rate = oil_rate + oil_fraction * delta_liquid;
        auto limited_water_rate = water_rate + (1.0 - oil_fraction) * delta_liquid;
        return {limited_oil_rate, limited_water_rate, limited, limited};
    }
    return {new_oil_rate, new_water_rate, limited, limited};
}

template<typename Scalar, typename IndexTraits>
std::tuple<Scalar, const std::string*, Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getRateWithGroupLimit_(Rate rate_type,
                       const Scalar new_rate,
                       const Scalar old_rate,
                       const std::string& gr_name_dont_limit) const
{
    const Scalar delta_rate = new_rate - old_rate;
    if (delta_rate > 0) {
        // It is required that the production rate for a given group is
        // is less than or equal to its target rate.
        // Then it only makes sense to check if the group target is exceeded
        //  if delta_rate > 0
        const auto& pairs = this->group_info_.getWellGroups(this->well_name_);
        Scalar limited_rate = new_rate;
        Scalar gr_target{}, new_gr_rate{}, efficiency{};
        const std::string* group_name = nullptr;
        for (const auto& [group_name_temp, efficiency_temp] : pairs) {
            // in stage 2 we don't want to limit the rate to the group
            // target we are trying to redistribute the gaslift within
            if (gr_name_dont_limit == group_name_temp) {
                continue;
            }

            auto gr_target_opt = this->group_info_.getTarget(rate_type, group_name_temp);
            if (gr_target_opt) {
                Scalar gr_target_temp = *gr_target_opt;
                Scalar gr_rate_temp = this->group_info_.getRate(rate_type, group_name_temp);
                if (gr_rate_temp > gr_target_temp) {
                    if (this->debug) {
                        debugInfoGroupRatesExceedTarget(rate_type, group_name_temp, gr_rate_temp, gr_target_temp);
                    }
                    group_name = &group_name_temp;
                    efficiency = efficiency_temp;
                    limited_rate = old_rate;
                    gr_target = gr_target_temp;
                    new_gr_rate = gr_rate_temp;
                    break;
                }
                Scalar new_gr_rate_temp = gr_rate_temp + efficiency_temp * delta_rate;
                if (new_gr_rate_temp > gr_target_temp) {
                    Scalar limited_rate_temp = old_rate + (gr_target_temp - gr_rate_temp) / efficiency_temp;
                    if (limited_rate_temp < limited_rate) {
                        group_name = &group_name_temp;
                        efficiency = efficiency_temp;
                        limited_rate = limited_rate_temp;
                        gr_target = gr_target_temp;
                        new_gr_rate = new_gr_rate_temp;
                    }
                }
            }
        }
        if (group_name) {
            if (this->debug) {
                const std::string msg = fmt::format("limiting {} rate from {} to {} to meet group target {} "
                                                    "for group {}. Computed group rate was: {}",
                                                    GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(rate_type),
                                                    new_rate,
                                                    limited_rate,
                                                    gr_target,
                                                    *group_name,
                                                    new_gr_rate);
                displayDebugMessage_(msg);
            }
            return {limited_rate, group_name, efficiency};
        }
    }
    return {new_rate, /*group_name =*/nullptr, /*efficiency dummy value*/ 0.0};
}

template<typename Scalar, typename IndexTraits>
std::pair<std::optional<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates>,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getInitialRatesWithLimit_() const
{
    std::optional<LimitedRates> limited_rates;
    Scalar initial_alq = this->orig_alq_;
    if (auto [rates, alq] = computeInitialWellRates_(); rates) {
        if (this->debug) {
            displayDebugMessage_("Maybe limiting initial rates before optimize loop..");
        }
        auto temp_rates = getLimitedRatesFromRates_(*rates);
        BasicRates old_rates = getWellStateRates_();
        limited_rates = updateRatesToGroupLimits_(old_rates, temp_rates);

        initial_alq = alq;
    }
    return {limited_rates, initial_alq};
}

template<typename Scalar, typename IndexTraits>
typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getLimitedRatesFromRates_(const BasicRates& rates) const
{
    auto [oil_rate, oil_limiting_target] = getOilRateWithLimit2_(rates);
    bool oil_is_limited = oil_limiting_target.has_value();
    auto [gas_rate, gas_is_limited] = getGasRateWithLimit_(rates);
    auto [water_rate, water_limiting_target] = getWaterRateWithLimit2_(rates);
   
    bool water_is_limited = water_limiting_target.has_value();
    return LimitedRates {oil_rate,
                         gas_rate,
                         water_rate,
                         oil_is_limited,
                         gas_is_limited,
                         water_is_limited,
                         rates.bhp_is_limited,
                         oil_limiting_target,
                         water_limiting_target};
}

template<typename Scalar, typename IndexTraits>
typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::BasicRates
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
getWellStateRates_() const
{
    const int well_index = this->well_state_.index(this->well_name_).value();
    const auto& ws = this->well_state_.well(well_index);
    const auto& wrate = ws.well_potentials;
    const auto& pu = this->well_state_.phaseUsageInfo();

    const auto oil_rate = pu.phaseIsActive(IndexTraits::oilPhaseIdx) ?
                         wrate[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] : Scalar{0.0};

    const auto gas_rate = pu.phaseIsActive(IndexTraits::gasPhaseIdx) ?
                          wrate[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] : Scalar{0.0};

    const auto water_rate = pu.phaseIsActive(IndexTraits::waterPhaseIdx) ?
                          wrate[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] : Scalar{0.0};
    if (this->debug) {
        const std::string msg = fmt::format("Initial surface rates: oil : {}, "
                                            "gas : {}, water : {}",
                                            oil_rate,
                                            gas_rate,
                                            water_rate);
        displayDebugMessage_(msg);
    }
    return BasicRates {oil_rate, gas_rate, water_rate, /*bhp_is_limited=*/false};
}

template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::
hasProductionControl_(Rate rate) const
{
    switch (rate) {
    case Rate::oil:
        return this->controls_.hasControl(Well::ProducerCMode::ORAT);
    case Rate::gas:
        return this->controls_.hasControl(Well::ProducerCMode::GRAT);
    case Rate::water:
        return this->controls_.hasControl(Well::ProducerCMode::WRAT);
    case Rate::liquid:
        return this->controls_.hasControl(Well::ProducerCMode::LRAT);
    default:
        // Need this to avoid compiler warning : control reaches end of non-void function
        throw std::runtime_error("This should not happen");
    }
}

template<typename Scalar, typename IndexTraits>
std::pair<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
increaseALQtoPositiveOilRate_(Scalar alq,
                              const LimitedRates& orig_rates) const
{
    Scalar temp_alq = alq;
    // use the copy constructor to only copy the rates
    BasicRates rates = orig_rates;
    while (true) {
        temp_alq += this->increment_;
        if (temp_alq > this->max_alq_)
            break;
        auto temp_rates = computeWellRatesWithALQ_(temp_alq);
        if (!temp_rates)
            break;
        alq = temp_alq;
        rates = *temp_rates;
        if (rates.oil > 0)
            break;
    }
    // TODO: what about group limits?
    return {getLimitedRatesFromRates_(rates), alq};
}

template<typename Scalar, typename IndexTraits>
std::pair<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
increaseALQtoMinALQ_(const Scalar orig_alq,
                    const LimitedRates& orig_rates) const
{
    auto min_alq = this->min_alq_;
    assert(min_alq >= 0);
    assert(orig_alq < min_alq);
    assert(min_alq < this->max_alq_);
    Scalar alq = orig_alq;
    LimitedRates rates = orig_rates;
    while (true) {
        Scalar temp_alq = alq + this->increment_;

        alq = temp_alq;
        if (temp_alq >= min_alq)
            break;

        auto temp_rates = computeLimitedWellRatesWithALQ_(temp_alq);

        if (temp_rates) {
            rates = *temp_rates;
            if (rates.limited())
                break;
        }
    }
    return std::make_pair(rates, alq);
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
logSuccess_(Scalar alq, const int iteration_idx)
{
    const std::string message = fmt::format("GLIFT, IT={}, WELL {} : {} ALQ from {} to {}",
                                            iteration_idx,
                                            this->well_name_,
                                            ((alq > this->orig_alq_) ? "increased" : "decreased"),
                                            this->orig_alq_,
                                            alq);
    this->deferred_logger_.debug(message);
}

template<typename Scalar, typename IndexTraits>
std::pair<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
maybeAdjustALQbeforeOptimizeLoop_(const LimitedRates& orig_rates,
                                  const Scalar orig_alq,
                                  const bool increase) const
{
    Scalar alq = orig_alq;
    LimitedRates rates = orig_rates;

    if (this->debug) {
        const std::string msg = fmt::format("initial ALQ: {}", alq);
        displayDebugMessage_(msg);
    }
    if (!increase) {
        // NOTE: Try to decrease ALQ down to a value where the groups
        // maximum alq target and the total gas + alq target is not violated
        std::tie(rates, alq) = reduceALQtoGroupAlqLimits_(alq, orig_rates);
        if (orig_rates.limited()) {
            // NOTE: Try to decrease ALQ down to a value where the well target is
            //   not exceeded.
            // NOTE: This may reduce ALQ below the minimum value set in WLIFTOPT
            //   item 5. However, this is OK since the rate target is met and there
            //   is no point in using a higher ALQ value then.
            auto [rates1, alq1] = reduceALQtoWellTarget_(alq, orig_rates);
            auto [rates2, alq2] = reduceALQtoGroupTarget(alq, orig_rates);
            if (alq1 < alq2) {
                alq = alq1;
                rates = rates1;
            } else {
                alq = alq2;
                rates = rates2;
            }
        }
    } else {
        if (orig_rates.oil < 0) {
            // Try to increase ALQ up to a value where oil_rate is positive
            std::tie(rates, alq) = increaseALQtoPositiveOilRate_(alq, rates);
        }
        if ((this->min_alq_ > 0) && (alq < this->min_alq_)) {
            // Try to increase ALQ up to the minimum limit without checking
            //   the economic gradient..
            std::tie(rates, alq) = increaseALQtoMinALQ_(alq, rates);
        }
    }
    if (orig_alq != alq) {
        if (this->debug) {
            const std::string msg = fmt::format("adjusted ALQ to: {}", alq);
            displayDebugMessage_(msg);
        }
    }
    return {rates, alq};
}

// Reduce ALQ to the lowest value greater than zero that still makes at
//   least one rate limited w.r.t. group targets, or reduce ALQ to zero if
//   such positive ALQ value cannot be found.
template<typename Scalar, typename IndexTraits>
std::pair<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
reduceALQtoGroupAlqLimits_(const Scalar orig_alq,
                           const LimitedRates& orig_rates) const
{
    Scalar alq = orig_alq;
    BasicRates rates {orig_rates};
    Scalar temp_alq = orig_alq;
    while (true) {
        if (temp_alq == 0)
            break;
        temp_alq -= this->increment_;
        if (temp_alq < 0)
            temp_alq = 0;
        auto new_rates = computeWellRatesWithALQ_(temp_alq);
        if (!new_rates)
            break;
        auto delta_alq = temp_alq - orig_alq;
        auto delta_gas_rate = new_rates->gas - orig_rates.gas;
        if (!checkGroupTotalRateExceeded(delta_alq, delta_gas_rate)) {
            break;
        }
        rates = *new_rates;
        alq = temp_alq;
    }
    if (alq == orig_alq) {
        return {orig_rates, orig_alq};
    } else {
        LimitedRates limited_rates = getLimitedRatesFromRates_(rates);
        return {limited_rates, alq};
    }
}
// Reduce ALQ to the lowest value greater than zero that still makes at
//   least one rate limited w.r.t. group targets, or reduce ALQ to zero if
//   such positive ALQ value cannot be found.
template<typename Scalar, typename IndexTraits>
std::pair<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
reduceALQtoGroupTarget(const Scalar orig_alq,
                      const LimitedRates& orig_rates) const
{
    bool stop_this_iteration = true;
    const std::vector<std::pair<std::string, Scalar>>& pairs = this->group_info_.getWellGroups(this->well_name_);
    for (const auto& pair /*<group_name, efficiency>*/ : pairs) {
        const auto& group_name = pair.first;
        if (!this->group_state_.has_production_control(group_name))
            continue;
        if (this->group_info_.hasAnyTarget(group_name)) {
            stop_this_iteration = false;
            displayDebugMessage_("Reducing ALQ to meet group target(s) before iteration starts.");
            break;
        }
    }
    Scalar alq = orig_alq;
    BasicRates rates {orig_rates};
    Scalar temp_alq = orig_alq;
    while (!stop_this_iteration) {
        if (temp_alq == 0)
            break;
        temp_alq -= this->increment_;
        if (temp_alq < 0)
            temp_alq = 0;
        auto new_rates = computeWellRatesWithALQ_(temp_alq);
        if (!new_rates)
            break;
        if (!checkGroupTargetsViolated(rates, *new_rates)) {
            break;
        }
        rates = *new_rates;
        alq = temp_alq;
    }
    if (alq == orig_alq) {
        return {orig_rates, orig_alq};
    } else {
        LimitedRates limited_rates = getLimitedRatesFromRates_(rates);
        return {limited_rates, alq};
    }
}

// Reduce ALQ to the lowest value greater than zero that still makes at
//   least one rate limited w.r.t. well targets, or reduce ALQ to zero if
//   such positive ALQ value cannot be found.
template<typename Scalar, typename IndexTraits>
std::pair<typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates,
          Scalar>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
reduceALQtoWellTarget_(const Scalar orig_alq,
                       const LimitedRates& rates) const
{
    // this method should only be called if "rates" is limited
    assert(rates.limited());
    if (this->debug) {
        displayDebugMessage_("Reducing ALQ to meet well targets before iteration starts..");
        debugShowLimitingTargets_(rates);
    }
    Scalar alq = orig_alq;
    Scalar temp_alq = alq;
    std::optional<LimitedRates> new_rates;
    while (true) {
        if (temp_alq == 0)
            break;
        temp_alq -= this->increment_;
        if (temp_alq < 0)
            temp_alq = 0;
        auto temp_rates = computeLimitedWellRatesWithALQ_(temp_alq);
        if (!temp_rates)
            break; // failed to compute BHP given THP limit and ALQ
        // keep iterating until no rate is limited
        if (!temp_rates->limited())
            break;
        alq = temp_alq;
        new_rates = temp_rates;
    }
    assert(alq <= orig_alq);
    if (this->debug) {
        if (alq < orig_alq) {
            // NOTE: ALQ may drop below zero before we are able to meet the target
            const std::string msg = fmt::format("Reduced ALQ from {} to {} to meet rate targets. Rates (new, old) : "
                                                "oil(({}, {}), gas({}, {}), water({}, {})",
                                                orig_alq,
                                                alq,
                                                new_rates->oil,
                                                rates.oil,
                                                new_rates->gas,
                                                rates.gas,
                                                new_rates->water,
                                                rates.water);
            displayDebugMessage_(msg);
        } else if (alq == orig_alq) {
            // We might not be able to reduce ALQ, for example if ALQ starts out at zero.
            const std::string msg = fmt::format("Not able to reduce ALQ {} further. ", orig_alq);
            displayDebugMessage_(msg);
        }
    }
    if (new_rates) {
        return {*new_rates, alq};
    } else {
        return {rates, orig_alq};
    }
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
template<typename Scalar, typename IndexTraits>
std::unique_ptr<GasLiftWellState<Scalar>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
runOptimizeLoop_(bool increase)
{
    OPM_TIMEFUNCTION();
    if (this->debug)
        debugShowProducerControlMode();
    std::unique_ptr<GasLiftWellState<Scalar>> ret_value; // nullptr initially
    auto [rates, cur_alq] = getInitialRatesWithLimit_();
    if (!rates)
        return ret_value;
    // if (this->debug) debugShowBhpAlqTable_();
    if (this->debug)
        debugShowAlqIncreaseDecreaseCounts_();
    if (this->debug)
        debugShowTargets_();
    bool success = false; // did we succeed to increase alq?
    bool alq_is_limited = false;
    LimitedRates new_rates = *rates;
    auto [temp_rates2, new_alq] = maybeAdjustALQbeforeOptimizeLoop_(*rates, cur_alq, increase);
    if (checkInitialALQmodified_(new_alq, this->orig_alq_)) {
        auto delta_alq = new_alq - cur_alq;
        new_rates = temp_rates2;
        cur_alq = new_alq;
        success = true;
        updateGroupRates_(*rates, new_rates, delta_alq);
    }

    OptimizeState state {*this, increase};
    auto temp_alq = cur_alq;
    if (checkThpControl_()) {
        if (this->debug)
            debugShowStartIteration_(temp_alq, increase, new_rates.oil);
    } else {
        // If the well is not under THP control, we can still use the previous
        //   initial adjustment of ALQ by using the well's THP limit to calculate
        //   BHP and then well rates from that.
        //   This is useful for example for wells under group control to reduce
        //   their gaslift. A typical case for this could be that a new well opens up.
        //   Then gaslift can be reduced while still keeping the group target.
        state.stop_iteration = true;
    }
    while (!state.stop_iteration && (++state.it <= this->max_iterations_)) {
        OPM_TIMEBLOCK(GasliftOptimizationLoopStage1);
        if (state.checkRatesViolated(new_rates))
            break;
        if (state.checkAlqOutsideLimits(temp_alq, new_rates.oil))
            break;
        std::optional<Scalar> alq_opt;
        std::tie(alq_opt, alq_is_limited) = state.addOrSubtractAlqIncrement(temp_alq);
        if (!alq_opt)
            break;
        auto delta_alq = *alq_opt - temp_alq;
        if (checkGroupALQrateExceeded(delta_alq))
            break;

        temp_alq = *alq_opt;
        if (this->debug)
            state.debugShowIterationInfo(temp_alq);
        rates = new_rates;
        auto temp_rates = computeLimitedWellRatesWithALQ_(temp_alq);
        if (!temp_rates)
            break;
        if (temp_rates->bhp_is_limited)
            state.stop_iteration = true;
        temp_rates = updateRatesToGroupLimits_(*rates, *temp_rates);

        auto delta_gas_rate = temp_rates->gas - rates->gas;
        if (checkGroupTotalRateExceeded(delta_alq, delta_gas_rate))
            break;

        /*        if (this->debug_abort_if_increase_and_gas_is_limited_) {
                    if (gas_is_limited && increase) {
                        // if gas is limited we do not want to increase
                        displayDebugMessage_(
                            "increasing ALQ and gas is limited -> aborting iteration");
                        break;
                    }
                }
        */
        auto gradient = state.calcEcoGradient(rates->oil, temp_rates->oil, rates->gas, temp_rates->gas);
        if (this->debug)
            debugCheckNegativeGradient_(gradient,
                                        cur_alq,
                                        temp_alq,
                                        rates->oil,
                                        temp_rates->oil,
                                        rates->gas,
                                        temp_rates->gas,
                                        increase);
        if (state.checkEcoGradient(gradient))
            break;
        cur_alq = temp_alq;
        success = true;
        new_rates = *temp_rates;
        updateGroupRates_(*rates, new_rates, delta_alq);
    }
    if (state.it > this->max_iterations_) {
        warnMaxIterationsExceeded_();
    }
    std::optional<bool> increase_opt;
    if (success) {
        this->well_state_.well(this->well_name_).alq_state.update_count(increase);
        increase_opt = increase;
    } else {
        increase_opt = std::nullopt;
    }
    ret_value = std::make_unique<GasLiftWellState<Scalar>>(new_rates.oil,
                                                           new_rates.oil_is_limited,
                                                           new_rates.gas,
                                                           new_rates.gas_is_limited,
                                                           cur_alq,
                                                           alq_is_limited,
                                                           new_rates.water,
                                                           new_rates.water_is_limited,
                                                           increase_opt);
    return ret_value;
}

template<typename Scalar, typename IndexTraits>
std::unique_ptr<GasLiftWellState<Scalar>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::runOptimize1_()
{
    OPM_TIMEFUNCTION();
    std::unique_ptr<GasLiftWellState<Scalar>> state;
    const int inc_count = this->well_state_.well(this->well_name_).alq_state.get_increment_count();
    const int dec_count = this->well_state_.well(this->well_name_).alq_state.get_decrement_count();
    if (dec_count == 0 && inc_count == 0) {
        state = tryIncreaseLiftGas_();
        if (!state || !(state->alqChanged())) {
            state = tryDecreaseLiftGas_();
        }
    } else if (dec_count == 0) {
        assert(inc_count > 0);
        state = tryIncreaseLiftGas_();
    } else if (inc_count == 0) {
        assert(dec_count > 0);
        state = tryDecreaseLiftGas_();
    }
    return state;
}

template<typename Scalar, typename IndexTraits>
std::unique_ptr<GasLiftWellState<Scalar>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::runOptimize2_()
{
    OPM_TIMEFUNCTION();
    std::unique_ptr<GasLiftWellState<Scalar>> state;
    state = tryIncreaseLiftGas_();
    if (!state || !(state->alqChanged())) {
        state = tryDecreaseLiftGas_();
    }
    return state;
}

template<typename Scalar, typename IndexTraits>
std::unique_ptr<GasLiftWellState<Scalar>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::tryDecreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/false);
}

template<typename Scalar, typename IndexTraits>
std::unique_ptr<GasLiftWellState<Scalar>>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::tryIncreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/true);
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
setAlqMinRate_(const GasLiftWell& well)
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

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
updateGroupRates_(const LimitedRates& rates,
                  const LimitedRates& new_rates,
                  Scalar delta_alq) const
{
    Scalar delta_oil = new_rates.oil - rates.oil;
    Scalar delta_gas = new_rates.gas - rates.gas;
    Scalar delta_water = new_rates.water - rates.water;
    const auto& pairs = this->group_info_.getWellGroups(this->well_name_);
    for (const auto& [group_name, efficiency] : pairs) {
        int idx = this->group_info_.getGroupIdx(group_name);
        // This will notify the optimize loop in BlackoilWellModel, see
        //   gasLiftOptimizationStage1() in BlackoilWellModel_impl.hpp
        // that this group_info needs to be synchronized to the other MPI ranks
        this->sync_groups_.insert(idx);
        this->group_info_.update(group_name,
                                 efficiency * delta_oil,
                                 efficiency * delta_gas,
                                 efficiency * delta_water,
                                 efficiency * delta_alq);
    }
}

template<typename Scalar, typename IndexTraits>
typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::LimitedRates
GasLiftSingleWellGeneric<Scalar, IndexTraits>::
updateRatesToGroupLimits_(const BasicRates& old_rates,
                          const LimitedRates& rates,
                          const std::string& gr_name) const
{
    LimitedRates new_rates = rates;
    auto [new_oil_rate, oil_is_limited] = getOilRateWithGroupLimit_(new_rates.oil, old_rates.oil, gr_name);
    auto mod_water_rate = new_rates.water;
    if (oil_is_limited) {
        new_rates.oil_limiting_target = Rate::oil;
        mod_water_rate *=  (new_oil_rate / new_rates.oil);
    }
    auto [new_gas_rate, gas_is_limited] = getGasRateWithGroupLimit_(new_rates.gas, old_rates.gas, gr_name);
    auto [new_water_rate, water_is_limited] = getWaterRateWithGroupLimit_(mod_water_rate, old_rates.water, gr_name);
    if (water_is_limited) {
        new_rates.water_limiting_target = Rate::water;
        new_oil_rate *=  (new_water_rate / new_rates.water);
    }
    auto [new_oil_rate2, new_water_rate2, oil_is_limited2, water_is_limited2]
        = getLiquidRateWithGroupLimit_(new_oil_rate, old_rates.oil, new_water_rate, old_rates.water, gr_name);
    if (oil_is_limited2) {
        new_rates.oil_limiting_target = Rate::liquid;
    }
    if (water_is_limited2) {
        new_rates.water_limiting_target = Rate::liquid;
    }
    new_rates.oil = new_oil_rate2;
    new_rates.gas = new_gas_rate;
    new_rates.water = new_water_rate2;
    new_rates.oil_is_limited = rates.oil_is_limited || oil_is_limited || oil_is_limited2;
    new_rates.gas_is_limited = rates.gas_is_limited || gas_is_limited;
    new_rates.water_is_limited = rates.water_is_limited || water_is_limited || water_is_limited2;
    if (oil_is_limited || oil_is_limited2 || gas_is_limited || water_is_limited || water_is_limited2) {
        new_rates.limit_type = LimitedRates::LimitType::group;
    }
    return new_rates;
}

// Called when we should use a fixed ALQ value
template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
updateWellStateAlqFixedValue_(const GasLiftWell& well)
{
    const auto& max_alq_optional = well.max_rate();
    if (max_alq_optional) {
        // According to WLIFTOPT, item 3:
        // If item 2 is NO, then item 3 is regarded as the fixed
        // lift gas injection rate for the well.
        auto new_alq = *max_alq_optional;
        this->well_state_.well(this->well_name_).alq_state.set(new_alq);
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
template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::
useFixedAlq_(const GasLiftWell& well)
{
    auto wliftopt_item2 = well.use_glo();
    if (wliftopt_item2) {
        return false;
    } else {
        displayDebugMessage_("WLIFTOPT item2 = NO. Skipping optimization.");
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

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::
debugInfoGroupRatesExceedTarget(Rate rate_type,
                                const std::string& gr_name,
                                Scalar rate,
                                Scalar target) const
{
    const std::string msg = fmt::format("{} rate for group {} exceeds target: "
                                        "rate = {}, target = {}, the old rate is kept.",
                                        GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(rate_type),
                                        gr_name,
                                        rate,
                                        target);
    displayDebugMessage_(msg);
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::warnMaxIterationsExceeded_()
{
    const std::string msg = fmt::format("Max iterations ({}) exceeded", this->max_iterations_);
    displayWarning_(msg);
}

/****************************************
 * Methods declared in OptimizeState
 ****************************************/

template<typename Scalar, typename IndexTraits>
std::pair<std::optional<Scalar>, bool>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::OptimizeState::
addOrSubtractAlqIncrement(Scalar alq)
{
    auto [alq_opt, limited] = this->parent.addOrSubtractAlqIncrement_(alq, this->increase);
    if (!alq_opt) {
        const std::string msg = fmt::format("iteration {}, alq = {} : not able to {} ALQ increment",
                                            this->it,
                                            alq,
                                            (this->increase ? "add" : "subtract"));
    }
    return {alq_opt, limited};
}

template<typename Scalar, typename IndexTraits>
Scalar
GasLiftSingleWellGeneric<Scalar, IndexTraits>::OptimizeState::
calcEcoGradient(Scalar oil_rate,
                Scalar new_oil_rate,
                Scalar gas_rate,
                Scalar new_gas_rate)
{
    return this->parent.calcEcoGradient_(oil_rate, new_oil_rate, gas_rate, new_gas_rate, this->increase);
}

// NOTE:  According to WLIFTOPT item 5 :
//   if min_rate() is negative, it means: allocate at least enough lift gas
//   to enable the well to flow
//  We will interpret this as (see discussion above GasLiftSingleWell()
//   in this file): Allocate at least the amount of lift gas needed to
//   get a positive oil production rate.
template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::OptimizeState::
checkAlqOutsideLimits(Scalar alq, [[maybe_unused]] Scalar oil_rate)
{
    std::ostringstream ss;
    bool result = false;

    if (this->increase) {
        if (alq >= this->parent.max_alq_) {
            ss << "ALQ >= " << this->parent.max_alq_ << " (max limit), "
               << "stopping iteration";
            result = true;
        } else { // checking the minimum limit...
            // NOTE: A negative min_alq_ means: allocate at least enough lift gas
            //  to enable the well to flow, see WLIFTOPT item 5.
            if (this->parent.min_alq_ < 0) {
                // - if oil rate is negative (i.e. the well is not flowing), continue to
                //    increase ALQ (according WLIFTOPT item 5) and try make the well
                //    flow.
                // - else if oil rate is already positive, there is no minimum
                //    limit for ALQ in this case
                result = false;
            } else {
                // NOTE: checking for a lower limit is not necessary
                //   when increasing alq. If ALQ was smaller than the minimum when
                //   we entered the runOptimizeLoop_() method,
                //   increaseALQtoMinALQ_() will ensure that ALQ >= min_alq
                assert(alq >= this->parent.min_alq_);
                result = false;
            }
        }
    } else { // we are decreasing lift gas
        if (alq == 0) {
            ss << "ALQ is zero, cannot decrease further. Stopping iteration. ";
        } else if (alq < 0) {
            ss << "Negative ALQ: " << alq << ". Stopping iteration. ";
        }
        // NOTE: A negative min_alq_ means: allocate at least enough lift gas
        //  to enable the well to flow, see WLIFTOPT item 5.
        if (this->parent.min_alq_ < 0) {
            // We know that the well is flowing (oil_rate > 0) since that was
            //  already checked in runOptimizeLoop_() by calling checkNegativeOilRate()
            assert(oil_rate >= 0);
            result = false;
        } else {
            if (alq <= this->parent.min_alq_) {
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
                //    checkRatesViolated().
                // - We also know that the rate limit was not exceeded since that was
                //    checked by checkRatesViolated()
                assert(oil_rate >= 0);
                ss << "ALQ <= " << this->parent.min_alq_ << " (min limit), "
                   << "stopping iteration";
                result = true;
            } else {
                // NOTE: checking for an upper limit should not be necessary
                // when decreasing alq.. so this is just to catch an
                // illegal state at an early point.
                if (this->parent.checkALQequal_(alq, this->parent.max_alq_)) {
                    return false;
                } else if (alq > this->parent.max_alq_) {
                    warn_("unexpected: alq above upper limit when trying to "
                          "decrease lift gas. aborting iteration.");
                    result = true;
                } else {
                    result = false;
                }
            }
        }
    }
    if (this->parent.debug) {
        const std::string msg = ss.str();
        if (!msg.empty())
            this->parent.displayDebugMessage_(msg);
    }
    return result;
}

template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::
checkGroupALQrateExceeded(Scalar delta_alq,
                          const std::string& gr_name_dont_limit) const
{
    const auto& pairs = group_info_.getWellGroups(well_name_);
    for (const auto& [group_name, efficiency] : pairs) {
        // in stage 2 we don't want to limit the rate to the group
        // target we are trying to redistribute the gaslift within
        if (gr_name_dont_limit == group_name)
            continue;
        auto max_alq_opt = group_info_.maxAlq(group_name);
        if (max_alq_opt) {
            Scalar alq = group_info_.alqRate(group_name) + efficiency * delta_alq;
            if (alq > *max_alq_opt) {
                if (this->debug) {
                    const std::string msg = fmt::format(
                        "Group {} : alq {} exceeds max_alq {}. Stopping iteration", group_name, alq, *max_alq_opt);
                    displayDebugMessage_(msg);
                }
                return true;
            }
        }
    }
    return false;
}

template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::
checkGroupTotalRateExceeded(Scalar delta_alq,
                            Scalar delta_gas_rate) const
{
    const auto& pairs = group_info_.getWellGroups(well_name_);
    for (const auto& [group_name, efficiency] : pairs) {
        auto max_total_rate_opt = group_info_.maxTotalGasRate(group_name);
        if (max_total_rate_opt) {
            Scalar alq = group_info_.alqRate(group_name) + efficiency * delta_alq;
            Scalar gas_rate = group_info_.gasRate(group_name) + efficiency * delta_gas_rate;

            if ((alq + gas_rate) > *max_total_rate_opt) {
                if (this->debug) {
                    const std::string msg
                        = fmt::format("Group {} : total gas rate {} exceeds max_total_gas_rate {}. Stopping iteration",
                                      group_name,
                                      alq + gas_rate,
                                      *max_total_rate_opt);
                    displayDebugMessage_(msg);
                }
                return true;
            }
        }
    }
    return false;
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
template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::OptimizeState::
checkEcoGradient(Scalar gradient)
{
    std::ostringstream ss;
    bool result = false;

    if (this->parent.debug) {
        ss << "checking gradient: " << gradient;
    }
    if (this->increase) {
        if (this->parent.debug)
            ss << " <= " << this->parent.eco_grad_ << " --> ";
        if (gradient <= this->parent.eco_grad_) {
            if (this->parent.debug)
                ss << "yes, stopping";
            result = true;
        } else {
            if (this->parent.debug)
                ss << "no, continue";
        }
    } else { // decreasing lift gas
        if (this->parent.debug)
            ss << " >= " << this->parent.eco_grad_ << " --> ";
        if (gradient >= this->parent.eco_grad_) {
            if (this->parent.debug)
                ss << "yes, stopping";
            result = true;
        } else {
            if (this->parent.debug)
                ss << "no, continue";
        }
    }
    if (this->parent.debug)
        this->parent.displayDebugMessage_(ss.str());
    return result;
}

template<typename Scalar, typename IndexTraits>
bool GasLiftSingleWellGeneric<Scalar, IndexTraits>::OptimizeState::
checkRatesViolated(const LimitedRates& rates) const
{
    if (!this->increase) {
        if (rates.oil < 0) {
            // The well is not flowing, and it will(?) not help to reduce lift
            // gas further. Note that this assumes that the oil rates drops with
            // decreasing lift gas.
            this->parent.displayDebugMessage_("Negative oil rate detected while descreasing "
                                              "lift gas. Stopping iteration.");
            return true;
        }
    }
    if (rates.limited()) {
        if (this->parent.debug) {
            const std::string well_or_group = rates.limit_type == LimitedRates::LimitType::well ? "well" : "group";
            std::string target_type;
            std::string rate_type;
            if (rates.oil_is_limited) {
                target_type = GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(*(rates.oil_limiting_target));
                rate_type = "oil";
            } else if (rates.gas_is_limited) {
                target_type = "gas";
                rate_type = "gas";
            } else if (rates.water_is_limited) {
                target_type = GasLiftGroupInfo<Scalar, IndexTraits>::rateToString(*(rates.water_limiting_target));
                rate_type = "water";
            }
            const std::string msg = fmt::format("iteration {} : {} rate was limited due to {} {} target. "
                                                "Stopping iteration",
                                                this->it,
                                                rate_type,
                                                well_or_group,
                                                target_type);
            this->parent.displayDebugMessage_(msg);
        }
        return true;
    }
    return false;
}

template<typename Scalar, typename IndexTraits>
void GasLiftSingleWellGeneric<Scalar, IndexTraits>::OptimizeState::
debugShowIterationInfo(Scalar alq)
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
template<typename Scalar, typename IndexTraits>
Scalar
GasLiftSingleWellGeneric<Scalar, IndexTraits>::OptimizeState::
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

/****************************************
 * Methods declared in BasicRates
 ****************************************/

template<typename Scalar, typename IndexTraits>
GasLiftSingleWellGeneric<Scalar, IndexTraits>::BasicRates::
BasicRates(const LimitedRates& rates)
{
    oil = rates.oil;
    gas = rates.gas;
    water = rates.water;
    bhp_is_limited = rates.bhp_is_limited;
}

template class GasLiftSingleWellGeneric<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class GasLiftSingleWellGeneric<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
