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

#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/GasLiftOpt.hpp>

#include <optional>
#include <string>

template<typename TypeTag>
GasLiftSingleWell<TypeTag>::
GasLiftSingleWell(
    const StdWell &std_well,
    const Simulator &ebos_simulator,
    const SummaryState &summary_state,
    DeferredLogger &deferred_logger,
    WellState &well_state
) :
    deferred_logger_{deferred_logger},
    ebos_simulator_{ebos_simulator},
    std_well_{std_well},
    summary_state_{summary_state},
    well_state_{well_state},
    ecl_well_{std_well_.wellEcl()},
    controls_{ecl_well_.productionControls(summary_state_)},
    num_phases_{well_state_.numPhases()},
    debug{false},  // extra debugging output
    debug_limit_increase_decrease_{false}
{
    const Schedule& schedule = this->ebos_simulator_.vanguard().schedule();
    const int report_step_idx = this->ebos_simulator_.episodeIndex();
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
    auto& gl_well = glo.well(this->well_name_);

    if(useFixedAlq_(gl_well)) {
        updateWellStateAlqFixedValue_(gl_well);
        this->optimize_ = false; // lift gas supply is fixed
    }
    else {
        setAlqMaxRate_(gl_well);
        this->optimize_ = true;
    }
    const auto& pu = std_well_.phaseUsage();
    this->oil_pos_ = pu.phase_pos[Oil];
    this->gas_pos_ = pu.phase_pos[Gas];
    this->water_pos_ = pu.phase_pos[Water];
    // get the alq value used for this well for the previous iteration (a
    //   nonlinear iteration in assemble() in BlackoilWellModel).
    //   If gas lift optimization has not been applied to this well yet, the
    //   default value is used.
    this->orig_alq_ = this->well_state_.getALQ(this->well_name_);
    if(this->optimize_) {
        setAlqMinRate_(gl_well);
        // NOTE: According to item 4 in WLIFTOPT, this value does not
        //    have to be positive.
        // TODO: Does it make sense to have a negative value?
        this->alpha_w_ = gl_well.weight_factor();
        if (this->alpha_w_ <= 0 ) {
            displayWarning_("Nonpositive value for alpha_w ignored");
            this->alpha_w_ = 1.0;
        }

        // NOTE: According to item 6 in WLIFTOPT:
        //   "If this value is greater than zero, the incremental gas rate will influence
        //    the calculation of the incremental gradient and may be used
        //    to discourage the allocation of lift gas to wells which produce more gas."
        // TODO: Does this mean that we should ignore this value if it
        //   is negative?
        this->alpha_g_ = gl_well.inc_weight_factor();

        // TODO: adhoc value.. Should we keep max_iterations_ as a safety measure
        //   or does it not make sense to have it?
        this->max_iterations_ = 1000;
    }
}

/****************************************
 * Public methods in alphabetical order
 ****************************************/

// NOTE: Used from GasLiftStage2
template<typename TypeTag>
std::optional<typename GasLiftSingleWell<TypeTag>::GradInfo>
GasLiftSingleWell<TypeTag>::
calcIncOrDecGradient(double oil_rate, double gas_rate, double alq, bool increase) const
{
    auto [new_alq_opt, alq_is_limited] = addOrSubtractAlqIncrement_(alq, increase);
    // TODO: What to do if ALQ is limited and new_alq != alq?
    if (!new_alq_opt)
        return std::nullopt;
    double new_alq = *new_alq_opt;
    if (auto bhp = computeBhpAtThpLimit_(new_alq)) {
        auto [new_bhp, bhp_is_limited] = getBhpWithLimit_(*bhp);
        // TODO: What to do if BHP is limited?
        std::vector<double> potentials(this->num_phases_, 0.0);
        computeWellRates_(new_bhp, potentials);
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

/* - At this point we know that this is a production well, and that its current
 * control mode is THP.
 *
 * - We would like to check if it is possible to
 *   1) increase the oil production by adding lift gas injection to the
 *   well, or if that is not possible, if we 2) should reduce the amount
 *   of lift gas injected due to a too small gain in oil production
 *   (with the current lift gas injection rate)
 * - For 1) above, we should not add lift gas if it would cause an oil
 *   rate target to be exceeded, and for 2) we should not reduce the
 *   amount of liftgas injected below the minimum lift gas injection
 *   rate.
 *
 *   NOTE: If reducing or adding lift-gas further would cause
 *     one of the well targets like ORAT, WRAT, GRAT, LRAT, CRAT, RESV, BHP,
 *     to become violated we should stop the lift gas optimization
 *     loop.. and then updateWellControls() will later (hopefully) switch the well's
 *     control mode from THP to the mode of the violated target.
 *
 * - Lift gas is added if it is economical viable depending on
 * the ratio of oil gained compared to the amount of liftgas added.
 *
 * - Lift gas supply may be limited.
 *
 * - The current value of liftgas for the well is stored in the WellState object.
 *
 * - It is assumed that the oil production rate is concave function F
 *   of the amount of lift gas, such that it increases initially due to the
 *   reduced density of the mixture in the tubing. However, as the
 *   lift gas supply is increased further, friction pressure losses in the
 *   tubing become more important, and the production rate peaks and
 *   then starts to decrease.
 *   Since lift gas injection has a price, e.g. compression costs can
 *   be expressed as a cost per unit rate of lift gas injection,
 *   it must be balanced against the value of the extra amount of
 *   oil produced. Thus there is a  "minimum economic gradient" of oil
 *   production rate versus lift gas injection rate, at which the
 *   value of the extra amount of oil produced by a small increase in
 *   the lift gas injection rate is equal to the cost of supplying the
 *   extra amount of lift gas. The optimum lift gas injection rate is then somewhat
 *   lower than the peak value.
 *
 *   Based on this assumption, we know that if the gradient (derivative) of F is
 *   positive, but greater than the economic gradient (assuming the
 *   economic gradient is positive), we should add
 *   lift gas. On the other hand, if the gradient of F is negative or
 *   if it is positive but smaller than the economic gradient, the amount
 *   of lift gas injected should be decreased.
 *
 */
template<typename TypeTag>
std::unique_ptr<typename GasLiftSingleWell<TypeTag>::GLiftWellState>
GasLiftSingleWell<TypeTag>::
runOptimize()
{
    std::unique_ptr<GLiftWellState> state;
    if (this->optimize_) {
        if (this->debug_limit_increase_decrease_) {
            state = runOptimize1_();
        }
        else {
            state = runOptimize2_();
        }
        if (state) {
            if (state->increase()) {
                double alq = state->alq();
                if (this->debug)
                    logSuccess_(alq);
                this->well_state_.setALQ(this->well_name_, alq);
            }
        }
    }
    return state;
}

template<typename TypeTag>
std::unique_ptr<typename GasLiftSingleWell<TypeTag>::GLiftWellState>
GasLiftSingleWell<TypeTag>::
runOptimize2_()
{
    std::unique_ptr<GLiftWellState> state;
    state = tryIncreaseLiftGas_();
    if (!state || !(state->alqChanged())) {
        state = tryDecreaseLiftGas_();
    }
    return state;
}

template<typename TypeTag>
std::unique_ptr<typename GasLiftSingleWell<TypeTag>::GLiftWellState>
GasLiftSingleWell<TypeTag>::
runOptimize1_()
{
    std::unique_ptr<GLiftWellState> state;
    auto inc_count = this->well_state_.gliftGetAlqIncreaseCount(this->well_name_);
    auto dec_count = this->well_state_.gliftGetAlqDecreaseCount(this->well_name_);
    if (dec_count == 0 && inc_count == 0) {
        state = tryIncreaseLiftGas_();
        if (!state && !state->alqChanged()) {
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


/****************************************
 * Private methods in alphabetical order
 ****************************************/

template<typename TypeTag>
std::pair<std::optional<double>, bool>
GasLiftSingleWell<TypeTag>::
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

template<typename TypeTag>
double
GasLiftSingleWell<TypeTag>::
calcEcoGradient_(
    double oil_rate, double new_oil_rate, double gas_rate,
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

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::
checkALQequal_(double alq1, double alq2) const
{
    return std::fabs(alq1-alq2) < (this->increment_*ALQ_EPSILON);
}

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::
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


template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::
checkWellRatesViolated_(
    std::vector<double> &potentials,
    const std::function<bool(double, double, const std::string &)> &callback,
    bool increase
)
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


template<typename TypeTag>
std::optional<double>
GasLiftSingleWell<TypeTag>::
computeBhpAtThpLimit_(double alq) const
{
    auto bhp_at_thp_limit = this->std_well_.computeBhpAtThpLimitProdWithAlq(
        this->ebos_simulator_,
        this->summary_state_,
        this->deferred_logger_,
        alq);
    if (bhp_at_thp_limit) {
        if (*bhp_at_thp_limit < this->controls_.bhp_limit) {
            const std::string msg = fmt::format(
                "Computed bhp ({}) from thp limit is below bhp limit ({}), (ALQ = {})."
                " Using bhp limit instead",
                *bhp_at_thp_limit, this->controls_.bhp_limit, alq);
            displayDebugMessage_(msg);
            bhp_at_thp_limit = this->controls_.bhp_limit;
        }
        //bhp_at_thp_limit = std::max(*bhp_at_thp_limit, this->controls_.bhp_limit);
    }
    else {
        const std::string msg = fmt::format(
            "Failed in getting converged bhp potential from thp limit (ALQ = {})", alq);
        displayDebugMessage_(msg);
    }
    return bhp_at_thp_limit;
}

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::
computeInitialWellRates_(std::vector<double> &potentials)
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

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
computeWellRates_(
    double bhp, std::vector<double> &potentials, bool debug_output) const
{
    // NOTE: If we do not clear the potentials here, it will accumulate
    //   the new potentials to the old values..
    std::fill(potentials.begin(), potentials.end(), 0.0);
    this->std_well_.computeWellRatesWithBhp(
        this->ebos_simulator_, bhp, potentials, this->deferred_logger_);
    if (debug_output) {
        const std::string msg = fmt::format("computed well potentials given bhp {}, "
            "oil: {}, gas: {}, water: {}", bhp,
            -potentials[this->oil_pos_], -potentials[this->gas_pos_],
            -potentials[this->water_pos_]);
        displayDebugMessage_(msg);
    }
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
debugCheckNegativeGradient_(double grad, double alq, double new_alq, double oil_rate,
    double new_oil_rate, double gas_rate, double new_gas_rate, bool increase) const
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

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
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

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
debugShowStartIteration_(double alq, bool increase, double oil_rate)
{
    const std::string msg =
        fmt::format("starting {} iteration, ALQ = {}, oilrate = {}",
            (increase ? "increase" : "decrease"),
            alq, oil_rate);
    displayDebugMessage_(msg);
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
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

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
debugShowAlqIncreaseDecreaseCounts_()
{
    auto inc_count = this->well_state_.gliftGetAlqIncreaseCount(this->well_name_);
    auto dec_count = this->well_state_.gliftGetAlqDecreaseCount(this->well_name_);
    const std::string msg =
        fmt::format("ALQ increase/decrease count : {}/{}", inc_count, dec_count);
    displayDebugMessage_(msg);

}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
displayDebugMessage_(const std::string &msg) const
{

    if (this->debug) {
        const std::string message = fmt::format(
            "  GLIFT (DEBUG) : Well {} : {}", this->well_name_, msg);
        this->deferred_logger_.info(message);
    }
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
displayWarning_(std::string msg)
{
    const std::string message = fmt::format(
        "GAS LIFT OPTIMIZATION, WELL {} : {}", this->well_name_, msg);
    this->deferred_logger_.warning("WARNING", message);
}

template<typename TypeTag>
std::pair<double, bool>
GasLiftSingleWell<TypeTag>::
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
template<typename TypeTag>
std::pair<double, bool>
GasLiftSingleWell<TypeTag>::
getGasRateWithLimit_(const std::vector<double> &potentials) const
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


// NOTE: If the computed oil rate is larger than the target
//   rate of the well, we reduce it to the target rate. This
//   will make the economic gradient smaller than it would be
//   if we did not reduce the rate, and it is less
//   likely that the current gas lift increment will be
//   accepted.
// TODO: If it still is accepted, we should ideally reduce the alq
//  also since we also reduced the rate. This might involve
//   some sort of iteration though..
template<typename TypeTag>
std::pair<double, bool>
GasLiftSingleWell<TypeTag>::
getOilRateWithLimit_(const std::vector<double> &potentials) const
{
    double new_rate = -potentials[this->oil_pos_];
    if (this->controls_.hasControl(Well::ProducerCMode::ORAT)) {
        auto target = this->controls_.oil_rate;
        if (new_rate > target) {
            const std::string msg = fmt::format("limiting oil rate to target: "
                "computed rate: {}, target: {}", new_rate, target);
            displayDebugMessage_(msg);
            new_rate = target;
            return { new_rate, /*limit=*/true};
        }
    }
    // TODO: What to do if both oil and lrat are violated?
    //      Currently oil rate violation takes precedence if both are
    //      violated, and the new rate is set to the oil rate target.
    if (this->controls_.hasControl(Well::ProducerCMode::LRAT)) {
        auto target = this->controls_.liquid_rate;
        auto oil_rate = -potentials[this->oil_pos_];
        auto water_rate = -potentials[this->water_pos_];
        auto liq_rate = oil_rate + water_rate;
        if (liq_rate > target) {
            new_rate = oil_rate * (target/liq_rate);
            const std::string msg = fmt::format(
                "limiting oil rate due to LRAT target {}: "
                "computed rate: {}, target: {}", target, oil_rate, new_rate);
            displayDebugMessage_(msg);
            return { new_rate, /*limit=*/true};
       }
    }
    return { new_rate, /*limit=*/false};
}

template<typename TypeTag>
std::tuple<double,double,bool,bool>
GasLiftSingleWell<TypeTag>::
getInitialRatesWithLimit_(const std::vector<double> &potentials)
{
    auto [oil_rate, oil_is_limited] = getOilRateWithLimit_(potentials);
    auto [gas_rate, gas_is_limited] = getGasRateWithLimit_(potentials);

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
    return std::make_tuple(oil_rate, gas_rate, oil_is_limited, gas_is_limited);
}

template<typename TypeTag>
std::tuple<double,double,bool,bool,double>
GasLiftSingleWell<TypeTag>::
increaseALQtoPositiveOilRate_(double alq, double oil_rate, double gas_rate,
    bool oil_is_limited, bool gas_is_limited, std::vector<double> &potentials)
{
    bool stop_iteration = false;
    double temp_alq = alq;
    while(!stop_iteration) {
        temp_alq += this->increment_;
        if (temp_alq > this->max_alq_) break;
        auto bhp_opt = computeBhpAtThpLimit_(temp_alq);
        if (!bhp_opt) break;
        alq = temp_alq;
        auto [bhp, bhp_is_limited] = getBhpWithLimit_(*bhp_opt);
        computeWellRates_(bhp, potentials);
        oil_rate = -potentials[this->oil_pos_];
        if (oil_rate > 0) break;
    }
    std::tie(oil_rate, oil_is_limited) = getOilRateWithLimit_(potentials);
    std::tie(gas_rate, gas_is_limited) = getGasRateWithLimit_(potentials);
    return std::make_tuple(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq);
}

template<typename TypeTag>
std::tuple<double,double,bool,bool,double>
GasLiftSingleWell<TypeTag>::
increaseALQtoMinALQ_(double alq, double oil_rate, double gas_rate,
    bool oil_is_limited, bool gas_is_limited, std::vector<double> &potentials)
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
        auto [bhp, bhp_is_limited] = getBhpWithLimit_(*bhp_opt);
        computeWellRates_(bhp, potentials);
        std::tie(oil_rate, oil_is_limited) = getOilRateWithLimit_(potentials);
        std::tie(gas_rate, gas_is_limited) = getGasRateWithLimit_(potentials);
        if (oil_is_limited || gas_is_limited) break;
    }
    return std::make_tuple(oil_rate, gas_rate, oil_is_limited, gas_is_limited, alq);
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
logSuccess_(double alq)
{
    const int iteration_idx =
        this->ebos_simulator_.model().newtonMethod().numIterations();
    const std::string message = fmt::format(
         "GLIFT, IT={}, WELL {} : {} ALQ from {} to {}",
         iteration_idx,
         this->well_name_,
         ((alq > this->orig_alq_) ? "increased" : "decreased"),
         this->orig_alq_, alq);
    this->deferred_logger_.info(message);
}

template<typename TypeTag>
std::tuple<double,double,bool,bool,double>
GasLiftSingleWell<TypeTag>::
reduceALQtoOilTarget_(double alq, double oil_rate, double gas_rate,
    bool oil_is_limited, bool gas_is_limited, std::vector<double> &potentials)
{
    displayDebugMessage_("Reducing ALQ to meet oil target before iteration starts..");
    double orig_oil_rate = oil_rate;
    double orig_alq = alq;
    // NOTE: This method should only be called if oil_is_limited, and hence
    //   we know that it has oil rate control
    assert(this->controls_.hasControl(Well::ProducerCMode::ORAT));
    auto target = this->controls_.oil_rate;
    bool stop_iteration = false;
    double temp_alq = alq;
    while(!stop_iteration) {
        temp_alq -= this->increment_;
        if (temp_alq <= 0) break;
        auto bhp_opt = computeBhpAtThpLimit_(temp_alq);
        if (!bhp_opt) break;
        alq = temp_alq;
        auto [bhp, bhp_is_limited] = getBhpWithLimit_(*bhp_opt);
        computeWellRates_(bhp, potentials);
        oil_rate = -potentials[this->oil_pos_];
        if (oil_rate < target) {
            break;
        }
    }
    std::tie(oil_rate, oil_is_limited) = getOilRateWithLimit_(potentials);
    std::tie(gas_rate, gas_is_limited) = getGasRateWithLimit_(potentials);
    if (this->debug) {
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
//  - return value: a new GLiftWellState or nullptr
//
template<typename TypeTag>
std::unique_ptr<typename GasLiftSingleWell<TypeTag>::GLiftWellState>
GasLiftSingleWell<TypeTag>::
runOptimizeLoop_(bool increase)
{
    std::vector<double> potentials(this->num_phases_, 0.0);
    std::unique_ptr<GLiftWellState> ret_value; // nullptr initially
    if (!computeInitialWellRates_(potentials)) return ret_value;
    bool alq_is_limited = false;
    bool oil_is_limited = false;
    bool gas_is_limited = false;
    double oil_rate, gas_rate;
    std::tie(oil_rate, gas_rate, oil_is_limited, gas_is_limited) =
        getInitialRatesWithLimit_(potentials);
    if (this->debug) debugShowBhpAlqTable_();
    if (this->debug) debugShowAlqIncreaseDecreaseCounts_();
    if (this->debug) debugShowTargets_();
    bool success = false;  // did we succeed to increase alq?
    auto cur_alq = this->orig_alq_;
    auto temp_alq = cur_alq;
    if (!increase && oil_is_limited) {
        // NOTE: Try to decrease ALQ down to a value where the oil target is
        //   not exceeded.
        // NOTE: This may reduce ALQ below the minimum value set in WLIFTOPT
        //   item 5. However, this is OK since the rate target is met and there
        //   is no point in using a higher ALQ value then.
        std::tie(oil_rate, gas_rate, oil_is_limited, gas_is_limited, temp_alq) =
            reduceALQtoOilTarget_(
                temp_alq, oil_rate, gas_rate, oil_is_limited, gas_is_limited, potentials);
    }
    else {
        if (increase && oil_rate < 0) {
            // NOTE: Try to increase ALQ up to a value where oil_rate is positive
            std::tie(oil_rate, gas_rate, oil_is_limited, gas_is_limited, temp_alq) =
               increaseALQtoPositiveOilRate_(temp_alq, oil_rate,
                   gas_rate, oil_is_limited, gas_is_limited, potentials);
        }
        if (increase && (this->min_alq_> 0) && (temp_alq < this->min_alq_)) {
            // NOTE: Try to increase ALQ up to the minimum limit without checking
            //   the economic gradient..
            std::tie(oil_rate, gas_rate, oil_is_limited, gas_is_limited, temp_alq) =
                increaseALQtoMinALQ_(temp_alq, oil_rate, gas_rate,
                    oil_is_limited, gas_is_limited, potentials);
        }
    }
    if (checkInitialALQmodified_(temp_alq, cur_alq)) {
        cur_alq = temp_alq;
        success = true;
    }
    OptimizeState state {*this, increase};
    if (this->debug) debugShowStartIteration_(temp_alq, increase, oil_rate);
    while (!state.stop_iteration && (++state.it <= this->max_iterations_)) {
        if (!increase && state.checkNegativeOilRate(oil_rate)) break;
        if (state.checkWellRatesViolated(potentials)) break;
        if (state.checkAlqOutsideLimits(temp_alq, oil_rate)) break;
        std::optional<double> alq_opt;
        std::tie(alq_opt, alq_is_limited)
            = state.addOrSubtractAlqIncrement(temp_alq);
        if (!alq_opt) break;
        temp_alq = *alq_opt;
        if (this->debug) state.debugShowIterationInfo(temp_alq);
        if (!state.computeBhpAtThpLimit(temp_alq)) break;
        // NOTE: if BHP is below limit, we set state.stop_iteration = true
        auto bhp = state.getBhpWithLimit();
        computeWellRates_(bhp, potentials);
        auto [new_oil_rate, new_oil_is_limited] = getOilRateWithLimit_(potentials);
/*        if (this->debug_abort_if_decrease_and_oil_is_limited_) {
            if (oil_is_limited && !increase) {
                // if oil is limited we do not want to decrease
                displayDebugMessage_(
                    "decreasing ALQ and oil is limited -> aborting iteration");
                break;
            }
            }
*/
        auto [new_gas_rate, new_gas_is_limited] = getGasRateWithLimit_(potentials);
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
        if (this->debug)
            debugCheckNegativeGradient_(
                gradient, cur_alq, temp_alq, oil_rate, new_oil_rate,
                gas_rate, new_gas_rate, increase);
        if (state.checkEcoGradient(gradient)) break;
        cur_alq = temp_alq;
        success = true;
        oil_rate = new_oil_rate;
        gas_rate = new_gas_rate;
        oil_is_limited = new_oil_is_limited;
        gas_is_limited = new_gas_is_limited;
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
    ret_value = std::make_unique<GLiftWellState>(oil_rate, oil_is_limited,
        gas_rate, gas_is_limited, cur_alq, alq_is_limited, increase_opt);
    return ret_value;
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
setAlqMaxRate_(const GasLiftOpt::Well &well)
{
    auto& max_alq_optional = well.max_rate();
    if (max_alq_optional) {
        // NOTE: To prevent extrapolation of the VFP tables, any value
        // entered here must not exceed the largest ALQ value in the well's VFP table.
        this->max_alq_ = *max_alq_optional;
    }
    else { // i.e. WLIFTOPT, item 3 has been defaulted
        // According to the manual for WLIFTOPT, item 3:
        //   The default value should be set to the largest ALQ
        //   value in the well's VFP table
        const auto& table = std_well_.vfp_properties_->getProd()->getTable(
                this->controls_.vfp_table_number);
        const auto& alq_values = table.getALQAxis();
        // Assume the alq_values are sorted in ascending order, so
        // the last item should be the largest value:
        this->max_alq_ = alq_values.back();
    }
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
setAlqMinRate_(const GasLiftOpt::Well &well)
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

template<typename TypeTag>
std::unique_ptr<typename GasLiftSingleWell<TypeTag>::GLiftWellState>
GasLiftSingleWell<TypeTag>::
tryDecreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/ false);
}

template<typename TypeTag>
std::unique_ptr<typename GasLiftSingleWell<TypeTag>::GLiftWellState>
GasLiftSingleWell<TypeTag>::
tryIncreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/ true);
}


// Called when we should use a fixed ALQ value
template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
updateWellStateAlqFixedValue_(const GasLiftOpt::Well &well)
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
template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::
useFixedAlq_(const GasLiftOpt::Well &well)
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

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
warnMaxIterationsExceeded_()
{
    const std::string msg = fmt::format(
        "Max iterations ({}) exceeded", this->max_iterations_);
    displayWarning_(msg);
}

/****************************************
 * Methods declared in OptimizeState
 ****************************************/

template<typename TypeTag>
std::pair<std::optional<double>, bool>
GasLiftSingleWell<TypeTag>::OptimizeState::
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

template<typename TypeTag>
double
GasLiftSingleWell<TypeTag>::OptimizeState::
calcEcoGradient(
       double oil_rate, double new_oil_rate, double gas_rate, double new_gas_rate)
{
    return this->parent.calcEcoGradient_(
        oil_rate, new_oil_rate, gas_rate, new_gas_rate, this->increase);
}

// NOTE:  According to WLIFTOPT item 5 :
//   if min_rate() is negative, it means: allocate at least enough lift gas
//   to enable the well to flow
//  We will interpret this as (see discussion above GasLiftSingleWell()
//   in this file): Allocate at least the amount of lift gas needed to
//   get a positive oil production rate.
template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::OptimizeState::
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
    if (this->parent.debug) {
        const std::string msg = ss.str();
        if (!msg.empty())
            this->parent.displayDebugMessage_(msg);
    }
    return result;
}

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::OptimizeState::
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
template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::OptimizeState::
checkEcoGradient(double gradient)
{
    std::ostringstream ss;
    bool result = false;

    if (this->parent.debug) {
        ss << "checking gradient: " << gradient;
    }
    if (this->increase) {
        if (this->parent.debug) ss << " <= " << this->parent.eco_grad_ << " --> ";
        if (gradient <= this->parent.eco_grad_) {
            if (this->parent.debug) ss << "yes, stopping";
            result = true;
        }
        else {
            if (this->parent.debug) ss << "no, continue";
        }
    }
    else {  // decreasing lift gas
        if (this->parent.debug) ss << " >= " << this->parent.eco_grad_ << " --> ";
        if (gradient >= this->parent.eco_grad_) {
            if (this->parent.debug) ss << "yes, stopping";
            result = true;
        }
        else {
            if (this->parent.debug) ss << "no, continue";
        }
    }
    if (this->parent.debug) this->parent.displayDebugMessage_(ss.str());
    return result;
}

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::OptimizeState::
checkRate(double rate, double limit, const std::string &rate_str) const
{
    if (limit < rate) {
        if (this->parent.debug) {
            const std::string msg = fmt::format(
                "iteration {} : {} rate {} exceeds target {}. Stopping iteration",
                this->it, rate_str, rate, limit);
            this->parent.displayDebugMessage_(msg);
        }
        return true;
    }
    return false;
}

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::OptimizeState::
checkWellRatesViolated(std::vector<double> &potentials)
{
    auto callback = [*this](double rate, double limit, const std::string &rate_str)
                    -> bool
                    { return this->checkRate(rate, limit, rate_str); };
    return this->parent.checkWellRatesViolated_(
              potentials, callback, this->increase);
}

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::OptimizeState::
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

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::OptimizeState::
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
template<typename TypeTag>
double
GasLiftSingleWell<TypeTag>::OptimizeState::
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
