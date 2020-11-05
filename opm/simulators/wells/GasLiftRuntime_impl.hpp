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
Opm::GasLiftRuntime<TypeTag>::
GasLiftRuntime(
    const StdWell &std_well,
    const Simulator &ebos_simulator,
    const SummaryState &summary_state,
    DeferredLogger &deferred_logger,
    const WellState &well_state,
    const Well::ProductionControls &controls
) :
    controls_{controls},
    deferred_logger_{deferred_logger},
    ebos_simulator_{ebos_simulator},
    potentials_(well_state.numPhases(),0.0),
    std_well_{std_well},
    summary_state_{summary_state},
    well_state_{well_state},
    debug{false}  // extra debugging output
{
    int well_index = this->std_well_.indexOfWell();
    const Well::ProducerCMode& control_mode
        = well_state_.currentProductionControls()[well_index];
    if (control_mode != Well::ProducerCMode::THP)
        throw std::logic_error("Bug in flow - invalid control mode detected\n");
    const Opm::Schedule& schedule = this->ebos_simulator_.vanguard().schedule();
    const int report_step_idx = this->ebos_simulator_.episodeIndex();
    auto ecl_well = this->std_well_.wellEcl();
    this->well_name_ = ecl_well.name();
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
    computeInitialWellRates_();

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

        const auto& pu = std_well_.phaseUsage();
        this->oil_pos_ = pu.phase_pos[Oil];
        this->gas_pos_ = pu.phase_pos[Gas];
        this->water_pos_ = pu.phase_pos[Water];
        this->new_alq_ = this->orig_alq_;
        // TODO: adhoc value.. Should we keep max_iterations_ as a safety measure
        //   or does it not make sense to have it?
        this->max_iterations_ = 1000;
    }
}

/****************************************
 * Methods in alphabetical order
 ****************************************/

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
computeInitialWellRates_()
{
    // get the alq value used for this well for the previous time step, or
    //   if gas lift optimization has not been applied to this well yet, the
    //   default value is used.
    this->orig_alq_ = this->well_state_.getALQ(this->well_name_);
    // NOTE: compute initial rates with current ALQ
    this->std_well_.computeWellRatesWithThpAlqProd(
        this->ebos_simulator_, this->summary_state_, this->deferred_logger_,
        this->potentials_, this->orig_alq_);

}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
computeWellRates_(double bhp, std::vector<double> &potentials)
{
    this->std_well_.computeWellRatesWithBhp(
        this->ebos_simulator_, bhp, potentials, this->deferred_logger_);

}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
debugShowIterationInfo_(OptimizeState &state, double alq)
{
    const std::string msg = fmt::format("iteration {}, ALQ = {}", state.it, alq);
    this->displayDebugMessage_(msg);
}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
debugShowStartIteration_(double alq, bool increase)
{
    const std::string msg =
        fmt::format("starting {} iteration, ALQ = {}, oilrate = {}",
            (increase ? "increase" : "decrease"),
            alq,
            -this->potentials_[this->oil_pos_]);
    this->displayDebugMessage_(msg);
}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
displayDebugMessage_(const std::string &msg)
{

    const std::string message = fmt::format(
        "  GLIFT (DEBUG) : Well {} : {}", this->well_name_, msg);
    this->deferred_logger_.info(message);
}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
displayWarning_(std::string msg)
{
    const std::string message = fmt::format(
        "GAS LIFT OPTIMIZATION, WELL {} : {}", this->well_name_, msg);
    this->deferred_logger_.warning("WARNING", message);
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
double
Opm::GasLiftRuntime<TypeTag>::
getGasRateWithLimit_(std::vector<double> &potentials)
{
    auto new_rate = -potentials[this->gas_pos_];
    if (this->controls_.hasControl(Well::ProducerCMode::GRAT)) {
        auto target = this->controls_.gas_rate;
        if (new_rate > target)
            new_rate = target;
    }
    return new_rate;
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
double
Opm::GasLiftRuntime<TypeTag>::
getOilRateWithLimit_(std::vector<double> &potentials)
{
    auto new_rate = -potentials[this->oil_pos_];
    if (this->controls_.hasControl(Well::ProducerCMode::ORAT)) {
        auto target = this->controls_.oil_rate;
        if (new_rate > target)
            new_rate = target;
    }
    if (this->controls_.hasControl(Well::ProducerCMode::LRAT)) {
        auto target = this->controls_.liquid_rate;
        auto oil_rate = -potentials[this->oil_pos_];
        auto water_rate = -potentials[this->water_pos_];
        auto liq_rate = oil_rate + water_rate;
        if (liq_rate > target) {
            new_rate = oil_rate * (target/liq_rate);
       }
    }
    return new_rate;
}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
logSuccess_()
{

    const std::string message = fmt::format(
         "GLIFT, WELL {} {} ALQ from {} to {}",
         this->well_name_,
         ((this->new_alq_ > this->orig_alq_) ? "increased" : "decreased"),
         this->orig_alq_,
         this->new_alq_);
    this->deferred_logger_.info(message);
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
void
Opm::GasLiftRuntime<TypeTag>::
runOptimize()
{
    if (this->optimize_) {
        if (!tryIncreaseLiftGas_()) {
            if (!tryDecreaseLiftGas_()) {
                return;
            }
        }
        logSuccess_();
        this->well_state_.setALQ(this->well_name_, this->new_alq_);
    }
    // NOTE: In addition to the new ALQ value, we also implicitly
    //    return this->potentials_
}

// INPUT:
// - increase (boolean) :
//   - true  : try increase the lift gas supply,
//   - false : try decrease lift gas supply.
//
// OUTPUT:
// - return value: success (true/false)
// - potentials_ : updated well potentials if success
// - new_alq_    : updated alq if success
template<typename TypeTag>
bool
Opm::GasLiftRuntime<TypeTag>::
runOptimizeLoop_(bool increase)
{
    auto cur_potentials = this->potentials_;  // make copy, since we may fail..
    auto oil_rate = -cur_potentials[this->oil_pos_];
    auto gas_rate = -cur_potentials[this->gas_pos_];
    bool success = false;  // did we succeed to increase alq?
    auto cur_alq = this->orig_alq_;
    if (cur_alq == 0 && !increase) // we don't decrease alq below zero
        return false;

    auto alq = cur_alq;
    OptimizeState state {*this, increase};
    if (this->debug) debugShowStartIteration_(alq, increase);
    while (!state.stop_iteration && (++state.it <= this->max_iterations_)) {
        if (state.checkWellRatesViolated(cur_potentials)) break;
        if (state.checkAlqOutsideLimits(alq, oil_rate)) break;
        alq = state.addOrSubtractAlqIncrement(alq);
        if(this->debug) debugShowIterationInfo_(state, alq);
        if (!state.computeBhpAtThpLimit(alq)) break;
        // NOTE: if BHP is below limit, we set state.stop_iteration = true
        auto bhp = state.getBhpWithLimit();
        computeWellRates_(bhp, cur_potentials);
        auto new_oil_rate = getOilRateWithLimit_(cur_potentials);
        auto new_gas_rate = getGasRateWithLimit_(cur_potentials);
        auto gradient = state.calcGradient(
            oil_rate, new_oil_rate, gas_rate, new_gas_rate);
        if (state.checkEcoGradient(gradient)) {
            if (state.it == 1) {
                break;
            }
            else {
                state.stop_iteration = true;
            }
        }
        cur_alq = alq;
        success = true;
        oil_rate = new_oil_rate;
        gas_rate = new_gas_rate;
    }
    if (state.it > this->max_iterations_) {
        warnMaxIterationsExceeded_();
    }
    if (success) {
        this->potentials_ = cur_potentials;
        this->new_alq_ = cur_alq;
    }
    return success;
}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
setAlqMaxRate_(const GasLiftOpt::Well &well)
{
    auto& max_alq_optional = well.max_rate();
    if (max_alq_optional) {
        // NOTE: To prevent extrapolation of the VFP tables, any value
        // entered here must not exceed the largest ALQ value in the well's VFP table.
        this->max_alq_ = *max_alq_optional;
    }
    else { // i.e. WLIFTOPT, item 3 has been defaulted
        // According to the Eclipse manual for WLIFTOPT, item 3:
        //   The default value should be set to the largest ALQ
        //   value in the well's VFP table
        const auto& table = *(std_well_.vfp_properties_->getProd()->getTable(
                this->controls_.vfp_table_number));
        const auto& alq_values = table.getALQAxis();
        // Assume the alq_values are sorted in ascending order, so
        // the last item should be the largest value:
        this->max_alq_ = alq_values.back();
    }
}

template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
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
bool
Opm::GasLiftRuntime<TypeTag>::
tryDecreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/ false);
}

template<typename TypeTag>
bool
Opm::GasLiftRuntime<TypeTag>::
tryIncreaseLiftGas_()
{
    return runOptimizeLoop_(/*increase=*/ true);
}


// Called when we should use a fixed ALQ value
template<typename TypeTag>
void
Opm::GasLiftRuntime<TypeTag>::
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
Opm::GasLiftRuntime<TypeTag>::
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
Opm::GasLiftRuntime<TypeTag>::
warnMaxIterationsExceeded_()
{
    std::ostringstream ss;
    ss << "Max iterations (" << this->max_iterations_ << ") exceeded in "
       << "gas lift optimization for well " << this->well_name_;
    deferred_logger_.warning("MAX_ITERATIONS_EXCEEDED", ss.str());
}

/****************************************
 * Methods declared in OptimizeState
 ****************************************/

template<typename TypeTag>
double
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
addOrSubtractAlqIncrement(double alq)
{
    if (this->increase) {
        alq += this->parent.increment_;
        // NOTE: if max_alq_ was defaulted in WCONPROD, item 3, it has
        //   already been set to the largest value in the VFP table in
        //   the contructor of GasLiftRuntime
        if (alq > this->parent.max_alq_) alq = this->parent.max_alq_;
    }
    else {
        alq -= this->parent.increment_;
        if (this->parent.min_alq_ > 0) {
            if (alq < this->parent.min_alq_) alq = this->parent.min_alq_;
        }
    }
    return alq;
}

template<typename TypeTag>
double
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
calcGradient(double oil_rate, double new_oil_rate, double gas_rate, double new_gas_rate)
{
    auto dqo = new_oil_rate - oil_rate;
    auto dqg = new_gas_rate - gas_rate;
    // TODO: Should we do any error checks on the calculation of the
    //   gradient?
    // NOTE: The eclipse techincal description (chapter 25) says:
    //   "The gas rate term in the denominator is subject to the
    //   constraint alpha_g_ * dqg >= 0.0"
    auto gradient = (this->parent.alpha_w_ * dqo) /
        (this->parent.increment_ + this->parent.alpha_g_*dqg);
    return gradient;
}

// NOTE:  According to WLIFTOPT item 5 :
//   if min_rate() is negative, it means: allocate at least enough lift gas
//   to enable the well to flow
//  We will interpret this as (see discussion above GasLiftRuntime()
//   in this file): Allocate at least the amount of lift gas needed to
//   get a positive oil production rate.
template<typename TypeTag>
bool
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
checkAlqOutsideLimits(double alq, double oil_rate)
{
    std::ostringstream ss;
    bool result = false;

    if (this->increase) {
        if (alq >= this->parent.max_alq_) {
            ss << "ALQ >= " << this->parent.max_alq_ << " (max limit), "
               << "stopping iteration";
            result = true;
        }
        else {
            // NOTE: A negative min_alq_ means: allocate at least enough lift gas
            //  to enable the well to flow, see WCONPROD item 5.
            if (this->parent.min_alq_ < 0) {
                result = false;
            }
            else {
                // NOTE: checking for a lower limit should not be necessary
                // when increasing alq.. so this is just to catch an
                // illegal state at an early point.
                if (alq < this->parent.min_alq_ ) {
                    warn_("unexpected: alq below lower limit when trying to "
                        "increase lift gas. aborting iteration.");
                    result = true;
                }
                else {
                    result = false;
                }
            }
        }
    }
    else { // we are decreasing lift gas
        // NOTE: A negative min_alq_ means: allocate at least enough lift gas
        //  to enable the well to flow, see WCONPROD item 5.
        if (this->parent.min_alq_ < 0) {
            // If the oil rate is already zero or negative (non-flowing well)
            // we assume we will not be able to increase it by decreasing the lift gas
            if ( oil_rate <= 0 ) {
                ss << "Oil rate ( " << oil_rate << " ) <= 0 when decreasing lift gas. "
                   << "We will not be able to make this well flowing by decreasing "
                   << "lift gas, stopping iteration.";
                result = true;
            }
            else {
                result = false;
            }
        }
        else {
            if (alq <= this->parent.min_alq_ ) {
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
    if (this->parent.debug) this->parent.displayDebugMessage_(ss.str());
    return result;
}

template<typename TypeTag>
bool
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
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
            if (this->parent.debug) ss << "true";
            result = true;
        }
        else {
            if (this->parent.debug) ss << "false";
        }
    }
    else {  // decreasing lift gas
        if (this->parent.debug) ss << " >= " << this->parent.eco_grad_ << " --> ";
        if (gradient >= this->parent.eco_grad_) {
            if (this->parent.debug) ss << "true";
            result = true;
        }
        else {
            if (this->parent.debug) ss << "false";
        }
    }
    if (this->parent.debug) this->parent.displayDebugMessage_(ss.str());
    return result;
}

template<typename TypeTag>
bool
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
checkRate(double rate, double limit, const std::string rate_str)
{

    if (limit < rate ) {
        if (this->parent.debug) {
            const std::string msg = fmt::format(
                "iteration {} : rate {} exceeds target rate {}. Stopping iteration",
                this->it, rate_str, rate, limit);
            this->parent.displayDebugMessage_(msg);
        }
        return true;
    }
    return false;
}

template<typename TypeTag>
bool
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
checkWellRatesViolated(std::vector<double> &potentials)
{
    if (this->parent.controls_.hasControl(Well::ProducerCMode::ORAT)) {
        auto oil_rate = -potentials[this->parent.oil_pos_];
        if (this->checkRate(oil_rate, this->parent.controls_.oil_rate, "oil"))
            return true;
    }
    if (this->parent.controls_.hasControl(Well::ProducerCMode::WRAT)) {
        auto water_rate = -potentials[this->parent.water_pos_];
        if (this->checkRate(water_rate, this->parent.controls_.water_rate, "water"))
            return true;
    }
    if (this->parent.controls_.hasControl(Well::ProducerCMode::GRAT)) {
        auto gas_rate = -potentials[this->parent.gas_pos_];
        if (this->checkRate(gas_rate, this->parent.controls_.gas_rate, "gas"))
            return true;
    }
    if (this->parent.controls_.hasControl(Well::ProducerCMode::LRAT)) {
        auto oil_rate = -potentials[this->parent.oil_pos_];
        auto water_rate = -potentials[this->parent.water_pos_];
        auto liq_rate = oil_rate + water_rate;
        if (this->checkRate(liq_rate, this->parent.controls_.liquid_rate, "liquid"))
            return true;
    }
    // TODO: Also check RESV, see checkIndividualContraints() in
    //   WellInterface_impl.hpp
    // TODO: Check group contraints?

    return false;
}



template<typename TypeTag>
bool
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
computeBhpAtThpLimit(double alq)
{
    auto bhp_at_thp_limit = this->parent.std_well_.computeBhpAtThpLimitProdWithAlq(
        this->parent.ebos_simulator_,
        this->parent.summary_state_,
        this->parent.deferred_logger_,
        alq);
    if (!bhp_at_thp_limit) {
        const std::string msg = fmt::format(
          "Failed in getting converged bhp potential for well {}",
          this->parent.well_name_);
        this->parent.deferred_logger_.warning(
            "FAILURE_GETTING_CONVERGED_POTENTIAL", msg);
        return false;
    }
    this->bhp = *bhp_at_thp_limit;
    return true;
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
Opm::GasLiftRuntime<TypeTag>::OptimizeState::
getBhpWithLimit()
{
    auto bhp_update = this->bhp;
    if (this->parent.controls_.hasControl(Well::ProducerCMode::BHP)) {
        auto limit = this->parent.controls_.bhp_limit;
        // TODO: is it possible that bhp falls below the limit when
        // adding lift gas? I.e. if this->increase == true..
        if (this->bhp < limit) {
            // TODO: we keep the current alq, but it should probably
            // be adjusted since we changed computed bhp. But how?
            bhp_update = limit;
            // Stop iteration, but first check the economic gradient
            //   with the bhp_update. If the gradient looks OK (see the
            //   main optimize loop) we keep the current ALQ value.
            this->stop_iteration = true;
        }
    }
    return bhp_update;
}
