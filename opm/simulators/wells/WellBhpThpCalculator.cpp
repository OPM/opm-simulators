/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <cassert>

static constexpr bool extraBhpAtThpLimitOutput = false;

namespace Opm
{

bool
WellBhpThpCalculator::wellHasTHPConstraints(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto controls = well_ecl.injectionControls(summaryState);
        if (controls.hasControl(Well::InjectorCMode::THP))
            return true;
    }

    if (well_ecl.isProducer()) {
        const auto controls = well_ecl.productionControls(summaryState);
        if (controls.hasControl(Well::ProducerCMode::THP))
            return true;
    }

    return false;
}

double WellBhpThpCalculator::getTHPConstraint(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto& controls = well_ecl.injectionControls(summaryState);
        return controls.thp_limit;
    }

    if (well_ecl.isProducer( )) {
        const auto& controls = well_ecl.productionControls(summaryState);
        return controls.thp_limit;
    }

    return 0.0;
}

double WellBhpThpCalculator::mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto& controls = well_ecl.injectionControls(summaryState);
        return controls.bhp_limit;
    }

    if (well_ecl.isProducer( )) {
        const auto& controls = well_ecl.productionControls(summaryState);
        return controls.bhp_limit;
    }

    return 0.0;
}

double WellBhpThpCalculator::calculateThpFromBhp(const std::vector<double>& rates,
                                                 const double bhp,
                                                 const double rho,
                                                 const double alq,
                                                 DeferredLogger& deferred_logger) const
{
    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const double aqua = rates[Water];
    const double liquid = rates[Oil];
    const double vapour = rates[Gas];

    // pick the density in the top layer
    double thp = 0.0;
    if (well_.isInjector()) {
        const int table_id = well_.wellEcl().vfp_table_number();
        const double vfp_ref_depth = well_.vfpProperties()->getInj()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());

        thp = well_.vfpProperties()->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
    }
    else if (well_.isProducer()) {
        const int table_id = well_.wellEcl().vfp_table_number();
        const double vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());

        thp = well_.vfpProperties()->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
    }
    else {
        OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
    }

    return thp;
}

std::optional<double>
WellBhpThpCalculator::
computeBhpAtThpLimitProd(const std::function<std::vector<double>(const double)>& frates,
                         const SummaryState& summary_state,
                         const double maxPerfPress,
                         const double rho,
                         const double alq_value,
                         DeferredLogger& deferred_logger) const
{
    // Given a VFP function returning bhp as a function of phase
    // rates and thp:
    //     fbhp(rates, thp),
    // a function extracting the particular flow rate used for VFP
    // lookups:
    //     flo(rates)
    // and the inflow function (assuming the reservoir is fixed):
    //     frates(bhp)
    // we want to solve the equation:
    //     fbhp(frates(bhp, thplimit)) - bhp = 0
    // for bhp.
    //
    // This may result in 0, 1 or 2 solutions. If two solutions,
    // the one corresponding to the lowest bhp (and therefore
    // highest rate) should be returned.

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    // Make the fbhp() function.
    const auto& controls = well_.wellEcl().productionControls(summary_state);
    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);
    const double vfp_ref_depth = table.getDatumDepth();
    const double thp_limit = this->getTHPConstraint(summary_state);
    const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());

    auto fbhp = [this, &controls, thp_limit, dp, alq_value](const std::vector<double>& rates) {
        assert(rates.size() == 3);
        const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
        const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());
        const bool use_vfpexp = well_.useVfpExplicit();
        return well_.vfpProperties()->getProd()
        ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], thp_limit, alq_value, wfr, gfr, use_vfpexp) - dp;
    };

    // Make the flo() function.
    auto flo = [&table](const std::vector<double>& rates) {
        return detail::getFlo(table, rates[Water], rates[Oil], rates[Gas]);
    };

    // Find the bhp-point where production becomes nonzero.
    auto fflo = [&flo, &frates](double bhp) { return flo(frates(bhp)); };
    auto bhp_max = this->bhpMax(fflo, controls.bhp_limit, maxPerfPress, table.getFloAxis().front(), deferred_logger);

    // could not solve for the bhp-point, we could not continue to find the bhp
    if (!bhp_max.has_value()) {
        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
                                "Robust bhp(thp) solve failed due to not being able to "
                                "find bhp-point where production becomes non-zero for well " + well_.name());
        return std::nullopt;
    }
    const std::array<double, 2> range {controls.bhp_limit, *bhp_max};
    return this->computeBhpAtThpLimit(frates, fbhp, range, deferred_logger);
}

std::optional<double>
WellBhpThpCalculator::
bhpMax(const std::function<double(const double)>& fflo,
       const double bhp_limit,
       const double maxPerfPress,
       const double vfp_flo_front,
       DeferredLogger& deferred_logger) const
{
    // Find the bhp-point where production becomes nonzero.
    double bhp_max = 0.0;
    double low = bhp_limit;
    double high = maxPerfPress + 1.0 * unit::barsa;
    double f_low = fflo(low);
    double f_high = fflo(high);
    if constexpr (extraBhpAtThpLimitOutput) {
        deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + well_.name() +
                              "  low = " + std::to_string(low) +
                              "  high = " + std::to_string(high) +
                              "  f(low) = " + std::to_string(f_low) +
                              "  f(high) = " + std::to_string(f_high));
    }
    int adjustments = 0;
    const int max_adjustments = 10;
    const double adjust_amount = 5.0 * unit::barsa;
    while (f_low * f_high > 0.0 && adjustments < max_adjustments) {
        // Same sign, adjust high to see if we can flip it.
        high += adjust_amount;
        f_high = fflo(high);
        ++adjustments;
    }
    if (f_low * f_high > 0.0) {
        if (f_low > 0.0) {
            // Even at the BHP limit, we are injecting.
            // There will be no solution here, return an
            // empty optional.
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
                                    "Robust bhp(thp) solve failed due to inoperability for well " + well_.name());
            return std::nullopt;
        } else {
            // Still producing, even at high bhp.
            assert(f_high < 0.0);
            bhp_max = high;
        }
    } else {
        // Bisect to find a bhp point where we produce, but
        // not a large amount ('eps' below).
        const double eps = 0.1 * std::fabs(vfp_flo_front);
        const int maxit = 50;
        int it = 0;
        while (std::fabs(f_low) > eps && it < maxit) {
            const double curr = 0.5*(low + high);
            const double f_curr = fflo(curr);
            if (f_curr * f_low > 0.0) {
                low = curr;
                f_low = f_curr;
            } else {
                high = curr;
                f_high = f_curr;
            }
            ++it;
        }
        if (it < maxit) {
            bhp_max = low;
        } else {
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
                                    "Bisect did not find the bhp-point where we produce for well " + well_.name());
            return std::nullopt;
        }
    }
    if constexpr (extraBhpAtThpLimitOutput) {
        deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + well_.name() +
                              "  low = " + std::to_string(low) +
                              "  high = " + std::to_string(high) +
                              "  f(low) = " + std::to_string(f_low) +
                              "  f(high) = " + std::to_string(f_high) +
                              "  bhp_max = " + std::to_string(bhp_max));
    }
    return bhp_max;
}

std::optional<double>
WellBhpThpCalculator::
computeBhpAtThpLimit(const std::function<std::vector<double>(const double)>& frates,
                     const std::function<double(const std::vector<double>)>& fbhp,
                     const std::array<double, 2>& range,
                     DeferredLogger& deferred_logger) const
{
    // Given a VFP function returning bhp as a function of phase
    // rates and thp:
    //     fbhp(rates, thp),
    // a function extracting the particular flow rate used for VFP
    // lookups:
    //     flo(rates)
    // and the inflow function (assuming the reservoir is fixed):
    //     frates(bhp)
    // we want to solve the equation:
    //     fbhp(frates(bhp, thplimit)) - bhp = 0
    // for bhp.
    //
    // This may result in 0, 1 or 2 solutions. If two solutions,
    // the one corresponding to the lowest bhp (and therefore
    // highest rate) should be returned.

    // Define the equation we want to solve.
    auto eq = [&fbhp, &frates](double bhp) {
        return fbhp(frates(bhp)) - bhp;
    };

    // Find appropriate brackets for the solution.
    std::optional<double> approximate_solution;
    double low, high;
    // trying to use bisect way to locate a bracket
    bool finding_bracket = this->bisectBracket(eq, range, low, high, approximate_solution, deferred_logger);

    // based on the origional design, if an approximate solution is suggested, we use this value directly
    // in the long run, we might change it
    if (approximate_solution.has_value()) {
        return *approximate_solution;
    }

    if (!finding_bracket) {
        deferred_logger.debug(" Trying the brute force search to bracket the bhp for last attempt ");
        finding_bracket = this->bruteForceBracket(eq, range, low, high, deferred_logger);
    }

    if (!finding_bracket) {
        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
                                "Robust bhp(thp) solve failed due to not being able to "
                                "bracket the bhp solution with the brute force search for " + well_.name());
        return std::nullopt;
    }

    // Solve for the proper solution in the given interval.
    const int max_iteration = 100;
    const double bhp_tolerance = 0.01 * unit::barsa;
    int iteration = 0;
    try {
        const double solved_bhp = RegulaFalsiBisection<ThrowOnError>::
            solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
        return solved_bhp;
    }
    catch (...) {
        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + well_.name());
        return std::nullopt;
    }
}

bool
WellBhpThpCalculator::
bisectBracket(const std::function<double(const double)>& eq,
              const std::array<double, 2>& range,
              double& low, double& high,
              std::optional<double>& approximate_solution,
              DeferredLogger& deferred_logger) const
{
    bool finding_bracket = false;
    low = range[0];
    high = range[1];

    double eq_high = eq(high);
    double eq_low = eq(low);
    const double eq_bhplimit = eq_low;
    if constexpr (extraBhpAtThpLimitOutput) {
        deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + well_.name() +
                              "  low = " + std::to_string(low) +
                              "  high = " + std::to_string(high) +
                              "  eq(low) = " + std::to_string(eq_low) +
                              "  eq(high) = " + std::to_string(eq_high));
    }
    if (eq_low * eq_high > 0.0) {
        // Failed to bracket the zero.
        // If this is due to having two solutions, bisect until bracketed.
        double abs_low = std::fabs(eq_low);
        double abs_high = std::fabs(eq_high);
        int bracket_attempts = 0;
        const int max_bracket_attempts = 20;
        double interval = high - low;
        const double min_interval = 1.0 * unit::barsa;
        while (eq_low * eq_high > 0.0 && bracket_attempts < max_bracket_attempts && interval > min_interval) {
            if (abs_high < abs_low) {
                low = 0.5 * (low + high);
                eq_low = eq(low);
                abs_low = std::fabs(eq_low);
            } else {
                high = 0.5 * (low + high);
                eq_high = eq(high);
                abs_high = std::fabs(eq_high);
            }
            ++bracket_attempts;
        }

        if (eq_low * eq_high <= 0.) {
            // We have a bracket!
            finding_bracket = true;
            // Now, see if (bhplimit, low) is a bracket in addition to (low, high).
            // If so, that is the bracket we shall use, choosing the solution with the
            // highest flow.
            if (eq_low * eq_bhplimit <= 0.0) {
                high = low;
                low = range[0];
            }
        } else { // eq_low * eq_high > 0.0
            // Still failed bracketing!
            const double limit = 0.1 * unit::barsa;
            if (std::min(abs_low, abs_high) < limit) {
                // Return the least bad solution if less off than 0.1 bar.
                deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_BRACKETING_FAILURE",
                                        "Robust bhp(thp) not solved precisely for well " + well_.name());
                approximate_solution = abs_low < abs_high ? low : high;
            } else {
                    // Return failure.
                deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_BRACKETING_FAILURE",
                                         "Robust bhp(thp) solve failed due to bracketing failure for well " +
                                         well_.name());
            }
        }
    } else {
        finding_bracket = true;
    }
    return finding_bracket;
}

bool
WellBhpThpCalculator::
bruteForceBracket(const std::function<double(const double)>& eq,
                  const std::array<double, 2>& range,
                  double& low, double& high,
                  DeferredLogger& deferred_logger)
{
    bool finding_bracket = false;
    low = range[0];
    high = range[1];
    const int sample_number = 100;
    const double interval = (high - low) / sample_number;
    double eq_low = eq(low);
    double eq_high;
    for (int i = 0; i < sample_number + 1; ++i) {
        high = range[0] + interval * i;
        eq_high = eq(high);
        if (eq_high * eq_low <= 0.) {
            finding_bracket = true;
            break;
        }
        low = high;
        eq_low = eq_high;
    }
    if (finding_bracket) {
        deferred_logger.debug(
                " brute force solve found low " + std::to_string(low) + " with eq_low " + std::to_string(eq_low) +
                " high " + std::to_string(high) + " with eq_high " + std::to_string(eq_high));
    }
    return finding_bracket;
}

} // namespace Opm
