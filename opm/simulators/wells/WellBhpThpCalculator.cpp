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
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <cassert>
#include <cmath>
#include <optional>

#include <fmt/format.h>

static constexpr bool extraBhpAtThpLimitOutput = false;
static constexpr bool extraThpFromBhpOutput = false;

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
                                                 const std::optional<double>& alq,
                                                 const double thp_limit,
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
    const int table_id = well_.wellEcl().vfp_table_number();
    if (well_.isInjector()) {
        assert(!alq.has_value());
        const double vfp_ref_depth = well_.vfpProperties()->getInj()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());
        auto thp_func =
            [this, table_id, aqua, liquid, vapour, dp](
                const double bhp_value, const double pressure_loss) {
                    return this->well_.vfpProperties()->getInj()->thp(
                           table_id, aqua, liquid, vapour, bhp_value + dp - pressure_loss);
                };
        thp = findThpFromBhpIteratively(thp_func, bhp, thp_limit, dp, deferred_logger);
    }
    else if (well_.isProducer()) {
        const double vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());
        auto thp_func =
            [this, table_id, aqua, liquid, vapour, dp, &alq]
                (const double bhp_value, const double pressure_loss) {
                    return this->well_.vfpProperties()->getProd()->thp(
                       table_id, aqua, liquid, vapour, bhp_value + dp - pressure_loss, alq.value());
                };
        thp = findThpFromBhpIteratively(thp_func, bhp, thp_limit, dp, deferred_logger);
    }
    else {
        OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
    }
    return thp;
}

double
WellBhpThpCalculator::
findThpFromBhpIteratively(
    const std::function<double(const double, const double)>& thp_func,
    const double bhp,
    const double thp_limit,
    const double dp,
    DeferredLogger& deferred_logger) const
{
    auto pressure_loss = getVfpBhpAdjustment(bhp + dp, thp_limit);
    auto thp = thp_func(bhp, pressure_loss);
    const double tolerance = 1e-5 * unit::barsa;
    bool do_iterate = true;
    int it = 1;
    int max_iterations = 50;
    while(do_iterate) {
        if (it > max_iterations) {
            break;
        }
        double thp_prev = thp;
        pressure_loss = getVfpBhpAdjustment(bhp + dp - pressure_loss, thp_prev);
        thp = thp_func(bhp, pressure_loss);
        auto error = std::fabs(thp-thp_prev);
        if (extraThpFromBhpOutput) {
            const std::string msg = fmt::format(
                "findThpFromBhpIteratively(): iteration {}, thp = {}, bhp = {}, "
                "pressure_loss = {}, error = {}", it, thp, bhp+dp-pressure_loss, pressure_loss, error);
            deferred_logger.debug(msg);
        }
        if (std::fabs(thp-thp_prev) < tolerance) {
            break;
        }
        it++;
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
                         const double thp_limit,
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
    const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());

    auto fbhp = [this, &controls, thp_limit, dp, alq_value](const std::vector<double>& rates) {
        assert(rates.size() == 3);
        const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
        const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());
        const bool use_vfpexp = well_.useVfpExplicit();
        const double bhp = well_.vfpProperties()->getProd()->bhp(controls.vfp_table_number,
                                                     rates[Water],
                                                     rates[Oil],
                                                     rates[Gas],
                                                     thp_limit,
                                                     alq_value,
                                                     wfr,
                                                     gfr,
                                                     use_vfpexp);
        return bhp - dp + getVfpBhpAdjustment(bhp, thp_limit);
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
computeBhpAtThpLimitInj(const std::function<std::vector<double>(const double)>& frates,
                        const SummaryState& summary_state,
                        const double rho,
                        const double flo_rel_tol,
                        const int max_iteration,
                        const bool throwOnError,
                        DeferredLogger& deferred_logger) const
{
    if (throwOnError) {
        return computeBhpAtThpLimitInjImpl<ThrowOnError>(frates, summary_state,
                                                         rho, flo_rel_tol,
                                                         max_iteration, deferred_logger);
    } else {
        return computeBhpAtThpLimitInjImpl<WarnAndContinueOnError>(frates, summary_state,
                                                                   rho, flo_rel_tol,
                                                                   max_iteration, deferred_logger);
    }
}

void WellBhpThpCalculator::updateThp(const double rho,
                                     const bool stop_or_zero_rate_target,
                                     const std::function<double()>& alq_value,
                                     const std::array<unsigned,3>& active,
                                     WellState& well_state,
                                     const SummaryState& summary_state,
                                     DeferredLogger& deferred_logger) const
{
    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;
    auto& ws = well_state.well(well_.indexOfWell());

    // When there is no vaild VFP table provided, we set the thp to be zero.
    if (!well_.isVFPActive(deferred_logger) || stop_or_zero_rate_target) {
        ws.thp = 0;
        return;
    }
    // For THP controlled wells, we know the thp value
    bool thp_controlled = well_.isInjector() ? ws.injection_cmode == Well::InjectorCMode::THP:
                                               ws.production_cmode == Well::ProducerCMode::THP;
    if (thp_controlled) {
        return;
    }

    // the well is under other control types, we calculate the thp based on bhp and rates
    std::vector<double> rates(3, 0.0);

    const PhaseUsage& pu = well_.phaseUsage();
    if (active[Water]) {
        rates[ Water ] = ws.surface_rates[pu.phase_pos[ Water ] ];
    }
    if (active[Oil]) {
        rates[ Oil ] = ws.surface_rates[pu.phase_pos[ Oil ] ];
    }
    if (active[Gas]) {
        rates[ Gas ] = ws.surface_rates[pu.phase_pos[ Gas ] ];
    }
    const std::optional<double> alq = this->well_.isProducer() ? std::optional<double>(alq_value()) : std::nullopt;
    const double thp_limit = well_.getTHPConstraint(summary_state);
    ws.thp = this->calculateThpFromBhp(rates, ws.bhp, rho, alq, thp_limit, deferred_logger);
}

template<class EvalWell>
EvalWell WellBhpThpCalculator::
calculateBhpFromThp(const WellState& well_state,
                    const std::vector<EvalWell>& rates,
                    const Well& well,
                    const SummaryState& summaryState,
                    const double rho,
                    DeferredLogger& deferred_logger) const
{
    // TODO: when well is under THP control, the BHP is dependent on the rates,
    // the well rates is also dependent on the BHP, so it might need to do some iteration.
    // However, when group control is involved, change of the rates might impacts other wells
    // so iterations on a higher level will be required. Some investigation might be needed when
    // we face problems under THP control.

    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    const EvalWell aqua = rates[Water];
    const EvalWell liquid = rates[Oil];
    const EvalWell vapour = rates[Gas];
    const double thp_limit = well_.getTHPConstraint(summaryState);
    double vfp_ref_depth;
    EvalWell bhp_tab;
    if (well_.isInjector() )
    {
        const auto& controls = well.injectionControls(summaryState);
        vfp_ref_depth = well_.vfpProperties()->getInj()->getTable(controls.vfp_table_number).getDatumDepth();
        bhp_tab = well_.vfpProperties()->getInj()->bhp(
               controls.vfp_table_number, aqua, liquid, vapour, thp_limit);
    }
    else if (well_.isProducer()) {
        const auto& controls = well.productionControls(summaryState);
        vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number).getDatumDepth();
        const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
        const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());
        const bool use_vfpexplicit = well_.useVfpExplicit();

        bhp_tab = well_.vfpProperties()->getProd()->bhp(controls.vfp_table_number,
                                                      aqua, liquid, vapour,
                                                      thp_limit,
                                                      well_.getALQ(well_state),
                                                      wfr, gfr, use_vfpexplicit);
    }
    else {
        OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER for well " + well_.name(), deferred_logger);
    }
    double bhp_tab_double_value;
    if constexpr (std::is_same_v<EvalWell, double>) {
        bhp_tab_double_value = bhp_tab;
    }
    else {  // EvalWell and bhp_tab is of type DenseAd::Evaluation<double,...,...>
        bhp_tab_double_value = bhp_tab.value();
    }
    const auto bhp_adjustment = getVfpBhpAdjustment(bhp_tab_double_value, thp_limit);
    const double dp_hydro = wellhelpers::computeHydrostaticCorrection(
            well_.refDepth(), vfp_ref_depth, rho, well_.gravity());
    return bhp_tab - dp_hydro + bhp_adjustment;
}

double
WellBhpThpCalculator::
calculateMinimumBhpFromThp(const WellState& well_state,
                           const Well& well,
                           const SummaryState& summaryState,
                           const double rho) const
{
    assert(well_.isProducer()); // only producers can go here for now

    const double thp_limit = well_.getTHPConstraint(summaryState);
    
    const auto& controls = well.productionControls(summaryState);
    const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
    const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());

    const double bhp_min = well_.vfpProperties()->getProd()->minimumBHP(controls.vfp_table_number,
                                                                        thp_limit, wfr, gfr, 
                                                                        well_.getALQ(well_state));
 
    const double vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number).getDatumDepth();
    const auto bhp_adjustment = getVfpBhpAdjustment(bhp_min, thp_limit);
    const double dp_hydro = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, 
                                                                      rho, well_.gravity());
    return bhp_min - dp_hydro + bhp_adjustment;
}

double WellBhpThpCalculator::getVfpBhpAdjustment(const double bhp_tab, const double thp_limit) const
{
    return well_.wellEcl().getWVFPDP().getPressureLoss(bhp_tab, thp_limit);
}

template<class ErrorPolicy>
std::optional<double>
WellBhpThpCalculator::
computeBhpAtThpLimitInjImpl(const std::function<std::vector<double>(const double)>& frates,
                            const SummaryState& summary_state,
                            const double rho,
                            const double flo_rel_tol,
                            const int max_iteration,
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
    // highest rate) is returned.
    //
    // In order to detect these situations, we will find piecewise
    // linear approximations both to the inverse of the frates
    // function and to the fbhp function.
    //
    // We first take the FLO sample points of the VFP curve, and
    // find the corresponding bhp values by solving the equation:
    //     flo(frates(bhp)) - flo_sample = 0
    // for bhp, for each flo_sample. The resulting (flo_sample,
    // bhp_sample) values give a piecewise linear approximation to
    // the true inverse inflow function, at the same flo values as
    // the VFP data.
    //
    // Then we extract a piecewise linear approximation from the
    // multilinear fbhp() by evaluating it at the flo_sample
    // points, with fractions given by the frates(bhp_sample)
    // values.
    //
    // When we have both piecewise linear curves defined on the
    // same flo_sample points, it is easy to distinguish between
    // the 0, 1 or 2 solution cases, and obtain the right interval
    // in which to solve for the solution we want (with highest
    // flow in case of 2 solutions).

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    // Make the fbhp() function.
    const auto& controls = well_.wellEcl().injectionControls(summary_state);
    const auto& table = well_.vfpProperties()->getInj()->getTable(controls.vfp_table_number);
    const double vfp_ref_depth = table.getDatumDepth();
    const double thp_limit = this->getTHPConstraint(summary_state);
    const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(),
                                                                vfp_ref_depth,
                                                                rho, well_.gravity());
    auto fbhp = [this, &controls, thp_limit, dp](const std::vector<double>& rates) {
        assert(rates.size() == 3);
        const auto bhp = well_.vfpProperties()->getInj()
                ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], thp_limit);
        return bhp - dp + getVfpBhpAdjustment(bhp, thp_limit);
    };

    // Make the flo() function.
    auto flo = [&table](const std::vector<double>& rates) {
        return detail::getFlo(table, rates[Water], rates[Oil], rates[Gas]);
    };

    // Get the flo samples, add extra samples at low rates and bhp
    // limit point if necessary.
    std::vector<double> flo_samples = table.getFloAxis();
    if (flo_samples[0] > 0.0) {
        const double f0 = flo_samples[0];
        flo_samples.insert(flo_samples.begin(), { f0/20.0, f0/10.0, f0/5.0, f0/2.0 });
    }
    const double flo_bhp_limit = flo(frates(controls.bhp_limit));
    if (flo_samples.back() < flo_bhp_limit) {
        flo_samples.push_back(flo_bhp_limit);
    }

    // Find bhp values for inflow relation corresponding to flo samples.
    std::vector<double> bhp_samples;
    for (double flo_sample : flo_samples) {
        if (flo_sample > flo_bhp_limit) {
            // We would have to go over the bhp limit to obtain a
            // flow of this magnitude. We associate all such flows
            // with simply the bhp limit. The first one
            // encountered is considered valid, the rest not. They
            // are therefore skipped.
            bhp_samples.push_back(controls.bhp_limit);
            break;
        }
        auto eq = [&flo, &frates, flo_sample](double bhp) {
            return flo(frates(bhp)) - flo_sample;
        };
        // TODO: replace hardcoded low/high limits.
        const double low = 10.0 * unit::barsa;
        const double high = 800.0 * unit::barsa;
        const double flo_tolerance = flo_rel_tol * std::fabs(flo_samples.back());
        int iteration = 0;
        try {
            const double solved_bhp = RegulaFalsiBisection<ErrorPolicy>::
                    solve(eq, low, high, max_iteration, flo_tolerance, iteration);
            bhp_samples.push_back(solved_bhp);
        }
        catch (...) {
            // Use previous value (or max value if at start) if we failed.
            bhp_samples.push_back(bhp_samples.empty() ? low : bhp_samples.back());
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_EXTRACT_SAMPLES",
                                    "Robust bhp(thp) solve failed extracting bhp values at flo samples for well " + well_.name());
        }
    }

    // Find bhp values for VFP relation corresponding to flo samples.
    const int num_samples = bhp_samples.size(); // Note that this can be smaller than flo_samples.size()
    std::vector<double> fbhp_samples(num_samples);
    for (int ii = 0; ii < num_samples; ++ii) {
        fbhp_samples[ii] = fbhp(frates(bhp_samples[ii]));
    }
    if constexpr (extraBhpAtThpLimitOutput) {
        std::string dbgmsg;
        dbgmsg += "flo: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(flo_samples[ii]);
        }
        dbgmsg += "\nbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(bhp_samples[ii]);
        }
        dbgmsg += "\nfbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(fbhp_samples[ii]);
        }
        OpmLog::debug(dbgmsg);
    }

    // Look for sign changes for the (fbhp_samples - bhp_samples) piecewise linear curve.
    // We only look at the valid
    int sign_change_index = -1;
    for (int ii = 0; ii < num_samples - 1; ++ii) {
        const double curr = fbhp_samples[ii] - bhp_samples[ii];
        const double next = fbhp_samples[ii + 1] - bhp_samples[ii + 1];
        if (curr * next < 0.0) {
            // Sign change in the [ii, ii + 1] interval.
            sign_change_index = ii; // May overwrite, thereby choosing the highest-flo solution.
        }
    }

    // Handle the no solution case.
    if (sign_change_index == -1) {
        return std::nullopt;
    }

    // Solve for the proper solution in the given interval.
    auto eq = [&fbhp, &frates](double bhp) {
        return fbhp(frates(bhp)) - bhp;
    };
    // TODO: replace hardcoded low/high limits.
    const double low = bhp_samples[sign_change_index + 1];
    const double high = bhp_samples[sign_change_index];
    const double bhp_tolerance = 0.01 * unit::barsa;
    int iteration = 0;
    if (low == high) {
        // We are in the high flow regime where the bhp_samples
        // are all equal to the bhp_limit.
        assert(low == controls.bhp_limit);
        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + well_.name());
        return std::nullopt;
    }
    try {
        const double solved_bhp = RegulaFalsiBisection<ErrorPolicy>::
                solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
        if constexpr (extraBhpAtThpLimitOutput) {
            OpmLog::debug("*****    " + well_.name() + "    solved_bhp = " + std::to_string(solved_bhp)
                          + "    flo_bhp_limit = " + std::to_string(flo_bhp_limit));
        }
        return solved_bhp;
    }
    catch (...) {
        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + well_.name());
        return std::nullopt;
    }
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
            interval = high - low;
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
    bool bracket_found = false;
    low = range[0];
    high = range[1];
    const int sample_number = 200;
    const double interval = (high - low) / sample_number;
    double eq_low = eq(low);
    double eq_high = 0.0;
    for (int i = 0; i < sample_number + 1; ++i) {
        high = range[0] + interval * i;
        eq_high = eq(high);
        if (eq_high * eq_low <= 0.) {
            bracket_found = true;
            break;
        }
        low = high;
        eq_low = eq_high;
    }
    if (bracket_found) {
        deferred_logger.debug(
                " brute force solve found low " + std::to_string(low) + " with eq_low " + std::to_string(eq_low) +
                " high " + std::to_string(high) + " with eq_high " + std::to_string(eq_high));
    }
    return bracket_found;
}

bool WellBhpThpCalculator::
isStableSolution(const WellState& well_state,
                 const Well& well,
                 const std::vector<double>& rates,
                 const SummaryState& summaryState,
                 DeferredLogger& deferred_logger) const
{
    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.
    assert(well_.isProducer()); // only valid for producers 

    static constexpr int Gas = BlackoilPhases::Vapour;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Water = BlackoilPhases::Aqua;

    const double aqua = rates[Water];
    const double liquid = rates[Oil];
    const double vapour = rates[Gas];
    const double thp = well_.getTHPConstraint(summaryState);

    const auto& controls = well.productionControls(summaryState);
    const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
    const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());

    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);
    const bool use_vfpexplicit = well_.useVfpExplicit();

    detail::VFPEvaluation bhp = detail::bhp(table, aqua, liquid, vapour, thp, well_.getALQ(well_state), wfr, gfr, use_vfpexplicit);

    bhp.value = bhp.value + getVfpBhpAdjustment(bhp.value, thp);

    if (bhp.dflo >= 0) {
        return true;
    } else {    // maybe check if ipr is available
        const auto ipr = getFloIPR(well_state, well, summaryState);
        return bhp.dflo + 1/ipr.second >= 0; 
    }                  
}

std::optional<double> WellBhpThpCalculator::
estimateStableBhp(const WellState& well_state,
                  const Well& well,
                  const std::vector<double>& rates,
                  const double rho,
                  const SummaryState& summaryState,
                  DeferredLogger& deferred_logger) const
{   
    // Given a *converged* well_state with ipr, estimate bhp of the stable solution 
    const auto& controls = well.productionControls(summaryState);
    const double thp = well_.getTHPConstraint(summaryState);
    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);

    const double aqua = rates[BlackoilPhases::Aqua];
    const double liquid = rates[BlackoilPhases::Liquid];
    const double vapour = rates[BlackoilPhases::Vapour];
    double flo = detail::getFlo(table, aqua, liquid, vapour);
    double wfr, gfr;
    if (well_.useVfpExplicit() || -flo < table.getFloAxis().front()) {
        wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
        gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());
    } else {
        wfr = detail::getWFR(table, aqua, liquid, vapour);  
        gfr = detail::getGFR(table, aqua, liquid, vapour);   
    }
    if (wfr <= 0.0 && gfr <= 0.0) {
        // warning message
    }
    auto ipr = getFloIPR(well_state, well, summaryState);
    
    const double vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number).getDatumDepth();

    const double dp_hydro = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, 
                                                                      rho, well_.gravity());
    if (ipr.first <= 0.0) {
        // error message
    }
    auto bhp_adjusted = [this, &thp, &dp_hydro](const double bhp) {
           return bhp - dp_hydro + getVfpBhpAdjustment(bhp, thp);
       };
    const auto retval = detail::intersectWithIPR(table, thp, wfr, gfr, well_.getALQ(well_state), ipr.first, ipr.second, bhp_adjusted);
    if (retval.has_value()) {
        // returned pair is (flo, bhp)
        return retval.value().second;
    } else {
        return std::nullopt;
    }
}    

std::pair<double, double> WellBhpThpCalculator::
getFloIPR(const WellState& well_state,
          const Well& well, 
          const SummaryState& summary_state) const 
{
    std::pair<double,double>retval(0.0, 0.0);
    const auto& controls = well.productionControls(summary_state);
    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);
    const auto& pu = well_.phaseUsage();
    const auto& ipr_a= well_state.well(well_.indexOfWell()).implicit_ipr_a;
    const auto& aqua_a = pu.phase_used[BlackoilPhases::Aqua]? ipr_a[pu.phase_pos[BlackoilPhases::Aqua]]:0.0;
    const auto& liquid_a = pu.phase_used[BlackoilPhases::Liquid]? ipr_a[pu.phase_pos[BlackoilPhases::Liquid]]:0.0;
    const auto& vapour_a = pu.phase_used[BlackoilPhases::Vapour]? ipr_a[pu.phase_pos[BlackoilPhases::Vapour]]:0.0;
    const auto& ipr_b= well_state.well(well_.indexOfWell()).implicit_ipr_b;
    const auto& aqua_b = pu.phase_used[BlackoilPhases::Aqua]? ipr_b[pu.phase_pos[BlackoilPhases::Aqua]]:0.0;
    const auto& liquid_b = pu.phase_used[BlackoilPhases::Liquid]? ipr_b[pu.phase_pos[BlackoilPhases::Liquid]]:0.0;
    const auto& vapour_b = pu.phase_used[BlackoilPhases::Vapour]? ipr_b[pu.phase_pos[BlackoilPhases::Vapour]]:0.0;
    return std::make_pair(detail::getFlo(table, aqua_a, liquid_a, vapour_a), 
                          detail::getFlo(table, aqua_b, liquid_b, vapour_b));
}

#define INSTANCE(...) \
template __VA_ARGS__ WellBhpThpCalculator:: \
calculateBhpFromThp<__VA_ARGS__>(const WellState&, \
                                 const std::vector<__VA_ARGS__>&, \
                                 const Well&, \
                                 const SummaryState&, \
                                 const double, \
                                 DeferredLogger&) const;

INSTANCE(double)
INSTANCE(DenseAd::Evaluation<double,3,0u>)
INSTANCE(DenseAd::Evaluation<double,4,0u>)
INSTANCE(DenseAd::Evaluation<double,5,0u>)
INSTANCE(DenseAd::Evaluation<double,6,0u>)
INSTANCE(DenseAd::Evaluation<double,7,0u>)
INSTANCE(DenseAd::Evaluation<double,8,0u>)
INSTANCE(DenseAd::Evaluation<double,9,0u>)
INSTANCE(DenseAd::Evaluation<double,10,0u>)
INSTANCE(DenseAd::Evaluation<double,-1,4u>)
INSTANCE(DenseAd::Evaluation<double,-1,5u>)
INSTANCE(DenseAd::Evaluation<double,-1,6u>)
INSTANCE(DenseAd::Evaluation<double,-1,7u>)
INSTANCE(DenseAd::Evaluation<double,-1,8u>)
INSTANCE(DenseAd::Evaluation<double,-1,9u>)
INSTANCE(DenseAd::Evaluation<double,-1,10u>)
INSTANCE(DenseAd::Evaluation<double,-1,11u>)
} // namespace Opm
