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

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

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

namespace Opm {

template<typename Scalar, typename IndexTraits>
bool WellBhpThpCalculator<Scalar, IndexTraits>::
wellHasTHPConstraints(const SummaryState& summaryState) const
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

template<typename Scalar, typename IndexTraits>
Scalar WellBhpThpCalculator<Scalar, IndexTraits>::
getTHPConstraint(const SummaryState& summaryState) const
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

template<typename Scalar, typename IndexTraits>
Scalar WellBhpThpCalculator<Scalar, IndexTraits>::
mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const
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

template<typename Scalar, typename IndexTraits>
Scalar WellBhpThpCalculator<Scalar, IndexTraits>::
calculateThpFromBhp(const std::vector<Scalar>& rates,
                    const Scalar bhp,
                    const Scalar rho,
                    const std::optional<Scalar>& alq,
                    const Scalar thp_limit,
                    DeferredLogger& deferred_logger) const
{
    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

    const Scalar aqua = rates[IndexTraits::waterPhaseIdx];
    const Scalar liquid = rates[IndexTraits::oilPhaseIdx];
    const Scalar vapour = rates[IndexTraits::gasPhaseIdx];

    // pick the density in the top layer
    Scalar thp = 0.0;
    const int table_id = well_.wellEcl().vfp_table_number();
    if (well_.isInjector()) {
        assert(!alq.has_value());
        const Scalar vfp_ref_depth = well_.vfpProperties()->getInj()->getTable(table_id).getDatumDepth();
        const Scalar dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());
        auto thp_func =
            [this, table_id, aqua, liquid, vapour, dp](
                const Scalar bhp_value, const Scalar pressure_loss) {
                    return this->well_.vfpProperties()->getInj()->thp(
                           table_id, aqua, liquid, vapour, bhp_value + dp - pressure_loss);
                };
        thp = findThpFromBhpIteratively(thp_func, bhp, thp_limit, dp, deferred_logger);
    }
    else if (well_.isProducer()) {
        const Scalar vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(table_id).getDatumDepth();
        const Scalar dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());
        const bool use_vfpexp = well_.useVfpExplicit();
        const Scalar wfr = well_.vfpProperties()->getExplicitWFR(table_id, well_.indexOfWell());
        const Scalar gfr = well_.vfpProperties()->getExplicitGFR(table_id, well_.indexOfWell());
        auto thp_func =
            [this, table_id, aqua, liquid, vapour, dp, &alq, wfr, gfr, use_vfpexp]
                (const Scalar bhp_value, const Scalar pressure_loss) {
                    return this->well_.vfpProperties()->getProd()->thp(
                       table_id, aqua, liquid, vapour, bhp_value + dp - pressure_loss, alq.value(), wfr, gfr, use_vfpexp);
                };
        thp = findThpFromBhpIteratively(thp_func, bhp, thp_limit, dp, deferred_logger);
    }
    else {
        OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
    }
    return thp;
}

template<typename Scalar, typename IndexTraits>
Scalar WellBhpThpCalculator<Scalar, IndexTraits>::
findThpFromBhpIteratively(const std::function<Scalar(const Scalar, const Scalar)>& thp_func,
                          const Scalar bhp,
                          const Scalar thp_limit,
                          const Scalar dp,
                          DeferredLogger& deferred_logger) const
{
    auto pressure_loss = getVfpBhpAdjustment(bhp + dp, thp_limit);
    auto thp = thp_func(bhp, pressure_loss);
    const Scalar tolerance = 1e-5 * unit::barsa;
    int it = 1;
    int max_iterations = 50;
    while (true) {
        if (it > max_iterations) {
            break;
        }
        Scalar thp_prev = thp;
        pressure_loss = getVfpBhpAdjustment(bhp + dp - pressure_loss, thp_prev);
        thp = thp_func(bhp, pressure_loss);
        auto error = std::fabs(thp - thp_prev);
        if (extraThpFromBhpOutput) {
            const std::string msg = fmt::format(
                "findThpFromBhpIteratively(): iteration {}, thp = {}, bhp = {}, "
                "pressure_loss = {}, error = {}", it, thp, bhp+dp-pressure_loss, pressure_loss, error);
            deferred_logger.debug(msg);
        }
        if (std::fabs(thp - thp_prev) < tolerance) {
            break;
        }
        it++;
    }
    return thp;
}

template<typename Scalar, typename IndexTraits>
std::optional<Scalar>
WellBhpThpCalculator<Scalar, IndexTraits>::
computeBhpAtThpLimitProd(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                         const SummaryState& summary_state,
                         const Scalar maxPerfPress,
                         const Scalar rho,
                         const Scalar alq_value,
                         const Scalar thp_limit,
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

    static constexpr int Water = IndexTraits::waterPhaseIdx;
    static constexpr int Oil = IndexTraits::oilPhaseIdx;
    static constexpr int Gas = IndexTraits::gasPhaseIdx;

    // Make the fbhp() function.
    const auto& controls = well_.wellEcl().productionControls(summary_state);
    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);
    const Scalar vfp_ref_depth = table.getDatumDepth();
    const Scalar dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(),
                                                                vfp_ref_depth,
                                                                rho,
                                                                well_.gravity());

    auto fbhp = [this, &controls, thp_limit, dp, alq_value](const std::vector<Scalar>& rates) {
        assert(rates.size() == 3);
        const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number,
                                                                well_.indexOfWell());
        const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number,
                                                                well_.indexOfWell());
        const bool use_vfpexp = well_.useVfpExplicit();
        const Scalar bhp = well_.vfpProperties()->getProd()->bhp(controls.vfp_table_number,
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
    auto flo = [&table](const std::vector<Scalar>& rates) {
        return detail::getFlo(table, rates[Water], rates[Oil], rates[Gas]);
    };

    // Find the bhp-point where production becomes nonzero.
    auto fflo = [&flo, &frates](Scalar bhp) { return flo(frates(bhp)); };
    auto bhp_max = this->bhpMax(fflo,
                                controls.bhp_limit,
                                maxPerfPress,
                                table.getFloAxis().front(),
                                deferred_logger);

    // could not solve for the bhp-point, we could not continue to find the bhp
    if (!bhp_max.has_value()) {
        deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
                              "Robust bhp(thp) solve failed due to not being able to "
                              "find bhp-point where production becomes non-zero for well " + well_.name());
        return std::nullopt;
    }
    const std::array<Scalar, 2> range {static_cast<Scalar>(controls.bhp_limit), *bhp_max};
    return this->computeBhpAtThpLimit(frates, fbhp, range, deferred_logger);
}

template<typename Scalar, typename IndexTraits>
std::optional<Scalar>
WellBhpThpCalculator<Scalar, IndexTraits>::
computeBhpAtThpLimitInj(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                        const SummaryState& summary_state,
                        const Scalar rho,
                        const Scalar flo_rel_tol,
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

template<typename Scalar, typename IndexTraits>
void WellBhpThpCalculator<Scalar, IndexTraits>::
updateThp(const Scalar rho,
          const std::function<Scalar()>& alq_value,
          WellState<Scalar, IndexTraits>& well_state,
          const SummaryState& summary_state,
          DeferredLogger& deferred_logger) const
{
    auto& ws = well_state.well(well_.indexOfWell());

    // When there is no vaild VFP table provided, we set the thp to be zero.
    if (!well_.isVFPActive(deferred_logger)) {
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
    std::vector<Scalar> rates(3, 0.0);

    const auto& pu = well_.phaseUsage();
    //TODO: the following code should be able to make a for loop
    if (pu.phaseIsActive(IndexTraits::waterPhaseIdx)) {
        const int water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        rates[IndexTraits::waterPhaseIdx] = ws.surface_rates[water_pos];
    }
    if (pu.phaseIsActive(IndexTraits::oilPhaseIdx)) {
        const int oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        rates[IndexTraits::oilPhaseIdx] = ws.surface_rates[oil_pos];
    }
    if (pu.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        const int gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        rates[IndexTraits::gasPhaseIdx] = ws.surface_rates[gas_pos];
    }
    const std::optional<Scalar> alq = this->well_.isProducer() ? std::optional<Scalar>(alq_value()) : std::nullopt;
    const Scalar thp_limit = well_.getTHPConstraint(summary_state);
    ws.thp = this->calculateThpFromBhp(rates, ws.bhp, rho, alq, thp_limit, deferred_logger);
}

template<typename Scalar, typename IndexTraits>
template<class EvalWell>
EvalWell WellBhpThpCalculator<Scalar, IndexTraits>::
calculateBhpFromThp(const WellState<Scalar, IndexTraits>& well_state,
                    const std::vector<EvalWell>& rates,
                    const Well& well,
                    const SummaryState& summaryState,
                    const Scalar rho,
                    DeferredLogger& deferred_logger) const
{
    // TODO: when well is under THP control, the BHP is dependent on the rates,
    // the well rates is also dependent on the BHP, so it might need to do some iteration.
    // However, when group control is involved, change of the rates might impacts other wells
    // so iterations on a higher level will be required. Some investigation might be needed when
    // we face problems under THP control.

    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

    static constexpr int Water = IndexTraits::waterPhaseIdx;
    static constexpr int Oil = IndexTraits::oilPhaseIdx;
    static constexpr int Gas = IndexTraits::gasPhaseIdx;

    const EvalWell aqua = rates[Water];
    const EvalWell liquid = rates[Oil];
    const EvalWell vapour = rates[Gas];
    const Scalar thp_limit = well_.getTHPConstraint(summaryState);
    Scalar vfp_ref_depth;
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
    Scalar bhp_tab_double_value;
    if constexpr (std::is_same_v<EvalWell, Scalar>) {
        bhp_tab_double_value = bhp_tab;
    }
    else {  // EvalWell and bhp_tab is of type DenseAd::Evaluation<double,...,...>
        bhp_tab_double_value = bhp_tab.value();
    }
    const auto bhp_adjustment = getVfpBhpAdjustment(bhp_tab_double_value, thp_limit);
    const Scalar dp_hydro = wellhelpers::computeHydrostaticCorrection(well_.refDepth(),
                                                                      vfp_ref_depth,
                                                                      rho,
                                                                      well_.gravity());
    return bhp_tab - dp_hydro + bhp_adjustment;
}

template<typename Scalar, typename IndexTraits>
Scalar WellBhpThpCalculator<Scalar, IndexTraits>::
calculateMinimumBhpFromThp(const WellState<Scalar, IndexTraits>& well_state,
                           const Well& well,
                           const SummaryState& summaryState,
                           const Scalar rho) const
{
    assert(well_.isProducer()); // only producers can go here for now

    const Scalar thp_limit = well_.getTHPConstraint(summaryState);

    const auto& controls = well.productionControls(summaryState);
    const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
    const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());

    const Scalar bhp_min = well_.vfpProperties()->getProd()->minimumBHP(controls.vfp_table_number,
                                                                        thp_limit, wfr, gfr,
                                                                        well_.getALQ(well_state));

    const Scalar vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number).getDatumDepth();
    const auto bhp_adjustment = getVfpBhpAdjustment(bhp_min, thp_limit);
    const Scalar dp_hydro = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth,
                                                                      rho, well_.gravity());
    return bhp_min - dp_hydro + bhp_adjustment;
}

template<typename Scalar, typename IndexTraits>
Scalar WellBhpThpCalculator<Scalar, IndexTraits>::
getVfpBhpAdjustment(const Scalar bhp_tab, const Scalar thp_limit) const
{
    return well_.wellEcl().getWVFPDP().getPressureLoss(bhp_tab, thp_limit);
}

template<typename Scalar, typename IndexTraits>
template<class ErrorPolicy>
std::optional<Scalar>
WellBhpThpCalculator<Scalar, IndexTraits>::
computeBhpAtThpLimitInjImpl(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                            const SummaryState& summary_state,
                            const Scalar rho,
                            const Scalar flo_rel_tol,
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

    static constexpr int Water = IndexTraits::waterPhaseIdx;
    static constexpr int Oil = IndexTraits::oilPhaseIdx;
    static constexpr int Gas = IndexTraits::gasPhaseIdx;

    // Make the fbhp() function.
    const auto& controls = well_.wellEcl().injectionControls(summary_state);
    const auto& table = well_.vfpProperties()->getInj()->getTable(controls.vfp_table_number);
    const Scalar vfp_ref_depth = table.getDatumDepth();
    const Scalar thp_limit = this->getTHPConstraint(summary_state);
    const Scalar dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(),
                                                                vfp_ref_depth,
                                                                rho, well_.gravity());
    auto fbhp = [this, &controls, thp_limit, dp](const std::vector<Scalar>& rates) {
        assert(rates.size() == 3);
        const auto bhp = well_.vfpProperties()->getInj()
                ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], thp_limit);
        return bhp - dp + getVfpBhpAdjustment(bhp, thp_limit);
    };

    // Make the flo() function.
    auto flo = [&table](const std::vector<Scalar>& rates) {
        return detail::getFlo(table, rates[Water], rates[Oil], rates[Gas]);
    };

    // Get the flo samples, add extra samples at low rates and bhp
    // limit point if necessary.
    std::vector<double> flo_samples = table.getFloAxis();
    if (flo_samples[0] > 0.0) {
        const double f0 = flo_samples[0];
        flo_samples.insert(flo_samples.begin(), { f0/20.0, f0/10.0, f0/5.0, f0/2.0 });
    }
    const Scalar flo_bhp_limit = flo(frates(controls.bhp_limit));
    if (flo_samples.back() < flo_bhp_limit) {
        flo_samples.push_back(flo_bhp_limit);
    }

    // Find bhp values for inflow relation corresponding to flo samples.
    std::vector<Scalar> bhp_samples;
    for (Scalar flo_sample : flo_samples) {
        if (flo_sample > flo_bhp_limit) {
            // We would have to go over the bhp limit to obtain a
            // flow of this magnitude. We associate all such flows
            // with simply the bhp limit. The first one
            // encountered is considered valid, the rest not. They
            // are therefore skipped.
            bhp_samples.push_back(controls.bhp_limit);
            break;
        }
        auto eq = [&flo, &frates, flo_sample](Scalar bhp) {
            return flo(frates(bhp)) - flo_sample;
        };
        // TODO: replace hardcoded low/high limits.
        const Scalar low = 10.0 * unit::barsa;
        const Scalar high = 800.0 * unit::barsa;
        const Scalar flo_tolerance = flo_rel_tol * std::fabs(flo_samples.back());
        try {
            int iteration = 0;
            const Scalar solved_bhp = RegulaFalsiBisection<ErrorPolicy>::
                    solve(eq, low, high, max_iteration, flo_tolerance, iteration);
            bhp_samples.push_back(solved_bhp);
        }
        catch (...) {
            // Use previous value (or max value if at start) if we failed.
            bhp_samples.push_back(bhp_samples.empty() ? low : bhp_samples.back());
            deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE_EXTRACT_SAMPLES",
                                  "Robust bhp(thp) solve failed extracting bhp values at flo samples for well " + well_.name());
        }
    }

    // Find bhp values for VFP relation corresponding to flo samples.
    const int num_samples = bhp_samples.size(); // Note that this can be smaller than flo_samples.size()
    std::vector<Scalar> fbhp_samples(num_samples);
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
        const Scalar curr = fbhp_samples[ii] - bhp_samples[ii];
        const Scalar next = fbhp_samples[ii + 1] - bhp_samples[ii + 1];
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
    auto eq = [&fbhp, &frates](Scalar bhp) {
        return fbhp(frates(bhp)) - bhp;
    };
    // TODO: replace hardcoded low/high limits.
    const Scalar low = bhp_samples[sign_change_index + 1];
    const Scalar high = bhp_samples[sign_change_index];
    const Scalar bhp_tolerance = 0.01 * unit::barsa;
    if (low == high) {
        // We are in the high flow regime where the bhp_samples
        // are all equal to the bhp_limit.
        assert(low == controls.bhp_limit);
        deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + well_.name());
        return std::nullopt;
    }
    try {
        int iteration = 0;
        const Scalar solved_bhp = RegulaFalsiBisection<ErrorPolicy>::
                solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
        if constexpr (extraBhpAtThpLimitOutput) {
            OpmLog::debug("*****    " + well_.name() + "    solved_bhp = " + std::to_string(solved_bhp)
                          + "    flo_bhp_limit = " + std::to_string(flo_bhp_limit));
        }
        return solved_bhp;
    }
    catch (...) {
        deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + well_.name());
        return std::nullopt;
    }
}

template<typename Scalar, typename IndexTraits>
std::optional<Scalar>
WellBhpThpCalculator<Scalar, IndexTraits>::
bhpMax(const std::function<Scalar(const Scalar)>& fflo,
       const Scalar bhp_limit,
       const Scalar maxPerfPress,
       const Scalar vfp_flo_front,
       DeferredLogger& deferred_logger) const
{
    // Find the bhp-point where production becomes nonzero.
    Scalar bhp_max = 0.0;
    Scalar low = bhp_limit;
    Scalar high = maxPerfPress + 1.0 * unit::barsa;
    Scalar f_low = fflo(low);
    Scalar f_high = fflo(high);
    if constexpr (extraBhpAtThpLimitOutput) {
        deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + well_.name() +
                              "  low = " + std::to_string(low) +
                              "  high = " + std::to_string(high) +
                              "  f(low) = " + std::to_string(f_low) +
                              "  f(high) = " + std::to_string(f_high));
    }
    int adjustments = 0;
    const int max_adjustments = 10;
    const Scalar adjust_amount = 5.0 * unit::barsa;
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
            deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
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
        const Scalar eps = 0.1 * std::fabs(vfp_flo_front);
        const int maxit = 50;
        int it = 0;
        while (std::fabs(f_low) > eps && it < maxit) {
            const Scalar curr = 0.5*(low + high);
            const Scalar f_curr = fflo(curr);
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
            deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
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

template<typename Scalar, typename IndexTraits>
std::optional<Scalar>
WellBhpThpCalculator<Scalar, IndexTraits>::
computeBhpAtThpLimit(const std::function<std::vector<Scalar>(const Scalar)>& frates,
                     const std::function<Scalar(const std::vector<Scalar>)>& fbhp,
                     const std::array<Scalar, 2>& range,
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
    auto eq = [&fbhp, &frates](Scalar bhp) {
        return fbhp(frates(bhp)) - bhp;
    };

    // Find appropriate brackets for the solution.
    std::optional<Scalar> approximate_solution;
    Scalar low, high;
    // trying to use bisect way to locate a bracket
    bool finding_bracket = this->bisectBracket(eq, range, low, high,
                                               approximate_solution, deferred_logger);

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
        deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE_INOPERABLE",
                                "Robust bhp(thp) solve failed due to not being able to "
                                "bracket the bhp solution with the brute force search for " + well_.name());
        return std::nullopt;
    }

    // Solve for the proper solution in the given interval.
    const int max_iteration = 100;
    const Scalar bhp_tolerance = 0.01 * unit::barsa;
    try {
        int iteration = 0;
        const Scalar solved_bhp = RegulaFalsiBisection<ThrowOnError>::
            solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
        return solved_bhp;
    }
    catch (...) {
        deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + well_.name());
        return std::nullopt;
    }
}

template<typename Scalar, typename IndexTraits>
bool WellBhpThpCalculator<Scalar, IndexTraits>::
bisectBracket(const std::function<Scalar(const Scalar)>& eq,
              const std::array<Scalar, 2>& range,
              Scalar& low, Scalar& high,
              std::optional<Scalar>& approximate_solution,
              DeferredLogger& deferred_logger) const
{
    bool finding_bracket = false;
    low = range[0];
    high = range[1];

    Scalar eq_high = eq(high);
    Scalar eq_low = eq(low);
    const Scalar eq_bhplimit = eq_low;
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
        Scalar abs_low = std::fabs(eq_low);
        Scalar abs_high = std::fabs(eq_high);
        int bracket_attempts = 0;
        const int max_bracket_attempts = 20;
        Scalar interval = high - low;
        const Scalar min_interval = 1.0 * unit::barsa;
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
            const Scalar limit = 0.1 * unit::barsa;
            if (std::min(abs_low, abs_high) < limit) {
                // Return the least bad solution if less off than 0.1 bar.
                deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE_BRACKETING_FAILURE",
                                      "Robust bhp(thp) not solved precisely for well " + well_.name());
                approximate_solution = abs_low < abs_high ? low : high;
            } else {
                    // Return failure.
                deferred_logger.debug("FAILED_ROBUST_BHP_THP_SOLVE_BRACKETING_FAILURE",
                                      "Robust bhp(thp) solve failed due to bracketing failure for well " +
                                         well_.name());
            }
        }
    } else {
        finding_bracket = true;
    }
    return finding_bracket;
}

template<typename Scalar, typename IndexTraits>
bool WellBhpThpCalculator<Scalar, IndexTraits>::
bruteForceBracket(const std::function<Scalar(const Scalar)>& eq,
                  const std::array<Scalar, 2>& range,
                  Scalar& low, Scalar& high,
                  DeferredLogger& deferred_logger)
{
    bool bracket_found = false;
    low = range[0];
    high = range[1];
    const int sample_number = 200;
    const Scalar interval = (high - low) / sample_number;
    Scalar eq_low = eq(low);
    Scalar eq_high = 0.0;
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

template<typename Scalar, typename IndexTraits>
bool WellBhpThpCalculator<Scalar, IndexTraits>::
isStableSolution(const WellState<Scalar, IndexTraits>& well_state,
                 const Well& well,
                 const std::vector<Scalar>& rates,
                 const SummaryState& summaryState) const
{
    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.
    assert(well_.isProducer()); // only valid for producers

    static constexpr int Water = IndexTraits::waterPhaseIdx;
    static constexpr int Oil = IndexTraits::oilPhaseIdx;
    static constexpr int Gas = IndexTraits::gasPhaseIdx;

    const Scalar aqua = rates[Water];
    const Scalar liquid = rates[Oil];
    const Scalar vapour = rates[Gas];
    const Scalar thp = well_.getTHPConstraint(summaryState);

    const auto& controls = well.productionControls(summaryState);
    const auto& wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
    const auto& gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());

    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);
    const bool use_vfpexplicit = well_.useVfpExplicit();

    auto bhp = VFPHelpers<double>::bhp(table, aqua, liquid, vapour, thp,
                                       well_.getALQ(well_state), wfr, gfr,
                                       use_vfpexplicit);

    if (bhp.dflo >= 0) {
        return true;
    } else {    // maybe check if ipr is available
        const auto ipr = getFloIPR(well_state, well, summaryState);
        return bhp.dflo + 1.0 / ipr.second >= 0;
    }
}

template<typename Scalar, typename IndexTraits>
std::optional<Scalar> WellBhpThpCalculator<Scalar, IndexTraits>::
estimateStableBhp(const WellState<Scalar, IndexTraits>& well_state,
                  const Well& well,
                  const std::vector<Scalar>& rates,
                  const Scalar rho,
                  const SummaryState& summaryState) const
{
    // Given a *converged* well_state with ipr, estimate bhp of the stable solution
    const auto& controls = well.productionControls(summaryState);
    const Scalar thp = well_.getTHPConstraint(summaryState);
    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);

    const Scalar aqua = rates[IndexTraits::waterPhaseIdx];
    const Scalar liquid = rates[IndexTraits::oilPhaseIdx];
    const Scalar vapour = rates[IndexTraits::gasPhaseIdx];
    Scalar flo = detail::getFlo(table, aqua, liquid, vapour);
    Scalar wfr, gfr;
    if (well_.useVfpExplicit() || -flo < table.getFloAxis().front()) {
        wfr =  well_.vfpProperties()->getExplicitWFR(controls.vfp_table_number, well_.indexOfWell());
        gfr = well_.vfpProperties()->getExplicitGFR(controls.vfp_table_number, well_.indexOfWell());
    } else {
        wfr = detail::getWFR(table, aqua, liquid, vapour);
        gfr = detail::getGFR(table, aqua, liquid, vapour);
    }

    auto ipr = getFloIPR(well_state, well, summaryState);

    const Scalar vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number).getDatumDepth();

    const Scalar dp_hydro = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth,
                                                                      rho, well_.gravity());
    auto bhp_adjusted = [this, &thp, &dp_hydro](const Scalar bhp) {
           return bhp - dp_hydro + getVfpBhpAdjustment(bhp, thp);
       };
    const auto retval = VFPHelpers<double>::intersectWithIPR(table, thp, wfr, gfr,
                                                             well_.getALQ(well_state),
                                                             ipr.first, ipr.second,
                                                             bhp_adjusted);
    if (retval.has_value()) {
        // returned pair is (flo, bhp)
        return retval.value().second;
    } else {
        return std::nullopt;
    }
}

template<typename Scalar, typename IndexTraits>
std::pair<Scalar, Scalar> WellBhpThpCalculator<Scalar, IndexTraits>::
getFloIPR(const WellState<Scalar, IndexTraits>& well_state,
          const Well& well,
          const SummaryState& summary_state) const
{
    // Convert ipr_a's and ipr_b's to our particular choice of FLO
    const auto& controls = well.productionControls(summary_state);
    const auto& table = well_.vfpProperties()->getProd()->getTable(controls.vfp_table_number);
    const auto& ipr_a = well_state.well(well_.indexOfWell()).implicit_ipr_a;
    const auto& pu = well_.phaseUsage();
    const Scalar& aqua_a = pu.phaseIsActive(IndexTraits::waterPhaseIdx) ?
                       ipr_a[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] : 0.0;
    const Scalar& liquid_a = pu.phaseIsActive(IndexTraits::oilPhaseIdx) ?
                           ipr_a[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] : 0.0;
    const Scalar& vapour_a = pu.phaseIsActive(IndexTraits::gasPhaseIdx) ?
                             ipr_a[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] : 0.0;
    const auto& ipr_b = well_state.well(well_.indexOfWell()).implicit_ipr_b;
    const Scalar& aqua_b = pu.phaseIsActive(IndexTraits::waterPhaseIdx) ?
                           ipr_b[pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)] : 0.0;
    const Scalar& liquid_b = pu.phaseIsActive(IndexTraits::oilPhaseIdx) ?
                           ipr_b[pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)] : 0.0;
    const Scalar& vapour_b = pu.phaseIsActive(IndexTraits::gasPhaseIdx) ?
                             ipr_b[pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)] : 0.0;
    // The getFlo helper is indended to pick one or add two of the phase rates (depending on FLO-type),
    // but we can equally use it to pick/add the corresponding ipr_a, ipr_b
    return std::make_pair(detail::getFlo(table, aqua_a, liquid_a, vapour_a),
                          detail::getFlo(table, aqua_b, liquid_b, vapour_b));
}

template<typename Scalar, typename IndexTraits>
bool
WellBhpThpCalculator<Scalar, IndexTraits>::
bruteForceBracketCommonTHP(const std::function<Scalar(const Scalar)>& eq,
                           const std::array<Scalar, 2>& range,
                           Scalar& low, Scalar& high,
                           std::optional<Scalar>& approximate_solution,
                           const Scalar& limit,
                           DeferredLogger& deferred_logger)
{
    bool bracket_found = false;
    low = range[0];
    high = range[1];
    const int sample_number = 300;
    const Scalar interval = (high - low) / sample_number;
    Scalar eq_low = eq(low);
    Scalar eq_high = 0.0;
    for (int i = 0; i < sample_number + 1; ++i) {
        high = range[0] + interval * i;
        eq_high = eq(high);
        if ((std::fabs(eq_high) < limit)) {
            approximate_solution = high;
            break;
        }
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

template<typename Scalar, typename IndexTraits>
bool
WellBhpThpCalculator<Scalar, IndexTraits>::
bruteForceBracketCommonTHP(const std::function<Scalar(const Scalar)>& eq,
                           Scalar& min_thp, Scalar& max_thp)
{
    bool bracket_found = false;
    constexpr int sample_number = 1000;
    constexpr Scalar interval = 1E5;
    Scalar eq_low = eq(min_thp);
    for (int i = 0; i < sample_number + 1; ++i) {
        max_thp = min_thp + interval * i;
        const Scalar eq_high = eq(max_thp);
        if (eq_high * eq_low <= 0.) {
            bracket_found = true;
            min_thp = max_thp - interval;
            break;
        }
        eq_low = eq_high;
    }
    return bracket_found;
}


#define INSTANTIATE(T,...)                                   \
    template __VA_ARGS__                                     \
    WellBhpThpCalculator<T, BlackOilDefaultFluidSystemIndices>::                                \
        calculateBhpFromThp(const WellState<T, BlackOilDefaultFluidSystemIndices>&,             \
                            const std::vector<__VA_ARGS__>&, \
                            const Well&,                     \
                            const SummaryState&,             \
                            const T,                         \
                            DeferredLogger&) const;

#define INSTANTIATE_TYPE(T)                      \
    template class WellBhpThpCalculator<T, BlackOilDefaultFluidSystemIndices>;      \
    INSTANTIATE(T,T)                             \
    INSTANTIATE(T,DenseAd::Evaluation<T,3,0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T,4,0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T,5,0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T,6,0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T,7,0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T,8,0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T,9,0u>)   \
    INSTANTIATE(T,DenseAd::Evaluation<T,10,0u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,4u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,5u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,6u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,7u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,8u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,9u>)  \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,10u>) \
    INSTANTIATE(T,DenseAd::Evaluation<T,-1,11u>)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
