/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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
#include <opm/simulators/wells/MultisegmentWellGeneric.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace Opm
{

template<typename Scalar>
MultisegmentWellGeneric<Scalar>::
MultisegmentWellGeneric(WellInterfaceGeneric& baseif)
    : baseif_(baseif)
    , segment_perforations_(numberOfSegments())
    , segment_inlets_(numberOfSegments())
    , segment_depth_diffs_(numberOfSegments(), 0.0)
    , perforation_segment_depth_diffs_(baseif_.numPerfs(), 0.0)
    , bhp_scaling_(1.0)
    , rate_scaling_(1.0)
    , bhp_control_scaling_(1.0)  
{
    // since we decide to use the WellSegments from the well parser. we can reuse a lot from it.
    // for other facilities needed but not available from parser, we need to process them here

    // initialize the segment_perforations_ and update perforation_segment_depth_diffs_
    const WellConnections& completion_set = baseif_.wellEcl().getConnections();
    // index of the perforation within wells struct
    // there might be some perforations not active, which causes the number of the perforations in
    // well_ecl_ and wells struct different
    // the current implementation is a temporary solution for now, it should be corrected from the parser
    // side
    int i_perf_wells = 0;
    baseif.perfDepth().resize(baseif_.numPerfs(), 0.);
    for (size_t perf = 0; perf < completion_set.size(); ++perf) {
        const Connection& connection = completion_set.get(perf);
        if (connection.state() == Connection::State::OPEN) {
            const int segment_index = segmentNumberToIndex(connection.segment());
            segment_perforations_[segment_index].push_back(i_perf_wells);
            baseif.perfDepth()[i_perf_wells] = connection.depth();
            const double segment_depth = segmentSet()[segment_index].depth();
            perforation_segment_depth_diffs_[i_perf_wells] = baseif.perfDepth()[i_perf_wells] - segment_depth;
            i_perf_wells++;
        }
    }

    // initialize the segment_inlets_
    for (int seg = 0; seg < numberOfSegments(); ++seg) {
        const Segment& segment = segmentSet()[seg];
        const int segment_number = segment.segmentNumber();
        const int outlet_segment_number = segment.outletSegment();
        if (outlet_segment_number > 0) {
            const int segment_index = segmentNumberToIndex(segment_number);
            const int outlet_segment_index = segmentNumberToIndex(outlet_segment_number);
            segment_inlets_[outlet_segment_index].push_back(segment_index);
        }
    }

    // calculating the depth difference between the segment and its oulet_segments
    // for the top segment, we will make its zero unless we find other purpose to use this value
    for (int seg = 1; seg < numberOfSegments(); ++seg) {
        const double segment_depth = segmentSet()[seg].depth();
        const int outlet_segment_number = segmentSet()[seg].outletSegment();
        const Segment& outlet_segment = segmentSet()[segmentNumberToIndex(outlet_segment_number)];
        const double outlet_depth = outlet_segment.depth();
        segment_depth_diffs_[seg] = segment_depth - outlet_depth;
    }
}

template<typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
scaleSegmentRatesWithWellRates(WellState& well_state) const
{
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    auto& segment_rates = segments.rates;
    for (int phase = 0; phase < baseif_.numPhases(); ++phase) {
        const double unscaled_top_seg_rate = segment_rates[phase];
        const double well_phase_rate = ws.surface_rates[phase];
        if (std::abs(unscaled_top_seg_rate) > 1e-12) {
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                segment_rates[baseif_.numPhases() * seg + phase] *= well_phase_rate / unscaled_top_seg_rate;
            }
        } else {
            // Due to various reasons, the well/top segment rate can be zero for this phase.
            // We can not scale this rate directly. The following approach is used to initialize the segment rates.
            double sumTw = 0;
            for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
                sumTw += baseif_.wellIndex()[perf];
            }

            // only handling this specific phase
            constexpr double num_single_phase = 1;
            std::vector<double> perforation_rates(num_single_phase * baseif_.numPerfs(), 0.0);
            const double perf_phaserate_scaled = ws.surface_rates[phase] / sumTw;
            for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
                perforation_rates[perf] = baseif_.wellIndex()[perf] * perf_phaserate_scaled;
            }

            std::vector<double> rates;
            WellState::calculateSegmentRates(segment_inlets_, segment_perforations_, perforation_rates, num_single_phase, 0, rates);
            for (int seg = 0; seg < numberOfSegments(); ++seg) {
                segment_rates[baseif_.numPhases() * seg + phase] = rates[seg];
            }
        }
    }
}

template <typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
scaleSegmentPressuresWithBhp(WellState& well_state) const
{
    auto& ws = well_state.well(baseif_.indexOfWell());
    auto& segments = ws.segments;
    segments.scale_pressure(ws.bhp);
}

template<typename Scalar>
const WellSegments&
MultisegmentWellGeneric<Scalar>::
segmentSet() const
{
    return baseif_.wellEcl().getSegments();
}

template <typename Scalar>
int
MultisegmentWellGeneric<Scalar>::
numberOfSegments() const
{
    return segmentSet().size();
}

template <typename Scalar>
WellSegments::CompPressureDrop
MultisegmentWellGeneric<Scalar>::
compPressureDrop() const
{
    return segmentSet().compPressureDrop();
}


template<typename Scalar>
int
MultisegmentWellGeneric<Scalar>::
segmentNumberToIndex(const int segment_number) const
{
    return segmentSet().segmentNumberToIndex(segment_number);
}

template<typename Scalar>
double
MultisegmentWellGeneric<Scalar>::
calculateThpFromBhp(const std::vector<double>& rates,
                    const double bhp,
                    const double rho,
                    DeferredLogger& deferred_logger) const
{
    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const double aqua = rates[Water];
    const double liquid = rates[Oil];
    const double vapour = rates[Gas];

    double thp = 0.0;
    if (baseif_.isInjector()) {
        const int table_id = baseif_.wellEcl().vfp_table_number();
        const double vfp_ref_depth = baseif_.vfpProperties()->getInj()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, rho, baseif_.gravity());

        thp = baseif_.vfpProperties()->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
    }
    else if (baseif_.isProducer()) {
        const int table_id = baseif_.wellEcl().vfp_table_number();
        const double alq = baseif_.wellEcl().alq_value();
        const double vfp_ref_depth = baseif_.vfpProperties()->getProd()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, rho, baseif_.gravity());

        thp = baseif_.vfpProperties()->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
    }
    else {
        OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
    }

    return thp;
}

template<typename Scalar>
void
MultisegmentWellGeneric<Scalar>::
detectOscillations(const std::vector<double>& measure_history,
                   const int it,
                   bool& oscillate,
                   bool& stagnate) const
{
    if ( it < 2 ) {
        oscillate = false;
        stagnate = false;
        return;
    }

    stagnate = true;
    const double F0 = measure_history[it];
    const double F1 = measure_history[it - 1];
    const double F2 = measure_history[it - 2];
    const double d1 = std::abs((F0 - F2) / F0);
    const double d2 = std::abs((F0 - F1) / F0);

    const double oscillaton_rel_tol = 0.2;
    oscillate = (d1 < oscillaton_rel_tol) && (oscillaton_rel_tol < d2);

    const double stagnation_rel_tol = 1.e-2;
    stagnate = std::abs((F1 - F2) / F2) <= stagnation_rel_tol;
}

template<typename Scalar>
std::optional<double>
MultisegmentWellGeneric<Scalar>::
computeBhpAtThpLimitInj(const std::function<std::vector<double>(const double)>& frates,
                        const SummaryState& summary_state,
                        const double rho,
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
    const auto& controls = baseif_.wellEcl().injectionControls(summary_state);
    const auto& table = baseif_.vfpProperties()->getInj()->getTable(controls.vfp_table_number);
    const double vfp_ref_depth = table.getDatumDepth();
    const double thp_limit = baseif_.getTHPConstraint(summary_state);
    const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, rho, baseif_.gravity());
    auto fbhp = [this, &controls, thp_limit, dp](const std::vector<double>& rates) {
        assert(rates.size() == 3);
        return baseif_.vfpProperties()->getInj()
                ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], thp_limit) - dp;
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
        const int max_iteration = 100;
        const double flo_tolerance = 0.05 * std::fabs(flo_samples.back());
        int iteration = 0;
        try {
            const double solved_bhp = RegulaFalsiBisection<WarnAndContinueOnError>::
                    solve(eq, low, high, max_iteration, flo_tolerance, iteration);
            bhp_samples.push_back(solved_bhp);
        }
        catch (...) {
            // Use previous value (or max value if at start) if we failed.
            bhp_samples.push_back(bhp_samples.empty() ? low : bhp_samples.back());
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_EXTRACT_SAMPLES",
                                    "Robust bhp(thp) solve failed extracting bhp values at flo samples for well " + baseif_.name());
        }
    }

    // Find bhp values for VFP relation corresponding to flo samples.
    const int num_samples = bhp_samples.size(); // Note that this can be smaller than flo_samples.size()
    std::vector<double> fbhp_samples(num_samples);
    for (int ii = 0; ii < num_samples; ++ii) {
        fbhp_samples[ii] = fbhp(frates(bhp_samples[ii]));
    }
// #define EXTRA_THP_DEBUGGING
#ifdef EXTRA_THP_DEBUGGING
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
#endif // EXTRA_THP_DEBUGGING

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
    const int max_iteration = 100;
    const double bhp_tolerance = 0.01 * unit::barsa;
    int iteration = 0;
    if (low == high) {
        // We are in the high flow regime where the bhp_samples
        // are all equal to the bhp_limit.
        assert(low == controls.bhp_limit);
        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + baseif_.name());
        return std::nullopt;
    }
    try {
        const double solved_bhp = RegulaFalsiBisection<WarnAndContinueOnError>::
                solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
#ifdef EXTRA_THP_DEBUGGING
        OpmLog::debug("*****    " + name() + "    solved_bhp = " + std::to_string(solved_bhp)
                      + "    flo_bhp_limit = " + std::to_string(flo_bhp_limit));
#endif // EXTRA_THP_DEBUGGING
        return solved_bhp;
    }
    catch (...) {
        deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                "Robust bhp(thp) solve failed for well " + baseif_.name());
        return std::nullopt;
    }
}

template<typename Scalar>
std::optional<double>
MultisegmentWellGeneric<Scalar>::
computeBhpAtThpLimitProdWithAlq(
                   const std::function<std::vector<double>(const double)>& frates,
                   const SummaryState& summary_state,
                   const double maxPerfPress,
                   const double rho,
                   DeferredLogger& deferred_logger,
                   double alq_value) const
{
   return baseif_.computeBhpAtThpLimitProdCommon(frates, summary_state, maxPerfPress, rho, alq_value, deferred_logger);
}

template<typename Scalar>
bool
MultisegmentWellGeneric<Scalar>::
frictionalPressureLossConsidered() const
{
    // HF- and HFA needs to consider frictional pressure loss
    return (segmentSet().compPressureDrop() != WellSegments::CompPressureDrop::H__);
}

template<typename Scalar>
bool
MultisegmentWellGeneric<Scalar>::
accelerationalPressureLossConsidered() const
{
    return (segmentSet().compPressureDrop() == WellSegments::CompPressureDrop::HFA);
}

template class MultisegmentWellGeneric<double>;

} // namespace Opm
