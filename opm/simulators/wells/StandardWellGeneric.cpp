/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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
#include <opm/simulators/wells/StandardWellGeneric.hpp>

#include <opm/common/utility/numeric/RootFinders.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>
#include <stdexcept>

namespace Opm
{

template<class Scalar>
StandardWellGeneric<Scalar>::
StandardWellGeneric(int Bhp,
                    const WellInterfaceGeneric& baseif)
    : baseif_(baseif)
    , perf_densities_(baseif_.numPerfs())
    , perf_pressure_diffs_(baseif_.numPerfs())
    , parallelB_(duneB_, baseif_.parallelWellInfo())
    , Bhp_(Bhp)
{
    duneB_.setBuildMode(OffDiagMatWell::row_wise);
    duneC_.setBuildMode(OffDiagMatWell::row_wise);
    invDuneD_.setBuildMode(DiagMatWell::row_wise);
}


template<class Scalar>
double
StandardWellGeneric<Scalar>::
relaxationFactorRate(const std::vector<double>& primary_variables,
                     const BVectorWell& dwells)
{
    double relaxation_factor = 1.0;
    static constexpr int WQTotal = 0;

    // For injector, we only check the total rates to avoid sign change of rates
    const double original_total_rate = primary_variables[WQTotal];
    const double newton_update = dwells[0][WQTotal];
    const double possible_update_total_rate = primary_variables[WQTotal] - newton_update;

    // 0.8 here is a experimental value, which remains to be optimized
    // if the original rate is zero or possible_update_total_rate is zero, relaxation_factor will
    // always be 1.0, more thoughts might be needed.
    if (original_total_rate * possible_update_total_rate < 0.) { // sign changed
        relaxation_factor = std::abs(original_total_rate / newton_update) * 0.8;
    }

    assert(relaxation_factor >= 0.0 && relaxation_factor <= 1.0);

    return relaxation_factor;
}

template<class Scalar>
double
StandardWellGeneric<Scalar>::
relaxationFactorFraction(const double old_value,
                         const double dx)
{
    assert(old_value >= 0. && old_value <= 1.0);

    double relaxation_factor = 1.;

    // updated values without relaxation factor
    const double possible_updated_value = old_value - dx;

    // 0.95 is an experimental value remains to be optimized
    if (possible_updated_value < 0.0) {
        relaxation_factor = std::abs(old_value / dx) * 0.95;
    } else if (possible_updated_value > 1.0) {
        relaxation_factor = std::abs((1. - old_value) / dx) * 0.95;
    }
    // if possible_updated_value is between 0. and 1.0, then relaxation_factor
    // remains to be one

    assert(relaxation_factor >= 0. && relaxation_factor <= 1.);

    return relaxation_factor;
}

template<class Scalar>
double
StandardWellGeneric<Scalar>::
calculateThpFromBhp(const WellState &well_state,
                    const std::vector<double>& rates,
                    const double bhp,
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
    if (baseif_.isInjector()) {
        const int table_id = baseif_.wellEcl().vfp_table_number();
        const double vfp_ref_depth = baseif_.vfpProperties()->getInj()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, getRho(), baseif_.gravity());

        thp = baseif_.vfpProperties()->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
    }
    else if (baseif_.isProducer()) {
        const int table_id = baseif_.wellEcl().vfp_table_number();
        const double alq = baseif_.getALQ(well_state);
        const double vfp_ref_depth = baseif_.vfpProperties()->getProd()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, getRho(), baseif_.gravity());

        thp = baseif_.vfpProperties()->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
    }
    else {
        OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
    }

    return thp;
}

template<class Scalar>
void
StandardWellGeneric<Scalar>::
computeConnectionPressureDelta()
{
    // Algorithm:

    // We'll assume the perforations are given in order from top to
    // bottom for each well.  By top and bottom we do not necessarily
    // mean in a geometric sense (depth), but in a topological sense:
    // the 'top' perforation is nearest to the surface topologically.
    // Our goal is to compute a pressure delta for each perforation.

    // 1. Compute pressure differences between perforations.
    //    dp_perf will contain the pressure difference between a
    //    perforation and the one above it, except for the first
    //    perforation for each well, for which it will be the
    //    difference to the reference (bhp) depth.

    const int nperf = baseif_.numPerfs();
    perf_pressure_diffs_.resize(nperf, 0.0);
    auto z_above = baseif_.parallelWellInfo().communicateAboveValues(baseif_.refDepth(), baseif_.perfDepth());

    for (int perf = 0; perf < nperf; ++perf) {
        const double dz = baseif_.perfDepth()[perf] - z_above[perf];
        perf_pressure_diffs_[perf] = dz * perf_densities_[perf] * baseif_.gravity();
    }

    // 2. Compute pressure differences to the reference point (bhp) by
    //    accumulating the already computed adjacent pressure
    //    differences, storing the result in dp_perf.
    //    This accumulation must be done per well.
    const auto beg = perf_pressure_diffs_.begin();
    const auto end = perf_pressure_diffs_.end();
    baseif_.parallelWellInfo().partialSumPerfValues(beg, end);
}

template<class Scalar>
std::optional<double>
StandardWellGeneric<Scalar>::
computeBhpAtThpLimitProdWithAlq(const std::function<std::vector<double>(const double)>& frates,
                                const SummaryState& summary_state,
                                DeferredLogger& deferred_logger,
                                double maxPerfPress,
                                double alq_value) const
{
    return baseif_.computeBhpAtThpLimitProdCommon(frates, summary_state, maxPerfPress, getRho(), alq_value, deferred_logger);
}

template<class Scalar>
std::optional<double>
StandardWellGeneric<Scalar>::
computeBhpAtThpLimitInj(const std::function<std::vector<double>(const double)>& frates,
                        const SummaryState& summary_state,
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
    const double dp = wellhelpers::computeHydrostaticCorrection(baseif_.refDepth(), vfp_ref_depth, getRho(), baseif_.gravity());
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
        const int max_iteration = 50;
        const double flo_tolerance = 1e-6 * std::fabs(flo_samples.back());
        int iteration = 0;
        try {
            const double solved_bhp = RegulaFalsiBisection<>::
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
    const int max_iteration = 50;
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
        const double solved_bhp = RegulaFalsiBisection<>::
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

template<class Scalar>
void
StandardWellGeneric<Scalar>::
checkConvergenceControlEq(const WellState& well_state,
                          ConvergenceReport& report,
                          DeferredLogger& deferred_logger,
                          const double max_residual_allowed) const
{
    double control_tolerance = 0.;
    using CR = ConvergenceReport;
    CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;

    const int well_index = baseif_.indexOfWell();
    const auto& ws = well_state.well(well_index);
    if (baseif_.wellIsStopped()) {
        ctrltype = CR::WellFailure::Type::ControlRate;
        control_tolerance = 1.e-6; // use smaller tolerance for zero control?
    }
    else if (baseif_.isInjector() )
    {
        auto current = ws.injection_cmode;
        switch(current) {
        case Well::InjectorCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = 1.e4; // 0.1 bar
            break;
        case Well::InjectorCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = 1.e3; // 0.01 bar
            break;
        case Well::InjectorCMode::RATE:
        case Well::InjectorCMode::RESV:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-4; //
            break;
        case Well::InjectorCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-6; //
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
        }
    }
    else if (baseif_.isProducer() )
    {
        auto current = ws.production_cmode;
        switch(current) {
        case Well::ProducerCMode::THP:
            ctrltype = CR::WellFailure::Type::ControlTHP;
            control_tolerance = 1.e4; // 0.1 bar
            break;
        case Well::ProducerCMode::BHP:
            ctrltype = CR::WellFailure::Type::ControlBHP;
            control_tolerance = 1.e3; // 0.01 bar
            break;
        case Well::ProducerCMode::ORAT:
        case Well::ProducerCMode::WRAT:
        case Well::ProducerCMode::GRAT:
        case Well::ProducerCMode::LRAT:
        case Well::ProducerCMode::RESV:
        case Well::ProducerCMode::CRAT:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-4; // smaller tolerance for rate control
            break;
        case Well::ProducerCMode::GRUP:
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-6; // smaller tolerance for rate control
            break;
        default:
            OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << baseif_.name(), deferred_logger);
        }
    }

    const double well_control_residual = std::abs(this->resWell_[0][Bhp_]);
    const int dummy_component = -1;
    if (std::isnan(well_control_residual)) {
        report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, baseif_.name()});
    } else if (well_control_residual > max_residual_allowed * 10.) {
        report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, baseif_.name()});
    } else if ( well_control_residual > control_tolerance) {
        report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, baseif_.name()});
    }
}

template<class Scalar>
void
StandardWellGeneric<Scalar>::
checkConvergencePolyMW(const std::vector<double>& res,
                       ConvergenceReport& report,
                       const double maxResidualAllowed) const
{
  if (baseif_.isInjector()) {
      //  checking the convergence of the perforation rates
      const double wat_vel_tol = 1.e-8;
      const int dummy_component = -1;
      using CR = ConvergenceReport;
      const auto wat_vel_failure_type = CR::WellFailure::Type::MassBalance;
      for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
          const double wat_vel_residual = res[Bhp_ + 1 + perf];
          if (std::isnan(wat_vel_residual)) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::NotANumber, dummy_component, baseif_.name()});
          } else if (wat_vel_residual > maxResidualAllowed * 10.) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::TooLarge, dummy_component, baseif_.name()});
          } else if (wat_vel_residual > wat_vel_tol) {
              report.setWellFailed({wat_vel_failure_type, CR::Severity::Normal, dummy_component, baseif_.name()});
          }
      }

      // checking the convergence of the skin pressure
      const double pskin_tol = 1000.; // 1000 pascal
      const auto pskin_failure_type = CR::WellFailure::Type::Pressure;
      for (int perf = 0; perf < baseif_.numPerfs(); ++perf) {
          const double pskin_residual = res[Bhp_ + 1 + perf + baseif_.numPerfs()];
          if (std::isnan(pskin_residual)) {
              report.setWellFailed({pskin_failure_type, CR::Severity::NotANumber, dummy_component, baseif_.name()});
          } else if (pskin_residual > maxResidualAllowed * 10.) {
              report.setWellFailed({pskin_failure_type, CR::Severity::TooLarge, dummy_component, baseif_.name()});
          } else if (pskin_residual > pskin_tol) {
              report.setWellFailed({pskin_failure_type, CR::Severity::Normal, dummy_component, baseif_.name()});
          }
      }
  }
}


template<class Scalar>
void
StandardWellGeneric<Scalar>::
getNumBlocks(unsigned int& numBlocks) const
{
    numBlocks = duneB_.nonzeroes();
}

template class StandardWellGeneric<double>;

}
