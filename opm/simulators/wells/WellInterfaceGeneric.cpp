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
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/common/utility/numeric/RootFinders.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace Opm
{

WellInterfaceGeneric::WellInterfaceGeneric(const Well& well,
                                           const ParallelWellInfo& pw_info,
                                           const int time_step,
                                           const int pvtRegionIdx,
                                           const int num_components,
                                           const int num_phases,
                                           const int index_of_well,
                                           const std::vector<PerforationData>& perf_data)
      : well_ecl_(well)
      , parallel_well_info_(pw_info)
      , current_step_(time_step)
      , pvtRegionIdx_(pvtRegionIdx)
      , num_components_(num_components)
      , number_of_phases_(num_phases)
      , index_of_well_(index_of_well)
      , perf_data_(&perf_data)
      , ipr_a_(num_components)
      , ipr_b_(num_components)
{
    assert(well.name()==pw_info.name());
    assert(std::is_sorted(perf_data.begin(), perf_data.end(),
                          [](const auto& perf1, const auto& perf2){
        return perf1.ecl_index < perf2.ecl_index;
    }));
    if (time_step < 0) {
        OPM_THROW(std::invalid_argument, "Negtive time step is used to construct WellInterface");
    }

    ref_depth_ = well.getRefDepth();

    // We do not want to count SHUT perforations here, so
    // it would be wrong to use wells.getConnections().size().
    number_of_perforations_ = perf_data.size();

    // perforations related
    {
        well_cells_.resize(number_of_perforations_);
        well_index_.resize(number_of_perforations_);
        saturation_table_number_.resize(number_of_perforations_);
        int perf = 0;
        for (const auto& pd : perf_data) {
            well_cells_[perf] = pd.cell_index;
            well_index_[perf] = pd.connection_transmissibility_factor;
            saturation_table_number_[perf] = pd.satnum_id;
            ++perf;
        }
    }

    // initialization of the completions mapping
    initCompletions();

    well_efficiency_factor_ = 1.0;

    this->wellStatus_ = Well::Status::OPEN;
    if (well.getStatus() == Well::Status::STOP) {
        this->wellStatus_ = Well::Status::STOP;
    }

    wsolvent_ = 0.0;

    well_control_log_.clear();
}

// Currently the VFP calculations requires three-phase input data, see
//  the documentation for keyword VFPPROD and its implementation in
//  VFPProdProperties.cpp. However, by setting the gas flow rate to a dummy
//  value in VFPPROD record 5 (GFR values) and supplying a dummy input value
//  for the gas rate to the methods in VFPProdProperties.cpp, we can extend
//  the VFP calculations to the two-phase oil-water case.
void WellInterfaceGeneric::adaptRatesForVFP(std::vector<double>& rates) const
{
    const auto& pu = this->phaseUsage();
    if (pu.num_phases == 2) {
        if (    pu.phase_used[BlackoilPhases::Aqua] == 1
             && pu.phase_used[BlackoilPhases::Liquid] == 1
             && pu.phase_used[BlackoilPhases::Vapour] == 0)
        {
            assert(rates.size() == 2);
            rates.push_back(0.0);  // set gas rate to zero
        }
        else {
            throw std::logic_error("Two-phase VFP calculation only "
                                   "supported for oil and water");
        }
    }
}

const std::vector<PerforationData>& WellInterfaceGeneric::perforationData() const
{
    return *perf_data_;
}

const std::string& WellInterfaceGeneric::name() const
{
    return well_ecl_.name();
}

bool WellInterfaceGeneric::isInjector() const
{
    return well_ecl_.isInjector();
}

bool WellInterfaceGeneric::isProducer() const
{
    return well_ecl_.isProducer();
}

int WellInterfaceGeneric::indexOfWell() const
{
    return index_of_well_;
}

bool WellInterfaceGeneric::getAllowCrossFlow() const
{
    return well_ecl_.getAllowCrossFlow();
}

const Well& WellInterfaceGeneric::wellEcl() const
{
    return well_ecl_;
}

const PhaseUsage& WellInterfaceGeneric::phaseUsage() const
{
    assert(phase_usage_ != nullptr);

    return *phase_usage_;
}

double WellInterfaceGeneric::wsolvent() const
{
    return wsolvent_;
}

double WellInterfaceGeneric::rsRvInj() const
{
    return well_ecl_.getInjectionProperties().rsRvInj;
}

bool WellInterfaceGeneric::wellHasTHPConstraints(const SummaryState& summaryState) const
{
    if (dynamic_thp_limit_) {
        return true;
    }

    if (well_ecl_.isInjector()) {
        const auto controls = well_ecl_.injectionControls(summaryState);
        if (controls.hasControl(Well::InjectorCMode::THP))
            return true;
    }

    if (well_ecl_.isProducer( )) {
        const auto controls = well_ecl_.productionControls(summaryState);
        if (controls.hasControl(Well::ProducerCMode::THP))
            return true;
    }

    return false;

}

double WellInterfaceGeneric::mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const
{
    if (well_ecl_.isInjector()) {
        const auto& controls = well_ecl_.injectionControls(summaryState);
        return controls.bhp_limit;
    }

    if (well_ecl_.isProducer( )) {
        const auto& controls = well_ecl_.productionControls(summaryState);
        return controls.bhp_limit;
    }

    return 0.0;
}

double WellInterfaceGeneric::getTHPConstraint(const SummaryState& summaryState) const
{
    if (dynamic_thp_limit_) {
        return *dynamic_thp_limit_;
    }
    if (well_ecl_.isInjector()) {
        const auto& controls = well_ecl_.injectionControls(summaryState);
        return controls.thp_limit;
    }

    if (well_ecl_.isProducer( )) {
        const auto& controls = well_ecl_.productionControls(summaryState);
        return controls.thp_limit;
    }

    return 0.0;
}

bool WellInterfaceGeneric::underPredictionMode() const
{
    return well_ecl_.predictionMode();
}

void WellInterfaceGeneric::initCompletions()
{
    assert(completions_.empty() );

    const WellConnections& connections = well_ecl_.getConnections();
    const std::size_t num_conns = connections.size();

    int num_active_connections = 0;
    auto my_next_perf = perf_data_->begin();
    for (std::size_t c = 0; c < num_conns; ++c) {
        if (my_next_perf == perf_data_->end())
        {
            break;
        }
        if (my_next_perf->ecl_index > c)
        {
            continue;
        }
        assert(my_next_perf->ecl_index == c);
        if (connections[c].state() == Connection::State::OPEN) {
            completions_[connections[c].complnum()].push_back(num_active_connections++);
        }
        ++my_next_perf;
    }
    assert(my_next_perf == perf_data_->end());
}

void WellInterfaceGeneric::closeCompletions(const WellTestState& wellTestState)
{
    const auto& connections = well_ecl_.getConnections();
    int perfIdx = 0;
    for (const auto& connection : connections) {
        if (connection.state() == Connection::State::OPEN) {
            if (wellTestState.completion_is_closed(name(), connection.complnum())) {
                this->well_index_[perfIdx] = 0.0;
            }
            perfIdx++;
        }
    }
}

void WellInterfaceGeneric::setVFPProperties(const VFPProperties* vfp_properties_arg)
{
    vfp_properties_ = vfp_properties_arg;
}

void WellInterfaceGeneric::setGuideRate(const GuideRate* guide_rate_arg)
{
    guide_rate_ = guide_rate_arg;
}

void WellInterfaceGeneric::setWellEfficiencyFactor(const double efficiency_factor)
{
    well_efficiency_factor_ = efficiency_factor;
}

void WellInterfaceGeneric::setRepRadiusPerfLength()
{
    const int nperf = number_of_perforations_;

    perf_rep_radius_.clear();
    perf_length_.clear();
    bore_diameters_.clear();

    perf_rep_radius_.reserve(nperf);
    perf_length_.reserve(nperf);
    bore_diameters_.reserve(nperf);

    const WellConnections& connections = well_ecl_.getConnections();
    const std::size_t num_conns = connections.size();
    int num_active_connections = 0;
    auto my_next_perf = perf_data_->begin();
    for (std::size_t c = 0; c < num_conns; ++c) {
        if (my_next_perf == perf_data_->end())
        {
            break;
        }
        if (my_next_perf->ecl_index > c)
        {
            continue;
        }
        assert(my_next_perf->ecl_index == c);
        const auto& connection = connections[c];
        if (connection.state() == Connection::State::OPEN) {
            double radius = connection.rw();
            double re = connection.re(); // area equivalent radius of the grid block
            double perf_length = connection.connectionLength(); // the length of the well perforation
            const double repR = std::sqrt(re * radius);
            perf_rep_radius_.push_back(repR);
            perf_length_.push_back(perf_length);
            bore_diameters_.push_back(2. * radius);
            num_active_connections++;
        }
        ++my_next_perf;
    }
    assert(my_next_perf == perf_data_->end());
    assert(num_active_connections == nperf);
}

void WellInterfaceGeneric::setWsolvent(const double wsolvent)
{
    wsolvent_ = wsolvent;
}

void WellInterfaceGeneric::setDynamicThpLimit(const double thp_limit)
{
    dynamic_thp_limit_ = thp_limit;
}

void WellInterfaceGeneric::updatePerforatedCell(std::vector<bool>& is_cell_perforated)
{

    for (int perf_idx = 0; perf_idx<number_of_perforations_; ++perf_idx) {
        is_cell_perforated[well_cells_[perf_idx]] = true;
    }
}

bool WellInterfaceGeneric::isVFPActive(DeferredLogger& deferred_logger) const
{
    // since the well_controls only handles the VFP number when THP constraint/target is there.
    // we need to get the table number through the parser, in case THP constraint/target is not there.
    // When THP control/limit is not active, if available VFP table is provided, we will still need to
    // update THP value. However, it will only used for output purpose.
    if (isProducer()) { // producer
        const int table_id = well_ecl_.vfp_table_number();
        if (table_id <= 0) {
            return false;
        } else {
            if (vfp_properties_->getProd()->hasTable(table_id)) {
                return true;
            } else {
                OPM_DEFLOG_THROW(std::runtime_error, "VFPPROD table " << std::to_string(table_id) << " is specified,"
                              << " for well " << name() << ", while we could not access it during simulation", deferred_logger);
            }
        }

    } else { // injector
        const int table_id = well_ecl_.vfp_table_number();
        if (table_id <= 0) {
            return false;
        } else {
            if (vfp_properties_->getInj()->hasTable(table_id)) {
                return true;
            } else {
                OPM_DEFLOG_THROW(std::runtime_error, "VFPINJ table " << std::to_string(table_id) << " is specified,"
                              << " for well " << name() << ", while we could not access it during simulation", deferred_logger);
            }
        }
    }
}

void WellInterfaceGeneric::updateWellTestStatePhysical(const double simulation_time,
                                                       const bool write_message_to_opmlog,
                                                       WellTestState& well_test_state,
                                                       DeferredLogger& deferred_logger) const
{
    if (!isOperableAndSolvable()) {
        if (well_test_state.well_is_closed(name())) {
            // Already closed, do nothing.
        } else {
            well_test_state.close_well(name(), WellTestConfig::Reason::PHYSICAL, simulation_time);
            if (write_message_to_opmlog) {
                const std::string action = well_ecl_.getAutomaticShutIn() ? "shut" : "stopped";
                const std::string msg = "Well " + name()
                    + " will be " + action + " as it can not operate under current reservoir conditions.";
                deferred_logger.info(msg);
            }
        }
    }
}

bool WellInterfaceGeneric::isOperableAndSolvable() const
{
    return operability_status_.isOperableAndSolvable();
}

bool WellInterfaceGeneric::useVfpExplicit() const
{
    const auto& wvfpexp = well_ecl_.getWVFPEXP();
    return ((wvfpexp.explicit_lookup() && !changedToOpenThisStep())|| operability_status_.use_vfpexplicit);
}

double WellInterfaceGeneric::getALQ(const WellState& well_state) const
{
    return well_state.getALQ(name());
}

void WellInterfaceGeneric::reportWellSwitching(const SingleWellState& ws, DeferredLogger& deferred_logger) const
{
    if (well_control_log_.empty())
        return;

    std::string msg = "    Well " + name()
        + " control mode changed from ";
    for (const std::string& from : well_control_log_) {
        msg += from + "->";
    }
    std::string to;
    if (isInjector()) {
        to = Well::InjectorCMode2String(ws.injection_cmode);
    } else {
        to = Well::ProducerCMode2String(ws.production_cmode);
    }
    msg += to;
    deferred_logger.info(msg);
}

std::optional<double>
WellInterfaceGeneric::
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
    if constexpr (extraBhpAtThpLimitProdOutput) {
        deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + this->name() +
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
                                    "Robust bhp(thp) solve failed due to inoperability for well " + this->name());
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
                                    "Bisect did not find the bhp-point where we produce for well " + this->name());
            return std::nullopt;
        }
    }
    if constexpr (extraBhpAtThpLimitProdOutput) {
        deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + this->name() +
                              "  low = " + std::to_string(low) +
                              "  high = " + std::to_string(high) +
                              "  f(low) = " + std::to_string(f_low) +
                              "  f(high) = " + std::to_string(f_high) +
                              "  bhp_max = " + std::to_string(bhp_max));
    }
    return bhp_max;
}


bool
WellInterfaceGeneric::
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
    if constexpr (extraBhpAtThpLimitProdOutput) {
        deferred_logger.debug("computeBhpAtThpLimitProd(): well = " + this->name() +
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
                                        "Robust bhp(thp) not solved precisely for well " + this->name());
                approximate_solution = abs_low < abs_high ? low : high;
            } else {
                    // Return failure.
                deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_BRACKETING_FAILURE",
                                         "Robust bhp(thp) solve failed due to bracketing failure for well " +
                                         this->name());
            }
        }
    } else {
        finding_bracket = true;
    }
    return finding_bracket;
}

bool
WellInterfaceGeneric::
bruteForceBracket(const std::function<double(const double)>& eq,
                  const std::array<double, 2>& range,
                  double& low, double& high,
                  DeferredLogger& deferred_logger) const
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

std::optional<double>
WellInterfaceGeneric::
computeBhpAtThpLimitProdCommon(const std::function<std::vector<double>(const double)>& frates,
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
    const auto& controls = this->wellEcl().productionControls(summary_state);
    const auto& table = this->vfpProperties()->getProd()->getTable(controls.vfp_table_number);
    const double vfp_ref_depth = table.getDatumDepth();
    const double thp_limit = this->getTHPConstraint(summary_state);
    const double dp = wellhelpers::computeHydrostaticCorrection(this->refDepth(), vfp_ref_depth, rho, this->gravity());

    auto fbhp = [this, &controls, thp_limit, dp, alq_value](const std::vector<double>& rates) {
        assert(rates.size() == 3);
        const auto& wfr =  this->vfpProperties()->getExplicitWFR(controls.vfp_table_number, this->indexOfWell());
        const auto& gfr = this->vfpProperties()->getExplicitGFR(controls.vfp_table_number, this->indexOfWell());
        const bool use_vfpexp = this->useVfpExplicit();
        return this->vfpProperties()->getProd()
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
                                "find bhp-point where production becomes non-zero for well " + this->name());
        return std::nullopt;
    }
    const std::array<double, 2> range {controls.bhp_limit, *bhp_max};
    return computeBhpAtThpLimitCommon(frates, fbhp, range, deferred_logger);
}

std::optional<double>
WellInterfaceGeneric::
computeBhpAtThpLimitCommon(const std::function<std::vector<double>(const double)>& frates,
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
                                "bracket the bhp solution with the brute force search for " + this->name());
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
                                "Robust bhp(thp) solve failed for well " + this->name());
        return std::nullopt;
    }
}
} // namespace Opm
