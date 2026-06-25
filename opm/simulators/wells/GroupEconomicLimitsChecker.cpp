/*
  Copyright 2023 Equinor ASA.

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
#include <opm/simulators/wells/GroupEconomicLimitsChecker.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <chrono>
#include <ctime>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace Opm {

std::string simTimeToString(const std::time_t start_time, const double sim_time)
{
    const auto start_timep = std::chrono::system_clock::from_time_t(start_time);
    const auto sim_duration = std::chrono::duration_cast<std::chrono::system_clock::duration>(
        std::chrono::duration<double>(sim_time)
    );
    const std::time_t cur_time = std::chrono::system_clock::to_time_t(start_timep + sim_duration);
    std::ostringstream ss;
    ss << std::put_time(std::localtime(&cur_time), "%d-%b-%Y");
    return ss.str();
}

template<typename Scalar, typename IndexTraits>
GroupEconomicLimitsChecker<Scalar, IndexTraits>::
GroupEconomicLimitsChecker(const BlackoilWellModelGeneric<Scalar, IndexTraits>& well_model,
                           WellTestState& well_test_state,
                           const Group& group,
                           const double simulation_time,
                           const int report_step_idx,
                           DeferredLogger& deferred_logger)
    : well_model_{well_model}
    , group_{group}
    , simulation_time_{simulation_time}
    , report_step_idx_{report_step_idx}
    , deferred_logger_{deferred_logger}
    , date_string_{simTimeToString(well_model.schedule().getStartTime(),simulation_time)}
    , unit_system_{well_model.eclipseState().getUnits()}
    , well_state_{well_model.wellState()}
    , well_test_state_{well_test_state}
    , schedule_{well_model.schedule()}
    , gecon_props_{schedule_[report_step_idx_].gecon().get_group_prop(
                   schedule_, well_model_.summaryState(), group_.name())}
{
    for (std::size_t i = 0; i < this->phase_idx_map_.size(); i++) {
        auto phase_idx = this->phase_idx_map_[i];
        this->phase_idx_reverse_map_[phase_idx] = static_cast<int>(i);
        auto phase_pos = this->well_model_.phaseUsage().canonicalToActivePhaseIdx(phase_idx);
        Scalar production_rate = this->well_model_.groupStateHelper().sumWellSurfaceRates(
            this->group_, phase_pos, /*isInjector*/false
        );
        this->production_rates_[i] = this->well_model_.comm().sum(production_rate);
    }
}

/****************************************
 * Public methods in alphabetical order
 ****************************************/

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
activateEndRun()
{
    displayDebugMessage("activate end run");
}

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
closeWells()
{
    closeWellsRecursive(this->group_);
}

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
doWorkOver(const RatioDetails& ratio_details)
{
    switch (this->gecon_props_.workover()) {
    case GroupEconProductionLimits::EconWorkover::NONE:
        break;
    case GroupEconProductionLimits::EconWorkover::WELL:
        this->closeWorstOffendingRatioWell(ratio_details);
        break;
    default:
        const auto workover_name = GroupEconProductionLimits::econWorkoverToString(this->gecon_props_.workover());
        throwNotImplementedError(fmt::format(
            "workover procedure {}",
            workover_name));
    }
}

template<typename Scalar, typename IndexTraits>
bool GroupEconomicLimitsChecker<Scalar, IndexTraits>::
endRun()
{
    if (this->gecon_props_.endRun()) {
        throwNotImplementedError("end run flag YES");
    }
    return false;
}

template<typename Scalar, typename IndexTraits>
std::optional<typename GroupEconomicLimitsChecker<Scalar, IndexTraits>::RatioDetails>
GroupEconomicLimitsChecker<Scalar, IndexTraits>::
GOR()
{
    const auto details = this->groupRatioDetails(RatioViolation::GOR);
    if (details.has_value()) {
        if (details->ratio > details->limit) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "GOR={} is greater than maximum: {}",
                    details->ratio, details->limit);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Gas/oil ratio = {:.2f} {} which is greater than the maximum economic value = {:.2f} {}",
                                details->ratio, details->limit, details->measure);
            return details;
        }
    }
    return std::nullopt;
}

template<typename Scalar, typename IndexTraits>
bool GroupEconomicLimitsChecker<Scalar, IndexTraits>::
minGasRate()
{
    auto phase_idx = this->phase_idx_reverse_map_[IndexTraits::gasPhaseIdx];
    auto gas_production_rate = this->production_rates_[phase_idx];
    if (this->debug_) {
        const std::string msg = fmt::format(
            "gecon: group: {}, gas_rate={}", this->group_.name(), gas_production_rate);
        displayDebugMessage(msg);
    }
    if (auto min_gas_rate = this->gecon_props_.minGasRate(); min_gas_rate) {
        if (gas_production_rate < *min_gas_rate) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "gas_rate={} is less than minimum: {}",
                    gas_production_rate, *min_gas_rate);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Gas rate = {:.2f} {} which is lower than the minimum economic value = {:.2f} {}",
                                gas_production_rate, *min_gas_rate, UnitSystem::measure::gas_surface_rate);
            return true;
        }
    }
    return false;
}

template<typename Scalar, typename IndexTraits>
bool GroupEconomicLimitsChecker<Scalar, IndexTraits>::
minOilRate()
{
    auto phase_idx = this->phase_idx_reverse_map_[IndexTraits::oilPhaseIdx];
    auto oil_production_rate = this->production_rates_[phase_idx];
    if (this->debug_) {
        const std::string msg = fmt::format(
            "oil_rate={}", oil_production_rate);
        displayDebugMessage(msg);
    }
    if (auto min_oil_rate = this->gecon_props_.minOilRate(); min_oil_rate) {
        if (oil_production_rate < *min_oil_rate) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "oil_rate={} is less than minimum: {}",
                    oil_production_rate, *min_oil_rate);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Oil rate = {:.2f} {} which is lower than the minimum economic value = {:.2f} {}",
                                oil_production_rate, *min_oil_rate, UnitSystem::measure::liquid_surface_rate);
            return true;
        }
    }
    return false;
}

template<typename Scalar, typename IndexTraits>
int GroupEconomicLimitsChecker<Scalar, IndexTraits>::
numProducersOpen()
{
    return 1;
}

template<typename Scalar, typename IndexTraits>
int GroupEconomicLimitsChecker<Scalar, IndexTraits>::
numProducersOpenInitially()
{
    return 1;
}

template <typename Scalar, typename IndexTraits>
std::optional<typename GroupEconomicLimitsChecker<Scalar, IndexTraits>::RatioDetails>
GroupEconomicLimitsChecker<Scalar, IndexTraits>::
ratioViolation()
{
    if (auto water_cut_details = this->waterCut(); water_cut_details.has_value())
    {
        return water_cut_details;
    }
    else if (auto gor_details = this->GOR(); gor_details.has_value())
    {
        return gor_details;
    }
    else if (auto wgr_details = this->WGR(); wgr_details.has_value())
    {
        return wgr_details;
    }

    return std::nullopt;
}

template<typename Scalar, typename IndexTraits>
std::optional<typename GroupEconomicLimitsChecker<Scalar, IndexTraits>::RatioDetails>
GroupEconomicLimitsChecker<Scalar, IndexTraits>::
waterCut()
{
    const auto details = this->groupRatioDetails(RatioViolation::WATER_CUT);
    if (details.has_value()) {
        if (details->ratio > details->limit) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "water_cut={} is greater than maximum: {}",
                    details->ratio, details->limit);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Water cut = {:.2f} {} which is greater than the maximum economic value = {:.2f} {}",
                                details->ratio, details->limit, details->measure);
            return details;
        }
    }
    return std::nullopt;
}

template<typename Scalar, typename IndexTraits>
std::optional<typename GroupEconomicLimitsChecker<Scalar, IndexTraits>::RatioDetails>
GroupEconomicLimitsChecker<Scalar, IndexTraits>::
WGR()
{
    const auto details = this->groupRatioDetails(RatioViolation::WGR);
    if (details.has_value()) {
        if (details->ratio > details->limit) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "WGR={} is greater than maximum: {}",
                    details->ratio, details->limit);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Water/gas ratio = {:.2f} {} which is greater than the maximum economic value = {:.2f} {}",
                                details->ratio, details->limit, details->measure);
            return details;
        }
    }
    return std::nullopt;
}

/****************************************
 * Private methods in alphabetical order
 ****************************************/

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
displayDebugMessage(const std::string& msg) const
{
    if (this->debug_) {
        const std::string msg2 = fmt::format(
            "GECON: group: {} : {}", this->group_.name(), msg);
        this->deferred_logger_.debug(msg2);
    }
}

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
addPrintMessage(const std::string& msg,
                const Scalar value,
                const Scalar limit,
                const UnitSystem::measure measure)
{
    const std::string header = fmt::format(
        "{}\nAt time = {:.2f} {} (date = {}): Group {} will close because: \n", this->message_separator(),
        this->unit_system_.from_si(UnitSystem::measure::time, this->simulation_time_),
        this->unit_system_.name(UnitSystem::measure::time),
        this->date_string_,
        this->group_.name()
    );
    const std::string measure_name(this->unit_system_.name(measure));
    const std::string message = fmt::format(fmt::runtime(msg),
                                         this->unit_system_.from_si(measure, value), measure_name,
                                         this->unit_system_.from_si(measure, limit), measure_name);

    this->message_ = header;
    this->message_ += message;
}

template<typename Scalar, typename IndexTraits>
bool GroupEconomicLimitsChecker<Scalar, IndexTraits>::
closeWellsRecursive(const Group& group, int level)
{
    bool wells_closed = false;

    if (this->debug_) {
        const std::string msg = fmt::format("closing wells recursive : group {} ", group.name());
        displayDebugMessage(msg);
    }
    for (const std::string& group_name : group.groups()) {
        auto next_group = this->schedule_.getGroup(group_name, this->report_step_idx_);
        wells_closed = wells_closed | closeWellsRecursive(next_group, level+1);
    }
    const auto indent = std::string(2*(level+1), ' ');
    if (level > 0) {
        const std::string msg = fmt::format("\n{}Closing group {}.", indent, group.name());
        this->message_ += msg;
    }

    if (this->debug_) {
        const std::string msg = fmt::format("closing wells recursive : group {} has {} wells",
                          group.name(), group.wells().size());
        displayDebugMessage(msg);
    }

    for (const std::string& well_name : group.wells()) {
        if (this->well_test_state_.well_is_closed(well_name)) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "well {} is already closed", well_name);
                displayDebugMessage(msg);
            }
        }
        else {
            wells_closed = true;
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "closing well {}", well_name);
                displayDebugMessage(msg);
            }
            const bool will_shut =
                this->schedule_.getWell(well_name, this->report_step_idx_).getAutomaticShutIn();
            const std::string msg = fmt::format("\n{} {} well {}", indent,
                                                will_shut ? "Shutting" : "Stopping",
                                                well_name);
            this->message_ += msg;

            // Only update the well_test_state_ on ranks that have the well in question.
            if (well_model_.hasLocalWell(well_name)) {
                this->well_test_state_.close_well(
                    well_name, WellTestConfig::Reason::GROUP, this->simulation_time_);
                this->well_model_.updateClosedWellsThisStep(well_name);
            }
        }
    }

    // If any wells were closed, output message at top level (group that hit constraint), on rank 0
    if (level == 0 && wells_closed && this->well_model_.comm().rank()==0) {
        this->message_ += ("\n" + this->message_separator() + "\n");
        this->deferred_logger_.info(this->message_);
    }
    return wells_closed;
}

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
throwNotImplementedError(const std::string& error) const
{
    const std::string msg = fmt::format("Group: {} : GECON : {} not implemented",
                                        this->group_.name(), error);
    OPM_DEFLOG_THROW(std::runtime_error, msg, this->deferred_logger_);
}

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
collectProducerWells(const Group& group, std::vector<std::string>& well_names) const
{
    for (const std::string& group_name : group.groups()) {
        const auto& sub_group = this->schedule_.getGroup(group_name, this->report_step_idx_);
        this->collectProducerWells(sub_group, well_names);
    }
    for (const std::string& well_name : group.wells()) {
        const auto& well_ecl = this->schedule_.getWell(well_name, this->report_step_idx_);
        if (well_ecl.isProducer()) {
            well_names.push_back(well_name);
        }
    }
}

template<typename Scalar, typename IndexTraits>
std::optional<Scalar> GroupEconomicLimitsChecker<Scalar, IndexTraits>::
computeWellRatio(const std::string& well_name,
                 const RatioViolation ratio_violation) const
{
    // The checks below avoid accessing WellState data for wells that are already closed,
    // not present on this rank, or not owned by this rank (important for parallel runs).
    if (this->well_test_state_.well_is_closed(well_name)) {
        return std::nullopt;
    }

    const auto well_index = this->well_state_.index(well_name);
    if (!well_index.has_value()) {
        return std::nullopt;
    }
    if (!this->well_state_.wellIsOwned(well_index.value(), well_name)) {
        return std::nullopt;
    }

    const auto& ws = this->well_state_.well(well_index.value());
    if (ws.status != Well::Status::OPEN) {
        return std::nullopt;
    }

    const auto& pu = this->well_model_.phaseUsage();
    constexpr Scalar big_value = std::numeric_limits<Scalar>::max();

    switch (ratio_violation) {
    case RatioViolation::WATER_CUT: {
        if (!pu.phaseIsActive(IndexTraits::oilPhaseIdx) ||
            !pu.phaseIsActive(IndexTraits::waterPhaseIdx))
        {
            return std::nullopt;
        }
        const auto oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        const auto water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        // Surface rates are negative for producers in WellState; flip sign here.
        const Scalar oil_rate = -ws.surface_rates[oil_pos];
        const Scalar water_rate = -ws.surface_rates[water_pos];
        const Scalar liquid = oil_rate + water_rate;
        if (liquid <= 0.0)      return Scalar{0};
        if (water_rate <= 0.0)  return Scalar{0};
        if (oil_rate   <= 0.0)  return Scalar{1};
        return water_rate / liquid;
    }
    case RatioViolation::GOR: {
        if (!pu.phaseIsActive(IndexTraits::oilPhaseIdx) ||
            !pu.phaseIsActive(IndexTraits::gasPhaseIdx))
        {
            return std::nullopt;
        }
        const auto oil_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        const auto gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        // Surface rates are negative for producers in WellState; flip sign here.
        const Scalar oil_rate = -ws.surface_rates[oil_pos];
        const Scalar gas_rate = -ws.surface_rates[gas_pos];
        if (gas_rate <= 0.0)    return Scalar{0};
        if (oil_rate <= 0.0)    return big_value;
        return gas_rate / oil_rate;
    }
    case RatioViolation::WGR: {
        if (!pu.phaseIsActive(IndexTraits::gasPhaseIdx) ||
            !pu.phaseIsActive(IndexTraits::waterPhaseIdx))
        {
            return std::nullopt;
        }
        const auto gas_pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        const auto water_pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        // Surface rates are negative for producers in WellState; flip sign here.
        const Scalar gas_rate = -ws.surface_rates[gas_pos];
        const Scalar water_rate = -ws.surface_rates[water_pos];
        if (water_rate <= 0.0)  return Scalar{0};
        if (gas_rate   <= 0.0)  return big_value;
        return water_rate / gas_rate;
    }
    case RatioViolation::NONE:
    default:
        return std::nullopt;
    }
}

template<typename Scalar, typename IndexTraits>
std::optional<typename GroupEconomicLimitsChecker<Scalar, IndexTraits>::RatioDetails>
GroupEconomicLimitsChecker<Scalar, IndexTraits>::
groupRatioDetails(const RatioViolation ratio_violation) const
{
    constexpr Scalar big_value = std::numeric_limits<Scalar>::max();
    switch (ratio_violation) {
    case RatioViolation::WGR: {
        const auto max_wgr = this->gecon_props_.maxWaterGasRatio();
        if (!max_wgr.has_value()) {
            return std::nullopt;
        }
        const auto water_phase_idx = this->phase_idx_reverse_map_.at(IndexTraits::waterPhaseIdx);
        const auto gas_phase_idx = this->phase_idx_reverse_map_.at(IndexTraits::gasPhaseIdx);
        const Scalar water_rate = this->production_rates_[water_phase_idx];
        const Scalar gas_rate = this->production_rates_[gas_phase_idx];
        const Scalar ratio = (water_rate <= 0.0) ? Scalar{0}
                           : (gas_rate <= 0.0)   ? Scalar{big_value}
                                                  : water_rate / gas_rate;
        // No dedicated water-gas-ratio measure exists; a water/gas ratio has the same
        // (liquid volume / gas volume) units as the oil-gas ratio.
        return RatioDetails{RatioViolation::WGR, "Water-gas ratio", ratio, static_cast<Scalar>(max_wgr.value()), UnitSystem::measure::oil_gas_ratio};
    }
    case RatioViolation::GOR: {
        const auto max_gor = this->gecon_props_.maxGasOilRatio();
        if (!max_gor.has_value()) {
            return std::nullopt;
        }
        const auto oil_phase_idx = this->phase_idx_reverse_map_.at(IndexTraits::oilPhaseIdx);
        const auto gas_phase_idx = this->phase_idx_reverse_map_.at(IndexTraits::gasPhaseIdx);
        const Scalar oil_rate = this->production_rates_[oil_phase_idx];
        const Scalar gas_rate = this->production_rates_[gas_phase_idx];
        const Scalar ratio = (gas_rate <= 0.0) ? Scalar{0}
                           : (oil_rate <= 0.0) ? Scalar{big_value}
                                               : gas_rate / oil_rate;
        return RatioDetails{RatioViolation::GOR, "Gas-oil ratio", ratio, static_cast<Scalar>(max_gor.value()), UnitSystem::measure::gas_oil_ratio};
    }
    case RatioViolation::WATER_CUT: {
        const auto max_water_cut = this->gecon_props_.maxWaterCut();
        if (!max_water_cut.has_value()) {
            return std::nullopt;
        }
        const auto oil_phase_idx = this->phase_idx_reverse_map_.at(IndexTraits::oilPhaseIdx);
        const auto water_phase_idx = this->phase_idx_reverse_map_.at(IndexTraits::waterPhaseIdx);
        const Scalar oil_rate = this->production_rates_[oil_phase_idx];
        const Scalar water_rate = this->production_rates_[water_phase_idx];
        const Scalar liquid_rate = oil_rate + water_rate;
        const Scalar ratio = (liquid_rate == 0.0) ? Scalar{0}
                           : (water_rate < 0.0)   ? Scalar{0}
                           : (oil_rate < 0.0)     ? Scalar{1}
                                                   : water_rate / liquid_rate;
        return RatioDetails{RatioViolation::WATER_CUT, "Water cut", ratio, static_cast<Scalar>(max_water_cut.value()), UnitSystem::measure::water_cut};
    }
    case RatioViolation::NONE:
    default:
        return std::nullopt;
    }
}

template<typename Scalar, typename IndexTraits>
void GroupEconomicLimitsChecker<Scalar, IndexTraits>::
closeWorstOffendingRatioWell(const RatioDetails& ratio_details)
{
    if (ratio_details.violation == RatioViolation::NONE) {
        // Should not happen: doWorkOver() is only called after a ratio
        // violation was reported, but guard against misuse.
        return;
    }

    std::vector<std::string> producer_wells;
    this->collectProducerWells(this->group_, producer_wells);

    // Each rank reports its locally-owned ratio, or a negative sentinel where the
    // ratio is not applicable. Ratios are non-negative, so the sentinel never wins
    // the MPI max-reduction; the worst ratio is then picked from the reduced array.
    constexpr Scalar not_applicable = std::numeric_limits<Scalar>::lowest();
    std::vector<Scalar> global_ratios;
    global_ratios.reserve(producer_wells.size());
    for (const std::string& well_name : producer_wells) {
        const auto ratio = this->computeWellRatio(well_name, ratio_details.violation);
        global_ratios.push_back(ratio.value_or(not_applicable));
    }

    std::string worst_well;
    Scalar worst_ratio = not_applicable;
    if (!global_ratios.empty()) {
        this->well_model_.comm().max(global_ratios.data(), global_ratios.size());
        for (std::size_t i = 0; i < global_ratios.size(); ++i) {
            if (global_ratios[i] > worst_ratio) {
                worst_ratio = global_ratios[i];
                worst_well = producer_wells[i];
            }
        }
    }

    if (worst_well.empty()) {
        // No candidate producer found (e.g. all already closed).
        if (this->debug_) {
            displayDebugMessage(
                "GECON WELL workover: no open producer well found to close.");
        }
        return;
    }


    // Perform the actual close on the owning rank(s).
    bool well_was_closed_locally = false;
    if (this->well_model_.hasLocalWell(worst_well)) {
        if (!this->well_test_state_.well_is_closed(worst_well)) {
            this->well_test_state_.close_well(
                worst_well, WellTestConfig::Reason::GROUP, this->simulation_time_);
            this->well_model_.updateClosedWellsThisStep(worst_well);
            well_was_closed_locally = true;
        }
    }

    if (this->well_model_.comm().max(static_cast<int>(well_was_closed_locally)) == 0) {
        // Nothing to do (well was already closed everywhere).
        return;
    }

    if (this->well_model_.comm().rank() == 0) {
        const Scalar group_ratio_display = this->unit_system_.from_si(ratio_details.measure, ratio_details.ratio);
        const Scalar limit_display = this->unit_system_.from_si(ratio_details.measure, ratio_details.limit);
        const std::string unit_name = this->unit_system_.name(ratio_details.measure);
        const bool will_shut =
            this->schedule_.getWell(worst_well, this->report_step_idx_).getAutomaticShutIn();

        const std::string msg = fmt::format(
            "{}\nAt time = {:.2f} {} (date = {}): Well {} will be {} because:\n"
            "  {} for group {} = {:.4e} {} exceeds the limit {:.4e} {}.\n{}",
            this->message_separator(),
            this->unit_system_.from_si(UnitSystem::measure::time, this->simulation_time_),
            this->unit_system_.name(UnitSystem::measure::time),
            this->date_string_,
            worst_well,
            will_shut ? "shut" : "stopped",
            ratio_details.ratio_type,
            this->group_.name(),
            group_ratio_display, unit_name,
            limit_display, unit_name,
            this->message_separator());
        this->deferred_logger_.info(msg);
    }
}

template class GroupEconomicLimitsChecker<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class GroupEconomicLimitsChecker<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
