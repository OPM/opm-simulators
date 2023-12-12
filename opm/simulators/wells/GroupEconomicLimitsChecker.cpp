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
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>

#include <fmt/format.h>

#include <ctime>
#include <chrono>
#include <sstream>
#include <iomanip>

namespace Opm {

std::string simTimeToString(const std::time_t start_time, const double sim_time) {
    const auto start_timep = std::chrono::system_clock::from_time_t(start_time);
    const auto sim_duration = std::chrono::duration_cast<std::chrono::system_clock::duration>(
        std::chrono::duration<double>(sim_time)
    );
    const std::time_t cur_time = std::chrono::system_clock::to_time_t(start_timep + sim_duration);
    std::ostringstream ss;
    ss << std::put_time(std::localtime(&cur_time), "%d-%b-%Y");
    return ss.str();
}

GroupEconomicLimitsChecker::
GroupEconomicLimitsChecker(
    const BlackoilWellModelGeneric &well_model,
    WellTestState &well_test_state,
    const Group &group,
    const double simulation_time,
    const int report_step_idx,
    DeferredLogger &deferred_logger
) :
    well_model_{well_model}
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
        auto phase_pos = this->well_model_.phaseUsage().phase_pos[phase_idx];
        double production_rate = WellGroupHelpers::sumWellSurfaceRates(
            this->group_, this->schedule_, this->well_state_,
            this->report_step_idx_, phase_pos, /*isInjector*/false);
        this->production_rates_[i] = this->well_model_.comm().sum(production_rate);
    }
}

/****************************************
 * Public methods in alphabetical order
 ****************************************/

void
GroupEconomicLimitsChecker::
activateEndRun()
{
    displayDebugMessage("activate end run");
}

void
GroupEconomicLimitsChecker::
closeWells()
{
    closeWellsRecursive(this->group_);
}

void
GroupEconomicLimitsChecker::
doWorkOver()
{
    if (this->gecon_props_.workover() != GroupEconProductionLimits::EconWorkover::NONE) {
        throwNotImplementedError("workover procedure");
    }
}

bool
GroupEconomicLimitsChecker::
endRun()
{
    if (this->gecon_props_.endRun()) {
        throwNotImplementedError("end run flag YES");
    }
    return false;
}

bool
GroupEconomicLimitsChecker::
GOR()
{
    auto oil_phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Liquid];
    auto gas_phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Vapour];
    auto oil_rate = this->production_rates_[oil_phase_idx];
    auto gas_rate = this->production_rates_[gas_phase_idx];
    double gor;
    if (gas_rate <= 0.0) {
        gor = 0.0;
    }
    else if (oil_rate <= 0.0) {
        gor = 1e100;
    }
    else {
        gor = gas_rate / oil_rate;
    }
    if (auto max_gor = this->gecon_props_.maxGasOilRatio(); max_gor) {
        if (gor > *max_gor) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "GOR={} is greater than maximum: {}",
                    gor, *max_gor);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Gas/oil ratio = {:.2f} {} which is greater than the minimum economic value = {:.2f} {}",
                                gor, *max_gor, UnitSystem::measure::gas_oil_ratio);
            return true;
        }
    }
    return false;
}

bool
GroupEconomicLimitsChecker::
minGasRate()
{
    auto phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Vapour];
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

bool
GroupEconomicLimitsChecker::
minOilRate()
{
    auto phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Liquid];
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

int
GroupEconomicLimitsChecker::
numProducersOpen()
{
    return 1;
}

int
GroupEconomicLimitsChecker::
numProducersOpenInitially()
{
    return 1;
}

bool
GroupEconomicLimitsChecker::
waterCut()
{
    auto oil_phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Liquid];
    auto water_phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Aqua];
    auto oil_rate = this->production_rates_[oil_phase_idx];
    auto water_rate = this->production_rates_[water_phase_idx];
    auto liquid_rate = oil_rate + water_rate;
    double water_cut;
    if (liquid_rate == 0.0) {
        water_cut = 0.0;
    }
    else {
        if (water_rate < 0.0) {
            water_cut = 0.0;
        }
        else if (oil_rate < 0.0) {
            water_cut = 1.0;
        }
        else {
            water_cut = water_rate / liquid_rate;
        }
    }
    if (auto max_water_cut = this->gecon_props_.maxWaterCut(); max_water_cut) {
        if (water_cut > *max_water_cut) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "water_cut={} is greater than maximum: {}",
                    water_cut, *max_water_cut);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Water cut = {:.2f} {} which is greater than the maximum economic value = {:.2f} {}",
                                water_cut, *max_water_cut, UnitSystem::measure::water_cut);
            return true;
        }
    }
    return false;
}

bool
GroupEconomicLimitsChecker::
WGR()
{
    auto water_phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Aqua];
    auto gas_phase_idx = this->phase_idx_reverse_map_[BlackoilPhases::Vapour];
    auto water_rate = this->production_rates_[water_phase_idx];
    auto gas_rate = this->production_rates_[gas_phase_idx];
    double wgr;
    if (water_rate <= 0.0) {
        wgr = 0.0;
    }
    else if (gas_rate <= 0.0) {
        wgr = 1e100;
    }
    else {
        wgr = water_rate / gas_rate;
    }
    if (auto max_wgr = this->gecon_props_.maxWaterGasRatio(); max_wgr) {
        if (wgr > *max_wgr) {
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "WGR={} is greater than maximum: {}",
                    wgr, *max_wgr);
                displayDebugMessage(msg);
            }
            addPrintMessage("  Water/gas ratio = {:.2f} {} which is greater than the maximum economic value = {:.2f} {}",
                                wgr, *max_wgr, UnitSystem::measure::gas_oil_ratio); // Same units
            return true;
        }
    }
    return false;
}

/****************************************
 * Private methods in alphabetical order
 ****************************************/

void
GroupEconomicLimitsChecker::
displayDebugMessage(const std::string &msg) const
{
    if (this->debug_) {
        const std::string msg2 = fmt::format(
            "GECON: group: {} : {}", this->group_.name(), msg);
        this->deferred_logger_.debug(msg2);
    }
}

void
GroupEconomicLimitsChecker::
addPrintMessage(const std::string &msg, const double value, const double limit, const UnitSystem::measure measure)
{
    const std::string header = fmt::format(
        "{}\nAt time = {:.2f} {} (date = {}): Group {} will close because: \n", this->message_separator(),
        this->unit_system_.from_si(UnitSystem::measure::time, this->simulation_time_),
        this->unit_system_.name(UnitSystem::measure::time),
        this->date_string_,
        this->group_.name()
    );
    const std::string measure_name(this->unit_system_.name(measure));
    const std::string message = fmt::format(msg,
                                         this->unit_system_.from_si(measure, value), measure_name,
                                         this->unit_system_.from_si(measure, limit), measure_name);

    this->message_ = header;
    this->message_ += message;
}

bool
GroupEconomicLimitsChecker::
closeWellsRecursive(Group group, int level)
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
            const std::string msg = fmt::format("\n{} Closing well {}", indent, well_name);
            this->message_ += msg;

            this->well_test_state_.close_well(
                well_name, WellTestConfig::Reason::GROUP, this->simulation_time_);
            this->well_model_.updateClosedWellsThisStep(well_name);
        }
    }

    // If any wells were closed, output message at top level (group that hit constraint), on rank 0
    if (level == 0 && wells_closed && this->well_model_.comm().rank()==0) {
        this->message_ += ("\n" + this->message_separator() + "\n");
        this->deferred_logger_.info(this->message_);
    }
    return wells_closed;
}

void
GroupEconomicLimitsChecker::
throwNotImplementedError(const std::string &error) const
{
    const std::string msg = fmt::format("Group: {} : GECON : {} not implemented", this->group_.name(), error);
    OPM_DEFLOG_THROW(std::runtime_error, msg, this->deferred_logger_);
}
} // namespace Opm
