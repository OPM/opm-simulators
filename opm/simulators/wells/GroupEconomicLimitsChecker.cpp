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
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <fmt/format.h>

namespace Opm {

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

}

void
GroupEconomicLimitsChecker::
closeWells()
{
    displayDebugMessage("closing wells..");
    closeWellsRecursive(this->group_);
}

void
GroupEconomicLimitsChecker::
doWorkOver()
{
    displayDebugMessage("do work over..");
}

bool
GroupEconomicLimitsChecker::
endRun()
{
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
displayDebugMessage(const std::string &msg)
{
    if (this->debug_) {
        const std::string msg2 = fmt::format(
            "GECON: group: {} : {}", this->group_.name(), msg);
        this->deferred_logger_.debug(msg2);
    }
}

void
GroupEconomicLimitsChecker::
closeWellsRecursive(Group group)
{
    if (this->debug_) {
        const std::string msg = fmt::format("closing wells recursive : group {} ", group.name());
        displayDebugMessage(msg);
    }
    for (const std::string& group_name : group.groups()) {
        auto next_group = this->schedule_.getGroup(group_name, this->report_step_idx_);
        closeWellsRecursive(next_group);
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
            if (this->debug_) {
                const std::string msg = fmt::format(
                    "closing well {}", well_name);
                displayDebugMessage(msg);
            }
            this->well_test_state_.close_well(
                well_name, WellTestConfig::Reason::ECONOMIC, this->simulation_time_);
            this->well_model_.updateClosedWellsThisStep(well_name);
        }
    }
}

} // namespace Opm
