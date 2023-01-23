/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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
#include <opm/simulators/wells/BlackoilWellModelRestart.hpp>

#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/output/data/Groups.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>

namespace {
    Opm::data::GuideRateValue::Item
    guideRateRestartItem(const Opm::GuideRateModel::Target target)
    {
        using Item = Opm::data::GuideRateValue::Item;
        using Target = Opm::GuideRateModel::Target;

        static const auto items = std::unordered_map<Target, Item> {
            { Target::OIL, Item::Oil   },
            { Target::GAS, Item::Gas   },
            { Target::WAT, Item::Water },
            { Target::RES, Item::ResV  },
        };

        auto i = items.find(target);
        return (i == items.end()) ? Item::NumItems : i->second;
    }

    Opm::GuideRate::GuideRateValue
    makeGuideRateValue(const Opm::data::GuideRateValue&  restart,
                       const Opm::GuideRateModel::Target target)
    {
        const auto item = guideRateRestartItem(target);

        if (! restart.has(item)) {
            return {};
        }

        return { 0.0, restart.get(item), target };
    }
} // Anonymous

namespace Opm {

void BlackoilWellModelRestart::
loadRestartConnectionData(const std::vector<data::Rates::opt>& phs,
                          const data::Well&                    rst_well,
                          const std::vector<PerforationData>&  old_perf_data,
                          SingleWellState&                     ws) const
{
    auto& perf_data        = ws.perf_data;
    auto  perf_pressure    = perf_data.pressure.begin();
    auto  perf_rates       = perf_data.rates.begin();
    auto  perf_phase_rates = perf_data.phase_rates.begin();

    for (const auto& pd : old_perf_data) {
        const auto& rst_connection = rst_well.connections[pd.ecl_index];

        *perf_pressure = rst_connection.pressure;       ++perf_pressure;
        *perf_rates    = rst_connection.reservoir_rate; ++perf_rates;

        for (const auto& phase : phs) {
            *perf_phase_rates = rst_connection.rates.get(phase);
            ++perf_phase_rates;
        }
    }
}

void BlackoilWellModelRestart::
loadRestartSegmentData(const std::string&                   well_name,
                       const std::vector<data::Rates::opt>& phs,
                       const data::Well&                    rst_well,
                       SingleWellState&                     ws) const
{
    const auto& segment_set = wellModel_.getWellEcl(well_name).getSegments();
    const auto& rst_segments = rst_well.segments;

    // \Note: Eventually we need to handle the situations that some segments are shut
    assert(0u + segment_set.size() == rst_segments.size());

    const auto np = phs.size();
    const auto pres_idx = data::SegmentPressures::Value::Pressure;

    auto& segments = ws.segments;
    auto& segment_pressure = segments.pressure;
    auto& segment_rates = segments.rates;
    for (const auto& [segNum, rst_segment] : rst_segments) {
        const int segment_index = segment_set.segmentNumberToIndex(segNum);

        // Recovering segment rates and pressure from the restart values
        segment_pressure[segment_index] = rst_segment.pressures[pres_idx];

        const auto& rst_segment_rates = rst_segment.rates;
        for (auto p = 0*np; p < np; ++p) {
            segment_rates[segment_index*np + p] = rst_segment_rates.get(phs[p]);
        }
    }
}

void BlackoilWellModelRestart::
loadRestartWellData(const std::string&                   well_name,
                    const bool                           handle_ms_well,
                    const std::vector<data::Rates::opt>& phs,
                    const data::Well&                    rst_well,
                    const std::vector<PerforationData>&  old_perf_data,
                    SingleWellState&                     ws) const
{
    const auto np = phs.size();

    ws.bhp = rst_well.bhp;
    ws.thp = rst_well.thp;
    ws.temperature = rst_well.temperature;

    if (rst_well.current_control.isProducer) {
        ws.production_cmode = rst_well.current_control.prod;
    }
    else {
        ws.injection_cmode = rst_well.current_control.inj;
    }

    for (auto i = 0*np; i < np; ++i) {
        assert( rst_well.rates.has( phs[ i ] ) );
        ws.surface_rates[i] = rst_well.rates.get(phs[i]);
    }

    this->loadRestartConnectionData(phs, rst_well, old_perf_data, ws);

    if (handle_ms_well && !rst_well.segments.empty()) {
        this->loadRestartSegmentData(well_name, phs, rst_well, ws);
    }
}

void BlackoilWellModelRestart::
loadRestartGroupData(const std::string&     group,
                     const data::GroupData& value,
                     GroupState& grpState) const
{
    using GPMode = Group::ProductionCMode;
    using GIMode = Group::InjectionCMode;

    const auto cpc = value.currentControl.currentProdConstraint;
    const auto cgi = value.currentControl.currentGasInjectionConstraint;
    const auto cwi = value.currentControl.currentWaterInjectionConstraint;

    if (cpc != GPMode::NONE) {
        grpState.production_control(group, cpc);
    }

    if (cgi != GIMode::NONE) {
        grpState.injection_control(group, Phase::GAS, cgi);
    }

    if (cwi != GIMode::NONE) {
        grpState.injection_control(group, Phase::WATER, cwi);
    }
}

void BlackoilWellModelRestart::
loadRestartGuideRates(const int                    report_step,
                      const GuideRateModel::Target target,
                      const data::Wells&           rst_wells,
                      GuideRate&                   guide_rate) const
{
    for (const auto& [well_name, rst_well] : rst_wells) {
        if (!wellModel_.hasWell(well_name) || wellModel_.getWellEcl(well_name).isInjector()) {
            continue;
        }

        guide_rate.init_grvalue_SI(report_step, well_name,
                                   makeGuideRateValue(rst_well.guide_rates, target));
    }
}

void BlackoilWellModelRestart::
loadRestartGuideRates(const int                                     report_step,
                      const GuideRateConfig&                        config,
                      const std::map<std::string, data::GroupData>& rst_groups,
                      GuideRate&                                    guide_rate) const
{
    const auto target = config.model().target();

    for (const auto& [group_name, rst_group] : rst_groups) {
        if (!config.has_production_group(group_name)) {
            continue;
        }

        const auto& group = config.production_group(group_name);
        if ((group.guide_rate > 0.0) || (group.target != Group::GuideRateProdTarget::FORM)) {
            continue;
        }

        guide_rate.init_grvalue_SI(report_step, group_name,
                                   makeGuideRateValue(rst_group.guideRates.production, target));
    }
}

void BlackoilWellModelRestart::
loadRestartData(const data::Wells&                 rst_wells,
                const data::GroupAndNetworkValues& grpNwrkValues,
                const bool                         handle_ms_well,
                WellState&                         well_state,
                GroupState&                        grpState) const
{
    using rt = data::Rates::opt;
    const auto& phases = wellModel_.phaseUsage();
    const auto np = phases.num_phases;

    std::vector<rt> phs(np);
    if (phases.phase_used[BlackoilPhases::Aqua]) {
        phs.at(phases.phase_pos[BlackoilPhases::Aqua]) = rt::wat;
    }

    if (phases.phase_used[BlackoilPhases::Liquid]) {
        phs.at( phases.phase_pos[BlackoilPhases::Liquid] ) = rt::oil;
    }

    if (phases.phase_used[BlackoilPhases::Vapour]) {
        phs.at( phases.phase_pos[BlackoilPhases::Vapour] ) = rt::gas;
    }

    for (auto well_index = 0*well_state.size();
         well_index < well_state.size();
         ++well_index)
    {
        const auto& well_name = well_state.name(well_index);

        this->loadRestartWellData(well_name, handle_ms_well, phs,
                                  rst_wells.at(well_name),
                                  wellModel_.perfData(well_index),
                                  well_state.well(well_index));
    }

    for (const auto& [group, value] : grpNwrkValues.groupData) {
        this->loadRestartGroupData(group, value, grpState);
    }
}

} // namespace Opm
