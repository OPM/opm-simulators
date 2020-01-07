/*
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_WELLANDGROUPSTATES_HEADER_INCLUDED
#define OPM_WELLANDGROUPSTATES_HEADER_INCLUDED

#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/SingleGroupState.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/simulators/wells/PerforationData.hpp>

#include <map>
#include <string>
#include <vector>

namespace Opm
{

/// Encapsulate all data needed to represent the state of all wells
/// and groups for persistence purposes, i.e. from one timestep to the
/// next or for restarting, and for accumulating rates to get correct
/// group rate limits.
/// In a parallel context, the well states will be those for local
/// wells only, while the group information will be global (and
/// communicated as needed).
template <int NumActivePhases>
class WellAndGroupStates
{
public:

    void init(const std::vector<double>& cell_pressures,
              const Schedule& schedule,
              const std::vector<Well>& wells,
              const int report_step,
              const WellAndGroupStates* prev_state,
              const PhaseUsage& phase_usage,
              const std::vector<std::vector<PerforationData>>& well_perf_data,
              const SummaryState& summary_state)
    {
        const int nw = wells.size();
        well_states_.resize(nw);
        well_names_.resize(nw);
        for (int w = 0; w < nw; ++w) {
            well_names_[w] = wells[w].name();
            initSingleWell(cell_pressures, schedule, wells[w], report_step, phase_usage, well_perf_data[w], summary_state, well_states_[w]);
        }
        static_cast<void>(prev_state);
        // TODO: deal with previous state.

        // TODO: initialize group state
    }

    data::Wells report(const PhaseUsage& phase_usage, const int* globalCellIdxMap) const
    {
        data::Wells dwells;
        const int nw = well_states_.size();
        for (int w = 0; w < nw; ++w) {
            if (well_states_[w].status == Well::Status::SHUT) {
                continue;
            }
            auto& dwell = dwells[well_names_[w]];
            reportSingleWell(phase_usage, globalCellIdxMap, well_states_[w], dwell);
        }
        return dwells;
    }

    using WellState = SingleWellState<NumActivePhases>;
    using GroupState = SingleGroupState<NumActivePhases>;

    std::vector<WellState>& wellStates()
    {
        return well_states_;
    }
    const std::vector<WellState>& wellStates() const
    {
        return well_states_;
    }

    std::map<std::string, GroupState>& groupStates()
    {
        return group_states_;
    }
    const std::map<std::string, GroupState>& groupStates() const
    {
        return group_states_;
    }


private:

    // -----------  Data members  -----------

    std::vector<WellState> well_states_;
    std::vector<std::string> well_names_;
    std::map<std::string, GroupState> group_states_;



    // -----------  Private functions  -----------

    static void initSingleWell(const std::vector<double>& cell_pressures,
                               const Schedule& schedule,
                               const Well& well,
                               const int report_step,
                               const PhaseUsage& phase_usage,
                               const std::vector<PerforationData>& perf_data,
                               const SummaryState& summary_state,
                               SingleWellState<NumActivePhases>& wstate)
    {
        assert(well.isInjector() || well.isProducer());

        wstate.status = well.getStatus();
        wstate.is_producer = well.isProducer();
        // At the moment, the following events are considered to be effective events
        // more events might join as effective events:
        // PRODUCTION_UPDATE, INJECTION_UPDATE, WELL_STATUS_CHANGE
        const uint64_t effective_events_mask =
            ScheduleEvents::WELL_STATUS_CHANGE
            + ScheduleEvents::PRODUCTION_UPDATE
            + ScheduleEvents::INJECTION_UPDATE;
        wstate.effective_events_occurred = schedule.hasWellGroupEvent(well.name(), effective_events_mask, report_step);

        const auto inj_controls = well.isInjector() ? well.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = well.isProducer() ? well.productionControls(summary_state) : Well::ProductionControls(0);

        if (well.isInjector()) {
            wstate.current_injection_control = inj_controls.cmode;
        } else {
            wstate.current_production_control = prod_controls.cmode;
        }

        const int num_perf_this_well = perf_data.size();
        wstate.connections.resize(num_perf_this_well);

        if ( num_perf_this_well == 0 ) {
            // No perforations of the well. Initialize pressures to zero.
            wstate.bhp = 0.0;
            wstate.thp = 0.0;
            return;
        }

        const bool is_bhp = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::BHP)
            : (prod_controls.cmode == Well::ProducerCMode::BHP);
        const double bhp_limit = well.isInjector() ? inj_controls.bhp_limit : prod_controls.bhp_limit;
        const bool is_grup = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::GRUP)
            : (prod_controls.cmode == Well::ProducerCMode::GRUP);

        const double inj_surf_rate = well.isInjector() ? inj_controls.surface_rate : 0.0; // To avoid a "maybe-uninitialized" warning.

        if (well.getStatus() == Well::Status::STOP) {
            // Stopped well:
            // 1. Rates: zero well rates.
            // 2. Bhp: assign bhp equal to bhp control, if
            //    applicable, otherwise assign equal to
            //    first perforation cell pressure.
            if (is_bhp) {
                wstate.bhp = bhp_limit;
            } else {
                const int first_cell = perf_data[0].cell_index;
                wstate.bhp = cell_pressures[first_cell];
            }
        } else if (is_grup) {
            // Well under group control.
            // 1. Rates: zero well rates.
            // 2. Bhp: initialize bhp to be a
            //    little above or below (depending on if
            //    the well is an injector or producer)
            //    pressure in first perforation cell.
            const int first_cell = perf_data[0].cell_index;
            const double safety_factor = well.isInjector() ? 1.01 : 0.99;
            wstate.bhp = safety_factor*cell_pressures[first_cell];
        } else {
            // Open well, under own control:
            // 1. Rates: initialize well rates to match
            //    controls if type is ORAT/GRAT/WRAT
            //    (producer) or RATE (injector).
            //    Otherwise, we cannot set the correct
            //    value here and initialize to zero rate.
            if (well.isInjector()) {
                if (inj_controls.cmode == Well::InjectorCMode::RATE) {
                    switch (inj_controls.injector_type) {
                    case Well::InjectorType::WATER:
                        assert(phase_usage.phase_used[BlackoilPhases::Aqua]);
                        wstate.surface_rates[phase_usage.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
                        break;
                    case Well::InjectorType::GAS:
                        assert(phase_usage.phase_used[BlackoilPhases::Vapour]);
                        wstate.surface_rates[phase_usage.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
                        break;
                    case Well::InjectorType::OIL:
                        assert(phase_usage.phase_used[BlackoilPhases::Liquid]);
                        wstate.surface_rates[phase_usage.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
                        break;
                    case Well::InjectorType::MULTI:
                        // Not currently handled, keep zero init.
                        break;
                    }
                } else {
                    // Keep zero init.
                }
            } else {
                assert(well.isProducer());
                // Note negative rates for producing wells.
                switch (prod_controls.cmode) {
                case Well::ProducerCMode::ORAT:
                    assert(phase_usage.phase_used[BlackoilPhases::Liquid]);
                    wstate.surface_rates[phase_usage.phase_pos[BlackoilPhases::Liquid]] = -prod_controls.oil_rate;
                    break;
                case Well::ProducerCMode::WRAT:
                    assert(phase_usage.phase_used[BlackoilPhases::Aqua]);
                    wstate.surface_rates[phase_usage.phase_pos[BlackoilPhases::Aqua]] = -prod_controls.water_rate;
                    break;
                case Well::ProducerCMode::GRAT:
                    assert(phase_usage.phase_used[BlackoilPhases::Vapour]);
                    wstate.surface_rates[phase_usage.phase_pos[BlackoilPhases::Vapour]] = -prod_controls.gas_rate;
                    break;
                default:
                    // Keep zero init.
                    break;
                }
            }
            // 2. Bhp: initialize bhp to be target pressure if
            //    bhp-controlled well, otherwise set to a
            //    little above or below (depending on if
            //    the well is an injector or producer)
            //    pressure in first perforation cell.
            if (is_bhp) {
                wstate.bhp = bhp_limit;
            } else {
                const int first_cell = perf_data[0].cell_index;
                const double safety_factor = well.isInjector() ? 1.01 : 0.99;
                wstate.bhp = safety_factor*cell_pressures[first_cell];
            }
        }

        // 3. Thp: assign thp equal to thp target/limit, if such a limit exists,
        //    otherwise keep it zero.
        const bool has_thp = well.isInjector() ? inj_controls.hasControl(Well::InjectorCMode::THP)
            : prod_controls.hasControl(Well::ProducerCMode::THP);
        const double thp_limit = well.isInjector() ? inj_controls.thp_limit : prod_controls.thp_limit;
        if (has_thp) {
            wstate.thp = thp_limit;
        }

        // Initialize multi-segment well parts.
        if (well.isMultiSegment()) {
            initMultiSegment(well, wstate);
        }
    } // initSingleWell()




    static void initMultiSegment(const Well& well, SingleWellState<NumActivePhases>& wstate)
    {
        const WellSegments& segment_set = well.getSegments();
        const int well_nseg = segment_set.size();
        wstate.segments.resize(well_nseg);

        // we need to know for each segment, how many perforation it has and how many segments using it as outlet_segment
        // that is why I think we should use a well model to initialize the WellState here
        std::vector<std::vector<int>> segment_perforations(well_nseg);
        std::vector<std::vector<int>> segment_inlets(well_nseg);
        {
            int n_activeperf = 0;
            const WellConnections& completion_set = well.getConnections();
            for (size_t perf = 0; perf < completion_set.size(); ++perf) {
                const Connection& connection = completion_set.get(perf);
                if (connection.state() == Connection::State::OPEN) {
                    const int segment_index = segment_set.segmentNumberToIndex(connection.segment());
                    segment_perforations[segment_index].push_back(n_activeperf);
                    n_activeperf++;
                }
            }

            for (int seg = 0; seg < well_nseg; ++seg) {
                const Segment& segment = segment_set[seg];
                const int segment_number = segment.segmentNumber();
                const int outlet_segment_number = segment.outletSegment();
                if (outlet_segment_number > 0) {
                    const int segment_index = segment_set.segmentNumberToIndex(segment_number);
                    const int outlet_segment_index = segment_set.segmentNumberToIndex(outlet_segment_number);
                    segment_inlets[outlet_segment_index].push_back(segment_index);
                }
            }
        }

        // Set segment numbers.
        for (int seg = 0; seg < well_nseg; ++seg) {
            wstate.segments[seg].segment_number = segment_set[seg].segmentNumber();
        }

        // Set segment pressures.
        wstate.segments[0].pressure = wstate.bhp;
        for (int seg = 1; seg < well_nseg; ++seg) {
            if (!segment_perforations[seg].empty()) {
                const int first_perf = segment_perforations[seg][0];
                wstate.segments[seg].pressure = wstate.connections[first_perf].pressure;
            } else {
                // using the outlet segment pressure // it needs the ordering is correct
                const int outlet_seg = segment_set[seg].outletSegment();
                wstate.segments[seg].pressure = wstate.segments[segment_set.segmentNumberToIndex(outlet_seg)].pressure;
            }
        }

        // Set segment surface_rates.
        // ...
        // TODO
        // ...
    }



    static void reportSingleWell(const PhaseUsage& pu,
                                 const int* globalCellIdxMap,
                                 const SingleWellState<NumActivePhases>& wstate,
                                 data::Well& dwell)
    {
        using rt = data::Rates::opt;

        dwell.bhp = wstate.bhp;
        dwell.thp = wstate.thp;
        dwell.temperature = wstate.temperature;

        if( pu.phase_used[BlackoilPhases::Aqua] ) {
            dwell.rates.set( rt::wat, wstate.surface_rates[ pu.phase_pos[BlackoilPhases::Aqua] ] );
        }

        if( pu.phase_used[BlackoilPhases::Liquid] ) {
            dwell.rates.set( rt::oil, wstate.surface_rates[ pu.phase_pos[BlackoilPhases::Liquid] ] );
        }

        if( pu.phase_used[BlackoilPhases::Vapour] ) {
            dwell.rates.set( rt::gas, wstate.surface_rates[ pu.phase_pos[BlackoilPhases::Vapour] ] );
        }

        // Make the connections vector.
        const int num_perf_well = wstate.connections.size();
        dwell.connections.resize(num_perf_well);
        for( int i = 0; i < num_perf_well; ++i ) {
            auto& connection = dwell.connections[ i ];
            // TODO
            // const auto active_index = this->well_perf_data_[well_index][i].cell_index;
            // connection.index = globalCellIdxMap[active_index];
            connection.pressure = wstate.connections[i].pressure;
            // TODO
            // connection.reservoir_rate = this->perfRates()[ itr.second[1] + i ];
        }

        // Make the segments map.
        for (const auto& seg : wstate.segments) {
            dwell.segments[seg.segment_number] = reportSegmentResults(pu, seg);
        }
    }



    static data::Segment
    reportSegmentResults(const PhaseUsage& pu, const typename SingleWellState<NumActivePhases>::Segment& seg)
    {
        auto seg_res = data::Segment{};

        seg_res.pressure = seg.pressure;

        if (pu.phase_used[BlackoilPhases::Aqua]) {
            seg_res.rates.set(data::Rates::opt::wat,
                              seg.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]]);
        }
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            seg_res.rates.set(data::Rates::opt::oil,
                              seg.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]]);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            seg_res.rates.set(data::Rates::opt::gas,
                              seg.surface_rates[pu.phase_pos[BlackoilPhases::Vapour]]);
        }

        seg_res.segNumber = seg.segment_number;

        return seg_res;
    }
};

} // namespace Opm

#endif // OPM_WELLANDGROUPSTATES_HEADER_INCLUDED
