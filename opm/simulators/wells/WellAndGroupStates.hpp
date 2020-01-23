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
        assert(wells.size() == well_perf_data.size());
        well_perf_data_ = well_perf_data;

        // Initialize all wells.
        const int nw = wells.size();
        well_states_.resize(nw);
        well_names_.resize(nw);
        for (int w = 0; w < nw; ++w) {
            well_names_[w] = wells[w].name();
            initSingleWell(cell_pressures, schedule, wells[w], report_step, phase_usage, well_perf_data_[w], summary_state, well_states_[w]);
        }

        // For wells present in a previous state, override initialization where appropriate.
        if (prev_state) {
            for (int w = 0; w < nw; ++w) {
                well_names_[w] = wells[w].name();
                auto it = std::find(prev_state->well_names_.begin(), prev_state->well_names_.end(), well_names_[w]);
                if (it != prev_state->well_names_.end()) {
                    const int w_in_prev = it - prev_state->well_names_.end();
                    initSingleWellFromPrevious(prev_state->well_states_[w_in_prev], well_states_[w]);
                }
            }
        }

        // Group states remain default-initialized.
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
            reportSingleWell(phase_usage, globalCellIdxMap, well_perf_data_[w], well_states_[w], dwell);
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

    std::vector<std::vector<PerforationData>> well_perf_data_;


    // -----------  Private functions  -----------

    template <unsigned long N>
    static double sum(const std::array<double, N>& x)
    {
        return std::accumulate(x.begin(), x.end(), 0.0);
    }

    static void initSingleWell(const std::vector<double>& cell_pressures,
                               const Schedule& schedule,
                               const Well& well,
                               const int report_step,
                               const PhaseUsage& phase_usage,
                               const std::vector<PerforationData>& perf_data,
                               const SummaryState& summary_state,
                               SingleWellState<NumActivePhases>& wstate)
    {
        const auto inj_controls = well.isInjector() ? well.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = well.isProducer() ? well.productionControls(summary_state) : Well::ProductionControls(0);

        // Initialize flag, statuses and current controls.
        initFlagsAndControls(schedule, well, report_step, inj_controls, prod_controls, wstate);

        // Initialize bhp, thp and well rates.
        initPressuresAndRates(cell_pressures, well, phase_usage, perf_data, inj_controls, prod_controls, wstate);

        // Initialize connection pressures and rates.
        initConnections(cell_pressures, perf_data, wstate);

        // Initialize multi-segment well parts (the 'segments' member).
        if (well.isMultiSegment()) {
            initMultiSegment(well, wstate);
        }
    }




    static void initSingleWellFromPrevious(const SingleWellState<NumActivePhases>& prev_wstate,
                                           SingleWellState<NumActivePhases>& wstate)
    {
        // Flags and statuses.
        wstate.status = prev_wstate.status;
        wstate.is_producer = prev_wstate.is_producer;
        wstate.effective_events_occurred = prev_wstate.effective_events_occurred;
        wstate.current_injection_control = prev_wstate.current_injection_control;
        wstate.current_production_control = prev_wstate.current_production_control;

        // Quantities.
        wstate.bhp = prev_wstate.bhp;
        wstate.thp = prev_wstate.thp;
        wstate.temperature = prev_wstate.temperature; // 20 degrees Celcius by default.
        wstate.surface_rates = prev_wstate.surface_rates;
        wstate.reservoir_rates = prev_wstate.reservoir_rates;
        wstate.dissolved_gas_rate = prev_wstate.dissolved_gas_rate;
        wstate.vaporized_oil_rate = prev_wstate.vaporized_oil_rate;
        wstate.potentials = prev_wstate.potentials;

        // Connections.
        if (wstate.connections.size() == prev_wstate.connections.size()) {
            wstate.connections = prev_wstate.connections;
        }

        // Segments.
        if (wstate.segments.size() == prev_wstate.segments.size()) {
            wstate.segments = prev_wstate.segments;
        }
    }




    static void initFlagsAndControls(const Schedule& schedule,
                                     const Well& well,
                                     const int report_step,
                                     const Well::InjectionControls& inj_controls,
                                     const Well::ProductionControls& prod_controls,
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


        if (well.isInjector()) {
            wstate.current_injection_control = inj_controls.cmode;
        } else {
            wstate.current_production_control = prod_controls.cmode;
        }
    }




    static void initPressuresAndRates(const std::vector<double>& cell_pressures,
                                      const Well& well,
                                      const PhaseUsage& phase_usage,
                                      const std::vector<PerforationData>& perf_data,
                                      const Well::InjectionControls& inj_controls,
                                      const Well::ProductionControls& prod_controls,
                                      SingleWellState<NumActivePhases>& wstate)
    {
        if (perf_data.empty()) {
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
    }




    static void initConnections(const std::vector<double>& cell_pressures,
                                const std::vector<PerforationData>& perf_data,
                                SingleWellState<NumActivePhases>& wstate)
    {
        const int num_conn = perf_data.size();
        wstate.connections.resize(num_conn);
        for (int conn = 0; conn < num_conn; ++conn) {
            auto& connection = wstate.connections[conn];
            connection.pressure = cell_pressures[perf_data[conn].cell_index];
            connection.surface_rates = wstate.surface_rates;
            for (double& q : connection.surface_rates) {
                q /= num_conn;
            }
        }
    }




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
        calculateSegmentRates(segment_inlets, segment_perforations, 0, wstate);
        assert(wstate.segments[0].surface_rates == wstate.surface_rates);
    }



    static void calculateSegmentRates(const std::vector<std::vector<int>>& segment_inlets,
                                      const std::vector<std::vector<int>>& segment_connections,
                                      const int segment,
                                      SingleWellState<NumActivePhases>& wstate)
    {
        // This assumes that the segment rates have been initialized to zero.

        // the rate of the segment equals to the sum of the contribution from the connections and inlet segment rates.
        // the first segment is always the top segment, its rates should be equal to the well rates.
        assert(segment_inlets.size() == segment_connections.size());
        const int np = wstate.connections.front().surface_rates.size();

        // contributions from the connections belong to this segment
        for (const int& perf : segment_connections[segment]) {
            for (int p = 0; p < np; ++p) {
                wstate.segments[segment].surface_rates[p] += wstate.connections[perf].surface_rates[p];
            }
        }
        for (const int& inlet_seg : segment_inlets[segment]) {
            calculateSegmentRates(segment_inlets, segment_connections, inlet_seg, wstate);
            for (int p = 0; p < np; ++p) {
                wstate.segments[segment].surface_rates[p] += wstate.segments[inlet_seg].surface_rates[p];
            }
        }
    }



    static void reportSingleWell(const PhaseUsage& pu,
                                 const int* globalCellIdxMap,
                                 const std::vector<PerforationData>& perf_data,
                                 const SingleWellState<NumActivePhases>& wstate,
                                 data::Well& dwell)
    {
        dwell.bhp = wstate.bhp;
        dwell.thp = wstate.thp;
        dwell.temperature = wstate.temperature;

        using rt = data::Rates::opt;
        std::vector<rt> phs(pu.num_phases);

        if( pu.phase_used[BlackoilPhases::Aqua] ) {
            const int apos = pu.phase_pos[BlackoilPhases::Aqua];
            phs[apos] = rt::wat;
            dwell.rates.set( rt::wat, wstate.surface_rates[apos] );
            dwell.rates.set( rt::reservoir_water, wstate.reservoir_rates[apos] );
            dwell.rates.set( rt::productivity_index_water, wstate.productivity_index[apos] );
            dwell.rates.set( rt::well_potential_water, wstate.potentials[apos] );
        }

        if( pu.phase_used[BlackoilPhases::Liquid] ) {
            const int lpos = pu.phase_pos[BlackoilPhases::Liquid];
            phs[lpos] = rt::oil;
            dwell.rates.set( rt::oil, wstate.surface_rates[lpos] );
            dwell.rates.set( rt::reservoir_oil, wstate.reservoir_rates[lpos] );
            dwell.rates.set( rt::productivity_index_oil, wstate.productivity_index[lpos] );
            dwell.rates.set( rt::well_potential_oil, wstate.potentials[lpos] );
        }

        if( pu.phase_used[BlackoilPhases::Vapour] ) {
            const int vpos = pu.phase_pos[BlackoilPhases::Vapour];
            phs[vpos] = rt::gas;
            dwell.rates.set( rt::gas, wstate.surface_rates[vpos] );
            dwell.rates.set( rt::reservoir_gas, wstate.reservoir_rates[vpos] );
            dwell.rates.set( rt::productivity_index_gas, wstate.productivity_index[vpos] );
            dwell.rates.set( rt::well_potential_gas, wstate.potentials[vpos] );
        }

        if ( pu.has_solvent ) {
            double solvent_rate = 0.0;
            for (const auto& conn : wstate.connections) {
                solvent_rate += conn.solvent_rate;
            }
            dwell.rates.set( rt::solvent, solvent_rate );
        }

        dwell.rates.set( rt::dissolved_gas, wstate.dissolved_gas_rate );
        dwell.rates.set( rt::vaporized_oil, wstate.vaporized_oil_rate );

        // Make the connections vector.
        assert(perf_data.size() == wstate.connections.size());
        const int num_conn = wstate.connections.size();
        dwell.connections.resize(num_conn);
        for(int i = 0; i < num_conn; ++i) {
            auto& connection = dwell.connections[i];
            const auto active_index = perf_data[i].cell_index;
            connection.index = globalCellIdxMap[active_index];
            connection.pressure = wstate.connections[i].pressure;
            connection.reservoir_rate = sum(wstate.connections[i].reservoir_rates);
            for (int ph = 0; ph < pu.num_phases; ++ph) {
                connection.rates.set(phs[ph], wstate.connections[i].surface_rates[ph]);
            }
        }

        // Make the segments map.
        for (const auto& seg : wstate.segments) {
            dwell.segments[seg.segment_number] = reportSegmentResults(pu, seg);
        }
    }



    static data::Segment
    reportSegmentResults(const PhaseUsage& pu, const typename SingleWellState<NumActivePhases>::Segment& seg)
    {
        using rt = data::Rates::opt;

        auto seg_res = data::Segment{};

        seg_res.pressure = seg.pressure;

        if (pu.phase_used[BlackoilPhases::Aqua]) {
            seg_res.rates.set(rt::wat, seg.surface_rates[pu.phase_pos[BlackoilPhases::Aqua]]);
        }
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            seg_res.rates.set(rt::oil, seg.surface_rates[pu.phase_pos[BlackoilPhases::Liquid]]);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            seg_res.rates.set(rt::gas, seg.surface_rates[pu.phase_pos[BlackoilPhases::Vapour]]);
        }

        seg_res.segNumber = seg.segment_number;

        return seg_res;
    }
};

} // namespace Opm

#endif // OPM_WELLANDGROUPSTATES_HEADER_INCLUDED
