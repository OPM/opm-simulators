/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2017 IRIS AS

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

#ifndef OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
#define OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED

#include <opm/simulators/wells/WellState.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <vector>
#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>

namespace Opm
{

    /// The state of a set of wells, tailored for use by the fully
    /// implicit blackoil simulator.
    class WellStateFullyImplicitBlackoil
        : public WellState
    {
        typedef WellState  BaseType;
    public:
        typedef BaseType :: WellMapType WellMapType;

        virtual ~WellStateFullyImplicitBlackoil() = default;

        // TODO: same definition with WellInterface, eventually they should go to a common header file.
        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;

        using BaseType :: wellRates;
        using BaseType :: bhp;
        using BaseType :: perfPress;
        using BaseType :: wellMap;
        using BaseType :: numWells;
        using BaseType :: numPhases;

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        void init(const std::vector<double>& cellPressures,
                  const Schedule& schedule,
                  const std::vector<Well>& wells_ecl,
                  const int report_step,
                  const WellStateFullyImplicitBlackoil* prevState,
                  const PhaseUsage& pu,
                  const std::vector<std::vector<PerforationData>>& well_perf_data,
                  const SummaryState& summary_state,
                  const int globalNumberOfWells)
        {
            // call init on base class
            BaseType :: init(cellPressures, wells_ecl, pu, well_perf_data, summary_state);

            globalIsInjectionGrup_.assign(globalNumberOfWells,0);
            globalIsProductionGrup_.assign(globalNumberOfWells,0);
            wellNameToGlobalIdx_.clear();

            const int nw = wells_ecl.size();

            if( nw == 0 ) return ;

            // Initialize perfphaserates_, which must be done here.
            const int np = pu.num_phases;

            int nperf = 0;
            for (const auto& wpd : well_perf_data) {
                nperf += wpd.size();
            }

            well_reservoir_rates_.resize(nw * np, 0.0);
            well_dissolved_gas_rates_.resize(nw, 0.0);
            well_vaporized_oil_rates_.resize(nw, 0.0);

            // checking whether some effective well control happens
            effective_events_occurred_.resize(nw, true);

            // a hack to make the resize() function used in RESTART related work
            if (!wells_ecl.empty() ) {
                // At the moment, the following events are considered to be effective events
                // more events might join as effective events
                // PRODUCTION_UPDATE, INJECTION_UPDATE, WELL_STATUS_CHANGE
                // 16 + 32 + 128
                const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                                                     + ScheduleEvents::PRODUCTION_UPDATE
                                                     + ScheduleEvents::INJECTION_UPDATE;
                for (int w = 0; w < nw; ++w) {
                    effective_events_occurred_[w]
                        = schedule.hasWellGroupEvent(wells_ecl[w].name(), effective_events_mask, report_step);
                }
            } // end of if (!well_ecl.empty() )

            // Ensure that we start out with zero rates by default.
            perfphaserates_.clear();
            perfphaserates_.resize(nperf * np, 0.0);

            perf_total_res_rates_.assign(nperf * np, 0.0);

            // these are only used to monitor the injectivity
            perf_water_throughput_.clear();
            perf_water_throughput_.resize(nperf, 0.0);
            perf_water_velocity_.clear();
            perf_water_velocity_.resize(nperf, 0.0);
            perf_skin_pressure_.clear();
            perf_skin_pressure_.resize(nperf, 0.0);

            int connpos = 0;
            first_perf_index_.resize(nw+1, 0);
            first_perf_index_[0] = connpos;
            for (int w = 0; w < nw; ++w) {
                // Initialize perfphaserates_ to well
                // rates divided by the number of perforations.
                const int num_perf_this_well = well_perf_data[w].size();
                for (int perf = connpos; perf < connpos + num_perf_this_well; ++perf) {
                    if (wells_ecl[w].getStatus() == Well::Status::OPEN) {
                        for (int p = 0; p < np; ++p) {
                            perfphaserates_[np*perf + p] = wellRates()[np*w + p] / double(num_perf_this_well);
                            // TODO: probably not good for gas rates
                            perf_total_res_rates_[np*perf + p] = perfphaserates_[np*perf + p];
                        }
                    }
                    perfPress()[perf] = cellPressures[well_perf_data[w][perf-connpos].cell_index];
                }
                connpos += num_perf_this_well;
                first_perf_index_[w+1] = connpos;
            }

            current_injection_controls_.resize(nw);
            current_production_controls_.resize(nw);

            perfRateSolvent_.clear();
            perfRateSolvent_.resize(nperf, 0.0);
            productivity_index_.resize(nw * np, 0.0);
            well_potentials_.resize(nw * np, 0.0);

            perfRatePolymer_.clear();
            perfRatePolymer_.resize(nperf, 0.0);

            perfRateBrine_.clear();
            perfRateBrine_.resize(nperf, 0.0);

            // intialize wells that have been there before
            // order may change so the mapping is based on the well name
            if (prevState && !prevState->wellMap().empty()) {
                connpos = 0;
                auto end = prevState->wellMap().end();
                for (int w = 0; w < nw; ++w) {
                    const Well& well = wells_ecl[w];
                    const int num_perf_this_well = well_perf_data[w].size();
                    auto it = prevState->wellMap().find(well.name());
                    if ( it != end )
                    {
                        const int oldIndex = (*it).second[ 0 ];
                        const int newIndex = w;

                        // bhp
                        bhp()[ newIndex ] = prevState->bhp()[ oldIndex ];

                        // thp
                        thp()[ newIndex ] = prevState->thp()[ oldIndex ];

                         // Currently this is taken care of by updateWellStateFromTarge. Maybe we should just remove the initialization and just use updateWellStateFromTarget
                         //if (effective_events_occurred_[w]) {
                         //   continue;
                         //}

                        // if there is no effective control event happens to the well, we use the current_injection/production_controls_ from prevState
                        // otherwise, we use the control specified in the deck
                        if (!effective_events_occurred_[w]) {
                            current_injection_controls_[ newIndex ] = prevState->currentInjectionControls()[ oldIndex ];
                            current_production_controls_[ newIndex ] = prevState->currentProductionControls()[ oldIndex ];
                        }

                        // wellrates
                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            wellRates()[ idx ] = prevState->wellRates()[ oldidx ];
                        }

                        // wellResrates
                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            wellReservoirRates()[ idx ] = prevState->wellReservoirRates()[ oldidx ];
                        }

                        // perfPhaseRates
                        const int oldPerf_idx_beg = (*it).second[ 1 ];
                        const int num_perf_old_well = (*it).second[ 2 ];
                        // copy perforation rates when the number of perforations is equal,
                        // otherwise initialize perfphaserates to well rates divided by the number of perforations.
                        if( num_perf_old_well == num_perf_this_well )
                        {
                            int old_perf_phase_idx = oldPerf_idx_beg *np;
                            for (int perf_phase_idx = connpos*np;
                                 perf_phase_idx < (connpos + num_perf_this_well)*np; ++perf_phase_idx, ++old_perf_phase_idx )
                            {
                                perfPhaseRates()[ perf_phase_idx ] = prevState->perfPhaseRates()[ old_perf_phase_idx ];
                                perfTotalResRates()[ perf_phase_idx ] = prevState->perfTotalResRates()[ old_perf_phase_idx ];
                            }
                        } else {
                            for (int perf = connpos; perf < connpos + num_perf_this_well; ++perf) {
                                for (int p = 0; p < np; ++p) {
                                    perfPhaseRates()[np*perf + p] = wellRates()[np*newIndex + p] / double(num_perf_this_well);
                                }
                            }
                        }
                        // perfPressures
                        if( num_perf_old_well == num_perf_this_well )
                        {
                            int oldPerf_idx = oldPerf_idx_beg;
                            for (int perf = connpos; perf < connpos + num_perf_this_well; ++perf, ++oldPerf_idx )
                            {
                                perfPress()[ perf ] = prevState->perfPress()[ oldPerf_idx ];
                            }
                        }
                        // perfSolventRates
                        if (pu.has_solvent) {
                            if( num_perf_old_well == num_perf_this_well )
                            {
                                int oldPerf_idx = oldPerf_idx_beg;
                                for (int perf = connpos; perf < connpos + num_perf_this_well; ++perf, ++oldPerf_idx )
                                {
                                    perfRateSolvent()[ perf ] = prevState->perfRateSolvent()[ oldPerf_idx ];
                                }
                            }
                        }

                        // polymer injectivity related
                        // here we did not consider the case that we close some perforation during the running
                        // and also, wells can be shut and re-opened
                        if (pu.has_polymermw) {
                            if( num_perf_old_well == num_perf_this_well )
                            {
                                int oldPerf_idx = oldPerf_idx_beg;
                                for (int perf = connpos; perf < connpos + num_perf_this_well; ++perf, ++oldPerf_idx )
                                {
                                    perf_water_throughput_[ perf ] = prevState->perfThroughput()[ oldPerf_idx ];
                                    perf_skin_pressure_[ perf ] = prevState->perfSkinPressure()[ oldPerf_idx ];
                                    perf_water_velocity_[ perf ] = prevState->perfWaterVelocity()[ oldPerf_idx ];
                                }
                            }
                        }
                    }

                    // If in the new step, there is no THP related target/limit anymore, its thp value should be
                    // set to zero.
                    const bool has_thp = well.isInjector() ? well.injectionControls(summary_state).hasControl(Well::InjectorCMode::THP)
                        : well.productionControls(summary_state).hasControl(Well::ProducerCMode::THP);
                    if (!has_thp) {
                        thp()[w] = 0.0;
                    }

                    // Increment connection position offset.
                    connpos += num_perf_this_well;
                }
            }

            {
                // we need to create a trival segment related values to avoid there will be some
                // multi-segment wells added later.
                nseg_ = nw;
                top_segment_index_.resize(nw);
                seg_number_.resize(nw);
                for (int w = 0; w < nw; ++w) {
                    top_segment_index_[w] = w;
                    seg_number_[w] = 1; // Top segment is segment #1
                }
                seg_press_ = bhp();
                seg_rates_ = wellRates();

                seg_pressdrop_.assign(nw, 0.);
                seg_pressdrop_hydorstatic_.assign(nw, 0.);
                seg_pressdrop_friction_.assign(nw, 0.);
                seg_pressdrop_acceleration_.assign(nw, 0.);
            }
        }


        void resize(const std::vector<Well>& wells_ecl,
                    const Schedule& schedule,
                    const bool handle_ms_well,
                    const size_t numCells,
                    const PhaseUsage& pu,
                    const std::vector<std::vector<PerforationData>>& well_perf_data,
                    const SummaryState& summary_state,
                    const int globalNumWells)
        {
            const std::vector<double> tmp(numCells, 0.0); // <- UGLY HACK to pass the size
            init(tmp, schedule, wells_ecl, 0, nullptr, pu, well_perf_data, summary_state, globalNumWells);

            if (handle_ms_well) {
                initWellStateMSWell(wells_ecl, pu, nullptr);
            }
        }

        /// One rate per phase and well connection.
        std::vector<double>& perfPhaseRates() { return perfphaserates_; }
        const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

        std::vector<double>& perfTotalResRates() { return perf_total_res_rates_; }
        const std::vector<double>& perfTotalResRates() const { return perf_total_res_rates_; }


        /// One current control per injecting well.
        std::vector<Opm::Well::InjectorCMode>& currentInjectionControls() { return current_injection_controls_; }
        const std::vector<Opm::Well::InjectorCMode>& currentInjectionControls() const { return current_injection_controls_; }


        /// One current control per producing well.
        std::vector<Well::ProducerCMode>& currentProductionControls() { return current_production_controls_; }
        const std::vector<Well::ProducerCMode>& currentProductionControls() const { return current_production_controls_; }

        bool hasProductionGroupControl(const std::string& groupName) const {
            return current_production_group_controls_.count(groupName) > 0;
        }

        bool hasInjectionGroupControl(const Opm::Phase& phase, const std::string& groupName) const {
            return current_injection_group_controls_.count(std::make_pair(phase, groupName)) > 0;
        }

        /// One current control per group.
        void setCurrentProductionGroupControl(const std::string& groupName, const Group::ProductionCMode& groupControl ) {
            current_production_group_controls_[groupName] = groupControl;
        }

        const Group::ProductionCMode& currentProductionGroupControl(const std::string& groupName) const {
            auto it = current_production_group_controls_.find(groupName);

            if (it == current_production_group_controls_.end())
                OPM_THROW(std::logic_error, "Could not find any control for production group " << groupName);

            return it->second;
        }

        /// One current control per group.
        void setCurrentInjectionGroupControl(const Opm::Phase& phase, const std::string& groupName, const Group::InjectionCMode& groupControl ) {
            current_injection_group_controls_[std::make_pair(phase, groupName)] = groupControl;
        }

        const Group::InjectionCMode& currentInjectionGroupControl(const Opm::Phase& phase, const std::string& groupName) const {
            auto it = current_injection_group_controls_.find(std::make_pair(phase, groupName));

            if (it == current_injection_group_controls_.end())
                OPM_THROW(std::logic_error, "Could not find any control for " << phase << " injection group " << groupName);

            return it->second;
        }
        
        void setCurrentWellRates(const std::string& wellName, const std::vector<double>& rates ) {
            well_rates[wellName] = rates;
        }

        const std::vector<double>& currentWellRates(const std::string& wellName) const {
            auto it = well_rates.find(wellName);

            if (it == well_rates.end())
                OPM_THROW(std::logic_error, "Could not find any rates for well  " << wellName);

            return it->second;
        }

        void setCurrentProductionGroupRates(const std::string& groupName, const std::vector<double>& rates ) {
            production_group_rates[groupName] = rates;
        }

        const std::vector<double>& currentProductionGroupRates(const std::string& groupName) const {
            auto it = production_group_rates.find(groupName);

            if (it == production_group_rates.end())
                OPM_THROW(std::logic_error, "Could not find any rates for productino group  " << groupName);

            return it->second;
        }
        
        void setCurrentProductionGroupReductionRates(const std::string& groupName, const std::vector<double>& target ) {
            production_group_reduction_rates[groupName] = target;
        }

        const std::vector<double>& currentProductionGroupReductionRates(const std::string& groupName) const {
            auto it = production_group_reduction_rates.find(groupName);

            if (it == production_group_reduction_rates.end())
                OPM_THROW(std::logic_error, "Could not find any reduction rates for production group  " << groupName);

            return it->second;
        }

        void setCurrentInjectionGroupReductionRates(const std::string& groupName, const std::vector<double>& target ) {
            injection_group_reduction_rates[groupName] = target;
        }

        const std::vector<double>& currentInjectionGroupReductionRates(const std::string& groupName) const {
            auto it = injection_group_reduction_rates.find(groupName);

            if (it == injection_group_reduction_rates.end())
                OPM_THROW(std::logic_error, "Could not find any reduction rates for injection group " << groupName);

            return it->second;
        }

        void setCurrentInjectionGroupReservoirRates(const std::string& groupName, const std::vector<double>& target ) {
            injection_group_reservoir_rates[groupName] = target;
        }

        const std::vector<double>& currentInjectionGroupReservoirRates(const std::string& groupName) const {
            auto it = injection_group_reservoir_rates.find(groupName);

            if (it == injection_group_reservoir_rates.end())
                OPM_THROW(std::logic_error, "Could not find any reservoir rates for injection group " << groupName);

            return it->second;
        }

        void setCurrentInjectionVREPRates(const std::string& groupName, const double& target ) {
            injection_group_vrep_rates[groupName] = target;
        }

        const double& currentInjectionVREPRates(const std::string& groupName) const {
            auto it = injection_group_vrep_rates.find(groupName);

            if (it == injection_group_vrep_rates.end())
                OPM_THROW(std::logic_error, "Could not find any VREP rates for group " << groupName);

            return it->second;
        }

        void setCurrentInjectionREINRates(const std::string& groupName, const std::vector<double>& target ) {
            injection_group_rein_rates[groupName] = target;
        }

        const std::vector<double>& currentInjectionREINRates(const std::string& groupName) const {
            auto it = injection_group_rein_rates.find(groupName);

            if (it == injection_group_rein_rates.end())
                OPM_THROW(std::logic_error, "Could not find any REIN rates for group " << groupName);

            return it->second;
        }

        void setCurrentGroupGratTargetFromSales(const std::string& groupName, const double& target ) {
            group_grat_target_from_sales[groupName] = target;
        }

        bool hasGroupGratTargetFromSales(const std::string& groupName) const {
            auto it = group_grat_target_from_sales.find(groupName);
            return it != group_grat_target_from_sales.end();
        }

        const double& currentGroupGratTargetFromSales(const std::string& groupName) const {
            auto it = group_grat_target_from_sales.find(groupName);

            if (it == group_grat_target_from_sales.end())
                OPM_THROW(std::logic_error, "Could not find any grat target from sales for group " << groupName);

            return it->second;
        }

        void setCurrentGroupInjectionPotentials(const std::string& groupName, const std::vector<double>& pot ) {
            injection_group_potentials[groupName] = pot;
        }

        const std::vector<double>& currentGroupInjectionPotentials(const std::string& groupName) const {
            auto it = injection_group_potentials.find(groupName);

            if (it == injection_group_potentials.end())
                OPM_THROW(std::logic_error, "Could not find any potentials for group " << groupName);

            return it->second;
        }



        data::Wells report(const PhaseUsage &pu, const int* globalCellIdxMap) const override
        {
            data::Wells res = WellState::report(pu, globalCellIdxMap);

            const int nw = this->numWells();
            if( nw == 0 ) return res;
            const int np = pu.num_phases;


            using rt = data::Rates::opt;
            std::vector< rt > phs( np );
            if( pu.phase_used[Water] ) {
                phs.at( pu.phase_pos[Water] ) = rt::wat;
            }

            if( pu.phase_used[Oil] ) {
                phs.at( pu.phase_pos[Oil] ) = rt::oil;
            }

            if( pu.phase_used[Gas] ) {
                phs.at( pu.phase_pos[Gas] ) = rt::gas;
            }

            /* this is a reference or example on **how** to convert from
             * WellState to something understood by opm-output. it is intended
             * to be properly implemented and maintained as a part of
             * simulators, as it relies on simulator internals, details and
             * representations.
             */

            for( const auto& wt : this->wellMap() ) {
                const auto w = wt.second[ 0 ];
                if (!this->open_for_output_[w])
                    continue;

                auto& well = res.at( wt.first );
                //well.injectionControl = static_cast<int>(this->currentInjectionControls()[ w ]);
                //well.productionControl = static_cast<int>(this->currentProductionControls()[ w ]);
                const int well_rate_index = w * pu.num_phases;

                if ( pu.phase_used[Water] ) {
                    well.rates.set( rt::reservoir_water, this->well_reservoir_rates_[well_rate_index + pu.phase_pos[Water]] );
                }

                if ( pu.phase_used[Oil] ) {
                    well.rates.set( rt::reservoir_oil, this->well_reservoir_rates_[well_rate_index + pu.phase_pos[Oil]] );
                }

                if ( pu.phase_used[Gas] ) {
                    well.rates.set( rt::reservoir_gas, this->well_reservoir_rates_[well_rate_index + pu.phase_pos[Gas]] );
                }

                if ( pu.phase_used[Water] ) {
                    well.rates.set( rt::productivity_index_water, this->productivity_index_[well_rate_index + pu.phase_pos[Water]] );
                }

                if ( pu.phase_used[Oil] ) {
                    well.rates.set( rt::productivity_index_oil, this->productivity_index_[well_rate_index + pu.phase_pos[Oil]] );
                }

                if ( pu.phase_used[Gas] ) {
                    well.rates.set( rt::productivity_index_gas, this->productivity_index_[well_rate_index + pu.phase_pos[Gas]] );
                }

                if ( pu.phase_used[Water] ) {
                    well.rates.set( rt::well_potential_water, this->well_potentials_[well_rate_index + pu.phase_pos[Water]] );
                }

                if ( pu.phase_used[Oil] ) {
                    well.rates.set( rt::well_potential_oil, this->well_potentials_[well_rate_index + pu.phase_pos[Oil]] );
                }

                if ( pu.phase_used[Gas] ) {
                    well.rates.set( rt::well_potential_gas, this->well_potentials_[well_rate_index + pu.phase_pos[Gas]] );
                }

                if ( pu.has_solvent ) {
                    well.rates.set( rt::solvent, solventWellRate(w) );
                }

                if ( pu.has_polymer ) {
                    well.rates.set( rt::polymer, polymerWellRate(w) );
                }

                if ( pu.has_brine ) {
                    well.rates.set( rt::brine, brineWellRate(w) );
                }

                well.rates.set( rt::dissolved_gas, this->well_dissolved_gas_rates_[w] );
                well.rates.set( rt::vaporized_oil, this->well_vaporized_oil_rates_[w] );

                {
                    auto& curr = well.current_control;

                    curr.prod = this->currentProductionControls()[w];
                    curr.inj  = this->currentInjectionControls() [w];
                }

                size_t local_comp_index = 0;
                for( auto& comp : well.connections) {
                    const auto rates = this->perfPhaseRates().begin()
                                     + (np * wt.second[ 1 ])
                                     + (np * local_comp_index);


                    for( int i = 0; i < np; ++i ) {
                        comp.rates.set( phs[ i ], *(rates + i) );
                    }
                    if ( pu.has_polymer ) {
                        comp.rates.set( rt::polymer, this->perfRatePolymer()[local_comp_index]);
                    }
                    if ( pu.has_brine ) {
                        comp.rates.set( rt::brine, this->perfRateBrine()[local_comp_index]);
                    }
                    if ( pu.has_solvent ) {
                        comp.rates.set( rt::solvent, this->perfRateSolvent()[local_comp_index]);
                    }

                    ++local_comp_index;
                }
                assert(local_comp_index == this->well_perf_data_[w].size());

                const auto nseg = this->numSegments(w);
                for (auto seg_ix = 0*nseg; seg_ix < nseg; ++seg_ix) {
                    const auto seg_no = this->segmentNumber(w, seg_ix);
                    well.segments[seg_no] =
                        this->reportSegmentResults(pu, w, seg_ix, seg_no);
                }
            }

            return res;
        }


        /// init the MS well related.
        void initWellStateMSWell(const std::vector<Well>& wells_ecl,
                                 const PhaseUsage& pu, const WellStateFullyImplicitBlackoil* prev_well_state)
        {
            // still using the order in wells
            const int nw = wells_ecl.size();
            if (nw == 0) {
                return;
            }

            top_segment_index_.clear();
            top_segment_index_.reserve(nw);
            seg_press_.clear();
            seg_press_.reserve(nw);
            seg_rates_.clear();
            seg_rates_.reserve(nw * numPhases());
            seg_number_.clear();

            nseg_ = 0;
            int connpos = 0;
            // in the init function, the well rates and perforation rates have been initialized or copied from prevState
            // what we do here, is to set the segment rates and perforation rates
            for (int w = 0; w < nw; ++w) {
                const auto& well_ecl = wells_ecl[w];
                int num_perf_this_well = well_perf_data_[w].size();
                top_segment_index_.push_back(nseg_);
                if ( !well_ecl.isMultiSegment() ) { // not multi-segment well
                    nseg_ += 1;
                    seg_number_.push_back(1); // Assign single segment (top) as number 1.
                    seg_press_.push_back(bhp()[w]);
                    const int np = numPhases();
                    for (int p = 0; p < np; ++p) {
                        seg_rates_.push_back(wellRates()[np * w + p]);
                    }
                } else { // it is a multi-segment well
                    const WellSegments& segment_set = well_ecl.getSegments();
                    // assuming the order of the perforations in well_ecl is the same with Wells
                    const WellConnections& completion_set = well_ecl.getConnections();
                    // number of segment for this single well
                    const int well_nseg = segment_set.size();
                    int n_activeperf = 0;
                    nseg_ += well_nseg;
                    for (auto segID = 0*well_nseg; segID < well_nseg; ++segID) {
                        this->seg_number_.push_back(segment_set[segID].segmentNumber());
                    }
                    // we need to know for each segment, how many perforation it has and how many segments using it as outlet_segment
                    // that is why I think we should use a well model to initialize the WellState here
                    std::vector<std::vector<int>> segment_perforations(well_nseg);
                    for (size_t perf = 0; perf < completion_set.size(); ++perf) {
                        const Connection& connection = completion_set.get(perf);
                        if (connection.state() == Connection::State::OPEN) {
                            const int segment_index = segment_set.segmentNumberToIndex(connection.segment());
                            segment_perforations[segment_index].push_back(n_activeperf);
                            n_activeperf++;
                        }
                    }

                    std::vector<std::vector<int>> segment_inlets(well_nseg);
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


                    // for the seg_rates_, now it becomes a recursive solution procedure.
                    {
                        const int np = numPhases();
                        const int start_perf = connpos;
                        const int start_perf_next_well = connpos + num_perf_this_well;

                        // make sure the information from wells_ecl consistent with wells
                        assert(n_activeperf == (start_perf_next_well - start_perf));

                        if (pu.phase_used[Gas]) {
                            const int gaspos = pu.phase_pos[Gas];
                            // scale the phase rates for Gas to avoid too bad initial guess for gas fraction
                            // it will probably benefit the standard well too, while it needs to be justified
                            // TODO: to see if this strategy can benefit StandardWell too
                            // TODO: it might cause big problem for gas rate control or if there is a gas rate limit
                            // maybe the best way is to initialize the fractions first then get the rates
                            for (int perf = 0; perf < n_activeperf; perf++) {
                                const int perf_pos = start_perf + perf;
                                perfPhaseRates()[np * perf_pos + gaspos] *= 100.;
                            }
                        }

                        const std::vector<double> perforation_rates(perfPhaseRates().begin() + np * start_perf,
                                                                    perfPhaseRates().begin() + np * start_perf_next_well); // the perforation rates for this well
                        std::vector<double> segment_rates;
                        calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, 0 /* top segment */, segment_rates);
                        std::copy(segment_rates.begin(), segment_rates.end(), std::back_inserter(seg_rates_));
                    }

                    // for the segment pressure, the segment pressure is the same with the first perforation belongs to the segment
                    // if there is no perforation associated with this segment, it uses the pressure from the outlet segment
                    // which requres the ordering is successful
                    // Not sure what is the best way to handle the initialization, hopefully, the bad initialization can be
                    // improved during the solveWellEq process
                    {
                        // top segment is always the first one, and its pressure is the well bhp
                        seg_press_.push_back(bhp()[w]);
                        const int top_segment = top_segment_index_[w];
                        const int start_perf = connpos;
                        for (int seg = 1; seg < well_nseg; ++seg) {
                            if ( !segment_perforations[seg].empty() ) {
                                const int first_perf = segment_perforations[seg][0];
                                seg_press_.push_back(perfPress()[start_perf + first_perf]);
                            } else {
                                // seg_press_.push_back(bhp); // may not be a good decision
                                // using the outlet segment pressure // it needs the ordering is correct
                                const int outlet_seg = segment_set[seg].outletSegment();
                                seg_press_.push_back(
                                    seg_press_[top_segment + segment_set.segmentNumberToIndex(outlet_seg)]);
                            }
                        }
                    }
                }
                connpos += num_perf_this_well;
            }
            assert(int(seg_press_.size()) == nseg_);
            assert(int(seg_rates_.size()) == nseg_ * numPhases() );

            seg_pressdrop_.assign(nseg_, 0.);
            seg_pressdrop_hydorstatic_.assign(nseg_, 0.);
            seg_pressdrop_friction_.assign(nseg_, 0.);
            seg_pressdrop_acceleration_.assign(nseg_, 0.);

            if (prev_well_state && !prev_well_state->wellMap().empty()) {
                // copying MS well related
                const auto& end = prev_well_state->wellMap().end();
                const int np = numPhases();
                for (int w = 0; w < nw; ++w) {
                    const auto& it = prev_well_state->wellMap().find( wells_ecl[w].name() );

                    if (it != end) { // the well is found in the prev_well_state
                        // TODO: the well with same name can change a lot, like they might not have same number of segments
                        // we need to handle that later.
                        // for now, we just copy them.
                        const int old_index_well = (*it).second[0];
                        const int new_index_well = w;
                        const int old_top_segment_index = prev_well_state->topSegmentIndex(old_index_well);
                        const int new_top_segmnet_index = topSegmentIndex(new_index_well);
                        int number_of_segment = 0;
                        // if it is the last well in list
                        if (new_index_well == int(top_segment_index_.size()) - 1) {
                            number_of_segment = nseg_ - new_top_segmnet_index;
                        } else {
                            number_of_segment = topSegmentIndex(new_index_well + 1) - new_top_segmnet_index;
                        }

                        for (int i = 0; i < number_of_segment * np; ++i) {
                            seg_rates_[new_top_segmnet_index * np + i] = prev_well_state->segRates()[old_top_segment_index * np + i];
                        }

                        for (int i = 0; i < number_of_segment; ++i) {
                            seg_press_[new_top_segmnet_index + i] = prev_well_state->segPress()[old_top_segment_index + i];
                        }
                    }
                }
            }
        }


        static void calculateSegmentRates(const std::vector<std::vector<int>>& segment_inlets, const std::vector<std::vector<int>>&segment_perforations,
                                          const std::vector<double>& perforation_rates, const int np, const int segment, std::vector<double>& segment_rates)
        {
            // the rate of the segment equals to the sum of the contribution from the perforations and inlet segment rates.
            // the first segment is always the top segment, its rates should be equal to the well rates.
            assert(segment_inlets.size() == segment_perforations.size());
            const int well_nseg = segment_inlets.size();
            if (segment == 0) { // beginning the calculation
                segment_rates.resize(np * well_nseg, 0.0);
            }
            // contributions from the perforations belong to this segment
            for (const int& perf : segment_perforations[segment]) {
                for (int p = 0; p < np; ++p) {
                    segment_rates[np * segment + p] += perforation_rates[np * perf + p];
                }
            }
            for (const int& inlet_seg : segment_inlets[segment]) {
                calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, inlet_seg, segment_rates);
                for (int p = 0; p < np; ++p) {
                    segment_rates[np * segment + p] += segment_rates[np * inlet_seg + p];
                }
            }
        }


        bool effectiveEventsOccurred(const int w) const {
            return effective_events_occurred_[w];
        }


        void setEffectiveEventsOccurred(const int w, const bool effective_events_occurred) {
            effective_events_occurred_[w] = effective_events_occurred;
        }

        const std::vector<int>& firstPerfIndex() const
        {
            return first_perf_index_;
        }

        /// One rate pr well connection.
        std::vector<double>& perfRateSolvent() { return perfRateSolvent_; }
        const std::vector<double>& perfRateSolvent() const { return perfRateSolvent_; }

        /// One rate pr well
        double solventWellRate(const int w) const {
            return std::accumulate(&perfRateSolvent_[0] + first_perf_index_[w], &perfRateSolvent_[0] + first_perf_index_[w+1], 0.0);
        }

        /// One rate pr well connection.
        std::vector<double>& perfRatePolymer() { return perfRatePolymer_; }
        const std::vector<double>& perfRatePolymer() const { return perfRatePolymer_; }

        /// One rate pr well
        double polymerWellRate(const int w) const {
            return std::accumulate(&perfRatePolymer_[0] + first_perf_index_[w], &perfRatePolymer_[0] + first_perf_index_[w+1], 0.0);
        }

        /// One rate pr well connection.
        std::vector<double>& perfRateBrine() { return perfRateBrine_; }
        const std::vector<double>& perfRateBrine() const { return perfRateBrine_; }

        /// One rate pr well
        double brineWellRate(const int w) const {
            return std::accumulate(&perfRateBrine_[0] + first_perf_index_[w], &perfRateBrine_[0] + first_perf_index_[w+1], 0.0);
        }

        std::vector<double>& wellReservoirRates()
        {
            return well_reservoir_rates_;
        }

        const std::vector<double>& wellReservoirRates() const
        {
            return well_reservoir_rates_;
        }

        std::vector<double>& wellDissolvedGasRates()
        {
            return well_dissolved_gas_rates_;
        }

        std::vector<double>& wellVaporizedOilRates()
        {
            return well_vaporized_oil_rates_;
        }

        const std::vector<double>& segRates() const
        {
            return seg_rates_;
        }

        std::vector<double>& segRates()
        {
            return seg_rates_;
        }

        const std::vector<double>& segPress() const
        {
            return seg_press_;
        }

        std::vector<double>& segPressDrop()
        {
            return seg_pressdrop_;
        }

        const std::vector<double>& segPressDrop() const
        {
            return seg_pressdrop_;
        }

        std::vector<double>& segPressDropFriction()
        {
            return seg_pressdrop_friction_;
        }

        const std::vector<double>& segPressDropFriction() const
        {
            return seg_pressdrop_friction_;
        }

        std::vector<double>& segPressDropHydroStatic()
        {
            return seg_pressdrop_hydorstatic_;
        }

        const std::vector<double>& segPressDropHydroStatic() const
        {
            return seg_pressdrop_hydorstatic_;
        }

        std::vector<double>& segPressDropAcceleration()
        {
            return seg_pressdrop_acceleration_;
        }

        const std::vector<double>& segPressDropAcceleration() const
        {
            return seg_pressdrop_acceleration_;
        }

        std::vector<double>& segPress()
        {
            return seg_press_;
        }

        int numSegment() const
        {
            return nseg_;
        }

        int topSegmentIndex(const int w) const
        {
            assert(w < int(top_segment_index_.size()) );

            return top_segment_index_[w];
        }

        std::vector<double>& productivityIndex() {
            return productivity_index_;
        }

        const std::vector<double>& productivityIndex() const {
            return productivity_index_;
        }

        std::vector<double>& wellPotentials() {
            return well_potentials_;
        }

        const std::vector<double>& wellPotentials() const {
            return well_potentials_;
        }

        std::vector<double>& perfThroughput() {
            return perf_water_throughput_;
        }

        const std::vector<double>& perfThroughput() const {
            return perf_water_throughput_;
        }

        std::vector<double>& perfSkinPressure() {
            return perf_skin_pressure_;
        }

        const std::vector<double>& perfSkinPressure() const {
            return perf_skin_pressure_;
        }

        std::vector<double>& perfWaterVelocity() {
            return perf_water_velocity_;
        }

        const std::vector<double>& perfWaterVelocity() const {
            return perf_water_velocity_;
        }

        virtual void shutWell(int well_index) override {
            WellState::shutWell(well_index);
            const int np = numPhases();
            for (int p = 0; p < np; ++p)
                this->well_reservoir_rates_[np * well_index + p] = 0;
        }

        template<class Comm>
        void communicateGroupRates(const Comm& comm) {
            // sum over all nodes
            for (auto& x : injection_group_rein_rates) {
                comm.sum(x.second.data(), x.second.size());
            }
            for (auto& x : injection_group_vrep_rates) {
                x.second = comm.sum(x.second);
            }
            for (auto& x : production_group_reduction_rates) {
                comm.sum(x.second.data(), x.second.size());
            }
            for (auto& x : injection_group_reduction_rates) {
                comm.sum(x.second.data(), x.second.size());
            }
            for (auto& x : injection_group_reservoir_rates) {
                comm.sum(x.second.data(), x.second.size());
            }
            for (auto& x : production_group_rates) {
                comm.sum(x.second.data(), x.second.size());
            }
            for (auto& x : well_rates) {
                comm.sum(x.second.data(), x.second.size());
            }
        }

        template<class Comm>
        void updateGlobalIsGrup(const Schedule& schedule, const int reportStepIdx, const Comm& comm)
        {
            std::fill(globalIsInjectionGrup_.begin(), globalIsInjectionGrup_.end(), 0);
            std::fill(globalIsProductionGrup_.begin(), globalIsProductionGrup_.end(), 0);
            int global_well_index = 0;
            const auto& end = wellMap().end();
            for (const auto& well : schedule.getWells(reportStepIdx)) {
                // Build global name->index map.
                wellNameToGlobalIdx_[well.name()] = global_well_index;

                // For wells on this process...
                const auto& it = wellMap().find( well.name());
                if (it != end) {
                    // ... set the GRUP/not GRUP states.
                    const int well_index = it->second[0];
                    if (!this->open_for_output_[well_index]) {
                        // Well is shut.
                        if (well.isInjector()) {
                            globalIsInjectionGrup_[global_well_index] = 0;
                        } else {
                            globalIsProductionGrup_[global_well_index] = 0;
                        }
                    } else {
                        if (well.isInjector()) {
                            globalIsInjectionGrup_[global_well_index] = (currentInjectionControls()[well_index] == Well::InjectorCMode::GRUP);
                        } else {
                            globalIsProductionGrup_[global_well_index] = (currentProductionControls()[well_index] == Well::ProducerCMode::GRUP);
                        }
                    }
                }
                ++global_well_index;
            }
            comm.sum(globalIsInjectionGrup_.data(), globalIsInjectionGrup_.size());
            comm.sum(globalIsProductionGrup_.data(), globalIsProductionGrup_.size());
        }

        bool isInjectionGrup(const std::string& name) const {

            auto it = wellNameToGlobalIdx_.find(name);

            if (it == wellNameToGlobalIdx_.end())
                OPM_THROW(std::logic_error, "Could not find global injection group for well" << name);

            return globalIsInjectionGrup_[it->second];
        }

        bool isProductionGrup(const std::string& name) const {

            auto it = wellNameToGlobalIdx_.find(name);

            if (it == wellNameToGlobalIdx_.end())
                OPM_THROW(std::logic_error, "Could not find global injection group for well" << name);

            return globalIsProductionGrup_[it->second];
        }

    private:
        std::vector<double> perfphaserates_;

        // vector with size number of wells +1.
        // iterate over all perforations of a given well
        // for (int perf = first_perf_index_[well_index]; perf < first_perf_index_[well_index+1]; ++perf)
        std::vector<int> first_perf_index_;

        // perforation total reservoir volume rates
        std::vector<double> perf_total_res_rates_;

        std::vector<Opm::Well::InjectorCMode> current_injection_controls_;
        std::vector<Well::ProducerCMode> current_production_controls_;

        // size of global number of wells
        std::vector<int> globalIsInjectionGrup_;
        std::vector<int> globalIsProductionGrup_;
        std::map<std::string, int> wellNameToGlobalIdx_;

        std::map<std::string, Group::ProductionCMode> current_production_group_controls_;
        std::map<std::pair<Opm::Phase, std::string>, Group::InjectionCMode> current_injection_group_controls_;

        std::map<std::string, std::vector<double>> well_rates;
        std::map<std::string, std::vector<double>> production_group_rates;
        std::map<std::string, std::vector<double>> production_group_reduction_rates;
        std::map<std::string, std::vector<double>> injection_group_reduction_rates;
        std::map<std::string, std::vector<double>> injection_group_reservoir_rates;
        std::map<std::string, std::vector<double>> injection_group_potentials;
        std::map<std::string, double> injection_group_vrep_rates;
        std::map<std::string, std::vector<double>> injection_group_rein_rates;
        std::map<std::string, double> group_grat_target_from_sales;

        std::vector<double> perfRateSolvent_;

        // only for output
        std::vector<double> perfRatePolymer_;
        std::vector<double> perfRateBrine_;

        // it is the throughput of water flow through the perforations
        // it is used as a measure of formation damage around well-bore due to particle deposition
        // it will only be used for injectors to check the injectivity
        std::vector<double> perf_water_throughput_;

        // skin pressure of peforation
        // it will only be used for injectors to check the injectivity
        std::vector<double> perf_skin_pressure_;

        // it will only be used for injectors to check the injectivity
        // water velocity of perforation
        std::vector<double> perf_water_velocity_;

        // phase rates under reservoir condition for wells
        // or voidage phase rates
        std::vector<double> well_reservoir_rates_;

        // dissolved gas rates or solution gas production rates
        // should be zero for injection wells
        std::vector<double> well_dissolved_gas_rates_;

        // vaporized oil rates or solution oil producation rates
        // should be zero for injection wells
        std::vector<double> well_vaporized_oil_rates_;

        // some events happens to the well, like this well is a new well
        // or new well control keywords happens
        // \Note: for now, only WCON* keywords, and well status change is considered
        std::vector<bool> effective_events_occurred_;

        // MS well related
        // for StandardWell, the number of segments will be one
        std::vector<double> seg_rates_;
        std::vector<double> seg_press_;
        // The following data are only recorded for output
        // pressure drop
        std::vector<double> seg_pressdrop_;
        // frictional pressure drop
        std::vector<double> seg_pressdrop_friction_;
        // hydrostatic pressure drop
        std::vector<double> seg_pressdrop_hydorstatic_;
        // accelerational pressure drop
        std::vector<double> seg_pressdrop_acceleration_;
        // the index of the top segments, which is used to locate the
        // multisegment well related information in WellState
        std::vector<int> top_segment_index_;
        int nseg_; // total number of the segments

        // Productivity Index
        std::vector<double> productivity_index_;

        // Well potentials
        std::vector<double> well_potentials_;

        /// Map segment index to segment number, mostly for MS wells.
        ///
        /// Segment number (one-based) of j-th segment in i-th well is
        /// \code
        ///    const auto top    = topSegmentIndex(i);
        ///    const auto seg_No = seg_number_[top + j];
        /// \end
        std::vector<int> seg_number_;

        ::Opm::data::Segment
        reportSegmentResults(const PhaseUsage& pu,
                             const int         well_id,
                             const int         seg_ix,
                             const int         seg_no) const
        {
            auto seg_res = ::Opm::data::Segment{};

            const auto seg_dof =
                this->topSegmentIndex(well_id) + seg_ix;

            const auto* rate =
                &this->segRates()[seg_dof * this->numPhases()];

            {
                using Value =::Opm::data::SegmentPressures::Value;
                auto& segpress = seg_res.pressures;
                segpress[Value::Pressure] = this->segPress()[seg_dof];
                segpress[Value::PDrop] = this->segPressDrop()[seg_dof];
                segpress[Value::PDropHydrostatic] = this->segPressDropHydroStatic()[seg_dof];
                segpress[Value::PDropFriction] = this->segPressDropFriction()[seg_dof];
                segpress[Value::PDropAccel] = this->segPressDropAcceleration()[seg_dof];
            }

            if (pu.phase_used[Water]) {
                seg_res.rates.set(data::Rates::opt::wat,
                                  rate[pu.phase_pos[Water]]);
            }

            if (pu.phase_used[Oil]) {
                seg_res.rates.set(data::Rates::opt::oil,
                                  rate[pu.phase_pos[Oil]]);
            }

            if (pu.phase_used[Gas]) {
                seg_res.rates.set(data::Rates::opt::gas,
                                  rate[pu.phase_pos[Gas]]);
            }

            seg_res.segNumber = seg_no;

            return seg_res;
        }

        int numSegments(const int well_id) const
        {
            const auto topseg = this->topSegmentIndex(well_id);

            return (well_id + 1 == this->numWells()) // Last well?
                ? (this->numSegment() - topseg)
                : (this->topSegmentIndex(well_id + 1) - topseg);
        }

        int segmentNumber(const int well_id, const int seg_id) const
        {
            const auto top_offset = this->topSegmentIndex(well_id);

            return this->seg_number_[top_offset + seg_id];
        }
    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
