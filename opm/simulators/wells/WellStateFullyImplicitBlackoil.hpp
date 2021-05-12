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
#include <opm/simulators/wells/ALQState.hpp>
#include <opm/simulators/wells/GlobalWellInfo.hpp>
#include <opm/simulators/wells/WellContainer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/simulators/wells/WellContainer.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

    /// The state of a set of wells, tailored for use by the fully
    /// implicit blackoil simulator.
    class WellStateFullyImplicitBlackoil
        : public WellState
    {
        typedef WellState  BaseType;
    public:
        static const uint64_t event_mask = ScheduleEvents::WELL_STATUS_CHANGE + ScheduleEvents::PRODUCTION_UPDATE + ScheduleEvents::INJECTION_UPDATE;
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
        using BaseType :: resetConnectionTransFactors;
        using BaseType :: updateStatus;

        explicit WellStateFullyImplicitBlackoil(const PhaseUsage& pu) :
            WellState(pu)
        {
        }


        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        void init(const std::vector<double>& cellPressures,
                  const Schedule& schedule,
                  const std::vector<Well>& wells_ecl,
                  const std::vector<ParallelWellInfo*>& parallel_well_info,
                  const int report_step,
                  const WellStateFullyImplicitBlackoil* prevState,
                  const std::vector<std::vector<PerforationData>>& well_perf_data,
                  const SummaryState& summary_state)
        {
            // call init on base class
            BaseType :: init(cellPressures, wells_ecl, parallel_well_info, well_perf_data, summary_state);
            this->global_well_info = std::make_optional<GlobalWellInfo>( schedule, report_step, wells_ecl );
            for (const auto& winfo: parallel_well_info)
            {
                well_rates.insert({winfo->name(), std::make_pair(winfo->isOwner(), std::vector<double>())});
            }

            const int nw = wells_ecl.size();

            if( nw == 0 ) return ;

            // Initialize perfphaserates_, which must be done here.
            const auto& pu = this->phaseUsage();
            const int np = pu.num_phases;

            int nperf = 0;
            for (const auto& wpd : well_perf_data) {
                nperf += wpd.size();
            }

            well_reservoir_rates_.resize(nw * this->numPhases(), 0.0);
            well_dissolved_gas_rates_.resize(nw, 0.0);
            well_vaporized_oil_rates_.resize(nw, 0.0);

            this->events_.clear();
            {
                const auto& wg_events = schedule[report_step].wellgroup_events();
                for (const auto& ecl_well : wells_ecl) {
                    const auto& wname = ecl_well.name();
                    if (wg_events.has(wname))
                        this->events_.add( wname, wg_events.at(wname) );
                    else
                        this->events_.add( wname, Events() );
                }
            }
            // Ensure that we start out with zero rates by default.
            perfphaserates_.clear();
            perfphaserates_.resize(nperf * this->numPhases(), 0.0);

            // these are only used to monitor the injectivity
            perf_water_throughput_.clear();
            perf_water_throughput_.resize(nperf, 0.0);
            perf_water_velocity_.clear();
            perf_water_velocity_.resize(nperf, 0.0);
            perf_skin_pressure_.clear();
            perf_skin_pressure_.resize(nperf, 0.0);

            num_perf_.resize(nw, 0);
            first_perf_index_.resize(nw, 0);
            first_perf_index_[0] = 0;
            for (int w = 0; w < nw; ++w) {
                // Initialize perfphaserates_ to well
                // rates divided by the number of perforations.
                const auto& wname = wells_ecl[w].name();
                const auto& well_info = this->wellMap().at(wname);
                const int connpos = well_info[1];
                const int num_perf_this_well = well_info[2];
                const int global_num_perf_this_well = parallel_well_info[w]->communication().sum(num_perf_this_well);
                auto& perf_press = this->perfPress(w);
                auto * phase_rates = &this->mutable_perfPhaseRates()[connpos * this->numPhases()];

                for (int perf = 0; perf < num_perf_this_well; ++perf) {
                    if (wells_ecl[w].getStatus() == Well::Status::OPEN) {
                        for (int p = 0; p < this->numPhases(); ++p) {
                            phase_rates[this->numPhases()*perf + p] = wellRates()[this->numPhases()*w + p] / double(global_num_perf_this_well);
                        }
                    }
                    perf_press[perf] = cellPressures[well_perf_data[w][perf].cell_index];
                }
                num_perf_[w] = num_perf_this_well;
                first_perf_index_[w] = connpos;
            }

            this->is_producer_.clear();
            for (int w = 0; w < nw; ++w) {
                const auto& ecl_well = wells_ecl[w];
                this->is_producer_.add( ecl_well.name(), ecl_well.isProducer());
            }

            current_injection_controls_.resize(nw, Well::InjectorCMode::CMODE_UNDEFINED);
            current_production_controls_.resize(nw, Well::ProducerCMode::CMODE_UNDEFINED);
            for (int w = 0; w < nw; ++w) {
                if (wells_ecl[w].isProducer()) {
                    const auto controls = wells_ecl[w].productionControls(summary_state);
                    currentProductionControl(w, controls.cmode);
                }
                else {
                    const auto controls = wells_ecl[w].injectionControls(summary_state);
                    currentInjectionControl(w, controls.cmode);
                }
            }

            perfRateSolvent_.clear();
            perfRateSolvent_.resize(nperf, 0.0);
            productivity_index_.resize(nw * this->numPhases(), 0.0);
            conn_productivity_index_.resize(nperf * this->numPhases(), 0.0);
            well_potentials_.resize(nw * this->numPhases(), 0.0);

            perfRatePolymer_.clear();
            perfRatePolymer_.resize(nperf, 0.0);

            perfRateBrine_.clear();
            perfRateBrine_.resize(nperf, 0.0);

            for (int w = 0; w < nw; ++w) {
                switch (wells_ecl[w].getStatus()) {
                case Well::Status::SHUT:
                    this->shutWell(w);
                    break;

                case Well::Status::STOP:
                    this->stopWell(w);
                    break;

                default:
                    this->openWell(w);
                    break;
                }
            }

            // intialize wells that have been there before
            // order may change so the mapping is based on the well name
            if (prevState && !prevState->wellMap().empty()) {
                auto end = prevState->wellMap().end();
                for (int w = 0; w < nw; ++w) {
                    const Well& well = wells_ecl[w];
                    if (well.getStatus() == Well::Status::SHUT) {
                        continue;
                    }

                    auto it = prevState->wellMap().find(well.name());
                    if (it != end)
                    {
                        const int newIndex = w;
                        const int oldIndex = it->second[ 0 ];
                        if (prevState->status_[oldIndex] == Well::Status::SHUT) {
                            // Well was shut in previous state, do not use its values.
                            continue;
                        }

                        if (is_producer_[newIndex] != prevState->is_producer_[oldIndex]) {
                            // Well changed to/from injector from/to producer, do not use its privious values.
                            continue;
                        }

                        // bhp
                        bhp()[ newIndex ] = prevState->bhp()[ oldIndex ];

                        // thp
                        thp()[ newIndex ] = prevState->thp()[ oldIndex ];

                        // If new target is set using WCONPROD, WCONINJE etc. we use the new control
                        if (!this->events_[w].hasEvent(WellStateFullyImplicitBlackoil::event_mask)) {
                            current_injection_controls_[ newIndex ] = prevState->currentInjectionControl( oldIndex );
                            current_production_controls_[ newIndex ] = prevState->currentProductionControl( oldIndex );
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

                        // Well potentials
                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            wellPotentials()[ idx ] = prevState->wellPotentials()[ oldidx ];
                        }

                        // perfPhaseRates
                        const int oldPerf_idx_beg = (*it).second[ 1 ];
                        const int num_perf_old_well = (*it).second[ 2 ];
                        const auto new_iter = this->wellMap().find(well.name());
                        if (new_iter == this->wellMap().end()) {
                            throw std::logic_error {
                                well.name() + " is not in internal well map - "
                                "Bug in WellStateFullyImplicitBlackoil"
                            };
                        }

                        const int connpos = new_iter->second[1];
                        const int num_perf_this_well = new_iter->second[2];

                        const int num_perf_changed = parallel_well_info[w]->communication()
                            .sum(static_cast<int>(num_perf_old_well != num_perf_this_well));
                        const bool global_num_perf_same = num_perf_changed == 0;


                        // copy perforation rates when the number of
                        // perforations is equal, otherwise initialize
                        // perfphaserates to well rates divided by the
                        // number of perforations.
                        if (global_num_perf_same)
                        {
                            const auto * src_rates = &prevState->perfPhaseRates()[oldPerf_idx_beg* np];
                            auto * target_rates = &this->mutable_perfPhaseRates()[connpos*np];
                            for (int perf_index = 0; perf_index < num_perf_this_well; perf_index++) {
                                for (int p = 0; p < np; p++) {
                                    target_rates[perf_index*np + p] = src_rates[perf_index*np + p];
                                }
                            }
                        } else {
                            const int global_num_perf_this_well = parallel_well_info[w]->communication().sum(num_perf_this_well);
                            auto * target_rates = &this->mutable_perfPhaseRates()[connpos*np];
                            for (int perf_index = 0; perf_index < num_perf_this_well; perf_index++) {
                                for (int p = 0; p < np; ++p) {
                                    target_rates[perf_index*np + p] = wellRates()[np*newIndex + p] / double(global_num_perf_this_well);
                                }
                            }
                        }

                        // perfPressures
                        if (global_num_perf_same)
                        {
                            auto& target_press = perfPress(w);
                            const auto& src_press = prevState->perfPress(well.name());
                            for (int perf = 0; perf < num_perf_this_well; ++perf)
                            {
                                target_press[perf] = src_press[perf];
                            }
                        }

                        // perfSolventRates
                        if (pu.has_solvent) {
                            if (global_num_perf_same)
                            {
                                int oldPerf_idx = oldPerf_idx_beg;
                                for (int perf = connpos; perf < connpos + num_perf_this_well; ++perf, ++oldPerf_idx )
                                {
                                    perfRateSolvent()[ perf ] = prevState->perfRateSolvent()[ oldPerf_idx ];
                                }
                            }
                        }

                        // polymer injectivity related
                        //
                        // here we did not consider the case that we close
                        // some perforation during the running and also,
                        // wells can be shut and re-opened
                        if (pu.has_polymermw) {
                            if (global_num_perf_same)
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

                        // Productivity index.
                        {
                            auto*       thisWellPI = &this     ->productivityIndex()[newIndex*np + 0];
                            const auto* thatWellPI = &prevState->productivityIndex()[oldIndex*np + 0];

                            for (int p = 0; p < np; ++p) {
                                thisWellPI[p] = thatWellPI[p];
                            }
                        }
                    }

                    // If in the new step, there is no THP related
                    // target/limit anymore, its thp value should be set to
                    // zero.
                    const bool has_thp = well.isInjector()
                        ? well.injectionControls (summary_state).hasControl(Well::InjectorCMode::THP)
                        : well.productionControls(summary_state).hasControl(Well::ProducerCMode::THP);

                    if (!has_thp) {
                        thp()[w] = 0.0;
                    }
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

            updateWellsDefaultALQ(wells_ecl);
            do_glift_optimization_ = true;
        }


        void resize(const std::vector<Well>& wells_ecl,
                    const std::vector<ParallelWellInfo*>& parallel_well_info,
                    const Schedule& schedule,
                    const bool handle_ms_well,
                    const size_t numCells,
                    const std::vector<std::vector<PerforationData>>& well_perf_data,
                    const SummaryState& summary_state)
        {
            const std::vector<double> tmp(numCells, 0.0); // <- UGLY HACK to pass the size
            init(tmp, schedule, wells_ecl, parallel_well_info, 0, nullptr, well_perf_data, summary_state);

            if (handle_ms_well) {
                initWellStateMSWell(wells_ecl, nullptr);
            }
        }

        /// One rate per phase and well connection.
        std::vector<double>& mutable_perfPhaseRates() { return perfphaserates_; }
        const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

        /// One current control per injecting well.
        std::vector<Opm::Well::InjectorCMode>& currentInjectionControls() { return current_injection_controls_; }
        const std::vector<Opm::Well::InjectorCMode>& currentInjectionControls() const { return current_injection_controls_; }
        void currentInjectionControl(std::size_t well_index, Well::InjectorCMode cmode) { current_injection_controls_[well_index] = cmode; }
        Well::InjectorCMode currentInjectionControl(std::size_t well_index) const { return current_injection_controls_[well_index]; }


        /// One current control per producing well.
        std::vector<Well::ProducerCMode>& currentProductionControls() { return current_production_controls_; }
        const std::vector<Well::ProducerCMode>& currentProductionControls() const { return current_production_controls_; }
        void currentProductionControl(std::size_t well_index, Well::ProducerCMode cmode) { current_production_controls_[well_index] = cmode; }
        Well::ProducerCMode currentProductionControl(std::size_t well_index) const { return current_production_controls_[well_index]; }

        void setCurrentWellRates(const std::string& wellName, const std::vector<double>& rates ) {
            well_rates[wellName].second = rates;
        }

        const std::vector<double>& currentWellRates(const std::string& wellName) const {
            auto it = well_rates.find(wellName);

            if (it == well_rates.end())
                OPM_THROW(std::logic_error, "Could not find any rates for well  " << wellName);

            return it->second.second;
        }

        bool hasWellRates(const std::string& wellName) const {
            return this->well_rates.find(wellName) != this->well_rates.end();
        }



        data::Wells
        report(const int* globalCellIdxMap,
               const std::function<bool(const int)>& wasDynamicallyClosed) const override
        {
            data::Wells res =
                WellState::report(globalCellIdxMap, wasDynamicallyClosed);

            const int nw = this->numWells();
            if (nw == 0) {
                return res;
            }

            const auto& pu = this->phaseUsage();
            const int np = pu.num_phases;

            using rt = data::Rates::opt;
            std::vector<rt> phs(np);
            if (pu.phase_used[Water]) {
                phs.at( pu.phase_pos[Water] ) = rt::wat;
            }

            if (pu.phase_used[Oil]) {
                phs.at( pu.phase_pos[Oil] ) = rt::oil;
            }

            if (pu.phase_used[Gas]) {
                phs.at( pu.phase_pos[Gas] ) = rt::gas;
            }

            // This is a reference or example on **how** to convert from
            // WellState to something understood by opm-common's output
            // layer.  It is intended to be properly implemented and
            // maintained as a part of simulators, as it relies on simulator
            // internals, details and representations.

            for (const auto& wt : this->wellMap()) {
                const auto w = wt.second[ 0 ];
                if (((this->status_[w] == Well::Status::SHUT) &&
                     ! wasDynamicallyClosed(w)) ||
                    ! this->parallel_well_info_[w]->isOwner())
                {
                    continue;
                }

                auto& well = res.at(wt.first);
                const int well_rate_index = w * pu.num_phases;

                if (pu.phase_used[Water]) {
                    const auto i = well_rate_index + pu.phase_pos[Water];
                    well.rates.set(rt::reservoir_water, this->well_reservoir_rates_[i]);
                    well.rates.set(rt::productivity_index_water, this->productivity_index_[i]);
                    well.rates.set(rt::well_potential_water, this->well_potentials_[i]);
                }

                if (pu.phase_used[Oil]) {
                    const auto i = well_rate_index + pu.phase_pos[Oil];
                    well.rates.set(rt::reservoir_oil, this->well_reservoir_rates_[i]);
                    well.rates.set(rt::productivity_index_oil, this->productivity_index_[i]);
                    well.rates.set(rt::well_potential_oil, this->well_potentials_[i]);
                }

                if (pu.phase_used[Gas]) {
                    const auto i = well_rate_index + pu.phase_pos[Gas];
                    well.rates.set(rt::reservoir_gas, this->well_reservoir_rates_[i]);
                    well.rates.set(rt::productivity_index_gas, this->productivity_index_[i]);
                    well.rates.set(rt::well_potential_gas, this->well_potentials_[i]);
                }

                if (pu.has_solvent || pu.has_zFraction) {
                    well.rates.set(rt::solvent, solventWellRate(w));
                }

                if (pu.has_polymer) {
                    well.rates.set(rt::polymer, polymerWellRate(w));
                }

                if (pu.has_brine) {
                    well.rates.set(rt::brine, brineWellRate(w));
                }

                if (is_producer_[w]) {
                    well.rates.set(rt::alq, getALQ(/*wellName=*/wt.first));
                }
                else {
                    well.rates.set(rt::alq, 0.0);
                }

                well.rates.set(rt::dissolved_gas, this->well_dissolved_gas_rates_[w]);
                well.rates.set(rt::vaporized_oil, this->well_vaporized_oil_rates_[w]);

                {
                    auto& curr = well.current_control;

                    curr.isProducer = this->is_producer_[w];
                    curr.prod = this->currentProductionControl(w);
                    curr.inj  = this->currentInjectionControl(w);
                }

                const auto nseg = this->numSegments(w);
                for (auto seg_ix = 0*nseg; seg_ix < nseg; ++seg_ix) {
                    const auto seg_no = this->segmentNumber(w, seg_ix);
                    well.segments[seg_no] =
                        this->reportSegmentResults(pu, w, seg_ix, seg_no);
                }
            }

            return res;
        }

        virtual void reportConnections(data::Well& well, const PhaseUsage &pu,
                                       const WellMapType::value_type& wt,
                                       const int* globalCellIdxMap) const override
        {
            using rt = data::Rates::opt;
            WellState::reportConnections(well, pu, wt, globalCellIdxMap);
            const int np = pu.num_phases;
            size_t local_comp_index = 0;
            std::vector< rt > phs( np );
            std::vector<rt> pi(np);
            if( pu.phase_used[Water] ) {
                phs.at( pu.phase_pos[Water] ) = rt::wat;
                pi .at( pu.phase_pos[Water] ) = rt::productivity_index_water;
            }

            if( pu.phase_used[Oil] ) {
                phs.at( pu.phase_pos[Oil] ) = rt::oil;
                pi .at( pu.phase_pos[Oil] ) = rt::productivity_index_oil;
            }

            if( pu.phase_used[Gas] ) {
                phs.at( pu.phase_pos[Gas] ) = rt::gas;
                pi .at( pu.phase_pos[Gas] ) = rt::productivity_index_gas;
            }
            for( auto& comp : well.connections) {
                const auto connPhaseOffset = np * (wt.second[1] + local_comp_index);

                const auto rates  = this->perfPhaseRates().begin() + connPhaseOffset;
                const auto connPI = this->connectionProductivityIndex().begin() + connPhaseOffset;

                for( int i = 0; i < np; ++i ) {
                    comp.rates.set( phs[ i ], *(rates  + i) );
                    comp.rates.set( pi [ i ], *(connPI + i) );
                }
                if ( pu.has_polymer ) {
                    comp.rates.set( rt::polymer, this->perfRatePolymer()[wt.second[1] + local_comp_index]);
                }
                if ( pu.has_brine ) {
                    comp.rates.set( rt::brine, this->perfRateBrine()[wt.second[1] + local_comp_index]);
                }
                if ( pu.has_solvent ) {
                    comp.rates.set( rt::solvent, this->perfRateSolvent()[wt.second[1] + local_comp_index]);
                }

                ++local_comp_index;
            }
            assert(local_comp_index == this->well_perf_data_[wt.second[0]].size());
        }

        /// init the MS well related.
        void initWellStateMSWell(const std::vector<Well>& wells_ecl,
                                 const WellStateFullyImplicitBlackoil* prev_well_state)
        {
            // still using the order in wells
            const int nw = wells_ecl.size();
            const auto& pu = this->phaseUsage();
            const int np = pu.num_phases;
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
            // in the init function, the well rates and perforation rates have been initialized or copied from prevState
            // what we do here, is to set the segment rates and perforation rates
            for (int w = 0; w < nw; ++w) {
                const auto& well_ecl = wells_ecl[w];
                const auto& wname = wells_ecl[w].name();
                const auto& well_info = this->wellMap().at(wname);
                const int connpos = well_info[1];
                const int num_perf_this_well = well_info[2];

                top_segment_index_.push_back(nseg_);
                if ( !well_ecl.isMultiSegment() ) { // not multi-segment well
                    nseg_ += 1;
                    seg_number_.push_back(1); // Assign single segment (top) as number 1.
                    seg_press_.push_back(bhp()[w]);
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
                        const int start_perf = connpos;

                        // make sure the information from wells_ecl consistent with wells
                        assert((n_activeperf == num_perf_this_well) &&
                               "Inconsistent number of reservoir connections in well");

                        if (pu.phase_used[Gas]) {
                            auto * perf_rates = &this->mutable_perfPhaseRates()[np * start_perf];
                            const int gaspos = pu.phase_pos[Gas];
                            // scale the phase rates for Gas to avoid too bad initial guess for gas fraction
                            // it will probably benefit the standard well too, while it needs to be justified
                            // TODO: to see if this strategy can benefit StandardWell too
                            // TODO: it might cause big problem for gas rate control or if there is a gas rate limit
                            // maybe the best way is to initialize the fractions first then get the rates
                            for (int perf = 0; perf < n_activeperf; perf++)
                                perf_rates[perf*np + gaspos] *= 100;
                        }

                        const auto * perf_rates = &perfPhaseRates()[np*start_perf];
                        std::vector<double> perforation_rates(perf_rates, perf_rates + num_perf_this_well*np);
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
                        const auto& perf_press = this->perfPress(w);
                        for (int seg = 1; seg < well_nseg; ++seg) {
                            if ( !segment_perforations[seg].empty() ) {
                                const int first_perf = segment_perforations[seg][0];
                                seg_press_.push_back(perf_press[first_perf]);
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
            }
            assert(int(seg_press_.size()) == nseg_);
            assert(int(seg_rates_.size()) == nseg_ * numPhases() );

            seg_pressdrop_.assign(nseg_, 0.);
            seg_pressdrop_hydorstatic_.assign(nseg_, 0.);
            seg_pressdrop_friction_.assign(nseg_, 0.);
            seg_pressdrop_acceleration_.assign(nseg_, 0.);

            if (prev_well_state && !prev_well_state->wellMap().empty()) {
                const auto& end = prev_well_state->wellMap().end();
                for (int w = 0; w < nw; ++w) {
                    const Well& well = wells_ecl[w];
                    if (well.getStatus() == Well::Status::SHUT)
                        continue;

                    const auto& it = prev_well_state->wellMap().find( wells_ecl[w].name() );
                    if (it != end) { // the well is found in the prev_well_state
                        // TODO: the well with same name can change a lot, like they might not have same number of segments
                        // we need to handle that later.
                        // for now, we just copy them.
                        const int old_index_well = (*it).second[0];
                        const int new_index_well = w;
                        if (prev_well_state->status_[old_index_well] == Well::Status::SHUT) {
                            continue;
                        }

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


        Events& events(std::size_t well_index) {
            return this->events_[well_index];
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
            return parallel_well_info_[w]->sumPerfValues(&perfRateSolvent_[0] + first_perf_index_[w],
                                                         &perfRateSolvent_[0] + first_perf_index_[w] + num_perf_[w]);
        }

        /// One rate pr well connection.
        std::vector<double>& perfRatePolymer() { return perfRatePolymer_; }
        const std::vector<double>& perfRatePolymer() const { return perfRatePolymer_; }

        /// One rate pr well
        double polymerWellRate(const int w) const {
            return parallel_well_info_[w]->sumPerfValues(&perfRatePolymer_[0] + first_perf_index_[w],
                                                         &perfRatePolymer_[0] + first_perf_index_[w] + num_perf_[w]);
        }

        /// One rate pr well connection.
        std::vector<double>& perfRateBrine() { return perfRateBrine_; }
        const std::vector<double>& perfRateBrine() const { return perfRateBrine_; }

        /// One rate pr well
        double brineWellRate(const int w) const {
            return parallel_well_info_[w]->sumPerfValues(&perfRateBrine_[0] + first_perf_index_[w],
                                                         &perfRateBrine_[0] + first_perf_index_[w] + num_perf_[w]);
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

        std::vector<double>& connectionProductivityIndex() {
            return this->conn_productivity_index_;
        }

        const std::vector<double>& connectionProductivityIndex() const {
            return this->conn_productivity_index_;
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

            auto* resv = &this->well_reservoir_rates_[np*well_index + 0];
            auto* wpi  = &this->productivity_index_[np*well_index + 0];

            for (int p = 0; p < np; ++p) {
                resv[p] = 0.0;
                wpi[p]  = 0.0;
            }

            const auto first = this->first_perf_index_[well_index]*np;
            const auto last  = first + this->num_perf_[well_index]*np;
            std::fill(this->conn_productivity_index_.begin() + first,
                      this->conn_productivity_index_.begin() + last, 0.0);
        }

        virtual void stopWell(int well_index) override {
            WellState::stopWell(well_index);
        }

        template<class Comm>
        void communicateGroupRates(const Comm& comm)
        {
            // Note that injection_group_vrep_rates is handled separate from
            // the forAllGroupData() function, since it contains single doubles,
            // not vectors.

            // Create a function that calls some function
            // for all the individual data items to simplify
            // the further code.
            auto iterateRatesContainer = [](auto& container, auto& func) {
                for (auto& x : container) {
                    if (x.second.first)
                    {
                        func(x.second.second);
                    }
                    else
                    {
                        // We might actually store non-zero values for
                        // distributed wells even if they are not owned.
                        std::vector<double> dummyRate;
                        dummyRate.assign(x.second.second.size(), 0);
                        func(dummyRate);
                    }
                }
            };


            // Compute the size of the data.
            std::size_t sz = 0;
            auto computeSize = [&sz](const auto& v) {
                sz += v.size();
            };
            iterateRatesContainer(this->well_rates, computeSize);
            sz += this->alq_state.pack_size();

            // Make a vector and collect all data into it.
            std::vector<double> data(sz);
            std::size_t pos = 0;
            auto collect = [&data, &pos](const auto& v) {
                for (const auto& x : v) {
                    data[pos++] = x;
                }
            };
            iterateRatesContainer(this->well_rates, collect);
            pos += this->alq_state.pack_data(&data[pos]);
            assert(pos == sz);

            // Communicate it with a single sum() call.
            comm.sum(data.data(), data.size());

            // Distribute the summed vector to the data items.
            pos = 0;
            auto distribute = [&data, &pos](auto& v) {
                for (auto& x : v) {
                    x = data[pos++];
                }
            };
            iterateRatesContainer(this->well_rates, distribute);
            pos += this->alq_state.unpack_data(&data[pos]);
            assert(pos == sz);
        }

        template<class Comm>
        void updateGlobalIsGrup(const Comm& comm)
        {
            this->global_well_info.value().update_group(this->status_.data(), this->current_injection_controls_, this->current_production_controls_);
            this->global_well_info.value().communicate(comm);
        }

        bool isInjectionGrup(const std::string& name) const {
            return this->global_well_info.value().in_injecting_group(name);
        }

        bool isProductionGrup(const std::string& name) const {
            return this->global_well_info.value().in_producing_group(name);
        }

        double getALQ( const std::string& name) const
        {
            return this->alq_state.get(name);
        }

        void setALQ( const std::string& name, double value)
        {
            this->alq_state.set(name, value);
        }

        bool gliftCheckAlqOscillation(const std::string &name) const {
            return this->alq_state.oscillation(name);
        }

        int gliftGetAlqDecreaseCount(const std::string &name) {
            return this->alq_state.get_decrement_count(name);
        }

        int gliftGetAlqIncreaseCount(const std::string &name) {
            return this->alq_state.get_increment_count(name);
        }

        void gliftUpdateAlqIncreaseCount(const std::string &name, bool increase) {
            this->alq_state.update_count(name, increase);
        }

        bool gliftOptimizationEnabled() const {
            return do_glift_optimization_;
        }

        void gliftTimeStepInit() {
            this->alq_state.reset_count();
            disableGliftOptimization();
        }

        void disableGliftOptimization() {
            do_glift_optimization_ = false;
        }

        void enableGliftOptimization() {
            do_glift_optimization_ = true;
        }

        int wellNameToGlobalIdx(const std::string &name) {
            return this->global_well_info.value().well_index(name);
        }

        std::string globalIdxToWellName(const int index) {
            return this->global_well_info.value().well_name(index);
        }


    private:
        std::vector<double> perfphaserates_;
        WellContainer<int> is_producer_;

        // vector with size number of wells +1.
        // iterate over all perforations of a given well
        // for (int perf = first_perf_index_[well_index]; perf < first_perf_index_[well_index] + num_perf_[well_index]; ++perf)
        std::vector<int> first_perf_index_;
        std::vector<int> num_perf_;
        std::vector<Opm::Well::InjectorCMode> current_injection_controls_;
        std::vector<Well::ProducerCMode> current_production_controls_;

        // Use of std::optional<> here is a technical crutch, the
        // WellStateFullyImplicitBlackoil class should be default constructible,
        // whereas the GlobalWellInfo is not.
        std::optional<GlobalWellInfo> global_well_info;
        std::map<std::string, std::pair<bool, std::vector<double>>> well_rates;

        ALQState alq_state;
        bool do_glift_optimization_;

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
        WellContainer<Events> events_;

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

        // Connection-level Productivity Index
        std::vector<double> conn_productivity_index_;

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

        // If the ALQ has changed since the previous report step,
        // reset current_alq and update default_alq. ALQ is used for
        // constant lift gas injection and for gas lift optimization
        // (THP controlled wells).
        //
        // NOTE: If a well is no longer used (e.g. it is shut down)
        // it is still kept in the maps "default_alq_" and "current_alq_". Since the
        // number of unused entries should be small (negligible memory
        // overhead) this is simpler than writing code to delete it.
        //
        void updateWellsDefaultALQ( const std::vector<Well>& wells_ecl )
        {
            const int nw = wells_ecl.size();
            for (int i = 0; i<nw; i++) {
                const Well &well = wells_ecl[i];
                if (well.isProducer()) {
                    // NOTE: This is the value set in item 12 of WCONPROD, or with WELTARG
                    auto alq = well.alq_value();
                    this->alq_state.update_default(well.name(), alq);
                }
            }
        }
    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
