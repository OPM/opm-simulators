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

#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <vector>
#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <array>

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
        void init(const Wells* wells, const std::vector<double>& cellPressures,
                  const std::vector<const Well*>& wells_ecl, const int report_step,
                  const WellStateFullyImplicitBlackoil* prevState, const PhaseUsage& pu)
        {
            // call init on base class
            BaseType :: init(wells, cellPressures);

            // if there are no well, do nothing in init
            if (wells == 0) {
                return;
            }

            const int nw = wells->number_of_wells;

            if( nw == 0 ) return ;

            // Initialize perfphaserates_, which must be done here.
            const int np = wells->number_of_phases;
            const int nperf = wells->well_connpos[nw];

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

                for (int w = 0; w <nw; ++w) {
                    const int nw_wells_ecl = wells_ecl.size();
                    int index_well_ecl = 0;
                    const std::string well_name(wells->name[w]);
                    for (; index_well_ecl < nw_wells_ecl; ++index_well_ecl) {
                        if (well_name == wells_ecl[index_well_ecl]->name()) {
                            break;
                        }
                    }

                    // It should be able to find in wells_ecl.
                    if (index_well_ecl == nw_wells_ecl) {
                        OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
                    }

                    const Well* well_ecl = wells_ecl[index_well_ecl];
                    effective_events_occurred_[w] = (well_ecl->hasEvent(effective_events_mask, report_step) );
                }
            } // end of if (!well_ecl.empty() )

            // Ensure that we start out with zero rates by default.
            perfphaserates_.clear();
            perfphaserates_.resize(nperf * np, 0.0);

            // these are only used to monitor the injectivity
            perf_water_throughput_.clear();
            perf_water_throughput_.resize(nperf, 0.0);
            perf_water_velocity_.clear();
            perf_water_velocity_.resize(nperf, 0.0);
            perf_skin_pressure_.clear();
            perf_skin_pressure_.resize(nperf, 0.0);

            for (int w = 0; w < nw; ++w) {
                assert((wells->type[w] == INJECTOR) || (wells->type[w] == PRODUCER));
                const WellControls* ctrl = wells->ctrls[w];

                if (well_controls_well_is_stopped(ctrl)) {
                    // Shut well: perfphaserates_ are all zero.
                } else {
                    const int num_perf_this_well = wells->well_connpos[w + 1] - wells->well_connpos[w];
                    // Open well: Initialize perfphaserates_ to well
                    // rates divided by the number of perforations.
                    for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                        for (int p = 0; p < np; ++p) {
                            perfphaserates_[np*perf + p] = wellRates()[np*w + p] / double(num_perf_this_well);
                        }
                        perfPress()[perf] = cellPressures[wells->well_cells[perf]];
                    }
                }
            }

            current_controls_.resize(nw);
            // The controls set in the Wells (specified in the DECK) are treated as default initial value
            for (int w = 0; w < nw; ++w) {
                current_controls_[w] = well_controls_get_current(wells->ctrls[w]);
            }
            perfRateSolvent_.clear();
            perfRateSolvent_.resize(nperf, 0.0);
            productivity_index_.resize(nw * np, 0.0);
            well_potentials_.resize(nw * np, 0.0);

            // intialize wells that have been there before
            // order may change so the mapping is based on the well name
            if(prevState && !prevState->wellMap().empty()) {
                typedef typename WellMapType :: const_iterator const_iterator;
                const_iterator end = prevState->wellMap().end();
                for (int w = 0; w < nw; ++w) {
                    const std::string name( wells->name[ w ] );
                    const_iterator it = prevState->wellMap().find( name );
                    if( it != end )
                    {
                        const int oldIndex = (*it).second[ 0 ];
                        const int newIndex = w;

                        // bhp
                        bhp()[ newIndex ] = prevState->bhp()[ oldIndex ];

                        // thp
                        thp()[ newIndex ] = prevState->thp()[ oldIndex ];

                        // if there is no effective control event happens to the well, we use the current_controls_ from prevState
                        // otherwise, we use the control specified in the deck
                        if (!effective_events_occurred_[w]) {
                            current_controls_[ newIndex ] = prevState->currentControls()[ oldIndex ];
                            // also change the one in the WellControls
                            well_controls_set_current(wells->ctrls[w], current_controls_[ newIndex ]);
                        }

                        // wellrates
                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            wellRates()[ idx ] = prevState->wellRates()[ oldidx ];
                        }

                        // perfPhaseRates
                        const int oldPerf_idx_beg = (*it).second[ 1 ];
                        const int num_perf_old_well = (*it).second[ 2 ];
                        const int num_perf_this_well = wells->well_connpos[newIndex + 1] - wells->well_connpos[newIndex];
                        // copy perforation rates when the number of perforations is equal,
                        // otherwise initialize perfphaserates to well rates divided by the number of perforations.
                        if( num_perf_old_well == num_perf_this_well )
                        {
                            int old_perf_phase_idx = oldPerf_idx_beg *np;
                            for (int perf_phase_idx = wells->well_connpos[ newIndex ]*np;
                                 perf_phase_idx < wells->well_connpos[ newIndex + 1]*np; ++perf_phase_idx, ++old_perf_phase_idx )
                            {
                                perfPhaseRates()[ perf_phase_idx ] = prevState->perfPhaseRates()[ old_perf_phase_idx ];
                            }
                        } else {
                            for (int perf = wells->well_connpos[newIndex]; perf < wells->well_connpos[newIndex + 1]; ++perf) {
                                for (int p = 0; p < np; ++p) {
                                    perfPhaseRates()[np*perf + p] = wellRates()[np*newIndex + p] / double(num_perf_this_well);
                                }
                            }
                        }
                        // perfPressures
                        if( num_perf_old_well == num_perf_this_well )
                        {
                            int oldPerf_idx = oldPerf_idx_beg;
                            for (int perf = wells->well_connpos[ newIndex ];
                                 perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
                            {
                                perfPress()[ perf ] = prevState->perfPress()[ oldPerf_idx ];
                            }
                        }
                        // perfSolventRates
                        if (pu.has_solvent) {
                            if( num_perf_old_well == num_perf_this_well )
                            {
                                int oldPerf_idx = oldPerf_idx_beg;
                                for (int perf = wells->well_connpos[ newIndex ];
                                     perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
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
                                for (int perf = wells->well_connpos[ newIndex ];
                                     perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
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
                    const WellControls* ctrl = wells->ctrls[w];
                    const int nwc = well_controls_get_num(ctrl);
                    int ctrl_index = 0;
                    for (; ctrl_index < nwc; ++ctrl_index) {
                        if (well_controls_iget_type(ctrl, ctrl_index) == THP) {
                            break;
                        }
                    }
                    // not finding any thp related control/limits
                    if (ctrl_index == nwc) {
                        thp()[w] = 0.;
                    }
                }
            }

            {
                // we need to create a trival segment related values to avoid there will be some
                // multi-segment wells added later.
                top_segment_index_.reserve(nw);
                for (int w = 0; w < nw; ++w) {
                    top_segment_index_.push_back(w);
                }
                segpress_ = bhp();
                segrates_ = wellRates();
            }
        }

        void resize(const Wells* wells, size_t numCells, const PhaseUsage& pu)
        {
            const std::vector<double> tmp(numCells, 0.0); // <- UGLY HACK to pass the size
            const std::vector<const Well*> wells_ecl;
            init(wells, tmp, wells_ecl, 0, nullptr, pu);
        }

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls.
        ///
        /// this method is only for flow_legacy!
        template <class PrevWellState>
        void initLegacy(const Wells* wells, const std::vector<double>& cellPressures , const PrevWellState& prevState, const PhaseUsage& pu)
        {
            // call init on base class
            BaseType :: init(wells, cellPressures);

            // if there are no well, do nothing in init
            if (wells == 0) {
                return;
            }

            const int nw = wells->number_of_wells;

            if( nw == 0 ) return ;

            // Initialize perfphaserates_, which must be done here.
            const int np = wells->number_of_phases;
            const int nperf = wells->well_connpos[nw];

            well_reservoir_rates_.resize(nw * np, 0.0);
            well_dissolved_gas_rates_.resize(nw, 0.0);
            well_vaporized_oil_rates_.resize(nw, 0.0);

            productivity_index_.resize(nw * np, 0.0);
            well_potentials_.resize(nw*np, 0.0);

            // Ensure that we start out with zero rates by default.
            perfphaserates_.clear();
            perfphaserates_.resize(nperf * np, 0.0);
            for (int w = 0; w < nw; ++w) {
                assert((wells->type[w] == INJECTOR) || (wells->type[w] == PRODUCER));
                const WellControls* ctrl = wells->ctrls[w];

                if (well_controls_well_is_stopped(ctrl)) {
                    // Shut well: perfphaserates_ are all zero.
                } else {
                    const int num_perf_this_well = wells->well_connpos[w + 1] - wells->well_connpos[w];
                    // Open well: Initialize perfphaserates_ to well
                    // rates divided by the number of perforations.
                    for (int perf = wells->well_connpos[w]; perf < wells->well_connpos[w + 1]; ++perf) {
                        for (int p = 0; p < np; ++p) {
                            perfphaserates_[np*perf + p] = wellRates()[np*w + p] / double(num_perf_this_well);
                        }
                        perfPress()[perf] = cellPressures[wells->well_cells[perf]];
                    }
                }
            }

            // Initialize current_controls_.
            // The controls set in the Wells object are treated as defaults,
            // and also used for initial values.
            current_controls_.resize(nw);
            for (int w = 0; w < nw; ++w) {
                current_controls_[w] = well_controls_get_current(wells->ctrls[w]);
            }

            perfRateSolvent_.clear();
            perfRateSolvent_.resize(nperf, 0.0);

            // intialize wells that have been there before
            // order may change so the mapping is based on the well name
            if( ! prevState.wellMap().empty() )
            {
                typedef typename WellMapType :: const_iterator const_iterator;
                const_iterator end = prevState.wellMap().end();
                for (int w = 0; w < nw; ++w) {
                    std::string name( wells->name[ w ] );
                    const_iterator it = prevState.wellMap().find( name );
                    if( it != end )
                    {
                        const int oldIndex = (*it).second[ 0 ];
                        const int newIndex = w;

                        // bhp
                        bhp()[ newIndex ] = prevState.bhp()[ oldIndex ];

                        // thp
                        thp()[ newIndex ] = prevState.thp()[ oldIndex ];

                        // wellrates
                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            wellRates()[ idx ] = prevState.wellRates()[ oldidx ];
                        }

                        // perfPhaseRates
                        const int oldPerf_idx_beg = (*it).second[ 1 ];
                        const int num_perf_old_well = (*it).second[ 2 ];
                        const int num_perf_this_well = wells->well_connpos[newIndex + 1] - wells->well_connpos[newIndex];
                        // copy perforation rates when the number of perforations is equal,
                        // otherwise initialize perfphaserates to well rates divided by the number of perforations.
                        if( num_perf_old_well == num_perf_this_well )
                        {
                            int old_perf_phase_idx = oldPerf_idx_beg *np;
                            for (int perf_phase_idx = wells->well_connpos[ newIndex ]*np;
                                 perf_phase_idx < wells->well_connpos[ newIndex + 1]*np; ++perf_phase_idx, ++old_perf_phase_idx )
                            {
                                perfPhaseRates()[ perf_phase_idx ] = prevState.perfPhaseRates()[ old_perf_phase_idx ];
                            }
                        } else {
                            for (int perf = wells->well_connpos[newIndex]; perf < wells->well_connpos[newIndex + 1]; ++perf) {
                                for (int p = 0; p < np; ++p) {
                                    perfPhaseRates()[np*perf + p] = wellRates()[np*newIndex + p] / double(num_perf_this_well);
                                }
                            }
                        }
                        // perfPressures
                        if( num_perf_old_well == num_perf_this_well )
                        {
                            int oldPerf_idx = oldPerf_idx_beg;
                            for (int perf = wells->well_connpos[ newIndex ];
                                 perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
                            {
                                perfPress()[ perf ] = prevState.perfPress()[ oldPerf_idx ];
                            }
                        }
                        // perfSolventRates
                        if (pu.has_solvent) {
                            if( num_perf_old_well == num_perf_this_well )
                            {
                                int oldPerf_idx = oldPerf_idx_beg;
                                for (int perf = wells->well_connpos[ newIndex ];
                                     perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
                                {
                                    perfRateSolvent()[ perf ] = prevState.perfRateSolvent()[ oldPerf_idx ];
                                }
                            }
                        }
                    }


                    // If in the new step, there is no THP related target/limit anymore, its thp value should be
                    // set to zero.
                    const WellControls* ctrl = wells->ctrls[w];
                    const int nwc = well_controls_get_num(ctrl);
                    int ctrl_index = 0;
                    for (; ctrl_index < nwc; ++ctrl_index) {
                        if (well_controls_iget_type(ctrl, ctrl_index) == THP) {
                            break;
                        }
                    }
                    // not finding any thp related control/limits
                    if (ctrl_index == nwc) {
                        thp()[w] = 0.;
                    }
                }
            }

            {
                // we need to create a trival segment related values to avoid there will be some
                // multi-segment wells added later.
                top_segment_index_.reserve(nw);
                for (int w = 0; w < nw; ++w) {
                    top_segment_index_.push_back(w);
                }
                segpress_ = bhp();
                segrates_ = wellRates();
            }
        }

        // this method is only for flow_legacy!
        template <class State, class PrevWellState>
        void initLegacy(const Wells* wells, const State& state, const PrevWellState& prevState, const PhaseUsage& pu)
        {
            initLegacy(wells, state.pressure(), prevState, pu);
        }

        // this method is only for flow_legacy!
        template <class State>
        void resizeLegacy(const Wells* wells, const State& state, const PhaseUsage& pu)
        {
            const WellStateFullyImplicitBlackoil dummy_state{}; // Init with an empty previous state only resizes
            initLegacy(wells, state, dummy_state, pu) ;
        }

        /// One rate per phase and well connection.
        std::vector<double>& perfPhaseRates() { return perfphaserates_; }
        const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

        /// One current control per well.
        std::vector<int>& currentControls() { return current_controls_; }
        const std::vector<int>& currentControls() const { return current_controls_; }

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

            if (pu.has_solvent) {
                // add solvent component
                for( int w = 0; w < nw; ++w ) {
                    res.at( wells_->name[ w ]).rates.set( rt::solvent, solventWellRate(w) );
                }
            }

            /* this is a reference or example on **how** to convert from
             * WellState to something understood by opm-output. it is intended
             * to be properly implemented and maintained as a part of
             * simulators, as it relies on simulator internals, details and
             * representations.
             */

            for( const auto& wt : this->wellMap() ) {
                const auto w = wt.second[ 0 ];
                auto& well = res.at( wt.first );
                well.control = this->currentControls()[ w ];

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

                well.rates.set( rt::dissolved_gas, this->well_dissolved_gas_rates_[w] );
                well.rates.set( rt::vaporized_oil, this->well_vaporized_oil_rates_[w] );

                int local_comp_index = 0;
                for( auto& comp : well.connections) {
                    const auto rates = this->perfPhaseRates().begin()
                                     + (np * wt.second[ 1 ])
                                     + (np * local_comp_index);
                    ++local_comp_index;

                    for( int i = 0; i < np; ++i ) {
                        comp.rates.set( phs[ i ], *(rates + i) );
                    }
                }
                assert(local_comp_index == this->wells_->well_connpos[ w + 1 ] - this->wells_->well_connpos[ w ]);
            }

            return res;
        }


        /// init the MS well related.
        template <typename PrevWellState>
        void initWellStateMSWell(const Wells* wells, const std::vector<const Well*>& wells_ecl,
                                 const int time_step, const PhaseUsage& pu, const PrevWellState& prev_well_state)
        {
            // still using the order in wells
            const int nw = wells->number_of_wells;
            if (nw == 0) {
                return;
            }

            top_segment_index_.clear();
            top_segment_index_.reserve(nw);
            segpress_.clear();
            segpress_.reserve(nw);
            segrates_.clear();
            segrates_.reserve(nw * numPhases());

            nseg_ = 0;
            // in the init function, the well rates and perforation rates have been initialized or copied from prevState
            // what we do here, is to set the segment rates and perforation rates
            for (int w = 0; w < nw; ++w) {
                const int nw_wells_ecl = wells_ecl.size();
                int index_well_ecl = 0;
                const std::string well_name(wells->name[w]);
                for (; index_well_ecl < nw_wells_ecl; ++index_well_ecl) {
                    if (well_name == wells_ecl[index_well_ecl]->name()) {
                        break;
                    }
                }

                // It should be able to find in wells_ecl.
                if (index_well_ecl == nw_wells_ecl) {
                    OPM_THROW(std::logic_error, "Could not find well " << well_name << " in wells_ecl ");
                }

                const Well* well_ecl = wells_ecl[index_well_ecl];
                top_segment_index_.push_back(nseg_);
                if ( !well_ecl->isMultiSegment(time_step) ) { // not multi-segment well
                    nseg_ += 1;
                    segpress_.push_back(bhp()[w]);
                    const int np = numPhases();
                    for (int p = 0; p < np; ++p) {
                        segrates_.push_back(wellRates()[np * w + p]);
                    }
                } else { // it is a multi-segment well
                    const WellSegments& segment_set = well_ecl->getWellSegments(time_step);
                    // assuming the order of the perforations in well_ecl is the same with Wells
                    const WellConnections& completion_set = well_ecl->getConnections(time_step);
                    // number of segment for this single well
                    const int well_nseg = segment_set.size();
                    const int nperf = completion_set.size();
                    nseg_ += well_nseg;
                    // we need to know for each segment, how many perforation it has and how many segments using it as outlet_segment
                    // that is why I think we should use a well model to initialize the WellState here
                    std::vector<std::vector<int>> segment_perforations(well_nseg);
                    for (int perf = 0; perf < nperf; ++perf) {
                        const Connection& connection = completion_set.get(perf);
                        const int segment_index = segment_set.segmentNumberToIndex(connection.segment());
                        segment_perforations[segment_index].push_back(perf);
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


                    // for the segrates_, now it becomes a recursive solution procedure.
                    {
                        const int np = numPhases();
                        const int start_perf = wells->well_connpos[w];
                        const int start_perf_next_well = wells->well_connpos[w + 1];
                        assert(nperf == (start_perf_next_well - start_perf)); // make sure the information from wells_ecl consistent with wells
                        if (pu.phase_used[Gas]) {
                            const int gaspos = pu.phase_pos[Gas];
                            // scale the phase rates for Gas to avoid too bad initial guess for gas fraction
                            // it will probably benefit the standard well too, while it needs to be justified
                            // TODO: to see if this strategy can benefit StandardWell too
                            // TODO: it might cause big problem for gas rate control or if there is a gas rate limit
                            // maybe the best way is to initialize the fractions first then get the rates
                            for (int perf = 0; perf < nperf; perf++) {
                                const int perf_pos = start_perf + perf;
                                perfPhaseRates()[np * perf_pos + gaspos] *= 100.;
                            }
                        }

                        const std::vector<double> perforation_rates(perfPhaseRates().begin() + np * start_perf,
                                                                    perfPhaseRates().begin() + np * start_perf_next_well); // the perforation rates for this well
                        std::vector<double> segment_rates;
                        calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, 0 /* top segment */, segment_rates);
                        std::copy(segment_rates.begin(), segment_rates.end(), std::back_inserter(segrates_));
                    }

                    // for the segment pressure, the segment pressure is the same with the first perforation belongs to the segment
                    // if there is no perforation associated with this segment, it uses the pressure from the outlet segment
                    // which requres the ordering is successful
                    // Not sure what is the best way to handle the initialization, hopefully, the bad initialization can be
                    // improved during the solveWellEq process
                    {
                        // top segment is always the first one, and its pressure is the well bhp
                        segpress_.push_back(bhp()[w]);
                        const int top_segment = top_segment_index_[w];
                        const int start_perf = wells->well_connpos[w];
                        for (int seg = 1; seg < well_nseg; ++seg) {
                            if ( !segment_perforations[seg].empty() ) {
                                const int first_perf = segment_perforations[seg][0];
                                segpress_.push_back(perfPress()[start_perf + first_perf]);
                            } else {
                                // segpress_.push_back(bhp); // may not be a good decision
                                // using the outlet segment pressure // it needs the ordering is correct
                                const int outlet_seg = segment_set[seg].outletSegment();
                                segpress_.push_back(segpress_[top_segment + segment_set.segmentNumberToIndex(outlet_seg)]);
                            }
                        }
                    }
                }
            }
            assert(int(segpress_.size()) == nseg_);
            assert(int(segrates_.size()) == nseg_ * numPhases() );

            if (!prev_well_state.wellMap().empty()) {
                // copying MS well related
                const auto& end = prev_well_state.wellMap().end();
                const int np = numPhases();
                for (int w = 0; w < nw; ++w) {
                    const std::string name( wells->name[w] );
                    const auto& it = prev_well_state.wellMap().find( name );

                    if (it != end) { // the well is found in the prev_well_state
                        // TODO: the well with same name can change a lot, like they might not have same number of segments
                        // we need to handle that later.
                        // for now, we just copy them.
                        const int old_index_well = (*it).second[0];
                        const int new_index_well = w;
                        const int old_top_segment_index = prev_well_state.topSegmentIndex(old_index_well);
                        const int new_top_segmnet_index = topSegmentIndex(new_index_well);
                        int number_of_segment = 0;
                        // if it is the last well in list
                        if (new_index_well == int(top_segment_index_.size()) - 1) {
                            number_of_segment = nseg_ - new_top_segmnet_index;
                        } else {
                            number_of_segment = topSegmentIndex(new_index_well + 1) - new_top_segmnet_index;
                        }

                        for (int i = 0; i < number_of_segment * np; ++i) {
                            segrates_[new_top_segmnet_index * np + i] = prev_well_state.segRates()[old_top_segment_index * np + i];
                        }

                        for (int i = 0; i < number_of_segment; ++i) {
                            segpress_[new_top_segmnet_index + i] = prev_well_state.segPress()[old_top_segment_index + i];
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


        /// One rate pr well connection.
        std::vector<double>& perfRateSolvent() { return perfRateSolvent_; }
        const std::vector<double>& perfRateSolvent() const { return perfRateSolvent_; }

        /// One rate pr well
        double solventWellRate(const int w) const {
            double solvent_well_rate = 0.0;
            for (int perf = wells_->well_connpos[w]; perf < wells_->well_connpos[w+1]; ++perf ) {
                solvent_well_rate += perfRateSolvent_[perf];
            }
            return solvent_well_rate;
        }

        std::vector<double>& wellReservoirRates()
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
            return segrates_;
        }

        std::vector<double>& segRates()
        {
            return segrates_;
        }

        const std::vector<double>& segPress() const
        {
            return segpress_;
        }

        std::vector<double>& segPress()
        {
            return segpress_;
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

    private:
        std::vector<double> perfphaserates_;
        std::vector<int> current_controls_;
        std::vector<double> perfRateSolvent_;

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
        std::vector<double> segrates_;
        std::vector<double> segpress_;
        // the index of the top segments, which is used to locate the
        // multisegment well related information in WellState
        std::vector<int> top_segment_index_;
        int nseg_; // total number of the segments

        // Productivity Index
        std::vector<double> productivity_index_;

        // Well potentials
        std::vector<double> well_potentials_;

    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
