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

#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
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

        using BaseType :: wellRates;
        using BaseType :: bhp;
        using BaseType :: perfPress;
        using BaseType :: wellMap;
        using BaseType :: numWells;
        using BaseType :: numPhases;

        template <class State, class PrevWellState>
        void init(const Wells* wells, const State& state, const PrevWellState& prevState, const PhaseUsage& pu)
        {
            init(wells, state.pressure(), prevState, pu);
        }

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        template <class PrevWellState>
        void init(const Wells* wells, const std::vector<double>& cellPressures , const PrevWellState& prevState, const PhaseUsage& pu)
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

            is_new_well_.resize(nw, true);

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
                        // this is not a new added well
                        is_new_well_[w] = false;

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
        }

        template <class State>
        void resize(const Wells* wells, const State& state, const PhaseUsage& pu) {
            const WellStateFullyImplicitBlackoil dummy_state{}; // Init with an empty previous state only resizes
            init(wells, state, dummy_state, pu) ;
        }

        /// One rate per phase and well connection.
        std::vector<double>& perfPhaseRates() { return perfphaserates_; }
        const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

        /// One current control per well.
        std::vector<int>& currentControls() { return current_controls_; }
        const std::vector<int>& currentControls() const { return current_controls_; }

        data::Wells report(const PhaseUsage &pu) const override {
            data::Wells res = WellState::report(pu);

            const int nw = this->numWells();
            if( nw == 0 ) return res;
            const int np = pu.num_phases;


            using rt = data::Rates::opt;
            std::vector< rt > phs( np );
            if( pu.phase_used[BlackoilPhases::Aqua] ) {
                phs.at( pu.phase_pos[BlackoilPhases::Aqua] ) = rt::wat;
            }

            if( pu.phase_used[BlackoilPhases::Liquid] ) {
                phs.at( pu.phase_pos[BlackoilPhases::Liquid] ) = rt::oil;
            }

            if( pu.phase_used[BlackoilPhases::Vapour] ) {
                phs.at( pu.phase_pos[BlackoilPhases::Vapour] ) = rt::gas;
            }

            if (pu.has_solvent) {
                // add solvent component
                for( int w = 0; w < nw; ++w ) {
                    using rt = data::Rates::opt;
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

                int local_comp_index = 0;
                for( auto& comp : well.completions ) {
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


        bool isNewWell(const int w) const {
            return is_new_well_[w];
        }


        void setNewWell(const int w, const bool is_new_well) {
            is_new_well_[w] = is_new_well;
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

    private:
        std::vector<double> perfphaserates_;
        std::vector<int> current_controls_;
        std::vector<double> perfRateSolvent_;

        // marking whether the well is just added
        // for newly added well, the current initialized rates from WellState
        // will have very wrong compositions for production wells, will mostly cause
        // problem with VFP interpolation
        std::vector<bool> is_new_well_;
    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
