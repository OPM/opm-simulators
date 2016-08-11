/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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
#include <opm/autodiff/BlackoilModelEnums.hpp>
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

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        template <class State, class PrevState>
        void init(const Wells* wells, const State& state, const PrevState& prevState)
        {
            // call init on base class
            BaseType :: init(wells, state);

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
                        perfPress()[perf] = state.pressure()[wells->well_cells[perf]];
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

            well_potentials_.clear();
            well_potentials_.resize(nperf * np, 0.0);

            well_solutions_.clear();
            well_solutions_.resize(nw * np, 0.0);
            std::vector<double> g = {1.0,1.0,0.01};
            for (int w = 0; w < nw; ++w) {
                const WellControls* wc = wells->ctrls[w];

                // The current control in the well state overrides
                // the current control set in the Wells struct, which
                // is instead treated as a default.
                const int current = current_controls_[w];
                const WellType& well_type = wells->type[w];

                switch (well_controls_iget_type(wc, current)) {
                case BHP:
                {
                    if (well_type == INJECTOR) {
                        for (int p = 0; p < np; ++p)  {
                            well_solutions_[w] += wellRates()[np*w + p] * wells->comp_frac[np*w + p];
                        }
                    } else {
                        for (int p = 0; p < np; ++p) {
                            well_solutions_[w] += g[p] * wellRates()[np*w + p];
                        }
                    }
                }
                    break;


                case RESERVOIR_RATE: // Intentional fall-through
                case SURFACE_RATE:
                {
                    wellSolutions()[w] = bhp()[w];

                }
                    break;
                }
                assert(np == 3);
                double total_rates = 0.0;
                for (int p = 0; p < np; ++p)  {
                    total_rates += g[p] * wellRates()[np*w + p];
                }

                //if(std::abs(total_rates) > 0) {
                //    wellSolutions()[nw + w] = g[Water] * wellRates()[np*w + Water] / total_rates; //wells->comp_frac[np*w + Water]; // Water;
                //    wellSolutions()[2*nw + w] = g[Gas] * wellRates()[np*w + Gas] / total_rates ; //wells->comp_frac[np*w + Gas]; //Gas
                //} else {
                wellSolutions()[nw + w] = wells->comp_frac[np*w + Water];
                wellSolutions()[2*nw + w] = wells->comp_frac[np*w + Gas];
                //}


            }


            // intialize wells that have been there before
            // order may change so the mapping is based on the well name
            if( ! prevState.wellMap().empty() )
            {
                typedef typename WellMapType :: const_iterator const_iterator;
                const_iterator end = prevState.wellMap().end();
                int nw_old = prevState.bhp().size();
                for (int w = 0; w < nw; ++w) {
                    std::string name( wells->name[ w ] );
                    const_iterator it = prevState.wellMap().find( name );
                    if( it != end )
                    {
                        const int oldIndex = (*it).second[ 0 ];
                        const int newIndex = w;

                        // bhp
                        bhp()[ newIndex ] = prevState.bhp()[ oldIndex ];

                        // wellrates
                        double total_well_rates = 0.0;
                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            total_well_rates += prevState.wellRates()[ oldidx ];
                        }

                        //if (std::abs(total_well_rates) > 0) {

                        for( int i=0, idx=newIndex*np, oldidx=oldIndex*np; i<np; ++i, ++idx, ++oldidx )
                        {
                            wellRates()[ idx ] = prevState.wellRates()[ oldidx ];
                        }

                        // wellSolutions

                        if (wells->type[w] == PRODUCER && std::abs(total_well_rates) > 0.0) {
                            for( int i=0;  i<np; ++i)
                            {
                                wellSolutions()[ i*nw + newIndex ] = prevState.wellSolutions()[i * nw_old + oldIndex ];
                            }
                        }




                        // perfPhaseRates
                        int oldPerf_idx = (*it).second[ 1 ];
                        const int num_perf_old_well = (*it).second[ 2 ];
                        const int num_perf_this_well = wells->well_connpos[newIndex + 1] - wells->well_connpos[newIndex];
                        // copy perforation rates when the number of perforations is equal,
                        // otherwise initialize perfphaserates to well rates divided by the number of perforations.
                        if( num_perf_old_well == num_perf_this_well )
                        {
                            int oldPerf = oldPerf_idx *np;
                            for (int perf = wells->well_connpos[ newIndex ]*np;
                                 perf < wells->well_connpos[ newIndex + 1]*np; ++perf, ++oldPerf )
                            {
                                perfPhaseRates()[ perf ] = prevState.perfPhaseRates()[ oldPerf ];
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
                            for (int perf = wells->well_connpos[ newIndex ];
                                 perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
                            {
                                perfPress()[ perf ] = prevState.perfPress()[ oldPerf_idx ];
                            }
                        }

                        // currentControls
                        const int old_control_index = prevState.currentControls()[ oldIndex ];
                        if (old_control_index < well_controls_get_num(wells->ctrls[w])) {
                            // If the set of controls have changed, this may not be identical
                            // to the last control, but it must be a valid control.
                            currentControls()[ newIndex ] = old_control_index;
                        }

                    }
                }
            }
        }

        template <class State>
        void resize(const Wells* wells, const State& state) {
            const WellStateFullyImplicitBlackoil dummy_state{}; // Init with an empty previous state only resizes
            init(wells, state, dummy_state) ;
        }


        /// One rate per phase and well connection.
        std::vector<double>& perfPhaseRates() { return perfphaserates_; }
        const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

        /// One current control per well.
        std::vector<int>& currentControls() { return current_controls_; }
        const std::vector<int>& currentControls() const { return current_controls_; }

        /// One rate per phase and well connection.
        std::vector<double>& wellPotentials() { return well_potentials_; }
        const std::vector<double>& wellPotentials() const { return well_potentials_; }

        /// One rate per phase and well connection.
        std::vector<double>& wellSolutions() { return well_solutions_; }
        const std::vector<double>& wellSolutions() const { return well_solutions_; }

        data::Wells report() const override {
            data::Wells res = WellState::report();

            const int nw = this->numWells();
            // If there are now wells numPhases throws a floating point
            // exception.
            const int np = nw ? this->numPhases() : -1;

            /* this is a reference or example on **how** to convert from
             * WellState to something understood by opm-output. it is intended
             * to be properly implemented and maintained as a part of
             * simulators, as it relies on simulator internals, details and
             * representations.
             */

            for( auto w = 0; w < nw; ++w ) {
                using rt = data::Rates::opt;
                std::map< size_t, data::Completion > completions;

                // completions aren't supported yet
                //const auto* begin = wells_->well_connpos + w;
                //const auto* end = wells_->well_connpos + w + 1;
                //for( auto* i = begin; i != end; ++i ) {
                //    const auto perfrate = this->perfPhaseRates().begin() + *i;
                //    data::Rates perfrates;
                //    perfrates.set( rt::wat, *(perfrate + 0) );
                //    perfrates.set( rt::oil, *(perfrate + 1) );
                //    perfrates.set( rt::gas, *(perfrate + 2) );

                //    const size_t active_index = wells_->well_cells[ *i ];

                //    completions.emplace( active_index,
                //        data::Completion{ active_index, perfrates } );
                //}

                const auto wellrate_index = np * w;
                const auto& wv = this->wellRates();

                data::Rates wellrates;
                if( np == 3 ) {
                    /* only write if 3-phase solution */
                    wellrates.set( rt::wat, wv[ wellrate_index + 0 ] );
                    wellrates.set( rt::oil, wv[ wellrate_index + 1 ] );
                    wellrates.set( rt::gas, wv[ wellrate_index + 2 ] );
                }

                const double bhp  = this->bhp()[ w ];
                const double thp  = this->thp()[ w ];

                res.emplace( wells_->name[ w ],
                    data::Well { wellrates, bhp, thp, std::move( completions ) } );
            }

            return res;
        }

    private:
        std::vector<double> perfphaserates_;
        std::vector<int> current_controls_;
        std::vector<double> well_potentials_;
        std::vector<double> well_solutions_;
    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
