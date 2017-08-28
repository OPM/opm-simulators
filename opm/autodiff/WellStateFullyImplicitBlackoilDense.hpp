/*
  Copyright 2016 IRIS

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

#ifndef OPM_WELLSTATEFULLYIMPLICITBLACKOILDENSE_HEADER_INCLUDED
#define OPM_WELLSTATEFULLYIMPLICITBLACKOILDENSE_HEADER_INCLUDED


#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <vector>
#include <cassert>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <array>
#include <cmath>

namespace Opm
{

    /// The state of a set of wells, tailored for use by the fully
    /// implicit blackoil simulator.
    class WellStateFullyImplicitBlackoilDense
        : public WellStateFullyImplicitBlackoil
    {
        typedef WellStateFullyImplicitBlackoil  BaseType;
    public:
        typedef BaseType :: WellMapType WellMapType;

        using BaseType :: wellRates;
        using BaseType :: bhp;
        using BaseType :: perfPress;
        using BaseType :: wellMap;
        using BaseType :: numWells;
        using BaseType :: numPhases;
        using BaseType :: perfPhaseRates;
        using BaseType :: currentControls;

        /// Allocate and initialize if wells is non-null.  Also tries
        /// to give useful initial values to the bhp(), wellRates()
        /// and perfPhaseRates() fields, depending on controls
        template <class PrevWellState>
        void init(const Wells* wells, const std::vector<double>& cellPressures, const PrevWellState& prevState, const PhaseUsage& pu)
        {

            // call init on base class
            BaseType :: init(wells, cellPressures, prevState);


            const int nw = wells->number_of_wells;
            if (nw == 0) {
                return;
            }
            const int nperf = wells->well_connpos[nw];
            perfRateSolvent_.clear();
            perfRateSolvent_.resize(nperf, 0.0);

            if (pu.has_solvent) {

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
                            const int newIndex = w;

                            // perfSolventRates
                            int oldPerf_idx = (*it).second[ 1 ];
                            const int num_perf_old_well = (*it).second[ 2 ];
                            const int num_perf_this_well = wells->well_connpos[newIndex + 1] - wells->well_connpos[newIndex];
                            if( num_perf_old_well == num_perf_this_well )
                            {
                                for (int perf = wells->well_connpos[ newIndex ];
                                     perf < wells->well_connpos[ newIndex + 1]; ++perf, ++oldPerf_idx )
                                {
                                    perfRateSolvent()[ perf ] = prevState.perfRateSolvent()[ oldPerf_idx ];
                                }
                            }
                        }
                    }
                }
            }
        }


        template <class State>
        void resize(const Wells* wells, const State& state, const PhaseUsage& pu ) {
            const WellStateFullyImplicitBlackoilDense dummy_state{}; // Init with an empty previous state only resizes            
            init(wells, state.pressure(), dummy_state, pu) ;
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


        data::Wells report(const PhaseUsage& pu) const override {
            data::Wells res = BaseType::report(pu);
            const int nw = WellState::numWells();
            // If there are now wells numPhases throws a floating point
            // exception.
            if (nw == 0) {
                return res;
            }
            if (pu.has_solvent) {
                // add solvent component
                for( int w = 0; w < nw; ++w ) {
                    using rt = data::Rates::opt;
                    res.at( wells_->name[ w ]).rates.set( rt::solvent, solventWellRate(w) );
                }
            }

            return res;
        }


    private:
        std::vector<double> perfRateSolvent_;

    };

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOILDENSE_HEADER_INCLUDED
