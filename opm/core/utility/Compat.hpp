/*
   Copyright 2016 Statoil ASA

   This file is part of the Open Porous Media project (OPM).

   OPM is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OPM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with OPM. If not, see <http://www.gnu.org/licenses/>.
   */

#ifndef OPM_CORE_COMPAT_HPP
#define OPM_CORE_COMPAT_HPP

#include <opm/common/data/SimulationDataContainer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/output/Cells.hpp>
#include <opm/output/Wells.hpp>

namespace Opm {

inline std::vector< double > destripe( const std::vector< double >& v,
                                       size_t stride,
                                       size_t offset ) {

    std::vector< double > dst( v.size() / stride );

    size_t di = 0;
    for( size_t i = offset; i < v.size(); i += stride ) {
        dst[ di++ ] = v[ i ];
    }

    return dst;
}





inline std::vector< double >& stripe( const std::vector< double >& v,
                                      size_t stride,
                                      size_t offset,
                                      std::vector< double >& dst ) {

    /* does little range checking etc; for future revisions */
    size_t vi = 0;
    for( size_t i = offset; i < dst.size(); i += stride ) {
        dst[ i ] = v[ vi++ ];
    }

    return dst;
}







inline data::Solution simToSolution( const SimulationDataContainer& reservoir,
                             PhaseUsage phases ) {
    using ds = data::Solution::key;

    data::Solution sol;
    sol.insert( ds::PRESSURE, reservoir.pressure() );
    sol.insert( ds::TEMP, reservoir.temperature() );

    const auto ph = reservoir.numPhases();
    const auto& sat = reservoir.saturation();

    const auto aqua = BlackoilPhases::Aqua;
    const auto vapour = BlackoilPhases::Vapour;

    if( phases.phase_used[ aqua ] ) {
        sol.insert( ds::SWAT, destripe( sat, ph, phases.phase_pos[ aqua ] ) );
    }

    if( phases.phase_used[ vapour ] ) {
        sol.insert( ds::SGAS, destripe( sat, ph, phases.phase_pos[ vapour ] ) );
    }

    if( reservoir.hasCellData( BlackoilState::GASOILRATIO ) ) {
        sol.insert( ds::RS, reservoir.getCellData( BlackoilState::GASOILRATIO ) );
    }

    if( reservoir.hasCellData( BlackoilState::RV ) ) {
        sol.insert( ds::RV, reservoir.getCellData( BlackoilState::RV ) );
    }

    sol.sdc = &reservoir;

    return sol;
}








inline void solutionToSim( const data::Solution& sol,
                          PhaseUsage phases,
                          SimulationDataContainer& state ) {
    using ds = data::Solution::key;

    const auto stride = phases.num_phases;
    if( sol.has( ds::SWAT ) ) {
        stripe( sol[ ds::SWAT ],
                stride,
                phases.phase_pos[ BlackoilPhases::Aqua ],
                state.saturation() );
    }

    if( sol.has( ds::SGAS ) ) {
        stripe( sol[ ds::SGAS ],
                stride,
                phases.phase_pos[ BlackoilPhases::Vapour ],
                state.saturation() );
    }

    if( sol.has( ds::PRESSURE ) ) {
        state.pressure() = sol[ ds::PRESSURE ];
    }

    if( sol.has( ds::TEMP ) ) {
        state.temperature() = sol[ ds::TEMP ];
    }

    if( sol.has( ds::RS ) ) {
        state.getCellData( "GASOILRATIO" ) = sol[ ds::RS ];
    }

    if( sol.has( ds::RV ) ) {
        state.getCellData( "RV" ) = sol[ ds::RV ];
    }
}















inline void wellsToState( const data::Wells& wells, WellState& state ) {
    state.bhp() = wells.bhp;
    state.temperature() = wells.temperature;
    state.wellRates() = wells.well_rate;
    state.perfPress() = wells.perf_pressure;
    state.perfRates() = wells.perf_rate;
}

}

#endif //OPM_CORE_COMPAT_HPP
