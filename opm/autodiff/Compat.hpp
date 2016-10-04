/*
   Copyright 2016 Statoil ASA
   2016 IRIS

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

#ifndef OPM_SIMULATORS_COMPAT_HPP
#define OPM_SIMULATORS_COMPAT_HPP

#include <opm/common/data/SimulationDataContainer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/autodiff/BlackoilSolventState.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>

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
    data::Solution sol;
    sol.insert( "PRESSURE", UnitSystem::measure::pressure, reservoir.pressure() , data::TargetType::RESTART_SOLUTION);
    sol.insert( "TEMP"    , UnitSystem::measure::temperature, reservoir.temperature() , data::TargetType::RESTART_SOLUTION );

    const auto ph = reservoir.numPhases();
    const auto& sat = reservoir.saturation();

    const auto aqua = BlackoilPhases::Aqua;
    const auto vapour = BlackoilPhases::Vapour;

    if( phases.phase_used[ aqua ] ) {
        sol.insert( "SWAT", UnitSystem::measure::identity, destripe( sat, ph, phases.phase_pos[ aqua ] ) , data::TargetType::RESTART_SOLUTION );
    }

    if( phases.phase_used[ vapour ] ) {
        sol.insert( "SGAS", UnitSystem::measure::identity, destripe( sat, ph, phases.phase_pos[ vapour ] ) , data::TargetType::RESTART_SOLUTION );
    }

    if( reservoir.hasCellData( BlackoilState::GASOILRATIO ) ) {
        sol.insert( "RS", UnitSystem::measure::identity, reservoir.getCellData( BlackoilState::GASOILRATIO ) , data::TargetType::RESTART_SOLUTION );
    }

    if( reservoir.hasCellData( BlackoilState::RV ) ) {
        sol.insert( "RV", UnitSystem::measure::identity, reservoir.getCellData( BlackoilState::RV ) , data::TargetType::RESTART_SOLUTION );
    }

    if (reservoir.hasCellData( BlackoilSolventState::SSOL)) {
        sol.insert( "SSOL", UnitSystem::measure::identity, reservoir.getCellData( BlackoilSolventState::SSOL ) , data::TargetType::RESTART_SOLUTION );
    }

    sol.sdc = &reservoir;

    return sol;
}








inline void solutionToSim( const data::Solution& sol,
                          PhaseUsage phases,
                          SimulationDataContainer& state ) {

    const auto stride = phases.num_phases;
    if( sol.has( "SWAT" ) ) {
        stripe( sol.data( "SWAT" ),
                stride,
                phases.phase_pos[ BlackoilPhases::Aqua ],
                state.saturation() );
    }

    if( sol.has( "SGAS" ) ) {
        stripe( sol.data( "SGAS" ),
                stride,
                phases.phase_pos[ BlackoilPhases::Vapour ],
                state.saturation() );
    }

    if( sol.has( "PRESSURE" ) ) {
        state.pressure() = sol.data( "PRESSURE" );
    }

    if( sol.has( "TEMP" ) ) {
        state.temperature() = sol.data( "TEMP" );
    }

    if( sol.has( "RS" ) ) {
        state.getCellData( "GASOILRATIO" ) = sol.data( "RS" );
    }

    if( sol.has( "RV" ) ) {
        state.getCellData( "RV" ) = sol.data( "RV" );
    }

    if (sol.has( "SSOL" ) ) {
        state.getCellData("SSOL") = sol.data("SSOL");
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

#endif //OPM_SIMULATORS_COMPAT_HPP
