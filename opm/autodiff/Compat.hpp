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

#include <algorithm>
#include <cassert>

#include <opm/common/data/SimulationDataContainer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
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













inline void wellsToState( const data::Wells& wells,
                          PhaseUsage phases,
                          WellStateFullyImplicitBlackoil& state ) {

    using rt = data::Rates::opt;

    const auto np = phases.num_phases;

    std::vector< rt > phs( np );
    if( phases.phase_used[BlackoilPhases::Aqua] ) {
        phs.at( phases.phase_pos[BlackoilPhases::Aqua] ) = rt::wat;
    }

    if( phases.phase_used[BlackoilPhases::Liquid] ) {
        phs.at( phases.phase_pos[BlackoilPhases::Liquid] ) = rt::oil;
    }

    if( phases.phase_used[BlackoilPhases::Vapour] ) {
        phs.at( phases.phase_pos[BlackoilPhases::Vapour] ) = rt::gas;
    }

    for( const auto& wm : state.wellMap() ) {
        const auto well_index = wm.second[ 0 ];
        const auto& well = wells.at( wm.first );

        state.bhp()[ well_index ] = well.bhp;
        state.temperature()[ well_index ] = well.temperature;
        state.currentControls()[ well_index ] = well.control;

        const auto wellrate_index = well_index * np;
        for( size_t i = 0; i < phs.size(); ++i ) {
            assert( well.rates.has( phs[ i ] ) );
            state.wellRates()[ wellrate_index + i ] = well.rates.get( phs[ i ] );
        }

        const auto perforation_pressure = []( const data::Completion& comp ) {
            return comp.pressure;
        };

        const auto perforation_reservoir_rate = []( const data::Completion& comp ) {
            return comp.reservoir_rate;
        };

        std::transform( well.completions.begin(),
                        well.completions.end(),
                        state.perfPress().begin() + wm.second[ 1 ],
                        perforation_pressure );

        std::transform( well.completions.begin(),
                        well.completions.end(),
                        state.perfRates().begin() + wm.second[ 1 ],
                        perforation_reservoir_rate );

        int local_comp_index = 0;
        for (const data::Completion& comp : well.completions) {
            const int global_comp_index = wm.second[1] + local_comp_index;
            for (int phase_index = 0; phase_index < np; ++phase_index) {
                state.perfPhaseRates()[global_comp_index*np + phase_index] = comp.rates.get(phs[phase_index]);
            }
            ++local_comp_index;
        }
    }
}

}

#endif //OPM_SIMULATORS_COMPAT_HPP
