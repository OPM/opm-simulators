/*
  Copyright 2016 Statoil ASA
  Copyright 2016 IRIS
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.

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

#include <opm/autodiff/Compat.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>

#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/common/data/SimulationDataContainer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>

namespace Opm {

std::vector< double > destripe( const std::vector< double >& v,
                                size_t stride,
                                size_t offset ) {

    std::vector< double > dst( v.size() / stride );

    size_t di = 0;
    for( size_t i = offset; i < v.size(); i += stride ) {
        dst[ di++ ] = v[ i ];
    }

    return dst;
}






std::vector< double >& stripe( const std::vector< double >& v,
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






data::Solution simToSolution( const SimulationDataContainer& reservoir,
                              const bool use_si_units,
                              PhaseUsage phases ) {

    // Set up unit system to use to suppress conversion if use_si_units is true.
    const UnitSystem::measure press_unit = use_si_units ? UnitSystem::measure::identity : UnitSystem::measure::pressure;
    const UnitSystem::measure temp_unit = use_si_units ? UnitSystem::measure::identity : UnitSystem::measure::temperature;
    const UnitSystem::measure rs_unit = use_si_units ? UnitSystem::measure::identity : UnitSystem::measure::gas_oil_ratio;
    const UnitSystem::measure rv_unit = use_si_units ? UnitSystem::measure::identity : UnitSystem::measure::oil_gas_ratio;

    data::Solution sol;
    sol.insert( "PRESSURE", press_unit, reservoir.pressure() , data::TargetType::RESTART_SOLUTION);
    sol.insert( "TEMP"    , temp_unit, reservoir.temperature() , data::TargetType::RESTART_SOLUTION );

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
        sol.insert( "RS", rs_unit, reservoir.getCellData( BlackoilState::GASOILRATIO ) , data::TargetType::RESTART_SOLUTION );
    }

    if( reservoir.hasCellData( BlackoilState::RV ) ) {
        sol.insert( "RV", rv_unit, reservoir.getCellData( BlackoilState::RV ) , data::TargetType::RESTART_SOLUTION );
    }

    if (phases.has_solvent) {
        sol.insert( "SSOL", UnitSystem::measure::identity, reservoir.getCellData( BlackoilState::SSOL ) , data::TargetType::RESTART_SOLUTION );
    }

    if (phases.has_polymer) {
        if (reservoir.hasCellData( PolymerBlackoilState::CONCENTRATION )) { // compatibility with legacy polymer
            sol.insert( "POLYMER", UnitSystem::measure::identity, reservoir.getCellData( PolymerBlackoilState::CONCENTRATION ) , data::TargetType::RESTART_SOLUTION );
        } else {
            sol.insert( "POLYMER", UnitSystem::measure::identity, reservoir.getCellData( BlackoilState::POLYMER ) , data::TargetType::RESTART_SOLUTION );
        }

    }

    return sol;
}






void solutionToSim( const RestartValue& restart_value,
                    PhaseUsage phases,
                    SimulationDataContainer& state ) {

    const auto& sol = restart_value.solution;
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

    for (size_t c = 0; c < state.numCells(); ++c) {
        double& so = state.saturation()[phases.num_phases*c + phases.phase_pos[ BlackoilPhases::Liquid ]];
        so = 1.0;
        if (phases.phase_used[ BlackoilPhases::Aqua]) {
            so -= state.saturation()[phases.num_phases*c + phases.phase_pos[ BlackoilPhases::Aqua ]];
        }
        if (phases.phase_used[ BlackoilPhases::Vapour]) {
            so -= state.saturation()[phases.num_phases*c + phases.phase_pos[ BlackoilPhases::Vapour ]];
        }
    }


    if( sol.has( "PRESSURE" ) ) {
        state.pressure() = sol.data( "PRESSURE" );
    }

    if( sol.has( "TEMP" ) ) {
        state.temperature() = sol.data( "TEMP" );
    }

    if( sol.has( "RS" ) ) {
        state.registerCellData("GASOILRATIO", 1);
        state.getCellData( "GASOILRATIO" ) = sol.data( "RS" );
    }

    if( sol.has( "RV" ) ) {
        state.registerCellData("RV", 1);
        state.getCellData( "RV" ) = sol.data( "RV" );
    }

    if ( sol.has( "SSOL" ) ) {
        state.registerCellData("SSOL", 1);
        state.getCellData("SSOL") = sol.data("SSOL");
    }

    if ( sol.has("SOMAX" ) ) {
        state.registerCellData("SOMAX", 1);
        state.getCellData("SOMAX") = sol.data("SOMAX");
    }

    if ( sol.has("PCSWM_OW" ) ) {
        state.registerCellData("PCSWMDC_OW", 1);
        state.getCellData("PCSWMDC_OW") = sol.data("PCSWM_OW");
    }

    if ( sol.has("KRNSW_OW" ) ) {
        state.registerCellData("KRNSWMDC_OW", 1);
        state.getCellData("KRNSWMDC_OW") = sol.data("KRNSW_OW");
    }

    if ( sol.has("PCSWM_GO" ) ) {
        state.registerCellData("PCSWMDC_GO", 1);
        state.getCellData("PCSWMDC_GO") = sol.data("PCSWM_GO");
    }

    if ( sol.has("KRNSW_GO" ) ) {
        state.registerCellData("KRNSWMDC_GO", 1);
        state.getCellData("KRNSWMDC_GO") = sol.data("KRNSW_GO");
    }
}






void wellsToState( const data::Wells& wells,
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

} // namespace Opm
