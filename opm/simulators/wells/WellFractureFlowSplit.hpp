/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_WELL_FRACTURE_FLOW_SPLIT_HPP_INCLUDED
#define OPM_WELL_FRACTURE_FLOW_SPLIT_HPP_INCLUDED

#include <opm/material/common/MathToolbox.hpp>

#include <cassert>

namespace Opm {

/// Split each injector connection's water flow between the matrix and a
/// fracture crossing the connection, based on the ratio of the regular
/// connection transmissibility to the fracture contribution registered via
/// WellInterfaceGeneric::addFracturePerforations().
///
/// Fills, per connection, the filtrate flow factor (matrix share) and the
/// fracture water rate, and accumulates the total fracture-routed rate of
/// the well (SingleWellState::frac_rate).  Wells without registered
/// fracture contributions are skipped.
template <class FluidSystem, class Scalar, class Simulator,
          class WellContainer, class WellState>
void updateWellFractureFlowSplit(const Simulator& simulator,
                                 const WellContainer& well_container,
                                 WellState& well_state,
                                 const WellState& nupcol_well_state)
{
    for (const auto& well : well_container) {
        if (!well->isInjector()) {
            continue;
        }

        const auto& fracture_indices = well->wellIndexFracture();
        if (fracture_indices.empty()) {
            continue; // no fracture crosses this well
        }

        const auto& wellstate_nupcol = nupcol_well_state.well(well->indexOfWell());
        auto& single_wellstate = well_state.well(well->indexOfWell());
        auto& perf_data = single_wellstate.perf_data;
        auto& filtrate_data = perf_data.filtrate_data;
        auto& fracture_data = perf_data.fracture_data;

        const int nperf = static_cast<int>(well->wellIndex().size());
        assert(static_cast<int>(fracture_indices.size()) == nperf);

        Scalar total_flow_fracture = 0.0;
        const int np = well_state.numPhases();

        auto obtain = [](const auto& value) { return getValue(value); };

        for (int perf = 0; perf < nperf; ++perf) {
            const auto cell_idx = well->perforationData()[perf].cell_index;
            const auto& intQuants =
                simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
            const auto trans_mult = simulator.problem()
                .template wellTransMultiplier<Scalar>(intQuants, cell_idx, obtain);

            const Scalar matrix_wi = well->wellIndex()[perf] * trans_mult;
            const Scalar perf_pressure = wellstate_nupcol.perf_data.pressure[perf];
            const Scalar fracture_wi = fracture_indices[perf].wellIndex(perf_pressure);
            const Scalar effective_wi = matrix_wi + fracture_wi;

            auto& fracture_rate = fracture_data.water_rate[perf];
            if (effective_wi > 0.0) {
                const Scalar matrix_frac = matrix_wi / effective_wi;
                filtrate_data.flow_factor[perf] = matrix_frac;
                fracture_rate = (1.0 - matrix_frac)
                    * perf_data.phase_rates[perf*np + FluidSystem::waterPhaseIdx];
                total_flow_fracture += fracture_rate;
            }
            else {
                filtrate_data.flow_factor[perf] = 0.0;
                fracture_rate = 0.0;
            }
        }

        single_wellstate.frac_rate = total_flow_fracture;
    }
}

} // namespace Opm

#endif // OPM_WELL_FRACTURE_FLOW_SPLIT_HPP_INCLUDED
