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
#include <opm/simulators/wells/InjectedPhase.hpp>

#include <algorithm>
#include <cassert>

namespace Opm {

/// Split each injector connection's flow of the injected phase (water, gas
/// or oil) between the matrix and a
/// fracture crossing the connection, based on the ratio of the regular
/// connection transmissibility to the fracture contribution registered via
/// WellInterfaceGeneric::addFracturePerforations().
///
/// Fills, per connection, the filtrate flow factor (matrix share) and the
/// fracture rate of the injected phase (surface volumes); the well total is
/// summed over all ranks when reported (SingleWellState::sum_fracture_rate).
/// The results of the previous step are reset first, since they carry over
/// to the next report step, so that wells without registered fracture
/// contributions, producers and MULTI injectors report no fracture flow.
template <class FluidSystem, class Scalar, class Simulator,
          class WellContainer, class WellState>
void updateWellFractureFlowSplit(const Simulator& simulator,
                                 const WellContainer& well_container,
                                 WellState& well_state,
                                 const WellState& nupcol_well_state)
{
    using IndexTraits = typename FluidSystem::IndexTraitsType;

    for (const auto& well : well_container) {
        auto& single_wellstate = well_state.well(well->indexOfWell());
        auto& perf_data = single_wellstate.perf_data;
        auto& filtrate_data = perf_data.filtrate_data;
        auto& fracture_data = perf_data.fracture_data;

        // fracture_data is only sized once a fracture model claims the well.
        std::ranges::fill(fracture_data.injected_rate, Scalar{0});
        std::ranges::fill(filtrate_data.flow_factor, Scalar{1});
        fracture_data.injected_phase = -1;

        if (!well->isInjector()) {
            continue;
        }

        // Set even without fractured connections on this rank: the owner
        // reports the total of a distributed well under this phase.
        const auto inj_type = well->wellEcl().injectorType();
        const auto inj_phase = injectedActivePhaseIdx(well->phaseUsage(), inj_type);
        if (!inj_phase.has_value()) {
            continue; // MULTI injector or injected phase not active
        }
        fracture_data.injected_phase =
            static_cast<int>(*injectedCanonicalPhaseIdx<IndexTraits>(inj_type));

        const auto& fracture_indices = well->wellIndexFracture();
        if (fracture_indices.empty()) {
            continue; // no fracture crosses this well's local connections
        }

        const auto& wellstate_nupcol = nupcol_well_state.well(well->indexOfWell());

        const int nperf = static_cast<int>(well->wellIndex().size());
        assert(static_cast<int>(fracture_indices.size()) == nperf);

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

            auto& fracture_rate = fracture_data.injected_rate[perf];
            if (effective_wi > 0.0) {
                const Scalar matrix_frac = matrix_wi / effective_wi;
                filtrate_data.flow_factor[perf] = matrix_frac;
                fracture_rate = (1.0 - matrix_frac)
                    * perf_data.phase_rates[perf*np + *inj_phase];
            }
            else {
                filtrate_data.flow_factor[perf] = 0.0;
                fracture_rate = 0.0;
            }
        }
    }
}

} // namespace Opm

#endif // OPM_WELL_FRACTURE_FLOW_SPLIT_HPP_INCLUDED
