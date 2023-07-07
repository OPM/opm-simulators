/*
  Copyright 2023 Equinor

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

#include <config.h>
#include <opm/simulators/wells/WellFilterCake.hpp>

#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

namespace Opm {

void WellFilterCake::
updateFiltrationParticleVolume(const WellInterfaceGeneric& well,
                               const double dt,
                               const std::size_t water_index,
                               const WellState& well_state,
                               std::vector<double>& filtration_particle_volume)
{
    if (!well.isInjector()) {
        return;
    }

    const auto injectorType = well.wellEcl().injectorType();
    if (injectorType != InjectorType::WATER) {
        return;
    }

    const double conc = well.wellEcl().getFilterConc();
    if (conc == 0.) {
        return;
    }

    // it is currently used for the filter cake modeling related to formation damage study
    auto& ws = well_state.well(well.indexOfWell());
    const auto& connection_rates = ws.perf_data.phase_rates;

    const std::size_t np = well_state.numPhases();
    for (int perf = 0; perf < well.numPerfs(); ++perf) {
        // not considering the production water
        const double water_rates = std::max(0., connection_rates[perf * np + water_index]);
        filtration_particle_volume[perf] += water_rates * conc * dt;
    }
}

} // namespace Opm
