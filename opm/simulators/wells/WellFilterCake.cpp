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

#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <stdexcept>

namespace Opm {

void WellFilterCake::
updateFiltrationParticleVolume(const WellInterfaceGeneric& well,
                               const double dt,
                               const double conc,
                               const std::size_t water_index,
                               const WellState& well_state)
{
    if (!well.isInjector()) {
        return;
    }

    if (filtration_particle_volume_.empty()) {
        const auto& ws = well_state.well(well.indexOfWell());
        filtration_particle_volume_.assign(ws.perf_data.size(), 0.); // initializing to be zero
    }

    const auto injectorType = well.wellEcl().injectorType();
    if (injectorType != InjectorType::WATER) {
        return;
    }

    // it is currently used for the filter cake modeling related to formation damage study
    auto& ws = well_state.well(well.indexOfWell());
    const auto& connection_rates = ws.perf_data.phase_rates;

    const std::size_t np = well_state.numPhases();
    for (int perf = 0; perf < well.numPerfs(); ++perf) {
        // not considering the production water
        const double water_rates = std::max(0., connection_rates[perf * np + water_index]);
        filtration_particle_volume_[perf] += water_rates * conc * dt;
    }
}

void WellFilterCake::
updateInjFCMult(const WellInterfaceGeneric& well,
                DeferredLogger& deferred_logger)
{
    if (inj_fc_multiplier_.empty()) {
        inj_fc_multiplier_.resize(well.numPerfs(), 1.0);
    }

    for (int perf = 0; perf < well.numPerfs(); ++perf) {
        const auto perf_ecl_index = well.perforationData()[perf].ecl_index;
        const auto& connections = well.wellEcl().getConnections();
        const auto& connection = connections[perf_ecl_index];
        if (well.isInjector() && connection.filterCakeActive()) {
            const auto& filter_cake = connection.getFilterCake();
            const double area = connection.getFilterCakeArea();
            const double poro = filter_cake.poro;
            const double perm = filter_cake.perm;
            const double rw = connection.getFilterCakeRadius();
            const auto cr0 = connection.r0();
            const auto crw = connection.rw();
            const double K = connection.Kh() / connection.connectionLength();
            const double factor = filter_cake.sf_multiplier;
            // the thickness of the filtration cake
            const double thickness = filtration_particle_volume_[perf] / (area * (1. - poro));

            double skin_factor = 0.;
            switch (filter_cake.geometry) {
                case FilterCake::FilterCakeGeometry::LINEAR: {
                    skin_factor = thickness / rw * K / perm * factor;
                    break;
                }
                case FilterCake::FilterCakeGeometry::RADIAL: {
                    const double rc = std::sqrt(rw * rw + 2. * rw * thickness);
                    skin_factor = K / perm * std::log(rc / rw) * factor;
                    break;
                }
                default:
                    const auto geometry =
                        FilterCake::filterCakeGeometryToString(filter_cake.geometry);
                    OPM_DEFLOG_THROW(std::runtime_error,
                                     fmt::format(" Invalid filtration cake geometry type ({}) for well {}",
                                                 geometry, well.name()),
                                     deferred_logger);
            }
            // the original skin factor for the connection
            const auto cskinfactor = connection.skinFactor();
            // compute a multiplier for the well connection transmissibility
            const auto denom = std::log(cr0 / std::min(crw, cr0)) + cskinfactor;
            const auto denom2 = denom + skin_factor;
            inj_fc_multiplier_[perf] = denom / denom2;
        } else {
            inj_fc_multiplier_[perf] = 1.0;
        }
    }
}

} // namespace Opm
