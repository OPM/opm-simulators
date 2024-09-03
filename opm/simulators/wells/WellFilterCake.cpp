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

template<class Scalar>
void WellFilterCake<Scalar>::
updateFiltrationParticleVolume(const WellInterfaceGeneric<Scalar>& well,
                               const double dt,
                               const Scalar conc,
                               const std::size_t water_index,
                               WellState<Scalar>& well_state)
{
    // Reset any filter cake cleaning multipliers (for all wells/perfs)
    const auto& connections = well.wellEcl().getConnections();
    for (auto& conn : connections) {
        if (conn.filterCakeActive()) {
            conn.getFilterCake().sf_multiplier = 1.0;
        }
    }


    if (!well.isInjector()) {
        return;
    }

    auto& ws = well_state.well(well.indexOfWell());
    const auto nperf = ws.perf_data.size();
    if (filtration_particle_volume_.empty()) {
        filtration_particle_volume_.assign(ws.perf_data.size(), 0.); // initializing to be zero (for each step..)
        prev_skin_factor_.assign(nperf, 0.);
        skin_factor_.assign(nperf, 0.);
        thickness_.assign(nperf, 0.);
        prev_thickness_.assign(nperf, 0.);
    }

    prev_skin_factor_ = skin_factor_;
    prev_thickness_ = thickness_;

    const auto injectorType = well.wellEcl().injectorType();
    if (injectorType != InjectorType::WATER) {
        return;
    }

    ws.filtrate_conc = conc;

    if (conc == 0.) {
        return;
    }

    const auto& connection_rates = ws.perf_data.phase_rates;

    const std::size_t np = well_state.numPhases();
    for (int perf = 0; perf < well.numPerfs(); ++perf) {
        // not considering the production water
        const Scalar water_rates = std::max(Scalar{0.}, connection_rates[perf * np + water_index]);
        const Scalar filtrate_rate = water_rates * conc;
        filtration_particle_volume_[perf] = filtrate_rate * dt;
        ws.perf_data.filtrate_data.rates[perf] = filtrate_rate;
        ws.perf_data.filtrate_data.total[perf] += filtration_particle_volume_[perf];
    }
}

template<class Scalar>
void WellFilterCake<Scalar>::
updateInjFCMult(const WellInterfaceGeneric<Scalar>& well,
                WellState<Scalar>& well_state,
                DeferredLogger& deferred_logger)
{
    if (inj_fc_multiplier_.empty()) {
        inj_fc_multiplier_.resize(well.numPerfs(), 1.0);
    }
    auto& ws = well_state.well(well.indexOfWell());
    auto& perf_data = ws.perf_data;

    for (int perf = 0; perf < well.numPerfs(); ++perf) {
        const auto perf_ecl_index = well.perforationData()[perf].ecl_index;
        const auto& connections = well.wellEcl().getConnections();
        const auto& connection = connections[perf_ecl_index];
        if (well.isInjector() && connection.filterCakeActive()) {
            const auto& filter_cake = connection.getFilterCake();
            const Scalar area = connection.getFilterCakeArea();
            const Scalar poro = filter_cake.poro;
            const Scalar perm = filter_cake.perm;
            const Scalar rw = connection.getFilterCakeRadius();
            const Scalar K = connection.Kh() / connection.connectionLength();
            const Scalar factor = filter_cake.sf_multiplier;
            // The thickness of the filtration cake due to particle deposition at current time step
            const Scalar delta_thickness = filtration_particle_volume_[perf] / (area * (1. - poro));
            auto& filtrate_data = perf_data.filtrate_data;
            filtrate_data.poro[perf] = poro;
            filtrate_data.perm[perf] = perm;
            filtrate_data.radius[perf] = connection.getFilterCakeRadius();
            filtrate_data.area_of_flow[perf] = connection.getFilterCakeArea();

            Scalar skin_factor = 0.;
            Scalar thickness = 0.;
            switch (filter_cake.geometry) {
                case FilterCake::FilterCakeGeometry::LINEAR: {
                    // Previuos thickness adjusted to give correct cleaning multiplier at start of time step
                    const auto prev_thickness = prev_thickness_[perf]*factor;
                    thickness = prev_thickness + delta_thickness;
                    skin_factor = thickness / rw * K / perm;
                    filtrate_data.thickness[perf] = thickness;
                    break;
                }
                case FilterCake::FilterCakeGeometry::RADIAL: {
                    auto prev_thickness = prev_thickness_[perf];
                    Scalar rc_prev = std::sqrt(rw * rw + 2. * rw * prev_thickness);
                    // Previuos thickness and rc adjusted to give correct cleaning multiplier at start of time step
                    if (factor != 1.0) {
                        rc_prev = rw*std::exp(factor*std::log(rc_prev/rw));
                        prev_thickness = (rc_prev * rc_prev - rw * rw) / (2. * rw);
                    }
                    thickness = prev_thickness + delta_thickness;

                    //const Scalar rc_prev = std::sqrt(rw * rw + 2. * rw * prev_thickness);
                    const Scalar rc = std::sqrt(rw * rw + 2. * rw * thickness);
                    filtrate_data.thickness[perf] = rc - rw;
                    skin_factor = K / perm * std::log(rc / rc_prev);
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
            skin_factor_[perf] = prev_skin_factor_[perf]*factor + skin_factor;
            filtrate_data.skin_factor[perf] = skin_factor_[perf];
            thickness_[perf] = thickness;

            const auto denom = connection.ctfProperties().peaceman_denom;
            const auto denom2 = denom + skin_factor_[perf];
            inj_fc_multiplier_[perf] = denom / denom2;
        } else {
            inj_fc_multiplier_[perf] = 1.0;
        }
    }
}

template class WellFilterCake<double>;

#if FLOW_INSTANTIATE_FLOAT
template class WellFilterCake<float>;
#endif

} // namespace Opm
