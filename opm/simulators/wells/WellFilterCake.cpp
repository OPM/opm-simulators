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

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <fmt/format.h>

#include <stdexcept>

namespace Opm {

template<typename Scalar, typename IndexTraits>
void WellFilterCake<Scalar, IndexTraits>::
updatePostStep(const WellInterfaceGeneric<Scalar, IndexTraits>& well,
               WellState<Scalar, IndexTraits>& well_state,
               const double dt,
               const Scalar conc,
               const std::size_t water_index,
               DeferredLogger& deferred_logger)
{
    if (! well.isInjector() || well.wellEcl().injectorType() != InjectorType::WATER)
        return;

    assert (conc > 0.);
    auto& ws = well_state.well(well.indexOfWell());
    const auto nperf = ws.perf_data.size();
    if (skin_factor_.empty()) {
        skin_factor_.assign(nperf, 0.);
        thickness_.assign(nperf, 0.);
    }

    ws.filtrate_conc = conc;

    updateSkinFactorsAndMultipliers(well, well_state, dt, water_index, deferred_logger);
}

template<typename Scalar, typename IndexTraits>
void WellFilterCake<Scalar, IndexTraits>::
updatePreStep(const WellInterfaceGeneric<Scalar, IndexTraits>& well, DeferredLogger& deferred_logger)
{
    // Apply cleaning and reset any filter cake cleaning multipliers (even if the well is producing at this time)
    applyCleaning(well, deferred_logger);
}

template<typename Scalar, typename IndexTraits>
void WellFilterCake<Scalar, IndexTraits>::
applyCleaning(const WellInterfaceGeneric<Scalar, IndexTraits>& well,
              DeferredLogger& deferred_logger)
{
    const auto& connections = well.wellEcl().getConnections();
    const auto nperf = well.numLocalPerfs();
    for (int perf = 0; perf < nperf; ++perf) {
        const auto perf_ecl_index = well.perforationData()[perf].ecl_index;
        const auto& connection = connections[perf_ecl_index];
        if (!connection.filterCakeActive())
            continue;

        const auto& filter_cake = connection.getFilterCake();
        const Scalar factor = filter_cake.sf_multiplier;
        if (factor == 1.0)
            continue;

        filter_cake.sf_multiplier = 1.0;
        skin_factor_[perf] *= factor;
        updateMultiplier(connection, perf);
        const Scalar rw = connection.getFilterCakeRadius();
        switch (filter_cake.geometry) {
            case FilterCake::FilterCakeGeometry::LINEAR:
            case FilterCake::FilterCakeGeometry::LINRAD: {
                // Previous thickness adjusted to give correct cleaning multiplier at start of time step
                thickness_[perf] *= factor;
                break;
            }
            case FilterCake::FilterCakeGeometry::RADIAL: {
                Scalar rc_prev = std::sqrt(rw * rw + 2. * rw * thickness_[perf]);
                // Previous thickness and rc adjusted to give correct cleaning multiplier at start of time step
                rc_prev = rw*std::exp(factor*std::log(rc_prev/rw));
                thickness_[perf] = (rc_prev * rc_prev - rw * rw) / (2. * rw);
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
    }
}


template<typename Scalar, typename IndexTraits>
void WellFilterCake<Scalar, IndexTraits>::
updateSkinFactorsAndMultipliers(const WellInterfaceGeneric<Scalar, IndexTraits>& well,
                  WellState<Scalar, IndexTraits>& well_state,
                  const double dt,
                  const std::size_t water_index,
                  DeferredLogger& deferred_logger)
{
    const auto nperf = well.numLocalPerfs();
    inj_fc_multiplier_.assign(nperf, 1.0);

    assert(well.isInjector());

    const auto& connections = well.wellEcl().getConnections();
    auto& ws = well_state.well(well.indexOfWell());
    auto& perf_data = ws.perf_data;
    assert (nperf == static_cast<int>(perf_data.size()));

    const auto& connection_rates = perf_data.phase_rates;
    const auto conc = ws.filtrate_conc;
    const std::size_t np = well_state.numPhases();

    for (int perf = 0; perf < nperf; ++perf) {
        const auto perf_ecl_index = well.perforationData()[perf].ecl_index;
        const auto& connection = connections[perf_ecl_index];
        if (!connection.filterCakeActive())
            continue;

        // not considering the production water
        const Scalar water_rates = std::max(Scalar{0.}, connection_rates[perf * np + water_index]);
        const Scalar filtrate_rate = water_rates * conc;
        const Scalar filtrate_particle_volume = filtrate_rate * dt;
        auto& filtrate_data = perf_data.filtrate_data;
        filtrate_data.rates[perf] = filtrate_rate;
        filtrate_data.total[perf] += filtrate_particle_volume;

        const auto& filter_cake = connection.getFilterCake();
        const Scalar area = connection.getFilterCakeArea();
        const Scalar poro = filter_cake.poro;
        const Scalar perm = filter_cake.perm;
        const Scalar rw = connection.getFilterCakeRadius();
        const Scalar K = connection.Kh() / connection.connectionLength();
        // The thickness of the filtration cake due to particle deposition at current time step
        const Scalar delta_thickness = filtrate_particle_volume / (area * (1. - poro));

        filtrate_data.poro[perf] = poro;
        filtrate_data.perm[perf] = perm;
        filtrate_data.radius[perf] = connection.getFilterCakeRadius();
        filtrate_data.area_of_flow[perf] = connection.getFilterCakeArea();

        Scalar delta_skin_factor = 0.;
        Scalar thickness = 0.;
        switch (filter_cake.geometry) {
            case FilterCake::FilterCakeGeometry::LINEAR: {
                thickness = thickness_[perf] + delta_thickness;
                filtrate_data.thickness[perf] = thickness;
                delta_skin_factor = delta_thickness / rw * K / perm;
                break;
            }
            case FilterCake::FilterCakeGeometry::RADIAL: {
                const auto prev_thickness = thickness_[perf];
                Scalar rc_prev = std::sqrt(rw * rw + 2. * rw * prev_thickness);
                thickness = prev_thickness + delta_thickness;

                const Scalar rc = std::sqrt(rw * rw + 2. * rw * thickness);
                filtrate_data.thickness[perf] = rc - rw;
                delta_skin_factor = K / perm * std::log(rc / rc_prev);
                break;
            }
            case FilterCake::FilterCakeGeometry::LINRAD: {
                const auto prev_thickness = thickness_[perf];
                Scalar rc_prev = rw + prev_thickness;
                thickness = thickness_[perf] + delta_thickness;
                filtrate_data.thickness[perf] = thickness;
                const Scalar rc = rw + thickness;
                delta_skin_factor = K / perm * std::log(rc / rc_prev);
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
        skin_factor_[perf] += delta_skin_factor;
        filtrate_data.skin_factor[perf] = skin_factor_[perf];
        thickness_[perf] = thickness;

        updateMultiplier(connection, perf);
    }
}

template<typename Scalar, typename IndexTraits>
template <class Conn>
void WellFilterCake<Scalar, IndexTraits>::
updateMultiplier(const Conn& connection, const int perf)
{
    const auto denom = connection.ctfProperties().peaceman_denom;
    const auto denom2 = denom + skin_factor_[perf];
    inj_fc_multiplier_[perf] = denom / denom2;
}

template class WellFilterCake<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class WellFilterCake<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
