/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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

#include <opm/simulators/wells/StandardWellConnections.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <variant>
#include <vector>

#include <fmt/format.h>

namespace {
    // Component and phase rates/mixtures at surface conditions are linked as
    //
    // [ componentMixture[oilpos]   ]   [ 1     Rv     0  ] [ phaseMixture[oilpos]  ]
    // [ componentMixture[gaspos]   ] = [ Rs     1    Rsw ] [ phaseMixture[gaspos]  ]
    // [ componentMixture[waterpos] ]   [ 0     Rvw    1  ] [ phaseMixture[waterpos ]
    //
    // The reapportion*Mixture() functions calculates phase rates/mixtures
    // in the oil/gas and gas/water systems, respectively, based on known
    // component rates/mixtures.

    template <typename Scalar, typename Properties>
    void reapportionGasOilMixture(std::string_view                              well,
                                  const typename std::vector<Scalar>::size_type gaspos,
                                  const typename std::vector<Scalar>::size_type oilpos,
                                  const typename std::vector<Scalar>::size_type perf,
                                  const Properties&                             props,
                                  const std::vector<Scalar>&                    componentMixture,
                                  std::vector<Scalar>&                          phaseMixture,
                                  Opm::DeferredLogger&                          deferred_logger)
    {
        const auto threshold = static_cast<Scalar>(1.0e-12);

        Scalar rs {0};
        if (!props.rsmax_perf.empty() && (componentMixture[oilpos] > threshold)) {
            rs = std::min(componentMixture[gaspos] / componentMixture[oilpos], props.rsmax_perf[perf]);
        }

        Scalar rv {0};
        if (!props.rvmax_perf.empty() && (componentMixture[gaspos] > threshold)) {
            rv = std::min(componentMixture[oilpos] / componentMixture[gaspos], props.rvmax_perf[perf]);
        }

        const Scalar denom = Scalar{1} - rs*rv;
        // too small denominator can cause small volrat and too big density
        if (! (denom > threshold)) {
            const auto msg = fmt::format(R"(Problematic denominator value {} for well {}.
    Connection density calculation with (Rs, Rv) = ({}, {}).
    Proceeding as if no dissolution or vaporisation for this connection)",
                                         denom, well, rs, rv);

            deferred_logger.debug(msg);
        }
        else {
            if (rs > Scalar{0}) {
                // Subtract gas in oil from gas mixture.
                phaseMixture[gaspos] =
                    (componentMixture[gaspos] - componentMixture[oilpos]*rs) / denom;
            }

            if (rv > Scalar{0}) {
                // Subtract oil in gas from oil mixture.
                phaseMixture[oilpos] =
                    (componentMixture[oilpos] - componentMixture[gaspos]*rv) / denom;
            }
        }
    }

    template <typename Scalar, typename Properties>
    void reapportionGasWaterMixture(std::string_view                              well,
                                    const typename std::vector<Scalar>::size_type gaspos,
                                    const typename std::vector<Scalar>::size_type waterpos,
                                    const typename std::vector<Scalar>::size_type perf,
                                    const Properties&                             props,
                                    const std::vector<Scalar>&                    componentMixture,
                                    std::vector<Scalar>&                          phaseMixture,
                                    Opm::DeferredLogger&                          deferred_logger)
    {
        const auto threshold = static_cast<Scalar>(1.0e-12);

        Scalar rvw {0};
        if (!props.rvwmax_perf.empty() && (componentMixture[gaspos] > threshold)) {
            rvw = std::min(componentMixture[waterpos] / componentMixture[gaspos], props.rvwmax_perf[perf]);
        }

        Scalar rsw {0};
        if (!props.rswmax_perf.empty() && (componentMixture[waterpos] > threshold)) {
            rsw = std::min(componentMixture[gaspos] / componentMixture[waterpos], props.rswmax_perf[perf]);
        }

        const Scalar d = Scalar{1} - rsw*rvw;
        // too small denominator can cause small volrat and too big density
        if (! (d > threshold)) {
            const auto msg = fmt::format(R"(Problematic denominator value {} for well {}.
    Connection density calculation with (Rsw, Rvw) = ({}, {}).
    Proceeding as if no dissolution or vaporisation for this connection)",
                                         d, well, rsw, rvw);

            deferred_logger.debug(msg);
        }
        else {
            if (rsw > Scalar{0}) {
                // Subtract gas in water from gas mixture
                phaseMixture[gaspos] =
                    (componentMixture[gaspos] - componentMixture[waterpos]*rsw) / d;
            }

            if (rvw > Scalar{0}) {
                // Subtract water in gas from water mixture
                phaseMixture[waterpos] =
                    (componentMixture[waterpos] - componentMixture[gaspos]*rvw) / d;
            }
        }
    }

} // Anonymous namespace

namespace Opm
{

template<typename FluidSystem, typename Indices>
StandardWellConnections<FluidSystem, Indices>::
StandardWellConnections(const WellInterfaceIndices<FluidSystem, Indices>& well)
    : well_(well)
    , perf_densities_(well.numLocalPerfs())
    , perf_pressure_diffs_(well.numLocalPerfs())
{
}

template<typename FluidSystem, typename Indices>
void StandardWellConnections<FluidSystem, Indices>::
computePressureDelta()
{
    // Algorithm:

    // We'll assume the perforations are given in order from top to
    // bottom for each well.  By top and bottom we do not necessarily
    // mean in a geometric sense (depth), but in a topological sense:
    // the 'top' perforation is nearest to the surface topologically.
    // Our goal is to compute a pressure delta for each perforation.

    // 1. Compute pressure differences between perforations.
    //    dp_perf will contain the pressure difference between a
    //    perforation and the one above it, except for the first
    //    perforation for each well, for which it will be the
    //    difference to the reference (bhp) depth.

    const int nperf = well_.numLocalPerfs();
    perf_pressure_diffs_.resize(nperf, 0.0);
    auto z_above = well_.parallelWellInfo().communicateAboveValues(well_.refDepth(), well_.perfDepth());

    for (int perf = 0; perf < nperf; ++perf) {
        const Scalar dz = well_.perfDepth()[perf] - z_above[perf];
        perf_pressure_diffs_[perf] = dz * perf_densities_[perf] * well_.gravity();
    }

    // 2. Compute pressure differences to the reference point (bhp) by
    //    accumulating the already computed adjacent pressure
    //    differences, storing the result in dp_perf.
    //    This accumulation must be done per well.
    const auto beg = perf_pressure_diffs_.begin();
    const auto end = perf_pressure_diffs_.end();
    well_.parallelWellInfo().partialSumPerfValues(beg, end);
}

template<typename FluidSystem, typename Indices>
void StandardWellConnections<FluidSystem, Indices>::
computeDensities(const std::vector<Scalar>& perfComponentRates,
                 const Properties& props,
                 DeferredLogger& deferred_logger)
{
    using Ix = typename std::vector<Scalar>::size_type;

    const auto activeGas   = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    const auto activeOil   = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
    const auto activeWater = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);

    // The *pos objects are unused unless the corresponding phase is active.
    // Therefore, it's safe to have -1 (== numeric_limits<size_t>::max()) be
    // their value in that case.
    const auto gaspos = activeGas
        ? static_cast<Ix>(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx))
        : static_cast<Ix>(-1);

    const auto oilpos = activeOil
        ? static_cast<Ix>(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx))
        : static_cast<Ix>(-1);

    const auto waterpos = activeWater
        ? static_cast<Ix>(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx))
        : static_cast<Ix>(-1);

    const int nperf    = this->well_.numLocalPerfs();
    const int num_comp = this->well_.numComponents();

    // 1. Compute the flow (in surface volume units for each component)
    //    exiting up the wellbore from each perforation, taking into account
    //    flow from lower in the well, and in/out-flow at each perforation.
    const auto q_out_perf = this->calculatePerforationOutflow(perfComponentRates);

    // 2. Compute the component mixture at each perforation as the absolute
    //    values of the surface rates divided by their sum.  Then compute
    //    volume ratios (formation factors) for each perforation.  Finally
    //    compute densities for the segments associated with each
    //    perforation.
    auto componentMixture = std::vector<Scalar>(num_comp, Scalar{0});
    auto phaseMixture     = std::vector<Scalar>(num_comp, Scalar{0});

    for (int perf = 0; perf < nperf; ++perf) {
        this->initialiseConnectionMixture(num_comp, perf,
                                          q_out_perf,
                                          phaseMixture,
                                          componentMixture);

        // Initial phase mixture for this perforation is the same as the
        // initialised component mixture, which in turn is possibly just
        // copied from the phase mixture of the previous perforation (i.e.,
        // "perf - 1").
        phaseMixture = componentMixture;

        // Separate components based on flowing rates.
        if (activeGas && activeOil) {
            reapportionGasOilMixture(this->well_.name(), gaspos, oilpos, perf,
                                     props, componentMixture, phaseMixture,
                                     deferred_logger);
        }

        if (activeGas && activeWater) {
            reapportionGasWaterMixture(this->well_.name(), gaspos, waterpos, perf,
                                       props, componentMixture, phaseMixture,
                                       deferred_logger);
        }

        // Compute connection level mixture density as a weighted average of
        // phase densities.
        const auto* const rho_s = &props.surf_dens_perf[perf*num_comp + 0];
        const auto* const b     = &props.b_perf        [perf*num_comp + 0];

        auto& rho = this->perf_densities_[perf];

        auto volrat = rho = Scalar{0};
        for (auto comp = 0*num_comp; comp < num_comp; ++comp) {
            rho    += componentMixture[comp] * rho_s[comp];
            volrat += phaseMixture    [comp] / b    [comp];
        }

        rho /= volrat;
    }
}

template<typename FluidSystem, typename Indices>
std::vector<typename FluidSystem::Scalar>
StandardWellConnections<FluidSystem, Indices>::
calculatePerforationOutflow(const std::vector<Scalar>& perfComponentRates) const
{
    const int nperf    = this->well_.numLocalPerfs();
    const int num_comp = this->well_.numComponents();

    auto q_out_perf = std::vector<Scalar>(nperf * num_comp, Scalar{0});

    // Component flow rates depend on the order of the perforations.  Thus,
    // we must use the global view of the well's perforation to get an
    // accurate picture of the component flow rates at the perforation
    // level.
    const auto& factory = this->well_.parallelWellInfo()
        .getGlobalPerfContainerFactory();

    auto global_q_out_perf = factory.createGlobal(q_out_perf, num_comp);

    const auto global_perf_comp_rates = factory
        .createGlobal(perfComponentRates, num_comp);

    // TODO: Investigate whether we should use the following techniques to
    // calcuate the composition of flows in the wellbore.  Iterate over well
    // perforations from bottom to top.
    for (int perf = factory.numGlobalPerfs() - 1; perf >= 0; --perf) {
        for (int component = 0; component < num_comp; ++component) {
            const auto index = perf*num_comp + component;
            auto& q_out = global_q_out_perf[index];

            // Initialise current perforation's component flow rate to that
            // of the perforation below.  Bottom perforation has zero flow
            // component rate.
            q_out = (perf == factory.numGlobalPerfs() - 1)
                ? Scalar{0}
                : global_q_out_perf[index + num_comp];

            // Subtract outflow through perforation.
            q_out -= global_perf_comp_rates[index];
        }
    }

    // Copy the data back to local view.
    factory.copyGlobalToLocal(global_q_out_perf, q_out_perf, num_comp);

    return q_out_perf;
}

template<typename FluidSystem, typename Indices>
void StandardWellConnections<FluidSystem, Indices>::
initialiseConnectionMixture(const int                  num_comp,
                            const int                  perf,
                            const std::vector<Scalar>& q_out_perf,
                            const std::vector<Scalar>& phaseMixture,
                            std::vector<Scalar>&       componentMixture) const
{
    // Find component mix.
    const auto tot_surf_rate =
        std::accumulate(q_out_perf.begin() + num_comp*(perf + 0),
                        q_out_perf.begin() + num_comp*(perf + 1), Scalar{0});

    if (tot_surf_rate != Scalar{0}) {
        const auto* const qo = &q_out_perf[perf*num_comp + 0];

        for (int component = 0; component < num_comp; ++component) {
            componentMixture[component] = std::abs(qo[component] / tot_surf_rate);
        }
    }
    else if (num_comp == 1) {
        componentMixture[num_comp - 1] = Scalar{1};
    }
    else {
        std::fill(componentMixture.begin(), componentMixture.end(), Scalar{0});

        // No flow => use fractions defined at well level for componentMixture.
        if (this->well_.isInjector()) {
            switch (this->well_.wellEcl().injectorType()) {
            case InjectorType::WATER:
                componentMixture[IndexTraits::waterCompIdx] = Scalar{1};
                break;

            case InjectorType::GAS:
                componentMixture[IndexTraits::gasCompIdx] = Scalar{1};
                break;

            case InjectorType::OIL:
                componentMixture[IndexTraits::oilCompIdx] = Scalar{1};
                break;

            case InjectorType::MULTI:
                // Not supported.
                // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                //                         "Multi phase injectors are not supported, requested for well " + name());
                break;
            }
        }
        else {
            assert(this->well_.isProducer());

            if (perf == 0) {
                // For the first perforation without flow we use the
                // preferred phase to decide the componentMixture initialization.

                switch (this->well_.wellEcl().getPreferredPhase()) {
                case Phase::OIL:
                    componentMixture[IndexTraits::oilCompIdx] = Scalar{1};
                    break;

                case Phase::GAS:
                    componentMixture[IndexTraits::gasCompIdx] = Scalar{1};
                    break;

                case Phase::WATER:
                    componentMixture[IndexTraits::waterCompIdx] = Scalar{1};
                    break;

                default:
                    // No others supported.
                    break;
                }
            }
            else {
                // For the rest of the perforations without flow we use the
                // componentMixture from the perforation above.
                componentMixture = phaseMixture;
            }
        }
    }
}

template<typename FluidSystem, typename Indices>
void StandardWellConnections<FluidSystem, Indices>::
computeDensitiesForStoppedProducer(const DensityPropertyFunctions& prop_func)
{
    const auto np = this->well_.numPhases();

    const auto modPhIx = [this, np]() {
        auto phIx = std::vector<int>(np);

        for (auto p = 0*np; p < np; ++p) {
            phIx[p] = this->well_.flowPhaseToModelPhaseIdx(p);
        }

        return phIx;
    }();

    auto mob = std::vector<Scalar>(np);
    auto rho = std::vector<Scalar>(np);

    const auto nperf = this->well_.numLocalPerfs();
    for (auto perf = 0*nperf; perf < nperf; ++perf) {
        const auto cell_idx = this->well_.cells()[perf];

        prop_func.mobility     (cell_idx, modPhIx, mob);
        prop_func.densityInCell(cell_idx, modPhIx, rho);

        auto& rho_i = this->perf_densities_[perf];

        auto mt = rho_i = Scalar{0};
        for (auto p = 0*np; p < np; ++p) {
            rho_i += mob[p] * rho[p];
            mt    += mob[p];
        }

        rho_i /= mt;
    }
}

template<typename FluidSystem, typename Indices>
typename StandardWellConnections<FluidSystem, Indices>::Properties
StandardWellConnections<FluidSystem, Indices>::
computePropertiesForPressures(const WellState<Scalar, IndexTraits>&         well_state,
                              const PressurePropertyFunctions& prop_func) const
{
    auto props = Properties{};

    const int nperf = well_.numLocalPerfs();

    props.b_perf        .resize(nperf * this->well_.numComponents());
    props.surf_dens_perf.resize(nperf * this->well_.numComponents());

    const auto& ws = well_state.well(this->well_.indexOfWell());

    const bool waterPresent = FluidSystem::phaseIsActive(IndexTraits::waterPhaseIdx);
    const bool oilPresent = FluidSystem::phaseIsActive(IndexTraits::oilPhaseIdx);
    const bool gasPresent = FluidSystem::phaseIsActive(IndexTraits::gasPhaseIdx);

    // Rs/Rv are used only if both oil and gas are present.
    if (oilPresent && gasPresent) {
        props.rsmax_perf.resize(nperf);
        props.rvmax_perf.resize(nperf);
    }

    // Rsw/Rvw are used only if both water and gas are present.
    if (waterPresent && gasPresent) {
        props.rvwmax_perf.resize(nperf);
        props.rswmax_perf.resize(nperf);
    }

    // Compute the average pressure in each well block
    const auto& perf_press = ws.perf_data.pressure;
    const auto  p_above    = this->well_.parallelWellInfo()
        .communicateAboveValues(ws.bhp, perf_press.data(), nperf);

    for (int perf = 0; perf < nperf; ++perf) {
        const int cell_idx = well_.cells()[perf];

        const Scalar p_avg = (perf_press[perf] + p_above[perf])/2;
        const Scalar temperature = prop_func.getTemperature(cell_idx, IndexTraits::oilPhaseIdx);
        const Scalar saltConcentration = prop_func.getSaltConcentration(cell_idx);
        const int region_idx = prop_func.pvtRegionIdx(cell_idx);

        if (waterPresent) {
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(IndexTraits::waterCompIdx);
            Scalar rsw = 0.0;
            if (FluidSystem::enableDissolvedGasInWater()) {
                // TODO support mutual solubility in water and oil
                assert(!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
                const int water_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
                const Scalar waterrate = std::abs(ws.surface_rates[water_pos]);
                props.rswmax_perf[perf] = FluidSystem::waterPvt().saturatedGasDissolutionFactor(region_idx, temperature, p_avg, saltConcentration);
                if (waterrate > 0) {
                    const int gas_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
                    const Scalar gasrate = std::abs(ws.surface_rates[gas_pos]) - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);
                    if (gasrate > 0) {
                        rsw = gasrate / waterrate;
                    }
                    rsw = std::min(rsw, props.rswmax_perf[perf]);
                }
            }

            props.b_perf[waterCompIdx + perf * well_.numComponents()] = FluidSystem::waterPvt()
                .inverseFormationVolumeFactor(region_idx, temperature, p_avg, rsw, saltConcentration);
        }

        if (gasPresent) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const int gaspos = gasCompIdx + perf * well_.numComponents();

            Scalar rvw = 0.0;
            Scalar rv = 0.0;
            if (oilPresent) {
                // in order to handle negative rates in producers
                const int oil_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
                const Scalar oilrate = std::abs(ws.surface_rates[oil_pos]);
                props.rvmax_perf[perf] = FluidSystem::gasPvt()
                    .saturatedOilVaporizationFactor(region_idx, temperature, p_avg);

                if (oilrate > 0) {
                    const int gas_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
                    const Scalar gasrate = std::abs(ws.surface_rates[gas_pos]) - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);
                    if (gasrate > 0) {
                        rv = oilrate / gasrate;
                    }
                    rv = std::min(rv, props.rvmax_perf[perf]);
                }
            }

            if (waterPresent) {
                // in order to handle negative rates in producers
                const int water_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
                const Scalar waterrate = std::abs(ws.surface_rates[water_pos]);
                props.rvwmax_perf[perf] = FluidSystem::gasPvt()
                    .saturatedWaterVaporizationFactor(region_idx, temperature, p_avg);

                if (waterrate > 0) {
                    const int gas_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
                    const Scalar gasrate = std::abs(ws.surface_rates[gas_pos])
                        - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);

                    if (gasrate > 0) {
                        rvw = waterrate / gasrate;
                    }

                    rvw = std::min(rvw, props.rvwmax_perf[perf]);
                }
            }

            props.b_perf[gaspos] = FluidSystem::gasPvt()
                .inverseFormationVolumeFactor(region_idx, temperature, p_avg, rv, rvw);
        }

        if (oilPresent) {
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            const int oilpos = oilCompIdx + perf * well_.numComponents();

            Scalar rs = 0.0;
            if (gasPresent) {
                props.rsmax_perf[perf] = FluidSystem::oilPvt()
                    .saturatedGasDissolutionFactor(region_idx, temperature, p_avg);

                const int gas_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
                const Scalar gasrate = std::abs(ws.surface_rates[gas_pos])
                    - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);

                if (gasrate > 0) {
                    const int oil_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
                    const Scalar oilrate = std::abs(ws.surface_rates[oil_pos]);
                    if (oilrate > 0) {
                        rs = gasrate / oilrate;
                    }
                    rs = std::min(rs, props.rsmax_perf[perf]);
                }
            }

            props.b_perf[oilpos] = FluidSystem::oilPvt()
                .inverseFormationVolumeFactor(region_idx, temperature, p_avg, rs);
        }

        // Surface density.
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            props.surf_dens_perf[well_.numComponents() * perf  + compIdx] =
                FluidSystem::referenceDensity( phaseIdx, region_idx );
        }

        // We use cell values for solvent injector
        if constexpr (Indices::enableSolvent) {
            props.b_perf[well_.numComponents() * perf + Indices::contiSolventEqIdx] =
                prop_func.solventInverseFormationVolumeFactor(cell_idx);

            props.surf_dens_perf[well_.numComponents() * perf + Indices::contiSolventEqIdx] =
                prop_func.solventRefDensity(cell_idx);
        }
    }

    return props;
}

template<typename FluidSystem, typename Indices>
std::vector<typename FluidSystem::Scalar>
StandardWellConnections<FluidSystem, Indices>::
copyInPerforationRates(const Properties&       props,
                       const PerfData<Scalar>& perf_data) const
{
    auto perfRates = std::vector<Scalar>(props.b_perf.size(), Scalar{0});

    const int nperf = this->well_.numLocalPerfs();
    const int np    = this->well_.numPhases();
    const int nc    = this->well_.numComponents();

    const auto srcIx = [this, np]() {
        auto ix = std::vector<int>(np);

        for (auto comp = 0; comp < np; ++comp) {
            ix[comp] = this->well_.modelCompIdxToFlowCompIdx(comp);
        }

        return ix;
    }();

    for (int perf = 0; perf < nperf; ++perf) {
        const auto* const src = &perf_data.phase_rates[perf*np + 0];
        auto* const       dst = &perfRates            [perf*nc + 0];

        for (int comp = 0; comp < np; ++comp) {
            dst[comp] = src[srcIx[comp]];
        }
    }

    if constexpr (Indices::enableSolvent) {
        for (int perf = 0, ix_s = 0*nc + Indices::contiSolventEqIdx;
             perf < nperf; ++perf, ix_s += nc)
        {
            perfRates[ix_s] = perf_data.solvent_rates[perf];
        }
    }

    return perfRates;
}

template<typename FluidSystem, typename Indices>
void StandardWellConnections<FluidSystem, Indices>::
computeProperties(const bool                      stopped_or_zero_rate_target,
                  const WellState<Scalar, IndexTraits>&        well_state,
                  const DensityPropertyFunctions& prop_func,
                  const Properties&               props,
                  DeferredLogger&                 deferred_logger)
{
    // Step 1: Compute mixture densities in well-bore.  Result values stored
    // in data member perf_densities_.
    if (stopped_or_zero_rate_target && this->well_.isProducer()) {
        this->computeDensitiesForStoppedProducer(prop_func);
    }
    else {
        // Injector or flowing producer.
        const auto perfRates = this->copyInPerforationRates
            (props, well_state.well(this->well_.indexOfWell()).perf_data);

        this->computeDensities(perfRates, props, deferred_logger);
    }

    // Step 2: Use those mixture densities to calculate pressure drops along
    // well-bore.  Result values stored in data member perf_pressure_diffs_.
    this->computePressureDelta();
}

template<typename FluidSystem, typename Indices>
typename StandardWellConnections<FluidSystem, Indices>::Eval
StandardWellConnections<FluidSystem, Indices>::
connectionRateBrine(Scalar& rate,
                    const Scalar vap_wat_rate,
                    const std::vector<EvalWell>& cq_s,
                    const std::variant<Scalar,EvalWell>& saltConcentration) const
{
    // TODO: the application of well efficiency factor has not been tested with an example yet
    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    // Correction salt rate; evaporated water does not contain salt
    EvalWell cq_s_sm = cq_s[waterCompIdx] - vap_wat_rate;
    if (well_.isInjector()) {
        cq_s_sm *= std::get<Scalar>(saltConcentration);
    } else {
        cq_s_sm *= std::get<EvalWell>(saltConcentration);
    }

    // Note. Efficiency factor is handled in the output layer
    rate = cq_s_sm.value();

    cq_s_sm *= well_.wellEfficiencyFactor();
    return well_.restrictEval(cq_s_sm);
}

template<typename FluidSystem, typename Indices>
typename StandardWellConnections<FluidSystem, Indices>::Eval
StandardWellConnections<FluidSystem, Indices>::
connectionRateFoam(const std::vector<EvalWell>& cq_s,
                    const std::variant<Scalar,EvalWell>& foamConcentration,
                    const Phase transportPhase,
                    DeferredLogger& deferred_logger) const
{
    // TODO: the application of well efficiency factor has not been tested with an example yet
    auto getFoamTransportIdx = [&deferred_logger,transportPhase] {
        switch (transportPhase) {
            case Phase::WATER: {
                return Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            }
            case Phase::GAS: {
                return Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            }
            case Phase::SOLVENT: {
                if constexpr (Indices::enableSolvent)
                    return Indices::contiSolventEqIdx;
                else
                    OPM_DEFLOG_THROW(std::runtime_error, "Foam transport phase is SOLVENT but SOLVENT is not activated.", deferred_logger);
            }
            default: {
                OPM_DEFLOG_THROW(std::runtime_error, "Foam transport phase must be GAS/WATER/SOLVENT.", deferred_logger);
            }
        }
    };

    EvalWell cq_s_foam = cq_s[getFoamTransportIdx()] * well_.wellEfficiencyFactor();;
    if (well_.isInjector()) {
        cq_s_foam *= std::get<Scalar>(foamConcentration);
    } else {
        cq_s_foam *= std::get<EvalWell>(foamConcentration);
    }

    return well_.restrictEval(cq_s_foam);
}

template<typename FluidSystem, typename Indices>
std::tuple<typename StandardWellConnections<FluidSystem, Indices>::Eval,
           typename StandardWellConnections<FluidSystem, Indices>::Eval,
           typename StandardWellConnections<FluidSystem, Indices>::Eval>
StandardWellConnections<FluidSystem, Indices>::
connectionRatesMICP(Scalar& rate_m,
                    Scalar& rate_o,
                    Scalar& rate_u,
                    const std::vector<EvalWell>& cq_s,
                    const std::variant<Scalar,EvalWell>& microbialConcentration,
                    const std::variant<Scalar,EvalWell>& oxygenConcentration,
                    const std::variant<Scalar,EvalWell>& ureaConcentration) const
{
    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    EvalWell cq_s_microbe = cq_s[waterCompIdx];
    if (well_.isInjector()) {
        cq_s_microbe *= std::get<Scalar>(microbialConcentration);
    } else {
        cq_s_microbe *= std::get<EvalWell>(microbialConcentration);
    }

    rate_m = cq_s_microbe.value();

    EvalWell cq_s_oxygen = cq_s[waterCompIdx];
    if (well_.isInjector()) {
        cq_s_oxygen *= std::get<Scalar>(oxygenConcentration);
    } else {
        cq_s_oxygen *= std::get<EvalWell>(oxygenConcentration);
    }

    rate_o = cq_s_oxygen.value();

    EvalWell cq_s_urea = cq_s[waterCompIdx];
    if (well_.isInjector()) {
        cq_s_urea *= std::get<Scalar>(ureaConcentration);
    } else {
        cq_s_urea *= std::get<EvalWell>(ureaConcentration);
    }

    rate_u = cq_s_urea.value();

    return {well_.restrictEval(cq_s_microbe),
            well_.restrictEval(cq_s_oxygen),
            well_.restrictEval(cq_s_urea)};
}

template<typename FluidSystem, typename Indices>
std::tuple<typename StandardWellConnections<FluidSystem, Indices>::Eval,
           typename StandardWellConnections<FluidSystem, Indices>::EvalWell>
StandardWellConnections<FluidSystem, Indices>::
connectionRatePolymer(Scalar& rate,
                      const std::vector<EvalWell>& cq_s,
                      const std::variant<Scalar,EvalWell>& polymerConcentration) const
{
    // TODO: the application of well efficiency factor has not been tested with an example yet
    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
    EvalWell cq_s_poly = cq_s[waterCompIdx];
    if (well_.isInjector()) {
        cq_s_poly *= std::get<Scalar>(polymerConcentration);
    } else {
        cq_s_poly *= std::get<EvalWell>(polymerConcentration);
    }
    // Note. Efficiency factor is handled in the output layer
    rate = cq_s_poly.value();

    cq_s_poly *= well_.wellEfficiencyFactor();

    return {well_.restrictEval(cq_s_poly), cq_s_poly};
}

template<typename FluidSystem, typename Indices>
std::tuple<typename StandardWellConnections<FluidSystem, Indices>::Eval,
           typename StandardWellConnections<FluidSystem, Indices>::EvalWell>
StandardWellConnections<FluidSystem, Indices>::
connectionRatezFraction(Scalar& rate,
                        const Scalar dis_gas_rate,
                        const std::vector<EvalWell>& cq_s,
                        const std::variant<Scalar, std::array<EvalWell,2>>& solventConcentration) const
{
    // TODO: the application of well efficiency factor has not been tested with an example yet
    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
    EvalWell cq_s_zfrac_effective = cq_s[gasCompIdx];
    if (well_.isInjector()) {
        cq_s_zfrac_effective *= std::get<Scalar>(solventConcentration);
    } else if (cq_s_zfrac_effective.value() != 0.0) {
        const Scalar dis_gas_frac = dis_gas_rate / cq_s_zfrac_effective.value();
        const auto& vol = std::get<std::array<EvalWell,2>>(solventConcentration);
        cq_s_zfrac_effective *= dis_gas_frac * vol[0] + (1.0 - dis_gas_frac) * vol[1];
    }

    rate = cq_s_zfrac_effective.value();

    cq_s_zfrac_effective *= well_.wellEfficiencyFactor();
    return {well_.restrictEval(cq_s_zfrac_effective), cq_s_zfrac_effective};
}

#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

INSTANTIATE_TYPE_INDICES(StandardWellConnections, double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(StandardWellConnections, float)
#endif

} // namespace Opm
