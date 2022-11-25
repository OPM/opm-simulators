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

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/WellInterfaceIndices.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <sstream>

namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
StandardWellConnections<FluidSystem,Indices,Scalar>::
StandardWellConnections(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well)
    : well_(well)
    , perf_densities_(well.numPerfs())
    , perf_pressure_diffs_(well.numPerfs())
{
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellConnections<FluidSystem,Indices,Scalar>::
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

    const int nperf = well_.numPerfs();
    perf_pressure_diffs_.resize(nperf, 0.0);
    auto z_above = well_.parallelWellInfo().communicateAboveValues(well_.refDepth(), well_.perfDepth());

    for (int perf = 0; perf < nperf; ++perf) {
        const double dz = well_.perfDepth()[perf] - z_above[perf];
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

template<class FluidSystem, class Indices, class Scalar>
void StandardWellConnections<FluidSystem,Indices,Scalar>::
computeDensities(const std::vector<Scalar>& perfComponentRates,
                 const std::vector<Scalar>& b_perf,
                 const std::vector<Scalar>& rsmax_perf,
                 const std::vector<Scalar>& rvmax_perf,
                 const std::vector<Scalar>& rvwmax_perf,
                 const std::vector<Scalar>& surf_dens_perf,
                 DeferredLogger& deferred_logger)
{
    // Verify that we have consistent input.
    const int nperf = well_.numPerfs();
    const int num_comp = well_.numComponents();

    // 1. Compute the flow (in surface volume units for each
    //    component) exiting up the wellbore from each perforation,
    //    taking into account flow from lower in the well, and
    //    in/out-flow at each perforation.
    std::vector<Scalar> q_out_perf((nperf)*num_comp, 0.0);

    // Step 1 depends on the order of the perforations. Hence we need to
    // do the modifications globally.
    // Create and get the global perforation information and do this sequentially
    // on each process

    const auto& factory = well_.parallelWellInfo().getGlobalPerfContainerFactory();
    auto global_q_out_perf = factory.createGlobal(q_out_perf, num_comp);
    auto global_perf_comp_rates = factory.createGlobal(perfComponentRates, num_comp);

    // TODO: investigate whether we should use the following techniques to calcuate the composition of flows in the wellbore
    // Iterate over well perforations from bottom to top.
    for (int perf = factory.numGlobalPerfs() - 1; perf >= 0; --perf) {
        for (int component = 0; component < num_comp; ++component) {
            auto index = perf * num_comp + component;
            if (perf == factory.numGlobalPerfs() - 1) {
                // This is the bottom perforation. No flow from below.
                global_q_out_perf[index] = 0.0;
            } else {
                // Set equal to flow from below.
                global_q_out_perf[index] = global_q_out_perf[index + num_comp];
            }
            // Subtract outflow through perforation.
            global_q_out_perf[index] -= global_perf_comp_rates[index];
        }
    }

    // Copy the data back to local view
    factory.copyGlobalToLocal(global_q_out_perf, q_out_perf, num_comp);

    // 2. Compute the component mix at each perforation as the
    //    absolute values of the surface rates divided by their sum.
    //    Then compute volume ratios (formation factors) for each perforation.
    //    Finally compute densities for the segments associated with each perforation.
    std::vector<Scalar> mix(num_comp,0.0);
    std::vector<Scalar> x(num_comp);
    std::vector<Scalar> surf_dens(num_comp);

    for (int perf = 0; perf < nperf; ++perf) {
        // Find component mix.
        const Scalar tot_surf_rate = std::accumulate(q_out_perf.begin() + num_comp*perf,
                                                     q_out_perf.begin() + num_comp*(perf+1), 0.0);
        if (tot_surf_rate != 0.0) {
            for (int component = 0; component < num_comp; ++component) {
                mix[component] = std::fabs(q_out_perf[perf*num_comp + component]/tot_surf_rate);
            }
        } else if (num_comp == 1) {
            mix[num_comp-1] = 1.0;
        } else {
            std::fill(mix.begin(), mix.end(), 0.0);
            // No flow => use well specified fractions for mix.
            if (well_.isInjector()) {
                switch (well_.wellEcl().injectorType()) {
                case InjectorType::WATER:
                    mix[FluidSystem::waterCompIdx] = 1.0;
                    break;
                case InjectorType::GAS:
                    mix[FluidSystem::gasCompIdx] = 1.0;
                    break;
                case InjectorType::OIL:
                    mix[FluidSystem::oilCompIdx] = 1.0;
                    break;
                case InjectorType::MULTI:
                    // Not supported.
                    // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                    //                         "Multi phase injectors are not supported, requested for well " + name());
                    break;
                }
            } else {
                assert(well_.isProducer());
                // For the frist perforation without flow we use the preferred phase to decide the mix initialization.
                if (perf == 0) { //
                    switch (well_.wellEcl().getPreferredPhase()) {
                    case Phase::OIL:
                        mix[FluidSystem::oilCompIdx] = 1.0;
                        break;
                    case Phase::GAS:
                        mix[FluidSystem::gasCompIdx] = 1.0;
                        break;
                    case Phase::WATER:
                        mix[FluidSystem::waterCompIdx] = 1.0;
                        break;
                    default:
                        // No others supported.
                        break;
                    }
                // For the rest of the perforation without flow we use mix from the above perforation.
                } else {
                    mix = x;
                }

            }
        }
        // Compute volume ratio.
        x = mix;

        // Subtract dissolved gas from oil phase and vapporized oil from gas phase and vaporized water from gas phase
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned oilpos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            Scalar rs = 0.0;
            Scalar rv = 0.0;
            if (!rsmax_perf.empty() && mix[oilpos] > 1e-12) {
                rs = std::min(mix[gaspos]/mix[oilpos], rsmax_perf[perf]);
            }
            if (!rvmax_perf.empty() && mix[gaspos] > 1e-12) {
                rv = std::min(mix[oilpos]/mix[gaspos], rvmax_perf[perf]);
            }
            const Scalar d = 1.0 - rs*rv;
            if (d <= 0.0) {
                std::ostringstream sstr;
                sstr << "Problematic d value " << d << " obtained for well " << well_.name()
                     << " during ccomputeConnectionDensities with rs " << rs
                     << ", rv " << rv
                     << " obtaining d " << d
                     << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                     << " for this connection.";
                deferred_logger.debug(sstr.str());
            } else {
                if (rs > 0.0) {
                    // Subtract gas in oil from gas mixture
                    x[gaspos] = (mix[gaspos] - mix[oilpos]*rs)/d;
                }
                if (rv > 0.0) {
                    // Subtract oil in gas from oil mixture
                    x[oilpos] = (mix[oilpos] - mix[gaspos]*rv)/d;
                }
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                //matrix system: (mix[oilpos] = q_os, x[oilpos] = bo*q_or, etc...)
                //┌             ┐   ┌                ┐  ┌           ┐
                //│mix[oilpos]  │   | 1     Rv     0 |  |x[oilpos]  |
                //│mix[gaspos]  │ = │ Rs    1      0 │  │x[gaspos]  │
                //│mix[waterpos]│   │ 0     Rvw    1 │  │x[waterpos │
                //└             ┘   └                ┘  └           ┘
                const unsigned waterpos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                Scalar rvw = 0.0;
                if (!rvwmax_perf.empty() && mix[gaspos] > 1e-12) {
                    rvw = std::min(mix[waterpos]/mix[gaspos], rvwmax_perf[perf]);
                }
                if (rvw > 0.0) {
                    // Subtract water in gas from water mixture
                    x[waterpos] = mix[waterpos] - x[gaspos] * rvw;
                }
            }
        } else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            //no oil
            const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned waterpos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            Scalar rvw = 0.0;
            if (!rvwmax_perf.empty() && mix[gaspos] > 1e-12) {
                rvw = std::min(mix[waterpos]/mix[gaspos], rvwmax_perf[perf]);
            }
            if (rvw > 0.0) {
               // Subtract water in gas from water mixture
               x[waterpos] = mix[waterpos] - mix[gaspos] * rvw;
            }
        }

        Scalar volrat = 0.0;
        for (int component = 0; component < num_comp; ++component) {
            volrat += x[component] / b_perf[perf*num_comp+ component];
        }
        for (int component = 0; component < num_comp; ++component) {
            surf_dens[component] = surf_dens_perf[perf*num_comp+ component];
        }

        // Compute segment density.
        perf_densities_[perf] = std::inner_product(surf_dens.begin(), surf_dens.end(), mix.begin(), 0.0) / volrat;
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellConnections<FluidSystem,Indices,Scalar>::
computePropertiesForPressures(const WellState& well_state,
                              const std::function<Scalar(int,int)>& getTemperature,
                              const std::function<Scalar(int)>& getSaltConcentration,
                              const std::function<int(int)>& pvtRegionIdx,
                              const std::function<Scalar(int)>& solventInverseFormationVolumeFactor,
                              const std::function<Scalar(int)>& solventRefDensity,
                              std::vector<Scalar>& b_perf,
                              std::vector<Scalar>& rsmax_perf,
                              std::vector<Scalar>& rvmax_perf,
                              std::vector<Scalar>& rvwmax_perf,
                              std::vector<Scalar>& surf_dens_perf) const
{
    const int nperf = well_.numPerfs();
    const PhaseUsage& pu = well_.phaseUsage();
    b_perf.resize(nperf * well_.numComponents());
    surf_dens_perf.resize(nperf * well_.numComponents());
    const auto& ws = well_state.well(well_.indexOfWell());

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const bool waterPresent = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
    const bool oilPresent = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
    const bool gasPresent = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);

    //rs and rv are only used if both oil and gas is present
    if (oilPresent && gasPresent) {
        rsmax_perf.resize(nperf);
        rvmax_perf.resize(nperf);
    }
      //rvw is only used if both water and gas is present
    if (waterPresent && gasPresent) {
        rvwmax_perf.resize(nperf);
    }

    // Compute the average pressure in each well block
    const auto& perf_press = ws.perf_data.pressure;
    auto p_above =  well_.parallelWellInfo().communicateAboveValues(ws.bhp,
                                                                    perf_press.data(),
                                                                    nperf);

    for (int perf = 0; perf < nperf; ++perf) {
        const int cell_idx = well_.cells()[perf];

        const double p_avg = (perf_press[perf] + p_above[perf])/2;
        const double temperature = getTemperature(cell_idx, FluidSystem::oilPhaseIdx);
        const double saltConcentration = getSaltConcentration(cell_idx);
        const int region_idx = pvtRegionIdx(cell_idx);

        if (waterPresent) {
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            b_perf[ waterCompIdx + perf * well_.numComponents()] =
            FluidSystem::waterPvt().inverseFormationVolumeFactor(region_idx, temperature, p_avg, saltConcentration);
        }

        if (gasPresent) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const int gaspos = gasCompIdx + perf * well_.numComponents();

            if (oilPresent && waterPresent) {
                const double oilrate = std::abs(ws.surface_rates[pu.phase_pos[Oil]]); //in order to handle negative rates in producers
                const double waterrate = std::abs(ws.surface_rates[pu.phase_pos[Water]]); //in order to handle negative rates in producers
                rvmax_perf[perf] = FluidSystem::gasPvt().saturatedOilVaporizationFactor(region_idx, temperature, p_avg);
                rvwmax_perf[perf] = FluidSystem::gasPvt().saturatedWaterVaporizationFactor(region_idx, temperature, p_avg);
                double rv = 0.0;
                double rvw = 0.0;
                if (oilrate > 0) {
                    const double gasrate = std::abs(ws.surface_rates[pu.phase_pos[Gas]]) - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);
                    if (gasrate > 0) {
                        rv = oilrate / gasrate;
                    }
                    rv = std::min(rv, rvmax_perf[perf]);
                }
                if (waterrate > 0) {
                    const double gasrate = std::abs(ws.surface_rates[pu.phase_pos[Gas]]) - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);
                    if (gasrate > 0) {
                        rvw = waterrate / gasrate;
                    }
                    rvw = std::min(rvw, rvwmax_perf[perf]);
                }
                if (rv > 0.0 || rvw > 0.0){
                    b_perf[gaspos] = FluidSystem::gasPvt().inverseFormationVolumeFactor(region_idx, temperature, p_avg, rv, rvw);
                }
                else {
                    b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(region_idx, temperature, p_avg);
                }
            } else if (oilPresent) {
                //no water
                const double oilrate = std::abs(ws.surface_rates[pu.phase_pos[Oil]]); //in order to handle negative rates in producers
                rvmax_perf[perf] = FluidSystem::gasPvt().saturatedOilVaporizationFactor(region_idx, temperature, p_avg);
                if (oilrate > 0) {
                    const double gasrate = std::abs(ws.surface_rates[pu.phase_pos[Gas]]) - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);
                    double rv = 0.0;
                    if (gasrate > 0) {
                        rv = oilrate / gasrate;
                    }
                    rv = std::min(rv, rvmax_perf[perf]);

                    b_perf[gaspos] = FluidSystem::gasPvt().inverseFormationVolumeFactor(region_idx, temperature, p_avg, rv, 0.0 /*Rvw*/);
                }
                else {
                    b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(region_idx, temperature, p_avg);
                }
            } else if (waterPresent) {
                //no oil
                const double waterrate = std::abs(ws.surface_rates[pu.phase_pos[Water]]); //in order to handle negative rates in producers
                rvwmax_perf[perf] = FluidSystem::gasPvt().saturatedWaterVaporizationFactor(region_idx, temperature, p_avg);
                if (waterrate > 0) {
                    const double gasrate = std::abs(ws.surface_rates[pu.phase_pos[Gas]]) - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);
                    double rvw = 0.0;
                    if (gasrate > 0) {
                        rvw = waterrate / gasrate;
                    }
                    rvw = std::min(rvw, rvwmax_perf[perf]);

                    b_perf[gaspos] = FluidSystem::gasPvt().inverseFormationVolumeFactor(region_idx, temperature, p_avg, 0.0 /*Rv*/, rvw);
                }
                else {
                    b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(region_idx, temperature, p_avg);
                }

            } else {
                b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(region_idx, temperature, p_avg);
            }
        }

        if (oilPresent) {
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            const int oilpos = oilCompIdx + perf * well_.numComponents();
            if (gasPresent) {
                rsmax_perf[perf] = FluidSystem::oilPvt().saturatedGasDissolutionFactor(region_idx, temperature, p_avg);
                const double gasrate = std::abs(ws.surface_rates[pu.phase_pos[Gas]]) - (Indices::enableSolvent ? ws.sum_solvent_rates() : 0.0);
                if (gasrate > 0) {
                    const double oilrate = std::abs(ws.surface_rates[pu.phase_pos[Oil]]);
                    double rs = 0.0;
                    if (oilrate > 0) {
                        rs = gasrate / oilrate;
                    }
                    rs = std::min(rs, rsmax_perf[perf]);
                    b_perf[oilpos] = FluidSystem::oilPvt().inverseFormationVolumeFactor(region_idx, temperature, p_avg, rs);
                } else {
                    b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(region_idx, temperature, p_avg);
                }
            } else {
                b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(region_idx, temperature, p_avg);
            }
        }

        // Surface density.
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            surf_dens_perf[well_.numComponents() * perf  + compIdx] = FluidSystem::referenceDensity( phaseIdx, region_idx );
        }

        // We use cell values for solvent injector
        if constexpr (Indices::enableSolvent) {
            b_perf[well_.numComponents() * perf + Indices::contiSolventEqIdx] = solventInverseFormationVolumeFactor(cell_idx);
            surf_dens_perf[well_.numComponents() * perf + Indices::contiSolventEqIdx] = solventRefDensity(cell_idx);
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void StandardWellConnections<FluidSystem,Indices,Scalar>::
computeProperties(const WellState& well_state,
                  const std::function<Scalar(int,int)>& invB,
                  const std::function<Scalar(int,int)>& mobility,
                  const std::function<Scalar(int)>& solventInverseFormationVolumeFactor,
                  const std::function<Scalar(int)>& solventMobility,
                  const std::vector<Scalar>& b_perf,
                  const std::vector<Scalar>& rsmax_perf,
                  const std::vector<Scalar>& rvmax_perf,
                  const std::vector<Scalar>& rvwmax_perf,
                  const std::vector<Scalar>& surf_dens_perf,
                  DeferredLogger& deferred_logger)
{
    // Compute densities
    const int nperf = well_.numPerfs();
    const int np = well_.numPhases();
    std::vector<double> perfRates(b_perf.size(),0.0);
    const auto& ws = well_state.well(well_.indexOfWell());
    const auto& perf_data = ws.perf_data;
    const auto& perf_rates_state = perf_data.phase_rates;

    for (int perf = 0; perf < nperf; ++perf) {
        for (int comp = 0; comp < np; ++comp) {
            perfRates[perf * well_.numComponents() + comp] =  perf_rates_state[perf * np + well_.ebosCompIdxToFlowCompIdx(comp)];
        }
    }

    if constexpr (Indices::enableSolvent) {
        const auto& solvent_perf_rates_state = perf_data.solvent_rates;
        for (int perf = 0; perf < nperf; ++perf) {
            perfRates[perf * well_.numComponents() + Indices::contiSolventEqIdx] = solvent_perf_rates_state[perf];
        }
    }

    // for producers where all perforations have zero rate we
    // approximate the perforation mixture using the mobility ratio
    // and weight the perforations using the well transmissibility.
    bool all_zero = std::all_of(perfRates.begin(), perfRates.end(),
                                [](double val) { return val == 0.0; });
    const auto& comm = well_.parallelWellInfo().communication();
    if (comm.size() > 1)
    {
        all_zero =  (comm.min(all_zero ? 1 : 0) == 1);
    }

    if (all_zero && well_.isProducer()) {
        double total_tw = 0;
        for (int perf = 0; perf < nperf; ++perf) {
            total_tw += well_.wellIndex()[perf];
        }
        if (comm.size() > 1)
        {
            total_tw = comm.sum(total_tw);
        }
        for (int perf = 0; perf < nperf; ++perf) {
            const int cell_idx = well_.cells()[perf];
            const double well_tw_fraction = well_.wellIndex()[perf] / total_tw;
            double total_mobility = 0.0;
            for (int p = 0; p < np; ++p) {
                int ebosPhaseIdx = well_.flowPhaseToEbosPhaseIdx(p);
                total_mobility += invB(cell_idx, ebosPhaseIdx) * mobility(cell_idx, ebosPhaseIdx);
            }
            if constexpr (Indices::enableSolvent) {
                total_mobility += solventInverseFormationVolumeFactor(cell_idx) * solventMobility(cell_idx);
            }
            for (int p = 0; p < np; ++p) {
                int ebosPhaseIdx = well_.flowPhaseToEbosPhaseIdx(p);
                perfRates[perf * well_.numComponents() + p] = well_tw_fraction * mobility(cell_idx, ebosPhaseIdx) / total_mobility;
            }
            if constexpr (Indices::enableSolvent) {
                perfRates[perf * well_.numComponents() + Indices::contiSolventEqIdx] = well_tw_fraction * solventInverseFormationVolumeFactor(cell_idx) / total_mobility;
            }
        }
    }

    this->computeDensities(perfRates, b_perf, rsmax_perf, rvmax_perf, rvwmax_perf, surf_dens_perf, deferred_logger);
    this->computePressureDelta();
}

#define INSTANCE(...) \
template class StandardWellConnections<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)

// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,1u,0u>)

}
