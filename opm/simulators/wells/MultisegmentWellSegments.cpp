/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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
#include <opm/simulators/wells/MultisegmentWellSegments.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/Schedule/MSW/AICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/material/densead/EvaluationFormat.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
MultisegmentWellSegments(const int numSegments,
                         WellInterfaceGeneric& well)
    : perforations_(numSegments)
    , perforation_depth_diffs_(well.numPerfs(), 0.0)
    , inlets_(well.wellEcl().getSegments().size())
    , depth_diffs_(numSegments, 0.0)
    , densities_(numSegments, 0.0)
    , mass_rates_(numSegments, 0.0)
    , viscosities_(numSegments, 0.0)
    , upwinding_segments_(numSegments, 0)
    , phase_densities_(numSegments, std::vector<EvalWell>(well.numComponents(), 0.0)) // number of phase here?
    , phase_fractions_(numSegments, std::vector<EvalWell>(well.numComponents(), 0.0)) // number of phase here?
    , phase_viscosities_(numSegments, std::vector<EvalWell>(well.numComponents(), 0.0)) // number of phase here?
    , well_(well)
{
    // since we decide to use the WellSegments from the well parser. we can reuse a lot from it.
    // for other facilities needed but not available from parser, we need to process them here

    // initialize the segment_perforations_ and update perforation_segment_depth_diffs_
    const WellConnections& completion_set = well_.wellEcl().getConnections();
    // index of the perforation within wells struct
    // there might be some perforations not active, which causes the number of the perforations in
    // well_ecl_ and wells struct different
    // the current implementation is a temporary solution for now, it should be corrected from the parser
    // side
    int i_perf_wells = 0;
    well.perfDepth().resize(well_.numPerfs(), 0.);
    const auto& segment_set = well_.wellEcl().getSegments();
    for (std::size_t perf = 0; perf < completion_set.size(); ++perf) {
        const Connection& connection = completion_set.get(perf);
        if (connection.state() == Connection::State::OPEN) {
            const int segment_index = segment_set.segmentNumberToIndex(connection.segment());
            if (segment_index == -1) {
                OPM_THROW(std::logic_error,
                          fmt::format("COMPSEGS: Well {} has connection in cell {}, {}, {} "
                                      "without associated segment.", well_.wellEcl().name(),
                                      connection.getI() + 1, connection.getJ() + 1,
                                      connection.getK() + 1));
            }
            perforations_[segment_index].push_back(i_perf_wells);
            well.perfDepth()[i_perf_wells] = connection.depth();
            const double segment_depth = segment_set[segment_index].depth();
            perforation_depth_diffs_[i_perf_wells] = well_.perfDepth()[i_perf_wells] - segment_depth;
            i_perf_wells++;
        }
    }

    // initialize the segment_inlets_
    for (const Segment& segment : segment_set) {
        const int segment_number = segment.segmentNumber();
        const int outlet_segment_number = segment.outletSegment();
        if (outlet_segment_number > 0) {
            const int segment_index = segment_set.segmentNumberToIndex(segment_number);
            const int outlet_segment_index = segment_set.segmentNumberToIndex(outlet_segment_number);
            inlets_[outlet_segment_index].push_back(segment_index);
        }
    }

    // calculating the depth difference between the segment and its oulet_segments
    // for the top segment, we will make its zero unless we find other purpose to use this value
    for (int seg = 1; seg < numSegments; ++seg) {
        const double segment_depth = segment_set[seg].depth();
        const int outlet_segment_number = segment_set[seg].outletSegment();
        const Segment& outlet_segment = segment_set[segment_set.segmentNumberToIndex(outlet_segment_number)];
        const double outlet_depth = outlet_segment.depth();
        depth_diffs_[seg] = segment_depth - outlet_depth;
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
computeFluidProperties(const EvalWell& temperature,
                       const EvalWell& saltConcentration,
                       const PrimaryVariables& primary_variables,
                       int pvt_region_index,
                       DeferredLogger& deferred_logger)
{
    std::vector<double> surf_dens(well_.numComponents());
    // Surface density.
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx)) {
            continue;
        }

        const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
        surf_dens[compIdx] = FluidSystem::referenceDensity( phaseIdx, pvt_region_index);
    }

    for (std::size_t seg = 0; seg < perforations_.size(); ++seg) {
        // the compostion of the components inside wellbore under surface condition
        std::vector<EvalWell> mix_s(well_.numComponents(), 0.0);
        for (int comp_idx = 0; comp_idx < well_.numComponents(); ++comp_idx) {
            mix_s[comp_idx] = primary_variables.surfaceVolumeFraction(seg, comp_idx);
        }

        std::vector<EvalWell> b(well_.numComponents(), 0.0);
        std::vector<EvalWell> visc(well_.numComponents(), 0.0);
        std::vector<EvalWell>& phase_densities = phase_densities_[seg];

        const EvalWell seg_pressure = primary_variables.getSegmentPressure(seg);
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            EvalWell rsw(0.0);
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            b[waterCompIdx] =
                FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rsw, saltConcentration);
            visc[waterCompIdx] =
                FluidSystem::waterPvt().viscosity(pvt_region_index, temperature, seg_pressure, rsw, saltConcentration);
            // TODO: double check here
            // TODO: should not we use phaseIndex here?
            phase_densities[waterCompIdx] = b[waterCompIdx] * surf_dens[waterCompIdx];
        }

        EvalWell rv(0.0);
        EvalWell rvw(0.0);
        // gas phase
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const EvalWell rvmax = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_index, temperature, seg_pressure);
                if (mix_s[oilCompIdx] > 0.0) {
                    if (mix_s[gasCompIdx] > 0.0) {
                        rv = mix_s[oilCompIdx] / mix_s[gasCompIdx];
                    }

                    if (rv > rvmax) {
                        rv = rvmax;
                    }
                    b[gasCompIdx] =
                        FluidSystem::gasPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rv, rvw);
                    visc[gasCompIdx] =
                        FluidSystem::gasPvt().viscosity(pvt_region_index, temperature, seg_pressure, rv, rvw);
                    phase_densities[gasCompIdx] = b[gasCompIdx] * surf_dens[gasCompIdx]
                                                + rv * b[gasCompIdx] * surf_dens[oilCompIdx];
                } else { // no oil exists
                    b[gasCompIdx] =
                        FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                    visc[gasCompIdx] =
                        FluidSystem::gasPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                    phase_densities[gasCompIdx] = b[gasCompIdx] * surf_dens[gasCompIdx];
                }
            } else { // no Liquid phase
                // it is the same with zero mix_s[Oil]
                b[gasCompIdx] =
                    FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                visc[gasCompIdx] =
                    FluidSystem::gasPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
            }
        }

        EvalWell rs(0.0);
        // oil phase
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const EvalWell rsmax = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_index, temperature, seg_pressure);
                if (mix_s[gasCompIdx] > 0.0) {
                    if (mix_s[oilCompIdx] > 0.0) {
                        rs = mix_s[gasCompIdx] / mix_s[oilCompIdx];
                    }

                    if (rs > rsmax) {
                        rs = rsmax;
                    }
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, rs);
                    visc[oilCompIdx] =
                        FluidSystem::oilPvt().viscosity(pvt_region_index, temperature, seg_pressure, rs);
                    phase_densities[oilCompIdx] = b[oilCompIdx] * surf_dens[oilCompIdx]
                                                + rs * b[oilCompIdx] * surf_dens[gasCompIdx];
                } else { // no oil exists
                    b[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                    visc[oilCompIdx] =
                        FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
                    phase_densities[oilCompIdx] = b[oilCompIdx] * surf_dens[oilCompIdx];
                }
            } else { // no Liquid phase
                // it is the same with zero mix_s[Oil]
                b[oilCompIdx] =
                    FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure);
                visc[oilCompIdx] =
                    FluidSystem::oilPvt().saturatedViscosity(pvt_region_index, temperature, seg_pressure);
            }
        }

        phase_viscosities_[seg] = visc;

        std::vector<EvalWell> mix(mix_s);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);

            const EvalWell d = 1.0 - rs * rv;
            if (d <= 0.0) {
                const std::string str =
                    fmt::format("Problematic d value {} obtained for well {} "
                                "during segment density calculations with rs {}, "
                                "rv {} and pressure {}. "
                                "Continue as if no dissolution (rs = 0) and "
                                "vaporization (rv = 0) for this connection.",
                                d, well_.name(), rv, seg_pressure);
                deferred_logger.debug(str);
            } else {
                if (rs > 0.0) {
                    mix[gasCompIdx] = (mix_s[gasCompIdx] - mix_s[oilCompIdx] * rs) / d;
                }
                if (rv > 0.0) {
                    mix[oilCompIdx] = (mix_s[oilCompIdx] - mix_s[gasCompIdx] * rv) / d;
                }
            }
        }

        EvalWell volrat(0.0);
        for (int comp_idx = 0; comp_idx < well_.numComponents(); ++comp_idx) {
            volrat += mix[comp_idx] / b[comp_idx];
        }

        viscosities_[seg] = 0.;
        // calculate the average viscosity
        for (int comp_idx = 0; comp_idx < well_.numComponents(); ++comp_idx) {
            const EvalWell fraction =  mix[comp_idx] / b[comp_idx] / volrat;
            // TODO: a little more work needs to be done to handle the negative fractions here
            phase_fractions_[seg][comp_idx] = fraction; // >= 0.0 ? fraction : 0.0;
            viscosities_[seg] += visc[comp_idx] * phase_fractions_[seg][comp_idx];
        }

        EvalWell density(0.0);
        for (int comp_idx = 0; comp_idx < well_.numComponents(); ++comp_idx) {
            density += surf_dens[comp_idx] * mix_s[comp_idx];
        }
        densities_[seg] = density / volrat;

        // calculate the mass rates
        mass_rates_[seg] = 0.;
        for (int comp_idx = 0; comp_idx < well_.numComponents(); ++comp_idx) {
            const int upwind_seg = upwinding_segments_[seg];
            const EvalWell rate = primary_variables.getSegmentRateUpwinding(seg,
                                                                            upwind_seg,
                                                                            comp_idx);
            mass_rates_[seg] += rate * surf_dens[comp_idx];
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
void MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
updateUpwindingSegments(const PrimaryVariables& primary_variables)
{
    for (std::size_t seg = 0; seg < perforations_.size(); ++seg) {
        // special treatment is needed for segment 0
        if (seg == 0) {
            // we are not supposed to have injecting producers and producing injectors
            assert(!(well_.isProducer() && primary_variables.eval(seg)[primary_variables.WQTotal] > 0.));
            assert(!(well_.isInjector() && primary_variables.eval(seg)[primary_variables.WQTotal] < 0.));
            upwinding_segments_[seg] = seg;
            continue;
        }

        // for other normal segments
        if (primary_variables.eval(seg)[primary_variables.WQTotal] <= 0.) {
            upwinding_segments_[seg] = seg;
        } else {
            const auto& segment_set = well_.wellEcl().getSegments();
            const int outlet_segment_index = segment_set.segmentNumberToIndex(segment_set[seg].outletSegment());
            upwinding_segments_[seg] = outlet_segment_index;
        }
    }
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
getHydroPressureLoss(const int seg,
                     const int seg_density) const
{
    return densities_[seg_density] * well_.gravity() * depth_diffs_[seg];
}

template<class FluidSystem, class Indices, class Scalar>
Scalar MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
getPressureDiffSegPerf(const int seg,
                       const int perf) const
{
    return well_.gravity() * densities_[seg].value() * perforation_depth_diffs_[perf];
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
getSurfaceVolume(const EvalWell& temperature,
                 const EvalWell& saltConcentration,
                 const PrimaryVariables& primary_variables,
                 const int pvt_region_index,
                 const int seg_idx) const
{
    const EvalWell seg_pressure = primary_variables.getSegmentPressure(seg_idx);

    std::vector<EvalWell> mix_s(well_.numComponents(), 0.0);
    for (int comp_idx = 0; comp_idx < well_.numComponents(); ++comp_idx) {
        mix_s[comp_idx] = primary_variables.surfaceVolumeFraction(seg_idx, comp_idx);
    }

    std::vector<EvalWell> b(well_.numComponents(), 0.);
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        EvalWell rsw(0.0);
        b[waterCompIdx] =
            FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index,
                                                                 temperature,
                                                                 seg_pressure,
                                                                 rsw,
                                                                 saltConcentration);
    }

    EvalWell rv(0.0);
    EvalWell rvw(0.0);
    // gas phase
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            EvalWell rvmax = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_index,
                                                                                  temperature,
                                                                                  seg_pressure);
            if (rvmax < 0.0) { // negative rvmax can happen if the seg_pressure is outside the range of the table
                rvmax = 0.0;
            }
            if (mix_s[oilCompIdx] > 0.0) {
                if (mix_s[gasCompIdx] > 0.0) {
                    rv = mix_s[oilCompIdx] / mix_s[gasCompIdx];
                }

                if (rv > rvmax) {
                    rv = rvmax;
                }
                b[gasCompIdx] =
                    FluidSystem::gasPvt().inverseFormationVolumeFactor(pvt_region_index,
                                                                       temperature,
                                                                       seg_pressure,
                                                                       rv,
                                                                       rvw);
            } else { // no oil exists
                b[gasCompIdx] =
                    FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                                temperature,
                                                                                seg_pressure);
            }
        } else { // no Liquid phase
            // it is the same with zero mix_s[Oil]
            b[gasCompIdx] =
                FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                            temperature,
                                                                            seg_pressure);
        }
    }

    EvalWell rs(0.0);
    // oil phase
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            EvalWell rsmax = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_index,
                                                                                 temperature,
                                                                                 seg_pressure);
            if (rsmax < 0.0) { // negative rsmax can happen if the seg_pressure is outside the range of the table
                rsmax = 0.0;
            }
            if (mix_s[gasCompIdx] > 0.0) {
                if (mix_s[oilCompIdx] > 0.0) {
                    rs = mix_s[gasCompIdx] / mix_s[oilCompIdx];
                }
                // std::cout << " rs " << rs.value() << " rsmax " << rsmax.value() << std::endl;

                if (rs > rsmax) {
                    rs = rsmax;
                }
                b[oilCompIdx] =
                    FluidSystem::oilPvt().inverseFormationVolumeFactor(pvt_region_index,
                                                                       temperature,
                                                                       seg_pressure,
                                                                       rs);
            } else { // no oil exists
                b[oilCompIdx] =
                    FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                                temperature,
                                                                                seg_pressure);
            }
        } else { // no gas phase
            // it is the same with zero mix_s[Gas]
            b[oilCompIdx] =
                FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_index,
                                                                            temperature,
                                                                            seg_pressure);
        }
    }

    std::vector<EvalWell> mix(mix_s);
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);

        const EvalWell d = 1.0 - rs * rv;
        if (d <= 0.0 || d > 1.0) {
            const std::string str =
                fmt::format("Problematic d value {} obtained for well {} "
                            "during conversion to surface volume with rs {}, "
                            "rv {} and pressure {}. "
                            "Continue as if no dissolution (rs = 0) and "
                            "vaporization (rv = 0) for this connection.",
                            d, well_.name(), rs, rv, seg_pressure);
            OpmLog::debug(str);
        } else {
            if (rs > 0.0) {
                mix[gasCompIdx] = (mix_s[gasCompIdx] - mix_s[oilCompIdx] * rs) / d;
            }
            if (rv > 0.0) {
                mix[oilCompIdx] = (mix_s[oilCompIdx] - mix_s[gasCompIdx] * rv) / d;
            }
        }
    }

    EvalWell vol_ratio(0.0);
    for (int comp_idx = 0; comp_idx < well_.numComponents(); ++comp_idx) {
        vol_ratio += mix[comp_idx] / b[comp_idx];
    }

    // We increase the segment volume with a factor 10 to stabilize the system.
    const double volume = well_.wellEcl().getSegments()[seg_idx].volume();

    return volume / vol_ratio;
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
getFrictionPressureLoss(const int seg, 
                        const bool extra_reverse_flow_derivatives /*false*/) const
{
    EvalWell mass_rate = mass_rates_[seg];
    const int seg_upwind = upwinding_segments_[seg];
    EvalWell density = densities_[seg_upwind];
    EvalWell visc = viscosities_[seg_upwind];
    // In the reverse flow case, we don't have enough slots for all derivatives, e.g.,
    // upwind pressure and flow. We amend this by a second function call optioin, where
    // only these remaining derivatives are considered.
    // For reference: the pressure equation assumes pressure/flow derivatives are given
    // at segment node while fraction derivatives are given at upwind node.

    if (seg != seg_upwind) {
        if (!extra_reverse_flow_derivatives){
            constexpr int WQTotal = Indices::numEq + PrimaryVariables::WQTotal;
            constexpr int SPres = Indices::numEq + PrimaryVariables::SPres;
            density.setDerivative(WQTotal, 0.0);
            density.setDerivative(SPres, 0.0);
            visc.setDerivative(WQTotal, 0.0);
            visc.setDerivative(SPres, 0.0);
        } else {
            if (PrimaryVariables::has_water){
                constexpr int WFrac = Indices::numEq + PrimaryVariables::WFrac;
                density.setDerivative(WFrac, 0.0);
                visc.setDerivative(WFrac, 0.0);
            }
            if (PrimaryVariables::has_gas){
                constexpr int GFrac = Indices::numEq + PrimaryVariables::GFrac;
                density.setDerivative(GFrac, 0.0);
                visc.setDerivative(GFrac, 0.0);
            }
            mass_rate.clearDerivatives();
        }
    }

    const auto& segment_set = well_.wellEcl().getSegments();
    const int outlet_segment_index = segment_set.segmentNumberToIndex(segment_set[seg].outletSegment());
    const double length = segment_set[seg].totalLength() - segment_set[outlet_segment_index].totalLength();
    assert(length > 0.);
    const double roughness = segment_set[seg].roughness();
    const double area = segment_set[seg].crossArea();
    const double diameter = segment_set[seg].internalDiameter();

    const double sign = mass_rate < 0. ? 1.0 : - 1.0;

    return sign * mswellhelpers::frictionPressureLoss(length, diameter, area, roughness, density, mass_rate, visc);
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
pressureDropSpiralICD(const int seg,
                      const bool extra_reverse_flow_derivatives /*false*/) const
{
    const auto& segment_set = well_.wellEcl().getSegments();
    const SICD& sicd = segment_set[seg].spiralICD();

    const int seg_upwind = upwinding_segments_[seg];
    const std::vector<EvalWell>& phase_fractions = phase_fractions_[seg_upwind];
    const std::vector<EvalWell>& phase_viscosities = phase_viscosities_[seg_upwind];

    EvalWell water_fraction = 0.;
    EvalWell water_viscosity = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        water_fraction = phase_fractions[water_pos];
        water_viscosity = phase_viscosities[water_pos];
    }

    EvalWell oil_fraction = 0.;
    EvalWell oil_viscosity = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const int oil_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        oil_fraction = phase_fractions[oil_pos];
        oil_viscosity = phase_viscosities[oil_pos];
    }

    EvalWell gas_fraction = 0.;
    EvalWell gas_viscosity = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        gas_fraction = phase_fractions[gas_pos];
        gas_viscosity = phase_viscosities[gas_pos];
    }

    EvalWell density = densities_[seg_upwind];
    EvalWell mass_rate = mass_rates_[seg];
    // In the reverse flow case, we don't have enough slots for all derivatives, e.g.,
    // upwind pressure and flow. We amend this by a second function call option, where 
    // only these remaining derivatives are considered.
    // For reference: the pressure equation assumes pressure/flow derivatives are given 
    // at segment node while fraction derivatives are given at upwind node.   
    if (seg != seg_upwind) {
        constexpr int nvar = FluidSystem::numPhases + 1;
        std::vector<bool> zero_mask(nvar, false);
        if (!extra_reverse_flow_derivatives){
            zero_mask[PrimaryVariables::WQTotal] = true;
            zero_mask[PrimaryVariables::SPres] = true;
        } else {
            if constexpr (PrimaryVariables::has_water){
                zero_mask[PrimaryVariables::WFrac] = true;
            }
            if constexpr (PrimaryVariables::has_gas){
                zero_mask[PrimaryVariables::GFrac] = true;
            }
            // mass_rate has no extra derivatives (they are organized as in equations)
            mass_rate.clearDerivatives();
        }
        for (int ii = 0; ii < nvar; ++ii) {
            if (zero_mask[ii]) {
                water_fraction.setDerivative(Indices::numEq + ii, 0.0);
                water_viscosity.setDerivative(Indices::numEq + ii, 0.0);
                oil_fraction.setDerivative(Indices::numEq + ii, 0.0);
                oil_viscosity.setDerivative(Indices::numEq + ii, 0.0);
                gas_fraction.setDerivative(Indices::numEq + ii, 0.0);
                gas_viscosity.setDerivative(Indices::numEq + ii, 0.0);
                density.setDerivative(Indices::numEq + ii, 0.0);
            }
        }
    }

    const EvalWell liquid_fraction = water_fraction + oil_fraction;

    // viscosity contribution from the liquid
    const EvalWell liquid_viscosity_fraction = liquid_fraction < 1.e-30 ? oil_fraction * oil_viscosity + water_fraction * water_viscosity :
            liquid_fraction * mswellhelpers::emulsionViscosity(water_fraction, water_viscosity, oil_fraction, oil_viscosity, sicd);

    const EvalWell mixture_viscosity = liquid_viscosity_fraction + gas_fraction * gas_viscosity;

    const EvalWell reservoir_rate = mass_rate / density;

    const EvalWell reservoir_rate_icd = reservoir_rate * sicd.scalingFactor();

    const double viscosity_cali = sicd.viscosityCalibration();

    using MathTool = MathToolbox<EvalWell>;

    const double density_cali = sicd.densityCalibration();
    // make sure we don't pass negative base to powers
    const EvalWell temp_value1 = density > 0.0 ? MathTool::pow(density / density_cali, 0.75) : 0.0;
    const EvalWell temp_value2 = mixture_viscosity > 0.0 ? MathTool::pow(mixture_viscosity / viscosity_cali, 0.25) : 0.0;

    // formulation before 2016, base_strength is used
    // const double base_strength = sicd.strength() / density_cali;
    // formulation since 2016, strength is used instead
    const double strength = sicd.strength();

    const double sign = reservoir_rate_icd <= 0. ? 1.0 : -1.0;

    return sign * temp_value1 * temp_value2 * strength * reservoir_rate_icd * reservoir_rate_icd;
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
pressureDropAutoICD(const int seg,
                    const UnitSystem& unit_system, 
                    const bool extra_reverse_flow_derivatives /*false*/) const
{
    const auto& segment_set = well_.wellEcl().getSegments();
    const AutoICD& aicd = segment_set[seg].autoICD();

    const int seg_upwind = upwinding_segments_[seg];
    const std::vector<EvalWell>& phase_fractions = phase_fractions_[seg_upwind];
    const std::vector<EvalWell>& phase_viscosities = phase_viscosities_[seg_upwind];
    const std::vector<EvalWell>& phase_densities = phase_densities_[seg_upwind];

    EvalWell water_fraction = 0.;
    EvalWell water_viscosity = 0.;
    EvalWell water_density = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        const int water_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        water_fraction = phase_fractions[water_pos];
        water_viscosity = phase_viscosities[water_pos];
        water_density = phase_densities[water_pos];
    }

    EvalWell oil_fraction = 0.;
    EvalWell oil_viscosity = 0.;
    EvalWell oil_density = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const int oil_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        oil_fraction = phase_fractions[oil_pos];
        oil_viscosity = phase_viscosities[oil_pos];
        oil_density = phase_densities[oil_pos];
    }

    EvalWell gas_fraction = 0.;
    EvalWell gas_viscosity = 0.;
    EvalWell gas_density = 0.;
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const int gas_pos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
        gas_fraction = phase_fractions[gas_pos];
        gas_viscosity = phase_viscosities[gas_pos];
        gas_density = phase_densities[gas_pos];
    }

    EvalWell density = densities_[seg_upwind];
    EvalWell mass_rate = mass_rates_[seg];
    // In the reverse flow case, we don't have enough slots for all derivatives, e.g.,
    // upwind pressure and flow. We amend this by a second function call option, where 
    // only these remaining derivatives are considered.
    // For reference: the pressure equation assumes pressure/flow derivatives are given 
    // at segment node while fraction derivatives are given at upwind node.   
    if (seg != seg_upwind) {
        constexpr int nvar = FluidSystem::numPhases + 1;
        std::vector<bool> zero_mask(nvar, false);
        if (!extra_reverse_flow_derivatives){
            zero_mask[PrimaryVariables::WQTotal] = true;
            zero_mask[PrimaryVariables::SPres] = true;
        } else {
            if (PrimaryVariables::has_water){
                zero_mask[PrimaryVariables::WFrac] = true;
            }
            if (PrimaryVariables::has_gas){
                zero_mask[PrimaryVariables::GFrac] = true;
            }
            // mass_rate has no extra derivatives (they are organized as in equations)
            mass_rate.clearDerivatives();
        }
        for (int ii = 0; ii < nvar; ++ii) {
            if (zero_mask[ii]) {
                water_fraction.setDerivative(Indices::numEq + ii, 0.0);
                water_viscosity.setDerivative(Indices::numEq + ii, 0.0);
                water_density.setDerivative(Indices::numEq + ii, 0.0);
                oil_fraction.setDerivative(Indices::numEq + ii, 0.0);
                oil_viscosity.setDerivative(Indices::numEq + ii, 0.0);
                oil_density.setDerivative(Indices::numEq + ii, 0.0);
                gas_fraction.setDerivative(Indices::numEq + ii, 0.0);
                gas_viscosity.setDerivative(Indices::numEq + ii, 0.0);
                gas_density.setDerivative(Indices::numEq + ii, 0.0);
                density.setDerivative(Indices::numEq + ii, 0.0);
            }
        }
    }

    using MathTool = MathToolbox<EvalWell>;
    // make sure we don't pass negative base to powers
    auto safe_pow = [](const auto& a, const double b) {
        return a > 0.0 ? MathTool::pow(a,b) : 0.0;
    };

    const EvalWell mixture_viscosity = safe_pow(water_fraction, aicd.waterViscExponent()) * water_viscosity
                                     + safe_pow(oil_fraction, aicd.oilViscExponent()) * oil_viscosity
                                     + safe_pow(gas_fraction, aicd.gasViscExponent()) * gas_viscosity;

    const EvalWell mixture_density = safe_pow(water_fraction, aicd.waterDensityExponent()) * water_density
                                   + safe_pow(oil_fraction, aicd.oilDensityExponent()) * oil_density
                                   + safe_pow(gas_fraction, aicd.gasDensityExponent()) * gas_density;

    const double rho_reference = aicd.densityCalibration();
    const double visc_reference = aicd.viscosityCalibration();
    const auto volume_rate_icd = mass_rate * aicd.scalingFactor() / mixture_density;
    const double sign = volume_rate_icd <= 0. ? 1.0 : -1.0;
    // convert 1 unit volume rate
    using M  = UnitSystem::measure;
    const double unit_volume_rate = unit_system.to_si(M::geometric_volume_rate, 1.);

    // TODO: we did not consider the maximum allowed rate here
    const auto result = sign / rho_reference * mixture_density * mixture_density
                      * safe_pow(visc_reference/mixture_viscosity, aicd.viscExponent())
                      * aicd.strength() * safe_pow( -sign * volume_rate_icd, aicd.flowRateExponent())
                      * std::pow(unit_volume_rate, (2. - aicd.flowRateExponent())) ;
    return result;
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
pressureDropValve(const int seg,
                  const SummaryState& summary_state,
                  const bool extra_reverse_flow_derivatives /*false*/) const
{
    const Opm::WellSegments& segment_set = well_.wellEcl().getSegments();
    const Valve& valve = segment_set[seg].valve();

    EvalWell mass_rate = mass_rates_[seg];
    const int seg_upwind = upwinding_segments_[seg];
    EvalWell visc = viscosities_[seg_upwind];
    EvalWell density = densities_[seg_upwind];
    // In the reverse flow case, we don't have enough slots for all derivatives, e.g.,
    // upwind pressure and flow. We amend this by a second function call optioin, where 
    // only these remaining derivatives are considered.
    // For reference: the pressure equation assumes pressure/flow derivatives are given 
    // at segment node while fraction derivatives are given at upwind node. 
    if (seg != seg_upwind) {
        if (!extra_reverse_flow_derivatives){
            constexpr int WQTotal = Indices::numEq + PrimaryVariables::WQTotal;
            constexpr int SPres = Indices::numEq + PrimaryVariables::SPres;
            density.setDerivative(WQTotal, 0.0);
            density.setDerivative(SPres, 0.0);
            visc.setDerivative(WQTotal, 0.0);
            visc.setDerivative(SPres, 0.0);
        } else {
            if (PrimaryVariables::has_water){
                constexpr int WFrac = Indices::numEq + PrimaryVariables::WFrac;
                density.setDerivative(WFrac, 0.0);
                visc.setDerivative(WFrac, 0.0);
            }
            if (PrimaryVariables::has_gas){
                constexpr int GFrac = Indices::numEq + PrimaryVariables::GFrac;
                density.setDerivative(GFrac, 0.0);
                visc.setDerivative(GFrac, 0.0);
            }
            mass_rate.clearDerivatives();
        }
    }

    const double additional_length = valve.pipeAdditionalLength();
    const double roughness = valve.pipeRoughness();
    const double diameter = valve.pipeDiameter();
    const double area = valve.pipeCrossArea();

    const EvalWell friction_pressure_loss =
        mswellhelpers::frictionPressureLoss(additional_length, diameter, area, roughness, density, mass_rate, visc);

    const ValveUDAEval uda_eval {summary_state, this->well_.name(), static_cast<std::size_t>(segment_set[seg].segmentNumber())};
    const double area_con = valve.conCrossArea(uda_eval);
    const double cv = valve.conFlowCoefficient();

    const EvalWell constriction_pressure_loss =
        mswellhelpers::valveContrictionPressureLoss(mass_rate, density, area_con, cv);

    const double sign = mass_rate <= 0. ? 1.0 : -1.0;
    return sign * (friction_pressure_loss + constriction_pressure_loss);
}

template<class FluidSystem, class Indices, class Scalar>
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
accelerationPressureLossContribution(const int seg,
                                     const double area,
                                     const bool extra_reverse_flow_derivatives /*false*/) const
{
    // Compute the *signed* velocity head for given segment (sign is positive for flow towards surface, i.e., negative rate) 
    // Optionally return derivatives for reversed flow case
    EvalWell mass_rate = mass_rates_[seg];
    const int seg_upwind = upwinding_segments_[seg];
    EvalWell density = densities_[seg_upwind];    
    if (seg != seg_upwind) {
        if (!extra_reverse_flow_derivatives){
            constexpr int WQTotal = Indices::numEq + PrimaryVariables::WQTotal;
            constexpr int SPres = Indices::numEq + PrimaryVariables::SPres;
            density.setDerivative(WQTotal, 0.0);
            density.setDerivative(SPres, 0.0);
        } else {
            if (PrimaryVariables::has_water){
                constexpr int WFrac = Indices::numEq + PrimaryVariables::WFrac;
                density.setDerivative(WFrac, 0.0);
            }
            if (PrimaryVariables::has_gas){
                constexpr int GFrac = Indices::numEq + PrimaryVariables::GFrac;
                density.setDerivative(GFrac, 0.0);
            }
            mass_rate.clearDerivatives();
        }
    }
    const double sign = mass_rate > 0 ? -1.0 : 1.0;
    return sign*mswellhelpers::velocityHead(area, mass_rate, density);
}                                     

template <class FluidSystem, class Indices, class Scalar>
void
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
copyPhaseDensities(const PhaseUsage& pu, SegmentState& segSol) const
{
    auto* rho = segSol.phase_density.data();

    const auto phaseMap = std::vector {
        std::pair { BlackoilPhases::Liquid, FluidSystem::oilPhaseIdx },
        std::pair { BlackoilPhases::Vapour, FluidSystem::gasPhaseIdx },
        std::pair { BlackoilPhases::Aqua  , FluidSystem::waterPhaseIdx },
    };

    // Densities stored in 'rho' as
    // [{ p0, p1, ..., (np - 1), mixture, mixture_with_exponents },
    //  { p0, p1, ..., (np - 1), mixture, mixture_with_exponents },
    //  ...
    //  { p0, p1, ..., (np - 1), mixture, mixture_with_exponents }]
    // Stride is np + 2.
    for (const auto& [boPhase, fsPhaseIdx] : phaseMap) {
        if (pu.phase_used[boPhase]) {
            this->copyPhaseDensities(fsPhaseIdx, pu.num_phases + 2,
                                     rho + pu.phase_pos[boPhase]);
        }
    }

    // Mixture densities.
    for (auto seg = 0*this->densities_.size(); seg < this->densities_.size(); ++seg) {
        const auto mixOffset = seg*(pu.num_phases + 2) + pu.num_phases;

        rho[mixOffset + 0] = this->mixtureDensity(seg);
        rho[mixOffset + 1] = this->mixtureDensityWithExponents(seg);
    }
}

template <class FluidSystem, class Indices, class Scalar>
void
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
copyPhaseDensities(const unsigned    phaseIdx,
                   const std::size_t stride,
                   double*           dens) const
{
    const auto compIdx = Indices::canonicalToActiveComponentIndex
        (FluidSystem::solventComponentIndex(phaseIdx));

    for (const auto& phase_density : this->phase_densities_) {
        *dens = phase_density[compIdx].value();
        dens += stride;
    }
}

template <class FluidSystem, class Indices, class Scalar>
double
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
mixtureDensity(const int seg) const
{
    auto mixDens = 0.0;

    const auto& rho = this->phase_densities_[seg];
    const auto& q   = this->phase_fractions_[seg];

    for (const auto& phIdx : {
            FluidSystem::oilPhaseIdx,
            FluidSystem::gasPhaseIdx,
            FluidSystem::waterPhaseIdx
        })
    {
        if (! FluidSystem::phaseIsActive(phIdx)) {
            continue;
        }

        const auto compIdx = Indices::
            canonicalToActiveComponentIndex(phIdx);

        mixDens += q[compIdx].value() * rho[compIdx].value();
    }

    return mixDens;
}

template <class FluidSystem, class Indices, class Scalar>
double
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
mixtureDensityWithExponents(const int seg) const
{
    if (const auto& segment = this->well_.wellEcl().getSegments()[seg];
        segment.isAICD())
    {
        return this->mixtureDensityWithExponents(segment.autoICD(), seg);
    }

    // No other segment type includes exponents of flowing fractions when
    // calculating the mixture/emulsion density.
    return this->mixtureDensity(seg);
}

template <class FluidSystem, class Indices, class Scalar>
double
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
mixtureDensityWithExponents(const AutoICD& aicd, const int seg) const
{
    auto mixDens = 0.0;

    const auto& rho = this->phase_densities_[seg];
    const auto& q   = this->phase_fractions_[seg];

    constexpr auto densityExponents = std::array {
        std::pair { FluidSystem::oilPhaseIdx  , &AutoICD::oilDensityExponent   },
        std::pair { FluidSystem::gasPhaseIdx  , &AutoICD::gasDensityExponent   },
        std::pair { FluidSystem::waterPhaseIdx, &AutoICD::waterDensityExponent },
    };

    for (const auto& [fsPhaseIdx, densityExponent] : densityExponents) {
        if (FluidSystem::phaseIsActive(fsPhaseIdx)) {
            const auto compIdx = Indices::
                canonicalToActiveComponentIndex(fsPhaseIdx);

            // exp = (aicd.*densityExponent)() in native syntax.
            const auto exp = std::invoke(densityExponent, aicd);

            mixDens += std::pow(q[compIdx].value(), exp) * rho[compIdx].value();
        }
    }

    return mixDens;
}

#define INSTANCE(...) \
template class MultisegmentWellSegments<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>,__VA_ARGS__,double>;

// One phase
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>)

// Two phase
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)
INSTANCE(BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>)
// Blackoil
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,true,2u,0u>)
INSTANCE(BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,0u,false,false,1u,0u>)
INSTANCE(BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)

INSTANCE(BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

} // namespace Opm
