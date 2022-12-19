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

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/MSWellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <fmt/format.h>
#include <sstream>

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
    for (size_t perf = 0; perf < completion_set.size(); ++perf) {
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

    for (size_t seg = 0; seg < perforations_.size(); ++seg) {
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
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            b[waterCompIdx] =
                FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index, temperature, seg_pressure, saltConcentration);
            visc[waterCompIdx] =
                FluidSystem::waterPvt().viscosity(pvt_region_index, temperature, seg_pressure, saltConcentration);
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
                std::ostringstream sstr;
                sstr << "Problematic d value " << d << " obtained for well " << well_.name()
                     << " during segment density calculations with rs " << rs
                     << ", rv " << rv << " and pressure " << seg_pressure
                     << " obtaining d " << d
                     << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                     << " for this connection.";
                deferred_logger.debug(sstr.str());
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
typename MultisegmentWellSegments<FluidSystem,Indices,Scalar>::EvalWell
MultisegmentWellSegments<FluidSystem,Indices,Scalar>::
getHydroPressureLoss(const int seg) const
{
    return densities_[seg] * well_.gravity() * depth_diffs_[seg];
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
        b[waterCompIdx] =
            FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_index,
                                                                 temperature,
                                                                 seg_pressure,
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
            std::ostringstream sstr;
            sstr << "Problematic d value " << d << " obtained for well " << well_.name()
                 << " during conversion to surface volume with rs " << rs
                 << ", rv " << rv << " and pressure " << seg_pressure
                 << " obtaining d " << d
                 << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                 << " for this connection.";
            OpmLog::debug(sstr.str());
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
getFrictionPressureLoss(const int seg) const
{
    const EvalWell mass_rate = mass_rates_[seg];
    const int seg_upwind = upwinding_segments_[seg];
    EvalWell density = densities_[seg_upwind];
    EvalWell visc = viscosities_[seg_upwind];
    // WARNING
    // We disregard the derivatives from the upwind density to make sure derivatives
    // wrt. to different segments dont get mixed.
    if (seg != seg_upwind) {
        density.clearDerivatives();
        visc.clearDerivatives();
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

} // namespace Opm
