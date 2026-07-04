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


#ifndef OPM_MULTISEGMENTWELL_SEGMENTS_HEADER_INCLUDED
#define OPM_MULTISEGMENTWELL_SEGMENTS_HEADER_INCLUDED

#include <opm/material/densead/EvaluationFormat.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace Opm {

    class  AutoICD;
    template<class Scalar> class SegmentState;
    class  UnitSystem;
    template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
    class  SummaryState;

} // namespace Opm

namespace Opm {

template<typename FluidSystem, typename Indices>
class MultisegmentWellSegments
{
    using PrimaryVariables = MultisegmentWellPrimaryVariables<FluidSystem,Indices>;
    using Scalar = typename FluidSystem::Scalar;
    using EvalWell = typename PrimaryVariables::EvalWell;
    using IndexTraits = typename FluidSystem::IndexTraitsType;

public:
    MultisegmentWellSegments(const int numSegments,
                             const ParallelWellInfo<Scalar>& parallel_well_info,
                             WellInterfaceGeneric<Scalar, IndexTraits>& well);

    //! \brief Compute per-segment fluid properties from the pre-computed segment fluid states.
    //!
    //! A BlackOilFluidState is maintained for each segment (see MultisegmentWell::
    //! updateSegmentFluidState) for all cases, thermal or not. Its stored inverse formation
    //! volume factors and dissolution/vaporization ratios are reused whenever they coincide
    //! with the values needed here (both oil and gas present in the segment, which is the
    //! common case); only the phase viscosities and a few single-phase corner cases require a
    //! fresh PVT evaluation. Also refreshes the cached segment volume ratios.
    template<class FluidState>
    void computeFluidProperties(const std::vector<FluidState>& segment_fluid_states,
                                const PrimaryVariables& primary_variables,
                                DeferredLogger& deferred_logger);

    //! \brief Refresh the cached volume ratio of a segment from its fluid state.
    template<class FluidState>
    void updateVolumeRatio(const int seg,
                           const FluidState& fluid_state,
                           const PrimaryVariables& primary_variables,
                           DeferredLogger& deferred_logger);

    //! \brief Update upwinding segments.
    void updateUpwindingSegments(const PrimaryVariables& primary_variables);

    EvalWell getHydroPressureLoss(const int seg,
                                  const int seg_side) const;

    //! Pressure difference between segment and perforation.
    Scalar getPressureDiffSegLocalPerf(const int seg,
                                       const int local_perf_index) const;

    EvalWell getFrictionPressureLoss(const int seg,
                                     const bool extra_reverse_flow_derivatives = false) const;

    // pressure drop for Spiral ICD segment (WSEGSICD)
    EvalWell pressureDropSpiralICD(const int seg,
                                   const bool extra_reverse_flow_derivatives = false) const;

    // pressure drop for Autonomous ICD segment (WSEGAICD)
    EvalWell pressureDropAutoICD(const int seg,
                                 const UnitSystem& unit_system,
                                 const bool extra_reverse_flow_derivatives = false) const;

    // pressure drop for sub-critical valve (WSEGVALV)
    EvalWell pressureDropValve(const int seg,
                               const SummaryState& st,
                               const bool extra_reverse_flow_derivatives = false) const;

    // pressure loss contribution due to acceleration
    EvalWell accelerationPressureLossContribution(const int seg,
                                                  const Scalar area,
                                                  const bool extra_reverse_flow_derivatives = false) const;

    const std::vector<std::vector<int>>& inlets() const
    {
        return inlets_;
    }

    const std::vector<int>& inlets(const int seg) const
    {
        return inlets_[seg];
    }

    const std::vector<std::vector<int>>& perforations() const
    {
        return perforations_;
    }

    int upwinding_segment(const int seg) const
    {
        return upwinding_segments_[seg];
    }

    Scalar getRefDensity() const
    {
        return densities_[0].value();
    }

    const EvalWell& density(const int seg) const
    {
        return densities_[seg];
    }

    //! \brief Segment volume ratio (reservoir volume per unit surface volume).
    //!
    //! Populated via updateVolumeRatio() so that getSegmentSurfaceVolume() can reuse it
    //! instead of recomputing the PVT quantities.
    const EvalWell& volumeRatio(const int seg) const
    {
        return volume_ratios_[seg];
    }

    Scalar local_perforation_depth_diff(const int local_perf_index) const
    {
        return local_perforation_depth_diffs_[local_perf_index];
    }

    void copyPhaseDensities(SegmentState<Scalar>& segSol) const;

private:
    // TODO: trying to use the information from the Well opm-parser as much
    // as possible, it will possibly be re-implemented later for efficiency reason.

    // the completions that is related to each segment
    // the completions's ids are their index in the vector well_index_, well_cell_
    // This is also assuming the order of the completions in Well is the same with
    // the order of the completions in wells.
    // it is for convenience reason. we can just calculate the information for segment once
    // then using it for all the perforations belonging to this segment
    std::vector<std::vector<int>> perforations_;

    // depth difference between the segment and the perforation
    // or in another way, the depth difference between the perforation and
    // the segment the perforation belongs to
    // This vector contains the depth differences for *all* perforations across all processes
    // that this well lies on, its size is well.wellEcl().getConnections().size(),
    // also it works with *global* perforation indices!
    std::vector<Scalar> local_perforation_depth_diffs_;

    // the inlet segments for each segment. It is for convenience and efficiency reason
    std::vector<std::vector<int>> inlets_;

    std::vector<Scalar> depth_diffs_;

    std::vector<Scalar> surface_densities_;

    // the densities of segment fluids
    // we should not have this member variable
    std::vector<EvalWell> densities_;

    // the segment volume ratio (reservoir volume per unit surface volume); taken from the
    // segment fluid state (see setVolumeRatio)
    std::vector<EvalWell> volume_ratios_;

    // the mass rate of the segments
    std::vector<EvalWell> mass_rates_;

    // the viscosity of the segments
    std::vector<EvalWell> viscosities_;

    // the upwinding segment for each segment based on the flow direction
    std::vector<int> upwinding_segments_;

    std::vector<std::vector<EvalWell>> phase_densities_;
    std::vector<std::vector<EvalWell>> phase_fractions_;
    std::vector<std::vector<EvalWell>> phase_viscosities_;

    WellInterfaceGeneric<Scalar, IndexTraits>& well_;

    void copyPhaseDensities(const unsigned    phaseIdx,
                            const std::size_t stride,
                            Scalar*           dens) const;

    Scalar mixtureDensity(const int seg) const;
    Scalar mixtureDensityWithExponents(const int seg) const;
    Scalar mixtureDensityWithExponents(const AutoICD& aicd, const int seg) const;

    // Per-segment phase quantities needed for the segment fluid properties: the inverse
    // formation volume factors and dissolution/vaporization ratios of the phases actually
    // present in the segment (a phase that is absent uses the saturated values of the other),
    // the reservoir-condition mixture composition and the volume ratio derived from them.
    struct PhaseState
    {
        explicit PhaseState(const std::size_t num_quantities)
            : b(num_quantities, 0.0)
            , mix(num_quantities, 0.0)
        {}

        std::vector<EvalWell> b;
        std::vector<EvalWell> mix;
        EvalWell rs{0.};
        EvalWell rv{0.};
        EvalWell vol_ratio{0.};
    };

    //! \brief Fill @p state from the segment fluid state and surface composition @p mix_s.
    //!
    //! The stored inverse formation volume factors and Rs/Rv of @p fluid_state are reused
    //! when both oil and gas are present in the segment (they coincide with the values
    //! needed here); the single-phase corner cases require their own PVT evaluations since
    //! the fluid state treats an absent phase differently than the property calculation
    //! (which uses the saturated values of the present phase).
    template<class FluidState>
    void phaseState(const FluidState& fluid_state,
                    const std::vector<EvalWell>& mix_s,
                    PhaseState& state,
                    DeferredLogger& deferred_logger) const;
};

template<typename FluidSystem, typename Indices>
template<class FluidState>
void MultisegmentWellSegments<FluidSystem,Indices>::
phaseState(const FluidState& fluid_state,
           const std::vector<EvalWell>& mix_s,
           PhaseState& state,
           DeferredLogger& deferred_logger) const
{
    const int pvt_region_index = well_.pvtRegionIdx();

    const bool waterActive = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
    const bool gasActive = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    const bool oilActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);

    const int waterActiveCompIdx = waterActive ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx) : -1;
    const int gasActiveCompIdx = gasActive ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx) : -1;
    const int oilActiveCompIdx = oilActive ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx) : -1;

    // Pressure and temperature are the same for all phases in the wellbore.
    const unsigned pressure_phase = waterActive ? FluidSystem::waterPhaseIdx
                                   : oilActive ? FluidSystem::oilPhaseIdx
                                               : FluidSystem::gasPhaseIdx;
    const EvalWell temperature = fluid_state.temperature(pressure_phase);
    const EvalWell seg_pressure = fluid_state.pressure(pressure_phase);

    state.rs = 0.;
    state.rv = 0.;
    state.vol_ratio = 0.;

    // water phase; the stored inverse formation volume factor was evaluated exactly this way
    if (waterActive) {
        state.b[waterActiveCompIdx] = fluid_state.invB(FluidSystem::waterPhaseIdx);
    }

    // gas phase
    if (gasActive) {
        const bool oil_exist = oilActive && mix_s[oilActiveCompIdx] > 0.0;
        if (oil_exist) {
            if (mix_s[gasActiveCompIdx] > 0.0) {
                // with free gas present the stored Rv is the capped (possibly saturated)
                // vaporization ratio needed here, and the stored inverse formation volume
                // factor was evaluated with it on the undersaturated table
                state.rv = fluid_state.Rv();
                state.b[gasActiveCompIdx] = fluid_state.invB(FluidSystem::gasPhaseIdx);
            } else {
                // dry gas; the fluid state stores the saturated Rv here, so evaluate anew
                const EvalWell rvw{0.};
                state.b[gasActiveCompIdx] = FluidSystem::gasPvt().inverseFormationVolumeFactor(
                        pvt_region_index, temperature, seg_pressure, state.rv, rvw);
            }
        } else { // no oil here
            state.b[gasActiveCompIdx] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(
                    pvt_region_index, temperature, seg_pressure);
        }
    }

    // oil phase
    if (oilActive) {
        const bool gas_exist = gasActive && mix_s[gasActiveCompIdx] > 0.0;
        if (gas_exist) {
            if (mix_s[oilActiveCompIdx] > 0.0) {
                state.rs = fluid_state.Rs();
                state.b[oilActiveCompIdx] = fluid_state.invB(FluidSystem::oilPhaseIdx);
            } else {
                // dead oil; the fluid state stores the saturated Rs here, so evaluate anew
                state.b[oilActiveCompIdx] = FluidSystem::oilPvt().inverseFormationVolumeFactor(
                        pvt_region_index, temperature, seg_pressure, state.rs);
            }
        } else { // no gas phase
            state.b[oilActiveCompIdx] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(
                    pvt_region_index, temperature, seg_pressure);
        }
    }

    state.mix = mix_s;
    if (oilActive && gasActive) {
        const EvalWell d = 1.0 - state.rs * state.rv;
        if (d <= 0.0) {
            const std::string str =
                fmt::format("Problematic d value {} obtained for well {} during segment density calculations "
                            "with rs {}, rv {} and pressure {}. Continue as if no dissolution (rs = 0) and "
                            "vaporization (rv = 0) for this connection.",
                            d, well_.name(), state.rs, state.rv, seg_pressure);
            deferred_logger.debug(str);
        } else {
            if (state.rs > 0.0) {
                state.mix[gasActiveCompIdx] = (mix_s[gasActiveCompIdx] - mix_s[oilActiveCompIdx] * state.rs) / d;
            }
            if (state.rv > 0.0) {
                state.mix[oilActiveCompIdx] = (mix_s[oilActiveCompIdx] - mix_s[gasActiveCompIdx] * state.rv) / d;
            }
        }
    }

    const int num_quantities = well_.numConservationQuantities();
    for (int comp_idx = 0; comp_idx < num_quantities; ++comp_idx) {
        state.vol_ratio += state.mix[comp_idx] / state.b[comp_idx];
    }
}

template<typename FluidSystem, typename Indices>
template<class FluidState>
void MultisegmentWellSegments<FluidSystem,Indices>::
updateVolumeRatio(const int seg,
                  const FluidState& fluid_state,
                  const PrimaryVariables& primary_variables,
                  DeferredLogger& deferred_logger)
{
    const int num_quantities = well_.numConservationQuantities();
    std::vector<EvalWell> mix_s(num_quantities, 0.0);
    for (int comp_idx = 0; comp_idx < num_quantities; ++comp_idx) {
        mix_s[comp_idx] = primary_variables.surfaceVolumeFraction(seg, comp_idx);
    }

    PhaseState state(num_quantities);
    phaseState(fluid_state, mix_s, state, deferred_logger);
    volume_ratios_[seg] = state.vol_ratio;
}

template<typename FluidSystem, typename Indices>
template<class FluidState>
void MultisegmentWellSegments<FluidSystem,Indices>::
computeFluidProperties(const std::vector<FluidState>& segment_fluid_states,
                       const PrimaryVariables& primary_variables,
                       DeferredLogger& deferred_logger)
{
    const int num_quantities = well_.numConservationQuantities();
    const int pvt_region_index = well_.pvtRegionIdx();

    const bool waterActive = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
    const bool gasActive = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    const bool oilActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);

    const int waterActiveCompIdx = waterActive ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx) : -1;
    const int gasActiveCompIdx = gasActive ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx) : -1;
    const int oilActiveCompIdx = oilActive ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx) : -1;

    std::vector<EvalWell> mix_s(num_quantities, 0.0);
    PhaseState state(num_quantities);

    for (std::size_t seg = 0; seg < perforations_.size(); ++seg) {
        const auto& fs = segment_fluid_states[seg];

        // Pressure and temperature are the same for all phases in the wellbore; the salt
        // concentration is zero when brine is inactive (the water PVT ignores it then).
        const unsigned pressure_phase = waterActive ? FluidSystem::waterPhaseIdx
                                       : oilActive ? FluidSystem::oilPhaseIdx
                                                   : FluidSystem::gasPhaseIdx;
        const EvalWell temperature = fs.temperature(pressure_phase);
        const EvalWell seg_pressure = fs.pressure(pressure_phase);
        const EvalWell saltConcentration = fs.saltConcentration();

        for (int comp_idx = 0; comp_idx < num_quantities; ++comp_idx) {
            mix_s[comp_idx] = primary_variables.surfaceVolumeFraction(seg, comp_idx);
        }

        phaseState(fs, mix_s, state, deferred_logger);
        volume_ratios_[seg] = state.vol_ratio;

        std::ranges::fill(phase_densities_[seg], 0.0);
        std::ranges::fill(phase_viscosities_[seg], 0.0);

        // water phase
        if (waterActive) {
            // rsw is only for interface usage
            const EvalWell rsw{0.};
            phase_viscosities_[seg][waterActiveCompIdx] = FluidSystem::waterPvt().viscosity(
                    pvt_region_index, temperature, seg_pressure, rsw, saltConcentration);
            phase_densities_[seg][waterActiveCompIdx] =
                    state.b[waterActiveCompIdx] * surface_densities_[waterActiveCompIdx];
        }

        // gas phase
        if (gasActive) {
            // rvw is only for interface usage
            const EvalWell rvw{0.};
            const bool oil_exist = oilActive && mix_s[oilActiveCompIdx] > 0.0;
            const EvalWell& b_g = state.b[gasActiveCompIdx];
            if (oil_exist) {
                phase_viscosities_[seg][gasActiveCompIdx] = FluidSystem::gasPvt().viscosity(
                        pvt_region_index, temperature, seg_pressure, state.rv, rvw);
                phase_densities_[seg][gasActiveCompIdx] = b_g * surface_densities_[gasActiveCompIdx]
                                                          + state.rv * b_g * surface_densities_[oilActiveCompIdx];
            } else { // no oil here
                phase_viscosities_[seg][gasActiveCompIdx] = FluidSystem::gasPvt().saturatedViscosity(
                        pvt_region_index, temperature, seg_pressure);
                phase_densities_[seg][gasActiveCompIdx] = b_g * surface_densities_[gasActiveCompIdx];
            }
        }

        // oil phase
        if (oilActive) {
            const bool gas_exist = gasActive && mix_s[gasActiveCompIdx] > 0.0;
            const EvalWell& b_o = state.b[oilActiveCompIdx];
            if (gas_exist) {
                phase_viscosities_[seg][oilActiveCompIdx] = FluidSystem::oilPvt().viscosity(
                        pvt_region_index, temperature, seg_pressure, state.rs);
                phase_densities_[seg][oilActiveCompIdx] = b_o * surface_densities_[oilActiveCompIdx]
                                                          + state.rs * b_o * surface_densities_[gasActiveCompIdx];
            } else { // no gas phase
                phase_viscosities_[seg][oilActiveCompIdx] = FluidSystem::oilPvt().saturatedViscosity(
                        pvt_region_index, temperature, seg_pressure);
                phase_densities_[seg][oilActiveCompIdx] = b_o * surface_densities_[oilActiveCompIdx];
            }
        }

        viscosities_[seg] = 0.;
        // calculate the average viscosity
        for (int comp_idx = 0; comp_idx < num_quantities; ++comp_idx) {
            const EvalWell fraction = state.mix[comp_idx] / state.b[comp_idx] / state.vol_ratio;
            // TODO: a little more work needs to be done to handle the negative fractions here
            phase_fractions_[seg][comp_idx] = fraction; // >= 0.0 ? fraction : 0.0;
            viscosities_[seg] += phase_viscosities_[seg][comp_idx] * phase_fractions_[seg][comp_idx];
        }

        EvalWell density(0.0);
        for (int comp_idx = 0; comp_idx < num_quantities; ++comp_idx) {
            density += surface_densities_[comp_idx] * mix_s[comp_idx];
        }
        densities_[seg] = density / state.vol_ratio;

        // calculate the mass rates
        mass_rates_[seg] = 0.;
        for (int comp_idx = 0; comp_idx < num_quantities; ++comp_idx) {
            const int upwind_seg = upwinding_segments_[seg];
            const EvalWell rate = primary_variables.getSegmentRateUpwinding(seg,
                                                                            upwind_seg,
                                                                            comp_idx);
            mass_rates_[seg] += rate * surface_densities_[comp_idx];
        }

    }
}

} // namespace Opm

#endif // OPM_MULTISEGMENTWELL_SEGMENTS_HEADER_INCLUDED
