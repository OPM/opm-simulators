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

#include <opm/simulators/wells/MultisegmentWellPrimaryVariables.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <array>
#include <cstddef>
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

    //! \brief The quantities taken from a segment fluid state.
    //!
    //! The fluid state type depends on the TypeTag, which this class does not have, so the
    //! caller extracts the values needed here. Pressure and temperature are the same for all
    //! phases in the wellbore, so a single value of each is stored. @p invB is indexed by
    //! canonical phase index.
    struct SegmentPvt
    {
        EvalWell pressure{0.};
        EvalWell temperature{0.};
        EvalWell saltConcentration{0.};
        std::array<EvalWell, FluidSystem::numPhases> invB{};
        EvalWell Rs{0.};
        EvalWell Rv{0.};
    };

    //! \brief Compute per-segment fluid properties from the pre-computed segment PVT values.
    //!
    //! Reuses the inverse formation volume factors and dissolution/vaporization ratios where
    //! applicable, and evaluates the remaining PVT properties. Also updates the cached segment
    //! volume ratios.
    void computeFluidProperties(const std::vector<SegmentPvt>& segment_pvt,
                                const PrimaryVariables& primary_variables,
                                DeferredLogger& deferred_logger);

    //! \brief Volume ratio of a segment from its PVT values and composition.
    EvalWell computeVolumeRatio(const int seg,
                                const SegmentPvt& pvt,
                                const PrimaryVariables& primary_variables,
                                DeferredLogger& deferred_logger) const;

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
    //! Set by computeFluidProperties(), like the other per-segment quantities.
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

    // Segment volume ratios (reservoir volume per unit surface volume), derived from the fluid
    // state and surface composition.
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

    // Intermediate PVT properties used to calculate segment fluid properties and volume ratios.
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

    //! \brief Fill @p state from the segment PVT values and surface composition @p mix_s.
    //!
    //! Reuses the PVT properties where they represent the requested composition, and
    //! evaluates the required saturated properties for single-phase cases.
    void calculatePhaseState(const SegmentPvt& pvt,
                             const std::vector<EvalWell>& mix_s,
                             PhaseState& state,
                             DeferredLogger& deferred_logger) const;
};

} // namespace Opm

#endif // OPM_MULTISEGMENTWELL_SEGMENTS_HEADER_INCLUDED
