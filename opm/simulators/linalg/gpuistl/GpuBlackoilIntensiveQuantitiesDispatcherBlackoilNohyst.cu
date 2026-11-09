/*
  Copyright 2026 Equinor ASA

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
/*!
 * \file
 *
 * GPU intensive-quantities dispatcher for three-phase black-oil decks with
 * live oil and dry gas (e.g. SPE1), run by \c flow with \c FlowProblemNohyst.
 */
#include <config.h>

// Tricks for older versions of ROCm that do not support alugrid or dune-fem
// (mirrors test_blackoilintensivequantities_gpu.cu).
#ifdef HAVE_DUNE_ALUGRID
#undef HAVE_DUNE_ALUGRID
#endif
#define HAVE_DUNE_ALUGRID 0
#ifdef HAVE_DUNE_FEM
#undef HAVE_DUNE_FEM
#endif
#define HAVE_DUNE_FEM 0

#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher_impl.hpp>

#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>

#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>

#include <opm/simulators/flow/FlowProblemNohystTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowBlackoilNohystTypeTag.hpp>

namespace Opm::gpuistl {

namespace dispatcher_detail {

struct BlackoilNohystDispatcherTraits
{
    using CpuTag = Properties::TTag::FlowProblemNohyst;
    using GpuTag = Properties::TTag::FlowBlackoilNohystGPU;

    using FluidSystemView = GetPropType<GpuTag, Properties::FluidSystem>;
    using GpuPrimaryVariables = GetPropType<GpuTag, Properties::PrimaryVariables>;
    using GpuIntensiveQuantities = BlackOilIntensiveQuantities<GpuTag>;
    using GpuFlowProblemBuf = GpuFlowProblemBufferFor<GpuTag>;

    static constexpr const char* name = "three-phase black-oil: live oil and dry gas";

    template <class Problem>
    static void validate(const Problem& problem)
    {
        assertNoDiffusionOrDispersion<GpuTag>();

        // The dispatcher assumes one entry per local grid DoF, and the overlap
        // and ghost cells of a parallel run have not been verified.
        if (problem.simulator().vanguard().grid().comm().size() > 1) {
            OPM_THROW(std::logic_error,
                      "GPU intensive quantities for three-phase black-oil are only "
                      "supported on a single MPI process");
        }

        rejectGeoModifier(problem);
        rejectUnsupportedFluidSystem();

        const auto& materialLawManager = problem.materialLawManager();
        if (materialLawManager.threePhaseApproach() != EclMultiplexerApproach::Default) {
            OPM_THROW(std::logic_error,
                      "GPU property evaluation requires the default three-phase relative "
                      "permeability model; STONE1 and STONE2 are not supported");
        }
        rejectUnsupportedSaturationFunctions(materialLawManager);

        rejectUnsupportedMixingControls(problem.simulator().vanguard().schedule());
    }

private:
    static void rejectUnsupportedFluidSystem()
    {
        using FluidSystem = GetPropType<CpuTag, Properties::FluidSystem>;

        if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
            || !FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
            || !FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            OPM_THROW(std::logic_error,
                      "GPU intensive quantities for black-oil require water, oil and gas "
                      "to be active");
        }
        if (FluidSystem::enableVaporizedOil()) {
            OPM_THROW(std::logic_error,
                      "VAPOIL (vaporized oil) is not supported on the GPU");
        }
        if (FluidSystem::enableDissolvedGasInWater()) {
            OPM_THROW(std::logic_error,
                      "DISGASW (gas dissolved in water) is not supported on the GPU");
        }
        if (FluidSystem::enableVaporizedWater()) {
            OPM_THROW(std::logic_error,
                      "VAPWAT (vaporized water) is not supported on the GPU");
        }

        // copy_to_gpu also checks the approaches, but only with the enum values.
        using Approaches = typename FluidSystemView::PvtApproachesType;
        static_assert(Approaches::oil == OilPvtApproach::LiveOil);
        static_assert(Approaches::gas == GasPvtApproach::DryGas);
        static_assert(Approaches::water == WaterPvtApproach::ConstantCompressibilityWater);
        if (FluidSystem::oilPvt().approach() != Approaches::oil) {
            OPM_THROW(std::logic_error,
                      "GPU intensive quantities for black-oil require live oil (PVTO); "
                      "other oil PVT models such as PVDO and PVCDO are not supported "
                      "on the GPU");
        }
        if (FluidSystem::gasPvt().approach() != Approaches::gas) {
            OPM_THROW(std::logic_error,
                      "GPU intensive quantities for black-oil require dry gas (PVDG); "
                      "other gas PVT models such as PVTG are not supported on the GPU");
        }
        if (FluidSystem::waterPvt().approach() != Approaches::water) {
            OPM_THROW(std::logic_error,
                      "GPU intensive quantities for black-oil require water PVT from PVTW; "
                      "other water PVT models are not supported on the GPU");
        }
    }

    /// Plain DRSDT is supported: the dispatcher refreshes the mixing limits on
    /// every call. DRSDTCON needs the convective-mixing quantities, which the
    /// GPU intensive quantities do not compute, and VAPPARS needs the maximum
    /// oil saturation history in a form that has not been verified on the GPU.
    static void rejectUnsupportedMixingControls(const Schedule& schedule)
    {
        for (std::size_t reportStep = 0; reportStep < schedule.size(); ++reportStep) {
            const auto& oilVaporization = schedule[reportStep].oilvap();
            if (oilVaporization.getType() == OilVaporizationProperties::OilVaporization::VAPPARS) {
                OPM_THROW(std::logic_error, "VAPPARS is not supported on the GPU");
            }
            if (oilVaporization.drsdtConvective()) {
                OPM_THROW(std::logic_error,
                          "DRSDTCON (convective DRSDT) is not supported on the GPU");
            }
        }
    }
};

} // namespace dispatcher_detail

template <>
struct GpuDispatcherTraits<Properties::TTag::FlowProblemNohyst>
    : dispatcher_detail::BlackoilNohystDispatcherTraits
{
};

// Used by flow for black-oil decks without HYSTER.
template class GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowProblemNohyst>;

} // namespace Opm::gpuistl
