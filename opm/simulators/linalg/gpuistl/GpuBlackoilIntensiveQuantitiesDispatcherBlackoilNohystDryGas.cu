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
 * GPU intensive-quantities dispatcher configuration for three-phase black-oil
 * decks with live oil and dry gas (e.g. SPE1), run by \c flow with
 * \c FlowProblemNohyst.
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

#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>

#include <opm/simulators/flow/FlowProblemNohystTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowBlackoilNohystTypeTag.hpp>

namespace Opm::gpuistl {

template <>
struct GpuDispatcherTraits<Properties::TTag::FlowProblemNohyst, LiveOilDryGasPvtApproaches>
{
    using CpuTag = Properties::TTag::FlowProblemNohyst;
    using GpuTag = Properties::TTag::FlowBlackoilNohystGPU;

    using FluidSystemView = GetPropType<GpuTag, Properties::FluidSystem>;
    using GpuPrimaryVariables = GetPropType<GpuTag, Properties::PrimaryVariables>;
    using GpuIntensiveQuantities = BlackOilIntensiveQuantities<GpuTag>;
    using GpuFlowProblemBuf = dispatcher_detail::GpuFlowProblemBufferFor<GpuTag>;

    using Approaches = typename FluidSystemView::PvtApproachesType;
    static_assert(Approaches::oil == OilPvtApproach::LiveOil);
    static_assert(Approaches::gas == GasPvtApproach::DryGas);
    static_assert(Approaches::water == WaterPvtApproach::ConstantCompressibilityWater);

    static constexpr const char* name = "three-phase black-oil: live oil and dry gas";

    template <class Problem>
    static void validate(const Problem& problem)
    {
        using FluidSystem = GetPropType<CpuTag, Properties::FluidSystem>;
        dispatcher_detail::assertNoDiffusionOrDispersion<GpuTag>();
        dispatcher_detail::validateBlackoilNohyst<FluidSystem>(problem);

        // PVDG with VAPOIL is a valid CPU combination, but no deck uses it.
        if (FluidSystem::enableVaporizedOil()) {
            OPM_THROW(std::logic_error,
                      "VAPOIL (vaporized oil) with dry gas (PVDG) is not supported on the GPU");
        }
    }
};

// Used by flow for black-oil decks without HYSTER.
template std::unique_ptr<
    GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowProblemNohyst>::Impl>
GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowProblemNohyst>::
    makeImpl<LiveOilDryGasPvtApproaches>(const Problem&);

} // namespace Opm::gpuistl
