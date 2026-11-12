// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * The GPU side of \c GpuBlackoilIntensiveQuantitiesDispatcher, shared by the
 * per-configuration translation units
 * (\c GpuBlackoilIntensiveQuantitiesDispatcher<Configuration>.cu).
 *
 * Only to be included from \c .cu files. Each of them specializes
 * \c GpuDispatcherTraits for its (CPU TypeTag, PVT approaches) pairs and
 * explicitly instantiates \c makeImpl for them. Keeping one configuration per
 * translation unit keeps the nvcc compile time and memory of each unit to a
 * single \c FlowProblemBlackoil instantiation.
 */
#ifndef OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_IMPL_HPP
#define OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_IMPL_HPP

#include <cuda_runtime.h>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/input/eclipse/Schedule/Events.hpp>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/GpuEclMaterialLawManager.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/material/thermal/GpuEclThermalLawManager.hpp>

#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::gpuistl {

/*!
 * \brief The GPU types and the deck validation of one configuration that the
 *        dispatcher supports, keyed on the CPU \c TypeTag of the simulator and
 *        the PVT approaches of the GPU fluid system.
 *
 * A specialization provides:
 *  - \c CpuTag: the CPU TypeTag that the GPU intensive-quantities prototype is
 *    converted from. Its intensive quantities must have the same module layout
 *    as those of \c GpuTag.
 *  - \c GpuTag: the TypeTag of the intensive quantities on the GPU. Its
 *    \c Problem is the \c GpuView-backed \c GpuFlowProblem.
 *  - \c FluidSystemView: the \c GpuView-backed fluid system of \c GpuTag.
 *  - \c GpuPrimaryVariables, \c GpuIntensiveQuantities.
 *  - \c GpuFlowProblemBuf: the owning (\c GpuBuffer) problem whose view is the
 *    \c Problem of \c GpuTag.
 *  - \c static void validate(const Problem&): rejects decks that the GPU
 *    types do not support, with an error that names the feature. The facade
 *    has already checked that the PVT approaches of the deck are \p PvtApproaches.
 *  - \c static constexpr const char* name: shown in the log at initialization.
 *
 * The specializations live in the \c .cu file that instantiates \c makeImpl
 * for that configuration.
 */
template <class CpuTypeTag, class PvtApproaches>
struct GpuDispatcherTraits;

namespace dispatcher_detail {

template <class ThermalLawManagerView>
struct ThermalLawManagerBuffer;

template <>
struct ThermalLawManagerBuffer<::Opm::NoThermalLawManager>
{
    using type = ::Opm::NoThermalLawManager;
};

template <class Scalar, class FluidSystem>
struct ThermalLawManagerBuffer<
    ::Opm::EclThermalLaw::GpuManager<Scalar, FluidSystem, GpuView, GpuView>>
{
    using type = ::Opm::EclThermalLaw::GpuManager<Scalar, FluidSystem, GpuBuffer, GpuView>;
};

template <class GpuProblemView>
struct GpuFlowProblemBuffer;

/// The owning \c GpuBuffer counterpart of a \c GpuView-backed \c GpuFlowProblem.
template <class Scalar, class MaterialLawManagerView, class ThermalLawManagerView>
struct GpuFlowProblemBuffer<
    ::Opm::GpuFlowProblem<Scalar, MaterialLawManagerView, GpuView, ThermalLawManagerView>>
{
    using MaterialLawManagerBuffer = ::Opm::EclMaterialLaw::GpuManager<
        typename MaterialLawManagerView::Traits,
        typename MaterialLawManagerView::GasOilLaw,
        typename MaterialLawManagerView::OilWaterLaw,
        GpuBuffer,
        typename MaterialLawManagerView::MaterialLaw>;

    using type = ::Opm::GpuFlowProblem<
        Scalar,
        MaterialLawManagerBuffer,
        GpuBuffer,
        typename ThermalLawManagerBuffer<ThermalLawManagerView>::type>;
};

/// The \c GpuBuffer-backed problem whose view is the \c Problem of \p GpuTag.
template <class GpuTag>
using GpuFlowProblemBufferFor =
    typename GpuFlowProblemBuffer<GetPropType<GpuTag, Properties::Problem>>::type;

// ---------------------------------------------------------------------------
// Validation shared by the configurations.
// ---------------------------------------------------------------------------

template <class Problem>
void rejectGeoModifier(const Problem& problem)
{
    const auto& schedule = problem.simulator().vanguard().schedule();
    for (std::size_t reportStep = 0; reportStep < schedule.size(); ++reportStep) {
        if (schedule[reportStep].events().hasEvent(ScheduleEvents::GEO_MODIFIER)) {
            OPM_THROW(std::logic_error,
                      "GPU intensive-quantities evaluation does not support GEO_MODIFIER");
        }
    }
}

/// The GPU material law manager evaluates piecewise-linear tables only, and
/// applies neither endpoint scaling nor hysteresis (OPM/opm-simulators#7644).
template <class MaterialLawManager>
void rejectUnsupportedSaturationFunctions(const MaterialLawManager& materialLawManager)
{
    if (!materialLawManager.satCurveIsAllPiecewiseLinear()) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation requires piecewise-linear saturation tables");
    }
    // enableEndPointScaling() is only set by ENKRVD. ENDSCALE, including the
    // capillary pressure scaling of SWATINIT, is enabled in the per-system
    // endpoint scaling configurations.
    const bool endPointScaling = materialLawManager.enableEndPointScaling()
        || materialLawManager.gasOilConfig().enableSatScaling()
        || materialLawManager.oilWaterConfig().enableSatScaling()
        || materialLawManager.gasWaterConfig().enableSatScaling();
    if (materialLawManager.enableHysteresis()
        || endPointScaling
        || materialLawManager.hasDirectionalRelperms()
        || materialLawManager.hasDirectionalImbnum()) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation does not support hysteresis, endpoint scaling "
                  "(ENDSCALE, ENKRVD, SWATINIT), or directional saturation functions");
    }
}

/// Validation shared by the three-phase black-oil configurations of
/// \c FlowProblemNohyst. The PVT approaches and VAPOIL are checked per
/// configuration.
template <class FluidSystem, class Problem>
void validateBlackoilNohyst(const Problem& problem)
{
    // The dispatcher assumes one entry per local grid DoF, and the overlap
    // and ghost cells of a parallel run have not been verified.
    if (problem.simulator().vanguard().grid().comm().size() > 1) {
        OPM_THROW(std::logic_error,
                  "GPU intensive quantities for three-phase black-oil are only "
                  "supported on a single MPI process");
    }

    rejectGeoModifier(problem);

    if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
        || !FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
        || !FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        OPM_THROW(std::logic_error,
                  "GPU intensive quantities for black-oil require water, oil and gas "
                  "to be active");
    }
    if (FluidSystem::enableDissolvedGasInWater()) {
        OPM_THROW(std::logic_error,
                  "DISGASW (gas dissolved in water) is not supported on the GPU");
    }
    if (FluidSystem::enableVaporizedWater()) {
        OPM_THROW(std::logic_error,
                  "VAPWAT (vaporized water) is not supported on the GPU");
    }

    const auto& materialLawManager = problem.materialLawManager();
    if (materialLawManager.threePhaseApproach() != EclMultiplexerApproach::Default) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation requires the default three-phase relative "
                  "permeability model; STONE1 and STONE2 are not supported");
    }
    rejectUnsupportedSaturationFunctions(materialLawManager);

    // Plain DRSDT is supported: the dispatcher refreshes the mixing limits on
    // every call. DRSDTCON needs the convective-mixing quantities, which the
    // GPU intensive quantities do not compute, and VAPPARS needs the maximum
    // oil saturation history in a form that has not been verified on the GPU.
    const auto& schedule = problem.simulator().vanguard().schedule();
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

/// Whether the kernel also updates the energy quantities. A namespace-scope
/// constant, because nvcc spells kernel template arguments out at global scope.
template <class GpuTag>
inline constexpr bool updateEnergyQuantities
    = getPropValue<GpuTag, Properties::EnergyModuleType>() == EnergyModules::FullyImplicitThermal;

template <class GpuTag>
constexpr void assertNoDiffusionOrDispersion()
{
    // The GPU intensive quantities intentionally do not implement these modules.
    static_assert(!getPropValue<GpuTag, Properties::EnableDiffusion>());
    static_assert(!getPropValue<GpuTag, Properties::EnableDispersion>());
}

/// The GPU problem data is indexed from zero over all grid DoFs. The checks
/// that do not need the problem are done by the facade, before the setup.
template <class ProblemT>
void validateNumDof(const ProblemT& problem, std::size_t numDof)
{
    if (numDof != static_cast<std::size_t>(problem.model().numGridDof())) {
        OPM_THROW(std::invalid_argument,
                  "GPU intensive-quantities dispatcher requires one entry per grid DoF");
    }
}

// The GPU kernel that runs the IQ update for each cell
template <bool updateEnergy, class GpuProblem, class PrimaryVariablesT, class IntensiveQuantitiesT>
__global__ void
dispatcherUpdateAllCellsKernel(GpuProblem problem,
                               GpuView<const PrimaryVariablesT> primaryVariables,
                               GpuView<IntensiveQuantitiesT> outIntensiveQuantities,
                               std::size_t numCells)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numCells) {
        return;
    }
    IntensiveQuantitiesT& iq = outIntensiveQuantities[i];
    iq.update(problem, primaryVariables[i], static_cast<unsigned>(i), 0);
    if constexpr (updateEnergy) {
        iq.updateEnergyQuantities_(problem, static_cast<unsigned>(i), 0u);
    }
}

// =============================================================================
// ImplFor: the GPU state and the kernel launch of one configuration
// =============================================================================
template <class CpuTypeTag, class PvtApproaches>
struct ImplFor final : GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Impl
{
    using Dispatcher = GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>;
    using Problem = typename Dispatcher::Problem;
    using PrimaryVariables = typename Dispatcher::PrimaryVariables;
    using IntensiveQuantities = typename Dispatcher::IntensiveQuantities;

    using Traits = GpuDispatcherTraits<CpuTypeTag, PvtApproaches>;
    using CpuTag = typename Traits::CpuTag;
    using GpuTag = typename Traits::GpuTag;
    using FluidSystemView = typename Traits::FluidSystemView;
    using GpuPrimaryVariables = typename Traits::GpuPrimaryVariables;
    using GpuIntensiveQuantities = typename Traits::GpuIntensiveQuantities;
    using GpuFlowProblemBuf = typename Traits::GpuFlowProblemBuf;
    using GpuFlowProblemView = decltype(::Opm::gpuistl::make_view(std::declval<GpuFlowProblemBuf&>()));

    static_assert(std::is_same_v<GpuFlowProblemView, GetPropType<GpuTag, Properties::Problem>>,
                  "the view of GpuFlowProblemBuf must be the Problem of the GPU TypeTag");
    static_assert(std::is_same_v<FluidSystemView, GetPropType<GpuTag, Properties::FluidSystem>>);
    static_assert(std::is_same_v<typename FluidSystemView::PvtApproachesType, PvtApproaches>,
                  "the GPU fluid system must have the PVT approaches of the configuration");

    using CpuFluidSystem = GetPropType<CpuTypeTag, Properties::FluidSystem>;
    using DynamicCpuFluidSystem
        = std::remove_reference_t<decltype(CpuFluidSystem::getNonStaticInstance())>;
    using FluidSystemBuffer
        = decltype(::Opm::gpuistl::copy_to_gpu<PvtApproaches>(
            std::declval<DynamicCpuFluidSystem&>()));
    using ManagedFluidSystemView
        = std::unique_ptr<FluidSystemView, GpuManagedDeleter<FluidSystemView>>;

    std::unique_ptr<GpuFlowProblemBuf> problemBuf;
    GpuFlowProblemView problemView{};
    std::unique_ptr<FluidSystemBuffer> fluidSystemBuffer;
    ManagedFluidSystemView managedFluidSystemView;
    std::optional<GpuIntensiveQuantities> prototype;
    // The constructor reads the mixing limits, so the first update() need not.
    bool mixingLimitsCurrent = true;

    explicit ImplFor(const Problem& cpuProblem)
    {
        // Build the GPU FlowProblem from the CPU FlowProblem (one-time setup).
        problemBuf = std::make_unique<GpuFlowProblemBuf>(cpuProblem);
        problemView = ::Opm::gpuistl::make_view(*problemBuf);

        // Place the FluidSystemView in unified memory so its device pointer
        // dereferences are valid both on host and device (mirrors the test).
        auto& dynamicCpuFluidSystem = CpuFluidSystem::getNonStaticInstance();
        // Keep the owning buffer with this dispatcher. A process-wide static
        // would retain the first deck's PVT tables across later simulations.
        fluidSystemBuffer = std::make_unique<FluidSystemBuffer>(
            ::Opm::gpuistl::copy_to_gpu<PvtApproaches>(dynamicCpuFluidSystem));
        auto fsView = ::Opm::gpuistl::make_view(*fluidSystemBuffer);

        managedFluidSystemView = make_gpu_managed_unique_ptr<FluidSystemView>(fsView);

        // Build a default IntensiveQuantities prototype for the GPU side.
        BlackOilIntensiveQuantities<CpuTag> cpuPrototype;
        prototype = cpuPrototype.template withOtherFluidSystem<GpuTag>(*managedFluidSystemView);
    }

    const char* name() const override
    {
        return Traits::name;
    }

    void update(const Problem& cpuProblem,
                const PrimaryVariables* const* cpuPriVars,
                IntensiveQuantities* const* outIQ,
                std::size_t numDof,
                unsigned timeIdx) override
    {
        validateNumDof(cpuProblem, numDof);

        if (!mixingLimitsCurrent) {
            // The mixing limits (DRSDT, DRVDT, VAPPARS) change between calls:
            // at the start of a time step, on a chopped time step and on restart.
            // Re-reading them on every call cannot miss any of these.
            problemBuf->refreshMixingLimits(cpuProblem, timeIdx);
        }
        mixingLimitsCurrent = false;

        // -------------------------------------------------------------------
        // 1. Convert CPU primary variables to the dispatcher's GPU primary
        //    variables (host-side).
        // -------------------------------------------------------------------
        std::vector<GpuPrimaryVariables> hostPriVars;
        hostPriVars.reserve(numDof);
        for (std::size_t i = 0; i < numDof; ++i) {
            // The BlackOilPrimaryVariables converting copy ctor handles the
            // TypeTag/storage difference.
            hostPriVars.emplace_back(*cpuPriVars[i]);
        }
        std::vector<GpuIntensiveQuantities> hostIQ(numDof, *prototype);

        // -------------------------------------------------------------------
        // 2. Upload host buffers to the device (host-to-device copies are
        //    embedded in the GpuBuffer ctors).
        // -------------------------------------------------------------------
        GpuBuffer<GpuPrimaryVariables> primaryVariablesBuffer(hostPriVars);
        GpuBuffer<GpuIntensiveQuantities> intensiveQuantitiesBuffer(hostIQ);

        // -------------------------------------------------------------------
        // 3. Launch the per-cell update kernel.
        // -------------------------------------------------------------------
        const unsigned blockSize = 64u;
        const unsigned gridSize = static_cast<unsigned>((numDof + blockSize - 1u) / blockSize);

        dispatcherUpdateAllCellsKernel<updateEnergyQuantities<GpuTag>>
            <<<gridSize, blockSize>>>(
                problemView,
                GpuView<const GpuPrimaryVariables>(
                    primaryVariablesBuffer.data(), primaryVariablesBuffer.size()),
                GpuView<GpuIntensiveQuantities>(
                    intensiveQuantitiesBuffer.data(), intensiveQuantitiesBuffer.size()),
                numDof);
        OPM_GPU_SAFE_CALL(cudaGetLastError());

        // -------------------------------------------------------------------
        // 4. Read the GPU IntensiveQuantities back to host memory.
        // -------------------------------------------------------------------
        intensiveQuantitiesBuffer.copyToHost(hostIQ.data(), numDof);

        // -------------------------------------------------------------------
        // 5. Field-by-field materialization onto the caller's CPU
        //    IntensiveQuantities. The supported GPU TypeTags have no diffusion
        //    or dispersion state, so the overlay covers the complete supported IQ.
        // -------------------------------------------------------------------
        for (std::size_t i = 0; i < numDof; ++i) {
            outIQ[i]->overlayBlackOilFieldsFrom(hostIQ[i]);
        }
    }
};

} // namespace dispatcher_detail

template <class CpuTypeTag>
template <class PvtApproaches>
std::unique_ptr<typename GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Impl>
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::makeImpl(const Problem& cpuProblem)
{
    // The first call comes from the initial solution, so an unsupported
    // deck is rejected before the first time step.
    GpuDispatcherTraits<CpuTypeTag, PvtApproaches>::validate(cpuProblem);
    return std::make_unique<dispatcher_detail::ImplFor<CpuTypeTag, PvtApproaches>>(cpuProblem);
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_IMPL_HPP
