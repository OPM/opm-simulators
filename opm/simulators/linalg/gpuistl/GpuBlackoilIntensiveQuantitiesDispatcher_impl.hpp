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
 * Implementation of \c GpuBlackoilIntensiveQuantitiesDispatcher, shared by the
 * per-configuration translation units
 * (\c GpuBlackoilIntensiveQuantitiesDispatcher<Configuration>.cu).
 *
 * Only to be included from \c .cu files. Each of them specializes
 * \c GpuDispatcherTraits for the CPU TypeTags it supports and explicitly
 * instantiates the dispatcher for them. Keeping one configuration per
 * translation unit keeps the nvcc compile time and memory of each unit to a
 * single \c FlowProblemBlackoil instantiation.
 */
#ifndef OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_IMPL_HPP
#define OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_IMPL_HPP

#include <cuda_runtime.h>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/input/eclipse/Schedule/Events.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

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
#include <format>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::gpuistl {

/*!
 * \brief The GPU types and the deck validation of one configuration that the
 *        dispatcher supports, keyed on the CPU \c TypeTag of the simulator.
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
 *    types do not support, with an error that names the feature.
 *  - \c static constexpr const char* name: shown in the log at initialization.
 *
 * The specializations live in the \c .cu file that instantiates the
 * dispatcher for that configuration.
 */
template <class CpuTypeTag>
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
    if (materialLawManager.enableHysteresis()
        || materialLawManager.enableEndPointScaling()
        || materialLawManager.hasDirectionalRelperms()
        || materialLawManager.hasDirectionalImbnum()) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation does not support hysteresis, endpoint scaling, "
                  "or directional saturation functions");
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

template <class ProblemT, class PrimaryVariablesT, class IntensiveQuantitiesT>
void validateDispatcherInputs(const ProblemT& problem,
                              const PrimaryVariablesT* const* primaryVariables,
                              IntensiveQuantitiesT* const* intensiveQuantities,
                              std::size_t numDof)
{
    if (primaryVariables == nullptr || intensiveQuantities == nullptr) {
        OPM_THROW(std::invalid_argument,
                  "GPU intensive-quantities dispatcher received a null pointer array");
    }

    if (numDof != static_cast<std::size_t>(problem.model().numGridDof())) {
        OPM_THROW(std::invalid_argument,
                  "GPU intensive-quantities dispatcher requires one entry per grid DoF");
    }

    for (std::size_t i = 0; i < numDof; ++i) {
        if (primaryVariables[i] == nullptr || intensiveQuantities[i] == nullptr) {
            OPM_THROW(std::invalid_argument,
                      "GPU intensive-quantities dispatcher received a null DoF entry");
        }
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

} // namespace dispatcher_detail

// =============================================================================
// Impl
// =============================================================================
template <class CpuTypeTag>
struct GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Impl {
    using Traits = GpuDispatcherTraits<CpuTypeTag>;
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

    using CpuFluidSystem = GetPropType<CpuTypeTag, Properties::FluidSystem>;
    using DynamicCpuFluidSystem
        = std::remove_reference_t<decltype(CpuFluidSystem::getNonStaticInstance())>;
    using FluidSystemBuffer
        = decltype(::Opm::gpuistl::copy_to_gpu<typename FluidSystemView::PvtApproachesType>(
            std::declval<DynamicCpuFluidSystem&>()));
    using ManagedFluidSystemView
        = std::unique_ptr<FluidSystemView, GpuManagedDeleter<FluidSystemView>>;

    bool initialized = false;
    std::unique_ptr<GpuFlowProblemBuf> problemBuf;
    GpuFlowProblemView problemView{};
    std::unique_ptr<FluidSystemBuffer> fluidSystemBuffer;
    ManagedFluidSystemView managedFluidSystemView;
    std::optional<GpuIntensiveQuantities> prototype;
};

template <class CpuTypeTag>
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::GpuBlackoilIntensiveQuantitiesDispatcher()
    : impl_(std::make_unique<Impl>())
{
}

template <class CpuTypeTag>
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::~GpuBlackoilIntensiveQuantitiesDispatcher()
    = default;

template <class CpuTypeTag>
void GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::update(
    const Problem& cpuProblem,
    const PrimaryVariables* const* cpuPriVars,
    IntensiveQuantities* const* outIQ,
    std::size_t numDof)
{
    using ImplT = Impl;
    // The kernel's template arguments are named through the public traits:
    // nvcc spells them out at global scope, where Impl is not accessible.
    using Traits = GpuDispatcherTraits<CpuTypeTag>;
    using GpuPrimaryVariables = typename Traits::GpuPrimaryVariables;
    using GpuIntensiveQuantities = typename Traits::GpuIntensiveQuantities;

    if (numDof == 0u) {
        return;
    }

    dispatcher_detail::validateDispatcherInputs(cpuProblem, cpuPriVars, outIQ, numDof);

    if (!impl_->initialized) {
        // The first call comes from the initial solution, so an unsupported
        // deck is rejected before the first time step.
        Traits::validate(cpuProblem);

        // Build the GPU FlowProblem from the CPU FlowProblem (one-time setup).
        impl_->problemBuf = std::make_unique<typename ImplT::GpuFlowProblemBuf>(cpuProblem);
        impl_->problemView = ::Opm::gpuistl::make_view(*impl_->problemBuf);

        // Place the FluidSystemView in unified memory so its device pointer
        // dereferences are valid both on host and device (mirrors the test).
        auto& dynamicCpuFluidSystem = ImplT::CpuFluidSystem::getNonStaticInstance();
        // Keep the owning buffer with this dispatcher. A process-wide static
        // would retain the first deck's PVT tables across later simulations.
        impl_->fluidSystemBuffer
            = std::make_unique<typename ImplT::FluidSystemBuffer>(
                ::Opm::gpuistl::copy_to_gpu<typename ImplT::FluidSystemView::PvtApproachesType>(
                    dynamicCpuFluidSystem));
        auto fsView = ::Opm::gpuistl::make_view(*impl_->fluidSystemBuffer);

        impl_->managedFluidSystemView
            = make_gpu_managed_unique_ptr<typename ImplT::FluidSystemView>(fsView);

        // Build a default IntensiveQuantities prototype for the GPU side.
        BlackOilIntensiveQuantities<typename ImplT::CpuTag> cpuPrototype;
        impl_->prototype = cpuPrototype.template withOtherFluidSystem<typename ImplT::GpuTag>(
            *impl_->managedFluidSystemView);

        impl_->initialized = true;
        OpmLog::info(std::format(
            "[GpuBlackoilIntensiveQuantitiesDispatcher] initialized for {} cells ({})",
            cpuProblem.model().numGridDof(), Traits::name));
    }
    else {
        // The mixing limits (DRSDT, DRVDT, VAPPARS) change between calls:
        // at the start of a time step, on a chopped time step and on restart.
        // Re-reading them on every call cannot miss any of these.
        impl_->problemBuf->refreshMixingLimits(cpuProblem, 0u);
    }

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
    std::vector<GpuIntensiveQuantities> hostIQ(numDof, *impl_->prototype);

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

    dispatcher_detail::dispatcherUpdateAllCellsKernel<
        dispatcher_detail::updateEnergyQuantities<typename Traits::GpuTag>>
        <<<gridSize, blockSize>>>(
            impl_->problemView,
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

} // namespace Opm::gpuistl

#endif // OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_IMPL_HPP
