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
 * The GPU side of \c GpuBlackoilNohystIntensiveQuantities, shared by the
 * per-configuration translation units
 * (\c GpuBlackoilNohystIntensiveQuantities<Configuration>.cu).
 *
 * Only to be included from \c .cu files. Each of them specializes
 * \c GpuDispatcherTraits for its (CPU TypeTag, PVT approaches) pairs and
 * explicitly instantiates \c makeImpl for them. Keeping one configuration per
 * translation unit keeps the nvcc compile time and memory of each unit to a
 * single \c FlowProblemBlackoil instantiation.
 */
#ifndef OPM_GPU_BLACKOIL_NOHYST_INTENSIVE_QUANTITIES_IMPL_HPP
#define OPM_GPU_BLACKOIL_NOHYST_INTENSIVE_QUANTITIES_IMPL_HPP

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
#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>
#include <opm/material/thermal/GpuEclThermalLawManager.hpp>

#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/discretization/common/tpfalinearizergpukernels.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBlackoilNohystIntensiveQuantities.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
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

/// The \c EclMaterialLaw::GpuManager of the material law of \p GpuTag. The
/// configurations call its static, host-only \c validate() on the CPU
/// manager. The converting constructor of the GPU manager runs the same check,
/// but only once the GPU problem is being built; calling it from \c validate
/// of the traits rejects the deck before anything is allocated on the GPU.
template <class GpuTag>
using GpuManagerFor = typename GetProp<GpuTag, Properties::MaterialLaw>::EclMaterialLawManager;

/// Validation shared by the three-phase black-oil configurations of
/// \c FlowProblemNohyst. The PVT approaches are checked by the facade, and
/// VAPOIL by \c BlackoilNohystDispatcherTraits.
template <class GpuTag, class FluidSystem, class Problem>
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
    GpuManagerFor<GpuTag>::validate(materialLawManager);

    // DRSDTCON needs the convective-mixing quantities, which the GPU
    // intensive quantities do not compute.
    const auto& schedule = problem.simulator().vanguard().schedule();
    for (std::size_t reportStep = 0; reportStep < schedule.size(); ++reportStep) {
        const auto& oilVaporization = schedule[reportStep].oilvap();
        if (oilVaporization.drsdtConvective()) {
            OPM_THROW(std::logic_error,
                      "DRSDTCON (convective DRSDT) is not supported on the GPU");
        }
    }
}

template <class GpuTag>
constexpr void assertNoDiffusionOrDispersion()
{
    // The GPU intensive quantities intentionally do not implement these modules.
    static_assert(!getPropValue<GpuTag, Properties::EnableDiffusion>());
    static_assert(!getPropValue<GpuTag, Properties::EnableDispersion>());
}

/// The \c GpuDispatcherTraits of the three-phase black-oil configurations of
/// \c FlowProblemNohyst, except \c name. They differ only in the GPU TypeTag,
/// whose fluid system fixes the PVT approaches: wet gas requires VAPOIL, dry gas
/// rejects it.
///
/// \c CpuTagT is a template parameter so that its properties are looked up only
/// on instantiation: this header is also included where \c FlowProblemNohyst is
/// only declared.
template <class GpuTagT, class CpuTagT = Properties::TTag::FlowProblemNohyst>
struct BlackoilNohystDispatcherTraits
{
    using CpuTag = CpuTagT;
    using GpuTag = GpuTagT;

    using FluidSystemView = GetPropType<GpuTag, Properties::FluidSystem>;
    using GpuPrimaryVariables = GetPropType<GpuTag, Properties::PrimaryVariables>;
    using GpuIntensiveQuantities = BlackOilIntensiveQuantities<GpuTag>;
    using GpuFlowProblemBuf = GpuFlowProblemBufferFor<GpuTag>;

    static constexpr bool wetGas
        = FluidSystemView::PvtApproachesType::gas == GasPvtApproach::WetGas;

    template <class Problem>
    static void validate(const Problem& problem)
    {
        using FluidSystem = GetPropType<CpuTag, Properties::FluidSystem>;
        assertNoDiffusionOrDispersion<GpuTag>();
        validateBlackoilNohyst<GpuTag, FluidSystem>(problem);

        // PVDG with VAPOIL is a valid CPU combination, but no deck uses it.
        // PVTG without VAPOIL is probably handled (Rv = 0 in every branch), but
        // no deck or test has verified it. Lift either check together with such
        // a test.
        if (FluidSystem::enableVaporizedOil() != wetGas) {
            const std::string message = wetGas
                ? "Wet gas (PVTG) without VAPOIL is not supported on the GPU"
                : "VAPOIL (vaporized oil) with dry gas (PVDG) is not supported on the GPU; "
                  "VAPOIL is supported with wet gas (PVTG)";
            OPM_THROW(std::logic_error, message);
        }
    }
};

// The GPU kernel that runs the IQ update for each cell
template <class GpuProblem, class PrimaryVariablesT, class IntensiveQuantitiesT>
__global__ void
dispatcherUpdateAllCellsKernel(GpuProblem problem,
                               IntensiveQuantitiesT prototype,
                               GpuView<const PrimaryVariablesT> primaryVariables,
                               GpuView<IntensiveQuantitiesT> outIntensiveQuantities,
                               std::size_t numCells)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numCells) {
        return;
    }
    prototype.update(problem, primaryVariables[i], static_cast<unsigned>(i), 0);
    outIntensiveQuantities[i] = prototype;
}

// =============================================================================
// ImplFor: the GPU state and the kernel launch of one configuration
// =============================================================================
template <class CpuTypeTag, class PvtApproaches>
struct ImplFor final : GpuBlackoilNohystIntensiveQuantities<CpuTypeTag>::Impl
{
    using Dispatcher = GpuBlackoilNohystIntensiveQuantities<CpuTypeTag>;
    using Problem = typename Dispatcher::Problem;
    using SolutionVector = typename Dispatcher::SolutionVector;
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
    using Scalar = GetPropType<CpuTypeTag, Properties::Scalar>;
    std::vector<Scalar> lastGasDissolution;
    std::vector<Scalar> lastOilVaporization;
    std::vector<Scalar> lastMaxOilSaturation;
    std::pair<Scalar, Scalar> vaporizationParameters{};
    std::vector<GpuPrimaryVariables> hostPriVars;
    std::vector<GpuIntensiveQuantities> hostIQ;
    GpuBuffer<GpuPrimaryVariables> devicePriVars;
    GpuBuffer<GpuIntensiveQuantities> deviceIQ;
    GpuBuffer<GpuIntensiveQuantities> previousDeviceIQ;
    GpuBuffer<Scalar> deviceVolumes;
    bool previousIqValid = false;


    explicit ImplFor(const Problem& cpuProblem)
    {
        // Build the GPU FlowProblem from the CPU FlowProblem (one-time setup).
        // The GPU material law manager is a snapshot of the CPU manager, so
        // every change to the scaled endpoints must come before this point:
        // applySwatinit() during equilibration and applyRestartSwatInit() on
        // restart. Both run in FlowProblemBlackoil::finishInit(), and the
        // first update() comes from Model::applyInitialSolution(), after it.
        // Hysteresis is the first state that will need a refresh after
        // construction (see GpuFlowProblem::refreshMixingLimits()).
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

        const auto numDof = cpuProblem.model().numGridDof();
        hostPriVars.resize(numDof);
        hostIQ.resize(numDof, *prototype);
        devicePriVars = GpuBuffer<GpuPrimaryVariables>(numDof);
        deviceIQ = GpuBuffer<GpuIntensiveQuantities>(numDof);
    }

    const char* name() const override
    {
        return Traits::name;
    }

    void invalidatePrevious() override
    {
        previousIqValid = false;
    }

    void advanceTimeLevel() override
    {
        if (previousDeviceIQ.size() != 0) {
            publishPrevious();
        }
    }

    void publishPrevious()
    {
        OPM_GPU_SAFE_CALL(cudaMemcpyAsync(previousDeviceIQ.data(), deviceIQ.data(),
                                         deviceIQ.size() * sizeof(GpuIntensiveQuantities),
                                         cudaMemcpyDeviceToDevice));
        previousIqValid = true;
    }

    void assemble(TpfaLinearizerGpuParams<CpuTypeTag>& buffers,
                  const Problem& cpuProblem, Scalar dt) override
    {
        const auto& model = cpuProblem.model();
        const auto numDof = deviceIQ.size();
        if (deviceVolumes.size() == 0) {
            std::vector<Scalar> volumes(numDof);
            for (std::size_t cell = 0; cell < numDof; ++cell) {
                volumes[cell] = model.dofTotalVolume(cell);
            }
            deviceVolumes = GpuBuffer<Scalar>(volumes);
            previousDeviceIQ = GpuBuffer<GpuIntensiveQuantities>(numDof);
        }

        if (model.enableStorageCache() && cpuProblem.recycleFirstIterationStorage()
            && cpuProblem.iterationContext().isFirstGlobalIteration()) {
            publishPrevious();
        } else if (!previousIqValid) {
#ifdef _OPENMP
#pragma omp parallel for if(numDof >= 16384)
#endif
            for (std::size_t cell = 0; cell < numDof; ++cell) {
                hostIQ[cell] = model.intensiveQuantityCache()[1][cell]
                    .template withOtherFluidSystem<GpuTag>(*managedFluidSystemView);
            }
            previousDeviceIQ.copyFromHost(hostIQ);
            previousIqValid = true;
        }

        const SimpleFIBlackOilModel<GpuTag, GpuView> modelView(
            make_view(deviceIQ), make_view(previousDeviceIQ), make_view(deviceVolumes));
        ThermalGasWaterFlowProblem<Scalar, GpuView> fluxProblem{};
        fluxProblem.moduleParams().saturationWeightedDensity =
            cpuProblem.moduleParams().saturationWeightedDensity;
        constexpr unsigned blockSize = 256;
        kernel_linearize<TpfaLinearizer<CpuTypeTag>, decltype(modelView),
                         BlackOilLocalResidualTPFA<GpuTag>>
            <<<(numDof + blockSize - 1) / blockSize, blockSize>>>(
                numDof, buffers.domainView(), buffers.neighborInfoView(),
                buffers.diagMatAddressView(), buffers.residualView(), modelView, dt, fluxProblem);
        OPM_GPU_SAFE_CALL(cudaGetLastError());
    }

    template <class Scalar, class ReadLimit>
    static void refreshLimit(GpuBuffer<Scalar>& device, std::vector<Scalar>& previous,
                             std::size_t count, ReadLimit readLimit)
    {
        bool changed = previous.size() != count;
        previous.resize(count);
#ifdef _OPENMP
#pragma omp parallel for reduction(||:changed) if(count >= 16384)
#endif
        for (std::size_t cell = 0; cell < count; ++cell) {
            const Scalar value = readLimit(static_cast<unsigned>(cell));
            changed |= previous[cell] != value;
            previous[cell] = value;
        }
        if (changed) {
            device.copyFromHost(previous);
        }
    }

    void update(const Problem& cpuProblem,
                const SolutionVector& solution,
                std::span<IntensiveQuantities> outIQ,
                unsigned timeIdx) override
    {
        const auto numDof = solution.size();
        const auto episode = std::max(0, cpuProblem.episodeIndex());
        const auto& oilVap = cpuProblem.simulator().vanguard().schedule()[episode].oilvap();
        const bool useVappars = oilVap.getType() == OilVaporizationProperties::OilVaporization::VAPPARS;
        const std::pair<Scalar, Scalar> parameters =
            useVappars ? std::pair<Scalar, Scalar>{oilVap.vap1(), oilVap.vap2()}
                       : std::pair<Scalar, Scalar>{};
        if (parameters != vaporizationParameters) {
            // The fluid-system view owns these scalars. Its PVT table views
            // still point at the original device allocations.
            managedFluidSystemView->setVapPars(parameters.first, parameters.second);
            vaporizationParameters = parameters;
        }
        if (!mixingLimitsCurrent) {
            if (useVappars || !lastMaxOilSaturation.empty()) {
                refreshLimit(problemBuf->maxOilSaturationStorage(), lastMaxOilSaturation,
                             numDof, [&](unsigned cell) { return cpuProblem.maxOilSaturation(cell); });
            }
            if (CpuFluidSystem::enableDissolvedGas()) {
                refreshLimit(problemBuf->maxGasDissolutionFactorStorage(), lastGasDissolution,
                             numDof, [&](unsigned cell) {
                                 return cpuProblem.maxGasDissolutionFactor(timeIdx, cell);
                             });
            }
            if (CpuFluidSystem::enableVaporizedOil()) {
                refreshLimit(problemBuf->maxOilVaporizationFactorStorage(), lastOilVaporization,
                             numDof, [&](unsigned cell) {
                                 return cpuProblem.maxOilVaporizationFactor(timeIdx, cell);
                             });
            }
        }
        mixingLimitsCurrent = false;
#ifdef _OPENMP
#pragma omp parallel for if(numDof >= 16384)
#endif
        for (std::size_t cell = 0; cell < numDof; ++cell) {
            hostPriVars[cell] = GpuPrimaryVariables(solution[cell]);
        }
        devicePriVars.copyFromHost(hostPriVars);

        constexpr unsigned blockSize = 64;
        const auto gridSize = (numDof + blockSize - 1) / blockSize;
        dispatcherUpdateAllCellsKernel<<<gridSize, blockSize>>>(
            problemView, *prototype,
            GpuView<const GpuPrimaryVariables>(devicePriVars.data(), numDof),
            make_view(deviceIQ), numDof);
        OPM_GPU_SAFE_CALL(cudaGetLastError());
        deviceIQ.copyToHost(hostIQ.data(), numDof);

#ifdef _OPENMP
#pragma omp parallel for if(numDof >= 16384)
#endif
        for (std::size_t cell = 0; cell < numDof; ++cell) {
            outIQ[cell].overlayBlackOilFieldsFrom(hostIQ[cell]);
        }
    }
};

} // namespace dispatcher_detail

template <class CpuTypeTag>
template <class PvtApproaches>
std::unique_ptr<typename GpuBlackoilNohystIntensiveQuantities<CpuTypeTag>::Impl>
GpuBlackoilNohystIntensiveQuantities<CpuTypeTag>::makeImpl(const Problem& cpuProblem)
{
    // The first call comes from the initial solution, so an unsupported
    // deck is rejected before the first time step.
    GpuDispatcherTraits<CpuTypeTag, PvtApproaches>::validate(cpuProblem);
    return std::make_unique<dispatcher_detail::ImplFor<CpuTypeTag, PvtApproaches>>(cpuProblem);
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_BLACKOIL_NOHYST_INTENSIVE_QUANTITIES_IMPL_HPP
