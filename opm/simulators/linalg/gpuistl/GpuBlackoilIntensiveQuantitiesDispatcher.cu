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

#include <cuda_runtime.h>

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/material/common/ResetLocale.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>

#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/material/fluidmatrixinteractions/GpuEclMaterialLawManager.hpp>
#include <opm/material/thermal/GpuEclThermalLawManager.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher.hpp>

#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowGasWaterEnergyTypeTags.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <format>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::gpuistl {

namespace {

using DispatcherCpuTag = Opm::Properties::TTag::FlowGasWaterEnergyCpuKernelBase;
using DispatcherGpuTag = Opm::Properties::TTag::FlowGasWaterEnergyDummyProblemGPU;

using DispatcherScalar = Opm::GetPropType<DispatcherCpuTag, Opm::Properties::Scalar>;

using DispatcherCpuFluidSystem = Opm::BlackOilFluidSystem<DispatcherScalar>;
using DispatcherFluidSystemView
    = Opm::BlackOilFluidSystemNonStatic<DispatcherScalar,
                                        Opm::BlackOilDefaultFluidSystemIndices,
                                        Opm::gpuistl::GpuView>;

using DispatcherGpuPrimaryVariables
    = Opm::BlackOilPrimaryVariables<DispatcherGpuTag, Opm::gpuistl::MiniVector>;
using DispatcherGpuIntensiveQuantities = Opm::BlackOilIntensiveQuantities<DispatcherGpuTag>;

using DispatcherCpuMaterialLawManager =
    typename Opm::GetProp<DispatcherCpuTag, Opm::Properties::MaterialLaw>::EclMaterialLawManager;
using DispatcherTraits = typename DispatcherCpuMaterialLawManager::MaterialLaw::Traits;
using DispatcherTwoPhaseTraits =
    Opm::TwoPhaseMaterialTraits<DispatcherScalar,
                                DispatcherTraits::wettingPhaseIdx,
                                DispatcherTraits::nonWettingPhaseIdx>;
using DispatcherGpuPiecewiseLinearParamsBuf =
    Opm::PiecewiseLinearTwoPhaseMaterialParams<DispatcherTwoPhaseTraits,
                                               Opm::gpuistl::GpuView<const DispatcherScalar>>;
using DispatcherGpuPiecewiseLinearLawBuf =
    Opm::PiecewiseLinearTwoPhaseMaterial<DispatcherTwoPhaseTraits,
                                         DispatcherGpuPiecewiseLinearParamsBuf>;
using DispatcherGpuMaterialLawParamsBuf =
    Opm::EclTwoPhaseMaterialParams<DispatcherTraits,
                                   DispatcherGpuPiecewiseLinearParamsBuf,
                                   DispatcherGpuPiecewiseLinearParamsBuf,
                                   DispatcherGpuPiecewiseLinearParamsBuf,
                                   Opm::gpuistl::ValueAsPointer>;
using DispatcherGpuMaterialLawBuf =
    Opm::EclTwoPhaseMaterial<DispatcherTraits,
                             DispatcherGpuPiecewiseLinearLawBuf,
                             DispatcherGpuPiecewiseLinearLawBuf,
                             DispatcherGpuPiecewiseLinearLawBuf,
                             DispatcherGpuMaterialLawParamsBuf>;
using DispatcherGpuManagerBuf =
    Opm::EclMaterialLaw::GpuManager<DispatcherTraits,
                                    DispatcherGpuPiecewiseLinearLawBuf,
                                    DispatcherGpuPiecewiseLinearLawBuf,
                                    Opm::gpuistl::GpuBuffer,
                                    DispatcherGpuMaterialLawBuf>;
using DispatcherGpuThermalManagerBuf =
    Opm::EclThermalLaw::GpuManager<DispatcherScalar,
                                   DispatcherFluidSystemView,
                                   Opm::gpuistl::GpuBuffer,
                                   Opm::gpuistl::GpuView>;
using DispatcherGpuFlowProblemBuf =
    Opm::GpuFlowProblem<DispatcherScalar,
                        DispatcherGpuManagerBuf,
                        Opm::gpuistl::GpuBuffer,
                        DispatcherGpuThermalManagerBuf>;
using DispatcherGpuFlowProblemView =
    decltype(Opm::gpuistl::make_view(std::declval<DispatcherGpuFlowProblemBuf&>()));

// The GPU view intentionally does not implement these modules
static_assert(!Opm::getPropValue<DispatcherGpuTag, Opm::Properties::EnableDiffusion>());
static_assert(!Opm::getPropValue<DispatcherGpuTag, Opm::Properties::EnableDispersion>());

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

template <class ProblemT>
void validateGpuPropertyInputs(const ProblemT& problem)
{
    const auto materialLawManager = problem.materialLawManager();
    if (materialLawManager == nullptr) {
        OPM_THROW(std::invalid_argument,
                  "GPU property evaluation requires a material-law manager");
    }
    if (materialLawManager->twoPhaseApproach() != Opm::EclTwoPhaseApproach::GasWater) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation requires the gas-water two-phase approach");
    }
    if (materialLawManager->hasOil() || !materialLawManager->hasGas()
        || !materialLawManager->hasWater()) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation requires water and gas without an active oil phase");
    }
    if (!materialLawManager->satCurveIsAllPiecewiseLinear()) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation requires piecewise-linear saturation tables");
    }
    if (materialLawManager->enableHysteresis()
        || materialLawManager->enableEndPointScaling()
        || materialLawManager->hasDirectionalRelperms()
        || materialLawManager->hasDirectionalImbnum()) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation does not support hysteresis, endpoint scaling, "
                  "or directional saturation functions");
    }

    const auto thermalLawManager = problem.thermalLawManager();
    if (thermalLawManager == nullptr) {
        OPM_THROW(std::invalid_argument,
                  "GPU property evaluation requires a thermal-law manager");
    }
    if (thermalLawManager->solidEnergyApproach()
            != Opm::EclSolidEnergyApproach::Specrock
        || thermalLawManager->thermalConductionApproach()
            != Opm::EclThermalConductionApproach::Thconr) {
        OPM_THROW(std::logic_error,
                  "GPU property evaluation requires SPECROCK solid energy and THCONR "
                  "thermal conduction");
    }
}

template <class GpuProblem, class PrimaryVariablesT, class IntensiveQuantitiesT>
__global__ void
dispatcherUpdateAllCellsKernel(GpuProblem problem,
                               Opm::gpuistl::GpuView<const PrimaryVariablesT> primaryVariables,
                               Opm::gpuistl::GpuView<IntensiveQuantitiesT> outIntensiveQuantities,
                               std::size_t numCells)
{
    const std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numCells) {
        return;
    }
    IntensiveQuantitiesT& iq = outIntensiveQuantities[i];
    iq.updateSaturations(primaryVariables[i], 0, Opm::LinearizationType{});
    iq.update(problem, primaryVariables[i], static_cast<unsigned>(i), 0);
    iq.updateEnergyQuantities_(problem, static_cast<unsigned>(i), 0u);
}

} // namespace

// =============================================================================
// Impl
// =============================================================================
template <class CpuTypeTag>
struct GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Impl {
    using DynamicCpuFluidSystem
        = std::remove_reference_t<decltype(DispatcherCpuFluidSystem::getNonStaticInstance())>;
    using FluidSystemBuffer
        = decltype(Opm::gpuistl::copy_to_gpu(std::declval<DynamicCpuFluidSystem&>()));
    using ManagedFluidSystemView
        = std::unique_ptr<DispatcherFluidSystemView,
                          Opm::gpuistl::GpuManagedDeleter<DispatcherFluidSystemView>>;

    bool initialized = false;
    std::unique_ptr<DispatcherGpuFlowProblemBuf> problemBuf;
    DispatcherGpuFlowProblemView problemView{};
    std::unique_ptr<FluidSystemBuffer> fluidSystemBuffer;
    ManagedFluidSystemView managedFluidSystemView;
    std::optional<DispatcherGpuIntensiveQuantities> prototype;
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
    if (numDof == 0u) {
        return;
    }

    validateDispatcherInputs(cpuProblem, cpuPriVars, outIQ, numDof);
    validateGpuPropertyInputs(cpuProblem);

    if (!impl_->initialized) {
        // Build the GPU FlowProblem from the CPU FlowProblem (one-time setup).
        impl_->problemBuf = std::make_unique<DispatcherGpuFlowProblemBuf>(cpuProblem);
        impl_->problemView = Opm::gpuistl::make_view(*impl_->problemBuf);

        // Place the FluidSystemView in unified memory so its device pointer
        // dereferences are valid both on host and device (mirrors the test).
        auto& dynamicCpuFluidSystem = DispatcherCpuFluidSystem::getNonStaticInstance();
        // Keep the owning buffer with this dispatcher. A process-wide static
        // would retain the first deck's PVT tables across later simulations.
        impl_->fluidSystemBuffer
            = std::make_unique<typename Impl::FluidSystemBuffer>(
                Opm::gpuistl::copy_to_gpu(dynamicCpuFluidSystem));
        auto fsView = Opm::gpuistl::make_view(*impl_->fluidSystemBuffer);

        impl_->managedFluidSystemView
            = Opm::gpuistl::make_gpu_managed_unique_ptr<DispatcherFluidSystemView>(fsView);

        // Build a default IntensiveQuantities prototype for the GPU side.
        Opm::BlackOilIntensiveQuantities<DispatcherCpuTag> cpuPrototype;
        impl_->prototype = cpuPrototype.template withOtherFluidSystem<DispatcherGpuTag>(
            *impl_->managedFluidSystemView);

        impl_->initialized = true;
        Opm::OpmLog::info(std::format(
            "[GpuBlackoilIntensiveQuantitiesDispatcher] initialized for {} cells",
            cpuProblem.model().numGridDof()));
    }

    // -------------------------------------------------------------------
    // 1. Convert CPU primary variables to the dispatcher's GPU primary
    //    variables (host-side).
    // -------------------------------------------------------------------
    std::vector<DispatcherGpuPrimaryVariables> hostPriVars;
    hostPriVars.reserve(numDof);
    for (std::size_t i = 0; i < numDof; ++i) {
        // The BlackOilPrimaryVariables converting copy ctor handles the
        // TypeTag/storage difference.
        hostPriVars.emplace_back(*cpuPriVars[i]);
    }
    std::vector<DispatcherGpuIntensiveQuantities> hostIQ(numDof, *impl_->prototype);

    // -------------------------------------------------------------------
    // 2. Upload host buffers to the device (host-to-device copies are
    //    embedded in the GpuBuffer ctors).
    // -------------------------------------------------------------------
    Opm::gpuistl::GpuBuffer<DispatcherGpuPrimaryVariables> primaryVariablesBuffer(hostPriVars);
    Opm::gpuistl::GpuBuffer<DispatcherGpuIntensiveQuantities> intensiveQuantitiesBuffer(hostIQ);

    // -------------------------------------------------------------------
    // 3. Launch the per-cell update kernel.
    // -------------------------------------------------------------------
    const unsigned blockSize = 64u;
    const unsigned gridSize = static_cast<unsigned>((numDof + blockSize - 1u) / blockSize);

    dispatcherUpdateAllCellsKernel<<<gridSize, blockSize>>>(
        impl_->problemView,
        Opm::gpuistl::GpuView<const DispatcherGpuPrimaryVariables>(
            primaryVariablesBuffer.data(), primaryVariablesBuffer.size()),
        Opm::gpuistl::GpuView<DispatcherGpuIntensiveQuantities>(
            intensiveQuantitiesBuffer.data(), intensiveQuantitiesBuffer.size()),
        numDof);
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    // -------------------------------------------------------------------
    // 4. Read the GPU IntensiveQuantities back to host memory.
    // -------------------------------------------------------------------
    OPM_GPU_SAFE_CALL(cudaMemcpy(hostIQ.data(),
                                 intensiveQuantitiesBuffer.data(),
                                 numDof * sizeof(DispatcherGpuIntensiveQuantities),
                                 cudaMemcpyDeviceToHost));

    // -------------------------------------------------------------------
    // 5. Field-by-field overlay onto the caller's CPU IntensiveQuantities.
    //    The caller is expected to have run the CPU update first to
    //    populate any fields the dispatcher does not overwrite (e.g.
    //    \c mobility_, which the GPU relperm path currently leaves at 0).
    // -------------------------------------------------------------------
    for (std::size_t i = 0; i < numDof; ++i) {
        outIQ[i]->overlayBlackOilFieldsFrom(hostIQ[i]);
    }
}

// =============================================================================
// Explicit instantiation for the user-facing CO2STORE simulation TypeTag,
// namely \c FlowGasWaterEnergyProblem (used by the Flow simulator binary
// when a CO2STORE deck with WATER+GAS+THERMAL is loaded). The GPU-specific
// property bindings live in \c GpuFlowGasWaterEnergyTypeTags.hpp.
// =============================================================================
template class GpuBlackoilIntensiveQuantitiesDispatcher<
    Opm::Properties::TTag::FlowGasWaterEnergyProblem>;

// =============================================================================
// Explicit instantiation for the GPU-assembly simulation TypeTag
// \c FlowGasWaterEnergyProblemGPU (used by the flow_gpu binary that runs both
// matrix assembly and intensive-quantities computation on the GPU).
// This TypeTag inherits all physics from \c FlowGasWaterEnergyProblem and
// therefore has the same CPU-side Problem / PrimaryVariables /
// IntensiveQuantities types; the GPU kernel dispatch internally still uses
// the \c FlowGasWaterEnergyKernelBaseGPU / \c FlowGasWaterEnergyDummyProblemGPU
// chain.
// =============================================================================
template class GpuBlackoilIntensiveQuantitiesDispatcher<
    Opm::Properties::TTag::FlowGasWaterEnergyProblemGPU>;

} // namespace Opm::gpuistl
