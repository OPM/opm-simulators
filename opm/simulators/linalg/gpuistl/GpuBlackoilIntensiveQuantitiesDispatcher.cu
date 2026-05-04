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
#include <opm/models/discretization/common/fvbaseelementcontextgpu.hh>
#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/GpuEclMaterialLawManager.hpp>
#include <opm/simulators/flow/GpuEclThermalLawManager.hpp>
#include <opm/simulators/flow/GpuFlowProblem.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher.hpp>

#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowGasWaterEnergyTypeTags.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <chrono>
#include <format>
#include <optional>
#include <vector>

namespace {

// Small helper to take a std::chrono::steady_clock timestamp.
inline auto dispatcherNow()
{
    return std::chrono::steady_clock::now();
}

inline double dispatcherElapsedMs(const std::chrono::steady_clock::time_point& start,
                                  const std::chrono::steady_clock::time_point& end)
{
    return std::chrono::duration<double, std::micro>(end - start).count()/ 1000.0;
}

} // namespace

namespace Opm::gpuistl {

namespace {

using DispatcherCpuTag = Opm::Properties::TTag::FlowGasWaterEnergyProblem;
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
    IntensiveQuantitiesT iq = outIntensiveQuantities[i];
    iq.updateSaturations(primaryVariables[i], 0, Opm::LinearizationType{});
    iq.update(problem, primaryVariables[i], static_cast<unsigned>(i), 0);
    iq.updateEnergyQuantities_(problem, static_cast<unsigned>(i), 0u);
    outIntensiveQuantities[i] = iq;
}

} // namespace

// =============================================================================
// Impl
// =============================================================================
template <class CpuTypeTag>
struct GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Impl {
    bool initialized = false;
    std::unique_ptr<DispatcherGpuFlowProblemBuf> problemBuf;
    DispatcherGpuFlowProblemView problemView{};
    DispatcherFluidSystemView* managedFluidSystemView = nullptr;
    std::optional<DispatcherGpuIntensiveQuantities> prototype;
    std::size_t callCount = 0;

    // Cumulative per-stage timings (milliseconds, host-side).
    double totalConvertMs = 0.0;
    double totalH2DMs = 0.0;
    double totalKernelMs = 0.0;
    double totalD2HMs = 0.0;
    double totalOverlayMs = 0.0;
    double totalDofsProcessed = 0.0;

    ~Impl()
    {
        if (managedFluidSystemView != nullptr) {
            managedFluidSystemView->~DispatcherFluidSystemView();
            (void)cudaFree(managedFluidSystemView);
            managedFluidSystemView = nullptr;
        }
        if (callCount > 0u) {
            const double totalMs =
                totalConvertMs + totalH2DMs + totalKernelMs + totalD2HMs + totalOverlayMs;
            const double avgConvertMs  = totalConvertMs  / static_cast<double>(callCount);
            const double avgH2DMs      = totalH2DMs      / static_cast<double>(callCount);
            const double avgKernelMs   = totalKernelMs   / static_cast<double>(callCount);
            const double avgD2HMs      = totalD2HMs      / static_cast<double>(callCount);
            const double avgOverlayMs  = totalOverlayMs  / static_cast<double>(callCount);
            const double avgTotalMs    = totalMs          / static_cast<double>(callCount);
            const double avgDofs       = totalDofsProcessed / static_cast<double>(callCount);
            Opm::OpmLog::info(std::format(
                "[GpuBlackoilIntensiveQuantitiesDispatcher] aggregate over {} call(s),"
                " {} dof-updates:\n"
                "  totals  : convert={:.3f}ms h2d={:.3f}ms kernel={:.3f}ms"
                " d2h={:.3f}ms overlay={:.3f}ms total={:.3f}ms\n"
                "  per-call: convert={:.3f}ms h2d={:.3f}ms kernel={:.3f}ms"
                " d2h={:.3f}ms overlay={:.3f}ms total={:.3f}ms avg_dofs={:.1f}",
                callCount,
                static_cast<std::size_t>(totalDofsProcessed),
                totalConvertMs, totalH2DMs, totalKernelMs, totalD2HMs, totalOverlayMs, totalMs,
                avgConvertMs,   avgH2DMs,   avgKernelMs,   avgD2HMs,   avgOverlayMs,   avgTotalMs,
                avgDofs));
        }
    }
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

    if (!impl_->initialized) {
        // Build the GPU FlowProblem from the CPU FlowProblem (one-time setup).
        impl_->problemBuf = std::make_unique<DispatcherGpuFlowProblemBuf>(cpuProblem);
        impl_->problemView = Opm::gpuistl::make_view(*impl_->problemBuf);

        // Place the FluidSystemView in unified memory so its device pointer
        // dereferences are valid both on host and device (mirrors the test).
        auto& dynamicCpuFluidSystem = DispatcherCpuFluidSystem::getNonStaticInstance();
        // The fluid-system buffer is owning device memory; leak it intentionally
        // for the lifetime of the process via a function-local static.
        using FsBufferType = decltype(Opm::gpuistl::copy_to_gpu(dynamicCpuFluidSystem));
        static FsBufferType s_fsBuffer = Opm::gpuistl::copy_to_gpu(dynamicCpuFluidSystem);
        auto fsView = Opm::gpuistl::make_view(s_fsBuffer);

        DispatcherFluidSystemView* managed = nullptr;
        OPM_GPU_SAFE_CALL(cudaMallocManaged(&managed, sizeof(DispatcherFluidSystemView)));
        new (managed) DispatcherFluidSystemView(fsView);
        impl_->managedFluidSystemView = managed;

        // Build a default IntensiveQuantities prototype for the GPU side.
        Opm::BlackOilIntensiveQuantities<DispatcherCpuTag> cpuPrototype;
        impl_->prototype = cpuPrototype.template withOtherFluidSystem<DispatcherGpuTag>(
            impl_->managedFluidSystemView);

        impl_->initialized = true;
        Opm::OpmLog::info(std::format(
            "[GpuBlackoilIntensiveQuantitiesDispatcher] initialized for {} cells",
            cpuProblem.model().numGridDof()));
    }

    // -------------------------------------------------------------------
    // 1. Convert CPU primary variables to the dispatcher's GPU primary
    //    variables (host-side).
    // -------------------------------------------------------------------
    const auto t0 = dispatcherNow();
    std::vector<DispatcherGpuPrimaryVariables> hostPriVars;
    hostPriVars.reserve(numDof);
    for (std::size_t i = 0; i < numDof; ++i) {
        // The BlackOilPrimaryVariables converting copy ctor handles the
        // TypeTag/storage difference.
        hostPriVars.emplace_back(*cpuPriVars[i]);
    }
    std::vector<DispatcherGpuIntensiveQuantities> hostIQ(numDof, *impl_->prototype);
    const auto t1 = dispatcherNow();

    // -------------------------------------------------------------------
    // 2. Upload host buffers to the device (host-to-device copies are
    //    embedded in the GpuBuffer ctors).
    // -------------------------------------------------------------------
    Opm::gpuistl::GpuBuffer<DispatcherGpuPrimaryVariables> primaryVariablesBuffer(hostPriVars);
    Opm::gpuistl::GpuBuffer<DispatcherGpuIntensiveQuantities> intensiveQuantitiesBuffer(hostIQ);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    const auto t2 = dispatcherNow();

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
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    const auto t3 = dispatcherNow();

    // -------------------------------------------------------------------
    // 4. Read the GPU IntensiveQuantities back to host memory.
    // -------------------------------------------------------------------
    OPM_GPU_SAFE_CALL(cudaMemcpy(hostIQ.data(),
                                 intensiveQuantitiesBuffer.data(),
                                 numDof * sizeof(DispatcherGpuIntensiveQuantities),
                                 cudaMemcpyDeviceToHost));
    const auto t4 = dispatcherNow();

    // -------------------------------------------------------------------
    // 5. Field-by-field overlay onto the caller's CPU IntensiveQuantities.
    //    The caller is expected to have run the CPU update first to
    //    populate any fields the dispatcher does not overwrite (e.g.
    //    \c mobility_, which the GPU relperm path currently leaves at 0).
    // -------------------------------------------------------------------
    for (std::size_t i = 0; i < numDof; ++i) {
        outIQ[i]->overlayBlackOilFieldsFrom(hostIQ[i]);
    }
    const auto t5 = dispatcherNow();

    // -------------------------------------------------------------------
    // Per-call timings.
    // -------------------------------------------------------------------
    const double convertMs = dispatcherElapsedMs(t0, t1);
    const double h2dMs     = dispatcherElapsedMs(t1, t2);
    const double kernelMs  = dispatcherElapsedMs(t2, t3);
    const double d2hMs     = dispatcherElapsedMs(t3, t4);
    const double overlayMs = dispatcherElapsedMs(t4, t5);

    impl_->totalConvertMs    += convertMs;
    impl_->totalH2DMs        += h2dMs;
    impl_->totalKernelMs     += kernelMs;
    impl_->totalD2HMs        += d2hMs;
    impl_->totalOverlayMs    += overlayMs;
    impl_->totalDofsProcessed += static_cast<double>(numDof);

    Opm::OpmLog::info(std::format(
        "[GpuBlackoilIntensiveQuantitiesDispatcher] call={} dofs={}"
        " convert={:.3f}ms h2d={:.3f}ms kernel={:.3f}ms d2h={:.3f}ms overlay={:.3f}ms total={:.3f}ms",
        impl_->callCount, numDof,
        convertMs, h2dMs, kernelMs, d2hMs, overlayMs,
        convertMs + h2dMs + kernelMs + d2hMs + overlayMs));

    ++impl_->callCount;
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
