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
#include <opm/input/eclipse/Schedule/Events.hpp>
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
#include <opm/simulators/linalg/gpuistl/GpuFlowGasWaterEnergyBridge.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowGasWaterEnergyContract.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <format>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::gpuistl {

namespace {

using DispatcherGpuTag =
    Opm::Properties::TTag::FlowGasWaterEnergyDeviceTypeTag<Opm::gpuistl::GpuView>;

template <class ProblemT, class SolutionVectorT>
void validateDispatcherInputs(const ProblemT& problem, const SolutionVectorT& solution)
{
    if (solution.size() != static_cast<std::size_t>(problem.model().numGridDof())) {
        OPM_THROW(std::invalid_argument,
                  "GPU intensive-quantities dispatcher requires one entry per grid DoF");
    }
}

template <class ProblemT>
void validateGpuPropertyInputs(const ProblemT& problem)
{
    // The GPU view intentionally does not implement these modules
    static_assert(!Opm::getPropValue<DispatcherGpuTag, Opm::Properties::EnableDiffusion>());
    static_assert(!Opm::getPropValue<DispatcherGpuTag, Opm::Properties::EnableDispersion>());

    const auto& schedule = problem.simulator().vanguard().schedule();
    for (std::size_t reportStep = 0; reportStep < schedule.size(); ++reportStep) {
        if (schedule[reportStep].events().hasEvent(Opm::ScheduleEvents::GEO_MODIFIER)) {
            OPM_THROW(std::logic_error,
                      "GPU intensive-quantities evaluation does not support GEO_MODIFIER");
        }
    }

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

// The GPU kernel that runs the IQ update for each cell
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
    iq.update(problem, primaryVariables[i], static_cast<unsigned>(i), 0);
    iq.updateEnergyQuantities_(problem, static_cast<unsigned>(i), 0u);
}

} // namespace

// =============================================================================
// Impl
// =============================================================================
template <class CpuTypeTag>
struct GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Impl {
    std::unique_ptr<Bridge> bridge;
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
    const SolutionVector& solution,
    unsigned timeIdx)
{
    if (solution.size() == 0u) {
        return;
    }

    validateDispatcherInputs(cpuProblem, solution);
    validateGpuPropertyInputs(cpuProblem);

    if (!impl_->bridge) {
        impl_->bridge = std::make_unique<Bridge>();
        Opm::OpmLog::info(std::format(
            "[GpuBlackoilIntensiveQuantitiesDispatcher] initialized for {} cells",
            cpuProblem.model().numGridDof()));
    }

    impl_->bridge->updatePrimaryVariables(cpuProblem, solution, timeIdx);

    const unsigned blockSize = 64u;
    const unsigned gridSize =
        static_cast<unsigned>((solution.size() + blockSize - 1u) / blockSize);

    dispatcherUpdateAllCellsKernel<<<gridSize, blockSize, 0, impl_->bridge->stream()>>>(
        impl_->bridge->flowProblemView(),
        impl_->bridge->primaryVariablesView(timeIdx),
        impl_->bridge->intensiveQuantitiesView(timeIdx),
        solution.size());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    impl_->bridge->recordPropertyReady(timeIdx);
}

template <class CpuTypeTag>
void GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::materializeHostIntensiveQuantities(
    unsigned timeIdx,
    IntensiveQuantities* const* destination,
    std::size_t numDof)
{
    if (!impl_->bridge) {
        OPM_THROW(std::logic_error, "GPU intensive-quantities bridge has not been initialized");
    }
    impl_->bridge->materializeHostIntensiveQuantities(timeIdx, destination, numDof);
}

template <class CpuTypeTag>
bool GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::hasDeviceModelView() const
{
    return impl_->bridge && impl_->bridge->hasModelView();
}

template <class CpuTypeTag>
const typename GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Bridge&
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::bridge() const
{
    if (!impl_->bridge) {
        OPM_THROW(std::logic_error, "GPU intensive-quantities bridge has not been initialized");
    }
    return *impl_->bridge;
}

template <class CpuTypeTag>
typename GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::Bridge&
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::bridge()
{
    if (!impl_->bridge) {
        OPM_THROW(std::logic_error, "GPU intensive-quantities bridge has not been initialized");
    }
    return *impl_->bridge;
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
