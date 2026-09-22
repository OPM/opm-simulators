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
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
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
#include <opm/simulators/linalg/gpuistl/GpuBlackoilNewtonUpdate.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBlackoilNewtonValidation.hpp>
#include <opm/simulators/flow/NonlinearSolver.hpp>

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

    if (!problem.usesDefaultRockCompaction()) {
        OPM_THROW(std::logic_error, "GPU properties do not support dynamic rock compaction");
    }
    const auto& schedule = problem.simulator().vanguard().schedule();
    for (std::size_t reportStep = 0; reportStep < schedule.size(); ++reportStep) {
        if (!schedule[reportStep].actions().empty()) {
            OPM_THROW(std::logic_error, "GPU properties do not support runtime schedule actions");
        }
        if (schedule[reportStep].oilvap().defined()) {
            OPM_THROW(std::logic_error,
                      "GPU properties do not support time-dependent VAPPARS/DRSDT/DRVDT/DRSDTCON inputs");
        }
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
    bool validatedBranches{false};
    std::array<bool, 3> reportedRoundoff{};
};

template <class CpuTypeTag>
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::GpuBlackoilIntensiveQuantitiesDispatcher()
    : impl_(std::make_unique<Impl>())
{
}

template <class CpuTypeTag>
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::~GpuBlackoilIntensiveQuantitiesDispatcher() = default;

template <class CpuTypeTag>
void GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::reportTransferCounters() const
{
    if (impl_->bridge) {
        const auto& c = impl_->bridge->transferCounters();
        OpmLog::info(std::format(
            "[GPU Newton transfers] updates={} pv_uploads={} pv_upload_bytes={} pv_downloads={} pv_download_bytes={} correction_downloads={} correction_download_bytes={} iq_downloads={} iq_download_bytes={} bridge_allocations={} static_upload_batches={} bridge_allocation_bytes={} static_upload_calls={} static_upload_bytes={} solver_correction_allocations={} correction_history_uploads={} correction_history_upload_bytes={}",
            c.successfulNewtonUpdates, c.primaryVariableUploads, c.primaryVariableUploadBytes,
            c.primaryVariableDownloads, c.primaryVariableDownloadBytes, c.correctionDownloads,
            c.correctionDownloadBytes, c.intensiveQuantityDownloads, c.intensiveQuantityDownloadBytes,
            c.ownedBufferAllocations, c.staticUploadBatches, c.ownedBufferAllocationBytes,
            c.staticUploadCalls, c.staticUploadBytes, c.solverCorrectionAllocations,
            c.correctionHistoryUploads, c.correctionHistoryUploadBytes));
    }
}

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
    if (!impl_->bridge || !impl_->bridge->initializedFor(solution.size())) {
        validateGpuPropertyInputs(cpuProblem);
    }

    if (!impl_->bridge) {
        impl_->bridge = std::make_unique<Bridge>();
        Opm::OpmLog::info(std::format(
            "[GpuBlackoilIntensiveQuantitiesDispatcher] initialized for {} cells",
            cpuProblem.model().numGridDof()));
    }

    impl_->bridge->updatePrimaryVariables(cpuProblem, solution, timeIdx);

    evaluateResident(timeIdx);
}

template <class CpuTypeTag>
bool GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::hasBridge() const
{
    return static_cast<bool>(impl_->bridge);
}

template <class CpuTypeTag>
void GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::evaluateResident(unsigned timeIdx)
{
    const auto numCells = bridge().primaryVariablesView(timeIdx).size();
    const unsigned blockSize = 64u;
    const unsigned gridSize =
        static_cast<unsigned>((numCells + blockSize - 1u) / blockSize);

    dispatcherUpdateAllCellsKernel<<<gridSize, blockSize, 0, impl_->bridge->stream()>>>(
        impl_->bridge->flowProblemView(),
        impl_->bridge->primaryVariablesView(timeIdx),
        impl_->bridge->intensiveQuantitiesView(timeIdx),
        numCells);
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    impl_->bridge->recordPropertyReady(timeIdx);
}

template <class CpuTypeTag>
unsigned GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::applyNewtonUpdate(
    const Problem& problem, const BlackoilNewtonParams<Scalar>& params,
    Scalar relaxation, bool useSOR, bool stabilize, bool validate)
{
    using Indices = GetPropType<CpuTypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<CpuTypeTag, Properties::FluidSystem>;
    using Correction = GetPropType<CpuTypeTag, Properties::GlobalEqVector>;
    auto& state = bridge();
    SolutionVector reference;
    std::vector<std::uint8_t> referenceSwitch;
    if (validate) {
        if (!impl_->validatedBranches) {
            validateBlackoilNewtonBranches<CpuTypeTag>(state, problem, problem.model().solution(0), params);
            impl_->validatedBranches = true;
        }
        reference = problem.model().solution(0);
        Correction correction(state.numDof()), previous(state.numDof());
        state.materializeHostCorrection(correction, "shadow validation");
        state.materializePreviousCorrection(previous, "shadow validation");
        referenceSwitch = state.materializeSwitchHistory();
        if (stabilize) {
            Opm::detail::stabilizeNonlinearUpdate(correction, previous, relaxation,
                useSOR ? NonlinearRelaxType::SOR : NonlinearRelaxType::Dampen);
        }
        for (unsigned cell = 0; cell < state.numDof(); ++cell) {
            const auto current = reference[cell];
            referenceSwitch[cell] = BlackOilNewtonUpdate<CpuTypeTag>::update(
                problem, FluidSystem{}, cell, reference[cell], current, correction[cell],
                params, referenceSwitch[cell] != 0);
        }
    }
    const auto before = state.transferCounters();
    launchBlackoilNewtonUpdate<typename Bridge::DeviceTypeTagPublic>(
        state, params, relaxation, useSOR, stabilize);
    const auto status = state.readUpdateStatus();
    if (status[0]) {
        OPM_THROW_PROBLEM(NumericalProblem, std::format(
            "GPU Newton update failed: cell {}, reason {}", status[0] - 1, status[1]));
    }
    if (validate) {
        SolutionVector candidate(state.numDof());
        state.materializeCandidatePrimaryVariables(candidate);
        const auto candidateSwitch = state.materializeCandidateSwitchHistory();
        const auto mismatch = [&](unsigned cell, const std::string& field) {
            OPM_THROW(std::runtime_error, std::format(
                "GPU Newton shadow mismatch: cell {}, time {} s, iteration {}, field {}",
                cell, problem.simulator().time(), state.transferCounters().successfulNewtonUpdates, field));
        };
        for (unsigned cell = 0; cell < state.numDof(); ++cell) {
            const auto& a = reference[cell];
            const auto& b = candidate[cell];
            if (a.primaryVarsMeaningWater() != b.primaryVarsMeaningWater()
                || a.primaryVarsMeaningGas() != b.primaryVarsMeaningGas()
                || a.primaryVarsMeaningPressure() != b.primaryVarsMeaningPressure()
                || a.primaryVarsMeaningBrine() != b.primaryVarsMeaningBrine()
                || a.primaryVarsMeaningSolvent() != b.primaryVarsMeaningSolvent()
                || a.pvtRegionIndex() != b.pvtRegionIndex()
                || a.capillaryPressureFactor() != b.capillaryPressureFactor()
                || a.pressureScale() != b.pressureScale()
                || referenceSwitch[cell] != candidateSwitch[cell]) {
                mismatch(cell, "metadata/switch history");
            }
            for (unsigned field = 0; field < Indices::numEq; ++field) {
                Scalar expected = a[field];
                Scalar actual = b[field];
                Scalar atol = 1e-12;
                if (field == Indices::pressureSwitchIdx) {
                    expected *= a.pressureScale();
                    actual *= b.pressureScale();
                    atol = 1e-4;
                } else if (field == Indices::temperatureIdx) {
                    atol = 1e-8;
                }
                if (actual != expected && !impl_->reportedRoundoff[field]) {
                    impl_->reportedRoundoff[field] = true;
                    OpmLog::info(std::format(
                        "[GPU Newton roundoff] time={} update={} cell={} field={} relaxation={} switched={} CPU={:.17g} GPU={:.17g} difference={:.17g}",
                        problem.simulator().time(), state.transferCounters().successfulNewtonUpdates,
                        cell, field, relaxation, referenceSwitch[cell], expected, actual, actual - expected));
                    if (field == Indices::waterSwitchIdx
                        && a.primaryVarsMeaningWater() == BlackOil::WaterMeaning::Rsw
                        && b.primaryVarsMeaningWater() == BlackOil::WaterMeaning::Rsw) {
                        diagnoseBlackoilNewtonRsw<CpuTypeTag>(state, a.pvtRegionIndex(),
                            a[Indices::temperatureIdx], a[Indices::pressureSwitchIdx] * a.pressureScale());
                    }
                }
                if (!std::isfinite(actual)
                    || std::abs(actual - expected) > atol + 1e-10 * std::abs(expected)) {
                    mismatch(cell, std::format("{} (CPU {}, GPU {})", field, expected, actual));
                }
            }
        }
    }
    state.commitDeviceUpdate();
    const auto& after = state.transferCounters();
    if (after.primaryVariableUploads != before.primaryVariableUploads
        || after.correctionDownloads != before.correctionDownloads
        || after.ownedBufferAllocations != before.ownedBufferAllocations) {
        OPM_THROW(std::logic_error, "Resident Newton update performed a forbidden transfer or allocation");
    }
    return status[2];
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
