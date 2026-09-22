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
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/flow/NonlinearSolver.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <format>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::gpuistl {


// status = {first failing cell + 1, failure code, switched count, reserved}.
// Failure codes: 1 = nonfinite correction; 2 = nonfinite resulting variable.
// Only candidate buffers are written: a failure never publishes partial PVs
// or switch history. Correction history stores the unrelaxed correction, as in
// detail::stabilizeNonlinearUpdate(), even when the relaxation factor is one.
template<class TypeTag, class Problem, class FluidSystem, class PrimaryVariables, class Scalar>
__global__ void blackoilNewtonUpdateKernel(
    Problem problem, const FluidSystem* fluidSystem,
    GpuView<const PrimaryVariables> current, GpuView<PrimaryVariables> candidate,
    Scalar* correction, Scalar* previousCorrection,
    GpuView<std::uint8_t> previousSwitch, GpuView<std::uint8_t> candidateSwitch,
    std::uint32_t* status, BlackoilNewtonParams<Scalar> params,
    Scalar relaxation, bool useSOR, bool stabilize)
{
    const std::size_t cell = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (cell >= current.size()) {
        return;
    }
    constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    MiniVector<Scalar, numEq> update;
    for (unsigned eq = 0; eq < numEq; ++eq) {
        const std::size_t index = cell * numEq + eq;
        const Scalar unrelaxed = correction[index];
        Scalar delta = unrelaxed;
        if (stabilize) {
            const Scalar old = previousCorrection[index];
            previousCorrection[index] = unrelaxed;
            if (relaxation != Scalar{1}) {
                delta *= relaxation;
                if (useSOR) {
                    delta += (Scalar{1} - relaxation) * old;
                }
            }
        }
        correction[index] = delta;
        update[eq] = delta;
        if (!std::isfinite(delta)) {
            if (atomicCAS(status, 0u, static_cast<unsigned>(cell + 1)) == 0u) {
                status[1] = 1u;
            }
            return;
        }
    }
    candidate[cell] = current[cell];
    const bool switched = BlackOilNewtonUpdate<TypeTag>::update(
        problem, *fluidSystem, static_cast<unsigned>(cell), candidate[cell],
        current[cell], update, params, previousSwitch[cell] != 0);
    for (unsigned eq = 0; eq < numEq; ++eq) {
        if (!std::isfinite(candidate[cell][eq])) {
            if (atomicCAS(status, 0u, static_cast<unsigned>(cell + 1)) == 0u) {
                status[1] = 2u;
            }
            return;
        }
    }
    candidateSwitch[cell] = switched;
    if (switched) {
        atomicAdd(status + 2, 1u);
    }
}

// The caller orders solver completion before this operation, then reads the
// compact status and commits the candidate only on success. There is no bulk
// transfer, allocation, or device-wide synchronization in this launch helper.
template<class DeviceTypeTag, class Bridge>
void launchBlackoilNewtonUpdate(Bridge& bridge,
                               const BlackoilNewtonParams<typename Bridge::Scalar>& params,
                               typename Bridge::Scalar relaxation,
                               bool useSOR, bool stabilize)
{
    bridge.clearUpdateStatus();
    const auto current = bridge.primaryVariablesView(0);
    if (current.size() == 0) {
        return;
    }
    constexpr unsigned blockSize = 128;
    const unsigned gridSize = static_cast<unsigned>((current.size() + blockSize - 1) / blockSize);
    blackoilNewtonUpdateKernel<DeviceTypeTag><<<gridSize, blockSize, 0, bridge.stream()>>>(
        bridge.flowProblemView(), &bridge.deviceFluidSystem(), current,
        bridge.scratchPrimaryVariablesView(), bridge.correction().data(),
        bridge.previousCorrection().data(), bridge.switchedLastIterationView(),
        bridge.scratchSwitchHistoryView(), bridge.updateStatusView().data(),
        params, relaxation, useSOR, stabilize);
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}



// Validation-only, bounded synthetic cases using the real deck's fluid and
// material laws. Every mutable device object belongs to this function. In
// particular neither resident primary variables nor correction/switch history
// is modified, and no commit operation is called.
template<class CpuTypeTag, class Bridge, class Problem, class Solution>
void validateBlackoilNewtonBranches(Bridge& bridge, const Problem& problem,
                                   const Solution& solution,
                                   const BlackoilNewtonParams<typename Bridge::Scalar>& params)
{
    using Scalar = typename Bridge::Scalar;
    using Indices = GetPropType<CpuTypeTag, Properties::Indices>;
    using CpuPV = GetPropType<CpuTypeTag, Properties::PrimaryVariables>;
    using DevicePV = typename Bridge::DevicePrimaryVariablesPublic;
    using FluidSystem = GetPropType<CpuTypeTag, Properties::FluidSystem>;
    using Correction = GetPropType<CpuTypeTag, Properties::GlobalEqVector>;
    using DeviceTag = typename Bridge::DeviceTypeTagPublic;
    constexpr unsigned numEq = getPropValue<CpuTypeTag, Properties::NumEq>();
    const unsigned count = static_cast<unsigned>(std::min<std::size_t>(8, solution.size()));
    if (count == 0) {
        return;
    }
    detail::GpuMemoryCounts transfers;
    const detail::ScopedGpuMemoryAccounting accounting(transfers);
    std::uint64_t downloadCalls = 0, downloadBytes = 0;
    unsigned batches = 0;
    for (unsigned scenario = 0; scenario < 4; ++scenario) {
        const bool waterVaporization = scenario == 0 || scenario == 2;
        if ((waterVaporization && !FluidSystem::enableVaporizedWater())
            || (!waterVaporization && !FluidSystem::enableDissolvedGasInWater())) {
            continue;
        }
        for (const bool useSOR : {false, true}) {
            for (const Scalar relaxation : {Scalar{0.65}, Scalar{1}}) {
                std::vector<CpuPV> current(count), reference(count);
                std::vector<DevicePV> deviceInput;
                deviceInput.reserve(count);
                std::vector<Scalar> raw(count * numEq), old(count * numEq);
                std::vector<std::uint8_t> history(count), expectedHistory(count);
                Correction correction(count), previous(count);
                for (unsigned cell = 0; cell < count; ++cell) {
                    auto& pv = current[cell];
                    pv = solution[cell];
                    const Scalar temperature = pv[Indices::temperatureIdx];
                    const Scalar pressure = pv[Indices::pressureSwitchIdx] * pv.pressureScale();
                    const auto region = pv.pvtRegionIndex();
                    pv.setPrimaryVarsMeaningGas(CpuPV::GasMeaning::Disabled);
                    pv.setPrimaryVarsMeaningPressure(scenario == 3
                        ? CpuPV::PressureMeaning::Pw : CpuPV::PressureMeaning::Pg);
                    Scalar compositionDelta = 0.002;
                    if (scenario < 2) {
                        pv.setPrimaryVarsMeaningWater(CpuPV::WaterMeaning::Sw);
                        pv[Indices::waterSwitchIdx] = scenario == 0 ? Scalar{-0.01} : Scalar{1.01};
                        if (scenario == 1) {
                            compositionDelta = -compositionDelta;
                        }
                    }
                    else if (scenario == 2) {
                        pv.setPrimaryVarsMeaningWater(CpuPV::WaterMeaning::Rvw);
                        const Scalar saturated = FluidSystem::gasPvt().saturatedWaterVaporizationFactor(
                            region, temperature, pressure, Scalar{0});
                        pv[Indices::waterSwitchIdx] = 2 * saturated + Scalar{1e-6};
                        compositionDelta = Scalar{0.001} * pv[Indices::waterSwitchIdx];
                    }
                    else {
                        pv.setPrimaryVarsMeaningWater(CpuPV::WaterMeaning::Rsw);
                        const Scalar saturated = std::min(
                            FluidSystem::waterPvt().saturatedGasDissolutionFactor(
                                region, temperature, pressure, Scalar{0}),
                            problem.maxGasDissolutionFactor(0, cell));
                        pv[Indices::waterSwitchIdx] = 2 * saturated + Scalar{1e-6};
                        compositionDelta = Scalar{0.001} * pv[Indices::waterSwitchIdx];
                    }
                    correction[cell] = Scalar{0};
                    correction[cell][Indices::pressureSwitchIdx] = Scalar{128} / pv.pressureScale();
                    correction[cell][Indices::temperatureIdx] = Scalar{0.125};
                    correction[cell][Indices::waterSwitchIdx] = compositionDelta;
                    for (unsigned eq = 0; eq < numEq; ++eq) {
                        raw[cell * numEq + eq] = correction[cell][eq];
                        old[cell * numEq + eq] = -Scalar{0.25} * correction[cell][eq];
                        previous[cell][eq] = old[cell * numEq + eq];
                    }
                    history[cell] = cell % 2;
                    deviceInput.emplace_back(pv);
                    reference[cell] = pv;
                }
                Opm::detail::stabilizeNonlinearUpdate(correction, previous, relaxation,
                    useSOR ? NonlinearRelaxType::SOR : NonlinearRelaxType::Dampen);
                unsigned switched = 0;
                const auto expectedMeaning = scenario == 0 ? CpuPV::WaterMeaning::Rvw
                    : scenario == 1 ? CpuPV::WaterMeaning::Rsw : CpuPV::WaterMeaning::Sw;
                for (unsigned cell = 0; cell < count; ++cell) {
                    expectedHistory[cell] = BlackOilNewtonUpdate<CpuTypeTag>::update(
                        problem, FluidSystem{}, cell, reference[cell], current[cell],
                        correction[cell], params, history[cell] != 0);
                    switched += expectedHistory[cell] != 0;
                    if (reference[cell].primaryVarsMeaningWater() != expectedMeaning) {
                        OPM_THROW(std::runtime_error, std::format(
                            "Synthetic Newton validation did not exercise transition {} at cell {}",
                            scenario, cell));
                    }
                }
                GpuBuffer<DevicePV> input(deviceInput), candidate(count);
                GpuBuffer<Scalar> delta(raw), previousDelta(old);
                GpuBuffer<std::uint8_t> previousSwitch(history), candidateSwitch(count);
                GpuBuffer<std::uint32_t> status(std::vector<std::uint32_t>(4, 0));
                blackoilNewtonUpdateKernel<DeviceTag><<<1, 32, 0, bridge.stream()>>>(
                    bridge.flowProblemView(), &bridge.deviceFluidSystem(),
                    GpuView<const DevicePV>(input.data(), count), make_view(candidate),
                    delta.data(), previousDelta.data(), make_view(previousSwitch),
                    make_view(candidateSwitch), status.data(), params, relaxation, useSOR, true);
                OPM_GPU_SAFE_CALL(cudaGetLastError());
                // These tiny explicit validation downloads have no place in the
                // ordinary resident path. Wait only for the validation stream.
                OPM_GPU_SAFE_CALL(cudaStreamSynchronize(bridge.stream()));
                std::vector<DevicePV> result(count);
                std::vector<Scalar> resultDelta(count * numEq), resultPrevious(count * numEq);
                std::vector<std::uint8_t> resultHistory(count);
                std::array<std::uint32_t, 4> resultStatus{};
                candidate.copyToHost(result);
                delta.copyToHost(resultDelta);
                previousDelta.copyToHost(resultPrevious);
                candidateSwitch.copyToHost(resultHistory);
                status.copyToHost(resultStatus.data(), resultStatus.size());
                downloadCalls += 5;
                downloadBytes += count * (sizeof(DevicePV) + 2 * numEq * sizeof(Scalar)
                                          + sizeof(std::uint8_t)) + sizeof(resultStatus);
                const auto fail = [&](unsigned cell, const char* field) {
                    OPM_THROW(std::runtime_error, std::format(
                        "Synthetic GPU Newton mismatch: transition {}, cell {}, relaxation {}, "
                        "mode {}, field {}", scenario, cell, relaxation, useSOR ? "SOR" : "Dampen", field));
                };
                if (resultStatus[0] || resultStatus[2] != switched) {
                    fail(resultStatus[0] ? resultStatus[0] - 1 : 0, "status");
                }
                for (unsigned cell = 0; cell < count; ++cell) {
                    const auto& actual = result[cell];
                    const auto& expected = reference[cell];
                    if (actual.primaryVarsMeaningWater() != expected.primaryVarsMeaningWater()
                        || actual.primaryVarsMeaningPressure() != expected.primaryVarsMeaningPressure()
                        || actual.primaryVarsMeaningGas() != expected.primaryVarsMeaningGas()
                        || actual.primaryVarsMeaningBrine() != expected.primaryVarsMeaningBrine()
                        || actual.primaryVarsMeaningSolvent() != expected.primaryVarsMeaningSolvent()
                        || actual.pvtRegionIndex() != expected.pvtRegionIndex()
                        || actual.capillaryPressureFactor() != expected.capillaryPressureFactor()
                        || actual.pressureScale() != expected.pressureScale()
                        || resultHistory[cell] != expectedHistory[cell]) {
                        fail(cell, "metadata/history");
                    }
                    for (unsigned eq = 0; eq < numEq; ++eq) {
                        Scalar value = actual[eq], ref = expected[eq];
                        const bool pressure = eq == Indices::pressureSwitchIdx;
                        if (pressure) {
                            value *= actual.pressureScale();
                            ref *= expected.pressureScale();
                        }
                        const Scalar atol = pressure ? Scalar{1e-4}
                            : eq == Indices::temperatureIdx ? Scalar{1e-8} : Scalar{1e-12};
                        if (!std::isfinite(value) || std::abs(value - ref) > atol + Scalar{1e-10} * std::abs(ref)) {
                            fail(cell, "primary variable");
                        }
                        const auto offset = cell * numEq + eq;
                        if (resultPrevious[offset] != raw[offset]) {
                            fail(cell, "unrelaxed correction history");
                        }
                        if (std::abs(resultDelta[offset] - correction[cell][eq])
                            > Scalar{1e-12} + Scalar{1e-14} * std::abs(correction[cell][eq])) {
                            fail(cell, "stabilized correction");
                        }
                    }
                }
                ++batches;
            }
        }
    }
    OpmLog::info(std::format(
        "[GPU Newton synthetic validation] {} batches passed (Sw/Rsw/Rvw, Dampen/SOR, "
        "relaxation 0.65/1); validation-only allocations={} bytes={}, uploads={} bytes={}, "
        "downloads={} bytes={}", batches, transfers.allocations, transfers.allocationBytes,
        transfers.hostToDeviceCalls, transfers.hostToDeviceBytes, downloadCalls, downloadBytes));
}


// Diagnostic only: sample public PVT stages and the elementary operations
// used by the low-temperature Spycher-Pruess fugacity calculation. These
// values localize a discrepancy; they do not replace any production formula.
template<class Pvt, class Scalar>
OPM_HOST_DEVICE std::array<Scalar, 20> newtonRswDiagnosticValues(
    const Pvt& pvt, unsigned region, Scalar temperature, Scalar pressure)
{
    using CO2 = typename Pvt::CO2;
    using Binary = typename Pvt::BinaryCoeffBrineCO2;
    std::array<Scalar, 20> out{};
    out[0] = temperature;
    out[1] = pressure;
    out[2] = pvt.salinity(region);
    out[3] = CO2::gasDensity(pvt.getParams(), temperature, pressure, true);
    out[4] = 1 / (out[3] / CO2::molarMass()) * 1e6;
    out[5] = pow(temperature, Scalar{1.5});
    out[6] = (out[4] + Scalar{27.8}) / out[4];
    out[7] = log(out[6]);
    out[8] = out[4] / (out[4] - Scalar{27.8});
    out[9] = log(out[8]);
    out[10] = (pressure / Scalar{1e5}) * out[4]
        / ((Opm::IdealGas<Scalar>::R * Scalar{10}) * temperature);
    out[11] = log(out[10]);
    out[12] = Binary::fugacityCoefficientCO2(pvt.getParams(), temperature, pressure,
                                           Scalar{0}, false, true, true);
    out[13] = Binary::fugacityCoefficientH2O(pvt.getParams(), temperature, pressure,
                                           Scalar{0}, false, true, true);
    Binary::calculateMoleFractions(pvt.getParams(), temperature, pressure, out[2], -1,
                                  out[14], out[15], pvt.getActivityModel(), true);
    out[16] = pvt.saturatedGasDissolutionFactor(region, temperature, pressure, Scalar{0});
    out[17] = pvt.waterReferenceDensity(region);
    out[18] = pvt.gasReferenceDensity(region);
    out[19] = pvt.getActivityModel();
    return out;
}

template<class FluidSystem, class Scalar>
__global__ void newtonRswDiagnosticKernel(const FluidSystem* fluidSystem, unsigned region,
                                        Scalar temperature, Scalar pressure, Scalar* out)
{
    const auto values = newtonRswDiagnosticValues(fluidSystem->waterPvt(), region,
                                                temperature, pressure);
    for (unsigned i = 0; i < values.size(); ++i) {
        out[i] = values[i];
    }
}

template<class CpuTypeTag, class Bridge>
void diagnoseBlackoilNewtonRsw(Bridge& bridge, unsigned region,
                              typename Bridge::Scalar temperature,
                              typename Bridge::Scalar pressure)
{
    using Scalar = typename Bridge::Scalar;
    using FluidSystem = GetPropType<CpuTypeTag, Properties::FluidSystem>;
    const auto& pvt = FluidSystem::waterPvt().template getRealPvt<WaterPvtApproach::BrineCo2>();
    const auto cpu = newtonRswDiagnosticValues(pvt, region, temperature, pressure);
    GpuBuffer<Scalar> output(cpu.size());
    newtonRswDiagnosticKernel<<<1, 1, 0, bridge.stream()>>>(
        &bridge.deviceFluidSystem(), region, temperature, pressure, output.data());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    OPM_GPU_SAFE_CALL(cudaStreamSynchronize(bridge.stream()));
    std::array<Scalar, 20> gpu{};
    output.copyToHost(gpu.data(), gpu.size());
    constexpr std::array<const char*, 20> names{
        "temperature", "pressure", "salinity", "CO2 density", "molar volume", "pow(T,1.5)",
        "log argument (V+b)/V", "log((V+b)/V)", "log argument V/(V-b)", "log(V/(V-b))",
        "log argument pV/RT", "log(pV/RT)", "CO2 fugacity", "H2O fugacity", "liquid CO2 mole fraction",
        "gas H2O mole fraction", "saturated Rsw", "brine reference density", "CO2 reference density", "activity model"};
    OpmLog::info(std::format(
        "[GPU Newton Rsw diagnostic] region={} validation-only allocations=1 bytes={} downloads=1 bytes={}",
        region, sizeof(gpu), sizeof(gpu)));
    for (unsigned i = 0; i < cpu.size(); ++i) {
        OpmLog::info(std::format(
            "[GPU Newton Rsw diagnostic] {} CPU={:.17g} GPU={:.17g} difference={:.17g}",
            names[i], cpu[i], gpu[i], gpu[i] - cpu[i]));
    }
}

// Keep the Newton kernel (including its PVT and material-law call graph) in
// this translation unit. CMake disables FP contraction only for this source;
// GPU property evaluation retains its existing compiler settings.
using NewtonDeviceTag = Properties::TTag::FlowGasWaterEnergyDeviceTypeTag<GpuView>;
using NewtonCpuTag = Properties::TTag::FlowGasWaterEnergyProblem;
using NewtonCpuAssemblyTag = Properties::TTag::FlowGasWaterEnergyProblemGPU;
using NewtonBridge = GpuFlowGasWaterEnergyBridge<NewtonCpuTag, NewtonDeviceTag>;
using NewtonAssemblyBridge = GpuFlowGasWaterEnergyBridge<NewtonCpuAssemblyTag, NewtonDeviceTag>;

template void launchBlackoilNewtonUpdate<NewtonDeviceTag, NewtonBridge>(
    NewtonBridge&, const BlackoilNewtonParams<NewtonBridge::Scalar>&,
    NewtonBridge::Scalar, bool, bool);
template void launchBlackoilNewtonUpdate<NewtonDeviceTag, NewtonAssemblyBridge>(
    NewtonAssemblyBridge&, const BlackoilNewtonParams<NewtonAssemblyBridge::Scalar>&,
    NewtonAssemblyBridge::Scalar, bool, bool);

template void validateBlackoilNewtonBranches<NewtonCpuTag, NewtonBridge,
    GetPropType<NewtonCpuTag, Properties::Problem>,
    GetPropType<NewtonCpuTag, Properties::SolutionVector>>(
    NewtonBridge&, const GetPropType<NewtonCpuTag, Properties::Problem>&,
    const GetPropType<NewtonCpuTag, Properties::SolutionVector>&,
    const BlackoilNewtonParams<NewtonBridge::Scalar>&);
template void validateBlackoilNewtonBranches<NewtonCpuAssemblyTag, NewtonAssemblyBridge,
    GetPropType<NewtonCpuAssemblyTag, Properties::Problem>,
    GetPropType<NewtonCpuAssemblyTag, Properties::SolutionVector>>(
    NewtonAssemblyBridge&, const GetPropType<NewtonCpuAssemblyTag, Properties::Problem>&,
    const GetPropType<NewtonCpuAssemblyTag, Properties::SolutionVector>&,
    const BlackoilNewtonParams<NewtonAssemblyBridge::Scalar>&);

template void diagnoseBlackoilNewtonRsw<NewtonCpuTag, NewtonBridge>(
    NewtonBridge&, unsigned, NewtonBridge::Scalar, NewtonBridge::Scalar);
template void diagnoseBlackoilNewtonRsw<NewtonCpuAssemblyTag, NewtonAssemblyBridge>(
    NewtonAssemblyBridge&, unsigned, NewtonAssemblyBridge::Scalar, NewtonAssemblyBridge::Scalar);

} // namespace Opm::gpuistl
