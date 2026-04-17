/*
  Copyright TODO ADD YEAR AND NAME OF AUTHOR

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

#define BOOST_TEST_MODULE TestGpuEclThermalLawManager

#include <boost/test/unit_test.hpp>

#include <cuda_runtime.h>

#include <opm/material/common/EnsureFinalized.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/material/thermal/EclSpecrockLaw.hpp>
#include <opm/material/thermal/EclSpecrockLawParams.hpp>
#include <opm/material/thermal/EclThconrLaw.hpp>
#include <opm/material/thermal/EclThconrLawParams.hpp>

#include <opm/simulators/flow/GpuEclThermalLawManager.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <array>
#include <vector>

namespace {

using Scalar = double;

// A minimal stand-in fluid system. The thermal-conduction law only ever
// asks the fluid system whether the gas phase is active and what its
// phase index is. Both are compile-time constants.
struct DummyFluidSystem {
    static constexpr int numPhases = 3;
    static constexpr int gasPhaseIdx = 2;
    OPM_HOST_DEVICE static constexpr bool phaseIsActive(int phaseIdx)
    {
        return phaseIdx == gasPhaseIdx;
    }
};

// A minimal stand-in fluid state with just enough surface area to feed
// EclSpecrockLaw / EclThconrLaw / their GPU counterparts.
template <class ScalarT>
struct DummyFluidState {
    using ValueType = ScalarT;
    ScalarT temperature_ = 0;
    ScalarT gasSaturation_ = 0;

    OPM_HOST_DEVICE ScalarT temperature(int /*phaseIdx*/) const { return temperature_; }
    OPM_HOST_DEVICE ScalarT saturation(int /*phaseIdx*/) const { return gasSaturation_; }
    OPM_HOST_DEVICE DummyFluidSystem fluidSystem() const { return DummyFluidSystem{}; }
};

using GpuSpecrockParamsBuf = Opm::EclSpecrockLawParams<Scalar, Opm::gpuistl::GpuView>;
using GpuSpecrockParamsView = Opm::EclSpecrockLawParams<Scalar, Opm::gpuistl::GpuView>;
using GpuSpecrockLawView = Opm::EclSpecrockLaw<Scalar, GpuSpecrockParamsView>;
using GpuThconrLaw = Opm::EclThconrLaw<Scalar, DummyFluidSystem>;

using ManagerCpu = Opm::EclThermalLaw::GpuManager<Scalar, DummyFluidSystem>;
using ManagerBuf
    = Opm::EclThermalLaw::GpuManager<Scalar, DummyFluidSystem, Opm::gpuistl::GpuBuffer, Opm::gpuistl::GpuView>;
using ManagerView
    = Opm::EclThermalLaw::GpuManager<Scalar, DummyFluidSystem, Opm::gpuistl::GpuView, Opm::gpuistl::GpuView>;

// Kernel: per-cell, evaluate solid internal energy from the law manager
// at the supplied temperature, and total thermal conductivity at the
// supplied gas saturation.
__global__ void evaluateThermalKernel(ManagerView view,
                                      const Scalar* temperatures,
                                      const Scalar* gasSaturations,
                                      Scalar* outSolidEnergy,
                                      Scalar* outConductivity,
                                      std::size_t numCells)
{
    const std::size_t i =
        static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numCells) {
        return;
    }
    DummyFluidState<Scalar> fs;
    fs.temperature_ = temperatures[i];
    fs.gasSaturation_ = gasSaturations[i];

    const auto solidEnergyParams = view.solidEnergyLawParams(static_cast<unsigned>(i));
    outSolidEnergy[i] =
        GpuSpecrockLawView::solidInternalEnergy(solidEnergyParams, fs);

    const auto thermalConductionParams =
        view.thermalConductionLawParams(static_cast<unsigned>(i));
    outConductivity[i] =
        GpuThconrLaw::thermalConductivity(thermalConductionParams, fs);
}

} // namespace

BOOST_AUTO_TEST_CASE(SpecrockAndThconrEvaluateMatchesCpu)
{
    // -- Build a synthetic two-region SPECROCK plus per-cell THCONR setup
    //    on the host. Cells 0..2 belong to region 0, cells 3..5 to region 1.
    constexpr std::size_t numCells = 6;

    const std::vector<Scalar> regionTemperatures0 = { 280.0, 320.0, 360.0, 400.0 };
    const std::vector<Scalar> regionEnergies0     = { 0.0,   4.0e7, 9.0e7, 1.5e8 };
    const std::vector<Scalar> regionTemperatures1 = { 280.0, 350.0, 420.0 };
    const std::vector<Scalar> regionEnergies1     = { 0.0,   3.5e7, 8.4e7 };

    // Build a host-side CPU manager directly from the basic constructor.
    Opm::EclSpecrockLawParams<Scalar> region0Params;
    region0Params.setSamples(regionTemperatures0, regionEnergies0);
    Opm::EclSpecrockLawParams<Scalar> region1Params;
    region1Params.setSamples(regionTemperatures1, regionEnergies1);

    std::vector<Opm::EclSpecrockLawParams<Scalar>> cpuSolidEnergyParams;
    cpuSolidEnergyParams.push_back(region0Params);
    cpuSolidEnergyParams.push_back(region1Params);

    std::vector<int> cpuElementToRegion(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        cpuElementToRegion[i] = (i < 3) ? 0 : 1;
    }

    std::vector<Opm::EclThconrLawParams<Scalar>> cpuThermalConductionParams(numCells);
    for (std::size_t i = 0; i < numCells; ++i) {
        cpuThermalConductionParams[i].setReferenceTotalThermalConductivity(
            100.0 + 5.0 * static_cast<Scalar>(i));
        cpuThermalConductionParams[i].setDTotalThermalConductivity_dSg(
            0.20 + 0.01 * static_cast<Scalar>(i));
        cpuThermalConductionParams[i].finalize();
    }

    ManagerCpu cpuManager(std::move(cpuSolidEnergyParams),
                          std::move(cpuElementToRegion),
                          std::move(cpuThermalConductionParams));

    // -- Upload to the GPU and create a non-owning view.
    auto bufManager = Opm::gpuistl::copy_to_gpu(cpuManager);
    auto viewManager = Opm::gpuistl::make_view(bufManager);

    // -- Pick a per-cell test temperature and gas saturation each.
    std::array<Scalar, numCells> hostTemperatures = {
        290.0, 310.0, 355.0, 285.0, 360.0, 410.0
    };
    std::array<Scalar, numCells> hostGasSaturations = {
        0.0, 0.1, 0.25, 0.5, 0.75, 0.99
    };

    // -- Run the kernel.
    Scalar* dT = nullptr;
    Scalar* dSg = nullptr;
    Scalar* dEnergy = nullptr;
    Scalar* dConductivity = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&dT, numCells * sizeof(Scalar)));
    OPM_GPU_SAFE_CALL(cudaMalloc(&dSg, numCells * sizeof(Scalar)));
    OPM_GPU_SAFE_CALL(cudaMalloc(&dEnergy, numCells * sizeof(Scalar)));
    OPM_GPU_SAFE_CALL(cudaMalloc(&dConductivity, numCells * sizeof(Scalar)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(dT, hostTemperatures.data(),
                                 numCells * sizeof(Scalar), cudaMemcpyHostToDevice));
    OPM_GPU_SAFE_CALL(cudaMemcpy(dSg, hostGasSaturations.data(),
                                 numCells * sizeof(Scalar), cudaMemcpyHostToDevice));

    evaluateThermalKernel<<<1, 32>>>(viewManager, dT, dSg, dEnergy, dConductivity, numCells);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());

    std::array<Scalar, numCells> hostEnergyResult {};
    std::array<Scalar, numCells> hostConductivityResult {};
    OPM_GPU_SAFE_CALL(cudaMemcpy(hostEnergyResult.data(), dEnergy,
                                 numCells * sizeof(Scalar), cudaMemcpyDeviceToHost));
    OPM_GPU_SAFE_CALL(cudaMemcpy(hostConductivityResult.data(), dConductivity,
                                 numCells * sizeof(Scalar), cudaMemcpyDeviceToHost));

    OPM_GPU_SAFE_CALL(cudaFree(dT));
    OPM_GPU_SAFE_CALL(cudaFree(dSg));
    OPM_GPU_SAFE_CALL(cudaFree(dEnergy));
    OPM_GPU_SAFE_CALL(cudaFree(dConductivity));

    // -- Compare against the CPU manager.
    for (std::size_t i = 0; i < numCells; ++i) {
        DummyFluidState<Scalar> fs;
        fs.temperature_ = hostTemperatures[i];
        fs.gasSaturation_ = hostGasSaturations[i];

        const auto& cpuSolidParams =
            cpuManager.solidEnergyLawParams(static_cast<unsigned>(i));
        const Scalar cpuEnergy =
            Opm::EclSpecrockLaw<Scalar>::solidInternalEnergy(cpuSolidParams, fs);
        BOOST_CHECK_CLOSE(cpuEnergy, hostEnergyResult[i], 1e-9);

        const auto& cpuConductionParams =
            cpuManager.thermalConductionLawParams(static_cast<unsigned>(i));
        const Scalar cpuConductivity =
            GpuThconrLaw::thermalConductivity(cpuConductionParams, fs);
        BOOST_CHECK_CLOSE(cpuConductivity, hostConductivityResult[i], 1e-9);
    }
}
