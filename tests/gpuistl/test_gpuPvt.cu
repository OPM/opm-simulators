#include <config.h>

#define BOOST_TEST_MODULE TestGpuPvt

#include <boost/test/unit_test.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/components/CO2Tables.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/SimpleHuDuanH2O.hpp>
#include <opm/material/components/BrineDynamic.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/fluidsystems/blackoilpvt/Co2GasPvt.hpp>
#include <opm/input/eclipse/EclipseState/Co2StoreConfig.hpp>

#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>

#include <cuda_runtime.h>
#include <vector>
#include <utility>
#include <cmath>

using Evaluation = Opm::DenseAd::Evaluation<double, 3>;
using GpuB = const Opm::gpuistl::GpuBuffer<double>;
using GpuV = Opm::gpuistl::GpuView<const double>;

using GpuTab = Opm::UniformTabulated2DFunction<double, GpuV>;

using GpuBufCo2Tables = Opm::CO2Tables<double, GpuB>;
using GpuViewCO2Tables = Opm::CO2Tables<double, GpuV>;
using GpuCO2 = Opm::CO2<double, GpuViewCO2Tables>;

using HuDuan = Opm::SimpleHuDuanH2O<double>;
using BrineDyn = Opm::BrineDynamic<double, HuDuan>;

using CpuBrine_CO2 = Opm::BinaryCoeff::Brine_CO2<double, HuDuan, Opm::CO2<double>>;
using GpuBrine_CO2 = Opm::BinaryCoeff::Brine_CO2<double, HuDuan, GpuCO2>;

using CpuCo2Pvt = Opm::Co2GasPvt<double>;
using GpuBufCo2Pvt = Opm::Co2GasPvt<double, GpuBufCo2Tables, GpuB>;
using GpuViewCo2Pvt = Opm::Co2GasPvt<double, GpuViewCO2Tables, GpuV>;

namespace {

/*
    This file contains unit tests for Pvt objects and function on the GPU, with additional helper classes
*/

const double ABS_TOL = 1e-6;

struct Fixture {
    Fixture(){
        temp = Evaluation(290.5);
        pressure = Evaluation(200000.0);

        gpuComputedResultOnCpu = 0.0;

        // move pvt evaluations to gpu
        OPM_GPU_SAFE_CALL(cudaMalloc(&gpuTemp, sizeof(Evaluation)));
        OPM_GPU_SAFE_CALL(cudaMemcpy(gpuTemp, &temp, sizeof(Evaluation), cudaMemcpyHostToDevice));
        OPM_GPU_SAFE_CALL(cudaMalloc(&gpuPressure, sizeof(Evaluation)));
        OPM_GPU_SAFE_CALL(cudaMemcpy(gpuPressure, &pressure, sizeof(Evaluation), cudaMemcpyHostToDevice));
    }
    ~Fixture(){
        OPM_GPU_SAFE_CALL(cudaFree(gpuTemp));
        OPM_GPU_SAFE_CALL(cudaFree(gpuPressure));
    }

    Evaluation temp; // [K]
    Evaluation pressure; // [Pa]
    Evaluation* gpuTemp; // [K]
    Evaluation* gpuPressure; // [Pa]

    double gpuComputedResultOnCpu;

    Opm::CO2Tables<double, std::vector<double>> co2Tables;
};

// Kernel to evaluate a 2D function on the GPU
__global__ void gpuEvaluateUniformTabulated2DFunction(GpuTab gpuTab, Evaluation* inputX, Evaluation* inputY, double* result) {
    *result = gpuTab.eval(*inputX, *inputY, true).value();
}

// Kernel using a CO2 object on the GPU
__global__ void gpuCO2GasViscosity(GpuViewCO2Tables gpuViewCo2Tables, Evaluation* temp, Evaluation* pressure, double* result) {
    *result = GpuCO2::gasViscosity<Evaluation>(gpuViewCo2Tables, *temp, *pressure, true).value();
}

// Kernel using a SimpleHuDuanH20 object on a GPU
__global__ void huDuanLiquidDensity(Evaluation* temp, Evaluation* pressure, double* result) {
    *result = HuDuan::liquidDensity<Evaluation>(*temp, *pressure, true).value();
}

// Kernel using a BrineDynamic object on a GPU
__global__ void brineDynamicLiquidEnthalpy(Evaluation* temp, Evaluation* pressure, Evaluation* salinity, double* result) {
    *result = BrineDyn::liquidEnthalpy<Evaluation>(*temp, *pressure, *salinity).value();
}

// Kernel using a Brine_CO2 object on a GPU
__global__ void brineCO2GasDiffCoeff(GpuViewCO2Tables co2tables, Evaluation* temp, Evaluation* pressure, double* result) {
    *result = GpuBrine_CO2::gasDiffCoeff<Evaluation, GpuViewCO2Tables>(co2tables, *temp, *pressure, true).value();
}

// Kernel using a Co2GasPvt object on a GPU
__global__ void co2GasPvtInternalEnergy(GpuViewCo2Pvt gpuViewCo2Pvt, Evaluation* temp, Evaluation* pressure, double* result) {
    *result = gpuViewCo2Pvt.internalEnergy<Evaluation>(1, *temp, *pressure, Evaluation(0.4), Evaluation(0.0)).value();
}

// Helper function to launch a kernel and retrieve the result on the CPU to reduce code duplicatoin
template <typename KernelFunc, typename... Args>
double launchKernelAndRetrieveResult(KernelFunc kernel, Args... args) {
    double* resultOnGpu;
    double gpuComputedResultOnCpu;

    // Allocate memory for the result on the GPU
    OPM_GPU_SAFE_CALL(cudaMalloc(&resultOnGpu, sizeof(double)));

    // Launch the kernel
    kernel<<<1, 1>>>(args..., resultOnGpu);

    // Check for any errors in kernel launch
    OPM_GPU_SAFE_CALL(cudaPeekAtLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    // Retrieve the result from the GPU to the CPU
    OPM_GPU_SAFE_CALL(cudaMemcpy(&gpuComputedResultOnCpu, resultOnGpu, sizeof(double), cudaMemcpyDeviceToHost));

    // Free allocated GPU memory
    OPM_GPU_SAFE_CALL(cudaFree(resultOnGpu));

    return gpuComputedResultOnCpu;
}

bool compareSignificantDigits(double a, double b, int significantDigits) {
    // Handle the case where both values are exactly equal
    if (a == b) {
        return true;
    }

    // Calculate the relative error
    double relativeError = std::abs(a - b) / std::max(std::abs(a), std::abs(b));

    // Compute the number of matching digits
    double digitsMatched = -std::log10(relativeError);

    // Return true if they match the required number of significant digits
    return digitsMatched >= significantDigits;
}

} // END EMPTY NAMESPACE

// Test case for evaluating a tabulated 2D function on both CPU and GPU
BOOST_FIXTURE_TEST_CASE(TestEvaluateUniformTabulated2DFunctionOnGpu, Fixture) {
    // Example tabulated data (2D)
    std::vector<std::vector<double>> tabData = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};

    // CPU-side function definition
    Opm::UniformTabulated2DFunction<double> cpuTab(1.0, 6.0, 3, 1.0, 6.0, 2, tabData);

    // Move data to GPU buffer and create a view for GPU operations
    Opm::UniformTabulated2DFunction<double, GpuB> gpuBufTab = Opm::gpuistl::move_to_gpu<double, GpuB>(cpuTab);
    GpuTab gpuViewTab = Opm::gpuistl::make_view<GpuV>(gpuBufTab);

    // Evaluation points on the CPU
    Evaluation a(2.3);
    Evaluation b(4.5);

    // Allocate GPU memory for the Evaluation inputs
    Evaluation* gpuA = nullptr;
    Evaluation* gpuB = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuA, sizeof(Evaluation)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuA, &a, sizeof(Evaluation), cudaMemcpyHostToDevice));
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuB, sizeof(Evaluation)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuB, &b, sizeof(Evaluation), cudaMemcpyHostToDevice));

    gpuComputedResultOnCpu = launchKernelAndRetrieveResult(gpuEvaluateUniformTabulated2DFunction, gpuViewTab, gpuA, gpuB);

    // Free allocated GPU memory
    OPM_GPU_SAFE_CALL(cudaFree(gpuA));
    OPM_GPU_SAFE_CALL(cudaFree(gpuB));

    // Verify that the CPU and GPU results match within a reasonable tolerance
    const double cpuComputedResult = cpuTab.eval(a, b, true).value();
    BOOST_CHECK(std::fabs(gpuComputedResultOnCpu - cpuComputedResult) < ABS_TOL);
}

// Test case evaluating CO2 pvt properties on CPU and GPU
BOOST_FIXTURE_TEST_CASE(TestUseCO2OnGpu, Fixture) {

    // use the CO2 tables to aquire the viscosity at 290[K] and 2e5[Pa]
    double viscosityReference = Opm::CO2<double, Opm::CO2Tables<double, std::vector<double>>>::gasViscosity<Evaluation>(co2Tables, temp, pressure, true).value();

    GpuBufCo2Tables gpuBufCo2Table = Opm::gpuistl::move_to_gpu<double, std::vector<double>, GpuB>(co2Tables);
    GpuViewCO2Tables gpuViewCo2Table = Opm::gpuistl::make_view<GpuV>(gpuBufCo2Table);

    gpuComputedResultOnCpu = launchKernelAndRetrieveResult(gpuCO2GasViscosity, gpuViewCo2Table, gpuTemp, gpuPressure);

    // Verify that the CPU and GPU results match within a reasonable tolerance
    BOOST_CHECK(std::fabs(gpuComputedResultOnCpu - viscosityReference) < ABS_TOL);
}

// Test case evaluating pvt values for SimpleHuDuanH20 on a GPU and CPU
BOOST_FIXTURE_TEST_CASE(TestUseH2OOnGpu, Fixture) {

    // use the CO2 tables to aquire the densityReference at 290[K] and 2e5[Pa]
    double densityReference = HuDuan::liquidDensity<Evaluation>(temp, pressure, true).value();

    gpuComputedResultOnCpu = launchKernelAndRetrieveResult(huDuanLiquidDensity, gpuTemp, gpuPressure);

    // Verify that the CPU and GPU results match within a reasonable tolerance
    BOOST_CHECK(std::fabs(gpuComputedResultOnCpu - densityReference) < ABS_TOL);
}

// Test case evaluating pvt values for BrineDynamic on a GPU and CPU
BOOST_FIXTURE_TEST_CASE(TestUseBrineDynamicOnGpu, Fixture) {
    Evaluation salinity(0.1); // [g/Kg]

    // use the CO2 tables to aquire the enthalpyReference at 290[K] and 2e5[Pa]
    double enthalpyReference = BrineDyn::liquidEnthalpy<Evaluation>(temp, pressure, salinity).value();

    // Allocate GPU memory for the Evaluation inputs
    Evaluation* gpuSalinity = nullptr;
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuSalinity, sizeof(Evaluation)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuSalinity, &salinity, sizeof(Evaluation), cudaMemcpyHostToDevice));

    gpuComputedResultOnCpu = launchKernelAndRetrieveResult(brineDynamicLiquidEnthalpy, gpuTemp, gpuPressure, gpuSalinity);

    // Verify that the CPU and GPU results match within a reasonable tolerance
    BOOST_CHECK(std::fabs(gpuComputedResultOnCpu - enthalpyReference) < ABS_TOL);
}

// Test case evaluating pvt values for BrineDynamic on a GPU and CPU
BOOST_FIXTURE_TEST_CASE(TestBrine_CO2OnGPU, Fixture) {

    // use the CO2 tables to aquire the gasDiffCoeffReference at 290[K] and 2e5[Pa]
    double gasDiffCoeffReference = CpuBrine_CO2::gasDiffCoeff<Evaluation>(co2Tables, temp, pressure, true).value();

    // use the CO2 tables to aquire the viscosity at 290[K] and 2e5[Pa]
    double viscosity = Opm::CO2<double, Opm::CO2Tables<double, std::vector<double>>>::gasViscosity<Evaluation>(co2Tables, temp, pressure, true).value();

    GpuBufCo2Tables gpuBufCo2Table = Opm::gpuistl::move_to_gpu<double, std::vector<double>, GpuB>(co2Tables);
    GpuViewCO2Tables gpuViewCo2Table = Opm::gpuistl::make_view<GpuV>(gpuBufCo2Table);

    gpuComputedResultOnCpu = launchKernelAndRetrieveResult(brineCO2GasDiffCoeff, gpuViewCo2Table, gpuTemp, gpuPressure);

    // Verify that the CPU and GPU results match within a reasonable tolerance
    BOOST_CHECK(std::fabs(gpuComputedResultOnCpu - gasDiffCoeffReference) < ABS_TOL);
}

// Test case evaluating pvt values for BrineDynamic on a GPU and CPU
BOOST_FIXTURE_TEST_CASE(TestCo2GasPvt, Fixture) {
    std::vector<double> salinities = {0.2, 0.3, 0.4};

    CpuCo2Pvt cpuCo2Pvt(salinities);
    double internalEnergyReference = cpuCo2Pvt.internalEnergy<Evaluation>(1, temp, pressure, Evaluation(0.4), Evaluation(0.0)).value();

    const GpuBufCo2Pvt gpuBufCo2Pvt = Opm::gpuistl::move_to_gpu<double, GpuBufCo2Tables, GpuB>(cpuCo2Pvt);
    const auto brineReferenceDensityCPUCopy = gpuBufCo2Pvt.getBrineReferenceDensity().asStdVector();
    const GpuViewCo2Pvt gpuViewCo2Pvt = Opm::gpuistl::make_view<GpuV, GpuViewCO2Tables>(gpuBufCo2Pvt);

    gpuComputedResultOnCpu = launchKernelAndRetrieveResult(co2GasPvtInternalEnergy, gpuViewCo2Pvt, gpuTemp, gpuPressure);

    // Verify that the CPU and GPU results match within a reasonable tolerance
    BOOST_CHECK(compareSignificantDigits(gpuComputedResultOnCpu, internalEnergyReference, 6));
}
