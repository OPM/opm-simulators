/*
  Copyright 2024 SINTEF AS
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

#define BOOST_TEST_MODULE TestGpuTwoPhaseMaterial

#include <cuda_runtime.h>

#include <boost/test/unit_test.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>

#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <array>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <iostream>



// these types are taken from SPE11C
const int numPhases = 3;

using Scalar = float;

const int staticSize = 0;

// using SPE11CEvaluation = Opm::DenseAd::Evaluation<Scalar, numPhases, staticSize>;

using ValueVector = std::vector<Scalar>;
using GPUBuffer = Opm::gpuistl::GpuBuffer<Scalar>;
using GPUView = Opm::gpuistl::GpuView<const Scalar>;

using SatOnlyFluidState = Opm::SimpleModularFluidState<Scalar,
                                                  numPhases,
                                                  3,
                                                  /*FluidSystem, just a placeholder for now*/ long,
                                                  /*storePressure=*/false,
                                                  /*storeTemperature=*/false,
                                                  /*storeComposition=*/false,
                                                  /*storeFugacity=*/false,
                                                  /*storeSaturation=*/true,
                                                  /*storeDensity=*/false,
                                                  /*storeViscosity=*/false,
                                                  /*storeEnthalpy=*/false>;

using TwoPhaseTraitsT = Opm::TwoPhaseMaterialTraits<Scalar, 1, 2>; // is this ordering reasonable?
using ThreePhaseTraitsT = Opm::ThreePhaseMaterialTraits<Scalar, 1, 2, 0>; // is this ordering reasonable? (i have also tried 1, 2, 3)

using CPUPLParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TwoPhaseTraitsT>;
using constGPUPiecewiseLinearParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TwoPhaseTraitsT, const GPUBuffer>;
// using GPUBufferParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TwoPhaseTraitsT, GPUBuffer>;
using ViewGPUPiecewiseLinearParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TwoPhaseTraitsT, GPUView>;
using CPUPLMaterialLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraitsT, CPUPLParams>;
using GPUPLTwoPhaseLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraitsT, ViewGPUPiecewiseLinearParams>;

using CPUTwoPhaseLawParams = Opm::EclTwoPhaseMaterialParams<ThreePhaseTraitsT, CPUPLParams, CPUPLParams, CPUPLParams>;
using CPUTwoPhaseLaw = Opm::EclTwoPhaseMaterial<ThreePhaseTraitsT, CPUPLMaterialLaw, CPUPLMaterialLaw, CPUPLMaterialLaw>;
using GPUTwoPhaseLawParams = Opm::EclTwoPhaseMaterialParams<ThreePhaseTraitsT, ViewGPUPiecewiseLinearParams, ViewGPUPiecewiseLinearParams, ViewGPUPiecewiseLinearParams, Opm::gpuistl::ViewPointer>;
using GPUTwoPhaseLaw = Opm::EclTwoPhaseMaterial<ThreePhaseTraitsT, GPUPLTwoPhaseLaw, GPUPLTwoPhaseLaw, GPUPLTwoPhaseLaw, Opm::gpuistl::ViewPointer>;

__global__ void gpuCapillaryPressure(GPUTwoPhaseLawParams params, Scalar* values){
    // createa a simplified fluidstate object on the GPU from inside a gpu kernel
    SatOnlyFluidState fluidState;
    fluidState.setSaturation(0, 0.6); // water
    fluidState.setSaturation(1, 0.0); // oil
    fluidState.setSaturation(2, 0.4); // gas

    // use the created fluidstate to create a materialLaw query
    GPUTwoPhaseLaw::capillaryPressures(values, params, fluidState);
}

BOOST_AUTO_TEST_CASE(TestSimpleInterpolation)
{
    CPUPLParams cpuParams;

    // create some vectors that will be used for linear interpolation
    ValueVector cx = {0.0, 0.5, 1.0};
    ValueVector cy = {0.0, 0.9, 1.0};

    // place the data of linear interpolation on the GPU in a GPUBuffer
    const GPUBuffer gx(cx);
    const GPUBuffer gy(cy);

    // Create a PiecewiseLinearTwoPhaseMaterial params object with the interpolation data on the CPU
    cpuParams.setPcnwSamples(cx, cy);
    cpuParams.setKrwSamples(cx, cy);
    cpuParams.setKrnSamples(cx, cy);
    cpuParams.finalize();
    std::shared_ptr<CPUPLParams> CPUPLParamsPointer = std::make_shared<CPUPLParams>(cpuParams);

    // do the same on the GPU
    constGPUPiecewiseLinearParams gpuBufferParams(gx, gy, gx, gy, gx, gy);

    // make a view of the GPU parameters
    ViewGPUPiecewiseLinearParams gpuViewParams = Opm::gpuistl::make_view<TwoPhaseTraitsT, const GPUBuffer, GPUView>(gpuBufferParams);
    auto gpuViewParamsPointer = Opm::gpuistl::ViewPointer<ViewGPUPiecewiseLinearParams>(gpuViewParams);

    // make a twoPhaseLaw with the parameters on the CPU
    auto cpuTwoPhaseLawParams = CPUTwoPhaseLawParams();
    cpuTwoPhaseLawParams.setGasOilParams(CPUPLParamsPointer);
    cpuTwoPhaseLawParams.setOilWaterParams(CPUPLParamsPointer);
    cpuTwoPhaseLawParams.setGasWaterParams(CPUPLParamsPointer);
    cpuTwoPhaseLawParams.setApproach(Opm::EclTwoPhaseApproach::GasWater);
    cpuTwoPhaseLawParams.finalize();

    // make a twoPhaseLaw with the parameters on the GPU
    auto gpuTwoPhaseLawParams = GPUTwoPhaseLawParams();
    //gpuTwoPhaseLawParams.setGasOilParams(&gpuViewParams);
    gpuTwoPhaseLawParams.setGasOilParams(gpuViewParamsPointer);
    gpuTwoPhaseLawParams.setOilWaterParams(gpuViewParamsPointer);
    gpuTwoPhaseLawParams.setGasWaterParams(gpuViewParamsPointer);
    gpuTwoPhaseLawParams.setApproach(Opm::EclTwoPhaseApproach::GasWater);
    gpuTwoPhaseLawParams.finalize();

    // Create test data and helper variables
    ValueVector testXs = {-1.0, 0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.1};
    const size_t VEC_BYTES = 3*sizeof(Scalar);
    Scalar* gpuValues;
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuValues, VEC_BYTES));

    for (Scalar x_i : testXs){

        // initialize some data for computation
        Scalar arr[3] = {x_i, x_i, x_i };
        Scalar cpuResult[3] = {x_i, x_i, x_i };
        Scalar gpuResult[3];

        // compute capillary pressure on the CPU
        SatOnlyFluidState fluidState;
        fluidState.setSaturation(0, 0.6); // water
        fluidState.setSaturation(1, 0.0); // oil
        fluidState.setSaturation(2, 0.4); // gas
        CPUTwoPhaseLaw::capillaryPressures(cpuResult, cpuTwoPhaseLawParams, fluidState);

        // compute capillary pressure on the GPU
        OPM_GPU_SAFE_CALL(cudaMemcpy(gpuValues, arr, VEC_BYTES, cudaMemcpyHostToDevice));
        gpuCapillaryPressure<<<1,1>>>(gpuTwoPhaseLawParams, gpuValues);
        OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
        OPM_GPU_SAFE_CALL(cudaMemcpy(gpuResult, gpuValues, VEC_BYTES, cudaMemcpyDeviceToHost));

        // check for each phase that the state of the input/output vector is the same
        printf("%f %f %f === %f %f %f\n", gpuResult[0], gpuResult[1], gpuResult[2], cpuResult[0], cpuResult[1], cpuResult[2]);
        BOOST_CHECK(gpuResult[0] == cpuResult[0]);
        BOOST_CHECK(gpuResult[1] == cpuResult[1]);
        BOOST_CHECK(gpuResult[2] == cpuResult[2]);
    }
    OPM_GPU_SAFE_CALL(cudaFree(gpuValues));
}
