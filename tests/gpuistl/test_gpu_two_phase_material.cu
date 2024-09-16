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

using SPE11CEvaluation = Opm::DenseAd::Evaluation<Scalar, numPhases, staticSize>;

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

using TraitsT = Opm::TwoPhaseMaterialTraits<Scalar, 1, 2>;
using CPUParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TraitsT>;
using constGPUPiecewiseLinearParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TraitsT, const GPUBuffer>;
// using GPUBufferParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TraitsT, GPUBuffer>;
using ViewGPUPiecewiseLinearParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TraitsT, GPUView>;
using CPUTPiecewiseLinearMaterialLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TraitsT, CPUParams>;
using GPUPLTwoPhaseLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TraitsT, ViewGPUPiecewiseLinearParams>;

using GPUTwoPhaseLawParams = Opm::EclTwoPhaseMaterialParams<TraitsT, ViewGPUPiecewiseLinearParams, ViewGPUPiecewiseLinearParams, ViewGPUPiecewiseLinearParams, false>;
using GPUTwoPhaseLaw = Opm::EclTwoPhaseMaterial<TraitsT, GPUPLTwoPhaseLaw, GPUPLTwoPhaseLaw, GPUPLTwoPhaseLaw, false>;


__global__ void gpuTwoPhaseSatPcnwWrapper(){
    SatOnlyFluidState fluidState;
    fluidState.setSaturation(0, 0.6); // water
    fluidState.setSaturation(0, 0.0); // oil
    fluidState.setSaturation(0, 0.4); // gas

    
}

BOOST_AUTO_TEST_CASE(TestSimpleInterpolation)
{


    CPUParams cpuParams;
    ViewGPUPiecewiseLinearParams gpuViewParams;

    ValueVector cx = {0.0, 0.5, 1.0};
    ValueVector cy = {0.0, 0.9, 1.0};
    const GPUBuffer gx(cx);
    const GPUBuffer gy(cy);

    cpuParams.setPcnwSamples(cx, cy);
    cpuParams.setKrwSamples(cx, cy);
    cpuParams.setKrnSamples(cx, cy);
    cpuParams.finalize();

    constGPUPiecewiseLinearParams gpuBufferParams(gx, gy, gx, gy, gx, gy);

    // make a view of the parameters
    gpuViewParams = Opm::gpuistl::make_view<TraitsT, const GPUBuffer, GPUView>(gpuBufferParams);

    // make a twoPhaseLaw with the parameters
    auto gpuTwoPhaseLawParams = GPUTwoPhaseLawParams();
    gpuTwoPhaseLawParams.setGasOilParams(&gpuViewParams);
    gpuTwoPhaseLawParams.setOilWaterParams(&gpuViewParams);
    gpuTwoPhaseLawParams.setGasWaterParams(&gpuViewParams);
    gpuTwoPhaseLawParams.setApproach(Opm::EclTwoPhaseApproach::GasWater);
    gpuTwoPhaseLawParams.finalize();

    ValueVector testXs = {-1.0, 0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.1};

    // for (Scalar x_i : testXs){
    //     auto cpuMadeAd = SPE11CEvaluation(x_i, 0);
    //     SPE11CEvaluation cpuInterpolatedEval = CPUTPiecewiseLinearMaterialLaw::twoPhaseSatPcnw<SPE11CEvaluation>(cpuParams, cpuMadeAd);

    //     SPE11CEvaluation* gpuAdInput;
    //     SPE11CEvaluation* gpuAdResOnGPU;
    //     SPE11CEvaluation gpuAdResOnCPU[1];

    //     OPM_GPU_SAFE_CALL(cudaMalloc(&gpuAdInput, sizeof(SPE11CEvaluation)));
    //     OPM_GPU_SAFE_CALL(cudaMalloc(&gpuAdResOnGPU, sizeof(SPE11CEvaluation)));

    //     OPM_GPU_SAFE_CALL(cudaMemcpy(gpuAdInput, &cpuMadeAd, sizeof(SPE11CEvaluation), cudaMemcpyHostToDevice));
    //     gpuTwoPhaseSatPcnwWrapper<<<1,1>>>(gpuViewParams, *gpuAdInput, gpuAdResOnGPU);
    //     OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    //     OPM_GPU_SAFE_CALL(cudaMemcpy(gpuAdResOnCPU, gpuAdResOnGPU, sizeof(SPE11CEvaluation), cudaMemcpyDeviceToHost));

    //     BOOST_CHECK(gpuAdResOnCPU[0] == cpuInterpolatedEval);
    // }
}
