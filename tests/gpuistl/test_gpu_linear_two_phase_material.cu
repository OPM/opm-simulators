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

#define BOOST_TEST_MODULE TestGpuLinearTwoPhaseMaterial

#include <cuda_runtime.h>

#include <boost/test/unit_test.hpp>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>

#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <array>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <iostream>

  // these types are taken from Norne
  using Scalar = float;
  using ValueVector = std::vector<Scalar>;
  using GPUBuffer = Opm::gpuistl::GpuBuffer<Scalar>;
  using GPUView = Opm::gpuistl::GpuView<Scalar>;

  using TraitsT = Opm::TwoPhaseMaterialTraits<Scalar, 1, 2>;
  using CPUParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TraitsT>;
  using GPUBufferParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TraitsT, GPUBuffer>;
  using GPUViewParams = Opm::PiecewiseLinearTwoPhaseMaterialParams<TraitsT, GPUView>;

  using CPUTwoPhaseMaterialLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TraitsT, CPUParams>;
  using GPUTwoPhaseViewMaterialLaw = Opm::PiecewiseLinearTwoPhaseMaterial<TraitsT, GPUViewParams>;
  using NorneEvaluation = Opm::DenseAd::Evaluation<Scalar, 3, 0u>;

__global__ void gpuTwoPhaseSatPcnwWrapper(GPUTwoPhaseViewMaterialLaw::Params params, NorneEvaluation* Sw, NorneEvaluation* res){
    *res = GPUTwoPhaseViewMaterialLaw::twoPhaseSatPcnw(params, *Sw);
}

BOOST_AUTO_TEST_CASE(TestSimpleInterpolation)
{
    ValueVector cx = {0.0, 0.5, 1.0};
    ValueVector cy = {0.0, 0.9, 1.0};
    const GPUBuffer gx(cx);
    const GPUBuffer gy(cy);

    CPUParams cpuParams;
    cpuParams.setPcnwSamples(cx, cy);
    cpuParams.setKrwSamples(cx, cy);
    cpuParams.setKrnSamples(cx, cy);
    cpuParams.finalize();

    GPUBufferParams gpuBufferParams = Opm::gpuistl::copy_to_gpu<GPUBuffer>(cpuParams);

    GPUViewParams gpuViewParams = Opm::gpuistl::make_view<GPUView>(gpuBufferParams);

    ValueVector testXs = {-1.0, 0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.1};

    NorneEvaluation* gpuAdInput;
    NorneEvaluation* gpuAdResOnGPU;
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuAdInput, sizeof(NorneEvaluation)));
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuAdResOnGPU, sizeof(NorneEvaluation)));

    for (Scalar x_i : testXs){
        auto cpuMadeAd = NorneEvaluation(x_i, 0);
        NorneEvaluation cpuInterpolatedEval = CPUTwoPhaseMaterialLaw::twoPhaseSatPcnw<NorneEvaluation>(cpuParams, cpuMadeAd);

        NorneEvaluation gpuAdResOnCPU;

        OPM_GPU_SAFE_CALL(cudaMemcpy(gpuAdInput, &cpuMadeAd, sizeof(NorneEvaluation), cudaMemcpyHostToDevice));
        gpuTwoPhaseSatPcnwWrapper<<<1,1>>>(gpuViewParams, gpuAdInput, gpuAdResOnGPU);
        OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
        OPM_GPU_SAFE_CALL(cudaMemcpy(&gpuAdResOnCPU, gpuAdResOnGPU, sizeof(NorneEvaluation), cudaMemcpyDeviceToHost));

        BOOST_CHECK(gpuAdResOnCPU == cpuInterpolatedEval);
    }

    OPM_GPU_SAFE_CALL(cudaFree(gpuAdInput));
    OPM_GPU_SAFE_CALL(cudaFree(gpuAdResOnGPU));
}
