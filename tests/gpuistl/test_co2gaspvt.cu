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
  You shouf have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <config.h>

#define BOOST_TEST_MODULE TestGpuAD

#include <boost/test/unit_test.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>

#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>

#include <cuda_runtime.h>
#include <vector>
#include <utility>

using Evaluation = Opm::DenseAd::Evaluation<float, 3>;
using GpuB = const Opm::gpuistl::GpuBuffer<double>;
using GpuV = Opm::gpuistl::GpuView<const double>;
using GpuTab = Opm::UniformTabulated2DFunction<double, GpuV>;

namespace{
__global__ void gpuEvaluateUniformTabulated2DFunction(GpuTab gpuTab, Evaluation* inputX, Evaluation* inputY, double* res){
    *res = gpuTab.eval(*inputX, *inputY, true).value();
}

} // END EMPTY NAMESPACE


BOOST_AUTO_TEST_CASE(TestInstantiateCO2Object)
{

    std::vector<std::vector<double>> tabData = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};

    Opm::UniformTabulated2DFunction<double> cpuTab(1.0, 6.0, 3, 1.0, 6.0, 2, tabData);

    Opm::UniformTabulated2DFunction<double, GpuB> gpuBufTab = Opm::gpuistl::move_to_gpu<double, GpuB>(cpuTab);
    GpuTab gpuViewTab = Opm::gpuistl::make_view<double, GpuB, GpuV>(gpuBufTab);

    Evaluation a(2.3);
    Evaluation b(4.5);

    double* resultOnGpu;
    double resultOnCpu;
    OPM_GPU_SAFE_CALL(cudaMalloc(&resultOnGpu, sizeof(double)));

    Evaluation* gpuA;
    Evaluation* gpuB;

    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuA, sizeof(Evaluation)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuA, &a, sizeof(Evaluation), cudaMemcpyHostToDevice));
    OPM_GPU_SAFE_CALL(cudaMalloc(&gpuB, sizeof(Evaluation)));
    OPM_GPU_SAFE_CALL(cudaMemcpy(gpuB, &b, sizeof(Evaluation), cudaMemcpyHostToDevice));

    gpuEvaluateUniformTabulated2DFunction<<<1,1>>>(gpuViewTab, gpuA, gpuB, resultOnGpu);
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaMemcpy(&resultOnCpu, resultOnGpu, sizeof(double), cudaMemcpyDeviceToHost));
    
    OPM_GPU_SAFE_CALL(cudaFree(resultOnGpu));
    OPM_GPU_SAFE_CALL(cudaFree(gpuA));
    OPM_GPU_SAFE_CALL(cudaFree(gpuB));


    printf("%f == %f\n", cpuTab.eval(a, b, true).value(), resultOnCpu);

    BOOST_CHECK(resultOnCpu == cpuTab.eval(a, b, true).value());
}
