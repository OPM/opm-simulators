/*
  Copyright 2025 Equinor ASA

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
#include <boost/test/tools/old/interface.hpp>
#include <config.h>
#include <stdexcept>

#define BOOST_TEST_MODULE TestDenseVector

#include <cuda.h>
#include <cuda_runtime.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/dense/DenseVector.hpp>
#include <opm/simulators/linalg/gpuistl/dense/FieldVector.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

namespace {
  __global__ void writeToGlobalFieldVector(Opm::gpuistl::dense::FieldVector<double, 10> *vec) {
    (*vec)[2] = 1.0;
  }

  __global__ void makeFieldVector(double* value) {
    Opm::gpuistl::dense::FieldVector<double, 10> vec;
    vec[0] = 10;
    *value = vec[0];
  }
}

BOOST_AUTO_TEST_CASE(TestDenseCreation) 
{
    auto someVec = Opm::gpuistl::dense::FieldVector<double, 10>();

    BOOST_CHECK_EQUAL(someVec.size(), 10);
    someVec[9] = 1.0;
    BOOST_CHECK_EQUAL(someVec[9], 1.0);
}


BOOST_AUTO_TEST_CASE(TestDenseVector) 
{
    Opm::gpuistl::dense::FieldVector<double, 10> vec;
    auto gpuVec = Opm::gpuistl::make_gpu_shared_ptr(vec);
    writeToGlobalFieldVector<<<1, 1>>>(gpuVec.get());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    auto vecFromGPU = Opm::gpuistl::copyFromGPU(gpuVec);
    BOOST_CHECK_EQUAL(vecFromGPU[2], 1.0);

    auto valueOnGPU = Opm::gpuistl::make_gpu_shared_ptr(0.0);
    makeFieldVector<<<1, 1>>>(valueOnGPU.get());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    auto value = Opm::gpuistl::copyFromGPU(valueOnGPU);
    BOOST_CHECK_EQUAL(value, 10.0);
}