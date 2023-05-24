/*
  Copyright 2022-2023 SINTEF AS

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

#define BOOST_TEST_MODULE TestCudaSafeCall
#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>

BOOST_AUTO_TEST_CASE(TestCudaMalloc)
{
    void* pointer;
    BOOST_CHECK_NO_THROW(OPM_CUDA_SAFE_CALL(cudaMalloc(&pointer, 1)););
}


BOOST_AUTO_TEST_CASE(TestThrows)
{
    // Just testing a subset here.
    std::vector<cudaError_t> errorCodes {{cudaErrorAddressOfConstant, cudaErrorAlreadyAcquired}};
    for (auto code : errorCodes) {
        BOOST_CHECK_THROW(OPM_CUDA_SAFE_CALL(code), std::exception);
    }
}
