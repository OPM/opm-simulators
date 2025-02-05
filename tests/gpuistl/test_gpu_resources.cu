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
#include <config.h>
#include <opm/simulators/linalg/gpuistl/gpu_resources.hpp>

#define BOOST_TEST_MODULE TestGPUResources

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>

namespace
{
__global__ void
dummyKernel()
{
}
} // namespace

BOOST_AUTO_TEST_CASE(TestStreams)
{
    using namespace Opm::gpuistl;
    BOOST_CHECK_NO_THROW(GPUStream());
    auto stream = GPUStream();
    BOOST_CHECK(stream.get() != nullptr);

    // Check that we can launch a kernel with said stream

    dummyKernel<<<1, 1, 0, stream.get()>>>();
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}

BOOST_AUTO_TEST_CASE(TestEvents)
{
    using namespace Opm::gpuistl;
    BOOST_CHECK_NO_THROW(GPUEvent());
    auto event = GPUEvent();
    BOOST_CHECK(event.get() != nullptr);

    // Check that we can record and wait for said event

    OPM_GPU_SAFE_CALL(cudaEventRecord(event.get()));
    OPM_GPU_SAFE_CALL(cudaEventSynchronize(event.get()));
    OPM_GPU_SAFE_CALL(cudaGetLastError());
}

BOOST_AUTO_TEST_CASE(TestGraph)
{
    using namespace Opm::gpuistl;
    BOOST_CHECK_NO_THROW(GPUGraph());
    auto graph = GPUGraph();
    BOOST_CHECK(graph.get() != nullptr);
}

BOOST_AUTO_TEST_CASE(TestGraphExec)
{
    using namespace Opm::gpuistl;
    BOOST_CHECK_NO_THROW(GPUGraphExec());
}
