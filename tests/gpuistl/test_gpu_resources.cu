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
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_constants.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_memcpy.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_resources.hpp>

#define BOOST_TEST_MODULE TestGPUResources

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>

#include <cstddef>
#include <initializer_list>
#include <vector>

namespace
{
__global__ void
dummyKernel()
{
}
__global__ void transform(int* values, std::size_t count, int factor, int offset)
{
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < count) {
        values[i] = factor * values[i] + offset;
    }
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
    GPUGraph graph;
    GPUGraphExec executable;
    // Give the owner a valid resource: destroying a null executable leaves a
    // runtime error that contaminates subsequent kernel error checks.
    OPM_GPU_SAFE_CALL(cudaGraphInstantiate(&executable.get(), graph.get(), nullptr, nullptr, 0));
    BOOST_CHECK(executable.get() != nullptr);
}

BOOST_AUTO_TEST_CASE(TestEventOrdersAsyncDeviceCopy)
{
    using namespace Opm::gpuistl;
    const std::vector<int> initial {3, -7, 11, 4, -2, 9, 5};
    GPUStream producer;
    GPUStream consumer;
    GPUEvent ready;
    GPUEvent complete;
    GpuBuffer<int> source(initial);
    GpuBuffer<int> destination(std::vector<int>(initial.size(), -99));

    // Exercise a separate consumer and the default-stream assembly boundary.
    for (const auto stream : {consumer.get(), detail::DEFAULT_STREAM}) {
        detail::gpuMemcpyHostToDeviceAsync(source.data(), initial.data(), initial.size(), producer.get());
        transform<<<1, 32, 0, producer.get()>>>(source.data(), initial.size(), 2, 1);
        OPM_GPU_SAFE_CALL(cudaGetLastError());
        recordEvent(ready, producer.get());
        waitForEventAsync(stream, ready);
        // Copy a subrange to check element counts, offsets and untouched guards.
        detail::gpuMemcpyDeviceToDeviceAsync(destination.data() + 1, source.data() + 2,
                                             initial.size() - 3, stream);
        transform<<<1, 32, 0, stream>>>(destination.data() + 1, initial.size() - 3, 3, -2);
        OPM_GPU_SAFE_CALL(cudaGetLastError());
        recordEvent(complete, stream);
        synchronizeEvent(complete);

        auto expected = std::vector<int>(initial.size(), -99);
        for (std::size_t i = 1; i < initial.size() - 2; ++i) {
            expected[i] = 3 * (2 * initial[i + 1] + 1) - 2;
        }
        const auto actual = destination.asStdVector();
        BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(), expected.end());
    }
    BOOST_CHECK_NO_THROW(detail::gpuMemcpyDeviceToDeviceAsync<int>(nullptr, nullptr, 0, producer.get()));
}

BOOST_AUTO_TEST_CASE(TestStreamDrainWithoutRecordedEvent)
{
    using namespace Opm::gpuistl;
    GPUStream stream;
    GpuBuffer<int> values(std::vector<int> {2, -4, 7});
    transform<<<1, 32, 0, stream.get()>>>(values.data(), values.size(), 5, 3);
    OPM_GPU_SAFE_CALL(cudaGetLastError());
    // The bridge must also drain work when an update exits before recording an event.
    synchronizeStreamNoThrow(stream.get());
    BOOST_CHECK_EQUAL(cudaStreamQuery(stream.get()), cudaSuccess);
    const auto actual = values.asStdVector();
    const std::vector<int> expected {13, -17, 38};
    BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(), expected.end());
}
