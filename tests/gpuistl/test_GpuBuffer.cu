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

#define BOOST_TEST_MODULE TestGpuBuffer

#include <boost/test/unit_test.hpp>
#include <cuda_runtime.h>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <type_traits>
#include <vector>

BOOST_AUTO_TEST_CASE(TestDocumentedUsage)
{
    auto someDataOnCPU = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});

    auto dataOnGPU = ::Opm::gpuistl::GpuBuffer<double>(someDataOnCPU);

    auto stdVectorOnCPU = dataOnGPU.asStdVector();

    BOOST_CHECK_EQUAL_COLLECTIONS(
        stdVectorOnCPU.begin(), stdVectorOnCPU.end(), someDataOnCPU.begin(), someDataOnCPU.end());
}

BOOST_AUTO_TEST_CASE(TestDefaultConstruction)
{
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>();
    BOOST_CHECK_EQUAL(0u, bufferOnGPU.size());
    BOOST_CHECK(bufferOnGPU.data() == nullptr);
}

BOOST_AUTO_TEST_CASE(TestConstructionSize)
{
    const size_t numberOfElements = 1234;
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(numberOfElements);
    BOOST_CHECK_EQUAL(numberOfElements, bufferOnGPU.size());
    BOOST_CHECK(bufferOnGPU.data() != nullptr);
}

BOOST_AUTO_TEST_CASE(TestCopyFromHostConstructor)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data.data(), data.size());
    BOOST_CHECK_EQUAL(data.size(), bufferOnGPU.size());
    std::vector<double> buffer(data.size(), 0.0);
    bufferOnGPU.copyToHost(buffer.data(), buffer.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(buffer.begin(), buffer.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(TestCopyFromHostConstructorWithGPUPointer)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data.data(), data.size());
    BOOST_CHECK_THROW(Opm::gpuistl::GpuBuffer<double>(bufferOnGPU.data(), data.size()), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestCopyFromHostFunction)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data.size());
    BOOST_CHECK_EQUAL(data.size(), bufferOnGPU.size());
    bufferOnGPU.copyFromHost(data.data(), data.size());
    std::vector<double> hostBuffer(data.size(), 0.0);
    bufferOnGPU.copyToHost(hostBuffer.data(), hostBuffer.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(hostBuffer.begin(), hostBuffer.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(TestCopyFromHostStdVector)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data.size());
    bufferOnGPU.copyFromHost(data);
    auto roundTrip = bufferOnGPU.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(roundTrip.begin(), roundTrip.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(TestCopyFromBvector)
{
    auto blockVector = Dune::BlockVector<Dune::FieldVector<double, 2>> {{{42, 43}, {44, 45}, {46, 47}}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(blockVector.dim());
    bufferOnGPU.copyFromHost(blockVector);
    std::vector<double> hostBuffer(blockVector.dim());
    bufferOnGPU.copyToHost(hostBuffer);

    BOOST_CHECK_EQUAL_COLLECTIONS(
        hostBuffer.begin(), hostBuffer.end(), &blockVector[0][0], &blockVector[0][0] + blockVector.dim());
}

BOOST_AUTO_TEST_CASE(TestCopyToBvector)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7, 8, 9}};
    auto blockVector = Dune::BlockVector<Dune::FieldVector<double, 3>>(3);
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data.data(), data.size());
    bufferOnGPU.copyToHost(blockVector);

    BOOST_CHECK_EQUAL_COLLECTIONS(data.begin(), data.end(), &blockVector[0][0], &blockVector[0][0] + blockVector.dim());
}

BOOST_AUTO_TEST_CASE(TestDataPointer)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7, 8, 9}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data.data(), data.size());

    std::vector<double> hostBuffer(data.size(), 0.0);
    OPM_GPU_SAFE_CALL(cudaMemcpy(hostBuffer.data(), bufferOnGPU.data(), sizeof(double) * data.size(), cudaMemcpyDeviceToHost));
    BOOST_CHECK_EQUAL_COLLECTIONS(data.begin(), data.end(), hostBuffer.begin(), hostBuffer.end());
}

BOOST_AUTO_TEST_CASE(TestCopyConstructor)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data);
    auto bufferCopy = Opm::gpuistl::GpuBuffer<double>(bufferOnGPU);

    BOOST_CHECK_EQUAL(bufferOnGPU.size(), bufferCopy.size());
    BOOST_CHECK(bufferOnGPU.data() != bufferCopy.data());

    auto original = bufferOnGPU.asStdVector();
    auto copied = bufferCopy.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(original.begin(), original.end(), copied.begin(), copied.end());
}

BOOST_AUTO_TEST_CASE(TestMoveConstructor)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data);
    auto* originalPointer = bufferOnGPU.data();

    auto moved = Opm::gpuistl::GpuBuffer<double>(std::move(bufferOnGPU));

    BOOST_CHECK_EQUAL(data.size(), moved.size());
    BOOST_CHECK_EQUAL(originalPointer, moved.data());
    BOOST_CHECK_EQUAL(0u, bufferOnGPU.size());
    BOOST_CHECK(bufferOnGPU.data() == nullptr);

    auto roundTrip = moved.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(roundTrip.begin(), roundTrip.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(TestMoveAssignment)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto source = Opm::gpuistl::GpuBuffer<double>(data);
    auto* originalPointer = source.data();

    auto destination = Opm::gpuistl::GpuBuffer<double>(3);
    destination = std::move(source);

    BOOST_CHECK_EQUAL(data.size(), destination.size());
    BOOST_CHECK_EQUAL(originalPointer, destination.data());
    BOOST_CHECK_EQUAL(0u, source.size());
    BOOST_CHECK(source.data() == nullptr);

    auto roundTrip = destination.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(roundTrip.begin(), roundTrip.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(TestResizeGrowPreservesData)
{
    std::vector<double> data {{1, 2, 3, 4}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data);

    bufferOnGPU.resize(7);
    BOOST_CHECK_EQUAL(7u, bufferOnGPU.size());

    // Prefix from before the resize must be preserved
    std::vector<double> hostBuffer(bufferOnGPU.size(), -1.0);
    bufferOnGPU.copyToHost(hostBuffer.data(), hostBuffer.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(hostBuffer.begin(), hostBuffer.begin() + data.size(), data.begin(), data.end());

    // Full new capacity must be usable as a real buffer of size 7
    std::vector<double> grownContent {{10, 20, 30, 40, 50, 60, 70}};
    bufferOnGPU.copyFromHost(grownContent);
    auto roundTrip = bufferOnGPU.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(roundTrip.begin(), roundTrip.end(), grownContent.begin(), grownContent.end());
}

BOOST_AUTO_TEST_CASE(TestResizeShrinkTruncatesData)
{
    std::vector<double> data {{1, 2, 3, 4, 5, 6, 7}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data);

    bufferOnGPU.resize(3);
    BOOST_CHECK_EQUAL(3u, bufferOnGPU.size());

    auto truncated = bufferOnGPU.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(truncated.begin(), truncated.end(), data.begin(), data.begin() + 3);
}

BOOST_AUTO_TEST_CASE(TestResizeFromEmpty)
{
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>();
    bufferOnGPU.resize(5);
    BOOST_CHECK_EQUAL(5u, bufferOnGPU.size());
    BOOST_CHECK(bufferOnGPU.data() != nullptr);
}

BOOST_AUTO_TEST_CASE(TestResizeSameSizeIsNoOp)
{
    std::vector<double> data {{1, 2, 3, 4}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data);
    auto* ptrBefore = bufferOnGPU.data();
    bufferOnGPU.resize(data.size());
    BOOST_CHECK_EQUAL(data.size(), bufferOnGPU.size());
    BOOST_CHECK_EQUAL(ptrBefore, bufferOnGPU.data());  // no realloc
    auto roundTrip = bufferOnGPU.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(roundTrip.begin(), roundTrip.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(TestResizeRejectsZero)
{
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(4);
    BOOST_CHECK_THROW(bufferOnGPU.resize(0), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestCopyTooManyElementsThrows)
{
    std::vector<double> data {{1, 2, 3, 4, 5}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(3);
    BOOST_CHECK_THROW(bufferOnGPU.copyFromHost(data.data(), data.size()), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestCopyToHostWrongSizeThrows)
{
    std::vector<double> data {{1, 2, 3, 4}};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<double>(data);
    std::vector<double> tooSmall(2);
    BOOST_CHECK_THROW(bufferOnGPU.copyToHost(tooSmall), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestBoolRoundTrip)
{
    std::vector<bool> data {true, false, true, true, false};
    auto bufferOnGPU = Opm::gpuistl::GpuBuffer<bool>(data);
    BOOST_CHECK_EQUAL(data.size(), bufferOnGPU.size());

    std::vector<bool> roundTrip(data.size());
    bufferOnGPU.copyToHost(roundTrip);
    BOOST_CHECK_EQUAL_COLLECTIONS(roundTrip.begin(), roundTrip.end(), data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(TestMakeView)
{
    // check creation of buffers and views for mutable buffers
    auto buf = std::vector<int>({1, 2, 3, 4, 5, 6});
    auto gpubuf = ::Opm::gpuistl::GpuBuffer<int>(buf);
    auto gpuview = ::Opm::gpuistl::make_view(gpubuf);
    bool gpuBufCreatedView = std::is_same_v<::Opm::gpuistl::GpuView<int>, decltype(gpuview)>;
    BOOST_CHECK(gpuBufCreatedView);

    auto gpubufOnCpu = gpubuf.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(gpubufOnCpu.begin(), gpubufOnCpu.end(), buf.begin(), buf.end());

    // check creation of buffers and views for const buffers
    const auto buf2 = std::vector<int>({2, 3, 4, 5, 6});
    const auto gpubuf2 = ::Opm::gpuistl::GpuBuffer<int>(buf2);
    auto gpuview2 = ::Opm::gpuistl::make_view(gpubuf2);
    bool gpuBufCreatedView2 = std::is_same_v<::Opm::gpuistl::GpuView<const int>, decltype(gpuview2)>;
    BOOST_CHECK(gpuBufCreatedView2);

    auto gpubufOnCpu2 = gpubuf2.asStdVector();
    BOOST_CHECK_EQUAL_COLLECTIONS(gpubufOnCpu2.begin(), gpubufOnCpu2.end(), buf2.begin(), buf2.end());
}
