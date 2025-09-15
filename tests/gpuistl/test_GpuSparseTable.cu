/*
  Copyright 2025 SINTEF AS

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

#define BOOST_TEST_MODULE TestGpuSparseTable

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/istlsparsematrixadapter.hh>

#include <dune/istl/bcrsmatrix.hh>

#include <boost/mpl/range_c.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/SparseTable.hpp>

#include <cstddef>
#include <random>
#include <type_traits>


__global__ void fetchSparseTableValuesInKernel(Opm::SparseTable<int, Opm::gpuistl::GpuView> tableView,
                                               std::array<int, 7>* output)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx == 0) {
        int offset = 0;
        for (size_t row = 0; row < 3; ++row) {
            auto rowData = tableView[row];
            for (auto e : rowData) {
                (*output)[offset++] = e;
            }
        }
        (*output)[offset++] = tableView.dataSize();
    }
}

BOOST_AUTO_TEST_CASE(TestUsingSparseTableInKernel)
{
    /*
        Sparse table storting the data
            R0 (idx=0)  R1 (idx=3)  R2 (idx=4)
            {1,2,3},    {4},        {5,6}
    */
    Opm::SparseTable<int> cpuTable({1, 2, 3, 4, 5, 6}, {0, 3, 4, 6});
    for (size_t row = 0; row < 3; ++row) {
        if (row == 0) {
            BOOST_CHECK_EQUAL(cpuTable[row].size(), 3);
            BOOST_CHECK_EQUAL(cpuTable[row][0], 1);
            BOOST_CHECK_EQUAL(cpuTable[row][1], 2);
            BOOST_CHECK_EQUAL(cpuTable[row][2], 3);
        } else if (row == 1) {
            BOOST_CHECK_EQUAL(cpuTable[row].size(), 1);
            BOOST_CHECK_EQUAL(cpuTable[row][0], 4);
        } else if (row == 2) {
            BOOST_CHECK_EQUAL(cpuTable[row].size(), 2);
            BOOST_CHECK_EQUAL(cpuTable[row][0], 5);
            BOOST_CHECK_EQUAL(cpuTable[row][1], 6);
        }
    }

    auto gpuBufferTable = Opm::gpuistl::copy_to_gpu<int>(cpuTable);
    auto gpuTableView = Opm::gpuistl::make_view(gpuBufferTable);

    auto valuesFetchedFromGpu = Opm::gpuistl::make_gpu_unique_ptr<std::array<int, 7>>({-1});

    fetchSparseTableValuesInKernel<<<1, 1>>>(gpuTableView, valuesFetchedFromGpu.get());

    auto hostValues = Opm::gpuistl::copyFromGPU<std::array<int, 7>>(valuesFetchedFromGpu);

    BOOST_CHECK_EQUAL(hostValues[0], 1);
    BOOST_CHECK_EQUAL(hostValues[1], 2);
    BOOST_CHECK_EQUAL(hostValues[2], 3);
    BOOST_CHECK_EQUAL(hostValues[3], 4);
    BOOST_CHECK_EQUAL(hostValues[4], 5);
    BOOST_CHECK_EQUAL(hostValues[5], 6);
    BOOST_CHECK_EQUAL(hostValues[6], 6);
}
