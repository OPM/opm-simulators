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
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <stdexcept>
namespace Opm::cuistl::detail
{

namespace
{   
    // TODO: combine the three following kernels into one using template arguments so that
    // no compile time optimization is lost while making the code more compact
    template <class T>
    __global__ void cuInvertDiagonalAndFlattenBlocksize1(
        T* matNonZeroValues, int* rowIndices, int* colIndices, size_t numberOfRows, T* vec)
    {
        const auto row = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 1;

        if (row < numberOfRows) {
            size_t nnzIdx = rowIndices[row];
            size_t nnzIdxLim = rowIndices[row + 1];

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a
            // zero element
            while (colIndices[nnzIdx] != row && nnzIdx <= nnzIdxLim) {
                ++nnzIdx;
            }

            // diagBlock points to the start of where the diagonal block is stored
            T* diagBlock = &matNonZeroValues[blocksize * blocksize * nnzIdx];
            // vecBlock points to the start of the block element in the vector where the inverse of the diagonal block
            // element should be stored
            T* vecBlock = &vec[blocksize * blocksize * row];

            vecBlock[0] = 1.0 / (diagBlock[0]);
        }
    }

    template <class T>
    __global__ void cuInvertDiagonalAndFlattenBlocksize2(
        T* matNonZeroValues, int* rowIndices, int* colIndices, size_t numberOfRows, T* vec)
    {
        const auto row = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 2;

        if (row < numberOfRows) {
            size_t nnzIdx = rowIndices[row];
            size_t nnzIdxLim = rowIndices[row + 1];

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a
            // zero element
            while (colIndices[nnzIdx] != row && nnzIdx <= nnzIdxLim) {
                ++nnzIdx;
            }

            // diagBlock points to the start of where the diagonal block is stored
            T* diagBlock = &matNonZeroValues[blocksize * blocksize * nnzIdx];
            // vecBlock points to the start of the block element in the vector where the inverse of the diagonal block
            // element should be stored
            T* vecBlock = &vec[blocksize * blocksize * row];

            // based on Dune implementation
            T detInv = 1.0 / (diagBlock[0] * diagBlock[3] - diagBlock[1] * diagBlock[2]);
            vecBlock[0] = diagBlock[3] * detInv;
            vecBlock[1] = -diagBlock[1] * detInv;
            vecBlock[2] = -diagBlock[2] * detInv;
            vecBlock[3] = diagBlock[0] * detInv;
        }
    }
    template <class T>
    __global__ void cuInvertDiagonalAndFlattenBlocksize3(
        T* matNonZeroValues, int* rowIndices, int* colIndices, size_t numberOfRows, T* vec)
    {
        const auto row = blockDim.x * blockIdx.x + threadIdx.x;
        const int blocksize = 3;

        if (row < numberOfRows) {
            size_t nnzIdx = rowIndices[row];
            size_t nnzIdxLim = rowIndices[row + 1];

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a
            // zero element
            while (colIndices[nnzIdx] != row && nnzIdx <= nnzIdxLim) {
                ++nnzIdx;
            }

            // diagBlock points to the start of where the diagonal block is stored
            T* diagBlock = &matNonZeroValues[blocksize * blocksize * nnzIdx];
            // vecBlock points to the start of the block element in the vector where the inverse of the diagonal block
            // element should be stored
            T* vecBlock = &vec[blocksize * blocksize * row];

            // based on Dune implementation
            T t4 = diagBlock[0] * diagBlock[4];
            T t6 = diagBlock[0] * diagBlock[5];
            T t8 = diagBlock[1] * diagBlock[3];
            T t10 = diagBlock[2] * diagBlock[3];
            T t12 = diagBlock[1] * diagBlock[6];
            T t14 = diagBlock[2] * diagBlock[6];

            T t17 = 1.0
                / (t4 * diagBlock[8] - t6 * diagBlock[7] - t8 * diagBlock[8] + t10 * diagBlock[7] + t12 * diagBlock[5]
                   - t14 * diagBlock[4]); // t17 is 1/determinant

            vecBlock[0] = (diagBlock[4] * diagBlock[8] - diagBlock[5] * diagBlock[7]) * t17;
            vecBlock[1] = -(diagBlock[1] * diagBlock[8] - diagBlock[2] * diagBlock[7]) * t17;
            vecBlock[2] = (diagBlock[1] * diagBlock[5] - diagBlock[2] * diagBlock[4]) * t17;
            vecBlock[3] = -(diagBlock[3] * diagBlock[8] - diagBlock[5] * diagBlock[6]) * t17;
            vecBlock[4] = (diagBlock[0] * diagBlock[8] - t14) * t17;
            vecBlock[5] = -(t6 - t10) * t17;
            vecBlock[6] = (diagBlock[3] * diagBlock[7] - diagBlock[4] * diagBlock[6]) * t17;
            vecBlock[7] = -(diagBlock[0] * diagBlock[7] - t12) * t17;
            vecBlock[8] = (t4 - t8) * t17;
        }
    }

    constexpr inline size_t getThreads([[maybe_unused]] size_t numberOfRows)
    {
        return 1024;
    }

    inline size_t getBlocks(size_t numberOfRows)
    {
        const auto threads = getThreads(numberOfRows);
        return (numberOfRows + threads - 1) / threads;
    }
} // namespace

template <class T>
void
invertDiagonalAndFlatten(T* mat, int* rowIndices, int* colIndices, size_t numberOfRows, size_t blocksize, T* vec)
{
    switch (blocksize) {
    case 1:
        cuInvertDiagonalAndFlattenBlocksize1<<<getBlocks(numberOfRows), getThreads(numberOfRows)>>>(
            mat, rowIndices, colIndices, numberOfRows, vec);
        break;
    case 2:
        cuInvertDiagonalAndFlattenBlocksize2<<<getBlocks(numberOfRows), getThreads(numberOfRows)>>>(
            mat, rowIndices, colIndices, numberOfRows, vec);
        break;
    case 3:
        cuInvertDiagonalAndFlattenBlocksize3<<<getBlocks(numberOfRows), getThreads(numberOfRows)>>>(
            mat, rowIndices, colIndices, numberOfRows, vec);
        break;
    default:
        // TODO: Figure out what is why it did not produce an error or any output in the output stream or the DBG file
        // when I forced this case to execute
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
        break;
    }
}

template void invertDiagonalAndFlatten(double*, int*, int*, size_t, size_t, double*);
template void invertDiagonalAndFlatten(float*, int*, int*, size_t, size_t, float*);

} // namespace Opm::cuistl::detail
