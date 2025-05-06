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
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/deviceBlockOperations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpuThreadUtils.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/JacKernels.hpp>
#include <stdexcept>

namespace Opm::gpuistl::detail::JAC
{
namespace
{
    template <class T, int blocksize>
    __global__ void
    cuInvertDiagonalAndFlatten(const T* matNonZeroValues, const int* rowIndices, const int* colIndices, size_t numberOfRows, T* vec)
    {
        const auto row = blockDim.x * blockIdx.x + threadIdx.x;

        if (row < numberOfRows) {
            size_t nnzIdx = rowIndices[row];
            size_t nnzIdxLim = rowIndices[row + 1];

            // this loop will cause some extra checks that we are within the limit in the case of the diagonal having a
            // zero element
            while (colIndices[nnzIdx] != row && nnzIdx <= nnzIdxLim) {
                ++nnzIdx;
            }

            // diagBlock points to the start of where the diagonal block is stored
            const T* diagBlock = &matNonZeroValues[blocksize * blocksize * nnzIdx];
            // vecBlock points to the start of the block element in the vector where the inverse of the diagonal block
            // element should be stored
            T* vecBlock = &vec[blocksize * blocksize * row];

            invBlockOutOfPlace<T, blocksize>(diagBlock, vecBlock);
        }
    }
} // namespace

template <class T, int blocksize>
void
invertDiagonalAndFlatten(const T* mat, const int* rowIndices, const int* colIndices, size_t numberOfRows, T* vec)
{
    if (blocksize <= 3) {
        int threadBlockSize
            = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuInvertDiagonalAndFlatten<T, blocksize>);
        int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(numberOfRows, threadBlockSize);
        cuInvertDiagonalAndFlatten<T, blocksize>
            <<<nThreadBlocks, threadBlockSize>>>(mat, rowIndices, colIndices, numberOfRows, vec);
    } else {
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
    }
}

#define INSTANTIATE_KERNEL_WRAPPERS(T, blocksize)                                                                      \
    template void invertDiagonalAndFlatten<T, blocksize>(const T*, const int*, const int*, size_t, T*);

INSTANTIATE_KERNEL_WRAPPERS(float, 1);
INSTANTIATE_KERNEL_WRAPPERS(float, 2);
INSTANTIATE_KERNEL_WRAPPERS(float, 3);
INSTANTIATE_KERNEL_WRAPPERS(float, 4);
INSTANTIATE_KERNEL_WRAPPERS(float, 5);
INSTANTIATE_KERNEL_WRAPPERS(float, 6);
INSTANTIATE_KERNEL_WRAPPERS(double, 1);
INSTANTIATE_KERNEL_WRAPPERS(double, 2);
INSTANTIATE_KERNEL_WRAPPERS(double, 3);
INSTANTIATE_KERNEL_WRAPPERS(double, 4);
INSTANTIATE_KERNEL_WRAPPERS(double, 5);
INSTANTIATE_KERNEL_WRAPPERS(double, 6);

} // namespace Opm::gpuistl::detail::JAC
