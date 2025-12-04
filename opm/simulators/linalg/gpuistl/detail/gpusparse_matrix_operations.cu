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
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/deviceBlockOperations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpuThreadUtils.hpp>
#include <stdexcept>

namespace Opm::gpuistl::detail
{
namespace
{
    template <class T, int blocksize>
    __global__ void cuMoveDataToReordered(
        const T* srcMatrix, const int* srcRowIndices, T* dstMatrix, int* dstRowIndices, int* indexConversion, size_t numberOfRows)
    {
        const auto srcRow = blockDim.x * blockIdx.x + threadIdx.x;
        if (srcRow < numberOfRows) {

            const auto dstRow = indexConversion[srcRow];

            for (int srcBlock = srcRowIndices[srcRow], dstBlock = dstRowIndices[dstRow];
                 srcBlock < srcRowIndices[srcRow + 1];
                 ++srcBlock, ++dstBlock) {
                for (int i = 0; i < blocksize; ++i) {
                    for (int j = 0; j < blocksize; ++j) {
                        dstMatrix[dstBlock * blocksize * blocksize + i * blocksize + j]
                            = srcMatrix[srcBlock * blocksize * blocksize + i * blocksize + j];
                    }
                }
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuMoveDataToReorderedSplit(const T* srcMatrix,
                                               const int* srcRowIndices,
                                               const int* srcColumnIndices,
                                               T* dstLowerMatrix,
                                               int* dstLowerRowIndices,
                                               T* dstUpperMatrix,
                                               int* dstUpperRowIndices,
                                               T* dstDiag,
                                               int* naturalToReordered,
                                               size_t numberOfRows)
    {
        const auto srcRow = blockDim.x * blockIdx.x + threadIdx.x;
        if (srcRow < numberOfRows) {

            const auto dstRow = naturalToReordered[srcRow];
            const auto rowStart = srcRowIndices[srcRow];
            const auto rowEnd = srcRowIndices[srcRow + 1];

            auto lowerBlock = dstLowerRowIndices[dstRow];
            auto upperBlock = dstUpperRowIndices[dstRow];

            for (int srcBlock = rowStart; srcBlock < rowEnd; srcBlock++) {
                int dstBlock;
                T* dstBuffer;

                if (srcColumnIndices[srcBlock] < srcRow) { // we are writing a value to the lower triangular matrix
                    dstBlock = lowerBlock;
                    ++lowerBlock;
                    dstBuffer = dstLowerMatrix;
                } else if (srcColumnIndices[srcBlock]
                           > srcRow) { // we are writing a value to the upper triangular matrix
                    dstBlock = upperBlock;
                    ++upperBlock;
                    dstBuffer = dstUpperMatrix;
                } else { // we are writing a value to the diagonal
                    dstBlock = dstRow;
                    dstBuffer = dstDiag;
                }
                for (int i = 0; i < blocksize; ++i) {
                    for (int j = 0; j < blocksize; ++j) {
                        dstBuffer[dstBlock * blocksize * blocksize + i * blocksize + j]
                            = srcMatrix[srcBlock * blocksize * blocksize + i * blocksize + j];
                    }
                }
            }
        }
    }

    template <class T>
    __global__ void cuComputeDiagPtrs(const int* rowIndices,
                                      const int* colIndices,
                                      int rows,
                                      T** diagPtrs,
                                      T* values,
                                      size_t blocksize)
    {
        const auto rawIdx = blockDim.x * blockIdx.x + threadIdx.x;
        if (rawIdx < rows) {
            const size_t nnzIdx = rowIndices[rawIdx];

            size_t diagIdx = nnzIdx;
            while (colIndices[diagIdx] != rawIdx) {
                ++diagIdx;
            }

            // Pointer arithmetic based on first element is sufficient here
            diagPtrs[rawIdx] = values + diagIdx*blocksize*blocksize;
        }
    }
} // namespace

template <class T, int blocksize>
void
copyMatDataToReordered(const T* srcMatrix,
                       const int* srcRowIndices,
                       T* dstMatrix,
                       int* dstRowIndices,
                       int* naturalToReordered,
                       size_t numberOfRows,
                       int thrBlockSize)
{
    int threadBlockSize
        = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuMoveDataToReordered<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(numberOfRows, threadBlockSize);
    cuMoveDataToReordered<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        srcMatrix, srcRowIndices, dstMatrix, dstRowIndices, naturalToReordered, numberOfRows);
}

template <class T, int blocksize>
void
copyMatDataToReorderedSplit(const T* srcMatrix,
                            const int* srcRowIndices,
                            const int* srcColumnIndices,
                            T* dstLowerMatrix,
                            int* dstLowerRowIndices,
                            T* dstUpperMatrix,
                            int* dstUpperRowIndices,
                            T* dstDiag,
                            int* naturalToReordered,
                            size_t numberOfRows,
                            int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
        cuMoveDataToReorderedSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(numberOfRows, threadBlockSize);
    cuMoveDataToReorderedSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(srcMatrix,
                                                                                 srcRowIndices,
                                                                                 srcColumnIndices,
                                                                                 dstLowerMatrix,
                                                                                 dstLowerRowIndices,
                                                                                 dstUpperMatrix,
                                                                                 dstUpperRowIndices,
                                                                                 dstDiag,
                                                                                 naturalToReordered,
                                                                                 numberOfRows);
}

template <class T>
GpuBuffer<T*>
getDiagPtrs(GpuSparseMatrixWrapper<T>& matrix)
{
    GpuBuffer<T*> diagPtrs(matrix.N());

    // Have cuda/hip automatically pick a reasonable block size for this kernel based on static analysis
    int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuComputeDiagPtrs<T>);
    // Calculate the number of blocks needed
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(matrix.N(), threadBlockSize);

    cuComputeDiagPtrs<T>
        <<<nThreadBlocks, threadBlockSize>>>(
            matrix->getRowIndices().data(), matrix->getColumnIndices().data(), matrix.N(), diagPtrs.data(), matrix->getNonZeroValues().data(), matrix.blockSize());
    return diagPtrs;
}

#define INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(T, blocksize)                                                                      \
    template void copyMatDataToReordered<T, blocksize>(const T*, const int*, T*, int*, int*, size_t, int);                         \
    template void copyMatDataToReorderedSplit<T, blocksize>(const T*, const int*, const int*, T*, int*, T*, int*, T*, int*, size_t, int);

INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(float, 1);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(float, 2);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(float, 3);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(float, 4);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(float, 5);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(float, 6);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(float, 7);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(double, 1);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(double, 2);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(double, 3);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(double, 4);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(double, 5);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(double, 6);
INSTANTIATE_KERNEL_WRAPPERS_WITH_SCALAR_AND_BLOCKSIZE(double, 7);

#define INSTANTIATE_KERNEL_WRAPPER_WITH_SCALAR(T) \
    template Opm::gpuistl::GpuBuffer<T*> Opm::gpuistl::detail::getDiagPtrs<T>(                           \
        Opm::gpuistl::GpuSparseMatrixWrapper<T>&);

INSTANTIATE_KERNEL_WRAPPER_WITH_SCALAR(float)
INSTANTIATE_KERNEL_WRAPPER_WITH_SCALAR(double)

} // namespace Opm::gpuistl::detail
