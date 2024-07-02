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
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/detail/deviceBlockOperations.hpp>
#include <opm/simulators/linalg/cuistl/detail/preconditionerKernels/DILUKernels.hpp>
#include <opm/simulators/linalg/cuistl/detail/gpuThreadUtils.hpp>
#include <stdexcept>
#include <config.h>

namespace Opm::gpuistl::detail::DILU
{
namespace
{

    template <class T, int blocksize>
    __global__ void gpuSolveLowerLevelSet(T* mat,
                                                int* rowIndices,
                                                int* colIndices,
                                                int* indexConversion,
                                                int startIdx,
                                                int rowsInLevelSet,
                                                const T* dInv,
                                                const T* d,
                                                T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {

            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            for (int block = nnzIdx; colIndices[block] < naturalRowIdx; ++block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <class T, int blocksize>
    __global__ void gpuSolveLowerLevelSetSplit(T* mat,
                                                     int* rowIndices,
                                                     int* colIndices,
                                                     int* indexConversion,
                                                     int startIdx,
                                                     int rowsInLevelSet,
                                                     const T* dInv,
                                                     const T* d,
                                                     T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {

            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            // TODO: removce the first condition in the for loop
            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }


    template <class T, int blocksize>
    __global__ void gpuSolveUpperLevelSet(T* mat,
                                                int* rowIndices,
                                                int* colIndices,
                                                int* indexConversion,
                                                int startIdx,
                                                int rowsInLevelSet,
                                                const T* dInv,
                                                T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize] = {0};
            for (int block = nnzIdxLim - 1; colIndices[block] > naturalRowIdx; --block) {
                const int col = colIndices[block];
                umv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mmv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <class T, int blocksize>
    __global__ void gpuSolveUpperLevelSetSplit(T* mat,
                                                     int* rowIndices,
                                                     int* colIndices,
                                                     int* indexConversion,
                                                     int startIdx,
                                                     int rowsInLevelSet,
                                                     const T* dInv,
                                                     T* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            T rhs[blocksize] = {0};
            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                umv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mmv<T, blocksize>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <class T, int blocksize>
    __global__ void gpuComputeDiluDiagonal(T* mat,
                                          int* rowIndices,
                                          int* colIndices,
                                          int* reorderedToNatural,
                                          int* naturalToReordered,
                                          const int startIdx,
                                          int rowsInLevelSet,
                                          T* dInv)
    {
        const auto reorderedRowIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const int naturalRowIdx = reorderedToNatural[reorderedRowIdx];
            const size_t nnzIdx = rowIndices[reorderedRowIdx];

            int diagIdx = nnzIdx;
            while (colIndices[diagIdx] != naturalRowIdx) {
                ++diagIdx;
            }

            T dInvTmp[blocksize * blocksize];
            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInvTmp[i * blocksize + j] = mat[diagIdx * blocksize * blocksize + i * blocksize + j];
                }
            }

            for (int block = nnzIdx; colIndices[block] < naturalRowIdx; ++block) {
                const int col = naturalToReordered[colIndices[block]];
                // find element with indices swapped
                // Binary search over block in the right row, [rowIndices[col], rowindices[col+1]-1] defines the range
                // we binary search over
                int left = rowIndices[col];
                int right = rowIndices[col + 1] - 1;
                int mid;

                while (left <= right) {
                    mid = left + (right - left) / 2; // overflow-safe average
                    const int col = colIndices[mid];

                    if (col == naturalRowIdx) {
                        break;
                    } else if (col < naturalRowIdx) {
                        left = mid + 1;
                    } else {
                        right = mid - 1;
                    }
                }

                const int symOpposite = mid;

                mmx2Subtraction<T, blocksize>(&mat[block * blocksize * blocksize],
                                              &dInv[col * blocksize * blocksize],
                                              &mat[symOpposite * blocksize * blocksize],
                                              dInvTmp);
            }

            invBlockInPlace<T, blocksize>(dInvTmp);

            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInv[reorderedRowIdx * blocksize * blocksize + i * blocksize + j] = dInvTmp[i * blocksize + j];
                }
            }
        }
    }

    template <class T, int blocksize>
    __global__ void gpuComputeDiluDiagonalSplit(T* reorderedLowerMat,
                                               int* lowerRowIndices,
                                               int* lowerColIndices,
                                               T* reorderedUpperMat,
                                               int* upperRowIndices,
                                               int* upperColIndices,
                                               T* diagonal,
                                               int* reorderedToNatural,
                                               int* naturalToReordered,
                                               const int startIdx,
                                               int rowsInLevelSet,
                                               T* dInv)
    {
        const auto reorderedRowIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const int naturalRowIdx = reorderedToNatural[reorderedRowIdx];
            const size_t lowerRowStart = lowerRowIndices[reorderedRowIdx];
            const size_t lowerRowEnd = lowerRowIndices[reorderedRowIdx + 1];

            T dInvTmp[blocksize * blocksize];
            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInvTmp[i * blocksize + j] = diagonal[reorderedRowIdx * blocksize * blocksize + i * blocksize + j];
                }
            }

            for (int block = lowerRowStart; block < lowerRowEnd; ++block) {
                const int col = naturalToReordered[lowerColIndices[block]];

                int symOppositeIdx = upperRowIndices[col];
                for (; symOppositeIdx < upperRowIndices[col + 1]; ++symOppositeIdx) {
                    if (naturalRowIdx == upperColIndices[symOppositeIdx]) {
                        break;
                    }
                }

                const int symOppositeBlock = symOppositeIdx;

                mmx2Subtraction<T, blocksize>(&reorderedLowerMat[block * blocksize * blocksize],
                                              &dInv[col * blocksize * blocksize],
                                              &reorderedUpperMat[symOppositeBlock * blocksize * blocksize],
                                              dInvTmp);
            }

            invBlockInPlace<T, blocksize>(dInvTmp);

            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInv[reorderedRowIdx * blocksize * blocksize + i * blocksize + j] = dInvTmp[i * blocksize + j];
                }
            }
        }
    }
} // EMPTY namespace

// perform the lower solve for all rows in the same level set
template <class T, int blocksize>
void
solveLowerLevelSet(T* reorderedMat,
                          int* rowIndices,
                          int* colIndices,
                          int* indexConversion,
                          int startIdx,
                          int rowsInLevelSet,
                          const T* dInv,
                          const T* d,
                          T* v,
                          int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getRecomendedThreadBlockSize(gpuSolveLowerLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    gpuSolveLowerLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, d, v);
}


template <class T, int blocksize>
void
solveLowerLevelSetSplit(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               const T* d,
                               T* v,
                               int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getRecomendedThreadBlockSize(gpuSolveLowerLevelSetSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    gpuSolveLowerLevelSetSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, d, v);
}
// perform the upper solve for all rows in the same level set
template <class T, int blocksize>
void
solveUpperLevelSet(T* reorderedMat,
                          int* rowIndices,
                          int* colIndices,
                          int* indexConversion,
                          int startIdx,
                          int rowsInLevelSet,
                          const T* dInv,
                          T* v,
                          int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getRecomendedThreadBlockSize(gpuSolveUpperLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    gpuSolveUpperLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, v);
}

template <class T, int blocksize>
void
solveUpperLevelSetSplit(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               T* v,
                               int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getRecomendedThreadBlockSize(gpuSolveUpperLevelSetSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    gpuSolveUpperLevelSetSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, v);
}

template <class T, int blocksize>
void
computeDiluDiagonal(T* reorderedMat,
                    int* rowIndices,
                    int* colIndices,
                    int* reorderedToNatural,
                    int* naturalToReordered,
                    const int startIdx,
                    int rowsInLevelSet,
                    T* dInv,
                    int thrBlockSize)
{
    if (blocksize <= 3) {
        int threadBlockSize = ::Opm::gpuistl::detail::getRecomendedThreadBlockSize(gpuComputeDiluDiagonal<T, blocksize>, thrBlockSize);
        int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
        gpuComputeDiluDiagonal<T, blocksize>
            <<<nThreadBlocks, threadBlockSize>>>(reorderedMat,
                                                                        rowIndices,
                                                                        colIndices,
                                                                        reorderedToNatural,
                                                                        naturalToReordered,
                                                                        startIdx,
                                                                        rowsInLevelSet,
                                                                        dInv);
    } else {
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
    }
}

template <class T, int blocksize>
void
computeDiluDiagonalSplit(T* reorderedLowerMat,
                         int* lowerRowIndices,
                         int* lowerColIndices,
                         T* reorderedUpperMat,
                         int* upperRowIndices,
                         int* upperColIndices,
                         T* diagonal,
                         int* reorderedToNatural,
                         int* naturalToReordered,
                         const int startIdx,
                         int rowsInLevelSet,
                         T* dInv,
                         int thrBlockSize)
{
    if (blocksize <= 3) {
        int threadBlockSize = ::Opm::gpuistl::detail::getRecomendedThreadBlockSize(gpuComputeDiluDiagonalSplit<T, blocksize>, thrBlockSize);
        int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
        gpuComputeDiluDiagonalSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(reorderedLowerMat,
                                                                                     lowerRowIndices,
                                                                                     lowerColIndices,
                                                                                     reorderedUpperMat,
                                                                                     upperRowIndices,
                                                                                     upperColIndices,
                                                                                     diagonal,
                                                                                     reorderedToNatural,
                                                                                     naturalToReordered,
                                                                                     startIdx,
                                                                                     rowsInLevelSet,
                                                                                     dInv);
    } else {
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
    }
}

#define INSTANTIATE_KERNEL_WRAPPERS(T, blocksize)                                                                      \
    template void computeDiluDiagonal<T, blocksize>(T*, int*, int*, int*, int*, const int, int, T*, int);                   \
    template void computeDiluDiagonalSplit<T, blocksize>(                                                              \
        T*, int*, int*, T*, int*, int*, T*, int*, int*, const int, int, T*, int);                                           \
    template void solveUpperLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, T*, int);               \
    template void solveLowerLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, const T*, T*, int);     \
    template void solveUpperLevelSetSplit<T, blocksize>(T*, int*, int*, int*, int, int, const T*, T*, int);          \
    template void solveLowerLevelSetSplit<T, blocksize>(T*, int*, int*, int*, int, int, const T*, const T*, T*, int);

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
} // namespace Opm::gpuistl::detail::DILU
