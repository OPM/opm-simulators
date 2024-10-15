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
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/DILUKernels.hpp>
#include <stdexcept>

namespace Opm::gpuistl::detail::DILU
{
namespace
{

    template <class T, int blocksize>
    __global__ void cuSolveLowerLevelSet(T* mat,
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

    template <int blocksize, class LinearSolverScalar, class MatrixScalar>
    __global__ void cuSolveLowerLevelSetSplit(MatrixScalar* mat,
                                              int* rowIndices,
                                              int* colIndices,
                                              int* indexConversion,
                                              int startIdx,
                                              int rowsInLevelSet,
                                              const MatrixScalar* dInv,
                                              const LinearSolverScalar* d,
                                              LinearSolverScalar* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {

            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            LinearSolverScalar rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            // TODO: removce the first condition in the for loop
            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                mmvMixedGeneral<blocksize, MatrixScalar, LinearSolverScalar, LinearSolverScalar, LinearSolverScalar>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mvMixedGeneral<blocksize, MatrixScalar, LinearSolverScalar, LinearSolverScalar, LinearSolverScalar>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }


    template <class T, int blocksize>
    __global__ void cuSolveUpperLevelSet(T* mat,
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

    template <int blocksize, class LinearSolverScalar, class MatrixScalar>
    __global__ void cuSolveUpperLevelSetSplit(MatrixScalar* mat,
                                              int* rowIndices,
                                              int* colIndices,
                                              int* indexConversion,
                                              int startIdx,
                                              int rowsInLevelSet,
                                              const MatrixScalar* dInv,
                                              LinearSolverScalar* v)
    {
        const auto reorderedRowIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const size_t nnzIdx = rowIndices[reorderedRowIdx];
            const size_t nnzIdxLim = rowIndices[reorderedRowIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedRowIdx];

            LinearSolverScalar rhs[blocksize] = {0};
            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                umvMixedGeneral<blocksize, MatrixScalar, LinearSolverScalar, LinearSolverScalar, LinearSolverScalar>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mmvMixedGeneral<blocksize, MatrixScalar, LinearSolverScalar, LinearSolverScalar, LinearSolverScalar>(&dInv[reorderedRowIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <class T, int blocksize>
    __global__ void cuComputeDiluDiagonal(T* mat,
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

    // TODO: rewrite such that during the factorization there is a dInv of InputScalar type that stores intermediate results
    // TOOD: The important part is to only cast after that is fully computed
    template <int blocksize, class InputScalar, class OutputScalar, bool copyResultToOtherMatrix>
    __global__ void cuComputeDiluDiagonalSplit(const InputScalar* srcReorderedLowerMat,
                                               int* lowerRowIndices,
                                               int* lowerColIndices,
                                               const InputScalar* srcReorderedUpperMat,
                                               int* upperRowIndices,
                                               int* upperColIndices,
                                               const InputScalar* srcDiagonal,
                                               int* reorderedToNatural,
                                               int* naturalToReordered,
                                               const int startIdx,
                                               int rowsInLevelSet,
                                               InputScalar* dInv,
                                               OutputScalar* dstDiag, // TODO: should this be diag or dInv?
                                               OutputScalar* dstLowerMat,
                                               OutputScalar* dstUpperMat)
    {
        const auto reorderedRowIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        if (reorderedRowIdx < rowsInLevelSet + startIdx) {
            const int naturalRowIdx = reorderedToNatural[reorderedRowIdx];
            const size_t lowerRowStart = lowerRowIndices[reorderedRowIdx];
            const size_t lowerRowEnd = lowerRowIndices[reorderedRowIdx + 1];

            InputScalar dInvTmp[blocksize * blocksize];
            for (int i = 0; i < blocksize; ++i) {
                for (int j = 0; j < blocksize; ++j) {
                    dInvTmp[i * blocksize + j] = srcDiagonal[reorderedRowIdx * blocksize * blocksize + i * blocksize + j];
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

                if constexpr (copyResultToOtherMatrix) {
                    // TODO: think long and hard about whether this performs only the wanted memory transfers
                    moveBlock<blocksize, InputScalar, OutputScalar>(&srcReorderedLowerMat[block * blocksize * blocksize], &dstLowerMat[block * blocksize * blocksize]);
                    moveBlock<blocksize, InputScalar, OutputScalar>(&srcReorderedUpperMat[symOppositeBlock * blocksize * blocksize], &dstUpperMat[symOppositeBlock * blocksize * blocksize]);
                }

                mmx2Subtraction<InputScalar, blocksize>(&srcReorderedLowerMat[block * blocksize * blocksize],
                                              &dInv[col * blocksize * blocksize],
                                              &srcReorderedUpperMat[symOppositeBlock * blocksize * blocksize],
                                              dInvTmp);
            }

            invBlockInPlace<InputScalar, blocksize>(dInvTmp);

            // for (int i = 0; i < blocksize; ++i) {
            //     for (int j = 0; j < blocksize; ++j) {
            //         dInv[reorderedRowIdx * blocksize * blocksize + i * blocksize + j] = dInvTmp[i * blocksize + j];
            //     }
            // }
            moveBlock<blocksize, InputScalar, InputScalar>(dInvTmp, &dInv[reorderedRowIdx * blocksize * blocksize]);
            if constexpr (copyResultToOtherMatrix) {
                moveBlock<blocksize, InputScalar, OutputScalar>(dInvTmp, &dstDiag[reorderedRowIdx * blocksize * blocksize]); // important!
            }
        }
    }
} // namespace

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
    int threadBlockSize
        = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveLowerLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveLowerLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, d, v);
}


template <int blocksize, class LinearSolverScalar, class MatrixScalar>
void
solveLowerLevelSetSplit(MatrixScalar* reorderedMat,
                        int* rowIndices,
                        int* colIndices,
                        int* indexConversion,
                        int startIdx,
                        int rowsInLevelSet,
                        const MatrixScalar* dInv,
                        const LinearSolverScalar* d,
                        LinearSolverScalar* v,
                        int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
        cuSolveLowerLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveLowerLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar><<<nThreadBlocks, threadBlockSize>>>(
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
    int threadBlockSize
        = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveUpperLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveUpperLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, v);
}

template <int blocksize, class LinearSolverScalar, class MatrixScalar>
void
solveUpperLevelSetSplit(MatrixScalar* reorderedMat,
                        int* rowIndices,
                        int* colIndices,
                        int* indexConversion,
                        int startIdx,
                        int rowsInLevelSet,
                        const MatrixScalar* dInv,
                        LinearSolverScalar* v,
                        int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
        cuSolveUpperLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveUpperLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar><<<nThreadBlocks, threadBlockSize>>>(
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
        int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
            cuComputeDiluDiagonal<T, blocksize>, thrBlockSize);
        int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
        cuComputeDiluDiagonal<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(reorderedMat,
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

template <int blocksize, class InputScalar, class OutputScalar, bool copyResultToOtherMatrix>
void
computeDiluDiagonalSplit(const InputScalar* srcReorderedLowerMat,
                         int* lowerRowIndices,
                         int* lowerColIndices,
                         const InputScalar* srcReorderedUpperMat,
                         int* upperRowIndices,
                         int* upperColIndices,
                         const InputScalar* srcDiagonal,
                         int* reorderedToNatural,
                         int* naturalToReordered,
                         const int startIdx,
                         int rowsInLevelSet,
                         InputScalar* dInv,
                         OutputScalar* dstDiag,
                         OutputScalar* dstLowerMat,
                         OutputScalar* dstUpperMat,
                         int thrBlockSize)
{
    if (blocksize <= 3) {
        int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
            cuComputeDiluDiagonalSplit<blocksize, InputScalar, OutputScalar, copyResultToOtherMatrix>, thrBlockSize);
        int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
        cuComputeDiluDiagonalSplit<blocksize, InputScalar, OutputScalar, copyResultToOtherMatrix><<<nThreadBlocks, threadBlockSize>>>(srcReorderedLowerMat,
                                                                                     lowerRowIndices,
                                                                                     lowerColIndices,
                                                                                     srcReorderedUpperMat,
                                                                                     upperRowIndices,
                                                                                     upperColIndices,
                                                                                     srcDiagonal,
                                                                                     reorderedToNatural,
                                                                                     naturalToReordered,
                                                                                     startIdx,
                                                                                     rowsInLevelSet,
                                                                                     dInv,
                                                                                     dstDiag,
                                                                                     dstLowerMat,
                                                                                     dstUpperMat);
    } else {
        OPM_THROW(std::invalid_argument, "Inverting diagonal is not implemented for blocksizes > 3");
    }
}

// TODO: format
#define INSTANTIATE_KERNEL_WRAPPERS(T, blocksize)                                                                      \
    template void computeDiluDiagonal<T, blocksize>(T*, int*, int*, int*, int*, const int, int, T*, int);              \
    template void computeDiluDiagonalSplit<blocksize, T, double, false>(                                                              \
        const T*, int*, int*, const T*, int*, int*, const T*, int*, int*, const int, int, T*, double*, double*, double*, int);                                      \
    template void computeDiluDiagonalSplit<blocksize, T, float, false>(                                                              \
        const T*, int*, int*, const T*, int*, int*, const T*, int*, int*, const int, int, T*, float*, float*, float*, int);                                      \
    template void computeDiluDiagonalSplit<blocksize, T, float, true>(                                                              \
        const T*, int*, int*, const T*, int*, int*, const T*, int*, int*, const int, int, T*, float*, float*, float*, int);                                      \
    template void computeDiluDiagonalSplit<blocksize, T, double, true>(                                                              \
        const T*, int*, int*, const T*, int*, int*, const T*, int*, int*, const int, int, T*, double*, double*, double*, int);                                      \
    template void solveUpperLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, T*, int);                 \
    template void solveLowerLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, const T*, T*, int);

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

#define INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(blocksize, LinearSolverScalar, MatrixScalar)                                 \
    template void solveUpperLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar>(                                \
        MatrixScalar*, int*, int*, int*, int, int, const MatrixScalar*, LinearSolverScalar*, int);                     \
    template void solveLowerLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar>(                                \
        MatrixScalar*, int*, int*, int*, int, int, const MatrixScalar*, const LinearSolverScalar*, LinearSolverScalar*, int);

INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(1, float, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(2, float, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(3, float, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(4, float, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(5, float, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(6, float, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(1, double, double);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(2, double, double);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(3, double, double);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(4, double, double);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(5, double, double);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(6, double, double);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(1, double, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(2, double, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(3, double, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(4, double, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(5, double, float);
INSTANTIATE_SOLVE_LEVEL_SET_SPLIT(6, double, float);

} // namespace Opm::gpuistl::detail::DILU
