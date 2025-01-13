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
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/ILU0Kernels.hpp>
#include <stdexcept>

/*
    The LU factorization and apply step is written based on the Dune implementations
*/

namespace Opm::gpuistl::detail::ILU0
{
namespace
{
    template <class T, int blocksize>
    __global__ void cuLUFactorization(T* srcMatrix,
                                      int* srcRowIndices,
                                      int* srcColumnIndices,
                                      int* naturalToReordered,
                                      int* reorderedToNatual,
                                      size_t rowsInLevelSet,
                                      int startIdx)
    {
        auto reorderedIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        constexpr int scalarsInBlock = blocksize * blocksize;

        if (reorderedIdx < rowsInLevelSet + startIdx) {
            int naturalIdx = reorderedToNatual[reorderedIdx];
            // for each element under the diagonal
            int endOfRowI = srcRowIndices[reorderedIdx + 1];
            int ij;
            for (ij = srcRowIndices[reorderedIdx];
                 ij < srcRowIndices[reorderedIdx + 1] && srcColumnIndices[ij] < naturalIdx;
                 ++ij) {
                // consider the block we are looking at to be A_ij
                // we want to set B_ij = A_ij * A_jj, where A_jj is a diagonal element on a previous row, meaning it is
                // already inverted find the A_jj block
                int jj; // index in bcsr of the A_jj element
                int j = naturalToReordered[srcColumnIndices[ij]];
                int startOfRowJ = srcRowIndices[j];
                int endOfRowJ = srcRowIndices[j + 1];
                for (int blockInRowJ = startOfRowJ; blockInRowJ < endOfRowJ; ++blockInRowJ) {
                    if (srcColumnIndices[blockInRowJ] == srcColumnIndices[ij]) {
                        jj = blockInRowJ;
                        break;
                    }
                }

                // the DUNE implementation is inplace, and I need to take care to make sure
                // this out of place implementation is equivalent
                mmOverlap<T, blocksize>(
                    &srcMatrix[ij * scalarsInBlock], &srcMatrix[jj * scalarsInBlock], &srcMatrix[ij * scalarsInBlock]);

                // we have now accounted for an element under the diagonal and the diagonal element above it
                // now iterate over the blocks that align in the
                int jk = jj + 1;
                int ik = ij + 1;
                // this code is NOT gpu friendly, thread divergence will probably slow this down significantly...
                while (ik != endOfRowI && jk != endOfRowJ) {
                    if (srcColumnIndices[ik] == srcColumnIndices[jk]) {
                        // A_jk = A_ij * A_jk
                        T tmp[scalarsInBlock] = {0};
                        mmNoOverlap<T, blocksize>(
                            &srcMatrix[ij * scalarsInBlock], &srcMatrix[jk * scalarsInBlock], tmp);
                        matrixSubtraction<T, blocksize>(&srcMatrix[ik * scalarsInBlock], tmp);
                        ++ik;
                        ++jk;
                    } else {
                        if (srcColumnIndices[ik] < srcColumnIndices[jk]) {
                            ++ik;
                        } else {
                            ++jk;
                        }
                    }
                }
            }
            invBlockInPlace<T, blocksize>(&srcMatrix[ij * scalarsInBlock]);
        }
    }

    // An index in the BCSR matrix can be in one of three states, above, on, or above the diagonal
    enum class POSITION_TYPE { UNDER_DIAG, ON_DIAG, ABOVE_DIAG };
    // this handles incrementation of an index on a row that is in a lower/upper split datastructure.
    // The point is that if we are in the lower or upper part we increment as usual.
    // if we reach the diagonal, we update the state, the index is not needed.
    // when we increment from the diagonal, we set it to be the first block in the upper part of that row
    __device__ __forceinline__ void
    incrementAcrossSplitStructure(int& idx, POSITION_TYPE& currentState, int lastLowerIdx, int firstUpperIdx)
    {
        if (currentState == POSITION_TYPE::UNDER_DIAG) {
            if (idx != lastLowerIdx - 1) { // we should increment to be on the diagonal
                ++idx;
            } else {
                currentState = POSITION_TYPE::ON_DIAG;
            }
        } else if (currentState == POSITION_TYPE::ON_DIAG) {
            idx = firstUpperIdx;
            currentState = POSITION_TYPE::ABOVE_DIAG;
        } else { // currentState = POSITION_TYPE::ABOVE_DIAG
            ++idx;
        }
    }

    template <int blocksize, class InputScalar, class OutputScalar, bool copyResultToOtherMatrix>
    __global__ void cuLUFactorizationSplit(InputScalar* srcReorderedLowerMat,
                                           int* lowerRowIndices,
                                           int* lowerColIndices,
                                           InputScalar* srcReorderedUpperMat,
                                           int* upperRowIndices,
                                           int* upperColIndices,
                                           InputScalar* srcDiagonal,
                                           OutputScalar* dstReorderedLowerMat,
                                           OutputScalar* dstReorderedUpperMat,
                                           OutputScalar* dstDiagonal,
                                           int* reorderedToNatural,
                                           int* naturalToReordered,
                                           const int startIdx,
                                           int rowsInLevelSet)
    {
        auto reorderedIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        constexpr int scalarsInBlock = blocksize * blocksize;

        if (reorderedIdx < rowsInLevelSet + startIdx) {
            // if (reorderedIdx < rowsInLevelSet){
            int naturalIdx = reorderedToNatural[reorderedIdx];
            // for each element under the diagonal
            int endOfRowILower = lowerRowIndices[reorderedIdx + 1];
            int endOfRowIUpper = upperRowIndices[reorderedIdx + 1];
            int startOfRowILower = lowerRowIndices[reorderedIdx];
            int startOfRowIUpper = upperRowIndices[reorderedIdx];
            for (int ij = startOfRowILower; ij < endOfRowILower; ++ij) {

                // consider the block we are looking at to be A_ij
                // we want to set B_ij = A_ij * A_jj, where A_jj is a diagonal element on a previous row, meaning it is
                // already inverted find the A_jj block
                int j = naturalToReordered[lowerColIndices[ij]];
                int endOfRowJ = upperRowIndices[j + 1];
                int startOfRowJ = upperRowIndices[j];

                // the DUNE implementation is inplace, and I need to take care to make sure
                // this out of place implementation is equivalent
                mmOverlap<InputScalar, blocksize>(&srcReorderedLowerMat[ij * scalarsInBlock],
                                                  &srcDiagonal[j * scalarsInBlock],
                                                  &srcReorderedLowerMat[ij * scalarsInBlock]);
                if (copyResultToOtherMatrix) {
                    moveBlock<blocksize, InputScalar, OutputScalar>(&srcReorderedLowerMat[ij * scalarsInBlock],
                                                                    &dstReorderedLowerMat[ij * scalarsInBlock]);
                }

                // we have now accounted for an element under the diagonal and the diagonal element above it
                // now iterate over the blocks that align in the
                int jk = startOfRowJ;
                POSITION_TYPE ikState = POSITION_TYPE::UNDER_DIAG; // it is guaranteed that the ij element we consider
                                                                   // is under the diagonal
                int ik = ij;
                incrementAcrossSplitStructure(ik, ikState, endOfRowILower, startOfRowIUpper);
                // this code is NOT gpu friendly, thread divergence will probably slow this down significantly...

                // checking for ik != endOfRowUpper is not sufficient alone because
                // the same block index might be found in the lower part of the matrix
                while (!(ik == endOfRowIUpper && ikState == POSITION_TYPE::ABOVE_DIAG) && jk != endOfRowJ) {

                    InputScalar* ikBlockPtr;
                    OutputScalar* dstIkBlockPtr;
                    int ikColumn;
                    if (ikState == POSITION_TYPE::UNDER_DIAG) {
                        ikBlockPtr = &srcReorderedLowerMat[ik * scalarsInBlock];
                        if (copyResultToOtherMatrix)
                            dstIkBlockPtr = &dstReorderedLowerMat[ik * scalarsInBlock];
                        ikColumn = lowerColIndices[ik];
                    } else if (ikState == POSITION_TYPE::ON_DIAG) {
                        ikBlockPtr = &srcDiagonal[reorderedIdx * scalarsInBlock];
                        if (copyResultToOtherMatrix)
                            dstIkBlockPtr = &dstDiagonal[reorderedIdx * scalarsInBlock];
                        ikColumn = naturalIdx;
                    } else { // ikState == POSITION_TYPE::ABOVE_DIAG
                        ikBlockPtr = &srcReorderedUpperMat[ik * scalarsInBlock];
                        if (copyResultToOtherMatrix)
                            dstIkBlockPtr = &dstReorderedUpperMat[ik * scalarsInBlock];
                        ikColumn = upperColIndices[ik];
                    }

                    if (ikColumn == upperColIndices[jk]) {
                        // A_jk = A_ij * A_jk
                        InputScalar tmp[scalarsInBlock] = {0};

                        mmNoOverlap<InputScalar, blocksize>(&srcReorderedLowerMat[ij * scalarsInBlock],
                                                            &srcReorderedUpperMat[jk * scalarsInBlock],
                                                            tmp);
                        matrixSubtraction<InputScalar, blocksize>(ikBlockPtr, tmp);
                        incrementAcrossSplitStructure(ik, ikState, endOfRowILower, startOfRowIUpper);
                        ++jk;
                    } else {
                        if (ikColumn < upperColIndices[jk]) {
                            incrementAcrossSplitStructure(ik, ikState, endOfRowILower, startOfRowIUpper);
                        } else {
                            ++jk;
                        }
                    }
                }
            }
            invBlockInPlace<InputScalar, blocksize>(&srcDiagonal[reorderedIdx * scalarsInBlock]);
            if (copyResultToOtherMatrix) {
                moveBlock<blocksize, InputScalar, OutputScalar>(&srcDiagonal[reorderedIdx * scalarsInBlock],
                                                                &dstDiagonal[reorderedIdx * scalarsInBlock]);

                // also move all values above the diagonal on this row
                for (int block = startOfRowIUpper; block < endOfRowIUpper; ++block) {
                    moveBlock<blocksize, InputScalar, OutputScalar>(&srcReorderedUpperMat[block * scalarsInBlock],
                                                                    &dstReorderedUpperMat[block * scalarsInBlock]);
                }
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuSolveLowerLevelSet(T* mat,
                                         int* rowIndices,
                                         int* colIndices,
                                         int* indexConversion,
                                         int startIdx,
                                         int rowsInLevelSet,
                                         const T* d,
                                         T* v)
    {
        auto reorderedIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedIdx < rowsInLevelSet + startIdx) {

            const size_t nnzIdx = rowIndices[reorderedIdx];
            const int naturalRowIdx = indexConversion[reorderedIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            int block = nnzIdx;
            for (; colIndices[block] < naturalRowIdx; ++block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }
            for (int i = 0; i < blocksize; ++i) {
                v[naturalRowIdx * blocksize + i] = rhs[i];
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuSolveUpperLevelSet(
        T* mat, int* rowIndices, int* colIndices, int* indexConversion, int startIdx, int rowsInLevelSet, T* v)
    {
        auto reorderedIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedIdx < rowsInLevelSet + startIdx) {

            const size_t nnzIdxLim = rowIndices[reorderedIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = v[naturalRowIdx * blocksize + i];
            }

            int block = nnzIdxLim - 1;
            for (; colIndices[block] > naturalRowIdx; --block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mv<T, blocksize>(&mat[block * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }

    template <int blocksize, class LinearSolverScalar, class MatrixScalar>
    __global__ void cuSolveLowerLevelSetSplit(MatrixScalar* mat,
                                              int* rowIndices,
                                              int* colIndices,
                                              int* indexConversion,
                                              int startIdx,
                                              int rowsInLevelSet,
                                              const LinearSolverScalar* d,
                                              LinearSolverScalar* v)
    {
        auto reorderedIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedIdx < rowsInLevelSet + startIdx) {
            const size_t nnzIdx = rowIndices[reorderedIdx];
            const size_t nnzIdxLim = rowIndices[reorderedIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedIdx];

            LinearSolverScalar rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                mmvMixedGeneral<blocksize, MatrixScalar, LinearSolverScalar, LinearSolverScalar, LinearSolverScalar>(
                    &mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }
            for (int i = 0; i < blocksize; ++i) {
                v[naturalRowIdx * blocksize + i] = LinearSolverScalar(rhs[i]);
            }
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
        auto reorderedIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedIdx < rowsInLevelSet + startIdx) {
            const size_t nnzIdx = rowIndices[reorderedIdx];
            const size_t nnzIdxLim = rowIndices[reorderedIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedIdx];

            LinearSolverScalar rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = v[naturalRowIdx * blocksize + i];
            }

            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                mmvMixedGeneral<blocksize, MatrixScalar, LinearSolverScalar, LinearSolverScalar, LinearSolverScalar>(
                    &mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mvMixedGeneral<blocksize, MatrixScalar, LinearSolverScalar, LinearSolverScalar, LinearSolverScalar>(
                &dInv[reorderedIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
        }
    }
} // namespace

template <class T, int blocksize>
void
solveLowerLevelSet(T* reorderedMat,
                   int* rowIndices,
                   int* colIndices,
                   int* indexConversion,
                   int startIdx,
                   int rowsInLevelSet,
                   const T* d,
                   T* v,
                   int thrBlockSize)
{
    int threadBlockSize
        = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveLowerLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveLowerLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, d, v);
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
                   T* v,
                   int thrBlockSize)
{
    int threadBlockSize
        = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveUpperLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveUpperLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, v);
}

template <int blocksize, class LinearSolverScalar, class MatrixScalar>
void
solveLowerLevelSetSplit(MatrixScalar* reorderedMat,
                        int* rowIndices,
                        int* colIndices,
                        int* indexConversion,
                        int startIdx,
                        int rowsInLevelSet,
                        const LinearSolverScalar* d,
                        LinearSolverScalar* v,
                        int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
        cuSolveLowerLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveLowerLevelSetSplit<blocksize, LinearSolverScalar, MatrixScalar><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, d, v);
}
// perform the upper solve for all rows in the same level set
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
LUFactorization(T* srcMatrix,
                int* srcRowIndices,
                int* srcColumnIndices,
                int* naturalToReordered,
                int* reorderedToNatual,
                size_t rowsInLevelSet,
                int startIdx,
                int thrBlockSize)
{
    int threadBlockSize
        = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(cuLUFactorization<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuLUFactorization<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        srcMatrix, srcRowIndices, srcColumnIndices, naturalToReordered, reorderedToNatual, rowsInLevelSet, startIdx);
}

template <int blocksize, class InputScalar, class OutputScalar, bool copyResultToOtherMatrix>
void
LUFactorizationSplit(InputScalar* srcReorderedLowerMat,
                     int* lowerRowIndices,
                     int* lowerColIndices,
                     InputScalar* reorderedUpperMat,
                     int* upperRowIndices,
                     int* upperColIndices,
                     InputScalar* srcDiagonal,
                     OutputScalar* dstReorderedLowerMat,
                     OutputScalar* dstReorderedUpperMat,
                     OutputScalar* dstDiagonal,
                     int* reorderedToNatural,
                     int* naturalToReordered,
                     const int startIdx,
                     int rowsInLevelSet,
                     int thrBlockSize)
{
    int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
        cuLUFactorizationSplit<blocksize, InputScalar, OutputScalar, copyResultToOtherMatrix>, thrBlockSize);
    int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuLUFactorizationSplit<blocksize, InputScalar, OutputScalar, copyResultToOtherMatrix>
        <<<nThreadBlocks, threadBlockSize>>>(srcReorderedLowerMat,
                                             lowerRowIndices,
                                             lowerColIndices,
                                             reorderedUpperMat,
                                             upperRowIndices,
                                             upperColIndices,
                                             srcDiagonal,
                                             dstReorderedLowerMat,
                                             dstReorderedUpperMat,
                                             dstDiagonal,
                                             reorderedToNatural,
                                             naturalToReordered,
                                             startIdx,
                                             rowsInLevelSet);
}

#define INSTANTIATE_KERNEL_WRAPPERS(T, blocksize)                                                                      \
    template void solveUpperLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, T*, int);                           \
    template void solveLowerLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, T*, int);                 \
    template void LUFactorization<T, blocksize>(T*, int*, int*, int*, int*, size_t, int, int);                         \
    template void LUFactorizationSplit<blocksize, T, float, true>(                                                     \
        T*, int*, int*, T*, int*, int*, T*, float*, float*, float*, int*, int*, const int, int, int);                  \
    template void LUFactorizationSplit<blocksize, T, double, true>(                                                    \
        T*, int*, int*, T*, int*, int*, T*, double*, double*, double*, int*, int*, const int, int, int);               \
    template void LUFactorizationSplit<blocksize, T, float, false>(                                                    \
        T*, int*, int*, T*, int*, int*, T*, float*, float*, float*, int*, int*, const int, int, int);                  \
    template void LUFactorizationSplit<blocksize, T, double, false>(                                                   \
        T*, int*, int*, T*, int*, int*, T*, double*, double*, double*, int*, int*, const int, int, int);

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

#define INSTANTIATE_MIXED_PRECISION_KERNEL_WRAPPERS(blocksize)                                                         \
    /* double preconditioner */                                                                                        \
    template void solveLowerLevelSetSplit<blocksize, double, double>(                                                  \
        double*, int*, int*, int*, int, int, const double*, double*, int);                                             \
    /* float matrix, double compute preconditioner */                                                                  \
    template void solveLowerLevelSetSplit<blocksize, double, float>(                                                   \
        float*, int*, int*, int*, int, int, const double*, double*, int);                                              \
    /* float preconditioner */                                                                                         \
    template void solveLowerLevelSetSplit<blocksize, float, float>(                                                    \
        float*, int*, int*, int*, int, int, const float*, float*, int);                                                \
                                                                                                                       \
    /* double preconditioner */                                                                                        \
    template void solveUpperLevelSetSplit<blocksize, double, double>(                                                  \
        double*, int*, int*, int*, int, int, const double*, double*, int);                                             \
    /* float matrix, double compute preconditioner */                                                                  \
    template void solveUpperLevelSetSplit<blocksize, double, float>(                                                   \
        float*, int*, int*, int*, int, int, const float*, double*, int);                                               \
    /* float preconditioner */                                                                                         \
    template void solveUpperLevelSetSplit<blocksize, float, float>(                                                    \
        float*, int*, int*, int*, int, int, const float*, float*, int);

INSTANTIATE_MIXED_PRECISION_KERNEL_WRAPPERS(1);
INSTANTIATE_MIXED_PRECISION_KERNEL_WRAPPERS(2);
INSTANTIATE_MIXED_PRECISION_KERNEL_WRAPPERS(3);
INSTANTIATE_MIXED_PRECISION_KERNEL_WRAPPERS(4);
INSTANTIATE_MIXED_PRECISION_KERNEL_WRAPPERS(5);
INSTANTIATE_MIXED_PRECISION_KERNEL_WRAPPERS(6);
} // namespace Opm::gpuistl::detail::ILU0
