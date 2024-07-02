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
#include <opm/simulators/linalg/cuistl/detail/gpuThreadUtils.hpp>
#include <opm/simulators/linalg/cuistl/detail/preconditionerKernels/ILU0Kernels.hpp>
#include <stdexcept>
#include <config.h>

namespace Opm::cuistl::detail::ILU0
{
namespace
{
    template <class T, int blocksize>
    __global__ void cuLUFactorization(T* srcMatrix, int* srcRowIndices, int* srcColumnIndices, int* naturalToReordered, int* reorderedToNatual, size_t rowsInLevelSet, int startIdx){
        auto reorderedIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        constexpr int scalarsInBlock = blocksize * blocksize;

        if (reorderedIdx < rowsInLevelSet + startIdx){
        // if (reorderedIdx < rowsInLevelSet){
            int naturalIdx = reorderedToNatual[reorderedIdx];
            // for each element under the diagonal
            int endOfRowI = srcRowIndices[reorderedIdx + 1];
            // printf("idx: %d->%d. In color [%d->%d]\n", naturalIdx, reorderedIdx, startIdx, startIdx + rowsInLevelSet);
            int ij;
            for (ij = srcRowIndices[reorderedIdx]; ij < srcRowIndices[reorderedIdx+1] && srcColumnIndices[ij] < naturalIdx; ++ij){
                // consider the block we are looking at to be A_ij
                // we want to set B_ij = A_ij * A_jj, where A_jj is a diagonal element on a previous row, meaning it is already inverted
                // find the A_jj block
                int jj; // index in bcsr of the A_jj element
                int j = naturalToReordered[srcColumnIndices[ij]];
                int startOfRowJ = srcRowIndices[j];
                int endOfRowJ = srcRowIndices[j + 1];
                for (int blockInRowJ = startOfRowJ; blockInRowJ < endOfRowJ; ++blockInRowJ){
                    if (srcColumnIndices[blockInRowJ] == srcColumnIndices[ij]){
                        jj = blockInRowJ;
                        break;
                    }
                }

                // the DUNE implementation is inplace, and I need to take care to make sure
                // this out of place implementation is equivalent
                mmOverlap<T, blocksize>(&srcMatrix[ij * scalarsInBlock], &srcMatrix[jj * scalarsInBlock], &srcMatrix[ij * scalarsInBlock]);

                // we have now accounted for an element under the diagonal and the diagonal element above it
                // now iterate over the blocks that align in the
                int jk = jj+1;
                int ik = ij+1;
                // this code is NOT gpu friendly, thread divergence will probably slow this down significantly...
                while (ik != endOfRowI && jk != endOfRowJ){
                    if (srcColumnIndices[ik] == srcColumnIndices[jk]){
                        // A_jk = A_ij * A_jk
                        T tmp[scalarsInBlock] = {0};
                        mmNoOverlap<T, blocksize>(&srcMatrix[ij * scalarsInBlock], &srcMatrix[jk * scalarsInBlock], tmp);
                        matrixSubtraction<T, blocksize>(&srcMatrix[ik * scalarsInBlock], tmp);
                        ++ik;
                        ++jk;
                    }
                    else{
                        if (srcColumnIndices[ik] < srcColumnIndices[jk]){
                            ++ik;
                        }
                        else{
                            ++jk;
                        }
                    }
                }
            }
            invBlockInPlace<T, blocksize>(&srcMatrix[ij * scalarsInBlock]);
        }
    }

    // An index in the BCSR matrix can be in one of three states, above, on, or above the diagonal
    enum class POSITION_TYPE {UNDER_DIAG, ON_DIAG, ABOVE_DIAG };
    // this handles incrementation of an index on a row that is in a lower/upper split datastructure.
    // The point is that if we are in the lower or upper part we increment as usual.
    // if we reach the diagonal, we update the state, the index is not needed.
    // when we increment from the diagonal, we set it to be the first block in the upper part of that row
    __device__ __forceinline__ void incrementAcrossSplitStructure(int &idx, POSITION_TYPE &currentState, int lastLowerIdx, int firstUpperIdx){
        if (currentState == POSITION_TYPE::UNDER_DIAG){
            if (idx != lastLowerIdx - 1){ // we should increment to be on the diagonal
                ++idx;
            }
            else{
                currentState = POSITION_TYPE::ON_DIAG;
            }
        }
        else if (currentState == POSITION_TYPE::ON_DIAG){
            idx = firstUpperIdx;
            currentState = POSITION_TYPE::ABOVE_DIAG;
        }
        else{ // currentState = POSITION_TYPE::ABOVE_DIAG
            ++idx;
        }
    }

    template <class T, int blocksize>
    __global__ void cuLUFactorizationSplit(T* reorderedLowerMat,
                                               int* lowerRowIndices,
                                               int* lowerColIndices,
                                               T* reorderedUpperMat,
                                               int* upperRowIndices,
                                               int* upperColIndices,
                                               T* diagonal,
                                               int* reorderedToNatural,
                                               int* naturalToReordered,
                                               const int startIdx,
                                               int rowsInLevelSet){
        auto reorderedIdx = startIdx + blockDim.x * blockIdx.x + threadIdx.x;
        constexpr int scalarsInBlock = blocksize * blocksize;

        if (reorderedIdx < rowsInLevelSet + startIdx){
        // if (reorderedIdx < rowsInLevelSet){
            int naturalIdx = reorderedToNatural[reorderedIdx];
            // for each element under the diagonal
            int endOfRowILower = lowerRowIndices[reorderedIdx + 1];
            int endOfRowIUpper = upperRowIndices[reorderedIdx + 1];
            int startOfRowILower = lowerRowIndices[reorderedIdx];
            int startOfRowIUpper = upperRowIndices[reorderedIdx];
            for (int ij = startOfRowILower; ij < endOfRowILower; ++ij){

                // consider the block we are looking at to be A_ij
                // we want to set B_ij = A_ij * A_jj, where A_jj is a diagonal element on a previous row, meaning it is already inverted
                // find the A_jj block
                int j = naturalToReordered[lowerColIndices[ij]];
                int endOfRowJ = upperRowIndices[j + 1];
                int startOfRowJ = upperRowIndices[j];

                // the DUNE implementation is inplace, and I need to take care to make sure
                // this out of place implementation is equivalent
                mmOverlap<T, blocksize>(&reorderedLowerMat[ij * scalarsInBlock], &diagonal[j * scalarsInBlock], &reorderedLowerMat[ij * scalarsInBlock]);

                // we have now accounted for an element under the diagonal and the diagonal element above it
                // now iterate over the blocks that align in the
                int jk = startOfRowJ;
                POSITION_TYPE ikState = POSITION_TYPE::UNDER_DIAG; // it is guaranteed that the ij element we consider is under the diagonal
                int ik = ij;
                incrementAcrossSplitStructure(ik, ikState, endOfRowILower, startOfRowIUpper);
                // this code is NOT gpu friendly, thread divergence will probably slow this down significantly...

                // checking for ik != endOfRowUpper is not sufficient alone because
                // the same block index might be found in the lower part of the matrix
                while (!(ik == endOfRowIUpper && ikState == POSITION_TYPE::ABOVE_DIAG) && jk != endOfRowJ){

                    T* ikBlockPtr;
                    int ikColumn;
                    if (ikState == POSITION_TYPE::UNDER_DIAG){
                        ikBlockPtr = &reorderedLowerMat[ik * scalarsInBlock];
                        ikColumn = lowerColIndices[ik];
                    }
                    else if(ikState == POSITION_TYPE::ON_DIAG){
                        ikBlockPtr = &diagonal[reorderedIdx*scalarsInBlock];
                        ikColumn = naturalIdx;
                    }
                    else{ // ikState == POSITION_TYPE::ABOVE_DIAG
                        ikBlockPtr = &reorderedUpperMat[ik * scalarsInBlock];
                        ikColumn = upperColIndices[ik];
                    }

                    // if (threadIdx.x == 0){
                    //     printf("ROW=%d, state: %d, col[ik]=%d, endOfRowILower: %d startOfRowIUpper: %d\n", reorderedIdx, ikState, ikColumn, endOfRowILower, startOfRowIUpper);
                    // }

                    if (ikColumn == upperColIndices[jk]){
                        // A_jk = A_ij * A_jk
                        T tmp[scalarsInBlock] = {0};

                        mmNoOverlap<T, blocksize>(&reorderedLowerMat[ij * scalarsInBlock], &reorderedUpperMat[jk * scalarsInBlock], tmp);
                        matrixSubtraction<T, blocksize>(ikBlockPtr, tmp);
                        incrementAcrossSplitStructure(ik, ikState, endOfRowILower, startOfRowIUpper);;
                        ++jk;
                    }
                    else{
                        if (ikColumn < upperColIndices[jk]){
                            incrementAcrossSplitStructure(ik, ikState, endOfRowILower, startOfRowIUpper);;
                        }
                        else{
                            ++jk;
                        }
                    }
                }
            }
            invBlockInPlace<T, blocksize>(&diagonal[reorderedIdx * scalarsInBlock]);
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
        if (reorderedIdx < rowsInLevelSet + startIdx){

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
            for (int i = 0; i < blocksize; ++i){
                v[naturalRowIdx*blocksize + i] = rhs[i];
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuSolveUpperLevelSet(T* mat,
                                                int* rowIndices,
                                                int* colIndices,
                                                int* indexConversion,
                                                int startIdx,
                                                int rowsInLevelSet,
                                                T* v)
    {
        auto reorderedIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedIdx < rowsInLevelSet + startIdx){

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

    template <class T, int blocksize>
    __global__ void cuSolveLowerLevelSetSplit(T* mat,
                                                     int* rowIndices,
                                                     int* colIndices,
                                                     int* indexConversion,
                                                     int startIdx,
                                                     int rowsInLevelSet,
                                                     const T* dInv,
                                                     const T* d,
                                                     T* v)
    {
        auto reorderedIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedIdx < rowsInLevelSet + startIdx){
            const size_t nnzIdx = rowIndices[reorderedIdx];
            const size_t nnzIdxLim = rowIndices[reorderedIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = d[naturalRowIdx * blocksize + i];
            }

            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }
            for (int i = 0; i < blocksize; ++i){
                v[naturalRowIdx*blocksize + i] = rhs[i];
            }
        }
    }

    template <class T, int blocksize>
    __global__ void cuSolveUpperLevelSetSplit(T* mat,
                                                     int* rowIndices,
                                                     int* colIndices,
                                                     int* indexConversion,
                                                     int startIdx,
                                                     int rowsInLevelSet,
                                                     const T* dInv,
                                                     T* v)
    {
        auto reorderedIdx = startIdx + (blockDim.x * blockIdx.x + threadIdx.x);
        if (reorderedIdx < rowsInLevelSet + startIdx){
            const size_t nnzIdx = rowIndices[reorderedIdx];
            const size_t nnzIdxLim = rowIndices[reorderedIdx + 1];
            const int naturalRowIdx = indexConversion[reorderedIdx];

            T rhs[blocksize];
            for (int i = 0; i < blocksize; i++) {
                rhs[i] = v[naturalRowIdx * blocksize + i];
            }

            for (int block = nnzIdx; block < nnzIdxLim; ++block) {
                const int col = colIndices[block];
                mmv<T, blocksize>(&mat[block * blocksize * blocksize], &v[col * blocksize], rhs);
            }

            mv<T, blocksize>(&dInv[reorderedIdx * blocksize * blocksize], rhs, &v[naturalRowIdx * blocksize]);
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
    int threadBlockSize = ::Opm::cuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveLowerLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::cuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
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
    int threadBlockSize = ::Opm::cuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveUpperLevelSet<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::cuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveUpperLevelSet<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, v);
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
    int threadBlockSize = ::Opm::cuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveLowerLevelSetSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::cuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveLowerLevelSetSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, d, v);
}
// perform the upper solve for all rows in the same level set
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
    int threadBlockSize = ::Opm::cuistl::detail::getCudaRecomendedThreadBlockSize(cuSolveUpperLevelSetSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::cuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuSolveUpperLevelSetSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(
        reorderedMat, rowIndices, colIndices, indexConversion, startIdx, rowsInLevelSet, dInv, v);
}

template <class T, int blocksize>
void LUFactorization(T* srcMatrix, int* srcRowIndices, int* srcColumnIndices, int* naturalToReordered, int* reorderedToNatual, size_t rowsInLevelSet, int startIdx, int thrBlockSize){
    int threadBlockSize = ::Opm::cuistl::detail::getCudaRecomendedThreadBlockSize(cuLUFactorization<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::cuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuLUFactorization<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(srcMatrix, srcRowIndices, srcColumnIndices, naturalToReordered, reorderedToNatual, rowsInLevelSet, startIdx);
}

template <class T, int blocksize>
void LUFactorizationSplit(T* reorderedLowerMat,
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
                         int thrBlockSize){
    int threadBlockSize = ::Opm::cuistl::detail::getCudaRecomendedThreadBlockSize(cuLUFactorizationSplit<T, blocksize>, thrBlockSize);
    int nThreadBlocks = ::Opm::cuistl::detail::getNumberOfBlocks(rowsInLevelSet, threadBlockSize);
    cuLUFactorizationSplit<T, blocksize><<<nThreadBlocks, threadBlockSize>>>(reorderedLowerMat,
                                                                                     lowerRowIndices,
                                                                                     lowerColIndices,
                                                                                     reorderedUpperMat,
                                                                                     upperRowIndices,
                                                                                     upperColIndices,
                                                                                     diagonal,
                                                                                     reorderedToNatural,
                                                                                     naturalToReordered,
                                                                                     startIdx,
                                                                                     rowsInLevelSet);
}

#define INSTANTIATE_KERNEL_WRAPPERS(T, blocksize)                                                                      \
    template void solveUpperLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, T*, int);             \
    template void solveLowerLevelSet<T, blocksize>(T*, int*, int*, int*, int, int, const T*, T*, int);     \
    template void LUFactorization<T, blocksize>(T*, int*, int*, int*, int*, size_t, int, int); \
    template void LUFactorizationSplit<T, blocksize>(T*, int*, int*, T*, int*, int*, T*, int*, int*, const int, int, int); \
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
} // namespace Opm::cuistl::detail::ILU0
