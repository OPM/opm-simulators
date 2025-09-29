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

#include <opm/simulators/linalg/gpuistl/detail/cpr_amg_operations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/deviceBlockOperations.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpuThreadUtils.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <cuda_runtime.h>

namespace Opm::gpuistl::detail
{

namespace
{
    // Kernel for calculating quasi-IMPES weights
    template <typename T, bool transpose, int blockSize>
    __global__ void quasiImpesWeightsKernel(const T* matrix,
                                            T* weights,
                                            const int* diagonalIndices,
                                            const int numberOfRows,
                                            const int pressureVarIndex)
    {
        const auto row = blockDim.x * blockIdx.x + threadIdx.x;

        if (row < numberOfRows) {
            const int diagIdx = diagonalIndices[row];
            const int blockOffset = diagIdx * blockSize * blockSize;
            const T* block = matrix + blockOffset;

            // Set up RHS with 1.0 at pressure index
            T rhs[blockSize] = {0};
            rhs[pressureVarIndex] = 1.0;

            // Storage for solution
            T bweights[blockSize] = {0};

            // Solve the system
            if constexpr (transpose) {
                // Solve using original matrix
                solveBlock<T, blockSize>(block, rhs, bweights);
            } else {
                // Create transposed block for solving
                T transposed[blockSize * blockSize];
                transposeBlock<T, blockSize>(block, transposed);
                solveBlock<T, blockSize>(transposed, rhs, bweights);
            }

            // Find maximum absolute value for normalization
            T invMaxAbs = abs(bweights[0]);
            for (int j = 1; j < blockSize; ++j) {
                invMaxAbs = max(invMaxAbs, abs(bweights[j]));
            }
            invMaxAbs = T(1.0) / invMaxAbs;

            // Normalize and store weights
            for (int j = 0; j < blockSize; ++j) {
                weights[row * blockSize + j] = bweights[j] * invMaxAbs;
            }
        }
    }

    // Kernel to calculate matrix entries for the coarse level - processes each row in parallel
    template <typename T, bool transpose>
    __global__ void calculateCoarseEntriesKernel(const T* fineNonZeroValues,
                                                 T* coarseNonZeroValues,
                                                 const T* weights,
                                                 const int* rowIndices,
                                                 const int* colIndices,
                                                 const int numberOfRows,
                                                 const int blockSize,
                                                 const int pressureVarIndex)
    {
        // Each thread processes one row of the matrix
        const auto row = blockDim.x * blockIdx.x + threadIdx.x;

        if (row < numberOfRows) {
            // Get start and end indices for this row
            const int start = rowIndices[row];
            const int end = rowIndices[row + 1];

            // Process all non-zeros in this row
            for (int i = start; i < end; i++) {
                const int col = colIndices[i];
                const int blockOffset = i * blockSize * blockSize;
                T matrixEl = 0.0;

                if constexpr (transpose) {
                    // Use column weight
                    const T* bw = weights + col * blockSize;
                    for (int j = 0; j < blockSize; ++j) {
                        matrixEl += fineNonZeroValues[blockOffset + pressureVarIndex * blockSize + j] * bw[j];
                    }
                } else {
                    // Use row weight
                    const T* bw = weights + row * blockSize;
                    for (int j = 0; j < blockSize; ++j) {
                        matrixEl += fineNonZeroValues[blockOffset + j * blockSize + pressureVarIndex] * bw[j];
                    }
                }

                coarseNonZeroValues[i] = matrixEl;
            }
        }
    }

    // Kernel to restrict a fine vector to a coarse vector
    template <typename T, bool transpose>
    __global__ void restrictVectorKernel(const T* fine,
                                         T* coarse,
                                         const T* weights,
                                         const int numberOfBlocks,
                                         const int blockSize,
                                         const int pressureVarIndex)
    {
        const auto blockIndex = blockDim.x * blockIdx.x + threadIdx.x;

        if (blockIndex < numberOfBlocks) {
            T rhsEl = 0.0;

            if constexpr (transpose) {
                // Just extract the pressure component
                rhsEl = fine[blockIndex * blockSize + pressureVarIndex];
            } else {
                // Weighted sum of components
                const T* bw = weights + blockIndex * blockSize;
                for (int i = 0; i < blockSize; ++i) {
                    rhsEl += fine[blockIndex * blockSize + i] * bw[i];
                }
            }

            coarse[blockIndex] = rhsEl;
        }
    }

    // Kernel to prolongate a coarse vector to a fine vector
    template <typename T, bool transpose>
    __global__ void prolongateVectorKernel(const T* coarse,
                                           T* fine,
                                           const T* weights,
                                           const int numberOfBlocks,
                                           const int blockSize,
                                           const int pressureVarIndex)
    {
        const auto blockIndex = blockDim.x * blockIdx.x + threadIdx.x;

        if (blockIndex < numberOfBlocks) {
            if constexpr (transpose) {
                // Distribute the coarse value using weights
                const T* bw = weights + blockIndex * blockSize;
                for (int i = 0; i < blockSize; ++i) {
                    fine[blockIndex * blockSize + i] = coarse[blockIndex] * bw[i];
                }
            } else {
                // Only update the pressure component
                fine[blockIndex * blockSize + pressureVarIndex] = coarse[blockIndex];
            }
        }
    }

} // anonymous namespace

template <typename T, bool transpose, int blocksize>
void
dispatchQuasiImpesWeights(const GpuSparseMatrixWrapper<T>& matrix,
                          std::size_t pressureVarIndex,
                          GpuVector<T>& weights,
                          const GpuVector<int>& diagonalIndices,
                          int numberOfRows)
{
    if (matrix.blockSize() != blocksize) {
        if constexpr (blocksize > 1) {
            dispatchQuasiImpesWeights<T, transpose, blocksize - 1>(
                matrix, pressureVarIndex, weights, diagonalIndices, numberOfRows);
        } else {
            throw std::runtime_error("Unsupported block size for getQuasiImpesWeights: " + 
                                   std::to_string(matrix.blockSize()) + ". Only block sizes 1-3 are supported.");
        }
    } else {
        // Launch kernel with the correct block size
        int threadBlockSize = getCudaRecomendedThreadBlockSize(quasiImpesWeightsKernel<T, transpose, blocksize>);
        int nThreadBlocks = getNumberOfBlocks(numberOfRows, threadBlockSize);
        quasiImpesWeightsKernel<T, transpose, blocksize><<<nThreadBlocks, threadBlockSize>>>(
            matrix.getNonZeroValues().data(),
            weights.data(),
            diagonalIndices.data(),
            numberOfRows,
            pressureVarIndex);
    }
}

// Implementation of getQuasiImpesWeights for GPU
template <typename T, bool transpose>
void
getQuasiImpesWeights(const GpuSparseMatrixWrapper<T>& matrix,
                     std::size_t pressureVarIndex,
                     GpuVector<T>& weights,
                     const GpuVector<int>& diagonalIndices)
{
    const int blockSize = matrix.blockSize();
    const int numberOfRows = matrix.N();

    // Ensure weights vector has the right size
    if (weights.dim() != numberOfRows * blockSize) {
        throw std::runtime_error("Weights vector has incorrect size");
    }

    // Dispatch based on block size, max block size is 3
    dispatchQuasiImpesWeights<T, transpose, 3>(matrix, pressureVarIndex, weights, diagonalIndices, numberOfRows);
}

template <typename T, bool transpose>
void
calculateCoarseEntries(const GpuSparseMatrixWrapper<T>& fineMatrix,
                       GpuSparseMatrixWrapper<T>& coarseMatrix,
                       const GpuVector<T>& weights,
                       std::size_t pressureVarIndex)
{
    const int blockSize = fineMatrix.blockSize();
    const int numberOfRows = fineMatrix.N();

    int threadBlockSize = getCudaRecomendedThreadBlockSize(calculateCoarseEntriesKernel<T, transpose>);
    int nThreadBlocks = getNumberOfBlocks(numberOfRows, threadBlockSize);

    calculateCoarseEntriesKernel<T, transpose>
        <<<nThreadBlocks, threadBlockSize>>>(fineMatrix.getNonZeroValues().data(),
                                             coarseMatrix.getNonZeroValues().data(),
                                             weights.data(),
                                             fineMatrix.getRowIndices().data(),
                                             fineMatrix.getColumnIndices().data(),
                                             fineMatrix.N(),
                                             fineMatrix.blockSize(),
                                             pressureVarIndex);
}

template <typename T, bool transpose>
void
restrictVector(const GpuVector<T>& fine,
               GpuVector<T>& coarse,
               const GpuVector<T>& weights,
               std::size_t pressureVarIndex)
{
    const int blockSize = fine.dim() / coarse.dim();
    const int numberOfBlocks = coarse.dim();

    int threadBlockSize = getCudaRecomendedThreadBlockSize(restrictVectorKernel<T, transpose>);
    int nThreadBlocks = getNumberOfBlocks(numberOfBlocks, threadBlockSize);

    restrictVectorKernel<T, transpose><<<nThreadBlocks, threadBlockSize>>>(
        fine.data(), coarse.data(), weights.data(), numberOfBlocks, blockSize, pressureVarIndex);
}

template <typename T, bool transpose>
void
prolongateVector(const GpuVector<T>& coarse,
                 GpuVector<T>& fine,
                 const GpuVector<T>& weights,
                 std::size_t pressureVarIndex)
{
    const int blockSize = fine.dim() / coarse.dim();
    const int numberOfBlocks = coarse.dim();

    int threadBlockSize = getCudaRecomendedThreadBlockSize(prolongateVectorKernel<T, transpose>);
    int nThreadBlocks = getNumberOfBlocks(numberOfBlocks, threadBlockSize);

    prolongateVectorKernel<T, transpose><<<nThreadBlocks, threadBlockSize>>>(
        coarse.data(), fine.data(), weights.data(), numberOfBlocks, blockSize, pressureVarIndex);
}

#define INSTANTIATE_CPR_AMG_FUNCTIONS(ScalarType, TransposeMode)                                                       \
    template void getQuasiImpesWeights<ScalarType, TransposeMode>(const GpuSparseMatrixWrapper<ScalarType>& matrix,           \
                                                                  std::size_t pressureVarIndex,                        \
                                                                  GpuVector<ScalarType>& weights,                      \
                                                                  const GpuVector<int>& diagonalIndices);              \
    template void calculateCoarseEntries<ScalarType, TransposeMode>(const GpuSparseMatrixWrapper<ScalarType>& fineMatrix,     \
                                                                    GpuSparseMatrixWrapper<ScalarType>& coarseMatrix,         \
                                                                    const GpuVector<ScalarType>& weights,              \
                                                                    std::size_t pressureVarIndex);                     \
    template void restrictVector<ScalarType, TransposeMode>(const GpuVector<ScalarType>& fine,                         \
                                                            GpuVector<ScalarType>& coarse,                             \
                                                            const GpuVector<ScalarType>& weights,                      \
                                                            std::size_t pressureVarIndex);                             \
    template void prolongateVector<ScalarType, TransposeMode>(const GpuVector<ScalarType>& coarse,                     \
                                                              GpuVector<ScalarType>& fine,                             \
                                                              const GpuVector<ScalarType>& weights,                    \
                                                              std::size_t pressureVarIndex);

INSTANTIATE_CPR_AMG_FUNCTIONS(double, false)
INSTANTIATE_CPR_AMG_FUNCTIONS(double, true)
INSTANTIATE_CPR_AMG_FUNCTIONS(float, false)
INSTANTIATE_CPR_AMG_FUNCTIONS(float, true)

#undef INSTANTIATE_CPR_AMG_FUNCTIONS

} // namespace Opm::gpuistl::detail
