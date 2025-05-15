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
template <typename T>
__global__ void quasiImpesWeightsKernel(const T* matrix,
                                       T* weights,
                                       const int* diagonalIndices,
                                       const int numberOfRows,
                                       const int blockSize,
                                       const int pressureVarIndex,
                                       const bool transpose)
{
    const auto row = blockDim.x * blockIdx.x + threadIdx.x;

    if (row < numberOfRows) {
        const int diagIdx = diagonalIndices[row];
        const int blockOffset = diagIdx * blockSize * blockSize;
        const T* block = matrix + blockOffset;

        // Set up RHS with 1.0 at pressure index
        T rhs[3] = {0}; // Maximum block size is 3
        rhs[pressureVarIndex] = 1.0;

        // Storage for solution
        T bweights[3] = {0}; // Maximum block size is 3

        // Solve the system
        if (transpose) {
            // Solve using original matrix
            solveBlock(block, rhs, bweights, blockSize);
        } else {
            // Create transposed block for solving
            T transposed[9]; // Maximum block size squared is 3x3=9
            transposeBlock(block, transposed, blockSize);
            solveBlock(transposed, rhs, bweights, blockSize);
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
template <typename T>
__global__ void calculateCoarseEntriesKernel(const T* fineNonZeroValues,
                                           T* coarseNonZeroValues,
                                           const T* weights,
                                           const int* rowIndices,
                                           const int* colIndices,
                                           const int numberOfRows,
                                           const int blockSize,
                                           const int pressureVarIndex,
                                           const bool transpose)
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

            if (transpose) {
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
template <typename T>
__global__ void restrictVectorKernel(const T* fine,
                                    T* coarse,
                                    const T* weights,
                                    const int numberOfBlocks,
                                    const int blockSize,
                                    const int pressureVarIndex,
                                    const bool transpose)
{
    const auto blockIndex = blockDim.x * blockIdx.x + threadIdx.x;

    if (blockIndex < numberOfBlocks) {
        T rhsEl = 0.0;

        if (transpose) {
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
template <typename T>
__global__ void prolongateVectorKernel(const T* coarse,
                                      T* fine,
                                      const T* weights,
                                      const int numberOfBlocks,
                                      const int blockSize,
                                      const int pressureVarIndex,
                                      const bool transpose)
{
    const auto blockIndex = blockDim.x * blockIdx.x + threadIdx.x;

    if (blockIndex < numberOfBlocks) {
        if (transpose) {
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

// Implementation of getQuasiImpesWeights for GPU
template <typename T>
void getQuasiImpesWeights(const GpuSparseMatrix<T>& matrix,
                          std::size_t pressureVarIndex,
                          bool transpose,
                          GpuVector<T>& weights,
                          const GpuVector<int>& diagonalIndices)
{
    const int blockSize = matrix.blockSize();
    const int numberOfRows = matrix.N();

    // Check that block size is supported (1-3)
    if (blockSize < 1 || blockSize > 3) {
        throw std::runtime_error("Unsupported block size for getQuasiImpesWeights: " +
                                std::to_string(blockSize) +
                                ". Only block sizes 1-3 are supported.");
    }

    // Ensure weights vector has the right size
    if (weights.dim() != numberOfRows * blockSize) {
        throw std::runtime_error("Weights vector has incorrect size");
    }

    // Initialize weights to zero
    weights = 0.0;

    // Calculate optimal thread block size for CUDA
    int threadBlockSize = getCudaRecomendedThreadBlockSize(quasiImpesWeightsKernel<T>);
    int nThreadBlocks = getNumberOfBlocks(numberOfRows, threadBlockSize);

    // Launch kernel
    quasiImpesWeightsKernel<<<nThreadBlocks, threadBlockSize>>>(
        matrix.getNonZeroValues().data(),
        weights.data(),
        diagonalIndices.data(),
        numberOfRows,
        blockSize,
        pressureVarIndex,
        transpose
    );
}

template <typename T>
void calculateCoarseEntries(const GpuSparseMatrix<T>& fineMatrix,
                           GpuSparseMatrix<T>& coarseMatrix,
                           const GpuVector<T>& weights,
                           std::size_t pressureVarIndex,
                           bool transpose)
{
    const int blockSize = fineMatrix.blockSize();
    const int numberOfRows = fineMatrix.N();

    int threadBlockSize = getCudaRecomendedThreadBlockSize(calculateCoarseEntriesKernel<T>);
    int nThreadBlocks = getNumberOfBlocks(numberOfRows, threadBlockSize);

    calculateCoarseEntriesKernel<<<nThreadBlocks, threadBlockSize>>>(
        fineMatrix.getNonZeroValues().data(),
        coarseMatrix.getNonZeroValues().data(),
        weights.data(),
        fineMatrix.getRowIndices().data(),
        fineMatrix.getColumnIndices().data(),
        fineMatrix.N(),
        fineMatrix.blockSize(),
        pressureVarIndex,
        transpose
    );
}

template <typename T>
void restrictVector(const GpuVector<T>& fine,
                   GpuVector<T>& coarse,
                   const GpuVector<T>& weights,
                   std::size_t pressureVarIndex,
                   bool transpose)
{
    const int blockSize = fine.dim() / coarse.dim();
    const int numberOfBlocks = coarse.dim();

    int threadBlockSize = getCudaRecomendedThreadBlockSize(restrictVectorKernel<T>);
    int nThreadBlocks = getNumberOfBlocks(numberOfBlocks, threadBlockSize);

    restrictVectorKernel<<<nThreadBlocks, threadBlockSize>>>(
        fine.data(),
        coarse.data(),
        weights.data(),
        numberOfBlocks,
        blockSize,
        pressureVarIndex,
        transpose
    );
}

template <typename T>
void prolongateVector(const GpuVector<T>& coarse,
                     GpuVector<T>& fine,
                     const GpuVector<T>& weights,
                     std::size_t pressureVarIndex,
                     bool transpose)
{
    const int blockSize = fine.dim() / coarse.dim();
    const int numberOfBlocks = coarse.dim();

    int threadBlockSize = getCudaRecomendedThreadBlockSize(prolongateVectorKernel<T>);
    int nThreadBlocks = getNumberOfBlocks(numberOfBlocks, threadBlockSize);

    prolongateVectorKernel<<<nThreadBlocks, threadBlockSize>>>(
        coarse.data(),
        fine.data(),
        weights.data(),
        numberOfBlocks,
        blockSize,
        pressureVarIndex,
        transpose
    );
}

// Explicit template instantiations
template void getQuasiImpesWeights(
    const GpuSparseMatrix<double>& matrix,
    std::size_t pressureVarIndex,
    bool transpose,
    GpuVector<double>& weights,
    const GpuVector<int>& diagonalIndices);

template void getQuasiImpesWeights(
    const GpuSparseMatrix<float>& matrix,
    std::size_t pressureVarIndex,
    bool transpose,
    GpuVector<float>& weights,
    const GpuVector<int>& diagonalIndices);

template void calculateCoarseEntries(const GpuSparseMatrix<double>& fineMatrix,
    GpuSparseMatrix<double>& coarseMatrix,
    const GpuVector<double>& weights,
    std::size_t pressureVarIndex,
    bool transpose);
template void calculateCoarseEntries(
    const GpuSparseMatrix<float>& fineMatrix,
    GpuSparseMatrix<float>& coarseMatrix,
    const GpuVector<float>& weights,
    std::size_t pressureVarIndex,
    bool transpose);

template void restrictVector(
    const GpuVector<double>& fine,
    GpuVector<double>& coarse,
    const GpuVector<double>& weights,
    std::size_t pressureVarIndex,
    bool transpose);
template void restrictVector(
    const GpuVector<float>& fine,
    GpuVector<float>& coarse,
    const GpuVector<float>& weights,
    std::size_t pressureVarIndex,
    bool transpose);

template void prolongateVector(
    const GpuVector<double>& coarse,
    GpuVector<double>& fine,
    const GpuVector<double>& weights,
    std::size_t pressureVarIndex,
    bool transpose);
template void prolongateVector(
    const GpuVector<float>& coarse,
    GpuVector<float>& fine,
    const GpuVector<float>& weights,
    std::size_t pressureVarIndex,
    bool transpose);
} // namespace Opm::gpuistl::detail