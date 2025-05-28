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
#include <iostream>
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/ILU_variants_helper_kernels.hpp>

namespace Opm::gpuistl::detail
{
    namespace {
        // precompute the diagonal indices to speedup the update
        // The matrix here is not reordered, so the i'th thread computes
        // the diagonal element index not the i'th row of the matrix,
        // That is what we would have done if there was no coloring.
        // Instead it represents what the j'th rows in the k'th levelset
        // correspons to via the indexConversion array.
        template <class T>
        __global__ void cuComputeDiagIndicesNoReorder(const T* mat,
                                                const int* rowIndices,
                                                const int* colIndices,
                                                const size_t* indexConversion,
                                                int rows,
                                                size_t* diagIndices)
        {
            const auto rawIdx = blockDim.x * blockIdx.x + threadIdx.x;
            if (rawIdx < rows) {
                const size_t realRowIdx = indexConversion[rawIdx];
                const size_t nnzIdx = rowIndices[realRowIdx];

                size_t diagIdx = nnzIdx;
                while (colIndices[diagIdx] != realRowIdx) {
                    ++diagIdx;
                }

                diagIndices[realRowIdx] = diagIdx;
            }
        }

        // In this case we have a reordered matrix, so the firsth i rows will
        // in the same levelset and contiguous in the matrix. For this reason
        // we convert from the reordered index (i) to the natural index (j) which
        // will give us the diagonal element index as the row and column indices match there
        template <class T>
        __global__ void cuComputeDiagIndices(const T* mat,
                             const int* rowIndices,
                             const int* colIndices,
                             const int* reorderedToNatural,
                             int rows,
                             size_t* diagIndices)
        {
            const auto rawIdx = blockDim.x * blockIdx.x + threadIdx.x;
            if (rawIdx < rows) {
                const int naturalRowIdx = reorderedToNatural[rawIdx];
                const size_t nnzIdx = rowIndices[rawIdx];

                size_t diagIdx = nnzIdx;
                while (colIndices[diagIdx] != naturalRowIdx) {
                    ++diagIdx;
                }

                // I choose to have this indexed in the same order as the matrix in hopes
                // of memory coalescing.
                diagIndices[rawIdx] = diagIdx;
            }
        }

    } // empty namespace


    template <class T>
    void computeDiagIndicesNoReorder(const T* mat,
                                         const int* rowIndices,
                                         const int* colIndices,
                                         const size_t* indexConversion,
                                         int rows,
                                         size_t* diagIndices)
    {
        // Have cuda/hip automatically pick a reasonable block size for this kernel based on static analysis
        int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
            cuComputeDiagIndicesNoReorder<T>
        );

        // Calculate the number of blocks needed
        int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rows, threadBlockSize);

        // launch kernel
        cuComputeDiagIndicesNoReorder<T><<<threadBlockSize, nThreadBlocks>>>(
            mat, rowIndices, colIndices, indexConversion, rows, diagIndices
        );
    }

    template <class T>
    void computeDiagIndices(const T* mat,
                            const int* rowIndices,
                            const int* colIndices,
                            const int* reorderedToNatural,
                            int rows,
                            size_t* diagIndices)
    {
        // Have cuda/hip automatically pick a reasonable block size for this kernel based on static analysis
        int threadBlockSize = ::Opm::gpuistl::detail::getCudaRecomendedThreadBlockSize(
            cuComputeDiagIndices<T>
        );

        // Calculate the number of blocks needed
        int nThreadBlocks = ::Opm::gpuistl::detail::getNumberOfBlocks(rows, threadBlockSize);

        // launch kernel
        cuComputeDiagIndices<T><<<threadBlockSize, nThreadBlocks>>>(
            mat, rowIndices, colIndices, reorderedToNatural, rows, diagIndices
        );
    }


} // namespace Opm::gpuistl::detail

#define DEFINE_KERNEL_WRAPPERS_FOR_TYPES(T) \
    template void Opm::gpuistl::detail::computeDiagIndicesNoReorder<T>(const T*, \
                                                                              const int*, \
                                                                              const int*, \
                                                                              const size_t*, \
                                                                              int, \
                                                                              size_t*); \
    template void Opm::gpuistl::detail::computeDiagIndices<T>(const T*, \
                                                              const int*, \
                                                              const int*, \
                                                              const int*, \
                                                              int, \
                                                              size_t*);

DEFINE_KERNEL_WRAPPERS_FOR_TYPES(float)
DEFINE_KERNEL_WRAPPERS_FOR_TYPES(double)
