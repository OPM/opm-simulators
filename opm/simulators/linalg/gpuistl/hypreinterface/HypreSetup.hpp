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

#ifndef OPM_HYPRE_SETUP_GPU_HPP
#define OPM_HYPRE_SETUP_GPU_HPP

#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreDataStructures.hpp>
#include <opm/simulators/linalg/hypreinterface/HypreErrorHandling.hpp>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/graph.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/repartition.hh>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrixWrapper.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#endif
#endif // HAVE_CUDA

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>

#include <algorithm>
#include <cstddef>

namespace Opm::gpuistl::HypreInterface
{

// GPU-specific helper functions
template <typename T, bool ForceLegacy>
linalg::HypreInterface::SparsityPattern
setupSparsityPatternFromGpuMatrix(const GpuSparseMatrixWrapper<T, ForceLegacy>& gpu_matrix,
                                  const linalg::HypreInterface::ParallelInfo& par_info,
                                  bool owner_first);

template <typename T, bool ForceLegacy>
std::vector<HYPRE_Int> computeRowIndexesWithMappingGpu(const GpuSparseMatrixWrapper<T, ForceLegacy>& gpu_matrix,
                                                       const std::vector<int>& local_dune_to_local_hypre);

/**
 * @brief Setup sparsity pattern from GPU matrix (GpuSparseMatrix)
 *
 * @param gpu_matrix The GPU matrix to analyze
 * @param par_info Parallel information structure
 * @param owner_first Whether all owned DOFs come first in the ordering
 * @return Sparsity pattern information
 */
template <typename T, bool ForceLegacy>
linalg::HypreInterface::SparsityPattern
setupSparsityPatternFromGpuMatrix(const GpuSparseMatrixWrapper<T, ForceLegacy>& gpu_matrix,
                                  const linalg::HypreInterface::ParallelInfo& par_info,
                                  bool owner_first)
{
    linalg::HypreInterface::SparsityPattern pattern;

    // Determine the size for cols array based on owner_first
    if (owner_first) {
        std::size_t cols_size = 0;
        // For owner_first=true case, we need to calculate how many owned entries there are
        auto host_row_ptrs = gpu_matrix.getRowIndices().asStdVector();
        for (int rind = 0; rind < static_cast<int>(gpu_matrix.N()); ++rind) {
            if (par_info.local_dune_to_local_hypre[rind] >= 0) {
                const int row_start = host_row_ptrs[rind];
                const int row_end = host_row_ptrs[rind + 1];
                cols_size += (row_end - row_start);
            }
        }
        pattern.nnz = cols_size;
    } else {
        // Full matrix space case: all entries (including ghost rows)
        pattern.nnz = gpu_matrix.nonzeroes();
    }

    // Setup host arrays
    pattern.ncols.resize(par_info.N_owned);
    pattern.rows.resize(par_info.N_owned);
    pattern.cols.resize(pattern.nnz);

    // Get row pointers and column indices from GPU matrix (one-time host copy during setup)
    auto host_row_ptrs = gpu_matrix.getRowIndices().asStdVector();
    auto host_col_indices = gpu_matrix.getColumnIndices().asStdVector();

    int pos = 0;
    for (int rind = 0; rind < static_cast<int>(gpu_matrix.N()); ++rind) {
        const int local_rowIdx = par_info.local_dune_to_local_hypre[rind];

        // For owner_first=true: skip ghost rows entirely
        // For owner_first=false: process all rows (owned + ghost)
        if (owner_first && local_rowIdx < 0) {
            continue;
        }

        const int row_start = host_row_ptrs[rind];
        const int row_end = host_row_ptrs[rind + 1];
        const int num_cols = row_end - row_start;

        if (local_rowIdx >= 0) {
            // This is an owned row - record its metadata
            const int global_rowIdx = par_info.local_dune_to_global_hypre[rind];
            pattern.rows[local_rowIdx] = global_rowIdx;
            pattern.ncols[local_rowIdx] = num_cols;
        }

        // Add column indices for this row
        for (int col_idx = row_start; col_idx < row_end; ++col_idx) {
            const int colIdx = host_col_indices[col_idx];
            const int global_colIdx = par_info.local_dune_to_global_hypre[colIdx];
            assert(global_colIdx >= 0);
            pattern.cols[pos++] = global_colIdx;
        }
    }

    return pattern;
}

/**
 * @brief Compute row indexes for GPU matrix with ownership mapping
 *
 * GPU-specialized version that uses CSR row pointers to create row_indexes
 * pointing directly into the FULL GPU matrix data for zero-copy operation.
 * Same principle as CPU version but leverages GPU matrix's CSR structure.
 *
 * @param gpu_matrix GPU matrix to analyze
 * @param local_dune_to_local_hypre Mapping from Dune indices to Hypre indices (-1 for ghost)
 * @return Vector containing row indexes that map to full GPU matrix data positions
 */
template <typename T, bool ForceLegacy>
std::vector<HYPRE_Int>
computeRowIndexesWithMappingGpu(const GpuSparseMatrixWrapper<T, ForceLegacy>& gpu_matrix,
                                const std::vector<int>& local_dune_to_local_hypre)
{
    const int N = std::count_if(
        local_dune_to_local_hypre.begin(), local_dune_to_local_hypre.end(), [](int val) { return val >= 0; });
    std::vector<HYPRE_Int> row_indexes(N);

    // Use pre-computed BSR row pointers (already contain row starting positions)
    auto host_row_ptrs = gpu_matrix.getRowIndices().asStdVector();

    // Map each owned Hypre row to its starting position in the FULL (including ghost) GPU matrix
    for (int dune_row_idx = 0; dune_row_idx < static_cast<int>(gpu_matrix.N()); ++dune_row_idx) {
        const int hypre_row_idx = local_dune_to_local_hypre[dune_row_idx];

        if (hypre_row_idx >= 0) {
            // This is an owned row - record where its data starts in the FULL GPU matrix
            // Use hypre_row_idx as index (maps to Hypre ordering up to N_owned)
            row_indexes[hypre_row_idx] = host_row_ptrs[dune_row_idx];
        }
        // Non-owned rows create natural gaps in the indexing
    }
    return row_indexes;
}

} // namespace Opm::gpuistl::HypreInterface

#endif // OPM_HYPRE_SETUP_GPU_HPP
