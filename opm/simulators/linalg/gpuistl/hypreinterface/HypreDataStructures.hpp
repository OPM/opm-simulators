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

#ifndef OPM_HYPRE_DATA_STRUCTURES_HPP
#define OPM_HYPRE_DATA_STRUCTURES_HPP

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>
#include <vector>

namespace Opm::gpuistl
{

/**
 * @brief Parallel domain decomposition information for HYPRE-Dune interface
 *
 * Contains all mappings and metadata to translate between Dune's distributed
 * indexing and HYPRE's global indexing in parallel MPI environments.
 */
 struct ParallelInfo {
    /** @brief Mapping from local Dune indices to local HYPRE indices
     *
     * Size: total local DOFs (owned + ghost)
     * Value: >=0 for owned DOFs (local HYPRE index), -1 for ghost/non-owned DOFs
     * Used to identify which DOFs are owned and their local HYPRE ordering
     */
    std::vector<int> local_dune_to_local_hypre;

    /** @brief Mapping from local Dune indices to global HYPRE indices
     *
     * Size: total local DOFs (owned + ghost)
     * Used for matrix assembly and vector operations with global indexing
     */
    std::vector<int> local_dune_to_global_hypre;

    /** @brief Mapping from local HYPRE indices to local Dune indices
     *
     * Size: N_owned (only owned DOFs)
     * Inverse mapping of local_dune_to_local_hypre for owned DOFs only
     * Used when transferring data from HYPRE back to Dune data structures
     */
    std::vector<int> local_hypre_to_local_dune;

    /** @brief Number of DOFs owned by this MPI process */
    HYPRE_Int N_owned;

    /** @brief Global index offset for this process's owned DOFs
     *
     * Starting global index for this process. Owned DOFs have global indices
     * in range [dof_offset, dof_offset + N_owned - 1].
     */
    HYPRE_Int dof_offset;

    /** @brief Whether owned DOFs appear first in local Dune ordering
     *
     * true: All owned DOFs have indices 0..N_owned-1, ghost DOFs follow
     * false: Owned and ghost DOFs are interleaved in local Dune ordering
     * Affects data layout optimization strategies and transfer efficiency
     */
    bool owner_first;
};

/**
 * @brief Compressed Sparse Row (CSR) sparsity pattern for HYPRE matrix assembly
 *
 * It represents only the owned rows in CSR-like
 * format with global HYPRE indexing for both rows and columns.
 */
struct SparsityPattern {
    /** @brief Non-zero entries per owned row (size: N_owned) */
    std::vector<HYPRE_Int> ncols;

    /** @brief Global row indices for owned rows (size: N_owned) */
    std::vector<HYPRE_BigInt> rows;

    /** @brief Global column indices in CSR format (size: nnz) */
    std::vector<HYPRE_BigInt> cols;

    /** @brief Number of non-zero entries in matrix */
    HYPRE_Int nnz;
};

/**
 * @brief Host arrays for HYPRE matrix and vector data transfers
 *
 * Pre-computed arrays to efficiently transfer data to/from HYPRE without
 * repeated calculations during solve operations.
 */
struct HypreHostDataArrays {
    /** @brief Pre-computed row start indexes for HYPRE_IJMatrixSetValues2
     *
     * For owner_first=true: prefix sum of ncols (contiguous data)
     * For owner_first=false: maps owned rows to positions in full matrix
     * Used to index into the full matrix data (including ghost rows)
     */
    std::vector<HYPRE_Int> row_indexes;

    /** @brief Global DOF indices for owned degrees of freedom
     *
     * Sequential indices starting from dof_offset, used for vector operations
     * Size: N_owned, contains [dof_offset, dof_offset+1, ..., dof_offset+N_owned-1]
     */
    std::vector<HYPRE_BigInt> indices;

    /** @brief Temporary buffer for vector values in non-owner-first ordering
     *
     * Only allocated when owner_first=false, to size N_owned. Used to collect owned
     * DOF values from scattered positions in the input vector into a contiguous array.
     * Size: N_owned when allocated, empty otherwise.
     */
    std::vector<HYPRE_Real> continuous_vector_values;
};

/**
 * @brief GPU device memory arrays for HYPRE operations with GPU backend
 *
 * These arrays mirror the host data but reside in GPU memory. All pointers are
 * managed manually using HYPRE's memory management functions.
 */
struct HypreDeviceDataArrays {
    /** @brief Mirrors host data arrays */
    HYPRE_Int* ncols_device = nullptr;
    HYPRE_BigInt* rows_device = nullptr;
    HYPRE_BigInt* cols_device = nullptr;
    HYPRE_Int* row_indexes_device = nullptr;
    HYPRE_BigInt* indices_device = nullptr;

    /** @brief Device buffer for vector operations
     * Used when input type and backend are different, for transferring values between host and device.
     * Size: N_owned when allocated, empty otherwise.
     */
    HYPRE_Real* vector_buffer_device = nullptr;

    /** @brief Device buffer for matrix values, only needed for CPU input + GPU backend
     *
     * Size: nnz when allocated, empty otherwise.
     */
    HYPRE_Real* matrix_buffer_device = nullptr;
};
} // namespace Opm::gpuistl

#endif // OPM_HYPRE_DATA_STRUCTURES_HPP
