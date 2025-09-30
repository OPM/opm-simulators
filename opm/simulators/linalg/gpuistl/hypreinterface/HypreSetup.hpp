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

#ifndef OPM_HYPRE_SETUP_HPP
#define OPM_HYPRE_SETUP_HPP

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_type_detection.hpp>
#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreDataStructures.hpp>
#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreErrorHandling.hpp>

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
#include <numeric>

namespace Opm::gpuistl::HypreInterface
{

// GPU-specific helper functions
template <typename T>
SparsityPattern setupSparsityPatternFromGpuMatrix(const GpuSparseMatrixWrapper<T>& gpu_matrix,
                                                  const ParallelInfo& par_info,
                                                  bool owner_first);

template <typename T>
std::vector<HYPRE_Int> computeRowIndexesWithMappingGpu(const GpuSparseMatrixWrapper<T>& gpu_matrix,
                                                       const std::vector<int>& local_dune_to_local_hypre);

// Serial helper functions
ParallelInfo setupHypreParallelInfoSerial(HYPRE_Int N);

// Parallel helper functions
template <typename CommType, typename MatrixType>
ParallelInfo setupHypreParallelInfoParallel(const CommType& comm, const MatrixType& matrix);

template <typename MatrixType>
SparsityPattern
setupSparsityPatternFromCpuMatrix(const MatrixType& matrix, const ParallelInfo& par_info, bool owner_first);

template <typename MatrixType>
std::vector<HYPRE_Int> computeRowIndexesWithMappingCpu(const MatrixType& matrix,
                                                       const std::vector<HYPRE_Int>& ncols,
                                                       const std::vector<int>& local_dune_to_local_hypre,
                                                       bool owner_first);

template <typename MatrixType>
std::vector<HYPRE_Int> computeRowIndexesWithMappingCpu(const MatrixType& matrix,
                                                       const std::vector<int>& local_dune_to_local_hypre);
/**
 * @brief Initialize the Hypre library and set memory/execution policy
 *
 * @param use_gpu_backend Whether to use GPU backend acceleration
 * @throws HypreError if initialization fails
 */
inline void
initialize([[maybe_unused]] bool use_gpu_backend)
{
    // Set memory location and execution policy
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
    if (use_gpu_backend) {
        OPM_HYPRE_SAFE_CALL(HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE));
        OPM_HYPRE_SAFE_CALL(HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE));
        // use hypre's SpGEMM instead of vendor implementation
        OPM_HYPRE_SAFE_CALL(HYPRE_SetSpGemmUseVendor(false));
        // use cuRand for PMIS
        OPM_HYPRE_SAFE_CALL(HYPRE_SetUseGpuRand(1));
        OPM_HYPRE_SAFE_CALL(HYPRE_DeviceInitialize());
        OPM_HYPRE_SAFE_CALL(HYPRE_PrintDeviceInfo());
    } else
#endif
    {
        OPM_HYPRE_SAFE_CALL(HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST));
        OPM_HYPRE_SAFE_CALL(HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST));
    }
}

/**
 * @brief Create Hypre BoomerAMG solver
 *
 * @return HYPRE_Solver The created solver handle
 * @throws HypreError if solver creation fails
 */
inline HYPRE_Solver
createAMGSolver()
{
    HYPRE_Solver solver;
    OPM_HYPRE_SAFE_CALL(HYPRE_BoomerAMGCreate(&solver));
    return solver;
}

/**
 * @brief Set solver parameters from property tree
 *
 * @param solver The solver handle
 * @param prm Property tree containing configuration parameters
 * @param use_gpu_backend Whether using GPU backend
 * @throws HypreError if parameter setting fails
 */
inline void
setSolverParameters(HYPRE_Solver solver, const PropertyTree& prm, bool use_gpu_backend)
{
    // Set parameters from property tree with defaults
    HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetPrintLevel(solver, prm.get<int>("print_level", 0)));
    HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetMaxIter(solver, prm.get<int>("max_iter", 1)));
    HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetStrongThreshold(solver, prm.get<double>("strong_threshold", 0.5)));
    HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggTruncFactor(solver, prm.get<double>("agg_trunc_factor", 0.3)));
    HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetInterpType(solver, prm.get<int>("interp_type", 6)));
    HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetMaxLevels(solver, prm.get<int>("max_levels", 15)));
    HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetTol(solver, prm.get<double>("tolerance", 0.0)));

    if (use_gpu_backend) {
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetRelaxType(solver, prm.get<int>("relax_type", 16)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetCoarsenType(solver, prm.get<int>("coarsen_type", 8)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggNumLevels(solver, prm.get<int>("agg_num_levels", 0)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggInterpType(solver, prm.get<int>("agg_interp_type", 6)));
        // Keep transpose to avoid SpMTV
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetKeepTranspose(solver, true));
    } else {
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetRelaxType(solver, prm.get<int>("relax_type", 13)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetCoarsenType(solver, prm.get<int>("coarsen_type", 10)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggNumLevels(solver, prm.get<int>("agg_num_levels", 1)));
        HYPRE_SAFE_CALL(HYPRE_BoomerAMGSetAggInterpType(solver, prm.get<int>("agg_interp_type", 4)));
    }
}

/**
 * @brief Create Hypre matrix
 *
 * @param N Number of rows/columns
 * @param dof_offset Offset for the matrix
 * @param comm The communicator
 * @return HYPRE_IJMatrix The created matrix handle
 * @throws HypreError if matrix creation fails
 */
template <typename CommType>
HYPRE_IJMatrix
createMatrix(HYPRE_Int N, HYPRE_Int dof_offset, const CommType& comm)
{
    HYPRE_IJMatrix matrix;
    MPI_Comm mpi_comm;
    if constexpr (std::is_same_v<CommType, Dune::Amg::SequentialInformation>) {
        mpi_comm = MPI_COMM_SELF;
    } else {
        mpi_comm = comm.communicator();
    }
    OPM_HYPRE_SAFE_CALL(
        HYPRE_IJMatrixCreate(mpi_comm, dof_offset, dof_offset + (N - 1), dof_offset, dof_offset + (N - 1), &matrix));
    OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixSetObjectType(matrix, HYPRE_PARCSR));
    OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixInitialize(matrix));
    return matrix;
}

/**
 * @brief Create Hypre vector
 *
 * @param N Size of the vector
 * @param dof_offset Offset for the vector
 * @param comm The communicator
 * @return HYPRE_IJVector The created vector handle
 * @throws HypreError if vector creation fails
 */
template <typename CommType>
HYPRE_IJVector
createVector(HYPRE_Int N, HYPRE_Int dof_offset, const CommType& comm)
{
    HYPRE_IJVector vector;
    MPI_Comm mpi_comm;
    if constexpr (std::is_same_v<CommType, Dune::Amg::SequentialInformation>) {
        mpi_comm = MPI_COMM_SELF;
    } else {
        mpi_comm = comm.communicator();
    }
    OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorCreate(mpi_comm, dof_offset, dof_offset + (N - 1), &vector));
    OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorSetObjectType(vector, HYPRE_PARCSR));
    OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorInitialize(vector));
    return vector;
}

/**
 * @brief Destroy Hypre solver
 *
 * @param solver The solver handle to destroy
 * @throws HypreError if solver destruction fails
 */
inline void
destroySolver(HYPRE_Solver solver)
{
    if (solver) {
        OPM_HYPRE_SAFE_CALL(HYPRE_BoomerAMGDestroy(solver));
    }
}

/**
 * @brief Destroy Hypre matrix
 *
 * @param matrix The matrix handle to destroy
 * @throws HypreError if matrix destruction fails
 */
inline void
destroyMatrix(HYPRE_IJMatrix matrix)
{
    if (matrix) {
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixDestroy(matrix));
    }
}

/**
 * @brief Destroy Hypre vector
 *
 * @param vector The vector handle to destroy
 * @throws HypreError if vector destruction fails
 */
inline void
destroyVector(HYPRE_IJVector vector)
{
    if (vector) {
        OPM_HYPRE_SAFE_CALL(HYPRE_IJVectorDestroy(vector));
    }
}

/**
 * @brief Setup parallel information for Hypre (automatically detects serial/parallel)
 *
 * @param comm The communication object (serial or parallel)
 * @param matrix The matrix to analyze
 * @return ParallelInfo structure containing mappings and offsets
 */
template <typename CommType, typename MatrixType>
ParallelInfo
setupHypreParallelInfo(const CommType& comm, const MatrixType& matrix)
{
    if constexpr (std::is_same_v<CommType, Dune::Amg::SequentialInformation>) {
        return setupHypreParallelInfoSerial(static_cast<HYPRE_Int>(matrix.N()));
    } else {
        return setupHypreParallelInfoParallel(comm, matrix);
    }
}

/**
 * @brief Setup parallel information for Hypre in serial case
 *
 * @param N Number of rows
 * @return ParallelInfo structure containing mappings and offsets
 */
inline ParallelInfo
setupHypreParallelInfoSerial(HYPRE_Int N)
{
    ParallelInfo info;
    info.N_owned = N;

    info.local_dune_to_local_hypre.resize(N);
    info.local_dune_to_global_hypre.resize(N);
    info.local_hypre_to_local_dune.resize(N);

    std::iota(info.local_dune_to_local_hypre.begin(), info.local_dune_to_local_hypre.end(), 0);
    std::iota(info.local_hypre_to_local_dune.begin(), info.local_hypre_to_local_dune.end(), 0);
    std::iota(info.local_dune_to_global_hypre.begin(), info.local_dune_to_global_hypre.end(), 0);

    info.dof_offset = 0;
    info.owner_first = true;

    return info;
}

/**
 * @brief Create mappings between Dune and HYPRE indexing for parallel decomposition
 *
 * This function interfaces between Dune's distributed matrix representation and HYPRE's
 * global indexing requirements. Note that Dune uses local indices with owner/ghost
 * classification, while HYPRE requires globally consistent indices for all DOFs.
 *
 * ## Index System Overview:
 *
 * **Dune Local Indices**: Each MPI process has local DOFs indexed 0..N_local-1, where:
 * - "Owner" DOFs: This process is responsible for their values/updates
 * - "Ghost" DOFs: Copies of DOFs owned by other processes (for stencil access)
 * - These can be interleaved or arranged in "owner-first" order where all owned DOFs come first
 *
 * **HYPRE Global Indices**: All DOFs across all processes have unique global indices:
 * - Process 0 owns global indices [0, N_owned_0-1]
 * - Process 1 owns global indices [N_owned_0, N_owned_0 + N_owned_1-1]
 * - Process k owns global indices [offset_k, offset_k + N_owned_k-1]
 *
 * ## Index Mappings Created:
 *
 * 1. **local_dune_to_local_hypre**: Maps Dune local index → HYPRE local index
 *    - Size: N_local (all DOFs)
 *    - Value: [0..N_owned-1] for owned DOFs, -1 for ghost DOFs
 *    - Purpose: Identify owned DOFs and their compact local ordering
 *
 * 2. **local_dune_to_global_hypre**: Maps Dune local index → HYPRE global index
 *    - Size: N_local (all DOFs)
 *    - Value: [offset..offset+N_owned-1] for owned, actual global index for ghost
 *    - Purpose: Matrix assembly and vector operations with global indexing
 *
 * 3. **local_hypre_to_local_dune**: Maps HYPRE local index → Dune local index
 *    - Size: N_owned (owned DOFs only)
 *    - Purpose: Transfer data from HYPRE back to Dune structures
 *
 * ## Algorithm Steps:
 *
 * 1. **Ownership Detection**: Scan Dune's index set to identify owner vs ghost DOFs
 * 2. **Local Reordering**: Create compact [0..N_owned-1] indexing for owned DOFs
 * 3. **Owner-First Detection**: Determine if all owned DOFs appear before ghost DOFs
 * 4. **Global Offset Calculation**: Coordinate with other processes to assign global ranges
 * 5. **Ghost Communication**: Exchange global indices for ghost DOFs
 *
 * ## Example (2 processes, 3 DOFs each):
 *
 * Process 0: Dune local [0,1,2] → owned=[0,1], ghost=[2 from P1]
 * - local_dune_to_local_hypre = [0, 1, -1]
 * - local_dune_to_global_hypre = [0, 1, 2] (after communication)
 * - local_hypre_to_local_dune = [0, 1]
 * - N_owned=2, dof_offset=0, owner_first=true
 *
 * Process 1: Dune local [0,1,2] → owned=[1,2], ghost=[0 from P0]
 * - local_dune_to_local_hypre = [-1, 0, 1]
 * - local_dune_to_global_hypre = [0, 2, 3] (after communication)
 * - local_hypre_to_local_dune = [1, 2]
 * - N_owned=2, dof_offset=2, owner_first=false
 *
 * @param comm Dune communication object with parallel index set
 * @param matrix Distributed matrix for consistency checks
 * @return ParallelInfo with all mappings, offsets, and metadata
 *
 * @note May modify comm if index set holes are detected (rare in OPM)
 */
template <typename CommType, typename MatrixType>
inline ParallelInfo
setupHypreParallelInfoParallel(const CommType& comm, const MatrixType& matrix)
{
    ParallelInfo info;
    const auto& collective_comm = comm.communicator();

    // Initialize mapping arrays to not owned (-1) state
    info.local_dune_to_local_hypre.resize(comm.indexSet().size(), -1);
    info.local_dune_to_global_hypre.resize(comm.indexSet().size(), -1);

    // Handle edge case: ensure index set covers all matrix rows
    if (!(matrix.N() == comm.indexSet().size())) {
          // in OPM this will likely not be trigged
          // ensure no holes in index sett
          const_cast<CommType&>(comm).buildGlobalLookup(matrix.N()); // need?
          Dune::Amg::MatrixGraph<MatrixType> graph(const_cast<MatrixType&>(matrix)); // do not know why not const ref is sufficient
          Dune::fillIndexSetHoles(graph, const_cast<CommType&>(comm));
          assert(matrix.N() == comm.indexSet().size());
    }

    // STEP 1: Ownership Detection
    // Scan Dune's index set to identify which DOFs this process owns
    // Note: iteration order in index set is NOT sequential by local index
    for (const auto& ind : comm.indexSet()) {
        int local_ind = ind.local().local();
        if (ind.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner) {
            // Mark as owned (temporarily use 1, will be replaced with proper local index)
            info.local_dune_to_local_hypre[local_ind] = 1;
        } else {
            // Mark as ghost/non-owned
            info.local_dune_to_local_hypre[local_ind] = -1;
        }
    }

    // STEP 2: Local Reordering & Owner-First Detection
    // Create compact [0..N_owned-1] local HYPRE indexing for owned DOFs
    // Simultaneously detect if owned DOFs appear before all ghost DOFs
    bool owner_first = true;
    bool visited_ghost = false;  // Have we seen any ghost DOF yet?
    std::size_t count = 0;      // Counter for owned DOFs

    for (std::size_t i = 0; i < info.local_dune_to_local_hypre.size(); ++i) {
        if (info.local_dune_to_local_hypre[i] < 0) {
            visited_ghost = true;
        } else {
            // This is an owned DOF - assign its local HYPRE index
            info.local_dune_to_local_hypre[i] = count;
            // Store the inverse mapping
            info.local_hypre_to_local_dune.push_back(i);

            // Check if we've seen ghost DOFs before this owner
            owner_first = owner_first && !visited_ghost;
            count += 1;
        }
    }

    // Owner first need other copying of data
    info.owner_first = owner_first;
    info.N_owned = count;

    // STEP 3: Global Offset Calculation
    // Coordinate with other processes to determine global index ranges
    // Each process owns a contiguous range of global indices
    std::vector<int> dof_counts_per_process(collective_comm.size());
    collective_comm.allgather(&info.N_owned, 1, dof_counts_per_process.data());

    // Calculate this process's global offset (sum of DOFs in processes with lower rank)
    info.dof_offset = std::accumulate(dof_counts_per_process.begin(),
                                      dof_counts_per_process.begin() + collective_comm.rank(),
                                      0);

    // STEP 4: Create Global Indices for Owned DOFs
    // Convert local HYPRE indices to global HYPRE indices by adding offset
    for (std::size_t i = 0; i < info.local_dune_to_local_hypre.size(); ++i) {
        if (info.local_dune_to_local_hypre[i] >= 0) {
            // Owned DOF: global index = local HYPRE index + this process's offset
            info.local_dune_to_global_hypre[i] = info.local_dune_to_local_hypre[i] + info.dof_offset;
        } else {
            info.local_dune_to_global_hypre[i] = -1;
        }
    }

    if (collective_comm.rank() > 0) {
        assert(info.dof_offset > 0);
    }

    // STEP 5: Exchange global indices for ghost DOFs
    // After this call, ghost DOFs will have their correct global indices
    comm.copyOwnerToAll(info.local_dune_to_global_hypre, info.local_dune_to_global_hypre);

    return info;
}

/**
 * @brief Setup sparsity pattern from matrix (automatically detects CPU/GPU type)
 *
 * @param matrix The matrix to analyze (CPU or GPU type)
 * @param par_info Parallel information structure
 * @param owner_first Whether all owned DOFs come first in the ordering
 * @return Sparsity pattern information
 */
template <typename MatrixType>
SparsityPattern
setupSparsityPattern(const MatrixType& matrix, const ParallelInfo& par_info, bool owner_first)
{
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
    if constexpr (is_gpu_type<MatrixType>::value) {
        return setupSparsityPatternFromGpuMatrix(matrix, par_info, owner_first);
    } else
#endif
    {
        return setupSparsityPatternFromCpuMatrix(matrix, par_info, owner_first);
    }
}

/**
 * @brief Setup sparsity pattern from CPU matrix (BCRSMatrix)
 *
 * @param matrix The matrix to analyze
 * @param par_info Parallel information structure
 * @param owner_first Whether all owned DOFs come first in the ordering
 * @return Sparsity pattern information
 */
template <typename MatrixType>
SparsityPattern
setupSparsityPatternFromCpuMatrix(const MatrixType& matrix, const ParallelInfo& par_info, bool owner_first)
{
    SparsityPattern pattern;

    // Determine the size for cols array based on owner_first
    if (owner_first) {
        std::size_t cols_size = 0;
        // For owner_first=true case, we need to calculate how many owned entries there are
        for (auto row = matrix.begin(); row != matrix.end(); ++row) {
            const int rowIdx = row.index();
            if (par_info.local_dune_to_local_hypre[rowIdx] >= 0) {
                cols_size += row->size();
            }
        }
        pattern.nnz = cols_size;
    } else {
        // Full matrix space case: all entries (including gaps)
        pattern.nnz = matrix.nonzeroes();
    }

    // Setup host arrays
    pattern.ncols.resize(par_info.N_owned);
    pattern.rows.resize(par_info.N_owned);
    pattern.cols.resize(pattern.nnz);

    int pos = 0;
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        const int rind = row.index();
        const int local_rowIdx = par_info.local_dune_to_local_hypre[rind];

        // For owner_first=true: skip ghost rows entirely
        // For owner_first=false: process all rows (owned + ghost)
        if (owner_first && local_rowIdx < 0) {
            continue;
        }

        if (local_rowIdx >= 0) {
            // This is an owned row - record its metadata
            const int global_rowIdx = par_info.local_dune_to_global_hypre[rind];
            pattern.rows[local_rowIdx] = global_rowIdx;
            pattern.ncols[local_rowIdx] = row->size();
        }

        // Add column indices for this row
        for (auto col = row->begin(); col != row->end(); ++col) {
            const int global_colIdx = par_info.local_dune_to_global_hypre[col.index()];
            assert(global_colIdx >= 0);
            pattern.cols[pos++] = global_colIdx;
        }
    }

    return pattern;
}

#if HYPRE_USING_CUDA || HYPRE_USING_HIP
/**
 * @brief Setup sparsity pattern from GPU matrix (GpuSparseMatrix)
 *
 * @param gpu_matrix The GPU matrix to analyze
 * @param par_info Parallel information structure
 * @param owner_first Whether all owned DOFs come first in the ordering
 * @return Sparsity pattern information
 */
template <typename T>
SparsityPattern
setupSparsityPatternFromGpuMatrix(const GpuSparseMatrixWrapper<T>& gpu_matrix,
                                  const ParallelInfo& par_info,
                                  bool owner_first)
{
    SparsityPattern pattern;

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
#endif // HYPRE_USING_CUDA || HYPRE_USING_HIP

/**
 * @brief Compute row indexes for HYPRE_IJMatrixSetValues2
 *
 * Chooses the appropriate row index computation method based on owner_first and matrix type.
 * - For owner_first=true: Simple prefix sum for contiguous data
 * - For owner_first=false: Mapping-based approach to handle gaps from ghost rows
 *
 * @param matrix The matrix to analyze
 * @param ncols Array containing number of columns per row
 * @param local_dune_to_local_hypre Mapping from Dune indices to Hypre indices
 * @param owner_first Whether all owned DOFs come first in the ordering
 * @return Vector containing row indexes for HYPRE_IJMatrixSetValues2
 */
template <typename MatrixType>
std::vector<HYPRE_Int>
computeRowIndexes(const MatrixType& matrix,
                  const std::vector<HYPRE_Int>& ncols,
                  const std::vector<int>& local_dune_to_local_hypre,
                  bool owner_first)
{
    if (owner_first) {
        // Simple contiguous case: prefix sum of ncols
        std::vector<HYPRE_Int> row_indexes(ncols.size());
        row_indexes[0] = 0;
        for (std::size_t i = 1; i < ncols.size(); ++i) {
            row_indexes[i] = row_indexes[i - 1] + ncols[i - 1];
        }
        return row_indexes;
    } else {
        // We need to compute the row indexes with mapping since we have gaps in the data
#if HYPRE_USING_CUDA || HYPRE_USING_HIP
        if constexpr (is_gpu_type<MatrixType>::value) {
            return computeRowIndexesWithMappingGpu(matrix, local_dune_to_local_hypre);
        } else
#endif // HYPRE_USING_CUDA || HYPRE_USING_HIP
        {
            return computeRowIndexesWithMappingCpu(matrix, local_dune_to_local_hypre);
        }
    }
}

/**
 * @brief Compute row indexes for CPU matrix with ownership mapping
 *
 * Creates row_indexes that point directly into the FULL matrix data, enabling using the original
  * matrix data without any data copying. This is necessary for some cases where ghost rows can come before all owned rows.
 *
 * Example with local_dune_to_local_hypre = [-1, 0, -1, 1]:
 *   Original matrix data (with gaps for ghost rows):
 *   Row 0 (ghost): {_, _}       at positions 0-1 (skipped)
 *   Row 1 (owned): {1.0, 2.0}   at positions 2-3
 *   Row 2 (ghost): {_}          at position 4 (skipped)
 *   Row 3 (owned): {3.0}        at position 5
 *   Resulting row_indexes: {2, 5}
 *   This tells HYPRE that Row 0 data starts at position 2, Row 1 data starts at position 5
 *
 * @param matrix Matrix to analyze
 * @param local_dune_to_local_hypre Mapping from Dune indices to Hypre indices (-1 for ghost)
 * @return Vector containing row indexes that map to full matrix data positions
 */
template <typename MatrixType>
std::vector<HYPRE_Int>
computeRowIndexesWithMappingCpu(const MatrixType& matrix, const std::vector<int>& local_dune_to_local_hypre)
{
    const int N = std::count_if(
        local_dune_to_local_hypre.begin(), local_dune_to_local_hypre.end(), [](int val) { return val >= 0; });
    std::vector<HYPRE_Int> row_indexes(N);
    int data_position = 0; // Current position in the FULL (including ghost) matrix data
    // Manually compute row starting positions by iterating through all matrix rows
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        const int dune_row_idx = row.index();
        const int hypre_row_idx = local_dune_to_local_hypre[dune_row_idx];

        if (hypre_row_idx >= 0) {
            // This is an owned row - record where its data starts in the FULL matrix
            // Use hypre_row_idx as index (maps to Hypre ordering up to N_owned)
            row_indexes[hypre_row_idx] = data_position;
        }
        // Always advance position (including ghost rows - this creates gaps for owned-only access)
        data_position += row->size();
    }
    return row_indexes;
}

#if HYPRE_USING_CUDA || HYPRE_USING_HIP
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
template <typename T>
std::vector<HYPRE_Int>
computeRowIndexesWithMappingGpu(const GpuSparseMatrixWrapper<T>& gpu_matrix, const std::vector<int>& local_dune_to_local_hypre)
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
#endif // HYPRE_USING_CUDA || HYPRE_USING_HIP

} // namespace Opm::gpuistl::HypreInterface

#endif // OPM_HYPRE_SETUP_HPP
