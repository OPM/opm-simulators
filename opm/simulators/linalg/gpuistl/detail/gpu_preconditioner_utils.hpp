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

#ifndef OPM_GPUISTL_DETAIL_GPU_PRECONDITIONER_UTILS_HEADER
#define OPM_GPUISTL_DETAIL_GPU_PRECONDITIONER_UTILS_HEADER

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>

#include <dune/istl/bcrsmatrix.hh>


namespace Opm::gpuistl::detail {

/**
 * @brief Utility functions for GPU preconditioner creation
 *
 * These functions are temporary workarounds needed because some GPU preconditioners
 * require CPU matrices for initial setup (e.g., graph coloring). They will be
 * removed once GPU preconditioners have native GPU constructors.
 */


/**
 * This function creates a CPU matrix from the operator holding a GPU matrix.
 *
 * This is a workaround for now since some of the GPU preconditioners need a
 * CPU matrix for the initial setup (graph coloring). The CPU matrix is only
 * used in the constructor, **not** in the update function or the apply function.
 */
template<class Operator, class BlockType>
Dune::BCRSMatrix<BlockType> makeCPUMatrix(const Operator& op) {
    // TODO: Make this more efficient. Maybe we can simply copy the memory areas directly?
    //       Do note that this function is anyway going away when we have a GPU
    //       constructor for the preconditioners, so it is not a priority.
    const auto& gpuMatrix = op.getmat();

    const auto nonZeros = gpuMatrix.getNonZeroValues().asStdVector();
    const auto rowIndices = gpuMatrix.getRowIndices().asStdVector();
    const auto columnIndices = gpuMatrix.getColumnIndices().asStdVector();

    const auto numberOfNonZeroes = gpuMatrix.nonzeroes();
    const auto N = gpuMatrix.N();

    Dune::BCRSMatrix<BlockType> matrix(N, N, numberOfNonZeroes, Dune::BCRSMatrix<BlockType>::row_wise);
    for (auto row = matrix.createbegin(); row != matrix.createend(); ++row) {
        for (auto j = rowIndices[row.index()]; j != rowIndices[row.index() + 1]; ++j) {
            const auto columnIndex = columnIndices[j];
            row.insert(columnIndex);
        }
    }

    for (std::size_t i = 0; i < N; ++i) {
        for (auto j = rowIndices[i]; j != rowIndices[i + 1]; ++j) {
            const auto columnIndex = columnIndices[j];
            // Now it gets a bit tricky, first we need to fetch the block matrix
            BlockType blockMatrix;
            constexpr static auto rows = BlockType::rows;

            for (std::size_t k = 0; k < rows; ++k) {
                for (std::size_t l = 0; l < rows; ++l) {
                    blockMatrix[k][l] = nonZeros[j * rows * rows + k * rows + l];
                }
            }
            matrix[i][columnIndex] = blockMatrix;
        }
    }

    return matrix;
}


} // namespace Opm::gpuistl::detail

#endif // OPM_GPUISTL_DETAIL_GPU_PRECONDITIONER_UTILS_HEADER
