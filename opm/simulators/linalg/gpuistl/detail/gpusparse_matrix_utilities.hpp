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
#ifndef OPM_GPUISTL_GPUSPARSE_MATRIX_UTILITIES_HPP
#define OPM_GPUISTL_GPUSPARSE_MATRIX_UTILITIES_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cusparse_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>

#include <fmt/core.h>
#include <memory>
#include <vector>

#include <cusparse.h>

namespace Opm::gpuistl::detail
{

/**
 * @brief Create RAII-managed cuSPARSE dense vector descriptor
 * @return unique_ptr with custom deleter for automatic cleanup
 */
inline auto
makeSafeVectorDescriptor()
{
    auto deleter = [](cusparseDnVecDescr_t* desc) {
        if (desc && *desc) {
            OPM_CUSPARSE_WARN_IF_ERROR(cusparseDestroyDnVec(*desc));
        }
        delete desc;
    };
    return std::unique_ptr<cusparseDnVecDescr_t, decltype(deleter)>(new cusparseDnVecDescr_t {nullptr}, deleter);
}

/**
 * @brief Create RAII-managed cuSPARSE sparse matrix descriptor
 * @return unique_ptr with custom deleter for automatic cleanup
 */
inline auto
makeSafeMatrixDescriptor()
{
    auto deleter = [](cusparseSpMatDescr_t* desc) {
        if (desc && *desc) {
            OPM_CUSPARSE_WARN_IF_ERROR(cusparseDestroySpMat(*desc));
        }
        delete desc;
    };
    return std::unique_ptr<cusparseSpMatDescr_t, decltype(deleter)>(new cusparseSpMatDescr_t {nullptr}, deleter);
}

/**
 * @brief Extract non-zero values from a DUNE matrix in the order expected by cuSPARSE
 *
 * This function iterates through the matrix and extracts the non-zero values in the
 * Block Compressed Sparse Row (BSR) format order that cuSPARSE expects.
 *
 * @tparam T the element type (float or double)
 * @tparam M the matrix type (assumed to be DUNE::BCRSMatrix compatible)
 * @param matrix the input matrix
 * @return vector containing the non-zero elements in BSR order
 */
template <class T, class M>
std::vector<T>
extractNonzeroValues(const M& matrix)
{
    const std::size_t blockSize = matrix[0][0].N();
    const std::size_t numberOfNonzeroBlocks = matrix.nonzeroes();
    const std::size_t numberOfNonzeroElements = blockSize * blockSize * numberOfNonzeroBlocks;

    std::vector<T> nonZeroElementsData(numberOfNonzeroElements);

    std::size_t ind = 0;
    for (auto& row : matrix) {
        for (const auto& element : row) {
            for (std::size_t c = 0; c < blockSize; ++c) {
                for (std::size_t d = 0; d < blockSize; ++d) {
                    nonZeroElementsData[ind++] = element[c][d];
                }
            }
        }
    }

    return nonZeroElementsData;
}

/**
 * @brief Extract sparsity pattern from a DUNE matrix
 *
 * This function extracts the row and column indices that define the sparsity pattern
 * of a DUNE matrix in the CSR format expected by cuSPARSE.
 *
 * @tparam M the matrix type (assumed to be DUNE::BCRSMatrix compatible)
 * @param matrix the input matrix
 * @return pair of (columnIndices, rowIndices) in CSR format
 */
template <class M>
std::pair<std::vector<int>, std::vector<int>>
extractSparsityPattern(const M& matrix)
{
    std::vector<int> columnIndices(matrix.nonzeroes());
    std::vector<int> rowIndices(matrix.N() + 1);

    rowIndices[0] = 0;
    std::size_t colIndex = 0;
    std::size_t rowIndex = 1;
    for (const auto& row : matrix) {
        for (auto columnIterator = row.begin(); columnIterator != row.end(); ++columnIterator) {
            columnIndices[colIndex++] = columnIterator.index();
        }
        rowIndices[rowIndex++] = colIndex;
    }

    return {columnIndices, rowIndices};
}

/**
 * @brief Find pointer to first matrix element for direct memory access
 *
 * This function finds a pointer to the first element in the matrix data
 * for cases where direct memory copying is desired.
 *
 * @tparam T the element type
 * @tparam M the matrix type
 * @param matrix the input matrix
 * @return pointer to first element, or nullptr if matrix is empty
 */
template <class T, class M>
T*
findFirstElementPointer(const M& matrix)
{
    // We must find the pointer to the first element in the matrix
    // Iterate until we find an element, we can get the blocksize from the element
    // TODO: Can this be done more cleanly in the DUNE api to access the raw data more directly?
    T* nonZeroElementsTmp = nullptr;
    for (auto rowIt = matrix.begin(); rowIt != matrix.end(); ++rowIt) {
        auto colIt = rowIt->begin();
        if (colIt != rowIt->end()) {
            nonZeroElementsTmp = const_cast<T*>(&((*colIt)[0][0]));
            break;
        }
    }
    OPM_ERROR_IF(nonZeroElementsTmp == nullptr, "No non-zero elements found in matrix");
    return nonZeroElementsTmp;
}

/**
 * @brief Validate matrix conversion parameters
 *
 * Performs sanity checks on the extracted sparsity pattern to ensure
 * it's valid for cuSPARSE operations.
 *
 * @tparam M the matrix type
 * @param matrix the original matrix
 * @param columnIndices the extracted column indices
 * @param rowIndices the extracted row indices
 */
template <class M>
void
validateMatrixConversion(const M& matrix, const std::vector<int>& columnIndices, const std::vector<int>& rowIndices)
{
    const size_t numberOfRows = matrix.N();
    const size_t numberOfNonzeroBlocks = matrix.nonzeroes();

    OPM_ERROR_IF(rowIndices[matrix.N()] != to_int(matrix.nonzeroes()),
                 "Error: Row indices do not sum to number of nonzeroes");
    OPM_ERROR_IF(rowIndices.size() != numberOfRows + 1, "Error: Row indices size does not match matrix dimensions");
    OPM_ERROR_IF(columnIndices.size() != numberOfNonzeroBlocks,
                 "Error: Column indices size does not match number of nonzero blocks");
}

/**
 * @brief Validate that two matrices have compatible dimensions
 *
 * @param nonzeroes1 number of nonzeros in first matrix
 * @param blockSize1 block size of first matrix
 * @param N1 number of rows in first matrix
 * @param nonzeroes2 number of nonzeros in second matrix
 * @param blockSize2 block size of second matrix
 * @param N2 number of rows in second matrix
 */
inline void
validateMatrixCompatibility(
    size_t nonzeroes1, size_t blockSize1, size_t N1, size_t nonzeroes2, size_t blockSize2, size_t N2)
{
    OPM_ERROR_IF(nonzeroes1 != nonzeroes2, "Matrix does not have the same number of non-zero elements.");
    OPM_ERROR_IF(blockSize1 != blockSize2, "Matrix does not have the same blocksize.");
    OPM_ERROR_IF(N1 != N2, "Matrix does not have the same number of rows.");
}

/**
 * @brief Validate vector-matrix size compatibility for operations
 *
 * @param vectorSize size of the vector
 * @param matrixBlockSize block size of the matrix
 * @param matrixColumns number of matrix columns
 */
inline void
validateVectorMatrixSizes(size_t vectorSize, size_t matrixBlockSize, size_t matrixColumns)
{
    const size_t expectedSize = matrixBlockSize * matrixColumns;
    if (vectorSize != expectedSize) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Size mismatch. Input vector has {} elements, while matrix expects {} elements.",
                              vectorSize,
                              expectedSize));
    }
}

/**
 * @brief Macro for generating template instantiations for DUNE matrix operations
 *
 * This macro generates the necessary explicit template instantiations for
 * fromMatrix and updateNonzeroValues functions with DUNE matrix types.
 */
#define INSTANTIATE_SPARSE_MATRIX_DUNE_OPERATIONS(CLASS_NAME, realtype, blockdim)                                      \
    template CLASS_NAME<realtype> CLASS_NAME<realtype>::fromMatrix(                                                    \
        const Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>&, bool);                               \
    template CLASS_NAME<realtype> CLASS_NAME<realtype>::fromMatrix(                                                    \
        const Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>&, bool);                                \
    template void CLASS_NAME<realtype>::updateNonzeroValues(                                                           \
        const Dune::BCRSMatrix<Dune::FieldMatrix<realtype, blockdim, blockdim>>&, bool);                               \
    template void CLASS_NAME<realtype>::updateNonzeroValues(                                                           \
        const Dune::BCRSMatrix<Opm::MatrixBlock<realtype, blockdim, blockdim>>&, bool)

/**
 * @brief Macro for generating template instantiations for a range of block sizes
 */
#define INSTANTIATE_FOR_TYPE_AND_CLASS(CLASS_NAME, T)                                                                  \
    INSTANTIATE_SPARSE_MATRIX_DUNE_OPERATIONS(CLASS_NAME, T, 1);                                                       \
    INSTANTIATE_SPARSE_MATRIX_DUNE_OPERATIONS(CLASS_NAME, T, 2);                                                       \
    INSTANTIATE_SPARSE_MATRIX_DUNE_OPERATIONS(CLASS_NAME, T, 3);                                                       \
    INSTANTIATE_SPARSE_MATRIX_DUNE_OPERATIONS(CLASS_NAME, T, 4);                                                       \
    INSTANTIATE_SPARSE_MATRIX_DUNE_OPERATIONS(CLASS_NAME, T, 5);                                                       \
    INSTANTIATE_SPARSE_MATRIX_DUNE_OPERATIONS(CLASS_NAME, T, 6)

} // namespace Opm::gpuistl::detail

#endif
