/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_CUISTL_CUSPARSE_MATRIX_OPERATIONS_HPP
#define OPM_CUISTL_CUSPARSE_MATRIX_OPERATIONS_HPP
#include <cstddef>
#include <vector>
namespace Opm::cuistl::detail
{

/**
 * @brief This function receives a matrix, and the inverse of the matrix containing only its diagonal is stored in d_vec
 * @param mat pointer to GPU memory containing nonzerovalues of the sparse matrix
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param numberOfRows Integer describing the number of rows in the matrix
 * @param[out] vec Pointer to the vector where the inverse of the diagonal matrix should be stored
 */
template <class T, int blocksize>
void invertDiagonalAndFlatten(T* mat, int* rowIndices, int* colIndices, size_t numberOfRows, T* vec);

/**
 * @brief Perform a lower solve on certain rows in a matrix that can safely be computed in parallel
 * @param reorderedMat pointer to GPU memory containing nonzerovalues of the sparse matrix. The matrix reordered such
 * that rows in the same level sets are contiguous
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param indexConversion Integer array containing mapping an index in the reordered matrix to its corresponding index
 * in the natural ordered matrix
 * @param startIdx Index of the first row of the matrix to be solve
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param dInv The diagonal matrix used by the Diagonal ILU preconditioner. Must be reordered in the same way as
 * reorderedMat
 * @param d Stores the defect
 * @param [out] v Will store the results of the lower solve
 */
template <class T, int blocksize>
void computeLowerSolveLevelSet(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               const T* d,
                               T* v);

/**
 * @brief Perform an upper solve on certain rows in a matrix that can safely be computed in parallel
 * @param reorderedMat pointer to GPU memory containing nonzerovalues of the sparse matrix. The matrix reordered such
 * that rows in the same level sets are contiguous
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param indexConversion Integer array containing mapping an index in the reordered matrix to its corresponding index
 * in the natural ordered matrix
 * @param startIdx Index of the first row of the matrix to be solve
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param dInv The diagonal matrix used by the Diagonal ILU preconditioner
 * @param [out] v Will store the results of the lower solve. To begin with it should store the output from the lower
 * solve
 */
template <class T, int blocksize>
void computeUpperSolveLevelSet(T* reorderedMat,
                               int* rowIndices,
                               int* colIndices,
                               int* indexConversion,
                               int startIdx,
                               int rowsInLevelSet,
                               const T* dInv,
                               T* v);

/**
 * @brief Computes the ILU0 of the diagonal elements of the reordered matrix and stores it in a reordered vector
 * containing the diagonal blocks
 * @param reorderedMat pointer to GPU memory containing nonzerovalues of the sparse matrix. The matrix reordered such
 * that rows in the same level sets are contiguous
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param reorderedToNatural Integer array containing mapping an index in the reordered matrix to its corresponding
 * index in the natural ordered matrix
 * @param naturalToreordered Integer array containing mapping an index in the reordered matrix to its corresponding
 * index in the natural ordered matrix
 * @param startIdx Index of the first row of the matrix to be solve
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param [out] dInv The diagonal matrix used by the Diagonal ILU preconditioner
 */
template <class T, int blocksize>
void computeDiluDiagonal(T* reorderedMat,
                         int* rowIndices,
                         int* colIndices,
                         int* reorderedToNatural,
                         int* naturalToReordered,
                         int startIdx,
                         int rowsInLevelSet,
                         T* dInv);

/**
 * @brief Reorders the elements of a matrix by copying them from one matrix to another using a permutation list
 * @param srcMatrix The source matrix we will copy data from
 * @param srcRowIndices Pointer to vector on GPU containing row indices for the source matrix compliant wiht bsr format
 * @param [out] dstMatrix The destination matrix that we copy data to
 * @param dstRowIndices Pointer to vector on GPU containing riw indices for the destination matrix compliant wiht bsr
 * format
 * @param naturalToReordered Permuation list that converts indices in the src matrix to the indices in the dst matrix
 * @param numberOfRows The number of rows in the matrices
 */
template <class T, int blocksize>
void copyMatDataToReordered(
    T* srcMatrix, int* srcRowIndices, T* dstMatrix, int* dstRowIndices, int* naturalToReordered, size_t numberOfRows);
} // namespace Opm::cuistl::detail
#endif
