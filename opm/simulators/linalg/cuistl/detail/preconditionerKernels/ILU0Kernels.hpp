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
#ifndef OPM_ILU0_KERNELS_HPP
#define OPM_ILU0_KERNELS_HPP
#include <cstddef>
#include <vector>
namespace Opm::cuistl::detail::ILU0
{

/**
 * @brief Perform a upper solve on certain rows in a matrix that can safely be computed in parallel
 * @param reorderedMat pointer to GPU memory containing nonzerovalues of the sparse matrix. The matrix reordered such
 * that rows in the same level sets are contiguous
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param indexConversion Integer array containing mapping an index in the reordered matrix to its corresponding index
 * in the natural ordered matrix
 * @param startIdx Index of the first row of the matrix to be solve
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param [out] v Will store the results of the upper solve
 * @param threadBlockSize The number of threads per threadblock. Leave as -1 if no blocksize is already chosen
 */
template <class T, int blocksize>
void solveUpperLevelSet(T* reorderedMat,
                           int* rowIndices,
                           int* colIndices,
                           int* indexConversion,
                           int startIdx,
                           int rowsInLevelSet,
                           T* v,
                           int threadBlockSize);

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
 * @param d Stores the defect
 * @param [out] v Will store the results of the lower solve
 * @param threadBlockSize The number of threads per threadblock. Leave as -1 if no blocksize is already chosen
 */
template <class T, int blocksize>
void solveLowerLevelSet(T* reorderedMat,
                           int* rowIndices,
                           int* colIndices,
                           int* indexConversion,
                           int startIdx,
                           int rowsInLevelSet,
                           const T* d,
                           T* v,
                           int threadBlockSize);

/**
 * @brief Perform an upper solve on certain rows in a matrix that can safely be computed in parallel
 * @param reorderedUpperMat pointer to GPU memory containing nonzerovalues of the sparse matrix. The matrix reordered
 * such that rows in the same level sets are contiguous. This matrix is assumed to be strictly upper triangular
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param indexConversion Integer array containing mapping an index in the reordered matrix to its corresponding index
 * in the natural ordered matrix
 * @param startIdx Index of the first row of the matrix to be solve
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param dInv The diagonal elements of the LU factorization. Store the inverse of the diagonal entries
 * @param [out] v Will store the results of the lower solve. To begin with it should store the output from the lower
 * solve
 * @param threadBlockSize The number of threads per threadblock. Leave as -1 if no blocksize is already chosen
 */
template <class T, int blocksize>
void solveUpperLevelSetSplit(T* reorderedMat,
                                int* rowIndices,
                                int* colIndices,
                                int* indexConversion,
                                int startIdx,
                                int rowsInLevelSet,
                                const T* dInv,
                                T* v,
                                int threadBlockSize);

/**
 * @brief Perform an lower solve on certain rows in a matrix that can safely be computed in parallel
 * @param reorderedLowerMat pointer to GPU memory containing nonzerovalues of the sparse matrix. The matrix reordered
 * such that rows in the same level sets are contiguous. This matrix is assumed to be strictly lower triangular
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param indexConversion Integer array containing mapping an index in the reordered matrix to its corresponding index
 * in the natural ordered matrix
 * @param startIdx Index of the first row of the matrix to be solve
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param dInv The diagonal elements of the LU factorization. Store the inverse of the diagonal entries
 * @param d Stores the defect
 * @param [out] v Will store the results of the lower solve. To begin with it should store the output from the lower
 * solve
 * @param threadBlockSize The number of threads per threadblock. Leave as -1 if no blocksize is already chosen
 */
template <class T, int blocksize>
void solveLowerLevelSetSplit(T* reorderedLowerMat,
                                int* rowIndices,
                                int* colIndices,
                                int* indexConversion,
                                int startIdx,
                                int rowsInLevelSet,
                                const T* dInv,
                                const T* d,
                                T* v,
                                int threadBlockSize);

/**
 * @brief Computes the ILU Factorization of the input bcsr matrix, which is stored in a reordered way. The diagonal
 * elements store the inverse of the diagonal entries
 * @param [out] reorderedMat pointer to GPU memory containing nonzerovalues of the sparse matrix. The matrix reordered
 * such that rows in the same level sets are contiguous
 * @param rowIndices Pointer to vector on GPU containing row indices compliant wiht bsr format
 * @param colIndices Pointer to vector on GPU containing col indices compliant wiht bsr format
 * @param naturalToreordered Integer array containing mapping an index in the reordered matrix to its corresponding
 * index in the natural ordered matrix
 * @param reorderedToNatural Integer array containing mapping an index in the reordered matrix to its corresponding
 * index in the natural ordered matrix
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param startIdx Index of the first row of the matrix to be solve
 * @param threadBlockSize The number of threads per threadblock. Leave as -1 if no blocksize is already chosen
 */
template <class T, int blocksize>
void LUFactorization(T* reorderedMat,
                     int* rowIndices,
                     int* columnIndices,
                     int* naturalToReordered,
                     int* reorderedToNatual,
                     size_t rowsInLevelSet,
                     int startIdx,
                     int threadBlockSize);

/**
 * @brief Computes the ILU0 factorization in-place of a bcsr matrix stored in a split format (lower, diagonal and upper
 * triangular part)
 * @param [out] reorderedLowerMat pointer to GPU memory containing nonzerovalues of the strictly lower triangular sparse
 * matrix. The matrix reordered such that rows in the same level sets are contiguous
 * @param lowerRowIndices Pointer to vector on GPU containing row indices of the lower matrix compliant wiht bsr format
 * @param lowerColIndices Pointer to vector on GPU containing col indices of the lower matrix compliant wiht bsr format
 * @param [out] reorderedUpperMat pointer to GPU memory containing nonzerovalues of the strictly upper triangular sparse
 * matrix. The matrix reordered such that rows in the same level sets are contiguous
 * @param upperRowIndices Pointer to vector on GPU containing row indices of the upper matrix compliant wiht bsr format
 * @param upperColIndices Pointer to vector on GPU containing col indices of the upper matrix compliant wiht bsr format
 * @param [out] diagonal The diagonal elements of the ILU0 factorization inverted
 * @param reorderedToNatural Integer array containing mapping an index in the reordered matrix to its corresponding
 * index in the natural ordered matrix
 * @param naturalToreordered Integer array containing mapping an index in the reordered matrix to its corresponding
 * index in the natural ordered matrix
 * @param startIdx Index of the first row of the matrix to be solve
 * @param rowsInLevelSet Number of rows in this level set, which number the amount of rows solved in parallel by this
 * function
 * @param threadBlockSize The number of threads per threadblock. Leave as -1 if no blocksize is already chosen
 */
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
                          int startIdx,
                          int rowsInLevelSet,
                          int threadBlockSize);

} // namespace Opm::cuistl::detail::ILU0
#endif
