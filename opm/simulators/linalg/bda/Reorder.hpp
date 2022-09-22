/*
  Copyright 2019 Equinor ASA

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

#ifndef REORDER_HPP
#define REORDER_HPP

#include <vector>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

namespace Opm
{
namespace Accelerator
{

#define MAX_COLORS 256

/// Give every node in the matrix a color so that no neighbouring nodes share a color
/// The color array must be allocated already
/// This function with throw an error if no coloring can be found within the given restrictions
/// This function does graph coloring based on random numbers
/// \param[in] rows            number of rows in the matrix
/// \param[in] CSRRowPointers  array of row pointers of the sparsity pattern stored in the CSR format
/// \param[in] CSRColIndices   array of column indices of the sparsity pattern stored in the CSR format
/// \param[in] CSCColPointers  array of column pointers of the sparsity pattern stored in the CSC format
/// \param[in] CSCRowIndices   array of row indices of the sparsity pattern stored in the CSC format
/// \param[inout] colors       output array containing the number of the color that each row is assigned to
/// \param[in] maxRowsPerColor the maximum number of rows that are allowed in one color (so: the maximum number of nodes per color)
/// \param[in] maxColsPerColor the maximum number of columns that the rows in a color are allowed to share (so: the maximum number of nodes that the nodes in one color may be connected to)
/// \return                    the number of colors needed for the coloring
template <unsigned int block_size>
int colorBlockedNodes(int rows, const int *CSRRowPointers, const int *CSRColIndices, const int *CSCColPointers, const int *CSCRowIndices, std::vector<int>& colors, int maxRowsPerColor, int maxColsPerColor);

/// Reorder the sparsity pattern of the matrix according to the mapping in toOrder and fromOrder
/// Also find mapping for nnz blocks
/// rmat must be allocated already
/// \param[in] mat                        matrix to be reordered
/// \param[out] reordermapping_nonzeroes  contains new index for every nnz block
/// \param[in] toOrder                    reorder pattern that lists for each index in the original order, to which index in the new order it should be moved
/// \param[in] fromOrder                  reorder pattern that lists for each index in the new order, from which index in the original order it was moved
/// \param[out] rmat                      reordered Matrix
void reorderBlockedMatrixByPattern(BlockedMatrix *mat, std::vector<int>& reordermapping_nonzeroes, int *toOrder, int *fromOrder, BlockedMatrix *rmat);

/// Write nnz blocks from mat to rmat, according to the mapping in reordermapping_nonzeroes
/// rmat must be allocated already
/// \param[in] mat                       matrix to be reordered
/// \param[in] reordermapping_nonzeroes  contains old index for every nnz block, so rmat_nnz[i] == mat_nnz[mapping[i]]
/// \param[inout] rmat                   reordered Matrix
void reorderNonzeroes(BlockedMatrix *mat, std::vector<int>& reordermapping_nonzeroes, BlockedMatrix *rmat);

/// Compute reorder mapping from the color that each node has received
/// The toOrder, fromOrder and iters arrays must be allocated already
/// \param[in] Nb              number of blocks in the vector
/// \param[in] colors          array containing the number of the color that each row is assigned to
/// \param[in] numColors       the total number of colors into which all rows have been divided
/// \param[inout] toOrder      reorder pattern that lists for each index in the original order, to which index in the new order it should be moved
/// \param[inout] fromOrder    reorder pattern that lists for each index in the new order, from which index in the original order it was moved
/// \param[inout] rowsPerColor array containing for each color the number of rows that it contains
void colorsToReordering(int Nb, std::vector<int>& colors, int numColors, int *toOrder, int *fromOrder, std::vector<int>& rowsPerColor);

/// Reorder a vector according to the mapping in fromOrder
/// The rVector array must be allocated already
/// \param[in] Nb            number of blocks in the vector
/// \param[in] vector        vector to be reordered
/// \param[in] toOrder       reorder pattern that lists for each index in the original order, to which index in the new order it should be moved
/// \param[in] fromOrder     reorder pattern that lists for each index in the new order, from which index in the original order it was moved
/// \param[inout] rVector    reordered vector
template <unsigned int block_size>
void reorderBlockedVectorByPattern(int Nb, double *vector, int *fromOrder, double *rVector);

/// Determine whether all rows that a certain row depends on are done already
/// \param[in] rowIndex      index of the row that needs to be checked for
/// \param[in] rowPointers   row pointers of the matrix that the row is in
/// \param[in] colIndices    column indices of the matrix that the row is in
/// \param[in] doneRows      array that for each row lists whether it is done or not
/// \return                  true iff all dependencies are done and if the result itself was not done yet
bool canBeStarted(const int rowIndex, const  int *rowPointers, const  int *colIndices, const std::vector<bool>& doneRows);

/// Find a level scheduling reordering for an input matrix
/// The toOrder and fromOrder arrays must be allocated already
/// \param[in] CSRColIndices    column indices array, obtained from storing the input matrix in the CSR format
/// \param[in] CSRRowPointers   row pointers array, obtained from storing the input matrix in the CSR format
/// \param[in] CSCRowIndices    row indices array, obtained from storing the input matrix in the CSC format
/// \param[in] CSCColPointers   column pointers array, obtained from storing the input matrix in the CSC format
/// \param[in] Nb               number of blockrows in the matrix
/// \param[out] numColors       a pointer to the number of colors needed for the level scheduling
/// \param[out] toOrder         the reorder pattern that was found, which lists for each index in the original order, to which index in the new order it should be moved
/// \param[out] fromOrder       the reorder pattern that was found, which lists for each index in the new order, from which index in the original order it was moved
/// \param[out] rowsPerColor    for each color, an array of all rowIndices in that color, this function uses emplace_back() to fill
void findLevelScheduling(int *CSRColIndices, int *CSRRowPointers, int *CSCRowIndices, int *CSCColPointers, int Nb, int *numColors, int *toOrder, int* fromOrder, std::vector<int>& rowsPerColor);

/// Find a graph coloring reordering for an input matrix
/// The toOrder and fromOrder arrays must be allocated already
/// \param[in] CSRColIndices    column indices of the input sparsity pattern stored in the CSR format
/// \param[in] CSRRowPointers   row pointers of the input sparsity pattern stored in the CSR format
/// \param[in] CSCRowIndices    row indices of the input sparsity pattern stored in the CSC format
/// \param[in] CSCColPointers   column pointers of the input sparsity pattern stored in the CSC format
/// \param[in] Nb               number of blockrows in the matrix
/// \param[in] maxRowsPerColor  the maximum number of rows that are allowed in one color (so: the maximum number of nodes per color)
/// \param[in] maxColsPerColor  the maximum number of columns that the rows in a color are allowed to share (so: the maximum number of nodes that the nodes in one color may be connected to)
/// \param[out] numColors       the number of colors used in the found graph coloring
/// \param[inout] toOrder       the reorder pattern that was found, which lists for each index in the original order, to which index in the new order it should be moved
/// \param[inout] fromOrder     the reorder pattern that was found, which lists for each index in the new order, from which index in the original order it was moved
/// \param[inout] rowsPerColor  for each used color, the number of rows assigned to that color, this function will resize()
template <unsigned int block_size>
void findGraphColoring(const int *CSRColIndices, const int *CSRRowPointers, const int *CSCRowIndices, const int *CSCColPointers, int Nb, int maxRowsPerColor, int maxColsPerColor, int *numColors, int *toOrder, int *fromOrder, std::vector<int>& rowsPerColor);

/// Convert a sparsity pattern stored in the CSR format to the CSC format
/// CSCRowIndices and CSCColPointers arrays must be allocated already
/// Based on the csr_tocsc() function from the scipy package from python, https://github.com/scipy/scipy/blob/master/scipy/sparse/sparsetools/csr.h
/// \param[in] CSRColIndices     column indices of the CSR representation of the pattern
/// \param[in] CSRRowPointers    row pointers of the CSR representation of the pattern
/// \param[inout] CSCRowIndices  row indices of the result CSC representation of the pattern
/// \param[inout] CSCColPointers column pointers of the result CSC representation of the pattern
/// \param[in] Nb                number of blockrows in the matrix
void csrPatternToCsc(int *CSRColIndices, int *CSRRowPointers, int *CSCRowIndices, int *CSCColPointers, int Nb);

} // namespace Accelerator
} // namespace Opm

#endif