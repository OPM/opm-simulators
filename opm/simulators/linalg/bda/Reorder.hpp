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

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

namespace bda
{

#define MAX_COLORS 256

/// Give every node in the matrix a color so that no neighbouring nodes share a color
/// The color array must be allocated already
/// This function with throw an error if no coloring can be found within the given restrictions
/// This function does graph coloring based on random numbers
/// \param[in] rows            number of rows in the matrix
/// \param[in] rowPointers     array of row pointers
/// \param[in] colIndices      array of column indices
/// \param[inout] colors       output array containing the number of the color that each row is assigned to
/// \param[in] maxRowsPerColor the maximum number of rows that are allowed in one color (so: the maximum number of nodes per color)
/// \param[in] maxColsPerColor the maximum number of columns that the rows in a color are allowed to share (so: the maximum number of nodes that the nodes in one color may be connected to)
/// \return                    the number of colors needed for the coloring
int colorBlockedNodes(int rows, const int *rowPointers, const int *colIndices, std::vector<int>& colors, int maxRowsPerColor, int maxColsPerColor);

/// Reorder the rows of the matrix according to the mapping in toOrder and fromOrder
/// rMat must be allocated already
/// \param[in] mat           matrix to be reordered
/// \param[in] toOrder       reorder pattern that lists for each index in the original order, to which index in the new order it should be moved
/// \param[in] fromOrder     reorder pattern that lists for each index in the new order, from which index in the original order it was moved
/// \param[inout] rMat       reordered Matrix 
template <unsigned int block_size>
void blocked_reorder_matrix_by_pattern(BlockedMatrix *mat, int *toOrder, int *fromOrder, BlockedMatrix *rMat);

/// Compute reorder mapping from the color that each node has received
/// The toOrder, fromOrder and iters arrays must be allocated already
/// \param[in] Nb            number of blocks in the vector
/// \param[in] colors        array containing the number of the color that each row is assigned to
/// \param[inout] toOrder    reorder pattern that lists for each index in the original order, to which index in the new order it should be moved
/// \param[inout] fromOrder  reorder pattern that lists for each index in the new order, from which index in the original order it was moved
/// \param[inout] iters      array containing for each color the number of rows that it contains
void colorsToReordering(int Nb, std::vector<int>& colors, int *toOrder, int *fromOrder, int *iters);

/// Reorder a vector according to the mapping in toOrder and fromOrder
/// The rVector array must be allocated already
/// \param[in] Nb            number of blocks in the vector
/// \param[in] vector        vector to be reordered
/// \param[in] toOrder       reorder pattern that lists for each index in the original order, to which index in the new order it should be moved
/// \param[in] fromOrder     reorder pattern that lists for each index in the new order, from which index in the original order it was moved
/// \param[inout] rVector    reordered vector
template <unsigned int block_size>
void blocked_reorder_vector_by_pattern(int Nb, double *vector, int *fromOrder, double *rVector);

/// Determine whether all rows that a certain row depends on are done already
/// \param[in] rowIndex      index of the row that needs to be checked for
/// \param[in] rowPointers   row pointers of the matrix that the row is in
/// \param[in] colIndices    column indices of the matrix that the row is in
/// \param[in] doneRows      array that for each row lists whether it is done or not
/// \return                  true iff all dependencies are done and if the result itself was not done yet
bool canBeStarted(int rowIndex, int *rowPointers, int *colIndices, std::vector<bool>& doneRows);

/// Find a level scheduling reordering for an input matrix
/// The toOrder and fromOrder arrays must be allocated already
/// \param[in] CSRColIndices  column indices array, obtained from storing the input matrix in the CSR format
/// \param[in] CSRRowPointers row pointers array, obtained from storing the input matrix in the CSR format
/// \param[in] CSCColIndices  row indices array, obtained from storing the input matrix in the CSC format
/// \param[in] CSCRowPointers column pointers array, obtained from storing the input matrix in the CSC format
/// \param[in] Nb             number of blockrows in the matrix
/// \param[out] iters         a pointer to the number of colors needed for the level scheduling
/// \param[inout] toOrder     the reorder pattern that was found, which lists for each index in the original order, to which index in the new order it should be moved
/// \param[inout] fromOrder   the reorder pattern that was found, which lists for each index in the new order, from which index in the original order it was moved
/// \return                   a pointer to an array that contains for each color, the number of rows that that color contains
int* findLevelScheduling(int *CSRColIndices, int *CSRRowPointers, int *CSCColIndices, int *CSCRowPointers, int Nb, int *iters, int *toOrder, int* fromOrder);

/// Find a graph coloring reordering for an input matrix
/// The toOrder and fromOrder arrays must be allocated already
/// \param[in] colIndices      column indices of the input matrix
/// \param[in] rowPointers     row pointers of the input matrix
/// \param[in] Nb              number of blockrows in the matrix
/// \param[in] maxRowsPerColor the maximum number of rows that are allowed in one color (so: the maximum number of nodes per color)
/// \param[in] maxColsPerColor the maximum number of columns that the rows in a color are allowed to share (so: the maximum number of nodes that the nodes in one color may be connected to)
/// \param[out] numColors      the number of colors used in the found graph coloring
/// \param[inout] toOrder      the reorder pattern that was found, which lists for each index in the original order, to which index in the new order it should be moved
/// \param[inout] fromOrder    the reorder pattern that was found, which lists for each index in the new order, from which index in the original order it was moved
/// \return                    a pointer to an array that contains for each color, the number of rows that that color contains
int* findGraphColoring(int *colIndices, int *rowPointers, int Nb, int maxRowsPerColor, int maxColsPerColor, int *numColors, int *toOrder, int* fromOrder);

/// Convert BCSR matrix to BCSC
/// Arrays for output matrix B must be allocated already
/// Based on the csr_tocsc() function from the scipy package from python, https://github.com/scipy/scipy/blob/master/scipy/sparse/sparsetools/csr.h
/// \param[in] Avals          non-zero values of the BCSR matrix
/// \param[in] Acols          column indices of the BCSR matrix
/// \param[in] Arows          row pointers of the BCSR matrix
/// \param[inout] Bvals       non-zero values of the result BCSC matrix
/// \param[inout] Bcols       row indices of the result BCSC matrix
/// \param[inout] Brows       column pointers of the result BCSC matrix
/// \param[in] Nb             number of blockrows in the matrix
template <unsigned int block_size>
void bcsr_to_bcsc(double *Avals, int *Acols, int *Arows, double *Bvals, int *Bcols, int *Brows, int Nb);

}

#endif