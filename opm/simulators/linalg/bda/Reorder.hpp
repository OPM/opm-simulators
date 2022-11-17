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

namespace Opm
{
namespace Accelerator
{

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