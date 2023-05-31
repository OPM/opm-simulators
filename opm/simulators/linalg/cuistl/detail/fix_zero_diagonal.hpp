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
#ifndef OPM_FIXZERODIAGONAL_HEADER_INCLUDED
#define OPM_FIXZERODIAGONAL_HEADER_INCLUDED

#include <limits>
#include <vector>

namespace Opm::cuistl::detail
{

/**
 * @brief makeMatrixWithNonzeroDiagonal creates a new matrix with the zero diagonal elements (when viewed as a matrix of
 * scalrars) set to replacementValue
 * @param matrix the matrix to replace
 * @param replacementValue the value to set in the diagonal elements that are zero
 * @return a new matrix with non non-zero diagonal elements.
 *
 * @note This modification is needed for the CuSparse implementation of ILU0. While the the algorithm operates on block
 * matrices, it still requires that the scalar matrix has no zero diagonal elements.
 */
template <class Matrix>
const Matrix
makeMatrixWithNonzeroDiagonal(const Matrix& matrix,
                              const typename Matrix::field_type replacementValue
                              = std::numeric_limits<typename Matrix::field_type>::epsilon())
{
    auto newMatrix = matrix;
    // TODO: [perf] Is this fast enough?
    for (size_t row = 0; row < newMatrix.N(); ++row) {
        for (size_t component = 0; component < Matrix::block_type::cols; ++component) {
            if (newMatrix[row][row][component][component] == 0) {
                newMatrix[row][row][component][component] = replacementValue;
            }
        }
    }

    return newMatrix;
}
} // namespace Opm::cuistl::detail

#endif
