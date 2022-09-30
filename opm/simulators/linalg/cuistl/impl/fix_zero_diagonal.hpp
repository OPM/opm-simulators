/*
  Copyright SINTEF AS 2022

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

//#define CUISTL_ASSUME_NON_ZERO_DIAGONAL 1
namespace Opm::cuistl::impl
{
template <class Matrix>
std::vector<typename Matrix::field_type>
fixZeroDiagonal(const Matrix& matrix,
                const typename Matrix::field_type replacementValue
                = std::numeric_limits<typename Matrix::field_type>::epsilon())
{
    using field_type = typename Matrix::field_type;
    std::vector<field_type> nonZeroes(matrix.nonzeroes() * Matrix::block_type::cols * Matrix::block_type::cols, 0.0);

    const auto dataPointer = static_cast<const field_type*>(&(matrix[0][0][0][0]));
    std::copy(dataPointer, dataPointer + nonZeroes.size(), nonZeroes.begin());

    // TODO: Is there a neater way of accessing the underlying CRS structure?
    size_t currentNonzeroPointer = 0u;
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        for (auto column = row->begin(); column != row->end(); ++column) {
            if (column.index() == row.index()) {
                for (int component = 0; component < Matrix::block_type::cols; ++component) {
                    const auto index = currentNonzeroPointer + Matrix::block_type::cols * component + component;
                    if (nonZeroes[index] == 0) {
                        nonZeroes[index] = replacementValue;
                    }
                }
            }
            currentNonzeroPointer += 1;
        }
    }

    return nonZeroes;
}

#ifdef CUISTL_ASSUME_NON_ZERO_DIAGONAL
template <class Matrix>
const Matrix&
#else
template <class Matrix>
const Matrix&
#endif
makeMatrixWithNonzeroDiagonal(const Matrix& matrix,
                              const typename Matrix::field_type replacementValue
                              = std::numeric_limits<typename Matrix::field_type>::epsilon())
{
#ifdef CUISTL_ASSUME_NON_ZERO_DIAGONAL
    return matrix;
#else
    auto& newMatrix = const_cast<Matrix&>(matrix);
    // TODO: Is this fast enough?
    for (int row = 0; row < newMatrix.N(); ++row) {
        for (int component = 0; component < Matrix::block_type::cols; ++component) {
            if (newMatrix[row][row][component][component] == 0) {
                newMatrix[row][row][component][component] = replacementValue;
            }
        }
    }

    return newMatrix;
#endif
}
} // namespace Opm::cuistl::impl

#endif
