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
#ifndef OPM_JAC_KERNELS_HPP
#define OPM_JAC_KERNELS_HPP
#include <cstddef>
#include <vector>
namespace Opm::gpuistl::detail::JAC
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
void invertDiagonalAndFlatten(const T* mat, const int* rowIndices, const int* colIndices, size_t numberOfRows, T* vec);

} // namespace Opm::gpuistl::detail::JAC
#endif
