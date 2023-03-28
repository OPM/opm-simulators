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
#ifndef OPM_CUISTL_VECTOR_OPERATIONS_HPP
#define OPM_CUISTL_VECTOR_OPERATIONS_HPP
#include <cstddef>
namespace Opm::cuistl::detail
{

/**
 * @brief setVectorValue sets every element of deviceData to value
 * @param deviceData pointer to GPU memory
 * @param numberOfElements number of elements to set to value
 * @param value the value to use
 */
template <class T>
void setVectorValue(T* deviceData, size_t numberOfElements, const T& value);

/**
 * @brief setZeroAtIndexSet sets deviceData to zero in the indices of contained in indices
 * @param deviceData the data to operate on
 * @param numberOfElements number of indices
 * @param indices the indices to use
 */
template <class T>
void setZeroAtIndexSet(T* deviceData, size_t numberOfElements, const int* indices);

/**
 * @brief innerProductAtIndices computes the inner product between deviceA[indices] and deviceB[indices]
 * @param deviceA data A
 * @param deviceB data B
 * @param buffer a buffer with number of elements equal to numberOfElements
 * @param numberOfElements number of indices
 * @param indices the indices to compute the inner product over
 * @return the result of the inner product
 *
 * @note This is equivalent to projecting the vectors to the indices contained in indices, then doing the inner product
 * of those indices.
 */
template <class T>
T innerProductAtIndices(const T* deviceA, const T* deviceB, T* buffer, size_t numberOfElements, const int* indices);
} // namespace Opm::cuistl::detail
#endif
