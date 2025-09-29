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
#ifndef OPM_GPU_CPR_AMG_OPERATIONS_HPP
#define OPM_GPU_CPR_AMG_OPERATIONS_HPP

#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <cstddef>

namespace Opm::gpuistl::detail
{

/**
 * @brief Calculates quasi-IMPES weights for CPR preconditioner on GPU
 * @param matrix The matrix used for weight calculation
 * @param pressureVarIndex The index of the pressure variable
 * @param weights Output vector to store the calculated weights
 * @param diagonalIndices The diagonal indices of the matrix
 * @tparam transpose Whether to use the transpose mode for calculation
 */
template <typename T, bool transpose>
void getQuasiImpesWeights(const GpuSparseMatrixWrapper<T>& matrix,
                          std::size_t pressureVarIndex,
                          GpuVector<T>& weights,
                          const GpuVector<int>& diagonalIndices);

/**
 * @brief Calculates the coarse level matrix entries based on the fine level matrix and weights
 * @param fineMatrix The fine level matrix
 * @param coarseMatrix The coarse level matrix to fill
 * @param weights The weights vector for pressure calculation
 * @param pressureVarIndex The index of the pressure variable
 * @tparam transpose Whether to use the transpose representation or not
 */
template <typename T, bool transpose>
void calculateCoarseEntries(const GpuSparseMatrixWrapper<T>& fineMatrix,
                            GpuSparseMatrixWrapper<T>& coarseMatrix,
                            const GpuVector<T>& weights,
                            std::size_t pressureVarIndex);


/**
 * @brief Restricts a fine level vector to a coarse level vector based on pressure index
 * @param fine The fine level vector
 * @param coarse The coarse level vector to fill
 * @param weights The weights vector for pressure calculation
 * @param pressureVarIndex The index of the pressure variable
 * @tparam transpose Whether to use the transpose representation or not
 */
template <typename T, bool transpose>
void restrictVector(const GpuVector<T>& fine,
                    GpuVector<T>& coarse,
                    const GpuVector<T>& weights,
                    std::size_t pressureVarIndex);

/**
 * @brief Prolongs a coarse level vector to a fine level vector based on pressure index
 * @param coarse The coarse level vector
 * @param fine The fine level vector to fill
 * @param weights The weights vector for pressure calculation
 * @param pressureVarIndex The index of the pressure variable
 * @tparam transpose Whether to use the transpose representation or not
 */
template <typename T, bool transpose>
void prolongateVector(const GpuVector<T>& coarse,
                      GpuVector<T>& fine,
                      const GpuVector<T>& weights,
                      std::size_t pressureVarIndex);

} // namespace Opm::gpuistl::detail

#endif // OPM_GPU_CPR_AMG_OPERATIONS_HPP
