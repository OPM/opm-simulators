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

#ifndef ILU_VARIANTS_HELPER_KERNELS_HPP
#define ILU_VARIANTS_HELPER_KERNELS_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include <opm/simulators/linalg/gpuistl/detail/gpuThreadUtils.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

namespace Opm::gpuistl::detail
{

/*
  The computeDiagIndices(NoReorder) functions are used to speed up the factorization time by
  not recomputing indices during each update kernel.
*/
/**
 * @brief Computes indices of diagonal elements for non-reordered GPU preconditioner.
 *
 * @tparam T Matrix element type
 * @param[in] rowIndices Row indices of the BCSR matrix
 * @param[in] colIndices Column indices of the BCSR matrix
 * @param[in] indexConversion Mapping from levelset-row index to matrix row
 * @param[in] rows Number of rows in the matrix
 * @param[out] diagIndices Array to store the indices of diagonal elements
 */
template <class T>
void computeDiagIndicesNoReorder(const int* rowIndices,
                                 const int* colIndices,
                                 const size_t* indexConversion,
                                 int rows,
                                 size_t* diagIndices);

/**
 * @brief Computes indices of diagonal elements for reordered GPU preconditioner.
 *
 * @tparam T Matrix element type
 * @param[in] rowIndices Row indices of the BCSR matrix
 * @param[in] colIndices Column indices of the BCSR matrix
 * @param[in] reorderedToNatural Table for reordering to natural numbering
 * @param[in] rows Number of rows in the matrix
 * @param[out] diagIndices Array to store the indices of diagonal elements
 */
template <class T>
void computeDiagIndices(const int* rowIndices,
                        const int* colIndices,
                        const int* reorderedToNatural,
                        int rows,
                        size_t* diagIndices);

} // namespace Opm::gpuistl::detail

#endif
