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

#ifndef OPM_GPUISTL_GPU_TYPE_DETECTION_HPP
#define OPM_GPUISTL_GPU_TYPE_DETECTION_HPP

#include <type_traits>

namespace Opm::gpuistl
{

// Forward declarations
template <typename T>
class GpuVector;
template <typename T>
class GpuSparseMatrixWrapper;
template <typename T>
class GpuSparseMatrixGeneric;

/**
 * @brief Type trait to detect if a type is a GPU type
 */
template <typename T>
struct is_gpu_type : std::false_type {
};

// Specializations for known GPU types
template <typename T>
struct is_gpu_type<GpuVector<T>> : std::true_type {
};

template <typename T>
struct is_gpu_type<GpuSparseMatrixWrapper<T>> : std::true_type {
};

template <typename T>
struct is_gpu_type<GpuSparseMatrixGeneric<T>> : std::true_type {
};

/**
 * @brief Helper variable template for easier usage
 */
template <typename T>
inline constexpr bool is_gpu_type_v = is_gpu_type<T>::value;

} // namespace Opm::gpuistl

#endif // OPM_GPUISTL_GPU_TYPE_DETECTION_HPP
