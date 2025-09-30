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

#ifndef OPM_IS_GPU_OPERATOR_HEADER
#define OPM_IS_GPU_OPERATOR_HEADER

#include <type_traits>

#if HAVE_CUDA // Avoid including GpuVector.hpp if CUDA 
              // is not enabled to avoid linking errors.
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#endif
namespace Opm {

/**
 * \brief Check if a given operator is a GPU operator.
 *
 * This is used to check if the operator is a GPU operator, which is used
 * in say the preconditioner factory to specialize the StandardPreconditioners class
 *
 * \tparam T The type of the operator to check.
 */
template <typename T>
struct is_gpu_operator
{
#if HAVE_CUDA
    // TODO: This can be done more thoroughly by checking if range and matrix also is a gpu operator, but this
    // works for now.
    static constexpr bool value
        = std::is_same_v<typename T::domain_type, Opm::gpuistl::GpuVector<typename T::field_type>>;
#else
    // If CUDA is not enabled, we assume that the operator is not a GPU operator.
    static constexpr bool value = false;
#endif
};

template <typename T>
static constexpr bool is_gpu_operator_v = is_gpu_operator<T>::value;


/**
 * \brief Check if a given operator is a GPU matrix.
 *
 * This is used to check if the matrix is a GPU matrix, which is used
 * in say the preconditioners to specialize
 *
 * \tparam T The type of the matrix to check.
 */
 template <typename T>
 struct is_gpu_matrix
 {
 #if HAVE_CUDA

     static constexpr bool value
         = std::is_same_v<T, Opm::gpuistl::GpuSparseMatrix<typename T::field_type>>
         || std::is_same_v<T, Opm::gpuistl::GpuSparseMatrixWrapper<typename T::field_type>>
         || std::is_same_v<T, Opm::gpuistl::GpuSparseMatrixGeneric<typename T::field_type>>;
 #else
     // If CUDA is not enabled, we assume that the matrix is not a GPU matrix.
     static constexpr bool value = false;
 #endif
 };
 
 template <typename T>
 static constexpr bool is_gpu_matrix_v = is_gpu_matrix<T>::value;

} // namespace Opm

#endif // OPM_IS_GPU_OPERATOR_HEADER
