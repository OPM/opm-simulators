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
#ifndef OPM_CUDA_CHECK_LAST_ERROR_HPP
#define OPM_CUDA_CHECK_LAST_ERROR_HPP
#include <cuda_runtime.h>
#include <fmt/core.h>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>

/**
 * @brief OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE checks the return type of cudaDeviceSynchronize(),
 * and throws an exception if cudaDeviceSynchronize() does not equal cudaSuccess.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cuda_check_last_error.hpp>
 *
 * void some_function() {
 *     OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE;
 * }
 * @endcode
 *
 * @note This can be used to debug the code, or simply make sure that no error has occured.
 * @note This is a rather heavy operation, so prefer to use only in Debug mode (see OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE_IF_DEBUG)
 */
#define OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE OPM_CUDA_SAFE_CALL(cudaDeviceSynchronize())

#ifdef NDEBUG
#define OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE_IF_DEBUG
#else

/**
 * @brief OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE_IF_DEBUG checks the return type of cudaDeviceSynchronize only if NDEBUG is not defined,
 * and throws an exception if cudaDeviceSynchronize() does not equal cudaSuccess.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
 *
 * void some_function() {
 *     OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE_IF_DEBUG;
 * }
 * @endcode
 *
 * @note This can be used to debug the code, or simply make sure that no error has occured.
 */
#define OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE_IF_DEBUG OPM_CUDA_CHECK_DEVICE_SYNCHRONIZE
#endif


/**
 * @brief OPM_CUDA_CHECK_LAST_ERROR checks the return type of cudaGetLastError(),
 * and throws an exception if cudaGetLastError() does not equal cudaSuccess.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cuda_check_last_error.hpp>
 *
 * void some_function() {
 *     OPM_CUDA_CHECK_LAST_ERROR;
 * }
 * @endcode
 *
 * @note This can be used to debug the code, or simply make sure that no error has occured.
 */
#define OPM_CUDA_CHECK_LAST_ERROR OPM_CUDA_SAFE_CALL(cudaGetLastError())

#ifdef NDEBUG
#define OPM_CUDA_CHECK_LAST_ERROR_IF_DEBUG
#else

/**
 * @brief OPM_CUDA_CHECK_LAST_ERROR_IF_DEBUG checks the return type of cudaGetLastError() only if NDEBUG is not defined,
 * and throws an exception if cudaGetLastError() does not equal cudaSuccess.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cuda_check_last_error.hpp>
 *
 * void some_function() {
 *     OPM_CUDA_CHECK_LAST_ERROR_IF_DEBUG;
 * }
 * @endcode
 *
 * @note This can be used to debug the code, or simply make sure that no error has occured.
 */
#define OPM_CUDA_CHECK_LAST_ERROR_IF_DEBUG OPM_CUDA_CHECK_LAST_ERROR
#endif

#endif
