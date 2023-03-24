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
#ifndef CUDA_SAFE_CALL_HPP
#define CUDA_SAFE_CALL_HPP
#include <cuda_runtime.h>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>

/**
 * @brief OPM_CUDA_SAFE_CALL checks the return type of the CUDA expression (function call) and throws an exception if it
 * does not equal cudaSuccess.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
 * #include <cuda_runtime.h>
 *
 * void some_function() {
 *     void* somePointer;
 *     OPM_CUDA_SAFE_CALL(cudaMalloc(&somePointer, 1));
 * }
 * @endcode
 *
 * @note This should be used for any call to cuSparse unless you have a good reason not to.
 */
#define OPM_CUDA_SAFE_CALL(expression)                                                                                 \
    do {                                                                                                               \
        cudaError_t error = expression;                                                                                \
        if (error != cudaSuccess) {                                                                                    \
            OPM_THROW(std::runtime_error,                                                                              \
                      fmt::format("CUDA expression did not execute correctly. Expression was: \n"                      \
                                  "    {}\n"                                                                           \
                                  "CUDA error was {}\n"                                                                \
                                  "in function {}, in {}, at line {}\n",                                               \
                                  #expression,                                                                         \
                                  cudaGetErrorString(error),                                                           \
                                  __func__,                                                                            \
                                  __FILE__,                                                                            \
                                  __LINE__));                                                                          \
        }                                                                                                              \
    } while (false)
#endif
