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
#ifndef CUDA_SAFE_CALL_HPP
#define CUDA_SAFE_CALL_HPP

#include <opm/common/ErrorMacros.hpp>

#define OPM_CUDA_SAFE_CALL(expression)                                                                                 \
    do {                                                                                                               \
        cudaError_t error = expression;                                                                                \
        if (error != cudaSuccess) {                                                                                    \
            OPM_THROW(std::runtime_error,                                                                              \
                      "CUDA expression did not execute correctly. Expression was: \n"                                  \
                          << "    " << #expression << "\n"                                                             \
                          << "CUDA error was " << cudaGetErrorString(error) << "\n"                                    \
                          << "in function " << __func__ << ", in " << __FILE__ << " at line " << __LINE__ << "\n");    \
        }                                                                                                              \
    } while (false)
#endif
