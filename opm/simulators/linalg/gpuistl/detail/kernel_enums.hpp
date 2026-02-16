/*
  Copyright 2024 Equinor ASA

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

#ifndef OPM_GPUISTL_KERNEL_ENUMS_HPP
#define OPM_GPUISTL_KERNEL_ENUMS_HPP

#include <cuda_runtime.h>

/*
    This file organizes a growing amount of different mixed precision options for the preconditioners.
*/

namespace Opm::gpuistl {
    // Mixed precision schemes used for storing the matrix in GPU memory
    enum class MatrixStorageMPScheme {
        DOUBLE_DIAG_DOUBLE_OFFDIAG = 0, // full precision should be default
        FLOAT_DIAG_FLOAT_OFFDIAG = 1,
        DOUBLE_DIAG_FLOAT_OFFDIAG = 2
    };

    namespace detail {
        bool isValidMatrixStorageMPScheme(int scheme);
    }

    inline MatrixStorageMPScheme makeMatrixStorageMPScheme(int scheme) {
        if (!detail::isValidMatrixStorageMPScheme(scheme)) {
            OPM_THROW(std::invalid_argument,
                      fmt::format(fmt::runtime("Invalid matrix storage mixed precision scheme chosen: {}.\n"
                                  "Valid Schemes:\n"
                                  "\t0: DOUBLE_DIAG_DOUBLE_OFFDIAG\n"
                                  "\t1: FLOAT_DIAG_FLOAT_OFFDIAG\n"
                                  "\t2: DOUBLE_DIAG_FLOAT_OFFDIAG"),
                                  scheme));
        }
        return static_cast<MatrixStorageMPScheme>(scheme);
    }

    namespace detail {

        __host__ __device__ constexpr bool storeDiagonalAsFloat(MatrixStorageMPScheme scheme) {
            switch (scheme) {
                case MatrixStorageMPScheme::DOUBLE_DIAG_DOUBLE_OFFDIAG:
                    return false;
                case MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG:
                    return true;
                case MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG:
                    return false;
                default:
                    return false;
            }
        }

        __host__ __device__ constexpr bool storeOffDiagonalAsFloat(MatrixStorageMPScheme scheme) {
            switch (scheme) {
                case MatrixStorageMPScheme::DOUBLE_DIAG_DOUBLE_OFFDIAG:
                    return false;
                case MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG:
                    return true;
                case MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG:
                    return true;
                default:
                    return false;
            }
        }

        // returns true if we use anything else that the the default double precision for everything
        __host__ __device__ constexpr bool usingMixedPrecision(MatrixStorageMPScheme scheme) {
            switch (scheme) {
                case MatrixStorageMPScheme::DOUBLE_DIAG_DOUBLE_OFFDIAG:
                    return false;
                case MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG:
                    return true;
                case MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG:
                    return true;
                default:
                    return false;
            }
        }

        inline bool isValidMatrixStorageMPScheme(int scheme) {
            switch (static_cast<MatrixStorageMPScheme>(scheme)) {
                case MatrixStorageMPScheme::DOUBLE_DIAG_DOUBLE_OFFDIAG:
                case MatrixStorageMPScheme::FLOAT_DIAG_FLOAT_OFFDIAG:
                case MatrixStorageMPScheme::DOUBLE_DIAG_FLOAT_OFFDIAG:
                    return true;
                default:
                    return false;
            }
        }
    }
}

#endif // OPM_GPUISTL_KERNEL_ENUMS_HPP
