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
#ifndef OPM_CUSPARSE_SAFE_CALL_HPP
#define OPM_CUSPARSE_SAFE_CALL_HPP
#include <cusparse.h>
#include <exception>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>

namespace Opm::cuistl::detail
{

#define CHECK_CUSPARSE_ERROR_TYPE(code, x)                                                                             \
    if (code == x) {                                                                                                   \
        return #x;                                                                                                     \
    }
/**
 * @brief getCusparseErrorMessage Converts an error code returned from a cusparse function a human readable string.
 * @param code an error code from a cusparse routine
 * @return a human readable string.
 */
inline std::string
getCusparseErrorMessage(int code)
{
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_SUCCESS);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_NOT_INITIALIZED);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_ALLOC_FAILED);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_INVALID_VALUE);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_ARCH_MISMATCH);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_MAPPING_ERROR);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_EXECUTION_FAILED);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_INTERNAL_ERROR);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_ZERO_PIVOT);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_NOT_SUPPORTED);
    CHECK_CUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_INSUFFICIENT_RESOURCES);
    return fmt::format("UNKNOWN CUSPARSE ERROR {}.", code);
}

#undef CHECK_CUSPARSE_ERROR_TYPE

/**
 * @brief cusparseSafeCall checks the return type of the CUSPARSE expression (function call) and throws an exception if
 * it does not equal CUSPARSE_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
 * #include <cublas_v2.h>
 *
 * void some_function() {
 *     cublasHandle_t cublasHandle;
 *     cudaSafeCall(cusparseCreate(&cusparseHandle), "cusparseCreate(&cusparseHandle)", __FILE__, __func__, __LINE__);
 * }
 * @endcode
 *
 * @note It is probably easier to use the macro OPM_CUBLAS_SAFE_CALL
 *
 * @todo Refactor to use std::source_location once we shift to C++20
 */
inline void
cusparseSafeCall(cusparseStatus_t error,
                 const std::string_view& expression,
                 const std::string_view& filename,
                 const std::string_view& functionName,
                 size_t lineNumber)
{
    if (error != CUSPARSE_STATUS_SUCCESS) {
        OPM_THROW(std::runtime_error,
                  fmt::format("cuSparse expression did not execute correctly. Expression was: \n\n"
                              "    {}\n\nin function {}, in {}, at line {}\n"
                              "CuSparse error code was: {}\n",
                              expression,
                              functionName,
                              filename,
                              lineNumber,
                              getCusparseErrorMessage(error)));
    }
}
} // namespace Opm::cuistl::detail



/**
 * @brief OPM_CUSPARSE_SAFE_CALL checks the return type of the cusparse expression (function call) and throws an
 * exception if it does not equal CUSPARSE_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
 * #include <cusparse.h>
 *
 * void some_function() {
 *     cusparseHandle_t cusparseHandle;
 *     OPM_CUSPARSE_SAFE_CALL(cusparseCreate(&cusparseHandle));
 * }
 * @endcode
 *
 * @note This should be used for any call to cuSparse unless you have a good reason not to.
 */
#define OPM_CUSPARSE_SAFE_CALL(expression)                                                                             \
    ::Opm::cuistl::detail::cusparseSafeCall(expression, #expression, __FILE__, __func__, __LINE__)
#endif // OPM_CUSPARSE_SAFE_CALL_HPP
