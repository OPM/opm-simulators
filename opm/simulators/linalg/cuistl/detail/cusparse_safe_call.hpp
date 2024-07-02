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
#ifndef OPM_GPUSPARSE_SAFE_CALL_HPP
#define OPM_GPUSPARSE_SAFE_CALL_HPP
#include <cusparse.h>
#include <exception>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

namespace Opm::gpuistl::detail
{

#define CHECK_GPUSPARSE_ERROR_TYPE(code, x)                                                                             \
    if (code == x) {                                                                                                   \
        return #x;                                                                                                     \
    }
/**
 * @brief getGpuSPARSEErrorCodeToString Converts an error code returned from a cu/hipSPARSE function a human readable string.
 * @param code an error code from a cu/hipSPARSE routine
 * @return a human readable string.
 */
inline std::string
getGpuSPARSEErrorCodeToString(int code)
{
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_SUCCESS);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_NOT_INITIALIZED);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_ALLOC_FAILED);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_INVALID_VALUE);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_ARCH_MISMATCH);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_MAPPING_ERROR);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_EXECUTION_FAILED);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_INTERNAL_ERROR);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_ZERO_PIVOT);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_NOT_SUPPORTED);
    CHECK_GPUSPARSE_ERROR_TYPE(code, CUSPARSE_STATUS_INSUFFICIENT_RESOURCES);
    return fmt::format("UNKNOWN CUSPARSE ERROR {}.", code);
}

#undef CHECK_GPUSPARSE_ERROR_TYPE
/**
 * @brief getGpuSPARSEErrorMessage generates the error message to display for a given error.
 *
 * @param error the error code from cublas
 * @param expression the expresison (say "cusparseCreate(&handle)")
 * @param filename the code file the error occured in (typically __FILE__)
 * @param functionName name of the function the error occured in (typically __func__)
 * @param lineNumber the line number the error occured in (typically __LINE__)
 *
 * @todo Refactor to use std::source_location once we shift to C++20
 *
 * @return An error message to be displayed.
 *
 * @note This function is mostly for internal use.
 */
inline std::string
getGpuSPARSEErrorMessage(cusparseStatus_t error,
                        const std::string_view& expression,
                        const std::string_view& filename,
                        const std::string_view& functionName,
                        size_t lineNumber)
{
    return fmt::format("cuSparse expression did not execute correctly. Expression was: \n\n"
                       "    {}\n\nin function {}, in {}, at line {}\n"
                       "CuSparse error code was: {}\n",
                       expression,
                       functionName,
                       filename,
                       lineNumber,
                       getGpuSPARSEErrorCodeToString(error));
}

/**
 * @brief gpuSPARSESafeCall checks the return type of the CUSPARSE expression (function call) and throws an exception if
 * it does not equal CUSPARSE_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
 * #include <cublas_v2.h>
 *
 * void some_function() {
 *     cusparseHandle_t cusparseHandle;
 *     gpuSPARSESafeCall(cusparseCreate(&cusparseHandle), "cusparseCreate(&cusparseHandle)", __FILE__, __func__,
 * __LINE__);
 * }
 * @endcode
 *
 * @note It is probably easier to use the macro OPM_CUBLAS_SAFE_CALL
 *
 * @todo Refactor to use std::source_location once we shift to C++20
 */
inline void
gpuSPARSESafeCall(cusparseStatus_t error,
                 const std::string_view& expression,
                 const std::string_view& filename,
                 const std::string_view& functionName,
                 size_t lineNumber)
{
    if (error != CUSPARSE_STATUS_SUCCESS) {
        OPM_THROW(std::runtime_error, getGpuSPARSEErrorMessage(error, expression, filename, functionName, lineNumber));
    }
}

/**
 * @brief cusparseWarnIfError checks the return type of the CUSPARSE expression (function call) and issues a warning if
 * it does not equal CUSPARSE_STATUS_SUCCESS.
 *
 * @param error the error code from cublas
 * @param expression the expresison (say "cublasCreate(&handle)")
 * @param filename the code file the error occured in (typically __FILE__)
 * @param functionName name of the function the error occured in (typically __func__)
 * @param lineNumber the line number the error occured in (typically __LINE__)
 *
 * @return the error sent in (for convenience).
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
 * #include <cublas_v2.h>
 *
 * void some_function() {
 *     cusparseHandle_t cusparseHandle;
 *     cusparseWarnIfError(cusparseCreate(&cusparseHandle), "cusparseCreate(&cusparseHandle)", __FILE__, __func__,
 * __LINE__);
 * }
 * @endcode
 *
 * @note It is probably easier to use the macro OPM_GPUSPARSE_WARN_IF_ERROR
 * @note Prefer the gpuSPARSESafeCall/OPM_GPUSPARSE_SAFE_CALL counterpart unless you really don't want to throw an
 * exception.
 * @todo Refactor to use std::source_location once we shift to C++20
 */
inline cusparseStatus_t
cusparseWarnIfError(cusparseStatus_t error,
                    const std::string_view& expression,
                    const std::string_view& filename,
                    const std::string_view& functionName,
                    size_t lineNumber)
{
    if (error != CUSPARSE_STATUS_SUCCESS) {
        OpmLog::warning(getGpuSPARSEErrorMessage(error, expression, filename, functionName, lineNumber));
    }

    return error;
}
} // namespace Opm::gpuistl::detail



/**
 * @brief OPM_GPUSPARSE_SAFE_CALL checks the return type of the cusparse expression (function call) and throws an
 * exception if it does not equal CUSPARSE_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
 * #include <cusparse.h>
 *
 * void some_function() {
 *     cusparseHandle_t cusparseHandle;
 *     OPM_GPUSPARSE_SAFE_CALL(cusparseCreate(&cusparseHandle));
 * }
 * @endcode
 *
 * @note This should be used for any call to cuSparse unless you have a good reason not to.
 */
#define OPM_GPUSPARSE_SAFE_CALL(expression)                                                                             \
    ::Opm::gpuistl::detail::gpuSPARSESafeCall(expression, #expression, __FILE__, __func__, __LINE__)

/**
 * @brief OPM_GPUSPARSE_WARN_IF_ERROR checks the return type of the cusparse expression (function call) and issues a
 * warning if it does not equal CUSPARSE_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>
 * #include <cusparse.h>
 *
 * void some_function() {
 *     cusparseHandle_t cusparseHandle;
 *     OPM_GPUSPARSE_WARN_IF_ERROR(cusparseCreate(&cusparseHandle));
 * }
 * @endcode
 *
 * @note Prefer the gpuSPARSESafeCall/OPM_CUBLAS_SAFE_CALL counterpart unless you really don't want to throw an
 * exception.
 */
#define OPM_GPUSPARSE_WARN_IF_ERROR(expression)                                                                         \
    ::Opm::gpuistl::detail::cusparseWarnIfError(expression, #expression, __FILE__, __func__, __LINE__)
#endif // OPM_GPUSPARSE_SAFE_CALL_HPP
