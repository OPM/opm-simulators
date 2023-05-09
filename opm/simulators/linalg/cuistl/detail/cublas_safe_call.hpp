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
#ifndef OPM_CUBLAS_SAFE_CALL_HPP
#define OPM_CUBLAS_SAFE_CALL_HPP
#include <cublas_v2.h>
#include <exception>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <string_view>




namespace Opm::cuistl::detail
{

#define CHECK_CUBLAS_ERROR_TYPE(code, x)                                                                               \
    if (code == x) {                                                                                                   \
        return #x;                                                                                                     \
    }

namespace
{
    /**
     * @brief getCublasErrorCodeToString Converts an error code returned from a cublas function a human readable string.
     * @param code an error code from a cublas routine
     * @return a human readable string.
     */
    inline std::string getCublasErrorCodeToString(int code)
    {
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_SUCCESS);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_NOT_INITIALIZED);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_ALLOC_FAILED);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_INVALID_VALUE);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_ARCH_MISMATCH);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_MAPPING_ERROR);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_EXECUTION_FAILED);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_INTERNAL_ERROR);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_NOT_SUPPORTED);
        CHECK_CUBLAS_ERROR_TYPE(code, CUBLAS_STATUS_LICENSE_ERROR);

        return fmt::format("UNKNOWN CUBLAS ERROR {}.", code);
    }

#undef CHECK_CUBLAS_ERROR_TYPE

} // namespace

/**
 * @brief getCublasErrorMessage generates the error message to display for a given error.
 *
 * @param error the error code from cublas
 * @param expression the expresison (say "cublasCreate(&handle)")
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
getCublasErrorMessage(cublasStatus_t error,
                      const std::string_view& expression,
                      const std::string_view& filename,
                      const std::string_view& functionName,
                      size_t lineNumber)
{
    return fmt::format("cuBLAS expression did not execute correctly. Expression was: \n\n"
                       "    {}\n\n"
                       "in function {}, in {}, at line {}.\n"
                       "CuBLAS error code was: {}\n",
                       expression,
                       functionName,
                       filename,
                       lineNumber,
                       getCublasErrorCodeToString(error));
}

/**
 * @brief cublasSafeCall checks the return type of the CUBLAS expression (function call) and throws an exception if it
 * does not equal CUBLAS_STATUS_SUCCESS.
 *
 * @param error the error code from cublas
 * @param expression the expresison (say "cublasCreate(&handle)")
 * @param filename the code file the error occured in (typically __FILE__)
 * @param functionName name of the function the error occured in (typically __func__)
 * @param lineNumber the line number the error occured in (typically __LINE__)
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
 * #include <cublas_v2.h>
 *
 * void some_function() {
 *     cublasHandle_t cublasHandle;
 *     cudaSafeCall(cublasCreate(&cublasHandle), "cublasCreate(&cublasHandle)", __FILE__, __func__, __LINE__);
 * }
 * @endcode
 *
 * @note It is probably easier to use the macro OPM_CUBLAS_SAFE_CALL
 *
 * @todo Refactor to use std::source_location once we shift to C++20
 */
inline void
cublasSafeCall(cublasStatus_t error,
               const std::string_view& expression,
               const std::string_view& filename,
               const std::string_view& functionName,
               size_t lineNumber)
{
    if (error != CUBLAS_STATUS_SUCCESS) {
        OPM_THROW(std::runtime_error, getCublasErrorMessage(error, expression, filename, functionName, lineNumber));
    }
}

/**
 * @brief cublasWarnIfError checks the return type of the CUBLAS expression (function call) and issues a warning if it
 * does not equal CUBLAS_STATUS_SUCCESS.
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
 * #include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
 * #include <cublas_v2.h>
 *
 * void some_function() {
 *     cublasHandle_t cublasHandle;
 *     cublasWarnIfError(cublasCreate(&cublasHandle), "cublasCreate(&cublasHandle)", __FILE__, __func__, __LINE__);
 * }
 * @endcode
 *
 * @note It is probably easier to use the macro OPM_CUBLAS_WARN_IF_ERROR
 * @note Prefer the cublasSafeCall/OPM_CUBLAS_SAFE_CALL counterpart unless you really don't want to throw an exception.
 *
 * @todo Refactor to use std::source_location once we shift to C++20
 */
inline cublasStatus_t
cublasWarnIfError(cublasStatus_t error,
                  const std::string_view& expression,
                  const std::string_view& filename,
                  const std::string_view& functionName,
                  size_t lineNumber)
{
    if (error != CUBLAS_STATUS_SUCCESS) {
        OpmLog::warning(getCublasErrorMessage(error, expression, filename, functionName, lineNumber));
    }

    return error;
}
} // namespace Opm::cuistl::detail

/**
 * @brief OPM_CUBLAS_SAFE_CALL checks the return type of the cublas expression (function call) and throws an exception
 * if it does not equal CUBLAS_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <cublas_v2.h>
 * #include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
 *
 * void some_function() {
 *     cublasHandle_t cublasHandle;
 *     OPM_CUBLAS_SAFE_CALL(cublasCreate(&cublasHandle));
 * }
 * @endcode
 *
 * @note This should be used for any call to cuBlas unless you have a good reason not to.
 */
#define OPM_CUBLAS_SAFE_CALL(expression)                                                                               \
    ::Opm::cuistl::detail::cublasSafeCall(expression, #expression, __FILE__, __func__, __LINE__)

/**
 * @brief OPM_CUBLAS_WARN_IF_ERROR checks the return type of the cublas expression (function call) and issues a warning
 * if it does not equal CUBLAS_STATUS_SUCCESS.
 *
 * Example usage:
 * @code{.cpp}
 * #include <cublas_v2.h>
 * #include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
 *
 * void some_function() {
 *     cublasHandle_t cublasHandle;
 *     OPM_CUBLAS_WARN_IF_ERROR(cublasCreate(&cublasHandle));
 * }
 * @endcode
 *
 * @note Prefer the cublasSafeCall/OPM_CUBLAS_SAFE_CALL counterpart unless you really don't want to throw an exception.
 */
#define OPM_CUBLAS_WARN_IF_ERROR(expression)                                                                           \
    ::Opm::cuistl::detail::cublasWarnIfError(expression, #expression, __FILE__, __func__, __LINE__)

#endif // OPM_CUBLAS_SAFE_CALL_HPP
