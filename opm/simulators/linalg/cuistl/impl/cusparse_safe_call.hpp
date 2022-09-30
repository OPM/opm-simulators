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
#ifndef CUSPARSE_SAFE_CALL_HPP
#define CUSPARSE_SAFE_CALL_HPP
#include <cusparse.h>
#include <exception>
#include <opm/common/ErrorMacros.hpp>
#define CHECK_CUSPARSE_ERROR_TYPE(code, x)                                                                             \
    if (code == x) {                                                                                                   \
        return #x;                                                                                                     \
    }
namespace
{
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
    return "UNKNOWN CUSPARSE ERROR " + std::to_string(code);
}
} // namespace
#undef CHECK_CUSPARSE_ERROR_TYPE
#define OPM_CUSPARSE_SAFE_CALL(expression)                                                                             \
    do {                                                                                                               \
        cusparseStatus_t error = expression;                                                                           \
        if (error != CUSPARSE_STATUS_SUCCESS) {                                                                        \
            OPM_THROW(std::runtime_error,                                                                              \
                      "cuSparse expression did not execute correctly. Expression was: \n\n"                            \
                          << "    " << #expression << "\n\n"                                                           \
                          << "in function " << __func__ << ", in " << __FILE__ << " at line " << __LINE__ << "\n"      \
                          << "CuSparse error code was: " << getCusparseErrorMessage(error) << "\n");                   \
        }                                                                                                              \
    } while (false)
#endif // CUSPARSE_SAFE_CALL_HPP
