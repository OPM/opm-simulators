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
#include <config.h>

#define BOOST_TEST_MODULE TestCublasSafeCall

#include <boost/test/unit_test.hpp>
#include <cublas_v2.h>
#include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>

BOOST_AUTO_TEST_CASE(TestCreateHandle)
{
    cublasHandle_t cublasHandle;
    BOOST_CHECK_NO_THROW(OPM_CUBLAS_SAFE_CALL(cublasCreate(&cublasHandle)););
}

BOOST_AUTO_TEST_CASE(TestThrows)
{
    std::vector<cublasStatus_t> errorCodes {{CUBLAS_STATUS_NOT_INITIALIZED,
                                             CUBLAS_STATUS_ALLOC_FAILED,
                                             CUBLAS_STATUS_INVALID_VALUE,
                                             CUBLAS_STATUS_ARCH_MISMATCH,
                                             CUBLAS_STATUS_MAPPING_ERROR,
                                             CUBLAS_STATUS_EXECUTION_FAILED,
                                             CUBLAS_STATUS_INTERNAL_ERROR,
                                             CUBLAS_STATUS_NOT_SUPPORTED,
                                             CUBLAS_STATUS_LICENSE_ERROR}};
    for (auto code : errorCodes) {
        BOOST_CHECK_THROW(OPM_CUBLAS_SAFE_CALL(code), std::exception);
    }
}
