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
#include <cublas_v2.h>
#include <opm/simulators/linalg/cuistl/detail/CuBlasHandle.hpp>
#include <opm/simulators/linalg/cuistl/detail/cublas_safe_call.hpp>
namespace Opm::cuistl::detail
{


CuBlasHandle::CuBlasHandle()
{
    OPM_CUBLAS_SAFE_CALL(cublasCreate(&m_handle));
}

CuBlasHandle::~CuBlasHandle()
{
    OPM_CUBLAS_WARN_IF_ERROR(cublasDestroy(m_handle));
}

cublasHandle_t
CuBlasHandle::get()
{
    return m_handle;
}

CuBlasHandle&
CuBlasHandle::getInstance()
{
    static CuBlasHandle instance;
    return instance;
}

} // namespace Opm::cuistl::detail
