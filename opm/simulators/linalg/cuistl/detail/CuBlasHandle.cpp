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
namespace Opm::gpuistl::detail
{


GpuBLASHandle::GpuBLASHandle()
{
    OPM_GPU_BLAS_SAFE_CALL(cublasCreate(&m_handle));
}

GpuBLASHandle::~GpuBLASHandle()
{
    OPM_CUBLAS_WARN_IF_ERROR(cublasDestroy(m_handle));
}

cublasHandle_t
GpuBLASHandle::get()
{
    return m_handle;
}

GpuBLASHandle&
GpuBLASHandle::getInstance()
{
    static GpuBLASHandle instance;
    return instance;
}

} // namespace Opm::gpuistl::detail
