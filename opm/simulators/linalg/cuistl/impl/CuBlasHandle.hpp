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
#ifndef CUBLASHANDLE_HPP
#define CUBLASHANDLE_HPP
#include <cublas_v2.h>
#include <memory>

namespace Opm::cuistl::impl
{
class CuBlasHandle
{
public:
    // This should not be copyable.
    CuBlasHandle(const CuBlasHandle&) = delete;
    CuBlasHandle& operator=(const CuBlasHandle&) = delete;

    ~CuBlasHandle();

    cublasHandle_t get();

    static CuBlasHandle& getInstance();

private:
    CuBlasHandle();
    cublasHandle_t m_handle;
};
} // namespace Opm::cuistl::impl
#endif // CuBlasHandle_HPP
