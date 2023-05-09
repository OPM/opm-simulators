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
#ifndef OPM_CUSPARSEHANDLE_HPP
#define OPM_CUSPARSEHANDLE_HPP
#include <cusparse.h>
#include <memory>

namespace Opm::cuistl::detail
{

/**
 * @brief The CuSparseHandle class provides a singleton for the simulator universal cuSparseHandle.
 *
 * Example use:
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/detail/CuSparseHandle.hpp>
 * void someFunction() {
 *     auto& cuSparseHandle = ::Opm::cuistl::detail::CuSparseHandle::getInstance();
 *     int cuSparseVersion = -1;
 *     OPM_CUSPARSE_SAFE_CALL(cusparseGetVersion(cuSparseHandle.get(), &cuSparseVersion));
 * }
 * @endcode
 */
class CuSparseHandle
{
public:
    // This should not be copyable.
    CuSparseHandle(const CuSparseHandle&) = delete;
    CuSparseHandle& operator=(const CuSparseHandle&) = delete;

    /**
     * Calls cuSparseDestroy on the handle
     */
    ~CuSparseHandle();

    /**
     * @brief get returns the underlying cuSparse handle (to be used in calls to cusparse)
     */
    cusparseHandle_t get();

    /**
     * @brief getInstance creates (if necessary) and returns the single unique instance of CuSparseHandle (singleton)
     */
    static CuSparseHandle& getInstance();

private:
    CuSparseHandle();
    cusparseHandle_t m_handle;
};
} // namespace Opm::cuistl::detail
#endif // OPM_CUSPARSEHANDLE_HPP
