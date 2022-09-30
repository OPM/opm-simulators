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
#ifndef CUSPARSERESOURCE_HPP
#define CUSPARSERESOURCE_HPP
#include <cusparse.h>
#include <functional>
#include <memory>
#include <type_traits>

namespace Opm::cuistl::impl
{

//! @brief wraps a CuSparseResource in proper RAII.
template <class T>
class CuSparseResource
{
public:
    using CreatorType = typename std::function<cusparseStatus_t(T*)>;
    using DeleterType = typename std::function<cusparseStatus_t(T)>;
    CuSparseResource(CreatorType creator, DeleterType deleter);
    CuSparseResource();
    ~CuSparseResource();

    // This should not be copyable.
    CuSparseResource(const CuSparseResource&) = delete;
    CuSparseResource& operator=(const CuSparseResource&) = delete;

    T get()
    {
        return m_resource;
    }

private:
    T m_resource;

    DeleterType m_deleter;
};

} // namespace Opm::cuistl::impl
#include <opm/simulators/linalg/cuistl/impl/CuSparseResource_impl.hpp>
#endif // CUSPARSERESOURCE_HPP
