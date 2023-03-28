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
#ifndef CUSPARSERESOURCE_HPP
#define CUSPARSERESOURCE_HPP
#include <cusparse.h>
#include <functional>
#include <memory>
#include <type_traits>

namespace Opm::cuistl::detail
{

/**
 * @brief The CuSparseResource class wraps a CuSparse resource in a proper RAII pattern
 *
 * Current we support the following types for T:
 *   - bsrilu02Info_t
 *   - bsrsv2Info_t
 *   - cusparseMatDescr_t
 *
 * More types are in principle supported by supplying a manual Creator and Destructor.
 *
 * In addition to acting as an easier-to-use smart_ptr specialized for these types, it
 * also adds error checking in the construction and deletion of these resources,
 * which a plain std::smart_ptr would not support out of the box. It also solves the
 * caveat of the pointer types of cuSparse resources not being exposed properly.
 *
 * Example usage:
 * @code{.cpp}
 * #include <opm/simulator/linalg/cuistl/detail/CuSparseResource.hpp>
 *
 * void someFunction() {
 *     auto resource = CuSparseResource<cuSparseMatDescr_t>();
 * }
 * @endcode
 */
template <class T>
class CuSparseResource
{
public:
    using CreatorType = typename std::function<cusparseStatus_t(T*)>;
    using DeleterType = typename std::function<cusparseStatus_t(T)>;

    /**
     * @brief CuSparseResource creates a new instance by calling creator, and will delete using deleter
     * @param creator a functor used to create an instance
     * @param deleter a functor used to delete the instance
     *
     * @note Using this constructor it is possible to add support for new types not already accounted for.
     */
    CuSparseResource(CreatorType creator, DeleterType deleter);

    /**
     * @brief CuSparseResource will automatically select the proper creator and deleter based on the type (and throw an exception if not available)
     */
    CuSparseResource();

    /**
     * Calls the deleter functor
     */
    ~CuSparseResource();

    // This should not be copyable.
    CuSparseResource(const CuSparseResource&) = delete;
    CuSparseResource& operator=(const CuSparseResource&) = delete;

    /**
     * @brief get returns the raw pointer to the resource.
     */
    T get()
    {
        return m_resource;
    }

private:
    T m_resource;

    DeleterType m_deleter;
};

} // namespace Opm::cuistl::impl
#include <opm/simulators/linalg/cuistl/detail/CuSparseResource_impl.hpp>
#endif // CUSPARSERESOURCE_HPP
