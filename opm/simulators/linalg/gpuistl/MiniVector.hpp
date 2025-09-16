/*
  Copyright 2025 Equinor ASA
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

#ifndef OPM_GPUISTL_MINIVECTOR_HPP
#define OPM_GPUISTL_MINIVECTOR_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

/**
 * @brief A small, fixed‑dimension MiniVector class backed by std::array that can be
 *        used in both host and CUDA device code.
 *
 * The implementation purposefully remains lightweight, containing only the
 * utilities required for element access, iteration, and initialization. It
 * avoids dynamic memory and leverages the compile‑time `Dimension` parameter
 * to enable full constexpr evaluation whenever possible.
 */

namespace Opm::gpuistl
{

/**
 * @tparam T         Element type (must be trivially copyable for use on the GPU).
 * @tparam Dimension Compile‑time dimension (must be strictly positive).
 */
template <class T, int Dimension>
class MiniVector
{
    static_assert(Dimension > 0, "Dimension must be positive");

public:
    //! Element type.
    using value_type = T;
    //! Index/size type.
    using size_type = std::size_t;
    //! Mutable element reference.
    using reference = value_type&;
    //! Immutable element reference.
    using const_reference = const value_type&;
    //! Mutable iterator.
    using iterator = typename std::array<T, Dimension>::iterator;
    //! Immutable iterator.
    using const_iterator = typename std::array<T, Dimension>::const_iterator;

    /**
     * @brief Default‑constructs the MiniVector; elements are value‑initialized.
     */
    OPM_HOST_DEVICE constexpr MiniVector() noexcept(std::is_nothrow_default_constructible<value_type>::value) = default;

    /**
     * @brief Uniform‑value constructor.
     * @param value Value used to initialize **all** components.
     */
    OPM_HOST_DEVICE constexpr explicit MiniVector(const value_type& value)
    {
        data_.fill(value);
    }

    /**
     * @brief Initializer‑list constructor.
     *
     * Enables natural brace initialization:
     *
     * @code
     * Opm::MiniVector<float, 3> v{1.0f, 2.0f, 3.0f};
     * @endcode
     *
     * @param init List containing exactly `Dimension` elements.
     * @throws std::runtime_error if `init.size() != Dimension` **(host only)**.
     */
    OPM_HOST_DEVICE MiniVector(std::initializer_list<value_type> init)
    {
        if (init.size() != Dimension) {
            OPM_THROW(std::runtime_error, "Opm::MiniVector – initializer‑list size mismatch");
        }
        std::copy_n(init.begin(), Dimension, data_.begin());
    }

    /**
     * @return Mutable reference to element `idx` (no bounds check).
     */
    OPM_HOST_DEVICE constexpr reference operator[](size_type idx) noexcept
    {
        return data_[idx];
    }

    /** @return Immutable reference to element `idx` (no bounds check). */
    OPM_HOST_DEVICE constexpr const_reference operator[](size_type idx) const noexcept
    {
        return data_[idx];
    }

    /**
     * @brief Safe element access with bounds checking (throws on host).
     * @throws std::out_of_range if `idx >= Dimension` **(host only)**.
     */
    OPM_HOST_DEVICE reference at(size_type idx)
    {
        if (idx >= Dimension) {
            OPM_THROW(std::out_of_range, "Opm::MiniVector::at – index out of range");
        }
        return data_[idx];
    }
    /** @copydoc at(size_type) */
    OPM_HOST_DEVICE const_reference at(size_type idx) const
    {
        if (idx >= Dimension) {
            OPM_THROW(std::out_of_range, "Opm::MiniVector::at – index out of range");
        }
        return data_[idx];
    }

    /** @return Iterator to first element (mutable). */
    OPM_HOST_DEVICE constexpr iterator begin() noexcept
    {
        return data_.begin();
    }
    /** @return Const iterator to first element. */
    OPM_HOST_DEVICE constexpr const_iterator begin() const noexcept
    {
        return data_.begin();
    }
    /** @return Const iterator to first element. */
    OPM_HOST_DEVICE constexpr const_iterator cbegin() const noexcept
    {
        return data_.cbegin();
    }

    /** @return One‑past‑the‑last iterator (mutable). */
    OPM_HOST_DEVICE constexpr iterator end() noexcept
    {
        return data_.end();
    }
    /** @return Const one‑past‑the‑last iterator. */
    OPM_HOST_DEVICE constexpr const_iterator end() const noexcept
    {
        return data_.end();
    }
    /** @return Const one‑past‑the‑last iterator. */
    OPM_HOST_DEVICE constexpr const_iterator cend() const noexcept
    {
        return data_.cend();
    }

    /** @return The fixed dimension of the MiniVector (compile‑time constant). */
    OPM_HOST_DEVICE static constexpr size_type size() noexcept
    {
        return Dimension;
    }


    /**
     * @brief Fill every component with the supplied value.
     * @param value The value to assign to each element.
     */
    OPM_HOST_DEVICE constexpr void fill(const value_type& value)
    {
        data_.fill(value);
    }

    /** @return `true` if all components compare equal. */
    OPM_HOST_DEVICE constexpr bool operator==(const MiniVector& other) const noexcept
    {
        return std::equal(data_.begin(), data_.end(), other.data_.begin());
    }

    /** @return `true` if any component differs. */
    OPM_HOST_DEVICE constexpr bool operator!=(const MiniVector& other) const noexcept
    {
        return !(*this == other);
    }

private:
    //! Element storage.
    std::array<value_type, Dimension> data_ {};
};

} // namespace Opm::gpuistl

#endif // OPM_GPUISTL_MINIVECTOR_HPP
