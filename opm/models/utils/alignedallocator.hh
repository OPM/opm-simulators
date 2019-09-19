// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief This is a stand-alone version of boost::alignment::aligned_allocator from Boost
 *        1.58
 *
 * The file has been modified to assume a C++-2011 compatible compiler on a POSIX
 * operating system to remove the boost dependencies which the original version
 * contained. The original copyright notice for this file is:
 *
<pre>
 (c) 2014 Glen Joseph Fernandes
 glenjofe at gmail dot com

 Distributed under the Boost Software
 License, Version 1.0.
 http://boost.org/LICENSE_1_0.txt
</pre>
*/
#ifndef EWOMS_ALIGNED_ALLOCATOR_HH
#define EWOMS_ALIGNED_ALLOCATOR_HH

#include <utility>
#include <memory>
#include <type_traits>
#include <cassert>

namespace Opm {

namespace detail {
constexpr inline bool is_alignment(std::size_t value) noexcept
{
    return (value > 0) && ((value & (value - 1)) == 0);
}

template<std::size_t N>
struct is_alignment_constant
    : std::integral_constant<bool, (N > 0) && ((N & (N - 1)) == 0)>
{};

template<std::size_t A, std::size_t B>
struct min_size
    : std::integral_constant<std::size_t, (A < B) ? A : B>
{ };

template<class T>
struct offset_object
{
    char offset;
    T object;
};

template<class T>
struct alignment_of
    : min_size<sizeof(T), sizeof(offset_object<T>) - sizeof(T)>::type
{};

template<std::size_t A, std::size_t B>
struct max_align
    : std::integral_constant<std::size_t,(A > B) ? A : B>
{};

template<class T>
struct max_count_of
    : std::integral_constant<std::size_t, ~static_cast<std::size_t>(0) / sizeof(T)>
{};

using std::addressof;
}

inline void* aligned_alloc(std::size_t alignment,
                           std::size_t size) noexcept
{
    assert(detail::is_alignment(alignment));
    if (alignment < sizeof(void*)) {
        alignment = sizeof(void*);
    }
    void* p;
    if (::posix_memalign(&p, alignment, size) != 0) {
        p = 0;
    }
    return p;
}

inline void aligned_free(void* ptr)
    noexcept
{
    ::free(ptr);
}


template<class T, std::size_t Alignment>
class aligned_allocator {
    static_assert(detail::is_alignment_constant<Alignment>::value, "Alignment must be powers of two!");

public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef void* void_pointer;
    typedef const void* const_void_pointer;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef T& reference;
    typedef const T& const_reference;

private:
    typedef detail::max_align<Alignment, detail::alignment_of<value_type>::value> MaxAlign;

public:
    template<class U>
    struct rebind {
        typedef aligned_allocator<U, Alignment> other;
    };

    aligned_allocator()
        noexcept = default;

    template<class U>
    aligned_allocator(const aligned_allocator<U,
                      Alignment>&) noexcept {
    }

    pointer address(reference value) const
        noexcept {
        return detail::addressof(value);
    }

    const_pointer address(const_reference value) const
        noexcept {
        return detail::addressof(value);
    }

    pointer allocate(size_type size,
                     const_void_pointer = 0) {
        void* p = aligned_alloc(MaxAlign::value,
                                sizeof(T) * size);
        if (!p && size > 0) {
            throw std::bad_alloc();
        }
        return static_cast<T*>(p);
    }

    void deallocate(pointer ptr, size_type) {
        aligned_free(ptr);
    }

    constexpr size_type max_size() const
        noexcept {
        return detail::max_count_of<T>::value;
    }

    template<class U, class... Args>
    void construct(U* ptr, Args&&... args) {
        void* p = ptr;
        ::new(p) U(std::forward<Args>(args)...);
    }

    template<class U>
    void construct(U* ptr) {
        void* p = ptr;
        ::new(p) U();
    }

    template<class U>
    void destroy(U* ptr) {
        (void)ptr;
        ptr->~U();
    }
};

template<std::size_t Alignment>
class aligned_allocator<void, Alignment> {
    static_assert(detail::is_alignment_constant<Alignment>::value,
                  "The specified alignment is not a power of two!");

public:
    typedef void value_type;
    typedef void* pointer;
    typedef const void* const_pointer;

    template<class U>
    struct rebind {
        typedef aligned_allocator<U, Alignment> other;
    };
};

template<class T1, class T2, std::size_t Alignment>
inline bool operator==(const aligned_allocator<T1,
                       Alignment>&, const aligned_allocator<T2,
                       Alignment>&) noexcept
{
    return true;
}

template<class T1, class T2, std::size_t Alignment>
inline bool operator!=(const aligned_allocator<T1,
                       Alignment>&, const aligned_allocator<T2,
                       Alignment>&) noexcept
{
    return false;
}
}

#endif
