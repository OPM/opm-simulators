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
 *
 * \copydoc Opm::ConditionalStorage
 */
#ifndef OPM_CONDITIONAL_STORAGE_HH
#define OPM_CONDITIONAL_STORAGE_HH

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>

#include <utility>

namespace Opm {
/*!
 * \ingroup Common
 *
 * \brief A simple class which only stores a given member attribute if a boolean
 *        condition is true
 *
 * If the condition is false, nothing is stored and an exception is thrown when trying to
 * access the object.
 */
template <bool cond, class T>
class ConditionalStorage
{
public:
    typedef T type;
    static constexpr bool condition = cond;

    ConditionalStorage()
    {}

    ConditionalStorage(const T& v)
        : data_(v)
    {}

    ConditionalStorage(T&& v)
        : data_(std::move(v))
    {}

    template <class ...Args>
    ConditionalStorage(Args... args)
        : data_(args...)
    {}

    ConditionalStorage(const ConditionalStorage& t)
        : data_(t.data_)
    {};

    ConditionalStorage(ConditionalStorage&& t)
        : data_(std::move(t.data_))
    {};

    ConditionalStorage& operator=(const ConditionalStorage& v)
    {
        data_ = v.data_;
        return *this;
    }

    ConditionalStorage& operator=(ConditionalStorage&& v)
    {
        data_ = std::move(v.data_);
        return *this;
    }

    const T& operator*() const
    { return data_; }
    T& operator*()
    { return data_; }

    const T* operator->() const
    { return &data_; }
    T* operator->()
    { return &data_; }

private:
    T data_;
};

template <class T>
class ConditionalStorage<false, T>
{
public:
    typedef T type;
    static constexpr bool condition = false;

    ConditionalStorage()
    {
        // ensure that T has a default constructor without actually calling it
        if (false) {
            T OPM_UNUSED dummy; // <- if the compiler bails out here, T does not have a default constructor
        }
    }

    ConditionalStorage(const T& v)
    {
        // ensure that T has a default constructor without actually calling it
        if (false) {
            T OPM_UNUSED dummy(v); // <- if the compiler bails out here, T does not have a copy constructor
        }
    }

    ConditionalStorage(const ConditionalStorage &)
    {
        // copying an empty conditional storage object does not do anything.
    };

    template <class ...Args>
    ConditionalStorage(Args... args)
    {
        // ensure that the arguments are valid without actually calling the constructor
        // of T
        if (false) {
            T OPM_UNUSED dummy(args...); // <- if the compiler bails out here, T does not have the requested constructor
        }
    }

    ConditionalStorage& operator=(const ConditionalStorage&)
    {
        // ensure that the stored object can actually be assined to but not actually do
        // anything
        if (false) {
            T *dummy;
            (*dummy) = (*dummy);   // <- if the compiler bails out here, T does not have an assignment operator
        }

        return *this;
    }

    const T& operator*() const
    { throw std::logic_error("data member deactivated"); }
    T& operator*()
    { throw std::logic_error("data member deactivated"); }

    const T* operator->() const
    { throw std::logic_error("data member deactivated"); }
    T* operator->()
    { throw std::logic_error("data member deactivated"); }
};

} // namespace Opm

#endif
