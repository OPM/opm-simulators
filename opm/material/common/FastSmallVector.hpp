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
 * \copydoc Opm::FastSmallVector
 */
#ifndef OPM_FAST_SMALL_VECTOR_HPP
#define OPM_FAST_SMALL_VECTOR_HPP

#include <array>
#include <cstring>

namespace Opm {

template <typename ValueType, unsigned N>
class FastSmallVector
{
public:
    FastSmallVector()
    {
        size_ = 0;
        dataPtr_ = smallBuf_.data();
    }

    explicit FastSmallVector(const size_t numElem)
    {
        init_(numElem);
    }

    FastSmallVector(const size_t numElem, const ValueType value)
    {
        init_(numElem);

        std::fill(dataPtr_, dataPtr_ + size_, value);
    }

    FastSmallVector(const FastSmallVector& other)
    {
        size_ = 0;
        dataPtr_ = smallBuf_.data();

        (*this) = other;
    }

    FastSmallVector(FastSmallVector&& other)
    {
        size_ = 0;
        dataPtr_ = smallBuf_.data();

        (*this) = std::move(other);
    }

    ~FastSmallVector()
    {
        if (dataPtr_ != smallBuf_.data())
            delete [] dataPtr_;
    }


    FastSmallVector& operator=(FastSmallVector&& other)
    {
        if (dataPtr_ != smallBuf_.data() && dataPtr_ != other.dataPtr_)
            delete [] dataPtr_;

        size_ = other.size_;
        if (size_ <= N) {
            smallBuf_ = std::move(other.smallBuf_);
            dataPtr_ = smallBuf_.data();
        }
        else
            dataPtr_ = other.dataPtr_;

        other.dataPtr_ = nullptr;
        other.size_ = 0;

        return (*this);
    }

    FastSmallVector& operator=(const FastSmallVector& other)
    {
        size_ = other.size_;

        if (size_ <= N) {
            smallBuf_ = other.smallBuf_;
            dataPtr_ = smallBuf_.data();
        }
        else if (dataPtr_ != other.dataPtr_) {
            if (dataPtr_ != smallBuf_.data())
                delete[] dataPtr_;
            dataPtr_ = new ValueType[size_];

            std::copy(other.dataPtr_, other.dataPtr_ + size_, dataPtr_);
        }

        return (*this);
    }

    ValueType& operator[](size_t idx)
    { return dataPtr_[idx]; }

    const ValueType& operator[](size_t idx) const
    { return dataPtr_[idx]; }

    size_t size() const
    { return size_; }

private:
    void init_(size_t numElem)
    {
        size_ = numElem;

        if (size_ > N)
            dataPtr_ = new ValueType[size_];
        else
            dataPtr_ = smallBuf_.data();
    }

    std::array<ValueType, N> smallBuf_;
    std::size_t size_;
    ValueType* dataPtr_;
};

} // namespace Opm

#endif // OPM_FAST_SMALL_VECTOR_HPP
