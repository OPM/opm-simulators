/*
  Copyright 2023 SINTEF AS

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
#ifndef OPM_CUISTL_SAFE_CONVERSION_HPP
#define OPM_CUISTL_SAFE_CONVERSION_HPP



#include <cstddef>
#include <fmt/format.h>
#include <limits>
#include <opm/common/ErrorMacros.hpp>
#include <type_traits>


/**
 * Provides various utilities for doing signed to unsigned conversion, unsigned to signed, 32 bits to 64 bits and 64
 * bits to 32 bits.
 *
 * The main use case within cuistl is that the cusparse library requires signed int for all its size parameters,
 * while Dune::BlockVector (and relatives) use unsigned size_t.
 */

namespace Opm::cuistl::detail
{

/**
 * @brief to_int converts a (on most relevant platforms) 64 bits unsigned size_t to a signed 32 bits signed  int
 * @param s the unsigned integer
 * @throw std::invalid_argument exception if s is out of range for an int
 * @return converted s to int if s is within the range of int
 *
 * @todo This can be done for more generic types, but then it is probably wise to wait for C++20's cmp-functions
 */
inline int
to_int(std::size_t s)
{
    static_assert(
        std::is_signed_v<int>,
        "Weird architecture or my understanding of the standard is flawed. Better have a look at this function.");
    static_assert(
        !std::is_signed_v<std::size_t>,
        "Weird architecture or my understanding of the standard is flawed. Better have a look at this function.");

    static_assert(
        sizeof(int) <= sizeof(std::size_t),
        "Weird architecture or my understanding of the standard is flawed. Better have a look at this function.");


    if (s > std::size_t(std::numeric_limits<int>::max())) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("Trying to convert {} to int, but it is out of range. Maximum possible int: {}. ",
                              s,
                              std::numeric_limits<int>::max()));
    }

    // We know it will be in range here:
    return int(s);
}

/**
 * @brief to_size_t converts a (on most relevant platforms) a 32 bit signed int to a 64 bits unsigned int
 * @param i the signed integer
 * @return converted i to size_t if it is a non-negative integer.
 *
 * @throw std::invalid_argument if i is negative.
 * @todo This can be done for more generic types, but then it is probably wise to wait for C++20's cmp-functions
 */
inline std::size_t
to_size_t(int i)
{
    static_assert(
        std::is_signed_v<int>,
        "Weird architecture or my understanding of the standard is flawed. Better have a look at this function.");
    static_assert(
        !std::is_signed_v<std::size_t>,
        "Weird architecture or my understanding of the standard is flawed. Better have a look at this function.");

    static_assert(
        sizeof(int) <= sizeof(std::size_t),
        "Weird architecture or my understanding of the standard is flawed. Better have a look at this function.");


    if (i < int(0)) {
        OPM_THROW(std::invalid_argument, fmt::format("Trying to convert the negative number {} to size_t.", i));
    }

    return std::size_t(i);
}
} // namespace Opm::cuistl::detail

#endif
