/*
  Copyright 2019 Equinor AS.

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
#ifndef SERIALIZATION_PACKERS_HPP
#define SERIALIZATION_PACKERS_HPP

#include <opm/common/utility/MemPacker.hpp>

#include <boost/date_time/gregorian/gregorian_types.hpp>

// Additional packers for serializers using the mempacker.

namespace Opm {
namespace Serialization {
namespace detail {

template<>
struct Packing<false,boost::gregorian::date>
{
    static std::size_t packSize(const boost::gregorian::date& data);

    static void pack(const boost::gregorian::date& data,
                     std::vector<char>& buffer, int& position);

    static void unpack(boost::gregorian::date& data,
                       std::vector<char>& buffer, int& position);
};

}

} // end namespace Serialization
} // end namespace Opm

#endif // SERIALIZATION_PACKERS_HPP
