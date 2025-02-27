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
#ifndef OPM_WELL_TRACER_RATE_HPP
#define OPM_WELL_TRACER_RATE_HPP

#include <string>

namespace Opm {

template<class Scalar>
struct WellTracerRate
{
    std::string name{};
    Scalar rate{};

    WellTracerRate() = default;

    WellTracerRate(const std::string& n, const Scalar r)
        : name(n), rate(r)
    {}

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(name);
        serializer(rate);
    }

    bool operator==(const WellTracerRate& that) const
    {
        return
               this->name == that.name
            && this->rate == that.rate;
    }
};

} // namespace Opm

#endif // OPM_WELL_TRACER_RATE_HPP
