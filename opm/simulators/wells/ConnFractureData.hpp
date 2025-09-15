/*
  Copyright 2023 Equinor ASA.


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

#ifndef OPM_CONNFRACTUREDATA_HPP
#define OPM_CONNFRACTUREDATA_HPP

#include <vector>

namespace Opm {

template<class Scalar>
struct ConnFractureData {

    void resize(std::size_t num_perf);

    template<class Serializer>
    void serializeOp(Serializer& serializer) {
        serializer(area);
        serializer(flux);
        serializer(height);
        serializer(length);
    }

    static ConnFractureData serializationTestObject();

    bool operator==(const ConnFractureData& rhs) const;

    std::vector<Scalar> area;
    std::vector<Scalar> flux;
    std::vector<Scalar> height;
    std::vector<Scalar> length;
};

}

#endif // OPM_CONNFILTRATEDATA_HPP
