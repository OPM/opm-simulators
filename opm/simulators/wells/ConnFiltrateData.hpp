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

#ifndef OPM_CONNFILTRATEDATA_HPP
#define OPM_CONNFILTRATEDATA_HPP

#include <vector>

namespace Opm {
    struct ConnFiltrateData {

        ConnFiltrateData() = default;

        void resize(std::size_t num_perf);

        template<class Serializer>
        void serializeOp(Serializer& serializer) {
            serializer(rates);
            serializer(total);
            serializer(skin_factor);
            serializer(thickness);
            serializer(perm);
            serializer(poro);
            serializer(radius);
            serializer(area_of_flow);
        }

        static ConnFiltrateData serializationTestObject();

        bool operator==(const ConnFiltrateData& rhs) const;

        std::vector<double> rates;
        std::vector<double> total;
        std::vector<double> skin_factor;
        std::vector<double> thickness;
        std::vector<double> perm;
        std::vector<double> poro;
        std::vector<double> radius;
        std::vector<double> area_of_flow;
    };
}

#endif //OPM_CONNFILTRATEDATA_HPP
