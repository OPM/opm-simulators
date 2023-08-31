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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/wells/ConnFiltrateData.hpp>

namespace Opm {

    void ConnFiltrateData::resize(std::size_t num_perf) {
        this->rates.resize(num_perf);
        this->total.resize(num_perf);
        this->skin_factor.resize(num_perf);
        this->thickness.resize(num_perf);
        this->perm.resize(num_perf);
        this->poro.resize(num_perf);
        this->radius.resize(num_perf);
        this->area_of_flow.resize(num_perf);
    }

    ConnFiltrateData ConnFiltrateData::serializationTestObject()
    {
        ConnFiltrateData result;
        result.rates = {8.};
        result.total = {100.};
        result.skin_factor = {0.5};
        result.thickness = {0.05};
        result.perm = {0.00001};
        result.poro = {0.3};
        result.radius = {0.05};
        result.area_of_flow = {0.7};
        return result;
    }

    bool ConnFiltrateData::operator==(const ConnFiltrateData& rhs) const
    {
        return this->rates == rhs.rates &&
               this->total == rhs.total &&
               this->skin_factor == rhs.skin_factor &&
               this->thickness == rhs.thickness &&
               this->perm == rhs.perm &&
               this->poro == rhs.poro &&
               this->radius == rhs.radius &&
               this->area_of_flow == rhs.area_of_flow;
    }
}