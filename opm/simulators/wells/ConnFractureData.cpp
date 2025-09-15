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

#include <opm/simulators/wells/ConnFractureData.hpp>

namespace Opm {

template<class Scalar>
void ConnFractureData<Scalar>::resize(std::size_t num_perf)
{
    this->area.resize(num_perf);
    this->flux.resize(num_perf);
    this->height.resize(num_perf);
    this->length.resize(num_perf);
}

template<class Scalar>
ConnFractureData<Scalar>
ConnFractureData<Scalar>::serializationTestObject()
{
    ConnFractureData result;
    result.area = {8.};
    result.flux = {100.};
    result.height = {0.5};
    result.length = {0.05};
    return result;
}

template<class Scalar>
bool ConnFractureData<Scalar>::operator==(const ConnFractureData& rhs) const
{
    return this->area == rhs.area &&
           this->flux == rhs.flux &&
           this->height == rhs.height &&
           this->length == rhs.length; 
}

template struct ConnFractureData<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct ConnFractureData<float>;
#endif

}
