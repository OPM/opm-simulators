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

#include <config.h>
#include <opm/simulators/flow/ExtboContainer.hpp>

namespace Opm {

template<class Scalar>
void ExtboContainer<Scalar>::
allocate(const unsigned bufferSize)
{
    X_volume_.resize(bufferSize, 0.0);
    Y_volume_.resize(bufferSize, 0.0);
    Z_fraction_.resize(bufferSize, 0.0);
    mFracOil_.resize(bufferSize, 0.0);
    mFracGas_.resize(bufferSize, 0.0);
    mFracCo2_.resize(bufferSize, 0.0);

    allocated_ = true;
}

template class ExtboContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class ExtboContainer<float>;
#endif

} // namespace Opm
