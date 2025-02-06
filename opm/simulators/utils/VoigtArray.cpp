/*
  Copyright 2025 Equinor ASA

  This file is part of the Open Porous Media Project (OPM).

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
#include <config.h>
#include <opm/simulators/utils/VoigtArray.hpp>

#include <dune/common/fvector.hh>

#include <algorithm>

namespace Opm {

template<class T>
template<class Array>
VoigtContainer<T>::
VoigtContainer(const Array& array)
{
    std::copy(array.begin(), array.end(), data_.begin());
}

template<class Scalar>
VoigtArray<Scalar>::
VoigtArray(const std::size_t size)
{
    this->resize(size);
}

template<class Scalar>
void VoigtArray<Scalar>::
resize(const std::size_t size)
{
    std::for_each(this->data_.begin(), this->data_.end(),
                  [size](auto& d) { d.resize(size); });
}

template<class Scalar>
void VoigtArray<Scalar>::
assign(const std::size_t i, const VoigtContainer<Scalar>& array)
{
    for (const auto idx : this->unique_indices) {
        (*this)[idx][i] = array[idx];
    }
}

template<class Scalar>
Scalar VoigtArray<Scalar>::
operator()(const VoigtIndex idx, const std::size_t i) const
{
    return (*this)[idx].at(i);
}

template<class Scalar>
Scalar& VoigtArray<Scalar>::
operator()(const VoigtIndex idx, const std::size_t i)
{
    return (*this)[idx].at(i);
}

#define INSTANTIATE_TYPE(T)                                                    \
    template class VoigtArray<T>;                                              \
    template VoigtContainer<T>::VoigtContainer(const std::array<T,6>&);        \
    template VoigtContainer<T>::VoigtContainer(const Dune::FieldVector<T,6>&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
