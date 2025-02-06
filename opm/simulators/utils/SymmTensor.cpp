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
#include <opm/simulators/utils/SymmTensor.hpp>

#include <algorithm>

namespace Opm {

template<class T>
void SymmTensor<T>::
operator+=(const T data)
{
    std::for_each(this->data_.begin(), this->data_.end(),
                  [&data](auto& v) { v += data; });
}

template<class T>
void SymmTensor<T>::
operator+=(const SymmTensor& data)
{
    auto it = data.data_.begin();
    std::for_each(this->data_.begin(), this->data_.end(),
                  [&it](auto& v) { v += *it++; });
}

template<class T>
void SymmTensor<T>::
operator*=(const T value)
{
    std::for_each(this->data_.begin(), this->data_.end(),
                  [&value](auto& v) { v *= value; });
}

template<class T>
void SymmTensor<T>::reset()
{
    this->data_.fill(T{0});
}

template<class T>
T SymmTensor<T>::
trace() const
{
    return (*this)[VoigtIndex::XX] +
           (*this)[VoigtIndex::YY] +
           (*this)[VoigtIndex::ZZ];
}

template<class T>
T SymmTensor<T>::traction(const Dune::FieldVector<T,3>& normal) const
{
    T traction = 0.0;
    constexpr auto& ind = Opm::SymmTensor<double>::diag_indices;
    for (std::size_t i = 0; i < 3; ++i){
        traction += (*this)[ind[i]] * normal[i] * normal[i];
    }
    traction += T{2} * (*this)[VoigtIndex::YZ] * normal[0] * normal[1]; // xy*nx*ny;
    traction += T{2} * (*this)[VoigtIndex::XZ] * normal[0] * normal[2]; // xz*nx*nz
    traction += T{2} * (*this)[VoigtIndex::XY] * normal[1] * normal[2]; // yz*ny*nz

    return traction;
}

template<class T>
SymmTensor<T>&
SymmTensor<T>::operator=(const T value)
{
    this->data_.fill(value);
    return *this;
}

template<class T1, class T2>
SymmTensor<T1> operator*(const T2 value, SymmTensor<T1> t1)
{
    t1 *= value;
    return t1;
}

template<class T1, class T2>
SymmTensor<T1> operator*(SymmTensor<T1> t1, const T2 value)
{
    t1 *= value;
    return t1;
}

template<class T>
SymmTensor<T> operator+(SymmTensor<T> t1, const SymmTensor<T>& t2)
{
    t1 += t2;
    return t1;
}

#define INSTANTIATE_OPS(T1, T2)                                  \
    template SymmTensor<T1> operator*(const T2, SymmTensor<T1>); \
    template SymmTensor<T1> operator*(SymmTensor<T1>, const T2);

#define INSTANTIATE_TYPE(T)                                                \
    template class SymmTensor<T>;                                          \
    INSTANTIATE_OPS(T, T)                                                  \
    INSTANTIATE_OPS(T, int)                                                \
    INSTANTIATE_OPS(T, unsigned)                                           \
    INSTANTIATE_OPS(T, std::size_t)                                        \
    template SymmTensor<T> operator+(SymmTensor<T>, const SymmTensor<T>&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
