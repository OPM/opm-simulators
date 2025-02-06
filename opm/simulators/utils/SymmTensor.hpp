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

#ifndef OPM_UTIL_SYMM_TENSOR_HPP
#define OPM_UTIL_SYMM_TENSOR_HPP

#include <dune/common/fvector.hh>

#include <opm/simulators/utils/VoigtArray.hpp>

#include <utility>

namespace Opm {

template<class T>
class SymmTensor : public VoigtContainer<T>
{
public:
    using field_type = T;

    SymmTensor() = default;
    SymmTensor(std::initializer_list<T> value)
        : VoigtContainer<T>(std::move(value))
    {}

    void operator+=(const T data);
    void operator+=(const SymmTensor<T>& data);

    void operator*=(const T data);

    SymmTensor<T>& operator=(const T value);

    void reset();

    T trace() const;
    T traction(const Dune::FieldVector<T,3>& normal) const;
};

template<class T1, class T2>
SymmTensor<T1> operator*(const T2 value, SymmTensor<T1> t1);
template<class T1, class T2>
SymmTensor<T1> operator*(SymmTensor<T1> t1, const T2 value);
template<class T>
SymmTensor<T> operator+(SymmTensor<T> t1, const SymmTensor<T>& t2);

} // namespace Opm

#endif // OPM_UTIL_SYMM_TENSOR_HPP
