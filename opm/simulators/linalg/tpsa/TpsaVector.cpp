// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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

#include <opm/simulators/linalg/tpsa/TpsaVector.hpp>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

#include <algorithm>

namespace Opm::Linear
{

//
// TpsaVector<ScalarT>::EntryProxy implementation
//
// The conversion type is spelled out rather than left as EqVector: in an
// out-of-line definition MSVC does not look the conversion-type-id up in the
// enclosing class template's scope, and rejects the short form with
// "C7753: ill-formed conversion-function-id". The qualified spelling means
// the same thing to every compiler.
template <class ScalarT>
TpsaVector<ScalarT>::EntryProxy::operator typename TpsaVector<ScalarT>::EqVector() const
{
    EqVector res;
    for (std::size_t eqIdx = 0; eqIdx < numTpsaEq; ++eqIdx) {
        res[eqIdx] = at_(eqIdx);
    }

    return res;
}

template <class ScalarT>
typename TpsaVector<ScalarT>::EntryProxy&
TpsaVector<ScalarT>::EntryProxy::operator=(Scalar value)
{
    for (std::size_t eqIdx = 0; eqIdx < numTpsaEq; ++eqIdx) {
        at_(eqIdx) = value;
    }

    return *this;
}

template <class ScalarT>
typename TpsaVector<ScalarT>::EntryProxy&
TpsaVector<ScalarT>::EntryProxy::operator=(const EqVector& value)
{
    for (std::size_t eqIdx = 0; eqIdx < numTpsaEq; ++eqIdx) {
        at_(eqIdx) = value[eqIdx];
    }

    return *this;
}

template <class ScalarT>
typename TpsaVector<ScalarT>::EntryProxy&
TpsaVector<ScalarT>::EntryProxy::operator=(const EntryProxy& other)
{
    return *this = static_cast<EqVector>(other);
}

template <class ScalarT>
typename TpsaVector<ScalarT>::EntryProxy&
TpsaVector<ScalarT>::EntryProxy::operator+=(const EqVector& value)
{
    for (std::size_t eqIdx = 0; eqIdx < numTpsaEq; ++eqIdx) {
        at_(eqIdx) += value[eqIdx];
    }

    return *this;
}

template <class ScalarT>
typename TpsaVector<ScalarT>::EntryProxy&
TpsaVector<ScalarT>::EntryProxy::operator-=(const EqVector& value)
{
    for (std::size_t eqIdx = 0; eqIdx < numTpsaEq; ++eqIdx) {
        at_(eqIdx) -= value[eqIdx];
    }

    return *this;
}

template <class ScalarT>
typename TpsaVector<ScalarT>::EntryProxy&
TpsaVector<ScalarT>::EntryProxy::operator*=(Scalar factor)
{
    for (std::size_t eqIdx = 0; eqIdx < numTpsaEq; ++eqIdx) {
        at_(eqIdx) *= factor;
    }

    return *this;
}

template <class ScalarT>
ScalarT&
TpsaVector<ScalarT>::EntryProxy::at_(std::size_t eqIdx) const
{
    using namespace Dune::Indices;
    switch (eqIdx) {
    case 0:
        return (*v_)[_0][i_][0];
    case 1:
        return (*v_)[_1][i_][0];
    case 2:
        return (*v_)[_2][i_][0];
    case 3:
        return (*v_)[_3][i_][0];
    case 4:
        return (*v_)[_3][i_][1];
    case 5:
        return (*v_)[_3][i_][2];
    default:
        return (*v_)[_4][i_][0];
    }
}

//
// TpsaVector<ScalarT> implementation
//
template <class ScalarT>
void
TpsaVector<ScalarT>::resize(std::size_t numDof)
{
    using namespace Dune::Indices;
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              v_[fieldIdx].resize(numDof);
                          });
    size_ = numDof;
}

template <class ScalarT>
TpsaVector<ScalarT>&
TpsaVector<ScalarT>::operator=(Scalar value)
{
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              v_[fieldIdx] = value;
                          });

    return *this;
}

template <class ScalarT>
TpsaVector<ScalarT>&
TpsaVector<ScalarT>::operator+=(const TpsaVector& other)
{
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              v_[fieldIdx] += other.v_[fieldIdx];
                          });

    return *this;
}

template <class ScalarT>
TpsaVector<ScalarT>&
TpsaVector<ScalarT>::operator-=(const TpsaVector& other)
{
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              v_[fieldIdx] -= other.v_[fieldIdx];
                          });

    return *this;
}

template <class ScalarT>
TpsaVector<ScalarT>&
TpsaVector<ScalarT>::operator*=(Scalar factor)
{
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              v_[fieldIdx] *= factor;
                          });

    return *this;
}

template <class ScalarT>
void
TpsaVector<ScalarT>::scaleFields(const std::array<Scalar, numTpsaFields>& factors)
{
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              v_[fieldIdx] *= factors[fieldIdx];
                          });
}

template <class ScalarT>
ScalarT
TpsaVector<ScalarT>::one_norm() const
{
    Scalar norm = 0.0;
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              norm += v_[fieldIdx].one_norm();
                          });

    return norm;
}

template <class ScalarT>
ScalarT
TpsaVector<ScalarT>::two_norm2() const
{
    Scalar norm = 0.0;
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              norm += v_[fieldIdx].two_norm2();
                          });

    return norm;
}

template <class ScalarT>
ScalarT
TpsaVector<ScalarT>::infinity_norm() const
{
    Scalar norm = 0.0;
    Dune::Hybrid::forEach(Dune::range(Dune::index_constant<numTpsaFields>{}),
                          [&](auto fieldIdx) {
                              norm = std::max(norm, v_[fieldIdx].infinity_norm());
                          });

    return norm;
}

template <class ScalarT>
typename TpsaVector<ScalarT>::EqVector
TpsaVector<ScalarT>::operator[](std::size_t dofIdx) const
{
    using namespace Dune::Indices;
    EqVector res;
    res[0] = v_[_0][dofIdx][0];
    res[1] = v_[_1][dofIdx][0];
    res[2] = v_[_2][dofIdx][0];
    res[3] = v_[_3][dofIdx][0];
    res[4] = v_[_3][dofIdx][1];
    res[5] = v_[_3][dofIdx][2];
    res[6] = v_[_4][dofIdx][0];

    return res;
}

template class TpsaVector<double>;

#if FLOW_INSTANTIATE_FLOAT
template class TpsaVector<float>;
#endif

} // namespace Opm::Linear
