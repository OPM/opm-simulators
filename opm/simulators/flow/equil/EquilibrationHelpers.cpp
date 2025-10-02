// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/simulators/flow/equil/EquilibrationHelpers_impl.hpp>

namespace Opm {
namespace EQUIL {

template<class Scalar> using FS = BlackOilFluidSystem<Scalar>;

#define INSTANTIATE_TYPE1(T, MatLaw) \
    template struct PcEq<FS<T>,MatLaw<T>>; \
    template T satFromPc<FS<T>,MatLaw<T>>(const MatLaw<T>&, \
                                          const int,const int, \
                                          const T,const bool); \
    template T satFromSumOfPcs<FS<T>,MatLaw<T>>(const MatLaw<T>&,    \
                                                const int,const int, \
                                                const int,const T); \
    template T satFromDepth<FS<T>,MatLaw<T>>(const MatLaw<T>&, \
                                             const T,const T, \
                                             const int,const int,const bool); \
    template bool isConstPc<FS<T>,MatLaw<T>>(const MatLaw<T>&,const int,const int);

#define INSTANTIATE_TYPE2(T) \
    template class EquilReg<T>; \
    template class Miscibility::PBVD<FS<T>>; \
    template class Miscibility::PDVD<FS<T>>; \
    template class Miscibility::RsVD<FS<T>>; \
    template class Miscibility::RsSatAtContact<FS<T>>; \
    template class Miscibility::RvSatAtContact<FS<T>>; \
    template class Miscibility::RvwSatAtContact<FS<T>>; \
    template class Miscibility::RvVD<FS<T>>; \
    template class Miscibility::RvwVD<FS<T>>;

#define INSTANTIATE_TYPE(T, ML1, ML2) \
INSTANTIATE_TYPE1(T, ML1) \
INSTANTIATE_TYPE1(T, ML2) \
INSTANTIATE_TYPE2(T)

template<class Scalar>
using MatLaw1 = EclMaterialLaw::Manager<ThreePhaseMaterialTraits<Scalar,0,1,2,true,true>>;
template<class Scalar>
using MatLaw2 = EclMaterialLaw::Manager<ThreePhaseMaterialTraits<Scalar,0,1,2,false,true>>;

INSTANTIATE_TYPE(double, MatLaw1, MatLaw2)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float, MatLaw1, MatLaw2)
#endif

} // namespace Equil
} // namespace Opm
