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

#include "config.h"
#include "equilibrationhelpers_impl.hh"

namespace Opm {
namespace EQUIL {
    
using MatLaw = EclMaterialLawManager<ThreePhaseMaterialTraits<double,0,1,2>>;
using FS = BlackOilFluidSystem<double>;

template struct PcEq<FS,MatLaw>;

template double satFromPc<FS,MatLaw>(const MatLaw&,const int,const int,
                                     const double,const bool);
template double satFromSumOfPcs<FS,MatLaw>(const MatLaw&,const int,const int,
                                           const int,const double);
template double satFromDepth<FS,MatLaw>(const MatLaw&,const double,const double,
                                        const int,const int,const bool);
template bool isConstPc<FS,MatLaw>(const MatLaw&,const int,const int);

namespace Miscibility {
    template class PBVD<FS>;
    template class PDVD<FS>;
    template class RsVD<FS>;
    template class RsSatAtContact<FS>;
    template class RvSatAtContact<FS>;
    template class RvwSatAtContact<FS>;
    template class RvVD<FS>;
    template class RvwVD<FS>;
}

} // namespace Equil
} // namespace Opm
