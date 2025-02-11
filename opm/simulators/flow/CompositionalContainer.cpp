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
#include <opm/simulators/flow/CompositionalContainer.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasFluidSystem.hpp>

namespace Opm {

template<class FluidSystem>
void CompositionalContainer<FluidSystem>::
allocate(const unsigned bufferSize,
         std::map<std::string, int>& rstKeywords)
{
    if (auto& zmf = rstKeywords["ZMF"]; zmf > 0) {
        zmf = 0;
        this->allocated_ = true;
        for (int i = 0; i < numComponents; ++i) {
            moleFractions_[i].resize(bufferSize, 0.0);
        }
    }

    if (auto& xmf = rstKeywords["XMF"]; xmf > 0 && FluidSystem::phaseIsActive(oilPhaseIdx)) {
        this->allocated_ = true;
        xmf = 0;
        for (int i = 0; i < numComponents; ++i) {
            phaseMoleFractions_[oilPhaseIdx][i].resize(bufferSize, 0.0);
        }
    }

    if (auto& ymf = rstKeywords["YMF"]; ymf > 0 && FluidSystem::phaseIsActive(gasPhaseIdx)) {
        this->allocated_ = true;
        ymf = 0;
        for (int i = 0; i < numComponents; ++i) {
            phaseMoleFractions_[gasPhaseIdx][i].resize(bufferSize, 0.0);
        }
    }
}

template<class T> using FS = BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>;

#define INSTANTIATE_TYPE(T) \
    template class CompositionalContainer<FS<T>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#define INSTANTIATE_COMP(NUM) \
    template<class T> using FS##NUM = GenericOilGasFluidSystem<T, NUM>; \
    template class CompositionalContainer<FS##NUM<double>>;

INSTANTIATE_COMP(0)
INSTANTIATE_COMP(2)
INSTANTIATE_COMP(3)
INSTANTIATE_COMP(4)
INSTANTIATE_COMP(5)
INSTANTIATE_COMP(6)
INSTANTIATE_COMP(7)

} // namespace Opm
