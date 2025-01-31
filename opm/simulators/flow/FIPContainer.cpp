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
#include <opm/simulators/flow/FIPContainer.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/GenericOilGasFluidSystem.hpp>

namespace Opm {

template<class FluidSystem>
bool
FIPContainer<FluidSystem>::
allocate(const std::size_t bufferSize,
         const SummaryConfig& summaryConfig,
         const bool forceAlloc)
{
    bool computeFip = false;
    bufferSize_ = bufferSize;
    for (const auto& phase : Inplace::phases()) {
        if (forceAlloc || summaryConfig.require3DField(Inplace::EclString(phase))) {
            this->add(phase);
            computeFip = true;
        }
        else {
            this->fip_[phase].clear();
        }
    }

    return computeFip;
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::add(const Inplace::Phase phase)
{
    this->fip_[phase].resize(bufferSize_, 0.0);
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignCo2InGas(const unsigned globalDofIdx, const Co2InGasInput& v)
{
    const Scalar massGas = (1.0 - v.xgW) * v.pv * v.rhog;
    if (!this->fip_[Inplace::Phase::CO2Mass].empty()) {
        this->fip_[Inplace::Phase::CO2Mass][globalDofIdx] = massGas * v.sg;
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhase].empty()) {
        this->fip_[Inplace::Phase::CO2MassInGasPhase][globalDofIdx] = massGas * v.sg;
    }

    if (!this->fip_[Inplace::Phase::CO2InGasPhaseInMob].empty()) {
        const Scalar imMobileGas = massGas / v.mM * std::min(v.sgcr , v.sg);
        this->fip_[Inplace::Phase::CO2InGasPhaseInMob][globalDofIdx] = imMobileGas;
    }

    if (!this->fip_[Inplace::Phase::CO2InGasPhaseMob].empty()) {
        const Scalar mobileGas = massGas / v.mM * std::max(Scalar{0.0}, v.sg - v.sgcr);
        this->fip_[Inplace::Phase::CO2InGasPhaseMob][globalDofIdx] = mobileGas;
    }

    if (!this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg].empty()) {
        if (v.sgcr >= v.sg) {
            const Scalar imMobileGasKrg = massGas / v.mM * v.sg;
            this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg][globalDofIdx] = imMobileGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg][globalDofIdx] = 0;
        }
    }

    if (!this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg].empty()) {
        if (v.sg > v.sgcr) {
            const Scalar mobileGasKrg = massGas / v.mM * v.sg;
            this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg][globalDofIdx] = mobileGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg][globalDofIdx] = 0;
        }
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseInMob].empty()) {
        const Scalar imMobileMassGas = massGas * std::min(v.sgcr , v.sg);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseInMob][globalDofIdx] = imMobileMassGas;
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMob].empty()) {
        const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, v.sg - v.sgcr);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseMob][globalDofIdx] = mobileMassGas;
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg].empty()) {
        if (v.sgcr >= v.sg) {
            const Scalar imMobileMassGasKrg = massGas * v.sg;
            this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg][globalDofIdx] = imMobileMassGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg][globalDofIdx] = 0;
        }
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg].empty()) {
        if (v.sg > v.sgcr) {
            const Scalar mobileMassGasKrg = massGas * v.sg;
            this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg][globalDofIdx] = mobileMassGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg][globalDofIdx] = 0;
        }
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumTrapped].empty()) {
        const Scalar imMobileMassGas = massGas * std::min(v.trappedGas, v.sg);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumTrapped][globalDofIdx] = imMobileMassGas;
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped].empty()) {
        const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, v.sg - v.trappedGas);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped][globalDofIdx] = mobileMassGas;
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped].empty()) {
        const Scalar imMobileMassGas = massGas * std::min(v.strandedGas, v.sg);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped][globalDofIdx] = imMobileMassGas;
    }

    if (!this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped].empty()) {
        const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, v.sg - v.strandedGas);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped][globalDofIdx] = mobileMassGas;
    }
}

template<class T> using FS = BlackOilFluidSystem<T,BlackOilDefaultIndexTraits>;

#define INSTANTIATE_TYPE(T) \
    template class FIPContainer<FS<T>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#define INSTANTIATE_COMP(NUM) \
    template<class T> using FS##NUM = GenericOilGasFluidSystem<T, NUM>; \
    template class FIPContainer<FS##NUM<double>>;

INSTANTIATE_COMP(0)
INSTANTIATE_COMP(2)
INSTANTIATE_COMP(3)
INSTANTIATE_COMP(4)
INSTANTIATE_COMP(5)
INSTANTIATE_COMP(6)
INSTANTIATE_COMP(7)

} // namespace Opm
