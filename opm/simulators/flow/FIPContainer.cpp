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

#include <opm/output/data/Solution.hpp>

#include <algorithm>

namespace Opm {

template<class FluidSystem>
bool
FIPContainer<FluidSystem>::
allocate(const std::size_t bufferSize,
         const SummaryConfig& summaryConfig,
         const bool forceAlloc,
         std::map<std::string, int>& rstKeywords)
{
    using namespace std::string_literals;

    const auto fipctrl = std::array {
        std::pair { "FIP"s , &OutputRestart::noPrefix  },
        std::pair { "SFIP"s, &OutputRestart::surface   },
        std::pair { "RFIP"s, &OutputRestart::reservoir },
    };

    this->outputRestart_.clearBits();

    for (const auto& [mnemonic, kind] : fipctrl) {
        if (auto fipPos = rstKeywords.find(mnemonic);
            fipPos != rstKeywords.end())
        {
            fipPos->second = 0;
            this->outputRestart_.*kind = true;
        }
    }

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
FIPContainer<FluidSystem>::
add(const Inplace::Phase phase)
{
    this->fip_[phase].resize(bufferSize_, 0.0);
}

template<class FluidSystem>
const std::vector<typename FIPContainer<FluidSystem>::Scalar>&
FIPContainer<FluidSystem>::
get(const Inplace::Phase phase) const
{
    return this->fip_.at(phase);
}

template<class FluidSystem>
bool
FIPContainer<FluidSystem>::
has(const Inplace::Phase phase) const
{
    const auto it = this->fip_.find(phase);
    return it != this->fip_.end() && !it->second.empty();
}

template<class FluidSystem>
bool
FIPContainer<FluidSystem>::
hasCo2InGas() const
{
    static const auto phases = std::array {
        Inplace::Phase::CO2InGasPhaseInMob,
        Inplace::Phase::CO2InGasPhaseMob,
        Inplace::Phase::CO2MassInGasPhaseInMob,
        Inplace::Phase::CO2MassInGasPhaseMob,
        Inplace::Phase::CO2Mass,
        Inplace::Phase::CO2MassInGasPhase,
        Inplace::Phase::CO2InGasPhaseInMobKrg,
        Inplace::Phase::CO2InGasPhaseMobKrg,
        Inplace::Phase::CO2MassInGasPhaseInMobKrg,
        Inplace::Phase::CO2MassInGasPhaseMobKrg,
        Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped,
        Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped,
        Inplace::Phase::CO2MassInGasPhaseMaximumTrapped,
        Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped,
    };

    return std::any_of(phases.begin(), phases.end(),
                       [this](const auto phase) { return has(phase); });
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignCo2InGas(const unsigned globalDofIdx, const Co2InGasInput& v)
{
    const Scalar massGas = (1.0 - v.xgW) * v.pv * v.rhog;
    if (this->has(Inplace::Phase::CO2Mass)) {
        this->fip_[Inplace::Phase::CO2Mass][globalDofIdx] = massGas * v.sg;
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhase)) {
        this->fip_[Inplace::Phase::CO2MassInGasPhase][globalDofIdx] = massGas * v.sg;
    }

    if (this->has(Inplace::Phase::CO2InGasPhaseInMob)) {
        const Scalar imMobileGas = massGas / v.mM * std::min(v.sgcr , v.sg);
        this->fip_[Inplace::Phase::CO2InGasPhaseInMob][globalDofIdx] = imMobileGas;
    }

    if (this->has(Inplace::Phase::CO2InGasPhaseMob)) {
        const Scalar mobileGas = massGas / v.mM * std::max(Scalar{0.0}, v.sg - v.sgcr);
        this->fip_[Inplace::Phase::CO2InGasPhaseMob][globalDofIdx] = mobileGas;
    }

    if (this->has(Inplace::Phase::CO2InGasPhaseInMobKrg)) {
        if (v.sgcr >= v.sg) {
            const Scalar imMobileGasKrg = massGas / v.mM * v.sg;
            this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg][globalDofIdx] = imMobileGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2InGasPhaseInMobKrg][globalDofIdx] = 0;
        }
    }

    if (this->has(Inplace::Phase::CO2InGasPhaseMobKrg)) {
        if (v.sg > v.sgcr) {
            const Scalar mobileGasKrg = massGas / v.mM * v.sg;
            this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg][globalDofIdx] = mobileGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2InGasPhaseMobKrg][globalDofIdx] = 0;
        }
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseInMob)) {
        const Scalar imMobileMassGas = massGas * std::min(v.sgcr , v.sg);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseInMob][globalDofIdx] = imMobileMassGas;
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseMob)) {
        const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, v.sg - v.sgcr);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseMob][globalDofIdx] = mobileMassGas;
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseInMobKrg)) {
        if (v.sgcr >= v.sg) {
            const Scalar imMobileMassGasKrg = massGas * v.sg;
            this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg][globalDofIdx] = imMobileMassGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2MassInGasPhaseInMobKrg][globalDofIdx] = 0;
        }
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseMobKrg)) {
        if (v.sg > v.sgcr) {
            const Scalar mobileMassGasKrg = massGas * v.sg;
            this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg][globalDofIdx] = mobileMassGasKrg;
        } else {
            this->fip_[Inplace::Phase::CO2MassInGasPhaseMobKrg][globalDofIdx] = 0;
        }
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseMaximumTrapped)) {
        const Scalar imMobileMassGas = massGas * std::min(v.trappedGas, v.sg);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumTrapped][globalDofIdx] = imMobileMassGas;
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped)) {
        const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, v.sg - v.trappedGas);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped][globalDofIdx] = mobileMassGas;
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped)) {
        const Scalar imMobileMassGas = massGas * std::min(v.strandedGas, v.sg);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped][globalDofIdx] = imMobileMassGas;
    }

    if (this->has(Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped)) {
        const Scalar mobileMassGas = massGas * std::max(Scalar{0.0}, v.sg - v.strandedGas);
        this->fip_[Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped][globalDofIdx] = mobileMassGas;
    }
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignGasWater(const unsigned  globalDofIdx,
               const std::array<Scalar, numPhases>& fip,
               const Scalar    gasInPlaceWater,
               const Scalar    waterInPlaceGas)
{
    if (this->has(Inplace::Phase::WaterInGasPhase)) {
        this->fip_[Inplace::Phase::WaterInGasPhase][globalDofIdx] = waterInPlaceGas;
    }

    if (this->has(Inplace::Phase::WaterInWaterPhase)) {
        this->fip_[Inplace::Phase::WaterInWaterPhase][globalDofIdx] = fip[waterPhaseIdx];
    }

    // For water+gas cases the gas in water is added to the GIPL value
    if (this->has(Inplace::Phase::GasInLiquidPhase) && !FluidSystem::phaseIsActive(oilPhaseIdx)) {
        this->fip_[Inplace::Phase::GasInLiquidPhase][globalDofIdx] = gasInPlaceWater;
    }

    // Add dissolved gas and vaporized water to total Fip
    if (this->has(Inplace::Phase::WATER)) {
        this->fip_[Inplace::Phase::WATER][globalDofIdx] += waterInPlaceGas;
    }

    if (this->has(Inplace::Phase::GAS)) {
        this->fip_[Inplace::Phase::GAS][globalDofIdx] += gasInPlaceWater;
    }
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignVolumesSurface(const unsigned globalDofIdx,
                     const std::array<Scalar, numPhases>& fip)
{
    if (FluidSystem::phaseIsActive(oilPhaseIdx) && this->has(Inplace::Phase::OIL)) {
        this->fip_[Inplace::Phase::OIL][globalDofIdx] = fip[oilPhaseIdx];
    }

    if (FluidSystem::phaseIsActive(oilPhaseIdx) && this->has(Inplace::Phase::OilInLiquidPhase)) {
        this->fip_[Inplace::Phase::OilInLiquidPhase][globalDofIdx] = fip[oilPhaseIdx];
    }

    if (FluidSystem::phaseIsActive(gasPhaseIdx) && this->has(Inplace::Phase::GAS)) {
        this->fip_[Inplace::Phase::GAS][globalDofIdx] = fip[gasPhaseIdx];
    }

    if (FluidSystem::phaseIsActive(gasPhaseIdx) && this->has(Inplace::Phase::GasInGasPhase)) {
        this->fip_[Inplace::Phase::GasInGasPhase][globalDofIdx] = fip[gasPhaseIdx];
    }

    if (FluidSystem::phaseIsActive(waterPhaseIdx) && this->has(Inplace::Phase::WATER)) {
        this->fip_[Inplace::Phase::WATER][globalDofIdx] = fip[waterPhaseIdx];
    }
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignVolumesReservoir(const unsigned    globalDofIdx,
                       const Scalar      saltConcentration,
                       const std::array<Scalar, numPhases>& fipr)
{
    if (FluidSystem::phaseIsActive(oilPhaseIdx) && this->has(Inplace::Phase::OilResVolume)) {
        this->fip_[Inplace::Phase::OilResVolume][globalDofIdx] = fipr[oilPhaseIdx];
    }

    if (FluidSystem::phaseIsActive(gasPhaseIdx) && this->has(Inplace::Phase::GasResVolume)) {
        this->fip_[Inplace::Phase::GasResVolume][globalDofIdx] = fipr[gasPhaseIdx];
    }

    if (FluidSystem::phaseIsActive(waterPhaseIdx) && this->has(Inplace::Phase::WaterResVolume)) {
        this->fip_[Inplace::Phase::WaterResVolume][globalDofIdx] = fipr[waterPhaseIdx];
    }

    if (FluidSystem::phaseIsActive(waterPhaseIdx) && this->has(Inplace::Phase::SALT)) {
        this->fip_[Inplace::Phase::SALT][globalDofIdx] =
            fipr[waterPhaseIdx] * saltConcentration;
    }
}

template<class FluidSystem>
bool
FIPContainer<FluidSystem>::
hasCo2InWater() const
{
    static const auto phases = std::array {
        Inplace::Phase::CO2InWaterPhase,
        Inplace::Phase::CO2MassInWaterPhase,
        Inplace::Phase::CO2Mass,
    };

    return std::any_of(phases.begin(), phases.end(),
                       [this](const auto phase) { return has(phase); });
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignCo2InWater(const unsigned globalDofIdx,
                 const Scalar   co2InWater,
                 const Scalar   mM)
{
    if (this->has(Inplace::Phase::CO2Mass)) {
        this->fip_[Inplace::Phase::CO2Mass][globalDofIdx] += co2InWater  * mM;
    }

    if (this->has(Inplace::Phase::CO2MassInWaterPhase)) {
        this->fip_[Inplace::Phase::CO2MassInWaterPhase][globalDofIdx] = co2InWater  * mM;
    }

    if (this->has(Inplace::Phase::CO2InWaterPhase)) {
        this->fip_[Inplace::Phase::CO2InWaterPhase][globalDofIdx] = co2InWater;
    }
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignOilGasDistribution(const unsigned globalDofIdx,
                         const Scalar   gasInPlaceLiquid,
                         const Scalar   oilInPlaceGas)
{
    if (this->has(Inplace::Phase::GasInLiquidPhase)) {
        this->fip_[Inplace::Phase::GasInLiquidPhase][globalDofIdx] = gasInPlaceLiquid;
    }

    if (this->has(Inplace::Phase::OilInGasPhase)) {
        this->fip_[Inplace::Phase::OilInGasPhase][globalDofIdx] = oilInPlaceGas;
    }

    // Add dissolved gas and vaporized oil to total Fip
    if (this->has(Inplace::Phase::OIL)) {
        this->fip_[Inplace::Phase::OIL][globalDofIdx] += oilInPlaceGas;
    }

    if (this->has(Inplace::Phase::GAS)) {
        this->fip_[Inplace::Phase::GAS][globalDofIdx] += gasInPlaceLiquid;
    }
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
outputRestart(data::Solution& sol)
{
    if (!this->outputRestart_) {
        return;
    }

    using namespace std::string_literals;

    using M = UnitSystem::measure;
    using FIPEntry = std::tuple<std::string, M, Inplace::Phase>;

    auto fipArrays = std::vector<FIPEntry> {};
    if (this->outputRestart_.surface) {
        fipArrays.insert(fipArrays.end(), {
                FIPEntry {"SFIPOIL"s, M::liquid_surface_volume, Inplace::Phase::OIL   },
                FIPEntry {"SFIPWAT"s, M::liquid_surface_volume, Inplace::Phase::WATER },
                FIPEntry {"SFIPGAS"s, M::gas_surface_volume,    Inplace::Phase::GAS   },
            });
    }

    if (this->outputRestart_.reservoir) {
        fipArrays.insert(fipArrays.end(), {
                FIPEntry {"RFIPOIL"s, M::volume, Inplace::Phase::OilResVolume   },
                FIPEntry {"RFIPWAT"s, M::volume, Inplace::Phase::WaterResVolume },
                FIPEntry {"RFIPGAS"s, M::volume, Inplace::Phase::GasResVolume   },
            });
    }

    if (this->outputRestart_.noPrefix && !this->outputRestart_.surface) {
        fipArrays.insert(fipArrays.end(), {
                FIPEntry { "FIPOIL"s, M::liquid_surface_volume, Inplace::Phase::OIL   },
                FIPEntry { "FIPWAT"s, M::liquid_surface_volume, Inplace::Phase::WATER },
                FIPEntry { "FIPGAS"s, M::gas_surface_volume,    Inplace::Phase::GAS   },
            });
    }

    for (const auto& [mnemonic, unit, phase] : fipArrays) {
        if (! this->fip_[phase].empty()) {
            sol.insert(mnemonic, unit, std::move(this->fip_[phase]),
                       data::TargetType::RESTART_SOLUTION);
        }
    }

    for (const auto& phase : Inplace::mixingPhases()) {
        if (! this->fip_[phase].empty()) {
            sol.insert(Inplace::EclString(phase),
                       UnitSystem::measure::volume,
                       std::move(this->fip_[phase]),
                       data::TargetType::SUMMARY);
        }
    }
}

template<class FluidSystem>
void
FIPContainer<FluidSystem>::
assignPoreVolume(const unsigned globalDofIdx,
                 const Scalar   value)
{
    this->fip_[Inplace::Phase::PoreVolume][globalDofIdx] = value;
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
