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

#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <tuple>

#include <fmt/format.h>

namespace Opm {

template<class FluidSystem>
void CompositionalContainer<FluidSystem>::
allocate(const unsigned bufferSize,
         std::map<std::string, int>& rstKeywords)
{
    if (auto& zmf = rstKeywords["ZMF"]; zmf > 0) {
        this->allocated_ = true;
        zmf = 0;
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

template<class FluidSystem>
void CompositionalContainer<FluidSystem>::
assignGasFractions(const unsigned globalDofIdx,
                   const AssignFunction& fractions)
{
    if (phaseMoleFractions_[gasPhaseIdx][0].empty()) {
        return;
    }

    std::for_each(phaseMoleFractions_[gasPhaseIdx].begin(),
                  phaseMoleFractions_[gasPhaseIdx].end(),
                  [globalDofIdx, &fractions, c = 0](auto& comp) mutable
                  { comp[globalDofIdx] = fractions(c++); });
}

template<class FluidSystem>
void CompositionalContainer<FluidSystem>::
assignMoleFractions(const unsigned globalDofIdx,
                    const AssignFunction& fractions)
{
    if (moleFractions_.empty()) {
        return;
    }

    std::for_each(moleFractions_.begin(), moleFractions_.end(),
                  [&fractions, globalDofIdx, c = 0](auto& comp) mutable
                  { comp[globalDofIdx] = fractions(c++); });
}

template<class FluidSystem>
void CompositionalContainer<FluidSystem>::
assignOilFractions(const unsigned globalDofIdx,
                   const AssignFunction& fractions)
{
    if (phaseMoleFractions_[oilPhaseIdx][0].empty()) {
        return;
    }

    std::for_each(phaseMoleFractions_[oilPhaseIdx].begin(),
                  phaseMoleFractions_[oilPhaseIdx].end(),
                  [globalDofIdx, &fractions, c = 0](auto& comp) mutable
                  { comp[globalDofIdx] = fractions(c++); });
}

template<class FluidSystem>
void CompositionalContainer<FluidSystem>::
outputRestart(data::Solution& sol,
              ScalarBuffer& oil_saturation)
{
    using DataEntry =
        std::tuple<std::string, UnitSystem::measure, std::vector<Scalar>&>;

    auto doInsert = [&sol](DataEntry&       entry,
                           const data::TargetType target)
    {
        if (std::get<2>(entry).empty()) {
            return;
        }

        sol.insert(std::get<std::string>(entry),
                   std::get<UnitSystem::measure>(entry),
                   std::move(std::get<2>(entry)),
                   target);
    };

    auto entries = std::vector<DataEntry>{};

    // ZMF
    if (!moleFractions_[0].empty()) {
        for (int i = 0; i < numComponents; ++i) {
            const auto name = fmt::format("ZMF{}", i + 1);  // Generate ZMF1, ZMF2, ...
            entries.emplace_back(name, UnitSystem::measure::identity, moleFractions_[i]);
        }
    }

    // XMF
    if (!phaseMoleFractions_[oilPhaseIdx][0].empty()) {
        for (int i = 0; i < numComponents; ++i) {
            const auto name = fmt::format("XMF{}", i + 1);  // Generate XMF1, XMF2, ...
            entries.emplace_back(name, UnitSystem::measure::identity,
                                 phaseMoleFractions_[oilPhaseIdx][i]);
        }
    }

    // YMF
    if (!phaseMoleFractions_[gasPhaseIdx][0].empty()) {
        for (int i = 0; i < numComponents; ++i) {
            const auto name = fmt::format("YMF{}", i + 1);  // Generate YMF1, YMF2, ...
            entries.emplace_back(name, UnitSystem::measure::identity,
                                 phaseMoleFractions_[gasPhaseIdx][i]);
        }
    }

    if (!oil_saturation.empty()) {
        entries.emplace_back("SOIL", UnitSystem::measure::identity, oil_saturation);
    }

    std::for_each(entries.begin(), entries.end(),
                  [&doInsert](auto& array)
                  { doInsert(array, data::TargetType::RESTART_SOLUTION); });

    this->allocated_ = false;
}

#define INSTANTIATE_COMP_THREEPHASE(NUM) \
    template<class T> using FS##NUM = GenericOilGasWaterFluidSystem<T, NUM, true>; \
    template class CompositionalContainer<FS##NUM<double>>;

#define INSTANTIATE_COMP_TWOPHASE(NUM) \
    template<class T> using GFS##NUM = GenericOilGasWaterFluidSystem<T, NUM, false>; \
    template class CompositionalContainer<GFS##NUM<double>>;

#define INSTANTIATE_COMP(NUM) \
    INSTANTIATE_COMP_THREEPHASE(NUM) \
    INSTANTIATE_COMP_TWOPHASE(NUM)

INSTANTIATE_COMP_THREEPHASE(0)  // \Note: to register the parameter ForceDisableFluidInPlaceOutput
INSTANTIATE_COMP(2)
INSTANTIATE_COMP(3)
INSTANTIATE_COMP(4)
INSTANTIATE_COMP(5)
INSTANTIATE_COMP(6)
INSTANTIATE_COMP(7)

} // namespace Opm
