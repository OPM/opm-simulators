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
/*!
 * \file
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_FIP_CONTAINER_HPP
#define OPM_FIP_CONTAINER_HPP

#include <opm/output/eclipse/Inplace.hpp>

#include <array>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>

namespace Opm::data {
class Solution;
}

namespace Opm {

class SummaryConfig;

template<class FluidSystem>
class FIPContainer {
public:
    using Scalar = typename FluidSystem::Scalar;
    using FIPMap = std::unordered_map<Inplace::Phase, std::vector<Scalar>>;

    static constexpr auto numPhases = FluidSystem::numPhases;
    static constexpr auto gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr auto oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

    bool allocate(const std::size_t bufferSize,
                  const SummaryConfig& summaryConfig,
                  const bool forceAlloc,
                  std::map<std::string, int>& rstKeywords);

    void add(const Inplace::Phase phase);

    struct Co2InGasInput
    {
        double pv;
        Scalar sg;
        Scalar sgcr;
        Scalar rhog;
        Scalar xgW;
        Scalar mM;
        Scalar trappedGas;
        Scalar strandedGas;
    };

    const std::vector<Scalar>& get(const Inplace::Phase phase) const;

    bool has(const Inplace::Phase phase) const;

    bool hasCo2InGas() const;
    void assignCo2InGas(const unsigned globalDofIdx,
                      const Co2InGasInput& v);

    bool hasCo2InWater() const;
    void assignCo2InWater(const unsigned globalDofIdx,
                          const Scalar   co2InWater,
                          const Scalar   mM);

    void assignGasWater(const unsigned  globalDofIdx,
                        const std::array<Scalar, numPhases>& fip,
                        const Scalar    gasInPlaceWater,
                        const Scalar    waterInPlaceGas);

    void assignOilGasDistribution(const unsigned globalDofIdx,
                                  const Scalar   gasInPlaceLiquid,
                                  const Scalar   oilInPlaceGas);

    void assignPoreVolume(const unsigned globalDofIdx,
                          const Scalar   value);

    void assignVolumesSurface(const unsigned globalDofIdx,
                              const std::array<Scalar, numPhases>& fip);

    void assignVolumesReservoir(const unsigned    globalDofIdx,
                                const Scalar      saltConcentration,
                                const std::array<Scalar, numPhases>& fipr);

    void outputRestart(data::Solution& sol);

private:
    FIPMap fip_{};
    std::size_t bufferSize_ = 0;

    struct OutputRestart
    {
        /// Whether or not run requests (surface condition) fluid-in-place
        /// restart file output using the 'FIP' mnemonic.
        bool noPrefix {false};

        /// Whether or not run requests surface condition fluid-in-place
        /// restart file output using the 'SFIP' mnemonic.
        bool surface {false};

        /// Whether or not run requests reservoir condition fluid-in-place
        /// restart file output using the 'RFIP' mnemonic.
        bool reservoir {false};

        void clearBits()
        {
            this->noPrefix = this->surface = this->reservoir = false;
        }

        explicit operator bool() const
        {
            return this->noPrefix || this->surface || this->reservoir;
        }
    } outputRestart_{};
};

} // namespace Opm

#endif // OPM_FIP_CONTAINER_HPP
