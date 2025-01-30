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
#include <unordered_map>
#include <vector>

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

    // Temporary constructor until we are ready to own the map
    explicit FIPContainer(FIPMap& fip)
        : fip_(fip) {}

    bool allocate(const std::size_t bufferSize,
                  const SummaryConfig& summaryConfig,
                  const bool forceAlloc);

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

    void assignCo2InGas(const unsigned globalDofIdx,
                      const Co2InGasInput& v);

    void assignGasWater(const unsigned  globalDofIdx,
                        const std::array<Scalar, numPhases>& fip,
                        const Scalar    gasInPlaceWater,
                        const Scalar    waterInPlaceGas);

    void assignVolumesSurface(const unsigned globalDofIdx,
                              const std::array<Scalar, numPhases>& fip);

private:
    FIPMap& fip_;
    std::size_t bufferSize_ = 0;
};

} // namespace Opm

#endif // OPM_FIP_CONTAINER_HPP
