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
#ifndef OPM_COMPOSITIONAL_CONTAINER_HPP
#define OPM_COMPOSITIONAL_CONTAINER_HPP

#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace Opm {

namespace data { class Solution; }

template<class FluidSystem>
class CompositionalContainer
{
    using Scalar = typename FluidSystem::Scalar;
    using ScalarBuffer = std::vector<Scalar>;

    static constexpr int numComponents = FluidSystem::numComponents;

    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;

public:
    void allocate(const unsigned bufferSize,
                  std::map<std::string, int>& rstKeywords);

    using AssignFunction = std::function<Scalar(const unsigned)>;

    void assignGasFractions(const unsigned globalDofIdx,
                            const AssignFunction& fractions);

    void assignMoleFractions(const unsigned globalDofIdx,
                             const AssignFunction& fractions);

    void assignOilFractions(const unsigned globalDofIdx,
                            const AssignFunction& fractions);

    void outputRestart(data::Solution& sol,
                       ScalarBuffer& oil_saturation);

    bool allocated() const
    { return allocated_; }

private:
    bool allocated_ = false;
    // total mole fractions for each component
    std::array<ScalarBuffer, numComponents> moleFractions_;
    // mole fractions for each component in each phase
    std::array<std::array<ScalarBuffer, numComponents>, numPhases> phaseMoleFractions_;
};

} // namespace Opm

#endif // OPM_COMPOSITIONAL_CONTAINER_HPP
