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
 * \brief Restart-output buffers specific to compositional simulations.
 */
#ifndef OPM_COMPOSITIONAL_CONTAINER_HPP
#define OPM_COMPOSITIONAL_CONTAINER_HPP

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <array>
#include <functional>
#include <map>
#include <optional>
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
    enum class RestartOutput { Disabled, Enabled };

    /// Allocate compositional restart fields requested by \p rstKeywords.
    /// PSAT is currently consumed only by restart output. Allocate it only
    /// for passes that will write restart data because filling the buffer
    /// requires a nonlinear solve for each single-phase cell.
    void allocate(const unsigned bufferSize,
                  std::map<std::string, int>& rstKeywords,
                  RestartOutput restartOutput);

    using AssignFunction = std::function<Scalar(const unsigned)>;

    void assignGasFractions(const unsigned globalDofIdx,
                            const AssignFunction& fractions);

    void assignMoleFractions(const unsigned globalDofIdx,
                             const AssignFunction& fractions);

    void assignOilFractions(const unsigned globalDofIdx,
                            const AssignFunction& fractions);

    void assignPhasePressures(const unsigned globalDofIdx,
                              const Scalar oilPressure,
                              const Scalar gasPressure);

    void assignSaturationPressure(const unsigned globalDofIdx,
                                  const Scalar psat);

    /// Return the oil-phase pressure for two hydrocarbon phases, or the
    /// bubble/dew pressure of the total composition for a single phase. Return
    /// std::nullopt if the solver cannot determine a saturation pressure.
    ///
    /// \p liquidFraction is the flash liquid fraction L: exactly one denotes
    /// liquid only, exactly zero vapour only, and other values denote two phases.
    /// Use L because computed phase saturations can contain round-off residuals.
    /// Plain-value arguments allow testing phase selection without a simulator.
    [[nodiscard]] static std::optional<Scalar>
    cellSaturationPressure(const Scalar liquidFraction,
                           const Scalar oilPressure,
                           const std::array<Scalar, numComponents>& moleFractions,
                           const Scalar temperature,
                           const CompositionalConfig::EOSType eosType);

    void assignVaporFraction(const unsigned globalDofIdx,
                             const Scalar vmf);

    void outputRestart(data::Solution& sol,
                       ScalarBuffer& oil_saturation);

    bool moleFractionsAllocated() const
    { return !moleFractions_[0].empty(); }

    bool gasFractionsAllocated() const
    { return !phaseMoleFractions_[gasPhaseIdx][0].empty(); }

    bool oilFractionsAllocated() const
    { return !phaseMoleFractions_[oilPhaseIdx][0].empty(); }

    bool phasePressuresAllocated() const
    { return !oilPressure_.empty() || !gasPressure_.empty(); }

    bool vaporFractionAllocated() const
    { return !vaporFraction_.empty(); }

    bool saturationPressureAllocated() const
    { return !saturationPressure_.empty(); }

    /// Whether PSAT was requested in this allocation pass. Unlike the buffer
    /// state above, this remains true on MPI ranks with no local cells.
    bool saturationPressureRequested() const
    { return saturationPressureRequested_; }

    bool allocated() const
    { return allocated_; }

private:
    bool allocated_ = false;
    // total mole fractions for each component
    std::array<ScalarBuffer, numComponents> moleFractions_;
    // mole fractions for each component in each phase
    std::array<std::array<ScalarBuffer, numComponents>, numPhases> phaseMoleFractions_;
    // phase pressures (POIL, PGAS)
    ScalarBuffer oilPressure_;
    ScalarBuffer gasPressure_;
    // saturation pressure (PSAT)
    ScalarBuffer saturationPressure_;
    bool saturationPressureRequested_ = false;
    // vapour mole fraction of the total mixture (VMF)
    ScalarBuffer vaporFraction_;
};

} // namespace Opm

#endif // OPM_COMPOSITIONAL_CONTAINER_HPP
