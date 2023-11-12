
/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>


namespace Opm
{

PhaseUsage phaseUsage(const Phases& phases)
{
    PhaseUsage pu;
    pu.phase_used.fill(0);

    // Discover phase usage.
    pu.phase_used[BlackoilPhases::Aqua] = phases.active(Phase::WATER);
    pu.phase_used[BlackoilPhases::Liquid] = phases.active(Phase::OIL);
    pu.phase_used[BlackoilPhases::Vapour] = phases.active(Phase::GAS);

    pu.num_phases = 0;
    int activePhaseIdx = -1;
    for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
        if (!pu.phase_used[phaseIdx]) {
            pu.phase_pos[phaseIdx] = -1;
        }
        else {
            pu.phase_pos[phaseIdx] = ++activePhaseIdx;
            pu.num_phases = activePhaseIdx+1;
        }
    }

    // Add solvent info
    pu.has_solvent = phases.active(Phase::SOLVENT);
    if (pu.has_solvent) {
        // this is quite a hack: even though solvent is not considered as in
        // MaxNumPhases and pu.num_phases because this would break a lot of
        // assumptions in old code, it is nevertheless an index to be translated
        // to. solvent and solvent are even larger hacks because not even this can be
        // done for them.
        pu.phase_pos[BlackoilPhases::Solvent] = ++activePhaseIdx;
    }
    else
        pu.phase_pos[BlackoilPhases::Solvent] = -1;

    // Add polymer info
    pu.has_polymer = phases.active(Phase::POLYMER);
    if (pu.has_polymer) {
        // this is quite a hack: even though polymer is not considered as in
        // MaxNumPhases and pu.num_phases because this would break a lot of
        // assumptions in old code, it is nevertheless an index to be translated
        // to. polymer and solvent are even larger hacks because not even this can be
        // done for them.
        pu.phase_pos[BlackoilPhases::Polymer] = ++activePhaseIdx;
    }
    else
        pu.phase_pos[BlackoilPhases::Polymer] = -1;

    // Add energy info
    pu.has_energy = phases.active(Phase::ENERGY);
    if (pu.has_energy) {
        // this is quite a hack: even though energy is not considered as in
        // MaxNumPhases and pu.num_phases because this would break a lot of
        // assumptions in old code, it is nevertheless an index to be translated
        // to. polymer and solvent are even larger hacks because not even this can be
        // done for them.
        pu.phase_pos[BlackoilPhases::Energy] = ++activePhaseIdx;
    }
    else
        pu.phase_pos[BlackoilPhases::Energy] = -1;

    // Add polymer molecular weight related
    pu.has_polymermw = phases.active(Phase::POLYMW);
    if (pu.has_polymermw) {
        if (!pu.has_polymer) {
            OPM_THROW(std::runtime_error, "pu.has_polymermw is true while pu.has_polymer is false");
        }
        pu.phase_pos[BlackoilPhases::PolymerMW] = ++activePhaseIdx;
    }
    else
        pu.phase_pos[BlackoilPhases::PolymerMW] = -1;

    // Add foam info
    pu.has_foam = phases.active(Phase::FOAM);
    if (pu.has_foam) {
        pu.phase_pos[BlackoilPhases::Foam] = ++activePhaseIdx;
    }
    else
        pu.phase_pos[BlackoilPhases::Foam] = -1;

    // Add brine info
    pu.has_brine = phases.active(Phase::BRINE);
    if (pu.has_brine) {
        pu.phase_pos[BlackoilPhases::Brine] = ++activePhaseIdx;
    }
    else
        pu.phase_pos[BlackoilPhases::Brine] = -1;

    // Add zFraction info
    pu.has_zFraction = phases.active(Phase::ZFRACTION);
    if (pu.has_zFraction) {
        pu.phase_pos[BlackoilPhases::ZFraction] = ++activePhaseIdx;
    }
    else
        pu.phase_pos[BlackoilPhases::ZFraction] = -1;


    return pu;
}

PhaseUsage phaseUsageFromDeck(const EclipseState& eclipseState)
{
    const auto& phases = eclipseState.runspec().phases();

    return phaseUsage(phases);
}

/// Looks at presence of WATER, OIL and GAS keywords in deck
/// to determine active phases.
PhaseUsage phaseUsageFromDeck(const Deck& deck)
{
    Runspec runspec( deck );
    const auto& phases = runspec.phases();

    return phaseUsage(phases);
}

}
