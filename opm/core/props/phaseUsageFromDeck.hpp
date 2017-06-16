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

#ifndef OPM_PHASEUSAGEFROMDECK_HEADER_INCLUDED
#define OPM_PHASEUSAGEFROMDECK_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Runspec.hpp>


namespace Opm
{

    /// Looks at presence of WATER, OIL and GAS keywords in state object
    /// to determine active phases.
    inline PhaseUsage phaseUsageFromDeck(const Opm::EclipseState& eclipseState)
    {
        PhaseUsage pu;
        std::fill(pu.phase_used, pu.phase_used + BlackoilPhases::MaxNumPhases, 0);

        const auto& phase = eclipseState.runspec().phases();
        // Discover phase usage.
        if (phase.active(Phase::WATER)) {
            pu.phase_used[BlackoilPhases::Aqua] = 1;
        }
        if (phase.active(Phase::OIL)) {
            pu.phase_used[BlackoilPhases::Liquid] = 1;
        }
        if (phase.active(Phase::GAS)) {
            pu.phase_used[BlackoilPhases::Vapour] = 1;
        }
        pu.num_phases = 0;
        for (int i = 0; i < BlackoilPhases::MaxNumPhases; ++i) {
            pu.phase_pos[i] = pu.num_phases;
            pu.num_phases += pu.phase_used[i];
        }

        // Only 2 or 3 phase systems handled.
        if (pu.num_phases < 2 || pu.num_phases > 3) {
            OPM_THROW(std::runtime_error, "Cannot handle cases with " << pu.num_phases << " phases.");
        }

        // We need oil systems, since we do not support the keywords needed for
        // water-gas systems.
        if (!pu.phase_used[BlackoilPhases::Liquid]) {
            OPM_THROW(std::runtime_error, "Cannot handle cases with no OIL, i.e. water-gas systems.");
        }

        // Add solvent info
        pu.has_solvent = false;
        if (phase.active(Phase::SOLVENT)) {
            pu.has_solvent = true;
        }

        // Add polyme info
        pu.has_polymer = false;
        if (phase.active(Phase::POLYMER)) {
            pu.has_polymer = true;
        }

        return pu;
    }

    /// Looks at presence of WATER, OIL and GAS keywords in deck
    /// to determine active phases.
    inline PhaseUsage phaseUsageFromDeck(const Opm::Deck& deck)
    {
        PhaseUsage pu;
        std::fill(pu.phase_used, pu.phase_used + BlackoilPhases::MaxNumPhases, 0);

        Runspec runspec( deck );
        const auto& phase = runspec.phases();

        // Discover phase usage.
        if (phase.active( Phase::WATER )) {
            pu.phase_used[BlackoilPhases::Aqua] = 1;
        }
        if (phase.active( Phase::OIL )) {
            pu.phase_used[BlackoilPhases::Liquid] = 1;
        }
        if (phase.active( Phase::GAS )) {
            pu.phase_used[BlackoilPhases::Vapour] = 1;
        }
        pu.num_phases = 0;
        for (int i = 0; i < BlackoilPhases::MaxNumPhases; ++i) {
            if (pu.phase_used[i]) {
                pu.phase_pos[i] = pu.num_phases;
                pu.num_phases += 1;
            }
            else {
                //Set to ridiculous value on purpose: should never be used
                pu.phase_pos[i] = 2000000000;
            }
        }

        // Only 2 or 3 phase systems handled.
        if (pu.num_phases < 2 || pu.num_phases > 3) {
            OPM_THROW(std::runtime_error, "Cannot handle cases with " << pu.num_phases << " phases.");
        }

        // We need oil systems, since we do not support the keywords needed for
        // water-gas systems.
        if (!pu.phase_used[BlackoilPhases::Liquid]) {
            OPM_THROW(std::runtime_error, "Cannot handle cases with no OIL, i.e. water-gas systems.");
        }

        // Add solvent info
        pu.has_solvent = false;
        if (phase.active(Phase::SOLVENT)) {
            pu.has_solvent = true;
        }

        // Add polyme info
        pu.has_polymer = false;
        if (phase.active(Phase::POLYMER)) {
            pu.has_polymer = true;
        }

        return pu;
    }

}

#endif // OPM_PHASEUSAGEFROMDECK_HEADER_INCLUDED
