
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
    /// Determine the active phases
    inline PhaseUsage phaseUsage(const Phases& phases) {
        PhaseUsage pu;
        std::fill(pu.phase_used, pu.phase_used + BlackoilPhases::MaxNumPhases + BlackoilPhases::NumCryptoPhases, 0);

        // Discover phase usage.
        pu.phase_used[BlackoilPhases::Aqua] = phases.active(Phase::WATER);
        pu.phase_used[BlackoilPhases::Liquid] = phases.active(Phase::OIL);
        pu.phase_used[BlackoilPhases::Vapour] = phases.active(Phase::GAS);

        pu.num_phases = 0;
        int numActivePhases = 0;
        for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
            if (!pu.phase_used[phaseIdx]) {
                pu.phase_pos[phaseIdx] = -1;
            }
            else {
                pu.phase_pos[phaseIdx] = numActivePhases;
                ++ numActivePhases;
                pu.num_phases = numActivePhases;
            }
        }

        // We need oil systems, since we do not support the keywords needed for
        // water-gas systems.
        if (!pu.phase_used[BlackoilPhases::Liquid] && !(pu.num_phases == 1)) {
            OPM_THROW(std::runtime_error, "Cannot handle cases with no OIL, i.e. water-gas systems.");
        }

        // Add solvent info
        pu.has_solvent = phases.active(Phase::SOLVENT);
        if (pu.has_solvent) {
            // this is quite a hack: even though solvent is not considered as in
            // MaxNumPhases and pu.num_phases because this would break a lot of
            // assumptions in old code, it is nevertheless an index to be translated
            // to. solvent and solvent are even larger hacks because not even this can be
            // done for them.
            pu.phase_pos[BlackoilPhases::Solvent] = numActivePhases;
            ++ numActivePhases;
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
            pu.phase_pos[BlackoilPhases::Polymer] = numActivePhases;
            ++ numActivePhases;
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
            pu.phase_pos[BlackoilPhases::Energy] = numActivePhases;
            ++ numActivePhases;
        }
        else
            pu.phase_pos[BlackoilPhases::Energy] = -1;

        // Add polymer molecular weight related
        pu.has_polymermw = phases.active(Phase::POLYMW);
        if (pu.has_polymermw) {
            if (!pu.has_polymer) {
                OPM_THROW(std::runtime_error, "pu.has_polymermw is true while pu.has_polymer is false");
            }
            pu.phase_pos[BlackoilPhases::PolymerMW] = numActivePhases;
            ++ numActivePhases;
        }
        else
            pu.phase_pos[BlackoilPhases::PolymerMW] = -1;

        // Add foam info
        pu.has_foam = phases.active(Phase::FOAM);
        if (pu.has_foam) {
            pu.phase_pos[BlackoilPhases::Foam] = numActivePhases;
            ++ numActivePhases;
        }
        else
            pu.phase_pos[BlackoilPhases::Foam] = -1;

        // Add salt info
        pu.has_salt = phases.active(Phase::SALTWATER);
        if (pu.has_salt) {
            pu.phase_pos[BlackoilPhases::Salt] = numActivePhases;
            ++ numActivePhases;
        }
        else
            pu.phase_pos[BlackoilPhases::Salt] = -1;

        return pu;
    }

    /// Looks at presence of WATER, OIL and GAS keywords in state object
    /// to determine active phases.
    inline PhaseUsage phaseUsageFromDeck(const Opm::EclipseState& eclipseState)
    {
        const auto& phases = eclipseState.runspec().phases();

        return phaseUsage(phases);
    }

    /// Looks at presence of WATER, OIL and GAS keywords in deck
    /// to determine active phases.
    inline PhaseUsage phaseUsageFromDeck(const Opm::Deck& deck)
    {
        Runspec runspec( deck );
        const auto& phases = runspec.phases();

        return phaseUsage(phases);
    }

}

#endif // OPM_PHASEUSAGEFROMDECK_HEADER_INCLUDED
