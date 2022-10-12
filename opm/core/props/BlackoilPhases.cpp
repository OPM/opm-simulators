/*
  Copyright 2021 Equinor ASA

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

#include <algorithm>
#include <vector>

#include <opm/core/props/BlackoilPhases.hpp>

Opm::PhaseUsage::PhaseUsage(std::vector<BlackoilPhases::PhaseIndex> phases)
{
    std::sort(phases.begin(), phases.end());
    this->phase_used.fill(0);
    this->phase_pos.fill(-1);

    std::size_t current_pos = 0;
    for (const auto& phase : phases) {
        this->phase_used[phase] = 1;
        this->phase_pos[phase] = current_pos;

        current_pos++;
    }

    this->num_phases = 0;
    if (this->phase_used[BlackoilPhases::Aqua])
        this->num_phases++;
    if (this->phase_used[BlackoilPhases::Liquid])
        this->num_phases++;
    if (this->phase_used[BlackoilPhases::Vapour])
        this->num_phases++;

    this->has_solvent   = this->phase_used[BlackoilPhases::Solvent];
    this->has_polymer   = this->phase_used[BlackoilPhases::Polymer];
    this->has_energy    = this->phase_used[BlackoilPhases::Energy];
    this->has_polymermw = this->phase_used[BlackoilPhases::PolymerMW];
    this->has_foam      = this->phase_used[BlackoilPhases::Foam];
    this->has_brine     = this->phase_used[BlackoilPhases::Brine];
    this->has_zFraction = this->phase_used[BlackoilPhases::ZFraction];
}
