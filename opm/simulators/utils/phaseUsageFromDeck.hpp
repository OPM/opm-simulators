
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

#include <opm/simulators/utils/BlackoilPhases.hpp>

namespace Opm
{
    class Deck;
    class EclipseState;
    class Phases;

    /// Determine the active phases
    PhaseUsage phaseUsage(const Phases& phases);

    /// Looks at presence of WATER, OIL and GAS keywords in state object
    /// to determine active phases.
    PhaseUsage phaseUsageFromDeck(const EclipseState& eclipseState);

    /// Looks at presence of WATER, OIL and GAS keywords in deck
    /// to determine active phases.
    PhaseUsage phaseUsageFromDeck(const Deck& deck);
}

#endif // OPM_PHASEUSAGEFROMDECK_HEADER_INCLUDED
