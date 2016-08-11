/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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
#ifndef OPM_DEFAULTBLACKOILSOLUTIONSTATE_HEADER_INCLUDED
#define OPM_DEFAULTBLACKOILSOLUTIONSTATE_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>

namespace Opm {
    /// Struct for containing iteration variables.
    struct DefaultBlackoilSolutionState
    {
        typedef AutoDiffBlock<double> ADB;
        explicit DefaultBlackoilSolutionState(const int np)
            : pressure  (    ADB::null())
            , temperature(   ADB::null())
            , saturation(np, ADB::null())
            , rs        (    ADB::null())
            , rv        (    ADB::null())
            , qs        (    ADB::null())
            , bhp       (    ADB::null())
            , wellVariables (    ADB::null())
            , canonical_phase_pressures(3, ADB::null())
        {
        }
        ADB              pressure;
        ADB              temperature;
        std::vector<ADB> saturation;
        ADB              rs;
        ADB              rv;
        ADB              qs;
        ADB              bhp;
        ADB wellVariables;
        // Below are quantities stored in the state for optimization purposes.
        std::vector<ADB> canonical_phase_pressures; // Always has 3 elements, even if only 2 phases active.
    };
} // namespace Opm

#endif // OPM_DEFAULTBLACKOILSOLUTIONSTATE_HEADER_INCLUDED
