/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILDEFS_HEADER_INCLUDED
#define OPM_BLACKOILDEFS_HEADER_INCLUDED

namespace Opm
{

    class BlackoilDefs
    {
    public:
        enum { numComponents = 3 };
        enum { numPhases = 3 };

        enum ComponentIndex { Water = 0, Gas = 1, Oil = 2 };
        enum PhaseIndex { Aqua = 0, Vapour = 1, Liquid = 2 };
    };

} // namespace Opm

#endif // OPM_BLACKOILDEFS_HEADER_INCLUDED
