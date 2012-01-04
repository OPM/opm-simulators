/*
  Copyright 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/fluid/blackoil/SinglePvtInterface.hpp>


namespace Opm
{

    SinglePvtInterface::SinglePvtInterface()
        : num_phases_(MaxNumPhases)
    {
        for (int i = 0; i < MaxNumPhases; ++i) {
            phase_pos_[i] = i;
        }
    }

    SinglePvtInterface::~SinglePvtInterface()
    {
    }

    void SinglePvtInterface::setPhaseConfiguration(const int num_phases, const int* phase_pos)
    {
        num_phases_ = num_phases;
        for (int i = 0; i < MaxNumPhases; ++i) {
            phase_pos_[i] = phase_pos[i];
        }
    }


} // namespace Opm
