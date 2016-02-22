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

#ifndef OPM_POLYMERSTATE_HEADER_INCLUDED
#define OPM_POLYMERSTATE_HEADER_INCLUDED


#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/grid.h>
#include <vector>

namespace Opm
{

    /// Simulator state for a two-phase simulator with polymer.
    class PolymerState : public SimulatorState
    {
    public:

        void init(int number_of_cells, int number_of_faces, int num_phases)
        {
            SimulatorState::init( number_of_cells , number_of_faces , num_phases );
            registerCellData("CONCENTRATION" , 1 , 0 );
            registerCellData("CMAX" , 1 , 0 );
        }

        const std::vector<double>& concentration() const {
            return getCellData("CONCENTRATION");
        }

        const std::vector<double>& maxconcentration() const {
            return getCellData("CMAX");
        }


        std::vector<double>& concentration()  {
            return getCellData("CONCENTRATION");
        }

        std::vector<double>& maxconcentration()  {
            return getCellData("CMAX");
        }

    };

} // namespace Opm




#endif // OPM_POLYMERSTATE_HEADER_INCLUDED
