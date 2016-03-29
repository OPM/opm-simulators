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


#include <opm/common/data/SimulationDataContainer.hpp>

#include <opm/polymer/PolymerState.hpp>

namespace Opm
{
    const std::string PolymerState::CONCENTRATION = "CONCENTRATION";
    const std::string PolymerState::CMAX = "CMAX";

    PolymerState::PolymerState(int number_of_cells, int number_of_faces, int num_phases) :
        SimulationDataContainer( number_of_cells , number_of_faces , num_phases )
    {
        registerCellData(CONCENTRATION , 1 , 0 );
        registerCellData(CMAX , 1 , 0 );
    }


} // namespace Opm




