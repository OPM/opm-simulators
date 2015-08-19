/*
  Copyright 2015 IRIS AS, Applied Mathematics.

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

#ifndef OPM_BLACKOILSOLVENTSTATE_HEADER_INCLUDED
#define OPM_BLACKOILSOLVENTSTATE_HEADER_INCLUDED

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/grid.h>
#include <vector>


namespace Opm
{

    /// Simulator state for blackoil simulator with solvent.
    /// We use the Blackoil state parameters.
    class BlackoilSolventState : public BlackoilState
    {
    public:
        void init(const UnstructuredGrid& g, int num_phases)
        {
            this->init(g.number_of_cells, g.number_of_faces, num_phases);
        }

        void init(int number_of_cells, int number_of_faces, int num_phases)
        {
            BlackoilState::init(number_of_cells, number_of_faces, num_phases);
            solventId_  = SimulatorState::registerCellData( "SSOL", 1 );
        }

        std::vector<double>& solvent_saturation()    { return cellData()[ solventId_ ]; }
        const std::vector<double>& solvent_saturation() const    { return cellData()[ solventId_ ]; }

    private:
        int solventId_;
    };

} // namespace Opm




#endif // OPM_BLACKOILSOLVENTSTATE_HEADER_INCLUDED
