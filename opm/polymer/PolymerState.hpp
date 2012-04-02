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


#include <opm/core/TwophaseState.hpp>
#include <opm/core/grid.h>
#include <vector>

namespace Opm
{

    /// Simulator state for a two-phase simulator with polymer.
    class PolymerState : public TwophaseState
    {
    public:
        void init(const UnstructuredGrid& g)
        {
            TwophaseState::init(g);
            concentration_.resize(g.number_of_cells, 0.0);
            cmax_.resize(g.number_of_cells, 0.0);
        }

        std::vector<double>& concentration() { return concentration_; }
        std::vector<double>& maxconcentration() { return cmax_; }

        const std::vector<double>& concentration() const { return concentration_; }
        const std::vector<double>& maxconcentration() const { return cmax_; }

    private:
        std::vector<double> concentration_;
        std::vector<double> cmax_;
    };

} // namespace Opm




#endif // OPM_POLYMERSTATE_HEADER_INCLUDED
