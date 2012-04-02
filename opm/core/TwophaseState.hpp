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

#ifndef OPM_TWOPHASESTATE_HEADER_INCLUDED
#define OPM_TWOPHASESTATE_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <vector>

namespace Opm
{

    /// Simulator state for a two-phase simulator.
    class TwophaseState
    {
    public:

        void init(const UnstructuredGrid& g)
        {
            press_.resize(g.number_of_cells, 0.0);
            fpress_.resize(g.number_of_faces, 0.0);
            flux_.resize(g.number_of_faces, 0.0);
            sat_.resize(2 * g.number_of_cells, 0.0);
            for (int cell = 0; cell < g.number_of_cells; ++cell) {
                sat_[2*cell + 1] = 1.0; // Defaulting oil saturations to 1.0.
            }
        }

        enum ExtremalSat { MinSat, MaxSat };

        void setWaterSat(const std::vector<int>& cells,
                         const Opm::IncompPropertiesInterface& props,
                         ExtremalSat es)
        {
            const int n = cells.size();
            std::vector<double> smin(2*n);
            std::vector<double> smax(2*n);
            props.satRange(n, &cells[0], &smin[0], &smax[0]);
            const double* svals = (es == MinSat) ? &smin[0] : &smax[0];
            for (int ci = 0; ci < n; ++ci) {
                const int cell = cells[ci];
                sat_[2*cell] = svals[2*ci];
                sat_[2*cell + 1] = 1.0 - sat_[2*cell];
            }
        }

        int numPhases() const
        {
            return 2;
        }

        std::vector<double>& pressure    () { return press_ ; }
        std::vector<double>& facepressure() { return fpress_; }
        std::vector<double>& faceflux    () { return flux_  ; }
        std::vector<double>& saturation  () { return sat_   ; }

        const std::vector<double>& pressure    () const { return press_ ; }
        const std::vector<double>& facepressure() const { return fpress_; }
        const std::vector<double>& faceflux    () const { return flux_  ; }
        const std::vector<double>& saturation  () const { return sat_   ; }

    private:
        std::vector<double> press_ ;
        std::vector<double> fpress_;
        std::vector<double> flux_  ;
        std::vector<double> sat_   ;
    };

} // namespace Opm


#endif // OPM_TWOPHASESTATE_HEADER_INCLUDED
