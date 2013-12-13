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

#ifndef OPM_BLACKOILSTATE_HEADER_INCLUDED
#define OPM_BLACKOILSTATE_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/simulator/SimulatorState.hpp>
#include <vector>

namespace Opm
{

    /// Simulator state for a blackoil simulator.
    class BlackoilState : public SimulatorState
    {
    public:
        virtual void init(const UnstructuredGrid& g, int num_phases);

        /// Set the first saturation to either its min or max value in
        /// the indicated cells. The second saturation value s2 is set
        /// to (1.0 - s1) for each cell. Any further saturation values
        /// are unchanged.
        void setFirstSat(const std::vector<int>& cells,
                         const Opm::BlackoilPropertiesInterface& props,
                         ExtremalSat es);

        virtual bool equals(const SimulatorState& other,
                            double epsilon = 1e-8) const;

        std::vector<double>& surfacevol  () { return surfvol_; }
        std::vector<double>& gasoilratio () { return gor_   ; }
        std::vector<double>& rv () {return rv_ ; }

        const std::vector<double>& surfacevol  () const { return surfvol_; }
        const std::vector<double>& gasoilratio () const { return gor_   ; }
        const std::vector<double>& rv () const {return rv_ ; }

    private:
        std::vector<double> surfvol_;
        std::vector<double> gor_   ;
        std::vector<double> rv_ ;
    };

} // namespace Opm


#endif // OPM_BLACKOILSTATE_HEADER_INCLUDED
