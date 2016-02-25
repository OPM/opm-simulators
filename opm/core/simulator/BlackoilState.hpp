/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2015 IRIS AS

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

#include <opm/common/data/SimulationDataContainer.hpp>

#include <opm/core/grid.h>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <vector>

namespace Opm
{

    /// Simulator state for a blackoil simulator.
    class BlackoilState : public SimulationDataContainer
    {
    public:
        static const std::string GASOILRATIO;
        static const std::string RV;
        static const std::string SURFACEVOL;

        BlackoilState(size_t num_cells , size_t num_faces, size_t num_phases);


        std::vector<double>& surfacevol  () { return getCellData("SURFACEVOL");  }
        std::vector<double>& gasoilratio () { return getCellData(GASOILRATIO); }
        std::vector<double>& rv ()          { return getCellData(RV);          }

        const std::vector<double>& surfacevol  () const { return getCellData("SURFACEVOL");  }
        const std::vector<double>& gasoilratio () const { return getCellData("GASOILRATIO"); }
        const std::vector<double>& rv ()          const { return getCellData(RV);          }

    };
} // namespace Opm


#endif // OPM_BLACKOILSTATE_HEADER_INCLUDED
