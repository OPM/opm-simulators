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

#ifndef OPM_ROCKCOMPRESSIBILITY_HEADER_INCLUDED
#define OPM_ROCKCOMPRESSIBILITY_HEADER_INCLUDED

#include <vector>

namespace Opm
{

    class EclipseGridParser;

    class RockCompressibility
    {
    public:
        /// Construct from input deck.
        RockCompressibility(const EclipseGridParser& deck);

        /// Porosity multiplier.
        double poroMult(double pressure);

        /// Rock compressibility = (d poro / d p)*(1 / poro).
        double rockComp(double pressure);

    private:
        std::vector<double> p_;
        std::vector<double> poromult_;
        double pref_;
        double rock_comp_;
    };

} // namespace Opm


#endif // OPM_ROCKCOMPRESSIBILITY_HEADER_INCLUDED
