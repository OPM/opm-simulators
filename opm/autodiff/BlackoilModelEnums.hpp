/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILMODELENUMS_HEADER_INCLUDED
#define OPM_BLACKOILMODELENUMS_HEADER_INCLUDED

#include <opm/autodiff/BlackoilPropsAdInterface.hpp>

namespace Opm
{

    enum Phases {
        Water        = BlackoilPropsAdInterface::Water,
        Oil          = BlackoilPropsAdInterface::Oil  ,
        Gas          = BlackoilPropsAdInterface::Gas  ,
        MaxNumPhases = BlackoilPropsAdInterface::MaxNumPhases
    };

    enum PrimalVariables {
        Sg = 0,
        RS = 1,
        RV = 2
    };


    enum CanonicalVariablePositions {
        Pressure = 0,
        Sw = 1,
        Xvar = 2,
        Qs = 3,
        Bhp = 4,
        Next // For extension.
    };

} // namespace Opm

#endif // OPM_BLACKOILMODELENUMS_HEADER_INCLUDED
