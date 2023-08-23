/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_PERFORATIONDATA_HEADER_INCLUDED
#define OPM_PERFORATIONDATA_HEADER_INCLUDED

#include <cstddef>

namespace Opm
{

/// Static data associated with a well perforation.
struct PerforationData
{
    int cell_index;
    double connection_transmissibility_factor;
    double connection_d_factor;
    int satnum_id;
    /// \brief The original index of the perforation in ECL Schedule
    std::size_t ecl_index;
};

struct PerforationRates
{
    double dis_gas = 0.0;
    double dis_gas_in_water = 0.0;
    double vap_oil = 0.0;
    double vap_wat = 0.0;
};

} // namespace Opm

#endif // OPM_PERFORATIONDATA_HEADER_INCLUDED
