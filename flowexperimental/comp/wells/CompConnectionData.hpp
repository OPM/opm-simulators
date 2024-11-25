/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_COMP_CONNECTION_DATA_HPP
#define OPM_COMP_CONNECTION_DATA_HPP

#include <cstddef>

namespace Opm {

template<typename Scalar>
struct CompConnectionData {
    std::size_t cell_index{};
    Scalar connection_transmissibility_factor{};
    Scalar connection_d_factor{};
    int satnum_id{};

    /// \brief The original index of the connection in ECL Schedule
    std::size_t ecl_index{};
};

} // namespace Opm

#endif // OPM_COMP_CONNECTION_DATA_HPP