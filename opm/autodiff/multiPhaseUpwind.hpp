/*
  Copyright 2015, 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil AS.

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

#ifndef OPM_MULTIPHASEUPWIND_HEADER_INCLUDED
#define OPM_MULTIPHASEUPWIND_HEADER_INCLUDED

#include <array>

namespace Opm
{
    /// Compute upwind directions for three-phase flow across a connection.
    ///
    /// @param[in]  head_diff           head differences by phase
    /// @param[in]  mob1                phase mobilities for first cell
    /// @param[in]  mob2                phase mobilities for second cell
    /// @param[in]  transmissibility    tranmissibility of connection
    /// @param[in]  flux                total volume flux across connection
    /// @return array containing, for each phase, 1.0 if flow in the
    ///         direction of the connection, -1.0 if flow in the opposite
    ///         direction.
    std::array<double, 3> connectionMultiPhaseUpwind(const std::array<double, 3>& head_diff,
                                                     const std::array<double, 3>& mob1,
                                                     const std::array<double, 3>& mob2,
                                                     const double transmissibility,
                                                     const double flux);
} // namespace Opm

#endif // OPM_MULTIPHASEUPWIND_HEADER_INCLUDED
