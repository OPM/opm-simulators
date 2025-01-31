/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_RUNTIME_PERFORATION_HPP_INCLUDED
#define OPM_RUNTIME_PERFORATION_HPP_INCLUDED

namespace Opm {

/// Simple model of a well connection created at runtime, possibly as a
/// result of a geo-mechanical fracturing process.
struct RuntimePerforation
{
    /// Active cell index, on current rank, that is dynamically perforated
    /// by a well.
    int cell{};

    /// Connection's transmissibility factor.
    double ctf{};

    /// Depth at which the new connection is created.
    double depth{};
};

} // namespace Opm

#endif // OPM_RUNTIME_PERFORATION_HPP_INCLUDED
