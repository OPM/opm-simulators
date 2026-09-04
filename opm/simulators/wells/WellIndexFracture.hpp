/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_WELL_INDEX_FRACTURE_HPP_INCLUDED
#define OPM_WELL_INDEX_FRACTURE_HPP_INCLUDED

namespace Opm {

/// Additional connection transmissibility contributed by a dynamically
/// created fracture crossing an existing well connection.  Populated by a
/// fracture model through WellInterfaceGeneric::addFracturePerforations();
/// zero-initialized entries contribute nothing.
struct WellIndexFracture
{
    /// Fracture connection transmissibility factor.
    double ctf{0.0};

    /// Connection pressure at which \c ctf was computed.
    double pressure{0.0};

    /// Reference values from the previous fracture solve; reserved for
    /// pressure interpolation of the connection factor between solves
    /// (not yet enabled -- needs validation).
    double ref_ctf{0.0};
    double ref_pressure{0.0};

    /// Fracture connection factor at the given connection pressure.
    ///
    /// The pressure argument is reserved for interpolating between
    /// (pressure, ctf) and (ref_pressure, ref_ctf); until that is
    /// validated the stored factor is returned unchanged.
    double wellIndex(const double /* connection_pressure */) const
    {
        return this->ctf;
    }
};

} // namespace Opm

#endif // OPM_WELL_INDEX_FRACTURE_HPP_INCLUDED
