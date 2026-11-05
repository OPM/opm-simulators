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
#include <limits>

namespace Opm {

/// Marks a perforation created at runtime against a degree of freedom that has
/// no schedule connection -- a fracture cell represented outside the grid.
/// Consumers that map a perforation back to its ECL connection must skip these.
inline constexpr std::size_t AUX_PERFORATION_ECL_INDEX =
    std::numeric_limits<std::size_t>::max();

/// Static data associated with a well perforation.
template<class Scalar>
struct PerforationData
{
  static PerforationData serializationTestObject();

  bool operator==(const PerforationData&) const;

    int cell_index{};
    Scalar connection_transmissibility_factor{};
    Scalar connection_d_factor{};
    int satnum_id{};
    /// \brief The original index of the perforation in ECL Schedule
    std::size_t ecl_index{};
    /// \brief Grid identity for this perforation: 0 = global grid,
    ///        > 0 = numeric LGR level (refined grid in CARFIN
    ///        declaration order).
    int grid_id{};
    /// \brief Linearised Cartesian cell index relative to the cell's
    ///        own grid origin.  For grid_id == 0 this is the global-grid
    ///        Cartesian index; for grid_id > 0 this is the LGR-local
    ///        Cartesian index.
    ///
    /// Together with grid_id, the pair (grid_id, global_index) uniquely
    /// identifies any cell in the deck -- including refined siblings
    /// under a common coarse parent, which would otherwise alias if
    /// identified via globalCellIdxMap[active_idx] alone (which returns
    /// the level-0 ancestor's Cartesian, shared by all siblings).
    std::size_t global_index{};
    /// \brief Segment number (1-based) of a multi-segment well's perforation that
    ///        has no schedule connection; 0 otherwise.
    int segment{};

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
      serializer(cell_index);
      serializer(connection_transmissibility_factor);
      serializer(connection_d_factor);
      serializer(satnum_id);
      serializer(ecl_index);
      serializer(grid_id);
      serializer(global_index);
      serializer(segment);
    }
};

template<class Scalar>
struct PerforationRates
{
    Scalar dis_gas = 0.0;
    Scalar dis_gas_in_water = 0.0;
    Scalar vap_oil = 0.0;
    Scalar vap_wat = 0.0;
};

} // namespace Opm

#endif // OPM_PERFORATIONDATA_HEADER_INCLUDED
