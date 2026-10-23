/*
  Copyright 2026 Equinor ASA

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

#ifndef OPM_GROUP_CONTROL_RETENTION_HPP
#define OPM_GROUP_CONTROL_RETENTION_HPP

#include <opm/input/eclipse/Schedule/ResCoup/GrupSlav.hpp>

#include <optional>

namespace Opm::GroupStateHelpers {

/// Whether a production group keeps its control mode when the groups without
/// active constraints are set to NONE at the end of a time step.
///
/// @param field_or_none The production control mode of the group is NONE or FLD, so
///   there is nothing to reset.
/// @param targeted A well under group control has the group as its target.
/// @param gas_lift_group The group is in the gas lift optimization, which needs a
///   control mode other than NONE and FLD.
/// @param rc_master_hierarchy The group is a reservoir coupling master group or one of
///   its ancestors, which distribute targets to the slaves.
/// @param master_target_filter Set when the group is a reservoir coupling slave group with a
///   target from the master: the GRUPSLAV flag for that target.  The master sends the target
///   in the first time step of a sync step only, and the control mode it gives the group
///   (not for SLAV, where the slave's own limit is the one used) has to survive the time steps
///   after it.  Not set when the group has no target from a master.
inline bool keepsProductionControl(const bool field_or_none,
                                   const bool targeted,
                                   const bool gas_lift_group,
                                   const bool rc_master_hierarchy,
                                   const std::optional<ReservoirCoupling::GrupSlav::FilterFlag>
                                       master_target_filter)
{
    using FilterFlag = ReservoirCoupling::GrupSlav::FilterFlag;
    const bool master_sets_control_mode =
        master_target_filter.has_value() && *master_target_filter != FilterFlag::SLAV;
    return targeted
        || field_or_none
        || gas_lift_group
        || rc_master_hierarchy
        || master_sets_control_mode;
}

} // namespace Opm::GroupStateHelpers

#endif // OPM_GROUP_CONTROL_RETENTION_HPP
