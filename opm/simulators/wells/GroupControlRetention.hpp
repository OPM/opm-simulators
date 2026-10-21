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

#include <opm/input/eclipse/Schedule/Group/Group.hpp>

namespace Opm::GroupStateHelpers {

/// Whether a production group keeps its control mode when the groups without
/// active constraints are set to NONE at the end of a time step.
///
/// @param cmode Production control mode of the group.
/// @param targeted A well under group control has the group as its target.
/// @param gas_lift_group The group is in the gas lift optimization, which needs a
///   control mode other than NONE and FLD.
/// @param rc_master_hierarchy The group is a reservoir coupling master group or one of
///   its ancestors, which distribute targets to the slaves.
/// @param rc_slave_master_target The group is a reservoir coupling slave group whose
///   control mode comes from the master, which sends its target only in the first time
///   step of a sync step, so the mode has to survive the time steps after it.
inline bool keepsProductionControl(const Group::ProductionCMode cmode,
                                   const bool targeted,
                                   const bool gas_lift_group,
                                   const bool rc_master_hierarchy,
                                   const bool rc_slave_master_target)
{
    return targeted
        || cmode == Group::ProductionCMode::NONE
        || cmode == Group::ProductionCMode::FLD
        || gas_lift_group
        || rc_master_hierarchy
        || rc_slave_master_target;
}

} // namespace Opm::GroupStateHelpers

#endif // OPM_GROUP_CONTROL_RETENTION_HPP
