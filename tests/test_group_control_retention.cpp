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

#include <config.h>

#define BOOST_TEST_MODULE GroupControlRetentionTests

#include <opm/simulators/wells/GroupControlRetention.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

using Opm::GroupStateHelpers::keepsProductionControl;

// A production group with a control mode other than NONE and FLD (say ORAT), that no well
// under group control has as its target, that is not in the gas lift optimization, not in a
// reservoir coupling master hierarchy and has no target from a master: it is reset to NONE.
BOOST_AUTO_TEST_CASE(untargeted_group_is_reset)
{
    BOOST_CHECK(!keepsProductionControl(/*field_or_none=*/false,
                                        /*targeted=*/false,
                                        /*gas_lift_group=*/false,
                                        /*rc_master_hierarchy=*/false,
                                        /*rc_slave_master_target=*/false));
}

// The same group, but a well under group control has it as its target: it keeps its control
// mode.
BOOST_AUTO_TEST_CASE(targeted_group_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/true,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/false));
}

// A group whose control mode is already NONE or FLD (field_or_none) has nothing to reset,
// whatever the other conditions are.
BOOST_AUTO_TEST_CASE(none_and_fld_have_nothing_to_reset)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/true,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/false));
}

// An untargeted group in the gas lift optimization keeps its control mode, as the
// optimization needs a mode other than NONE and FLD.
BOOST_AUTO_TEST_CASE(gas_lift_group_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/true,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/false));
}

// An untargeted reservoir coupling master group, or an ancestor of one, keeps its control
// mode: it distributes targets to the slaves.
BOOST_AUTO_TEST_CASE(master_hierarchy_group_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/true,
                                       /*rc_slave_master_target=*/false));
}

// An untargeted reservoir coupling slave group with a target from the master keeps the
// master's control mode.  Its wells can all be on THP control, but the master sends the
// target in the first time step of a sync step only, and the mode must stay for the later
// time steps, or the slave stops at "cmode mismatch between slave and master".
BOOST_AUTO_TEST_CASE(slave_group_with_master_target_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/true));
}
