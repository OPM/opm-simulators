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

BOOST_AUTO_TEST_CASE(untargeted_group_is_reset)
{
    BOOST_CHECK(!keepsProductionControl(/*field_or_none=*/false,
                                        /*targeted=*/false,
                                        /*gas_lift_group=*/false,
                                        /*rc_master_hierarchy=*/false,
                                        /*rc_slave_master_target=*/false));
}

BOOST_AUTO_TEST_CASE(targeted_group_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/true,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/false));
}

BOOST_AUTO_TEST_CASE(none_and_fld_have_nothing_to_reset)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/true,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/false));
}

BOOST_AUTO_TEST_CASE(gas_lift_group_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/true,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/false));
}

BOOST_AUTO_TEST_CASE(master_hierarchy_group_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/true,
                                       /*rc_slave_master_target=*/false));
}

// A slave group whose control mode is the master's has no wells under group control
// when they are all on THP control, but the mode must stay for the later time steps
// of the sync step, or the slave stops at "cmode mismatch between slave and master".
BOOST_AUTO_TEST_CASE(slave_group_with_master_target_keeps_its_control)
{
    BOOST_CHECK(keepsProductionControl(/*field_or_none=*/false,
                                       /*targeted=*/false,
                                       /*gas_lift_group=*/false,
                                       /*rc_master_hierarchy=*/false,
                                       /*rc_slave_master_target=*/true));
}
