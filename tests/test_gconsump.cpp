// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!\
 * \brief Test suite for GCONSUMP efficiency factor handling
 *
 * Tests that hierarchical efficiency factors are correctly applied to GCONSUMP calculations
 * following the REIN pattern where groups are not affected by their own efficiency factors,
 * but parent groups apply child efficiency factors when accumulating child contributions.
 */

#define BOOST_TEST_MODULE GConsumpEfficiencyFactor
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include "ToleranceAndUnitFixture.hpp"

#include <boost/test/unit_test.hpp>

using namespace Opm;

// Using shared ToleranceAndUnitFixture from ToleranceAndUnitFixture.hpp
//! \brief Test fixture for group constraint tests to avoid code duplication
//!
//! Tests can customize behavior via configuration options before calling setup().
struct GConsumpEfficiencyFactorTestFixture : ToleranceAndUnitFixture
{
    void setup(const std::string& deck_file = "GCONSUMP.DATA")
    {
        // Load the test deck with GCONSUMP and GEFAC specification
        Parser parser;
        auto deck = parser.parseFile(deck_file);
        EclipseState eclipseState{deck};
        schedule_ = Schedule{deck, eclipseState};
        report_step_ = 0;
        num_phases_ = 3;
        group_state_ = GroupState<double>{num_phases_};
        summary_state_ = SummaryState{schedule_.getStartTime()};
        group_state_.update_gconsump(schedule_, report_step_, summary_state_);
    }

    Schedule& schedule() { return schedule_; }
    int report_step() { return report_step_; }
    const ScheduleState& scheduleState() { return schedule()[report_step()]; }
    const GConSump& gconsump() { return scheduleState().gconsump(); }
    GroupState<double>& group_state() { return group_state_; }
    SummaryState& summary_state() { return summary_state_; }

private:
    Schedule schedule_;
    int report_step_;
    std::size_t num_phases_;
    GroupState<double> group_state_;
    SummaryState summary_state_;
};

BOOST_FIXTURE_TEST_SUITE(GConsumpEfficiencyFactorTests, GConsumpEfficiencyFactorTestFixture)

//! Test GCONSUMP efficiency factor implementation following REIN pattern
//!
//! Verifies that groups are not affected by their own efficiency factors,
//! but parent groups apply child efficiency factors when accumulating contributions.
//!
//! Group hierarchy:
//! \code
//!   FIELD (GEFAC=1.0)
//!   └── PLAT-A (GEFAC=0.9)
//!       └── SUB-B (GEFAC=0.8, GCONSUMP=100/50)
//!           └── WELL-C
//! \endcode
BOOST_AUTO_TEST_CASE(TestGconsumpEfficiencyFactorBug)
{
    setup();

    const auto& group_sub_b = schedule().getGroup("SUB-B", report_step());
    const auto& group_plat_a = schedule().getGroup("PLAT-A", report_step());
    const auto& group_field = schedule().getGroup("FIELD", report_step());

    const auto& sub_b_efficiency = group_sub_b.getGroupEfficiencyFactor();
    const auto& plat_a_efficiency = group_plat_a.getGroupEfficiencyFactor();
    const auto& field_efficiency = group_field.getGroupEfficiencyFactor();

    // Validate efficiency factors are as expected
    checkClose(sub_b_efficiency, 0.8);
    checkClose(plat_a_efficiency, 0.9);
    checkClose(field_efficiency, 1.0);

    // Note: Cumulative efficiency calculation is for reference only
    const double cumulative_efficiency = sub_b_efficiency * plat_a_efficiency;
    checkClose(cumulative_efficiency, 0.72);

    // Get the efficiency-adjusted rates from update_gconsump()
    const auto& sub_b_rates = group_state().gconsump_rates("SUB-B");
    const auto& plat_a_rates = group_state().gconsump_rates("PLAT-A");
    const auto& field_rates = group_state().gconsump_rates("FIELD");

    // Raw GCONSUMP values for SUB-B from deck
    const auto& group_gc = gconsump().get("SUB-B", summary_state());
    const double raw_consumption = metric_rate(group_gc.consumption_rate);
    const double raw_import = metric_rate(group_gc.import_rate);

    // Eclipse-compliant expected values (correct efficiency pattern)
    const double expected_sub_b_consumption = raw_consumption;    // 100.0
    const double expected_sub_b_import = raw_import;              // 50.0
    const double expected_plat_a_consumption = raw_consumption * sub_b_efficiency;  // 80.0
    const double expected_plat_a_import = raw_import * sub_b_efficiency;            // 40.0
    const double expected_field_consumption = expected_plat_a_consumption * plat_a_efficiency;  // 72.0
    const double expected_field_import = expected_plat_a_import * plat_a_efficiency;            // 36.0

    // Convert actual SI rates to metric for direct comparison with expected metric values
    const double actual_sub_b_consumption_metric = metric_rate(sub_b_rates.first);
    const double actual_sub_b_import_metric = metric_rate(sub_b_rates.second);
    const double actual_plat_a_consumption_metric = metric_rate(plat_a_rates.first);
    const double actual_plat_a_import_metric = metric_rate(plat_a_rates.second);
    const double actual_field_consumption_metric = metric_rate(field_rates.first);
    const double actual_field_import_metric = metric_rate(field_rates.second);

    // Test the corrected GCONSUMP efficiency factor implementation
    checkRate(actual_sub_b_consumption_metric, expected_sub_b_consumption);
    checkRate(actual_sub_b_import_metric, expected_sub_b_import);
    checkRate(actual_plat_a_consumption_metric, expected_plat_a_consumption);
    checkRate(actual_plat_a_import_metric, expected_plat_a_import);
    checkRate(actual_field_consumption_metric, expected_field_consumption);
    checkRate(actual_field_import_metric, expected_field_import);
}

//! Test for groups without GCONSUMP - should return zero rates
//!
//! Group hierarchy:
//! \code
//!   FIELD (GEFAC=1.0)
//!   └── PLAT-A (GEFAC=0.9)
//!       └── SUB-B (GEFAC=0.8, GCONSUMP=100/50)
//!           └── WELL-C (no GCONSUMP - should return 0/0)
//! \endcode
BOOST_AUTO_TEST_CASE(TestGconsumpZeroRates)
{
    setup();

    // WELL-C has no GCONSUMP specified - should return zero rates
    const auto& well_c_rates = group_state().gconsump_rates("WELL-C");

    BOOST_CHECK_EQUAL(well_c_rates.first, 0.0);
    BOOST_CHECK_EQUAL(well_c_rates.second, 0.0);
}

//! Test complex GCONSUMP hierarchy with multiple levels and sibling groups
//!
//! Tests realistic scenarios with GCONSUMP at multiple hierarchy levels,
//! sibling groups, and mixed scenarios (groups with/without GCONSUMP).
//!
//! Group hierarchy:
//! \code
//!   FIELD (GEFAC=1.0, GCONSUMP=20/10)
//!   └── PLAT-A (GEFAC=0.9, GCONSUMP=30/15)
//!       ├── SUB-B (GEFAC=0.8, GCONSUMP=100/50)
//!       │   └── WELL-C
//!       └── SUB-D (GEFAC=0.7, no GCONSUMP)
//!           └── WELL-E
//! \endcode
BOOST_AUTO_TEST_CASE(TestGconsumpComplexHierarchy)
{
    setup("GCONSUMP_COMPLEX.DATA");

    // Get all group rates
    const auto& sub_b_rates = group_state().gconsump_rates("SUB-B");
    const auto& sub_d_rates = group_state().gconsump_rates("SUB-D");
    const auto& plat_a_rates = group_state().gconsump_rates("PLAT-A");
    const auto& field_rates = group_state().gconsump_rates("FIELD");

    // Convert to metric for comparison
    const auto sub_b_consumption_metric = metric_rate(sub_b_rates.first);
    const auto sub_b_import_metric = metric_rate(sub_b_rates.second);
    const auto sub_d_consumption_metric = metric_rate(sub_d_rates.first);
    const auto sub_d_import_metric = metric_rate(sub_d_rates.second);
    const auto plat_a_consumption_metric = metric_rate(plat_a_rates.first);
    const auto plat_a_import_metric = metric_rate(plat_a_rates.second);
    const auto field_consumption_metric = metric_rate(field_rates.first);
    const auto field_import_metric = metric_rate(field_rates.second);

    // Expected values (following REIN efficiency pattern)
    const double expected_sub_b_consumption = 100.0;    // Own GCONSUMP only
    const double expected_sub_b_import = 50.0;

    const double expected_sub_d_consumption = 0.0;      // No GCONSUMP
    const double expected_sub_d_import = 0.0;

    // PLAT-A: own 30/15 + children (SUB-B: 100×0.8, SUB-D: 0×0.7) = 30+80+0/15+40+0 = 110/55
    const double expected_plat_a_consumption = 30.0 + (100.0 * 0.8) + (0.0 * 0.7);  // 110.0
    const double expected_plat_a_import = 15.0 + (50.0 * 0.8) + (0.0 * 0.7);        // 55.0

    // FIELD: own 20/10 + children (PLAT-A: 110×0.9) = 20+99/10+49.5 = 119/59.5
    const double expected_field_consumption = 20.0 + (110.0 * 0.9);  // 119.0
    const double expected_field_import = 10.0 + (55.0 * 0.9);        // 59.5

    // Test leaf group with own GCONSUMP
    checkRate(sub_b_consumption_metric, expected_sub_b_consumption);
    checkRate(sub_b_import_metric, expected_sub_b_import);

    // Test group without GCONSUMP
    checkRate(sub_d_consumption_metric, expected_sub_d_consumption);
    checkRate(sub_d_import_metric, expected_sub_d_import);

    // Test intermediate group (with both own GCONSUMP and children)
    checkRate(plat_a_consumption_metric, expected_plat_a_consumption);
    checkRate(plat_a_import_metric, expected_plat_a_import);

    // Test field-level accumulation (most complex)
    checkRate(field_consumption_metric, expected_field_consumption);
    checkRate(field_import_metric, expected_field_import);
}

BOOST_AUTO_TEST_SUITE_END()
