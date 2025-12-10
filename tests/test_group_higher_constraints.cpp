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
/*!
 * \file
 * \brief Test for efficiency factor handling in checkGroupConstraintsProd()
 *
 * This test verifies the behavior of GroupStateHelper::checkGroupConstraintsProd()
 * when checking higher-level group constraints with efficiency factors (GEFAC).
 *
 * Group hierarchy used in the test:
 * \code
 *                    FIELD (ORAT = 10000)
 *                      |
 *                    PLAT (GEFAC=0.9, FLD, has guide rate)
 *                    /   \
 *                   /     \
 *   MANI (GEFAC=0.8)       MANI2 (GEFAC=0.7)
 *   ORAT control           FLD control
 *        |                      |
 *     WELL-A                 WELL-B
 * \endcode
 *
 * The test exposes a bug where the "addback" calculation at local_reduction_level
 * uses the accumulated efficiency factor (E_MANI * E_PLAT = 0.72) instead of
 * the partial efficiency factor (E_MANI = 0.8). This causes incorrect constraint
 * violation detection when MANI's production should actually be within limits.
 *
 * Key concepts tested:
 * - local_reduction_level: Set at PLAT because it has a guide rate
 * - Reduction rates: PLAT's reduction = E_MANI * MANI_rate (groups with individual control)
 * - Addback: Should restore MANI's contribution using E_MANI, not accumulated efficiency
 * - Guide rate fractions: How MANI and MANI2 share PLAT's allocation
 */

#define BOOST_TEST_MODULE GroupHigherConstraints
#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>

#include <boost/test/unit_test.hpp>

namespace Opm::Properties::TTag {
    struct TestEfficiencyFactorTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

//! \brief Test fixture with unit conversion helpers and tolerance checks
struct ToleranceAndUnitFixture
{
    //! Relative tolerance for exact mathematical relationships (e.g., efficiency factors)
    static constexpr double tight_tol = 1e-12;

    //! Relative tolerance for rate values with minor numerical differences
    static constexpr double rate_tol = 1e-8;

    //! Relative tolerance for algorithm results (e.g., scale factors)
    static constexpr double algo_tol = 1e-3;

    //! Convert SI rate (m³/s) to metric rate (SM3/day)
    static double metric_rate(double si_rate)
    {
        using namespace Opm::unit;
        return convert::to(si_rate, cubic(meter) / day);
    }

    //! Convert metric rate (SM3/day) to SI rate (m³/s)
    static double si_rate(double metric_rate)
    {
        using namespace Opm::unit;
        return convert::from(metric_rate, cubic(meter) / day);
    }

    //! Check that two values are close within tight tolerance (for exact relationships)
    static void checkClose(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, tight_tol);
    }

    //! Check that two rate values are close
    static void checkRate(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, rate_tol);
    }

    //! Check that two algorithm results are close (looser tolerance)
    static void checkAlgo(double actual, double expected)
    {
        BOOST_CHECK_CLOSE(actual, expected, algo_tol);
    }
};

} // anonymous namespace

BOOST_FIXTURE_TEST_SUITE(GroupHigherConstraintsTests, ToleranceAndUnitFixture)

BOOST_AUTO_TEST_CASE(TestGroupHigherConstraints)
{
    using TypeTag = Opm::Properties::TTag::TestEfficiencyFactorTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;

    const std::string filename = "GROUP_HIGHER_CONSTRAINTS.DATA";
    auto simulator = Opm::initSimulator<TypeTag>(filename.data());

    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(Opm::unit::day);
    simulator->model().newtonMethod().setIterationIndex(0);

    WellModel& well_model = simulator->problem().wellModel();
    const auto& schedule = simulator->vanguard().schedule();
    const int report_step_idx = 0;
    const auto& pu = well_model.phaseUsage();

    const int oil_phase_pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);

    well_model.beginReportStep(report_step_idx);

    // We want predictable rates that exceed FIELD's target to trigger constraint violation:
    //   WELL-A: 10000 SM3/day oil (negative for production)
    //   WELL-B: 4000 SM3/day oil (negative for production)
    // Total at FIELD level = E_PLAT * (E_MANI * 10000 + E_MANI2 * 4000)
    //                      = 0.9 * (0.8 * 10000 + 0.7 * 4000) = 0.9 * 10800 = 9720 SM3/day
    // This is close to FIELD's ORAT of 10000, so MANI's share will be constrained.
    const double well_a_oil_rate = si_rate(-10000.0);  // m³/s (negative = production)
    const double well_b_oil_rate = si_rate(-4000.0);   // m³/s

    auto& ws = well_model.wellState();

    // Set WELL-A rates
    {
        auto well_idx = ws.index("WELL-A");
        BOOST_REQUIRE(well_idx.has_value());
        auto& well_a = ws.well(well_idx.value());
        std::fill(well_a.surface_rates.begin(), well_a.surface_rates.end(), 0.0);
        well_a.surface_rates[oil_phase_pos] = well_a_oil_rate;
        well_a.production_cmode = Opm::Well::ProducerCMode::GRUP;  // Under group control
    }

    // Set WELL-B rates
    {
        auto well_idx = ws.index("WELL-B");
        BOOST_REQUIRE(well_idx.has_value());
        auto& well_b = ws.well(well_idx.value());
        std::fill(well_b.surface_rates.begin(), well_b.surface_rates.end(), 0.0);
        well_b.surface_rates[oil_phase_pos] = well_b_oil_rate;
        well_b.production_cmode = Opm::Well::ProducerCMode::GRUP;  // Under group control
    }

    // After setting well control modes, we need to update the global well info
    // so that isProductionGrup() returns true for wells under group control.
    // This is needed for FractionCalculator::guideRateSum() to include wells.
    std::vector<Opm::Well::Status> well_status(ws.size(), Opm::Well::Status::OPEN);
    const auto& comm = simulator->vanguard().grid().comm();
    ws.updateGlobalIsGrup(comm, well_status);

    // Set up guide rates for wells to control how rates are distributed.
    // Without this, localFraction() falls back to potentials or returns 1.0.
    // We use guide rates proportional to well rates for simplicity.
    //
    // Guide rate values are in output units (SM3/day for metric), not SI units.
    // The init_grvalue() method stores values directly without conversion.
    auto& guide_rate = well_model.guideRate();
    const double sim_time = 0.0;

    // Set guide rates for wells (values in SM3/day, same as their production rates)
    guide_rate.init_grvalue(report_step_idx, "WELL-A",
        Opm::GuideRate::GuideRateValue{sim_time, 10000.0, Opm::GuideRateModel::Target::OIL});
    guide_rate.init_grvalue(report_step_idx, "WELL-B",
        Opm::GuideRate::GuideRateValue{sim_time, 4000.0, Opm::GuideRateModel::Target::OIL});

    // Set guide rate for PLAT group - this is needed to trigger local_reduction_level
    // at PLAT level in checkGroupConstraintsProd(). Without this, the algorithm
    // doesn't set local_reduction_level and the addback logic at line checkGroupConstraintsProd() is skipped.
    guide_rate.init_grvalue(report_step_idx, "PLAT",
        Opm::GuideRate::GuideRateValue{sim_time, 9000.0, Opm::GuideRateModel::Target::OIL});

    // The checkGroupConstraintsProd function uses production_rates() to compute
    // current_rate_fraction when checking guide rate violations. We need to
    // populate these for the groups in the chain.
    auto& gs = well_model.groupState();
    const int num_phases = pu.numActivePhases();

    // Get group objects and their efficiency factors
    const auto& mani_group = schedule.getGroup("MANI", report_step_idx);
    const auto& mani2_group = schedule.getGroup("MANI2", report_step_idx);
    const auto& plat_group = schedule.getGroup("PLAT", report_step_idx);
    const auto& field_group = schedule.getGroup("FIELD", report_step_idx);

    const double E_MANI = mani_group.getGroupEfficiencyFactor();
    const double E_MANI2 = mani2_group.getGroupEfficiencyFactor();
    const double E_PLAT = plat_group.getGroupEfficiencyFactor();

    // Set production rates for groups (these are the "current" group rates)
    // MANI produces 10000 SM3/day, but with efficiency factor this shows up as:
    // MANI: 10000 SM3/day (direct well rate)
    // PLAT: E_MANI * 10000 + E_MANI2 * 4000 = 8000 + 2800 = 10800 SM3/day
    // FIELD: E_PLAT * 10800 = 0.9 * 10800 = 9720 SM3/day
    std::vector<double> mani_prod_rates(num_phases, 0.0);
    mani_prod_rates[oil_phase_pos] = -well_a_oil_rate;  // 10000 SM3/day (positive)

    std::vector<double> mani2_prod_rates(num_phases, 0.0);
    mani2_prod_rates[oil_phase_pos] = -well_b_oil_rate;  // 4000 SM3/day (positive)

    std::vector<double> plat_prod_rates(num_phases, 0.0);
    plat_prod_rates[oil_phase_pos] = E_MANI * (-well_a_oil_rate) + E_MANI2 * (-well_b_oil_rate);  // 10800 SM3/day

    std::vector<double> field_prod_rates(num_phases, 0.0);
    field_prod_rates[oil_phase_pos] = E_PLAT * plat_prod_rates[oil_phase_pos];  // 9720 SM3/day

    gs.update_production_rates("MANI", mani_prod_rates);
    gs.update_production_rates("MANI2", mani2_prod_rates);
    gs.update_production_rates("PLAT", plat_prod_rates);
    gs.update_production_rates("FIELD", field_prod_rates);

    // MANI has individual ORAT control (this triggers checkGroupConstraintsProd)
    gs.production_control("MANI", Opm::Group::ProductionCMode::ORAT);

    // MANI2 has FLD control (follows parent)
    gs.production_control("MANI2", Opm::Group::ProductionCMode::FLD);

    // PLAT has FLD control (follows parent FIELD)
    gs.production_control("PLAT", Opm::Group::ProductionCMode::FLD);

    // FIELD has ORAT control
    gs.production_control("FIELD", Opm::Group::ProductionCMode::ORAT);

    Opm::DeferredLogger deferred_logger;
    auto& gsh = well_model.groupStateHelper();

    // Update the groupControlledWells count based on well control modes
    gsh.updateGroupControlledWells(/*is_production_group=*/true, Opm::Phase::OIL, deferred_logger);

    gsh.updateGroupTargetReduction(field_group, /*is_injector=*/false);

    // ========================================================================
    // VERIFY: Check the test setup
    // ========================================================================

    // Verify efficiency factors from schedule
    checkClose(E_MANI, 0.8);
    checkClose(E_MANI2, 0.7);
    checkClose(E_PLAT, 0.9);
    const double accumulated_eff = E_MANI * E_PLAT;
    checkClose(accumulated_eff, 0.72);

    // Verify well rates in WellState
    checkRate(metric_rate(ws.well("WELL-A").surface_rates[oil_phase_pos]), -10000.0);
    checkRate(metric_rate(ws.well("WELL-B").surface_rates[oil_phase_pos]), -4000.0);

    // Verify wells are marked as under group control
    BOOST_CHECK(ws.isProductionGrup("WELL-A"));
    BOOST_CHECK(ws.isProductionGrup("WELL-B"));

    // Verify guide rates are set
    BOOST_CHECK(guide_rate.has("WELL-A"));
    BOOST_CHECK(guide_rate.has("WELL-B"));
    BOOST_CHECK(guide_rate.has("PLAT"));

    // Verify groupControlledWells count (1 well per leaf group)
    BOOST_CHECK_EQUAL(gsh.groupControlledWells("MANI", "", true, Opm::Phase::OIL), 1);
    BOOST_CHECK_EQUAL(gsh.groupControlledWells("MANI2", "", true, Opm::Phase::OIL), 1);
    BOOST_CHECK_EQUAL(gsh.groupControlledWells("PLAT", "", true, Opm::Phase::OIL), 1);

    // Verify sumWellSurfaceRates for each group
    // MANI: direct well rate = 10000 SM3/day
    // MANI2: direct well rate = 4000 SM3/day
    // PLAT: E_MANI * 10000 + E_MANI2 * 4000 = 0.8 * 10000 + 0.7 * 4000 = 10800 SM3/day
    const double sum_mani = gsh.sumWellSurfaceRates(mani_group, oil_phase_pos, /*injector=*/false);
    const double sum_mani2 = gsh.sumWellSurfaceRates(mani2_group, oil_phase_pos, /*injector=*/false);
    const double sum_plat = gsh.sumWellSurfaceRates(plat_group, oil_phase_pos, /*injector=*/false);
    checkRate(metric_rate(sum_mani), 10000.0);
    checkRate(metric_rate(sum_mani2), 4000.0);
    checkRate(metric_rate(sum_plat), E_MANI * 10000.0 + E_MANI2 * 4000.0);  // 10800

    // Verify FIELD's ORAT target
    const auto& field_group_ctrl = field_group.productionControls(gsh.summaryState());
    checkRate(metric_rate(field_group_ctrl.oil_target), 10000.0);

    // Verify group reduction rates
    // FIELD: 0 (PLAT has guide rate, so doesn't contribute to FIELD's reduction)
    // PLAT: E_MANI * 10000 = 8000 (MANI has individual ORAT control)
    // MANI, MANI2: 0 (leaf groups with no sub-groups having individual control)
    const auto& field_red = gs.production_reduction_rates("FIELD");
    const auto& plat_red = gs.production_reduction_rates("PLAT");
    BOOST_CHECK_SMALL(metric_rate(field_red[oil_phase_pos]), tight_tol);
    checkRate(metric_rate(plat_red[oil_phase_pos]), E_MANI * 10000.0);  // 8000

    // Verify guide rate fraction for MANI
    // MANI's guide rate = E_MANI * well_guide_rate = 0.8 * 10000 = 8000
    // MANI2's guide rate = E_MANI2 * well_guide_rate = 0.7 * 4000 = 2800
    // Total = 10800, MANI fraction = 8000/10800 = 20/27
    const double expected_mani_fraction = (E_MANI * 10000.0) / (E_MANI * 10000.0 + E_MANI2 * 4000.0);
    checkRate(expected_mani_fraction, 20.0 / 27.0);  // 0.740740740...

    // Prepare rates array for MANI (negative for production)
    // This represents the "current rate" that will be compared against the target
    std::vector<double> mani_rates(num_phases, 0.0);
    mani_rates[oil_phase_pos] = well_a_oil_rate;  // -10000 SM3/day in SI units

    // Prepare resv_coeff (reservoir volume coefficients) - use 1.0 for simplicity
    std::vector<double> resv_coeff(num_phases, 1.0);

    // Initial efficiency factor is E_MANI (will accumulate E_PLAT during recursion)
    const double initial_efficiency = E_MANI;

    // Call checkGroupConstraintsProd
    // Arguments: name, parent, group, rates, efficiency_factor, resv_coeff, check_guide_rate, deferred_logger
    auto [constraint_violated, scale] = gsh.checkGroupConstraintsProd(
        "MANI",                    // name (the group/well being checked)
        "PLAT",                    // parent (immediate parent group name)
        plat_group,                // group (PLAT - will recurse to FIELD since PLAT has FLD)
        mani_rates.data(),         // rates (MANI's production rates)
        initial_efficiency,        // efficiency_factor (E_MANI, will accumulate E_PLAT)
        resv_coeff,                // resv_coeff
        /*check_guide_rate=*/true, // check_guide_rate
        deferred_logger            // deferred_logger
    );

    // ========================================================================
    // EXPECTED RESULT (with current buggy implementation)
    // ========================================================================
    //
    // The checkGroupConstraintsProd algorithm proceeds as follows:
    //
    // Setup:
    //   - Chain: [FIELD, PLAT, MANI]
    //   - efficiency_factor = E_MANI * E_PLAT = 0.8 * 0.9 = 0.72
    //   - current_rate_available = 10000 SM3/day (MANI's production)
    //   - orig_target = 10000 SM3/day (FIELD's ORAT)
    //   - local_reduction_level = 1 (PLAT has guide rate)
    //
    // Loop iteration ii=0 (FIELD):
    //   target = 10000
    //   local_reduction("FIELD") = 0 (PLAT has guide rate)
    //   target = 10000 - 0 = 10000
    //   local_fraction("PLAT") = 1.0 (only child of FIELD)
    //   target = 10000 * 1.0 = 10000
    //
    // Loop iteration ii=1 (PLAT):
    //   target = 10000
    //   local_reduction("PLAT") = E_MANI * 10000 = 8000 (MANI has ORAT control)
    //   target = 10000 - 8000 = 2000
    //
    //   BUG: Addback uses accumulated efficiency (0.72) instead of E_MANI (0.8):
    //   addback = current_rate * efficiency_factor = 10000 * 0.72 = 7200
    //   target = 2000 + 7200 = 9200
    //
    //   (Correct addback would be: 10000 * 0.8 = 8000, giving target = 10000)
    //
    //   local_fraction("MANI") = 8000/10800 = 0.740741
    //   target = 9200 * 0.740741 = 6814.81
    //
    // Final calculation (line 413):
    //   target_rate_available = target / efficiency_factor
    //                        = (184000/27) / 0.72 = 9465.0205761... SM3/day
    //
    // Scale factor:
    //   scale = target_rate_available / current_rate_available
    //         = 9465.0205761... / 10000 = 0.94650205761...
    //
    // Because scale < 1, the buggy implementation incorrectly constrains MANI,
    // reducing its production even though the correct calculation would give
    // scale = 1.0288 (no constraint violation).
    //
    // Exact value: 4600000 / 486 / 10000 = 2300/243 / 5 = 460/243
    const double expected_scale_with_bug = 460.0 / 486.0;  // 0.9465020576...
    BOOST_CHECK(constraint_violated);  // Bug causes false constraint violation
    checkAlgo(scale, expected_scale_with_bug);

    // Derived values for verification
    const double current_rate = 10000.0;  // SM3/day
    const double target_rate_available_with_bug = scale * current_rate;
    checkAlgo(target_rate_available_with_bug, 4600000.0 / 486.0);  // 9465.0205761...
}

BOOST_AUTO_TEST_SUITE_END()
