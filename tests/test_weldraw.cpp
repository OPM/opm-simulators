// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2026 Equinor ASA

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
 * \brief Tests for the WELDRAW keyword (maximum allowable drawdown).
 *
 * The drawdown limit Dmax is converted into a maximum production rate
 * for the selected phase,
 *
 *     Qmax = Dmax * sum_j (Tw_j * M_j),
 *
 * which is then imposed as an additional LRAT (or GRAT) constraint.  The
 * tests below cover the three properties that the implementation has to
 * get right: that the derived rate is proportional to the drawdown limit,
 * that it only ever tightens the well's controls, and that it is frozen
 * once the NUPCOL iteration window has passed.
 */

#define BOOST_TEST_MODULE Weldraw
#include "SimulatorFixture.hpp"
#include "ToleranceAndUnitFixture.hpp"

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <boost/test/unit_test.hpp>

#include <memory>
#include <string>

namespace Opm::Properties::TTag {
    struct TestWeldrawTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

struct WeldrawTestFixture : ToleranceAndUnitFixture
{
    using TypeTag = Opm::Properties::TTag::TestWeldrawTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;

    //! Boot the simulator on the WELDRAW deck and bring the well model up to
    //! the point where well controls are updated.
    void setup(const std::string& data_file = "WELDRAW.DATA")
    {
        simulator_ = Opm::initSimulator<TypeTag>(data_file.data(), "test_weldraw");

        simulator_->model().applyInitialSolution();
        simulator_->setEpisodeIndex(-1);
        simulator_->setEpisodeLength(0.0);
        simulator_->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
        simulator_->setTimeStepSize(Opm::unit::day);
        simulator_->problem().resetIterationForNewTimestep();

        wellModel().beginReportStep(report_step_idx_);
        wellModel().beginTimeStep();
    }

    //! One pass of the control update, i.e. the call chain that recomputes
    //! the WELDRAW rate limits.
    void updateWellControls()
    {
        auto logger_guard = wellModel().groupStateHelper().pushLogger();
        auto& deferred_logger = wellModel().groupStateHelper().deferredLogger();
        wellModel().calculateExplicitQuantities();
        wellModel().prepareTimeStep(deferred_logger);
        wellModel().updateWellControls(deferred_logger);
    }

    WellModel& wellModel() { return simulator_->problem().wellModel(); }
    auto& wellState() { return wellModel().wellState(); }

    const Opm::SummaryState& summaryState() const
    { return simulator_->vanguard().summaryState(); }

    int nupcol() const
    { return simulator_->vanguard().schedule()[report_step_idx_].nupcol(); }

    //! The state of a well by name.
    auto& ws(const std::string& well_name)
    {
        const auto index = wellState().index(well_name);
        BOOST_REQUIRE_MESSAGE(index.has_value(), "No well state for " << well_name);
        return wellState().well(*index);
    }

    //! The production controls a well is actually solved against, i.e. with
    //! the drawdown limit applied.
    Opm::Well::ProductionControls controlsWithWeldraw(const std::string& well_name)
    {
        return wellModel().getWell(well_name)
            .productionControlsWithWeldraw(summaryState(), ws(well_name));
    }

    //! The drawdown-derived maximum rate of a well; fails if none was set.
    double weldrawMaxRate(const std::string& well_name)
    {
        const auto& rate = ws(well_name).weldraw_max_rate;
        BOOST_REQUIRE_MESSAGE(rate.has_value(),
                              "No WELDRAW rate limit for " << well_name);
        return *rate;
    }

    std::unique_ptr<Opm::Simulator<TypeTag>> simulator_;
    int report_step_idx_{0};
};

} // anonymous namespace

BOOST_FIXTURE_TEST_SUITE(WeldrawTests, WeldrawTestFixture)

// The derived rate limit is proportional to the drawdown limit: PROD-B has
// twice the drawdown of PROD-A and is otherwise identical, so its rate limit
// must be exactly twice as large.  This pins down the Qmax = Dmax * sum(Tw*M)
// relation without duplicating the mobility calculation in the test.
BOOST_AUTO_TEST_CASE(RateLimitScalesWithDrawdown)
{
    setup();
    updateWellControls();

    const auto rate_a = weldrawMaxRate("PROD-A");
    const auto rate_b = weldrawMaxRate("PROD-B");

    BOOST_CHECK_GT(rate_a, 0.0);
    checkClose(rate_b, 2.0 * rate_a);
}

// A well whose deck targets are looser than the drawdown limit is put on an
// LRAT constraint carrying exactly the derived rate.
BOOST_AUTO_TEST_CASE(DrawdownLimitImposedAsLiquidRate)
{
    setup();
    updateWellControls();

    const auto controls = controlsWithWeldraw("PROD-A");

    BOOST_CHECK(controls.hasControl(Opm::Well::ProducerCMode::LRAT));
    checkClose(controls.liquid_rate, weldrawMaxRate("PROD-A"));

    // The well is reported as being under the drawdown limit.
    BOOST_REQUIRE(ws("PROD-A").weldraw_cmode.has_value());
    BOOST_CHECK(*ws("PROD-A").weldraw_cmode == Opm::Well::ProducerCMode::LRAT);
}

// PROD-C's deck LRAT target is far below what its 400 bar drawdown limit
// allows, so the limit must not loosen it.
BOOST_AUTO_TEST_CASE(TighterDeckTargetIsNotLoosened)
{
    setup();
    updateWellControls();

    const auto deck_target = si_rate(5.0);
    BOOST_CHECK_LT(deck_target, weldrawMaxRate("PROD-C"));

    const auto controls = controlsWithWeldraw("PROD-C");
    BOOST_CHECK(controls.hasControl(Opm::Well::ProducerCMode::LRAT));
    checkRate(metric_rate(controls.liquid_rate), 5.0);

    // ... and the well is not reported as drawdown limited.
    BOOST_CHECK(!ws("PROD-C").weldraw_cmode.has_value());
}

// The rate limit is recomputed while inside the NUPCOL window and frozen
// afterwards, mirroring how group targets behave.  A poisoned value is used
// to distinguish "recomputed" from "left alone".
BOOST_AUTO_TEST_CASE(RateLimitFrozenAfterNupcol)
{
    setup();
    updateWellControls();

    const auto physical_rate = weldrawMaxRate("PROD-A");
    constexpr double poison = 12345.0;

    // Still within NUPCOL: the poisoned value is overwritten.
    ws("PROD-A").weldraw_max_rate = poison;
    updateWellControls();
    checkClose(weldrawMaxRate("PROD-A"), physical_rate);

    // Step past the NUPCOL window.
    for (int i = 0; i < nupcol(); ++i) {
        simulator_->problem().advanceIteration();
    }
    BOOST_REQUIRE(!simulator_->problem().iterationContext().withinNupcol(nupcol()));

    // Beyond NUPCOL the stored limit must survive untouched.
    ws("PROD-A").weldraw_max_rate = poison;
    updateWellControls();
    BOOST_CHECK_EQUAL(weldrawMaxRate("PROD-A"), poison);
}

BOOST_AUTO_TEST_SUITE_END()
