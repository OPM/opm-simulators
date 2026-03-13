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
 * \brief Test for efficiency factor handling, addback and guide rate fractions in group constraint checking
 *
 *
 * Group hierarchy used in the test:
 * \code
 *                          FIELD (ORAT = 10000)
 *                          /           \
 *                         /             \
 *         PLAT (GEFAC=0.9)              PLAT-2 (GEFAC=0.85)
 *              /      \                       |
 *             /        \                   WELL-C
 *   MANI (0.8)       MANI-2 (0.75)
 *       /   \            |
 *  WELL-A  WELL-A2    WELL-B
 * \endcode
 *
 * Guide rates (explicit from GCONPROD):
 *   MANI:   8000, MANI-2: 1500, PLAT: 9500, PLAT-2: 2550
 *
 * Fractions:
 *   MANI at PLAT:    8000 / 9500 ≈ 0.842
 *   PLAT at FIELD:   9500 / 12050 ≈ 0.788
 *
 */

#define BOOST_TEST_MODULE GroupHigherConstraints
#include "SimulatorFixture.hpp"
#include "ToleranceAndUnitFixture.hpp"

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/FractionCalculator.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>

#include <boost/test/unit_test.hpp>

namespace Opm::Properties::TTag {
    struct TestEfficiencyFactorTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

// Using shared ToleranceAndUnitFixture from ToleranceAndUnitFixture.hpp

//! \brief Test fixture for group constraint tests to avoid code duplication
//!
//! Tests can customize behavior via configuration options before calling setup().
struct GroupConstraintTestFixture : ToleranceAndUnitFixture
{
    using TypeTag = Opm::Properties::TTag::TestEfficiencyFactorTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;

    void commonSetupForAllTests(
        bool plat_has_guide_rate,
        Opm::Group::ProductionCMode mani_control,
        Opm::Well::ProducerCMode well_a_control,
        Opm::Well::ProducerCMode well_a2_control,
        Opm::Group::ProductionCMode mani2_control = Opm::Group::ProductionCMode::FLD,
        const std::string& data_file = "GROUP_HIGHER_CONSTRAINTS.DATA",
        double well_a_rate_metric = 6000.0,
        double well_a2_rate_metric = 4000.0)
    {
        plat_has_guide_rate_ = plat_has_guide_rate;
        mani_control_ = mani_control;
        mani2_control_ = mani2_control;
        well_a_control_ = well_a_control;
        well_a2_control_ = well_a2_control;
        data_file_ = data_file;
        well_a_rate_metric_ = well_a_rate_metric;
        well_a2_rate_metric_ = well_a2_rate_metric;

        setupSimulator_();
        setupWellRates_();
        setupGroupState_();
        setupGuideRates_();
        finalizeSetup_();
    }

    // ========================================================================
    // Accessors
    // ========================================================================

    double wellARateMetric() const { return well_a_rate_metric_; }
    double wellA2RateMetric() const { return well_a2_rate_metric_; }
    double wellBRateMetric() const { return well_b_rate_metric_; }
    double wellCRateMetric() const { return well_c_rate_metric_; }
    double wellARateSI() const { return si_rate(-well_a_rate_metric_); }
    double wellA2RateSI() const { return si_rate(-well_a2_rate_metric_); }
    double wellBRateSI() const { return si_rate(-well_b_rate_metric_); }
    double wellCRateSI() const { return si_rate(-well_c_rate_metric_); }
    WellModel& wellModel() { return simulator_->problem().wellModel(); }
    const Opm::Schedule& schedule() const { return simulator_->vanguard().schedule(); }
    auto& wellState() { return wellModel().wellState(); }
    auto& groupState() { return wellModel().groupState(); }
    Opm::GuideRate& guideRate() { return wellModel().guideRate(); }
    auto& groupStateHelper() { return wellModel().groupStateHelper(); }
    int reportStepIdx() const { return report_step_idx_; }
    int oilPhasePos() const { return oil_phase_pos_; }
    Opm::DeferredLogger& deferredLogger() { return deferred_logger_; }

    const Opm::Group& maniGroup() const { return schedule().getGroup("MANI", report_step_idx_); }
    const Opm::Group& mani2Group() const { return schedule().getGroup("MANI-2", report_step_idx_); }
    const Opm::Group& platGroup() const { return schedule().getGroup("PLAT", report_step_idx_); }
    const Opm::Group& plat2Group() const { return schedule().getGroup("PLAT-2", report_step_idx_); }
    const Opm::Group& fieldGroup() const { return schedule().getGroup("FIELD", report_step_idx_); }

    //! Efficiency factors from schedule
    double eMani() const { return maniGroup().getGroupEfficiencyFactor(); }
    double eMani2() const { return mani2Group().getGroupEfficiencyFactor(); }
    double ePlat() const { return platGroup().getGroupEfficiencyFactor(); }
    double ePlat2() const { return plat2Group().getGroupEfficiencyFactor(); }

    //! Get MANI's total rate (SI units, positive for production)
    double maniTotalRate() const { return -wellARateSI() - wellA2RateSI(); }

    //! Get MANI's total rate (metric, positive for production)
    double maniTotalRateMetric() const { return wellARateMetric() + wellA2RateMetric(); }

    //! Get MANI-2's total rate (metric, positive for production)
    double mani2TotalRateMetric() const { return wellBRateMetric(); }

    //! Get PLAT's total rate (metric, positive for production)
    double platTotalRateMetric() const { return eMani() * maniTotalRateMetric() + eMani2() * mani2TotalRateMetric(); }

    //! Guide rates (from GCONPROD via guideRate object)
    double maniGuideRateMetric() {
        const Opm::GuideRate::RateVector zero_rates{0.0, 0.0, 0.0};
        return guideRate().get("MANI", Opm::GuideRateModel::Target::OIL, zero_rates);
    }
    double mani2GuideRateMetric() {
        const Opm::GuideRate::RateVector zero_rates{0.0, 0.0, 0.0};
        return guideRate().get("MANI-2", Opm::GuideRateModel::Target::OIL, zero_rates);
    }
    double platGuideRateMetric() {
        const Opm::GuideRate::RateVector zero_rates{0.0, 0.0, 0.0};
        return guideRate().get("PLAT", Opm::GuideRateModel::Target::OIL, zero_rates);
    }
    double plat2GuideRateMetric() {
        const Opm::GuideRate::RateVector zero_rates{0.0, 0.0, 0.0};
        return guideRate().get("PLAT-2", Opm::GuideRateModel::Target::OIL, zero_rates);
    }

    // ========================================================================
    // Helper methods for common test operations
    // ========================================================================

    //! Prepare rates array for MANI group (for checkGroupConstraintsProd)
    std::vector<double> maniRatesArray() const
    {
        std::vector<double> rates(num_phases_, 0.0);
        rates[oil_phase_pos_] = wellARateSI() + wellA2RateSI();  // negative
        return rates;
    }

    //! Prepare rates array for WELL-A (for getWellGroupTargetProducer)
    std::vector<double> wellARatesArray() const
    {
        std::vector<double> rates(num_phases_, 0.0);
        rates[oil_phase_pos_] = wellARateSI();  // negative
        return rates;
    }

    //! Get default resv_coeff (all 1.0)
    std::vector<double> resvCoeff() const
    {
        return std::vector<double>(num_phases_, 1.0);
    }
private:
    //! Whether PLAT group has a guide rate (affects local_reduction_level)
    bool plat_has_guide_rate_{true};

    //! Control mode for WELL-A (GRUP or ORAT)
    Opm::Well::ProducerCMode well_a_control_{Opm::Well::ProducerCMode::GRUP};

    //! Control mode for WELL-A2 (GRUP or ORAT)
    Opm::Well::ProducerCMode well_a2_control_{Opm::Well::ProducerCMode::GRUP};

    //! Control mode for MANI group (ORAT or FLD)
    Opm::Group::ProductionCMode mani_control_{Opm::Group::ProductionCMode::ORAT};

    //! Control mode for MANI-2 group (FLD or ORAT)
    Opm::Group::ProductionCMode mani2_control_{Opm::Group::ProductionCMode::FLD};

    //! DATA file to load
    std::string data_file_{"GROUP_HIGHER_CONSTRAINTS.DATA"};

    // Well rates (SM3/day, metric units for readability)
    double well_a_rate_metric_{6000.0};
    double well_a2_rate_metric_{4000.0};
    double well_b_rate_metric_{2000.0};
    double well_c_rate_metric_{3000.0};

    std::unique_ptr<Opm::Simulator<TypeTag>> simulator_;

    int report_step_idx_{0};
    int num_phases_{0};
    int oil_phase_pos_{0};

    Opm::DeferredLogger deferred_logger_;

    // ========================================================================
    // Private setup methods
    // ========================================================================

    void setupSimulator_()
    {
        simulator_ = Opm::initSimulator<TypeTag>(data_file_.data());

        simulator_->model().applyInitialSolution();
        simulator_->setEpisodeIndex(-1);
        simulator_->setEpisodeLength(0.0);
        simulator_->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
        simulator_->setTimeStepSize(Opm::unit::day);
        // Reset iteration context so code using problem().iterationContext() sees first iteration
        simulator_->problem().resetIterationForNewTimestep();

        wellModel().beginReportStep(report_step_idx_);

        const auto& pu = wellModel().phaseUsage();
        num_phases_ = pu.numActivePhases();
        oil_phase_pos_ = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    }

    void setupWellRates_()
    {
        auto& ws = wellModel().wellState();

        // For simplicity, assume the same well potentials as rates
        // Set WELL-A: surface_rates (negative) and well_potentials (positive)
        {
            auto well_idx = ws.index("WELL-A");
            BOOST_REQUIRE(well_idx.has_value());
            auto& well = ws.well(well_idx.value());
            std::ranges::fill(well.surface_rates, 0.0);
            std::ranges::fill(well.well_potentials, 0.0);
            well.surface_rates[oil_phase_pos_] = wellARateSI();
            well.well_potentials[oil_phase_pos_] = -wellARateSI();  // positive
            well.production_cmode = well_a_control_;
        }

        // Set WELL-A2
        {
            auto well_idx = ws.index("WELL-A2");
            BOOST_REQUIRE(well_idx.has_value());
            auto& well = ws.well(well_idx.value());
            std::ranges::fill(well.surface_rates, 0.0);
            std::ranges::fill(well.well_potentials, 0.0);
            well.surface_rates[oil_phase_pos_] = wellA2RateSI();
            well.well_potentials[oil_phase_pos_] = -wellA2RateSI();  // positive
            well.production_cmode = well_a2_control_;
        }

        // Set WELL-B
        {
            auto well_idx = ws.index("WELL-B");
            BOOST_REQUIRE(well_idx.has_value());
            auto& well = ws.well(well_idx.value());
            std::ranges::fill(well.surface_rates, 0.0);
            std::ranges::fill(well.well_potentials, 0.0);
            well.surface_rates[oil_phase_pos_] = wellBRateSI();
            well.well_potentials[oil_phase_pos_] = -wellBRateSI();  // positive
            well.production_cmode = Opm::Well::ProducerCMode::GRUP;
        }

        // Set WELL-C
        {
            auto well_idx = ws.index("WELL-C");
            BOOST_REQUIRE(well_idx.has_value());
            auto& well = ws.well(well_idx.value());
            std::ranges::fill(well.surface_rates, 0.0);
            std::ranges::fill(well.well_potentials, 0.0);
            well.surface_rates[oil_phase_pos_] = wellCRateSI();
            well.well_potentials[oil_phase_pos_] = -wellCRateSI();  // positive
            well.production_cmode = Opm::Well::ProducerCMode::GRUP;
        }

        // Update global well info
        std::vector<Opm::Well::Status> well_status(ws.size(), Opm::Well::Status::OPEN);
        const auto& comm = simulator_->vanguard().grid().comm();
        ws.updateGlobalIsGrup(comm, well_status);
    }

    void setupGuideRates_()
    {
        // Use updateGuideRates() to compute guide rates from:
        // - Explicit values in GCONPROD for MANI (8000), MANI-2 (1500), PLAT (9500), PLAT-2 (2550)
        // - GUIDERAT formula for wells (guide_rate = oil_potential)
        //
        // The explicit GCONPROD values ensure predictable guide rate fractions:
        //   MANI fraction at PLAT = 8000 / 9500 ≈ 0.842
        //   PLAT fraction at FIELD = 9500 / 12050 ≈ 0.788
        const double sim_time = 0.0;
        wellModel().updateGuideRates(report_step_idx_, sim_time);

        // If test requires PLAT to NOT have a guide rate (to test local_reduction_level=0),
        // erase it after updateGuideRates() has set it from GCONPROD
        if (!plat_has_guide_rate_) {
            guideRate().erase("PLAT");
        }
    }

    void setupGroupState_()
    {
        auto& gs = wellModel().groupState();
        auto& gsh = wellModel().groupStateHelper();

        // Compute group production rates from well surface rates
        // This applies efficiency factors at each level (see sumWellPhaseRates):
        //   MANI:   sum of wells = -10000 (WELL-A + WELL-A2)
        //   MANI-2: sum of wells = -2000 (WELL-B)
        //   PLAT:   E_MANI * MANI + E_MANI2 * MANI-2 = 0.8 * -10000 + 0.75 * -2000 = -9500
        //   PLAT-2: sum of wells = -3000 (WELL-C, no child groups)
        //   FIELD:  E_PLAT * PLAT + E_PLAT2 * PLAT-2 = 0.9 * -9500 + 0.85 * -3000 = -11100
        gsh.updateGroupProductionRates(fieldGroup());

        // Set group controls
        gs.production_control("MANI", mani_control_);
        gs.production_control("MANI-2", mani2_control_);
        gs.production_control("PLAT", Opm::Group::ProductionCMode::FLD);
        gs.production_control("PLAT-2", Opm::Group::ProductionCMode::FLD);
        gs.production_control("FIELD", Opm::Group::ProductionCMode::ORAT);
    }

    void finalizeSetup_()
    {
        auto& gsh = wellModel().groupStateHelper();
        const auto& field_group = schedule().getGroup("FIELD", report_step_idx_);

        gsh.updateGroupControlledWells(/*is_production_group=*/true, Opm::Phase::OIL);
        gsh.updateGroupTargetReduction(field_group, /*is_injector=*/false);
    }
};

} // anonymous namespace

BOOST_FIXTURE_TEST_SUITE(GroupHigherConstraintsTests, GroupConstraintTestFixture)

//! \brief Test checkGroupConstraintsProd() for group MANI with a guide rate at PLAT level between MANI and FIELD
//!
//! When there's a guide rate at PLAT level between MANI and FIELD, local_reduction_level is 1
//! and the addback happens at PLAT level. This test verifies that an addback efficiency factor of 0.8 is
//! applied to PLAT (local_reduction_level = 1) and not 0.8 * 0.9 = 0.72. (original version of
//! checkGroupConstraintsProd() had this bug)
//! \code
//!   (GCW = Group Controlled Wells, IDV = Individual control)
//!
//!   FIELD (ORAT = 10000, GCW=2)
//!   ├── PLAT (E=0.9, FLD, GUIDERATE=9500, GCW=1, local_reduction_level=1)
//!   │   ├── MANI (E=0.8, IDV, ORAT, GUIDERATE=8000, GCW=2)
//!   │   │   ├── WELL-A (RATE=6000, POT=6000, GRUP)
//!   │   │   └── WELL-A2 (RATE=4000, POT=4000, GRUP)
//!   │   └── MANI-2 (E=0.75, FLD, GUIDERATE=1500, GCW=1)
//!   │       └── WELL-B (RATE=2000, POT=2000, GRUP)
//!   └── PLAT-2 (E=0.85, FLD, GUIDERATE=2550, GCW=1)
//!       └── WELL-C (RATE=3000, POT=3000, GRUP)
//! \endcode
BOOST_AUTO_TEST_CASE(TestGroupHigherConstraintsWithGuideRate)
{
    // Configuration: Guide rate for PLAT (local_reduction_level = 1)
    commonSetupForAllTests(/*plat_has_guide_rate=*/true,
        /*mani_control=*/Opm::Group::ProductionCMode::ORAT,
        /*well_a_control=*/Opm::Well::ProducerCMode::GRUP,
        /*well_a2_control=*/Opm::Well::ProducerCMode::GRUP);

    // Verify the test setup
    auto& gsh = groupStateHelper();
    auto& ws = wellState();
    auto& gs = groupState();

    // Verify efficiency factors from schedule
    checkClose(eMani(), 0.8);
    checkClose(eMani2(), 0.75);
    checkClose(ePlat(), 0.9);
    checkClose(ePlat2(), 0.85);
    const double accumulated_eff = eMani() * ePlat();
    checkClose(accumulated_eff, 0.72);

    // Verify well rates in WellState
    checkRate(metric_rate(ws.well("WELL-A").surface_rates[oilPhasePos()]), -wellARateMetric()); // -6000.0
    checkRate(metric_rate(ws.well("WELL-A2").surface_rates[oilPhasePos()]), -wellA2RateMetric()); // -4000.0
    checkRate(metric_rate(ws.well("WELL-B").surface_rates[oilPhasePos()]), -wellBRateMetric()); // -2000.0
    checkRate(metric_rate(ws.well("WELL-C").surface_rates[oilPhasePos()]), -wellCRateMetric()); // -3000.0

    // Verify well potentials (stored in SI, convert to metric for comparison)
    checkRate(metric_rate(ws.well("WELL-A").well_potentials[oilPhasePos()]), wellARateMetric()); // 6000.0
    checkRate(metric_rate(ws.well("WELL-A2").well_potentials[oilPhasePos()]), wellA2RateMetric()); // 4000.0
    checkRate(metric_rate(ws.well("WELL-B").well_potentials[oilPhasePos()]), wellBRateMetric()); // 2000.0
    checkRate(metric_rate(ws.well("WELL-C").well_potentials[oilPhasePos()]), wellCRateMetric()); // 3000.0

    // Verify wells are marked as under group control
    BOOST_CHECK(ws.isProductionGrup("WELL-A"));
    BOOST_CHECK(ws.isProductionGrup("WELL-A2"));
    BOOST_CHECK(ws.isProductionGrup("WELL-B"));
    BOOST_CHECK(ws.isProductionGrup("WELL-C"));

    // Verify guide rates are set from GCONPROD via updateGuideRates()
    // Wells use GUIDERAT formula (guide_rate = oil_potential), groups use explicit GCONPROD values
    BOOST_CHECK(guideRate().has("WELL-A"));
    BOOST_CHECK(guideRate().has("WELL-A2"));
    BOOST_CHECK(guideRate().has("WELL-B"));
    BOOST_CHECK(guideRate().has("WELL-C"));
    BOOST_CHECK(guideRate().has("MANI"));
    BOOST_CHECK(guideRate().has("MANI-2"));
    BOOST_CHECK(guideRate().has("PLAT"));
    BOOST_CHECK(guideRate().has("PLAT-2"));

    // Verify group guide rates (explicit values from GCONPROD)
    const Opm::GuideRate::RateVector zero_rates{0.0, 0.0, 0.0};
    checkClose(guideRate().get("MANI", Opm::GuideRateModel::Target::OIL, zero_rates), maniGuideRateMetric());
    checkClose(guideRate().get("MANI-2", Opm::GuideRateModel::Target::OIL, zero_rates), mani2GuideRateMetric());
    checkClose(guideRate().get("PLAT", Opm::GuideRateModel::Target::OIL, zero_rates), platGuideRateMetric());
    checkClose(guideRate().get("PLAT-2", Opm::GuideRateModel::Target::OIL, zero_rates), plat2GuideRateMetric());

    // Verify local fractions using FractionCalculator
    // FractionCalculator only includes groups/wells that:
    //   1. Have FLD or NONE control (for groups) or GRUP control (for wells)
    //   2. Have groupControlledWells > 0
    {
        using FCalc = Opm::GroupStateHelpers::FractionCalculator<double, IndexTraits>;
        FCalc fcalc{schedule(), gsh, gsh.summaryState(), reportStepIdx(),
                    &guideRate(), Opm::GuideRateModel::Target::OIL,
                    /*is_producer=*/true, /*injection_phase=*/Opm::Phase::OIL};

        auto localFraction = [&](const std::string& name) {
            return fcalc.localFraction(name, /*always_included_child=*/"");
        };

        // Well fractions within MANI (both wells under GRUP control):
        //   WELL-A guide rate = oil_potential = 6000
        //   WELL-A2 guide rate = oil_potential = 4000
        //   WELL-A fraction = 6000 / (6000 + 4000) = 0.6
        //   WELL-A2 fraction = 4000 / (6000 + 4000) = 0.4
        const double well_a_frac = localFraction("WELL-A");
        const double well_a2_frac = localFraction("WELL-A2");
        checkClose(well_a_frac, wellARateMetric() / (wellARateMetric() + wellA2RateMetric())); // 0.6
        checkClose(well_a2_frac, wellA2RateMetric() / (wellARateMetric() + wellA2RateMetric())); // 0.4
        checkClose(well_a_frac + well_a2_frac, 1.0);

        // WELL-B fraction within MANI-2 (only well):
        const double well_b_frac = localFraction("WELL-B");
        checkClose(well_b_frac, 1.0);  // Only well in MANI-2

        // Group fractions within PLAT:
        //   MANI has ORAT control (not FLD), so it's excluded from parent's sum
        //   MANI-2 has FLD control and 1 group-controlled well (WELL-B) -> included
        //   Since MANI-2 is the only FLD child of PLAT, both return 1.0
        const double mani_frac = localFraction("MANI");
        const double mani2_frac = localFraction("MANI-2");
        checkClose(mani_frac, 1.0);   // MANI excluded (ORAT control), only MANI-2 active
        checkClose(mani2_frac, 1.0);  // Only active group at PLAT level

        // Group fractions at FIELD level:
        //   PLAT has FLD control and 1 group-controlled well (via MANI-2) -> included
        //   PLAT-2 has FLD control and 1 group-controlled well (WELL-C) -> included
        //   Both compete for rate allocation: PLAT=9500, PLAT-2=2550
        const double plat_frac = localFraction("PLAT");
        const double plat2_frac = localFraction("PLAT-2");
        // PLAT fraction = PLAT guide rate / (PLAT guide rate + PLAT-2 guide rate) = 9500/12050 = 0.788
        checkClose(plat_frac, platGuideRateMetric() / (platGuideRateMetric() + plat2GuideRateMetric()));
        // PLAT-2 fraction = PLAT-2 guide rate / (PLAT guide rate + PLAT-2 guide rate) = 2550/12050 = 0.212
        checkClose(plat2_frac, plat2GuideRateMetric() / (platGuideRateMetric() + plat2GuideRateMetric()));
        checkClose(plat_frac + plat2_frac, 1.0);

        // WELL-C fraction within PLAT-2 (only well):
        const double well_c_frac = localFraction("WELL-C");
        checkClose(well_c_frac, 1.0);  // Only well in PLAT-2
    }

    // Verify groupControlledWells count
    {
        const auto always_included_child = std::string("");
        const auto is_production_group = true;
        const auto injection_phase = Opm::Phase::OIL;
        BOOST_CHECK_EQUAL(gsh.groupControlledWells(
            "MANI", always_included_child, is_production_group, injection_phase), 2);
        BOOST_CHECK_EQUAL(gsh.groupControlledWells(
            "MANI-2", always_included_child, is_production_group, injection_phase), 1);
        BOOST_CHECK_EQUAL(gsh.groupControlledWells(
            "PLAT-2", always_included_child, is_production_group, injection_phase), 1);
        BOOST_CHECK_EQUAL(gsh.groupControlledWells(
            "PLAT", always_included_child, is_production_group, injection_phase), 1);  // via MANI-2
        BOOST_CHECK_EQUAL(gsh.groupControlledWells(
            "FIELD", always_included_child, is_production_group, injection_phase), 2);  // PLAT + PLAT-2
    }
    // Test sumWellSurfaceRates for each group
    {
        const double sum_mani = gsh.sumWellSurfaceRates(maniGroup(), oilPhasePos(), /*injector=*/false);
        const double sum_mani2 = gsh.sumWellSurfaceRates(mani2Group(), oilPhasePos(), /*injector=*/false);
        const double sum_plat2 = gsh.sumWellSurfaceRates(plat2Group(), oilPhasePos(), /*injector=*/false);
        const double sum_plat = gsh.sumWellSurfaceRates(platGroup(), oilPhasePos(), /*injector=*/false);
        // WELL-A + WELL-A2 = 6000.0 + 4000.0 = 10000.0
        checkRate(metric_rate(sum_mani), wellARateMetric() + wellA2RateMetric()); // 10000.0
        checkRate(metric_rate(sum_mani2), wellBRateMetric());  // WELL-B = 2000.0
        checkRate(metric_rate(sum_plat2), wellCRateMetric());  // WELL-C = 3000.0
        checkRate(metric_rate(sum_plat), platTotalRateMetric());  // E_MANI*10000 + E_MANI2*2000 = 9500.0
    }
    // Test FIELD's ORAT target
    {
        const auto& field_group_ctrl = fieldGroup().productionControls(gsh.summaryState());
        checkRate(metric_rate(field_group_ctrl.oil_target), 10000.0);
    }

    // Test production control modes from schedule (GCONPROD)
    {
        const auto& field_ctrl = fieldGroup().productionControls(gsh.summaryState());
        BOOST_CHECK(field_ctrl.cmode == Opm::Group::ProductionCMode::ORAT);
        const auto& plat_ctrl = platGroup().productionControls(gsh.summaryState());
        BOOST_CHECK(plat_ctrl.cmode == Opm::Group::ProductionCMode::FLD);
        const auto& plat2_ctrl = plat2Group().productionControls(gsh.summaryState());
        BOOST_CHECK(plat2_ctrl.cmode == Opm::Group::ProductionCMode::FLD);
        const auto& mani_ctrl = maniGroup().productionControls(gsh.summaryState());
        BOOST_CHECK(mani_ctrl.cmode == Opm::Group::ProductionCMode::ORAT);
        const auto& mani2_ctrl = mani2Group().productionControls(gsh.summaryState());
        BOOST_CHECK(mani2_ctrl.cmode == Opm::Group::ProductionCMode::FLD);
    }

    // Test production control modes from GroupState (set by setupGroupState())
    {
        BOOST_CHECK(gs.production_control("FIELD") == Opm::Group::ProductionCMode::ORAT);
        BOOST_CHECK(gs.production_control("PLAT") == Opm::Group::ProductionCMode::FLD);
        BOOST_CHECK(gs.production_control("PLAT-2") == Opm::Group::ProductionCMode::FLD);
        BOOST_CHECK(gs.production_control("MANI") == Opm::Group::ProductionCMode::ORAT);
        BOOST_CHECK(gs.production_control("MANI-2") == Opm::Group::ProductionCMode::FLD);
    }

    // Verify group reduction rates
    // PLAT: E_MANI * 10000 = 8000 (MANI has individual ORAT control, MANI-2 is FLD so doesn't add)
    // FIELD: 0 (PLAT has guide rate AND group-controlled wells via MANI-2, so reduction doesn't propagate)
    // Note: With MANI-2 under PLAT, PLAT now has group-controlled wells, so the else branch in
    //       updateGroupTargetReductionRecursive_ is taken. Since PLAT has a guide rate, its
    //       reduction is NOT accumulated to FIELD.
    {
        const auto& field_red = gs.production_reduction_rates("FIELD");
        checkRate(metric_rate(field_red[oilPhasePos()]), 0.0); // 0 because PLAT has guide rate and group-controlled wells
        const auto& plat_red = gs.production_reduction_rates("PLAT");
        checkRate(metric_rate(plat_red[oilPhasePos()]), eMani() * maniTotalRateMetric()); // 8000.0
        const auto& plat2_red = gs.production_reduction_rates("PLAT-2");
        checkRate(metric_rate(plat2_red[oilPhasePos()]), 0.0);
        const auto& mani_red = gs.production_reduction_rates("MANI");
        checkRate(metric_rate(mani_red[oilPhasePos()]), 0.0); // 0.0 since both WELL-A and WELL-A2 are under group control
        const auto& mani2_red = gs.production_reduction_rates("MANI-2");
        checkRate(metric_rate(mani2_red[oilPhasePos()]), 0.0); // 0.0 since WELL-B is under group control
    }


    // Call checkGroupConstraintsProd
    auto mani_rates = maniRatesArray();
    auto resv_coeff = resvCoeff();
    const double initial_efficiency = eMani();

    auto [constraint_violated, scale] = gsh.checkGroupConstraintsProd(
        "MANI",
        "PLAT",
        platGroup(),
        mani_rates.data(),
        initial_efficiency,
        resv_coeff,
        /*check_guide_rate=*/true
    );

    // ========================================================================
    // EXPECTED RESULT (with PLAT having guide rate -> local_reduction_level = 1)
    // ========================================================================
    //
    // The checkGroupConstraintsProd algorithm proceeds as follows:
    //
    // Setup:
    //   - Chain: [FIELD, PLAT, MANI]
    //   - efficiency_factor = E_MANI * E_PLAT = 0.8 * 0.9 = 0.72
    //   - current_rate_available = 10000 SM3/day (MANI's production)
    //   - orig_target = 10000 SM3/day (FIELD's ORAT)
    //   - local_reduction_level = 1 (PLAT is under FLD control AND has guide rate AND num_gr_ctrl > 0)
    //
    // Loop iteration ii=0 (FIELD):
    //   - reduction_fn("FIELD") = 0 (FIELD has no local reduction)
    //   - target = 10000 - 0 = 10000
    //   - No addback (ii != local_reduction_level)
    //   - PLAT fraction at FIELD = 9500 / (9500 + 2550) ≈ 0.788
    //   - target = 10000 * 0.788 ≈ 7884.65
    //
    // Loop iteration ii=1 (PLAT):
    //   - Reduction condition: local_reduction_level==1, ii==1, reduction applied
    //   - reduction_fn("PLAT") = E_MANI * 10000 = 8000 SM3/day (MANI has ORAT control)
    //   - target = 7884.65 - 8000 ≈ -115.35
    //
    //   - Addback (ii == local_reduction_level==1):
    //   - addback_efficiency = E_MANI = 0.8
    //   - target += current_rate * 0.8 = -115.35 + 10000 * 0.8 = -115.35 + 8000 ≈ 7884.65
    //
    //   - MANI fraction at PLAT = 8000 / (8000 + 1500) ≈ 0.842
    //   - target = 7884.65 * 0.842 ≈ 6640.07
    //
    // Final calculation:
    //   target_rate_available = target / efficiency_factor
    //                        = 6640.07 / 0.72 ≈ 9222.32
    //
    // Scale factor:
    //   scale = target_rate_available / current_rate_available
    //         = 9222.32 / 10000 ≈ 0.922
    //
    // The constraint is violated (scale < 1.0) because MANI's share of FIELD's target
    // (after applying PLAT and MANI fractions) is less than its current production.
    //
    // Expected: scale ≈ 0.922 (constraint violated)
    const double expected_scale = 0.92208390963577669;
    BOOST_CHECK(constraint_violated);  // Constraint IS violated with new hierarchy
    checkAlgo(scale, expected_scale);

    // Derived values for verification
    const double current_rate = 10000.0;  // SM3/day
    const double target_rate_available = scale * current_rate;
    checkAlgo(target_rate_available, expected_scale * current_rate);
}

//! \brief Test checkGroupConstraintsProd() for group MANI with no guide rate at PLAT level between MANI and FIELD
//!
//! In this test, since PLAT has no guide rate and not zero group-controlled wells (it has 1 group-controlled well,
//! namely WELL-B via MANI-2), local_reduction_level stays at 0 and the addback happens at FIELD level instead of
//! at PLAT level.
//! \code
//!   FIELD (ORAT = 10000, local_reduction_level=0)
//!   ├── PLAT (E=0.9, FLD, NO GUIDERATE, GCW=1)
//!   │   ├── MANI (E=0.8, IDV, ORAT, GUIDERATE=8000, GCW=2)
//!   │   │   ├── WELL-A (RATE=6000, POT=6000, GRUP)
//!   │   │   └── WELL-A2 (RATE=4000, POT=4000, GRUP)
//!   │   └── MANI-2 (E=0.75, FLD, GUIDERATE=1500, GCW=1)
//!   │       └── WELL-B (RATE=2000, POT=2000, GRUP)
//!   └── PLAT-2 (E=0.85, FLD, GUIDERATE=2550, GCW=1)
//!       └── WELL-C (RATE=3000, POT=3000, GRUP)
//! \endcode
BOOST_AUTO_TEST_CASE(TestGroupHigherConstraintsProdWithoutGuideRate)
{
    // Configuration: No guide rate for PLAT (local_reduction_level = 0)
    commonSetupForAllTests(
        /*plat_has_guide_rate=*/false,
        /*mani_control       =*/Opm::Group::ProductionCMode::ORAT,
        /*well_a_control     =*/Opm::Well::ProducerCMode::GRUP,
        /*well_a2_control    =*/Opm::Well::ProducerCMode::GRUP
    );

    auto& gsh = groupStateHelper();
    auto& gs = groupState();

    // Verify PLAT has no guide rate
    BOOST_CHECK(!guideRate().has("PLAT"));

    // Verify reduction rates
    // Without guide rate at PLAT, MANI's reduction bubbles up to FIELD
    // But PLAT still has group-controlled wells via MANI-2 (FLD control)
    // -> So PLAT goes to else branch in updateGroupTargetReductionRecursive_()
    // -> Since PLAT has no guide rate, MANI's reduction propagates to FIELD
    const auto& field_red = gs.production_reduction_rates("FIELD");
    // FIELD reduction = E_PLAT * E_MANI * 10000 (MANI has individual control)
    const double expected_field_reduction = ePlat() * eMani() * maniTotalRateMetric();
    checkRate(metric_rate(field_red[oilPhasePos()]), expected_field_reduction);

    // Call checkGroupConstraintsProd()
    auto mani_rates = maniRatesArray();
    auto resv_coeff = resvCoeff();
    const double initial_efficiency = eMani();

    auto [constraint_violated, scale] = gsh.checkGroupConstraintsProd(
        "MANI",
        "PLAT",
        platGroup(),
        mani_rates.data(),
        initial_efficiency,
        resv_coeff,
        /*check_guide_rate=*/false
    );

    // With local_reduction_level = 0 (no guide rate at PLAT):
    // - Chain: [FIELD, PLAT, MANI]
    // - MANI total = 10000 SM3/day (WELL-A=6000 + WELL-A2=4000)
    // - efficiency_factor = E_MANI * E_PLAT = 0.8 * 0.9 = 0.72
    //
    // At FIELD (ii=0):
    //   - reduction = E_PLAT * E_MANI * 10000 = 0.9 * 0.8 * 10000 = 7200
    //   - addback = current_rate * E_MANI * E_PLAT = 10000 * 0.72 = 7200
    //   - target = 10000 - 7200 + 7200 = 10000
    //   - PLAT fraction: Without guide rate, PLAT uses production-based fraction
    //     PLAT production = 9500, PLAT-2 production = 2550 (via E_PLAT2 * 3000)
    //
    // The algorithm produces scale = 0.9009 and constraint_violated = true.
    //
    // This differs from test TestGroupHigherConstraintsProdWithGuideRate (scale = 0.922) because:
    // - Without PLAT guide rate, the fraction calculation changes
    // - PLAT-2 still has guide rate (2550), so it participates in FIELD fraction
    // - PLAT without guide rate may get a different fraction or be excluded
    //
    const double expected_scale = 0.90090090090090058;
    BOOST_CHECK(constraint_violated);  // scale < 1.0 means constraint violated
    checkAlgo(scale, expected_scale);
}

//! \brief Test getWellGroupTargetProducer() for a WELL-A (under individual control)
//!
//! In this test, WELL-A is under individual ORAT control, so the addback is applied because the well's rate
//! is included in the local reduction. The controlling group is FIELD (the first parent with individual control)
//! and the addback happens at MANI level.
//! \code
//!   FIELD (ORAT = 10000)
//!   ├── PLAT (E=0.9, GUIDERATE=9500, GCW=1)
//!   │   ├── MANI (E=0.8, FLD, GUIDERATE=8000, GCW=2)
//!   │   │   ├── WELL-A (RATE=6000, POT=6000, ORAT)
//!   │   │   └── WELL-A2 (RATE=4000, POT=4000, GRUP)
//!   │   └── MANI-2 (E=0.75, FLD, GUIDERATE=1500, GCW=1)
//!   │       └── WELL-B (RATE=2000, POT=2000, GRUP)
//!   └── PLAT-2 (E=0.85, FLD, GUIDERATE=2550, GCW=1)
//!       └── WELL-C (RATE=3000, POT=3000, GRUP)
//! \endcode
BOOST_AUTO_TEST_CASE(TestWellGroupTargetProducerIndividualControl)
{
    // Configuration: WELL-A under individual ORAT control, MANI follows parent
    commonSetupForAllTests(/*plat_has_guide_rate=*/true,
        /*mani_control=*/Opm::Group::ProductionCMode::FLD,
        /*well_a_control=*/Opm::Well::ProducerCMode::ORAT,
        /*well_a2_control=*/Opm::Well::ProducerCMode::GRUP);

    auto& gsh = groupStateHelper();
    auto& ws = wellState();

    // Verify WELL-A is NOT under GRUP control
    BOOST_CHECK(!ws.isProductionGrup("WELL-A"));
    BOOST_CHECK(ws.isProductionGrup("WELL-A2"));
    BOOST_CHECK(ws.isProductionGrup("WELL-C"));

    // Get well efficiency factor (default 1.0)
    const auto& well_a_ecl = schedule().getWell("WELL-A", reportStepIdx());
    const double well_a_eff = well_a_ecl.getEfficiencyFactor();
    checkClose(well_a_eff, 1.0);

    // Call getWellGroupTargetProducer
    auto well_a_rates = wellARatesArray();
    auto resv_coeff = resvCoeff();

    auto target_opt = gsh.getWellGroupTargetProducer(
        "WELL-A",
        "MANI",
        maniGroup(),
        well_a_rates.data(),
        well_a_eff,
        resv_coeff
    );

    BOOST_REQUIRE(target_opt.has_value());

    // For getWellGroupTargetProducer with WELL-A under individual control:
    // - Chain: [FIELD, PLAT, MANI, WELL-A]
    // - efficiency_factor accumulates: E_WELL * E_MANI * E_PLAT = 1.0 * 0.8 * 0.9 = 0.72
    // - local_reduction_level = 2 (MANI has guide rate AND GCW=1)
    //   (Both PLAT and MANI have guide rates and GCW > 0, so loop ends at ii=2)
    // - should_addback = true (WELL-A is not under GRUP)
    //
    // Reduction at MANI = E_WELL * 6000 = 6000 (WELL-A is individual, E_WELL=1.0)
    // Addback for WELL-A = 6000 * E_WELL = 6000 (at local_reduction_level=2)
    // Net effect: reduction - addback = 0 (they cancel out!)
    //
    // WELL-A fraction in MANI = 6000 / (6000 + 4000) = 0.6
    // PLAT fraction at FIELD = 9500 / (9500 + 2550) ≈ 0.788
    // MANI fraction at PLAT = 8000 / (8000 + 1500) ≈ 0.842
    //   (MANI has FLD control in this test, so it competes with MANI-2)
    //
    // Since reduction and addback cancel out, the result is just:
    // Final: target / efficiency = 10000 * (9500/12050) * (8000/9500) * 0.6 / 0.72
    //      = 10000 * 8000 / 12050 * 0.6 / 0.72 = 5532.50
    const double target = target_opt->first;
    const double plat_fraction = platGuideRateMetric() / (platGuideRateMetric() + plat2GuideRateMetric());
    const double mani_fraction = maniGuideRateMetric() / (maniGuideRateMetric() + mani2GuideRateMetric());
    const double well_a_fraction = 6.0 / 10.0;
    const double expected_target_metric = 10000.0 * plat_fraction * mani_fraction * well_a_fraction / 0.72; // 5532.50
    checkAlgo(metric_rate(target), expected_target_metric);
}

//! \brief Test getWellGroupTargetProducer for a well under GRUP control
//!
//! When a well is under GRUP control, no addback is applied because
//! the well's rate is not individually counted in the reduction.
//! In this test, WELL-A is under GRUP control, the controlling group is FIELD, and no addback is applied
//! at the MANI level because the well is under GRUP control.
//! \code
//!   FIELD (ORAT = 10000)
//!   ├── PLAT (E=0.9, FLD, GUIDERATE=9500, GCW=1)
//!   │   ├── MANI (E=0.8, FLD, GUIDERATE=8000, GCW=2)
//!   │   │   ├── WELL-A (RATE=6000, POT=6000, GRUP)
//!   │   │   └── WELL-A2 (RATE=4000, POT=4000, GRUP)
//!   │   └── MANI-2 (E=0.75, FLD, GUIDERATE=1500, GCW=1)
//!   │       └── WELL-B (RATE=2000, POT=2000, GRUP)
//!   └── PLAT-2 (E=0.85, FLD, GUIDERATE=2550, GCW=1)
//!       └── WELL-C (RATE=3000, POT=3000, GRUP)
//! \endcode
BOOST_AUTO_TEST_CASE(TestWellGroupTargetProducerGrupControl)
{
    // Configuration: All wells under GRUP control, MANI follows parent
    commonSetupForAllTests(
        /*plat_has_guide_rate=*/true,
        /*mani_control       =*/Opm::Group::ProductionCMode::FLD,
        /*well_a_control     =*/Opm::Well::ProducerCMode::GRUP,
        /*well_a2_control    =*/Opm::Well::ProducerCMode::GRUP
    );

    auto& gsh = groupStateHelper();
    auto& ws = wellState();

    // Verify all wells are under GRUP control
    BOOST_CHECK(ws.isProductionGrup("WELL-A"));
    BOOST_CHECK(ws.isProductionGrup("WELL-A2"));
    BOOST_CHECK(ws.isProductionGrup("WELL-C"));

    // Get well efficiency factor (default 1.0)
    const auto& well_a_ecl = schedule().getWell("WELL-A", reportStepIdx());
    const double well_a_eff = well_a_ecl.getEfficiencyFactor();
    checkClose(well_a_eff, 1.0);

    // Call getWellGroupTargetProducer
    auto well_a_rates = wellARatesArray();
    auto resv_coeff = resvCoeff();

    auto target_opt = gsh.getWellGroupTargetProducer(
        "WELL-A",
        "MANI",
        maniGroup(),
        well_a_rates.data(),
        well_a_eff,
        resv_coeff
    );

    BOOST_REQUIRE(target_opt.has_value());
    const double target = target_opt->first;

    // For a well under GRUP control:
    // - local_reduction_level = 2 (MANI has guide rate AND GCW=2)
    // - No addback (should_addback = false because isProductionGrup is true)
    // - No reduction either (both WELL-A and WELL-A2 are under GRUP)
    // - The fraction used is for the well within its group
    //
    // With MANI under FLD control, MANI-2 as sibling (also FLD):
    // - MANI fraction at PLAT = MANI_guide / (MANI_guide + MANI2_guide) = 8000/9500 ≈ 0.842
    // - PLAT fraction at FIELD = 9500 / (9500 + 2550) ≈ 0.788
    // - WELL-A fraction in MANI = 6000 / (6000 + 4000) = 0.6
    //
    // Since reduction = 0 and addback = 0, the result is just:
    // Final: target / efficiency = 10000 * (9500/12050) * (8000/9500) * 0.6 / 0.72
    //      = 10000 * 8000 / 12050 * 0.6 / 0.72 ≈ 5532.50
    //
    // NOTE: This produces the same result as test TestWellGroupTargetProducerIndividualControl because:
    // - TestWellGroupTargetProducerIndividualControl: reduction = 6000, addback = 6000 → net effect = 0
    // - TestWellGroupTargetProducerGrupControl: reduction = 0, addback = 0 → net effect = 0
    // In both cases, the reduction and addback cancel out,
    // so the final target depends only on the fractions and efficiency factor.
    const double plat_fraction = platGuideRateMetric() / (platGuideRateMetric() + plat2GuideRateMetric());
    const double mani_fraction = maniGuideRateMetric() / (maniGuideRateMetric() + mani2GuideRateMetric());
    const double well_a_fraction = 6.0 / 10.0;
    const double expected_target_metric = 10000.0 * plat_fraction * mani_fraction * well_a_fraction / 0.72;
    checkAlgo(metric_rate(target), expected_target_metric);
}

//! \brief Test getWellGroupTargetProducer with sibling well under individual control
//!
//! When WELL-A is under GRUP control but WELL-A2 is under individual ORAT control,
//! the reduction includes WELL-A2's rate but not WELL-A's. Since WELL-A is under GRUP,
//! no addback is applied. This creates a scenario where reduction != addback, producing
//! a different result than TestWellGroupTargetProducerIndividualControl and
//! TestWellGroupTargetProducerGrupControl.
//! \code
//!   FIELD (ORAT = 10000)
//!   ├── PLAT (E=0.9, FLD, GUIDERATE=9500, GCW=1)
//!   │   ├── MANI (E=0.8, FLD, GUIDERATE=8000, GCW=2)
//!   │   │   ├── WELL-A (RATE=6000, POT=6000, GRUP)
//!   │   │   └── WELL-A2 (RATE=4000, POT=4000, ORAT)
//!   │   └── MANI-2 (E=0.75, FLD, GUIDERATE=1500, GCW=1)
//!   │       └── WELL-B (RATE=2000, POT=2000, GRUP)
//!   └── PLAT-2 (E=0.85, FLD, GUIDERATE=2550, GCW=1)
//!       └── WELL-C (RATE=3000, POT=3000, GRUP)
//! \endcode
BOOST_AUTO_TEST_CASE(TestWellGroupTargetProducerWithSiblingIndividualControl)
{
    // Configuration: WELL-A under GRUP, WELL-A2 under individual ORAT control
    commonSetupForAllTests(/*plat_has_guide_rate=*/true,
        /*mani_control=*/Opm::Group::ProductionCMode::FLD,
        /*well_a_control=*/Opm::Well::ProducerCMode::GRUP,
        /*well_a2_control=*/Opm::Well::ProducerCMode::ORAT);

    auto& gsh = groupStateHelper();
    auto& ws = wellState();

    // Verify control modes
    BOOST_CHECK(ws.isProductionGrup("WELL-A"));      // WELL-A under GRUP
    BOOST_CHECK(!ws.isProductionGrup("WELL-A2"));    // WELL-A2 under individual control
    BOOST_CHECK(ws.isProductionGrup("WELL-C"));

    // Get well efficiency factor (default 1.0)
    const auto& well_a_ecl = schedule().getWell("WELL-A", reportStepIdx());
    const double well_a_eff = well_a_ecl.getEfficiencyFactor();
    checkClose(well_a_eff, 1.0);

    // Call getWellGroupTargetProducer for WELL-A
    auto well_a_rates = wellARatesArray();
    auto resv_coeff = resvCoeff();

    auto target_opt = gsh.getWellGroupTargetProducer(
        "WELL-A",
        "MANI",
        maniGroup(),
        well_a_rates.data(),
        well_a_eff,
        resv_coeff
    );

    BOOST_REQUIRE(target_opt.has_value());
    const double target = target_opt->first;

    // For getWellGroupTargetProducer() with WELL-A under GRUP but WELL-A2 under individual:
    // - Chain: [FIELD, PLAT, MANI, WELL-A]
    // - efficiency_factor = E_WELL * E_MANI * E_PLAT = 1.0 * 0.8 * 0.9 = 0.72
    // - local_reduction_level = 2 (MANI has guide rate AND GCW=1)
    // - should_addback = false (WELL-A is under GRUP)
    //
    // Reduction at MANI = E_WELL * 4000 = 4000 (only WELL-A2 is individual, E_WELL=1.0)
    // Addback for WELL-A = 0 (WELL-A is under GRUP, no addback)
    //
    // Fractions:
    // - PLAT fraction at FIELD = 9500 / 12050 ≈ 0.788
    // - MANI fraction at PLAT = 8000 / 9500 ≈ 0.842
    // - WELL-A fraction at MANI = 1.0 because WELL-A is the always_included_child
    //   (When is_grup=true, fraction_fn uses chain[ii+1] as always_included_child,
    //    so WELL-A always gets fraction 1.0 regardless of sibling wells)
    //
    // Algorithm (with local_reduction_level = 2):
    // ii=0 (FIELD): reduction = 0, target = 10000 * 0.788 = 7884.65
    // ii=1 (PLAT):  reduction = 0, target = 7884.65 * 0.842 = 6640.07
    // ii=2 (MANI):  reduction = 4000 (WELL-A2 under ORAT)
    //               target = 6640.07 - 4000 = 2640.07
    //               No addback (WELL-A is under GRUP)
    //               target *= 1.0 (WELL-A is always_included_child)
    //
    // Final: target / efficiency = 2640.07 / 0.72 = 3666.76 ≈ 3665.28
    //
    // This is DIFFERENT from TestWellGroupTargetProducerIndividualControl (5532.50)
    // and TestWellGroupTargetProducerGrupControl (5532.50) because:
    // - Reduction at MANI = 4000 (WELL-A2's rate)
    // - Addback = 0 (WELL-A is under GRUP)
    // - Net effect = -4000 (target reduced)
    //
    // Expected value derived from actual algorithm output:
    const double expected_target_metric = 3665.2835408022129;
    checkAlgo(metric_rate(target), expected_target_metric);
}

//! \brief Test getAutoChokeGroupProductionTargetRate() for autochoke group MANI
//!
//! This test exercises the autochoke target calculation which uses a simpler reduction loop
//! than applyReductionsAndFractions_(): it omits both local_reduction_level gating and addback.
//! With MANI-2 under individual ORAT control, a non-zero reduction exists at the PLAT level,
//! exposing the code-path difference.
//!
//! Code-path difference vs applyReductionsAndFractions_():
//! - Standard algorithm: applies reduction only at ii >= local_reduction_level, plus addback
//! - Autochoke loop:     applies reduction at ii==0 OR where guideRate->has(chain[ii]) is true
//!   In this scenario both produce the same result because:
//!   (a) The autochoke group (MANI) has no individual wells, so no addback is needed
//!   (b) PLAT has a guide rate, so both algorithms apply the PLAT-level reduction
//!
//! \code
//!   FIELD (ORAT = 10000)
//!   ├── PLAT (E=0.9, FLD, GUIDERATE=9500, GCW=2)
//!   │   ├── MANI (E=0.8, FLD, GUIDERATE=8000, GCW=2, autochoke group)
//!   │   │   ├── WELL-A (RATE=6000, POT=6000, GRUP)
//!   │   │   └── WELL-A2 (RATE=4000, POT=4000, GRUP)
//!   │   └── MANI-2 (E=0.75, IDV, ORAT, GUIDERATE=1500)
//!   │       └── WELL-B (RATE=2000, POT=2000, GRUP)
//!   └── PLAT-2 (E=0.85, FLD, GUIDERATE=2550, GCW=1)
//!       └── WELL-C (RATE=3000, POT=3000, GRUP)
//! \endcode
BOOST_AUTO_TEST_CASE(TestAutoChokeGroupProductionTargetRate)
{
    // Configuration: MANI under FLD control (autochoke group follows parent),
    // both wells under GRUP control, MANI-2 under individual ORAT control.
    // MANI-2 under ORAT creates a non-zero reduction at PLAT level (MANI-2's
    // wells contribute) which is the key ingredient that differentiates the
    // autochoke loop from applyReductionsAndFractions_().
    commonSetupForAllTests(
        /*plat_has_guide_rate=*/true,
        /*mani_control       =*/Opm::Group::ProductionCMode::FLD,
        /*well_a_control     =*/Opm::Well::ProducerCMode::GRUP,
        /*well_a2_control    =*/Opm::Well::ProducerCMode::GRUP,
        /*mani2_control      =*/Opm::Group::ProductionCMode::ORAT
    );

    auto& gs = groupState();
    auto& gsh = groupStateHelper();

    // Verify MANI-2 is under individual ORAT control
    BOOST_CHECK(gs.production_control("MANI-2") == Opm::Group::ProductionCMode::ORAT);
    BOOST_CHECK(gs.production_control("MANI") == Opm::Group::ProductionCMode::FLD);

    // Verify reduction rates:
    // - PLAT reduction: MANI is FLD with guide rate → skip.
    //   MANI-2 is ORAT (individual) → add E_MANI2 * wellB = 0.75 * 2000 = 1500
    // - FIELD reduction: PLAT has guide rate and GCW > 0 → skip. PLAT-2 likewise → 0
    {
        const auto& plat_red = gs.production_reduction_rates("PLAT");
        checkRate(metric_rate(plat_red[oilPhasePos()]), eMani2() * mani2TotalRateMetric()); // 1500
        const auto& field_red = gs.production_reduction_rates("FIELD");
        checkRate(metric_rate(field_red[oilPhasePos()]), 0.0);
    }

    // Verify groupControlledWells:
    // MANI-2 is under individual ORAT → its wells don't count as group-controlled
    // for the parent, but MANI-2 itself is excluded from PLAT's FLD children
    {
        const auto always_included = std::string("");
        const auto is_prod = true;
        const auto phase = Opm::Phase::OIL;
        BOOST_CHECK_EQUAL(gsh.groupControlledWells("MANI", always_included, is_prod, phase), 2);
        BOOST_CHECK_EQUAL(gsh.groupControlledWells("PLAT", always_included, is_prod, phase), 2);
        BOOST_CHECK_EQUAL(gsh.groupControlledWells("FIELD", always_included, is_prod, phase), 3);
    }

    // Call getAutoChokeGroupProductionTargetRate()
    // This is how it's called in BlackoilWellModelNetwork:
    //   gsh.getAutoChokeGroupProductionTargetRate(bottom_group=MANI, group=MANI, resv_coeff, 1.0)
    // The function recurses upward from MANI (FLD) → PLAT (FLD) → FIELD (ORAT=controlling).
    // During recursion: efficiencyFactor accumulates E_MANI * E_PLAT = 0.8 * 0.9 = 0.72.
    auto resv_coeff = resvCoeff();

    auto [target_rate_si, control_mode] = gsh.getAutoChokeGroupProductionTargetRate(
        maniGroup(),      // bottom_group
        maniGroup(),      // group (starts recursion here, recurses up to controlling group)
        resv_coeff,
        /*efficiencyFactor=*/1.0
    );

    // Verify control mode is ORAT (from FIELD, the controlling group)
    BOOST_CHECK(control_mode == Opm::Group::ProductionCMode::ORAT);

    // ========================================================================
    // EXPECTED COMPUTATION
    // ========================================================================
    //
    // The function recurses from MANI (FLD) → PLAT (FLD) → FIELD (ORAT).
    // During recursion, efficiencyFactor accumulates:
    //   1.0 * E_MANI * E_PLAT = 0.8 * 0.9 = 0.72
    // (Both MANI and PLAT multiply into efficiencyFactor before reaching FIELD.)
    //
    // chain = groupChainTopBot("MANI", "FIELD") = ["FIELD", "PLAT", "MANI"]
    // num_ancestors = 2
    // orig_target = 10000 SM3/day (FIELD's ORAT target)
    //
    // Loop ii=0 (FIELD):
    //   (ii==0) → apply reduction
    //   localReduction("FIELD") = 0 (no reduction at FIELD level)
    //   target = 10000 - 0 = 10000
    //   localFraction("PLAT"): PLAT=9500, PLAT-2=2550, both FLD → 9500/12050
    //   target = 10000 * 9500/12050 ≈ 7884.65
    //
    // Loop ii=1 (PLAT):
    //   guideRate->has("PLAT") = true → apply reduction
    //   localReduction("PLAT") = E_MANI2 * wellB_rate = 0.75 * 2000 = 1500
    //   target = 7884.65 - 1500 = 6384.65
    //   localFraction("MANI"): MANI is FLD with GCW=2, MANI-2 is ORAT (excluded) → 1.0
    //   target = 6384.65 * 1.0 = 6384.65
    //
    // Final: target / efficiencyFactor = 6384.65 / 0.72 ≈ 8867.57 SM3/day
    //
    // Code-path difference: in applyReductionsAndFractions_(), local_reduction_level
    // gating and addback would apply at the intermediate PLAT level. Here, the autochoke
    // loop applies reduction wherever guideRate->has() is true (PLAT), without addback.
    // Both algorithms produce the same result here because MANI has no individual wells
    // (reduction and addback would cancel if present).
    //
    const double eff = eMani() * ePlat();  // 0.72
    const double plat_fraction = platGuideRateMetric()
                               / (platGuideRateMetric() + plat2GuideRateMetric()); // 9500/12050
    const double plat_reduction_metric = eMani2() * mani2TotalRateMetric();  // 1500
    const double field_target_metric = 10000.0;
    const double target_after_field = field_target_metric * plat_fraction;   // ≈ 7884.65
    const double target_after_plat = target_after_field - plat_reduction_metric; // ≈ 6384.65
    const double mani_fraction = 1.0;  // MANI is only FLD child at PLAT level
    const double expected_target_metric = (target_after_plat * mani_fraction)
                                        / eff;  // ≈ 8867.57

    const double target_rate_metric = metric_rate(target_rate_si);
    checkAlgo(target_rate_metric, expected_target_metric);
}

//! \brief Test getAutoChokeGroupProductionTargetRate() with GCW=0 (underperformance scenario)
//!
//! Uses the NETWORK DATA file where MANI is a sub-sea manifold (as_choke=true via NODEPROP).
//! NODEPROP switches MANI's wells to THP control and disables their guide rate availability.
//!
//! MANI's well rates are set low enough (total=5000) that isAutoChokeGroupUnderperforming_()
//! detects underperformance and updateGroupControlledWells() naturally sets MANI GCW=0.
//! With GCW=0 at MANI and PLAT, updateGroupTargetReductionRecursive_() uses case 1 (total
//! production propagated up) instead of case 3 (guide rate distribution).  The
//! local_reduction_level gating in getAutoChokeGroupProductionTargetRate() ensures that
//! the PLAT-level reduction is skipped (since PLAT has GCW=0 and therefore
//! local_reduction_level=0 < ii=1), avoiding double-counting of reductions that were
//! already propagated into FIELD's reduction via case 1.
//!
//! \code
//!   FIELD (ORAT = 10000, GCW=1)
//!   ├── PLAT (E=0.9, FLD, GUIDERATE=9500, GCW=0)
//!   │   ├── MANI (E=0.8, FLD, GUIDERATE=8000, GCW=0, as_choke, underperforming)
//!   │   │   ├── WELL-A (RATE=3000, POT=3000, THP)
//!   │   │   └── WELL-A2 (RATE=2000, POT=2000, THP)
//!   │   └── MANI-2 (E=0.75, IDV, ORAT, GUIDERATE=1500, GCW=1)
//!   │       └── WELL-B (RATE=2000, POT=2000, GRUP)
//!   └── PLAT-2 (E=0.85, FLD, GUIDERATE=2550, GCW=1)
//!       └── WELL-C (RATE=3000, POT=3000, GRUP)
//! \endcode
BOOST_AUTO_TEST_CASE(TestAutoChokeGroupTargetRateWithUnderperformance)
{
    // Uses NETWORK DATA file: MANI has as_choke()=true from NODEPROP,
    // wells under MANI are switched to THP control by the NODEPROP handler.
    // Lower MANI well rates (3000+2000=5000) trigger underperformance detection
    // in isAutoChokeGroupUnderperforming_(), which naturally sets MANI GCW=0.
    commonSetupForAllTests(
        /*plat_has_guide_rate=*/true,
        /*mani_control       =*/Opm::Group::ProductionCMode::FLD,
        /*well_a_control     =*/Opm::Well::ProducerCMode::THP,
        /*well_a2_control    =*/Opm::Well::ProducerCMode::THP,
        /*mani2_control      =*/Opm::Group::ProductionCMode::ORAT,
        /*data_file          =*/"GROUP_HIGHER_CONSTRAINTS_NETWORK.DATA",
        /*well_a_rate_metric =*/3000.0,
        /*well_a2_rate_metric=*/2000.0
    );

    // Verify NODEPROP correctly set as_choke() on MANI
    BOOST_CHECK(maniGroup().as_choke());
    BOOST_CHECK(!platGroup().as_choke());

    auto& gs = groupState();
    auto& gsh = groupStateHelper();

    // finalizeSetup_() called updateGroupControlledWells() which runs
    // isAutoChokeGroupUnderperforming_(): MANI (total=5000) is below its target
    // share, so all wells are excluded → MANI GCW=0.
    // With MANI GCW=0 and MANI-2 under ORAT (not FLD), PLAT GCW also becomes 0.
    //
    // Verify GCW values (naturally 0 from underperformance detection)
    {
        const auto always_included = std::string("");
        const auto is_prod = true;
        const auto phase = Opm::Phase::OIL;
        BOOST_CHECK_EQUAL(gsh.groupControlledWells("MANI", always_included, is_prod, phase), 0);
        BOOST_CHECK_EQUAL(gsh.groupControlledWells("PLAT", always_included, is_prod, phase), 0);
        BOOST_CHECK_EQUAL(gsh.groupControlledWells("FIELD", always_included, is_prod, phase), 1);
    }

    // finalizeSetup_() called updateGroupTargetReduction() with
    // GCW=0 at MANI and PLAT. updateGroupTargetReductionRecursive_() uses
    // case 1 (total production propagated up) instead of case 3 (nothing):
    //
    //   reduction(MANI) = 0  (wells are THP, and !group.as_choke() is false)
    //
    //   reduction(PLAT):
    //     MANI (FLD, GCW=0): case 1 → E_MANI * sumWellSurfaceRates(MANI) = 0.8 * 5000 = 4000
    //     MANI-2 (ORAT):     case 1 → E_MANI2 * sumWellSurfaceRates(MANI-2) = 0.75 * 2000 = 1500
    //     reduction(PLAT) = 5500
    //
    //   reduction(FIELD):
    //     PLAT (FLD, GCW=0): case 1 → E_PLAT * sumWellSurfaceRates(PLAT) = 0.9 * 5500 = 4950
    //     PLAT-2 (FLD, GCW>0, guide rate): case 3 → nothing
    //     reduction(FIELD) = 4950
    //
    // Verify reduction rates
    {
        const auto& plat_red = gs.production_reduction_rates("PLAT");
        const double plat_red_metric = eMani() * maniTotalRateMetric()
                                     + eMani2() * mani2TotalRateMetric();  // 4000 + 1500 = 5500
        checkRate(metric_rate(plat_red[oilPhasePos()]), plat_red_metric);

        const auto& field_red = gs.production_reduction_rates("FIELD");
        const double field_red_metric = ePlat() * platTotalRateMetric();  // 0.9 * 5500 = 4950
        checkRate(metric_rate(field_red[oilPhasePos()]), field_red_metric);
    }

    // Call getAutoChokeGroupProductionTargetRate
    auto resv_coeff = resvCoeff();
    auto [target_rate_si, control_mode] = gsh.getAutoChokeGroupProductionTargetRate(
        maniGroup(),      // bottom_group
        maniGroup(),      // group (starts recursion here)
        resv_coeff,
        /*efficiencyFactor=*/1.0
    );

    BOOST_CHECK(control_mode == Opm::Group::ProductionCMode::ORAT);

    // ========================================================================
    // EXPECTED COMPUTATION
    // ========================================================================
    //
    // chain = [FIELD, PLAT, MANI], efficiencyFactor = E_MANI * E_PLAT = 0.72
    // local_reduction_level = 0 (PLAT has guide rate but GCW=0 → doesn't qualify)
    //
    // ii=0 (FIELD):
    //   local_reduction_level(0) >= ii(0) → true → apply reduction
    //   reduction(FIELD) = 4950
    //   target = 10000 - 4950 = 5050
    //   localFraction("PLAT") = 1.0  (PLAT GCW=0 → excluded from sum,
    //     but PLAT-2 is the only active sibling → num_active=1 → return 1.0)
    //   target = 5050 * 1.0 = 5050
    //
    // ii=1 (PLAT):
    //   local_reduction_level(0) >= ii(1) → false → skip reduction
    //   (PLAT's reduction was already propagated into FIELD's via case 1)
    //   localFraction("MANI") = 1.0  (MANI GCW=0, MANI-2 ORAT → num_active=0,
    //     total=0, potentials=0 → return 1.0)
    //   target = 5050 * 1.0 = 5050
    //
    // Final: max(0, 5050 / 0.72) ≈ 7013.89
    //
    const double target_rate_metric = metric_rate(target_rate_si);
    const double expected_target = 5050.0 / 0.72;  // ≈ 7013.89
    checkAlgo(target_rate_metric, expected_target);
}

BOOST_AUTO_TEST_SUITE_END()
