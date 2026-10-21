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
 * \brief Tests for the reservoir coupling master's distribution of group targets
 *
 * A master group none of whose slave wells are under group control is excluded
 * from its parent's guide-rate distribution, and its siblings share the rest of
 * the parent's target.  The excluded group is sent the target it would get if
 * it returned to group control.  These tests check that for a real well model
 * and a real ReservoirCouplingMaster, with the data the slaves would have
 * reported set directly on the master; no slave process is started.
 *
 * Group hierarchy (RC_MASTER_GROUP_EXCLUSION.DATA):
 * \code
 *              FIELD (ORAT 1000, LRAT 2000, WATER RATE 1000, GAS RATE 1000)
 *          /       |       |       \
 *        M1       M2       M3       M4
 *   (GR 100)  (GR 200)  (GR 300)  (item 8 = NO, own limits 50)
 * \endcode
 * M1 and M2 belong to slave S1, M3 and M4 to slave S2.
 * RC_MASTER_GROUP_EXCLUSION_NESTED.DATA puts a plain group P between FIELD and
 * M1, and RC_MASTER_GROUP_EXCLUSION_SIBLINGS.DATA puts both M1 and M2 under P.
 */

#define BOOST_TEST_MODULE RescoupConstraints
#include "../SimulatorFixture.hpp"
#include "../ToleranceAndUnitFixture.hpp"

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/utils/readDeck.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/FractionCalculator.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/rescoup/RescoupConstraintsCalculator.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace Opm::Properties::TTag {
    struct TestRescoupConstraintsTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

namespace {

//! \brief Like Opm::initSimulator(), but reads the deck with low parsing strictness
//!
//! The reservoir coupling keywords (SLAVES, GRUPMAST) are still on the list of
//! keywords flow rejects by default, so a coupled run needs
//! --parsing-strictness=low, and so does this test's master deck.
template <class TypeTag>
std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>>
initMasterSimulator(const std::string& filename)
{
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    Opm::setupSimulatorParameters<TypeTag>(filename, "test_rescoup_constraints");
    auto& params = Opm::FlowGenericVanguard::modelParams_;
    Opm::readDeck(Opm::FlowGenericVanguard::comm(), filename,
                  params.eclState_, params.eclSchedule_, params.udqState_,
                  params.actionState_, params.wtestState_, params.eclSummaryConfig_,
                  /*python=*/nullptr, /*parsingStrictness=*/"low",
                  /*actionParsingStrictness=*/"normal", /*inputSkipMode=*/"100",
                  /*initFromRestart=*/false, /*checkDeck=*/false, /*keepKeywords=*/false,
                  /*outputInterval=*/{}, /*slaveMode=*/false);
    return std::make_unique<Simulator>();
}

const std::string flat_deck{"RC_MASTER_GROUP_EXCLUSION.DATA"};
const std::string nested_deck{"RC_MASTER_GROUP_EXCLUSION_NESTED.DATA"};
const std::string sibling_deck{"RC_MASTER_GROUP_EXCLUSION_SIBLINGS.DATA"};

//! What the slave reports for one of its groups, in SM3/day
struct SlaveGroupReport
{
    double oil_rate{0.0};
    int num_group_controlled_producers{1};
    double water_injection_rate{0.0};
    double gas_injection_rate{0.0};
    int num_group_controlled_water_injectors{1};
    int num_group_controlled_gas_injectors{1};
};

//! \brief A well model built from a master deck, with a reservoir coupling master
//!   attached whose slave data is set by the test.
struct RescoupConstraintsFixture : ToleranceAndUnitFixture
{
    using TypeTag = Opm::Properties::TTag::TestRescoupConstraintsTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using Master = Opm::ReservoirCouplingMaster<double>;
    using Calculator = Opm::RescoupConstraintsCalculator<double, IndexTraits>;
    using FractionCalculator = Opm::GroupStateHelpers::FractionCalculator<double, IndexTraits>;
    using RcPhase = Opm::ReservoirCoupling::Phase;

    //! \brief Build the well model and the master, and update the group data
    //! \param data_file The deck: flat_deck, nested_deck or sibling_deck
    //! \param groups_without_group_controlled_producers Master groups whose slave
    //!   reports no producers under group control
    //! \param groups_without_group_controlled_water_injectors Master groups whose slave
    //!   reports no water injectors under group control
    //! \param uncoupled_slaves Slaves that are not coupled: never activated, and
    //!   reporting zeros, as the master assumes for a slave that has not
    //!   activated yet or has ended.  At least one slave must be coupled.
    void commonSetupForAllTests(
        const std::string& data_file = flat_deck,
        const std::vector<std::string>& groups_without_group_controlled_producers = {},
        const std::vector<std::string>& groups_without_group_controlled_water_injectors = {},
        const std::set<std::string>& uncoupled_slaves = {})
    {
        for (const auto& group : groups_without_group_controlled_producers) {
            slave_group_reports_.at(group).num_group_controlled_producers = 0;
        }
        for (const auto& group : groups_without_group_controlled_water_injectors) {
            slave_group_reports_.at(group).num_group_controlled_water_injectors = 0;
        }
        uncoupled_slaves_ = uncoupled_slaves;
        simulator_ = initMasterSimulator<TypeTag>(data_file);
        setupMaster_();
        wellModel().setReservoirCouplingMaster(master_.get());

        simulator_->model().applyInitialSolution();
        simulator_->setEpisodeIndex(-1);
        simulator_->setEpisodeLength(0.0);
        simulator_->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
        simulator_->setTimeStepSize(Opm::unit::day);
        simulator_->problem().resetIterationForNewTimestep();

        wellModel().beginReportStep(report_step_idx_);
        updateGroupData_();
    }

    //! \brief Send the slave reports to the master again, as the slaves do at the
    //!   next sync step, and update the group data from them
    //!
    //! Change the reports with slaveGroupReport() first.
    void resendSlaveData()
    {
        setSlaveData_();
        updateGroupData_();
    }

    // ------------------------------------------------------------------------
    // Accessors
    // ------------------------------------------------------------------------

    WellModel& wellModel() { return simulator_->problem().wellModel(); }
    const Opm::Schedule& schedule() const { return simulator_->vanguard().schedule(); }
    auto& groupState() { return wellModel().groupState(); }
    auto& groupStateHelper() { return wellModel().groupStateHelper(); }
    Opm::GuideRate& guideRate() { return wellModel().guideRate(); }
    Master& master() { return *master_; }
    //! What the slave reports for the master group \p group
    SlaveGroupReport& slaveGroupReport(const std::string& group)
    {
        return slave_group_reports_.at(group);
    }
    int reportStepIdx() const { return report_step_idx_; }

    int phasePos(int canonical_phase_idx)
    {
        return wellModel().phaseUsage().canonicalToActivePhaseIdx(canonical_phase_idx);
    }

    //! The constraints the master would send, by master group name
    struct Result
    {
        std::map<std::string, Opm::ReservoirCoupling::ProductionGroupConstraints<double>> production;
        std::map<std::pair<std::string, RcPhase>, double> injection;  //!< targets in SM3/day
        double productionTarget(const std::string& group) const
        {
            return metric_rate(production.at(group).target);
        }
        double injectionTarget(const std::string& group, RcPhase phase) const
        {
            return injection.at({group, phase});
        }
        //! The limit for the rate type \p mode in SM3/day, or -1 if the group
        //! hierarchy has no limit for it
        double productionLimit(const std::string& group, Opm::Group::ProductionCMode mode) const
        {
            using Mode = Opm::Group::ProductionCMode;
            const auto& pc = production.at(group);
            double limit = -1.0;
            switch (mode) {
            case Mode::ORAT: limit = pc.oil_limit; break;
            case Mode::WRAT: limit = pc.water_limit; break;
            case Mode::GRAT: limit = pc.gas_limit; break;
            case Mode::LRAT: limit = pc.liquid_limit; break;
            case Mode::RESV: limit = pc.resv_limit; break;
            default: BOOST_FAIL("No limit for this control mode");
            }
            return limit < 0.0 ? limit : metric_rate(limit);
        }
    };

    Result calculateConstraints()
    {
        auto guard = groupStateHelper().pushLogger();
        Calculator calculator{wellModel().guideRateHandler(), groupStateHelper()};
        const auto constraints = calculator.calculateMasterGroupConstraints();
        Result result;
        for (std::size_t slave_idx = 0; slave_idx < master().numSlaves(); ++slave_idx) {
            const auto& groups = master().getMasterGroupNamesForSlave(slave_idx);
            for (const auto& pc : constraints.production_constraints[slave_idx]) {
                result.production.emplace(groups[pc.group_name_idx], pc);
            }
            for (const auto& it : constraints.injection_targets[slave_idx]) {
                result.injection.emplace(std::pair{groups[it.group_name_idx], it.phase},
                                         metric_rate(it.target));
            }
        }
        return result;
    }

    FractionCalculator productionFractionCalculator()
    {
        auto& gsh = groupStateHelper();
        return FractionCalculator{schedule(), gsh, gsh.summaryState(), reportStepIdx(),
                                  &guideRate(), Opm::GuideRateModel::Target::OIL,
                                  /*is_producer=*/true, /*injection_phase=*/Opm::Phase::OIL};
    }

    FractionCalculator injectionFractionCalculator(Opm::Phase phase)
    {
        auto& gsh = groupStateHelper();
        return FractionCalculator{schedule(), gsh, gsh.summaryState(), reportStepIdx(),
                                  &guideRate(), gsh.getInjectionGuideTargetMode(phase),
                                  /*is_producer=*/false, phase};
    }

private:
    //! The master groups of each slave, in the order the slave reports them
    const std::map<std::string, std::vector<std::string>> master_groups_of_slave_{
        {"S1", {"M1", "M2"}},
        {"S2", {"M3", "M4"}},
    };

    //! The slave reports, by master group
    std::map<std::string, SlaveGroupReport> slave_group_reports_{
        {"M1", {150.0, 1, 150.0, 150.0, 1, 1}},
        {"M2", {200.0, 1, 200.0, 200.0, 1, 1}},
        {"M3", {300.0, 1, 300.0, 300.0, 1, 1}},
        {"M4", {50.0, 1, 50.0, 50.0, 1, 1}},
    };

    std::set<std::string> uncoupled_slaves_;
    std::unique_ptr<Opm::Simulator<TypeTag>> simulator_;
    std::unique_ptr<Master> master_;
    int report_step_idx_{0};

    // Does for the master what spawning the slaves and receiving their first data
    // would have done.
    void setupMaster_()
    {
        master_ = std::make_unique<Master>(Opm::FlowGenericVanguard::comm(), schedule(),
                                           /*argc=*/0, /*argv=*/nullptr);
        for (const auto& [slave_name, groups] : master_groups_of_slave_) {
            master_->addSlaveName(slave_name);
            std::map<std::string, std::size_t> group_order;
            for (std::size_t i = 0; i < groups.size(); ++i) {
                master_->getMasterGroupToSlaveNameMap()[groups[i]] = slave_name;
                group_order[groups[i]] = i;
            }
            master_->getSlaveNameToMasterGroupsMap()[slave_name] = groups;
            master_->updateMasterGroupNameOrderMap(slave_name, group_order);
        }
        master_->rebuildSlaveIdxToMasterGroupsVector();
        master_->initTimeStepping();
        master_->initStartOfReportStep(report_step_idx_);

        for (const auto& [slave_name, groups] : master_groups_of_slave_) {
            if (uncoupled_slaves_.count(slave_name) == 0) {
                master_->setSlaveActivated(master_->getSlaveIdx(slave_name));
            }
        }
        setSlaveData_();
    }

    // Sets the data the slaves report.  A slave that is not coupled reports
    // nothing, and the master then assumes zero data for it
    // (see ReservoirCouplingMasterReportStep::receiveProductionDataFromSlaves()).
    void setSlaveData_()
    {
        for (const auto& [slave_name, groups] : master_groups_of_slave_) {
            std::vector<Master::SlaveGroupProductionData> production_data(groups.size());
            std::vector<Master::SlaveGroupInjectionData> injection_data(groups.size());
            if (uncoupled_slaves_.count(slave_name) > 0) {
                master_->setSlaveGroupProductionData(slave_name, std::move(production_data));
                master_->setSlaveGroupInjectionData(slave_name, std::move(injection_data));
                continue;
            }
            for (std::size_t i = 0; i < groups.size(); ++i) {
                const auto& report = slave_group_reports_.at(groups[i]);
                auto& prod = production_data[i];
                prod.potentials[RcPhase::Oil] = si_rate(report.oil_rate);
                prod.surface_rates[RcPhase::Oil] = si_rate(report.oil_rate);
                prod.network_surface_rates[RcPhase::Oil] = si_rate(report.oil_rate);
                prod.reservoir_rates[RcPhase::Oil] = si_rate(report.oil_rate);
                prod.num_group_controlled_producers = report.num_group_controlled_producers;
                auto& inj = injection_data[i];
                inj.surface_rates[RcPhase::Water] = si_rate(report.water_injection_rate);
                inj.surface_rates[RcPhase::Gas] = si_rate(report.gas_injection_rate);
                inj.reservoir_rates[RcPhase::Water] = si_rate(report.water_injection_rate);
                inj.reservoir_rates[RcPhase::Gas] = si_rate(report.gas_injection_rate);
                inj.num_group_controlled_injectors[static_cast<std::size_t>(RcPhase::Water)]
                    = report.num_group_controlled_water_injectors;
                inj.num_group_controlled_injectors[static_cast<std::size_t>(RcPhase::Gas)]
                    = report.num_group_controlled_gas_injectors;
            }
            master_->setSlaveGroupProductionData(slave_name, std::move(production_data));
            master_->setSlaveGroupInjectionData(slave_name, std::move(injection_data));
        }
    }

    // As at the start of a time step: the group rates, with the master groups'
    // rates taken from the slave reports, then the guide rates.
    void updateGroupData_()
    {
        wellModel().updateAndCommunicateGroupData(report_step_idx_,
                                                  /*update_wellgrouptarget=*/false);
        wellModel().updateGuideRates(report_step_idx_, /*sim_time=*/0.0);
    }
};

} // anonymous namespace

BOOST_FIXTURE_TEST_SUITE(RescoupConstraintsTests, RescoupConstraintsFixture)

//! \brief The fixture gives the master groups the reported rates and the deck's guide rates
BOOST_AUTO_TEST_CASE(FixtureSetup)
{
    commonSetupForAllTests();

    auto& gsh = groupStateHelper();
    for (const std::string name : {"M1", "M2", "M3", "M4"}) {
        BOOST_CHECK(gsh.rescoup().isMasterGroup(name));
    }
    BOOST_CHECK(master().slaveIsCoupled(0));
    BOOST_CHECK(master().slaveIsCoupled(1));
    BOOST_CHECK_EQUAL(master().getSlaveGroupNumGroupControlledProducers("M1"), 1);

    const int oil = phasePos(IndexTraits::oilPhaseIdx);
    checkRate(metric_rate(groupState().production_rates("M1")[oil]), 150.0);
    checkRate(metric_rate(groupState().production_rates("M4")[oil]), 50.0);
    checkRate(metric_rate(groupState().production_rates("FIELD")[oil]), 700.0);

    const Opm::GuideRate::RateVector zero_rates{0.0, 0.0, 0.0};
    checkClose(guideRate().get("M1", Opm::GuideRateModel::Target::OIL, zero_rates), 100.0);
    checkClose(guideRate().get("M2", Opm::GuideRateModel::Target::OIL, zero_rates), 200.0);
    checkClose(guideRate().get("M3", Opm::GuideRateModel::Target::OIL, zero_rates), 300.0);
}

//! \brief With every master group under group control, M1-M3 share FIELD's
//!   target less M4's rate by their guide rates
//!
//! M4 is not available for higher-level control, so its rate of 50 is a reduction
//! on FIELD.  M1, M2, M3 share 1000 - 50 = 950 by 100 : 200 : 300.
//!
//! FIELD's LRAT limit of 2000 is shared the same way: the slaves report no water
//! production, so the liquid rates equal the oil rates, and so do the guide rates
//! for liquid.  M4 has no LRAT limit of its own, and takes no share of FIELD's.
BOOST_AUTO_TEST_CASE(NoExclusion)
{
    using Mode = Opm::Group::ProductionCMode;
    commonSetupForAllTests();
    const auto result = calculateConstraints();

    checkRate(result.productionTarget("M1"), 950.0 * 100.0 / 600.0);
    checkRate(result.productionTarget("M2"), 950.0 * 200.0 / 600.0);
    checkRate(result.productionTarget("M3"), 950.0 * 300.0 / 600.0);
    checkRate(result.productionTarget("M4"), 50.0);
    for (const std::string name : {"M1", "M2", "M3", "M4"}) {
        BOOST_CHECK_EQUAL(master().effectiveGCW(name), 1);
    }

    // FIELD's LRAT limit less M4's liquid rate, 2000 - 50 = 1950, is shared by the
    // guide rates of GCONPROD item 9.  These are oil guide rates (item 10), which
    // are scaled to liquid by each group's liquid / oil rate ratio, 1 here.
    checkRate(result.productionLimit("M1", Mode::LRAT), 1950.0 * 100.0 / 600.0);
    checkRate(result.productionLimit("M2", Mode::LRAT), 1950.0 * 200.0 / 600.0);
    checkRate(result.productionLimit("M3", Mode::LRAT), 1950.0 * 300.0 / 600.0);
    // M4 has no LRAT limit of its own, and takes no share of FIELD's (GCONPROD
    // item 8 = NO)
    BOOST_CHECK_EQUAL(result.productionLimit("M4", Mode::LRAT), -1.0);
    // The limit for the active control mode is the target, and the hierarchy has
    // no limits for the other rate types
    checkRate(result.productionLimit("M2", Mode::ORAT), result.productionTarget("M2"));
    for (const auto mode : {Mode::WRAT, Mode::GRAT, Mode::RESV}) {
        BOOST_CHECK_EQUAL(result.productionLimit("M2", mode), -1.0);
    }
}

//! \brief M1 has no producers under group control: it is excluded, and M2, M3
//!   share the rest of FIELD's target
//!
//! The reduction on FIELD is the sum of the rates of the groups that are not in
//! the guide-rate distribution: M4's 50 (item 8 = NO) plus the excluded M1's 150.
//! The master assumes that M1 keeps producing its reported 150, so M2 and M3 only
//! need to produce the rest, 1000 - 200 = 800, which they share by 200 : 300.
//! M1 is sent its share if it rejoined: (1000 - 50) * 100 / (100 + 200 + 300).
//!
//! The LRAT limits are recalculated with the targets: M2 and M3 share
//! 2000 - 200 = 1800, and M1 is sent (2000 - 50) * 100 / 600.  Without the
//! exclusion, M2 and M3 would have 650 (= 1950 * 200 / 600) and 975
//! (= 1950 * 300 / 600), see NoExclusion.
BOOST_AUTO_TEST_CASE(ExcludeAndRedistributeProduction)
{
    using Mode = Opm::Group::ProductionCMode;
    commonSetupForAllTests(flat_deck, /*groups_without_group_controlled_producers=*/{"M1"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveGCW("M1"), 0);
    BOOST_CHECK_EQUAL(master().effectiveGCW("M2"), 1);
    BOOST_CHECK_EQUAL(master().effectiveGCW("M3"), 1);

    const int oil = phasePos(IndexTraits::oilPhaseIdx);
    checkRate(metric_rate(groupState().production_reduction_rates("FIELD")[oil]), 200.0);

    checkRate(result.productionTarget("M2"), 800.0 * 200.0 / 500.0);
    checkRate(result.productionTarget("M3"), 800.0 * 300.0 / 500.0);
    checkRate(result.productionTarget("M1"), 950.0 * 100.0 / 600.0);
    checkRate(result.productionTarget("M4"), 50.0);

    checkRate(result.productionLimit("M2", Mode::LRAT), 1800.0 * 200.0 / 500.0);
    checkRate(result.productionLimit("M3", Mode::LRAT), 1800.0 * 300.0 / 500.0);
    checkRate(result.productionLimit("M1", Mode::LRAT), 1950.0 * 100.0 / 600.0);
    BOOST_CHECK_EQUAL(result.productionLimit("M4", Mode::LRAT), -1.0);
}

//! \brief The fraction fallback: an excluded master group counts for its own
//!   fraction, and not for its siblings'
BOOST_AUTO_TEST_CASE(FractionFallback)
{
    commonSetupForAllTests(flat_deck, /*groups_without_group_controlled_producers=*/{"M1"});
    calculateConstraints();
    BOOST_REQUIRE_EQUAL(master().effectiveGCW("M1"), 0);

    auto fcalc = productionFractionCalculator();
    // M1's own target: M1 is the always-included child
    checkClose(fcalc.localFraction("M1", "M1"), 100.0 / 600.0);
    checkClose(fcalc.fraction("M1", "FIELD", /*always_include_this=*/true), 100.0 / 600.0);
    // M2's target, with M2 as the always-included child: M1 is not in the
    // guide-rate sum.  The fallback only counts the always-included child itself
    // and its ancestors, and M1 is neither, so M2 does not bring M1 back.
    checkClose(fcalc.localFraction("M2", "M2"), 200.0 / 500.0);
    // Without an always-included child the fraction is the same
    checkClose(fcalc.localFraction("M2", ""), 200.0 / 500.0);
    checkClose(fcalc.fraction("M3", "FIELD", /*always_include_this=*/true), 300.0 / 500.0);
}

//! \brief The exclusion step skips a master group not available for higher-level
//!   control (item 8 = NO)
//!
//! M4 never takes part in the guide-rate distribution: its rate is a reduction on
//! FIELD in any case.  When it reports no producers under group control, the
//! exclusion step leaves it alone: its effective GCW is not set to 0, the targets
//! are not recalculated, and the targets of the others do not change.
BOOST_AUTO_TEST_CASE(ItemEightNoIsSkippedByExclusion)
{
    commonSetupForAllTests(flat_deck, /*groups_without_group_controlled_producers=*/{"M4"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveGCW("M4"), 1);
    checkRate(result.productionTarget("M1"), 950.0 * 100.0 / 600.0);
    checkRate(result.productionTarget("M2"), 950.0 * 200.0 / 600.0);
    checkRate(result.productionTarget("M3"), 950.0 * 300.0 / 600.0);
    checkRate(result.productionTarget("M4"), 50.0);
}

//! \brief Injection is excluded per phase: M1 has no water injectors under group
//!   control, but has gas injectors under group control
//!
//! Water: the reduction on FIELD is M4's 50 plus M1's 150, so M2 and M3 share
//! 1000 - 200 = 800 by 200 : 300, and M1 is sent (1000 - 50) * 100 / 600.
//! Gas: nothing is excluded, so M1, M2, M3 share 950 by 100 : 200 : 300.
BOOST_AUTO_TEST_CASE(ExcludeAndRedistributeInjectionPerPhase)
{
    commonSetupForAllTests(flat_deck,
                           /*groups_without_group_controlled_producers=*/{},
                           /*groups_without_group_controlled_water_injectors=*/{"M1"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveInjectionGCW("M1", Opm::Phase::WATER), 0);
    BOOST_CHECK_EQUAL(master().effectiveInjectionGCW("M1", Opm::Phase::GAS), 1);

    checkRate(result.injectionTarget("M2", RcPhase::Water), 800.0 * 200.0 / 500.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 800.0 * 300.0 / 500.0);
    checkRate(result.injectionTarget("M1", RcPhase::Water), 950.0 * 100.0 / 600.0);

    checkRate(result.injectionTarget("M1", RcPhase::Gas), 950.0 * 100.0 / 600.0);
    checkRate(result.injectionTarget("M2", RcPhase::Gas), 950.0 * 200.0 / 600.0);
    checkRate(result.injectionTarget("M3", RcPhase::Gas), 950.0 * 300.0 / 600.0);

    // Production is untouched
    checkRate(result.productionTarget("M2"), 950.0 * 200.0 / 600.0);
}

//! \brief The nested case: M1 is the only child of a plain group P
//!
//! RC_MASTER_GROUP_EXCLUSION_NESTED.DATA puts P between FIELD and M1.  When M1 is
//! excluded P has no group-controlled wells either, and only the fallback in
//! FractionCalculator::groupControlledWells() keeps M1 and P in the fractions for
//! M1's own target.  The targets must be the same
//! as for the flat tree.
BOOST_AUTO_TEST_CASE(NestedExclusion)
{
    commonSetupForAllTests(nested_deck, /*groups_without_group_controlled_producers=*/{"M1"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveGCW("M1"), 0);
    checkRate(result.productionTarget("M2"), 800.0 * 200.0 / 500.0);
    checkRate(result.productionTarget("M3"), 800.0 * 300.0 / 500.0);
    checkRate(result.productionTarget("M1"), 950.0 * 100.0 / 600.0);

    auto fcalc = productionFractionCalculator();
    checkClose(fcalc.fraction("M1", "FIELD", /*always_include_this=*/true), 100.0 / 600.0);
    checkClose(fcalc.localFraction("P", "M1"), 100.0 / 600.0);
    // Without M1 as the always-included child, P has no group-controlled wells
    // and drops out of FIELD's guide-rate sum
    checkClose(fcalc.localFraction("M2", ""), 200.0 / 500.0);
}

//! \brief The fraction fallback (FractionCalculator::groupControlledWells()) for
//!   injection: an excluded master group counts for its own fraction of the water
//!   target, and not for its siblings'
//!
//! The injection counterpart of FractionFallback.  M1 has no water injectors
//! under group control, but has gas injectors under group control.
BOOST_AUTO_TEST_CASE(InjectionFractionFallback)
{
    // Make M1 have no water injectors under group control
    commonSetupForAllTests(flat_deck,
                           /*groups_without_group_controlled_producers=*/{},
                           /*groups_without_group_controlled_water_injectors=*/{"M1"});
    calculateConstraints();
    BOOST_REQUIRE_EQUAL(master().effectiveInjectionGCW("M1", Opm::Phase::WATER), 0);

    auto water = injectionFractionCalculator(Opm::Phase::WATER);
    // M1's own target: M1 is the always-included child
    checkClose(water.localFraction("M1", "M1"), 100.0 / 600.0);
    checkClose(water.fraction("M1", "FIELD", /*always_include_this=*/true), 100.0 / 600.0);
    // Without M1 as the always-included child, M1 has no group-controlled
    // injectors, and so no guide rate
    checkClose(water.localFraction("M1", ""), 0.0);
    // M2's target: FIELD's guide-rate sum is M2's 200 plus M3's 300.  M1 is
    // excluded, and the fallback only counts the always-included child and its
    // ancestors, so M2 does not bring M1 back; M4 is not available for
    // higher-level control.  So M2's fraction is 200 / 500, with M2 as the
    // always-included child or without one.
    checkClose(water.localFraction("M2", "M2"), 200.0 / 500.0);
    checkClose(water.localFraction("M2", ""), 200.0 / 500.0);

    // Gas: nothing is excluded
    auto gas = injectionFractionCalculator(Opm::Phase::GAS);
    checkClose(gas.localFraction("M1", "M1"), 100.0 / 600.0);
    checkClose(gas.localFraction("M2", ""), 200.0 / 600.0);
}

//! \brief The exclusion step skips a master group not available for higher-level
//!   injection control (item 8 = NO)
//!
//! The injection counterpart of ItemEightNoIsSkippedByExclusion.  M4 reports no
//! water injectors under group control; its effective injection GCW is left at
//! 1, and the water targets are those without exclusion.
BOOST_AUTO_TEST_CASE(ItemEightNoIsSkippedByInjectionExclusion)
{
    // Make M4 have no water injectors under group control
    commonSetupForAllTests(flat_deck,
                           /*groups_without_group_controlled_producers=*/{},
                           /*groups_without_group_controlled_water_injectors=*/{"M4"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveInjectionGCW("M4", Opm::Phase::WATER), 1);
    checkRate(result.injectionTarget("M1", RcPhase::Water), 950.0 * 100.0 / 600.0);
    checkRate(result.injectionTarget("M2", RcPhase::Water), 950.0 * 200.0 / 600.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 950.0 * 300.0 / 600.0);
    checkRate(result.injectionTarget("M4", RcPhase::Water), 50.0);
}

//! \brief The nested case for injection: M1, the only child of P, has no water
//!   injectors under group control
//!
//! M1's slave injectors inject 150 on their own limits (see
//! RescoupConstraintsFixture::slave_group_reports_).
//! The injection counterpart of NestedExclusion.  The water targets are the same
//! as for the flat tree, and production is untouched.
BOOST_AUTO_TEST_CASE(NestedInjectionExclusion)
{
    // Make M1 have no water injectors under group control
    commonSetupForAllTests(nested_deck,
                           /*groups_without_group_controlled_producers=*/{},
                           /*groups_without_group_controlled_water_injectors=*/{"M1"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveInjectionGCW("M1", Opm::Phase::WATER), 0);
    checkRate(result.injectionTarget("M2", RcPhase::Water), 800.0 * 200.0 / 500.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 800.0 * 300.0 / 500.0);
    checkRate(result.injectionTarget("M1", RcPhase::Water), 950.0 * 100.0 / 600.0);
    checkRate(result.productionTarget("M2"), 950.0 * 200.0 / 600.0);

    auto water = injectionFractionCalculator(Opm::Phase::WATER);
    // M1's own target: P's guide rate is M1's 100, counted by the fallback, of
    // FIELD's sum 100 + 200 + 300
    checkClose(water.localFraction("P", "M1"), 100.0 / 600.0);
    // M1 is P's only child, so it gets all of P's target
    checkClose(water.localFraction("M1", "M1"), 1.0);
    checkClose(water.localFraction("M2", ""), 200.0 / 500.0);
}

//! \brief The nested case with siblings: M1 and M2 are the children of a plain
//!   group P, and M1 has no producers and no water injectors under group control
//!
//! M1's slave wells produce and inject 150 on their own limits (see
//! RescoupConstraintsFixture::slave_group_reports_).  P has no guide rate of its own, so
//! M1's rate is a reduction on P that passes up to FIELD, and P's guide rate is
//! the sum of its children's.  The targets are the same as for the flat tree.
BOOST_AUTO_TEST_CASE(NestedExclusionWithSibling)
{
    // Make M1 have no producers and no water injectors under group control
    commonSetupForAllTests(sibling_deck,
                           /*groups_without_group_controlled_producers=*/{"M1"},
                           /*groups_without_group_controlled_water_injectors=*/{"M1"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveGCW("M1"), 0);
    const int oil = phasePos(IndexTraits::oilPhaseIdx);
    checkRate(metric_rate(groupState().production_reduction_rates("FIELD")[oil]), 200.0);
    checkRate(result.productionTarget("M2"), 800.0 * 200.0 / 500.0);
    checkRate(result.productionTarget("M3"), 800.0 * 300.0 / 500.0);
    checkRate(result.productionTarget("M1"), 950.0 * 100.0 / 600.0);
    checkRate(result.injectionTarget("M2", RcPhase::Water), 800.0 * 200.0 / 500.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 800.0 * 300.0 / 500.0);
    checkRate(result.injectionTarget("M1", RcPhase::Water), 950.0 * 100.0 / 600.0);

    auto production = productionFractionCalculator();
    auto water = injectionFractionCalculator(Opm::Phase::WATER);
    for (auto* fcalc : {&production, &water}) {
        // M1's own target.  P has no guide rate of its own, so its guide rate is
        // the sum of its children's (FractionCalculator::guideRate()), with M1
        // counted as the always-included child: 100 + 200 = 300.  FIELD's sum is
        // P's 300 plus M3's 300 = 600; M4 is not available for higher-level
        // control.
        checkClose(fcalc->localFraction("P", "M1"), 300.0 / 600.0);
        checkClose(fcalc->localFraction("M1", "M1"), 100.0 / 300.0);
        // M2's target: P counts only M2, so FIELD's sum is 200 + 300, and M2 is
        // the only child left in P's sum
        checkClose(fcalc->localFraction("P", ""), 200.0 / 500.0);
        checkClose(fcalc->localFraction("M2", ""), 1.0);
    }
}

//! \brief Two master groups excluded at once: M1 and M2 have no producers and no
//!   water injectors under group control
//!
//! The reduction on FIELD is 50 (from M4) + 150 (from M1) + 200 (from M2) = 400,
//! the rates M1 and M2 produce on their own limits (see
//! RescoupConstraintsFixture::slave_group_reports_) included, and M3, the only group left in
//! the distribution, is sent the rest, 600.  Each excluded group is sent its share
//! if it alone rejoined: M1 shares 1000 - 400 + 150 with M3 by 100 : 300, M2
//! shares 1000 - 400 + 200 with M3 by 200 : 300.
BOOST_AUTO_TEST_CASE(TwoGroupsExcluded)
{
    // Make M1 and M2 have no producers and no water injectors under group control
    commonSetupForAllTests(flat_deck,
                           /*groups_without_group_controlled_producers=*/{"M1", "M2"},
                           /*groups_without_group_controlled_water_injectors=*/{"M1", "M2"});
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveGCW("M1"), 0);
    BOOST_CHECK_EQUAL(master().effectiveGCW("M2"), 0);
    const int oil = phasePos(IndexTraits::oilPhaseIdx);
    checkRate(metric_rate(groupState().production_reduction_rates("FIELD")[oil]), 400.0);

    checkRate(result.productionTarget("M3"), 600.0);
    checkRate(result.productionTarget("M1"), 750.0 * 100.0 / 400.0);
    checkRate(result.productionTarget("M2"), 800.0 * 200.0 / 500.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 600.0);
    checkRate(result.injectionTarget("M1", RcPhase::Water), 750.0 * 100.0 / 400.0);
    checkRate(result.injectionTarget("M2", RcPhase::Water), 800.0 * 200.0 / 500.0);
}

//! \brief Every master group available for higher-level control is excluded:
//!   M1, M2 and M3 have no producers and no water injectors under group control
//!
//! No group is left in FIELD's guide-rate distribution, and the reduction on FIELD
//! is the rate of all four groups, 700.  Each excluded group is sent its share if
//! it alone rejoined, as in TwoGroupsExcluded: with no other group in the sum, that
//! is all of FIELD's target the others leave, 1000 - 700 + its own rate.  The
//! three targets add up to more than FIELD's target; each holds only if its group
//! is the only one that returns to group control.
BOOST_AUTO_TEST_CASE(AllEligibleGroupsExcluded)
{
    using Mode = Opm::Group::ProductionCMode;
    // Make M1, M2 and M3 have no producers and no water injectors under group
    // control
    commonSetupForAllTests(flat_deck,
                           /*groups_without_group_controlled_producers=*/{"M1", "M2", "M3"},
                           /*groups_without_group_controlled_water_injectors=*/{"M1", "M2", "M3"});
    const auto result = calculateConstraints();

    for (const std::string name : {"M1", "M2", "M3"}) {
        BOOST_CHECK_EQUAL(master().effectiveGCW(name), 0);
        BOOST_CHECK_EQUAL(master().effectiveInjectionGCW(name, Opm::Phase::WATER), 0);
    }
    BOOST_CHECK_EQUAL(groupState().number_of_wells_under_group_control("FIELD"), 0);
    const int oil = phasePos(IndexTraits::oilPhaseIdx);
    checkRate(metric_rate(groupState().production_reduction_rates("FIELD")[oil]), 700.0);

    // FIELD's ORAT target less the reduction, 1000 - 700 = 300, plus the group's
    // own rate, which is in the reduction and is added back for its own target.
    // The group is alone in FIELD's guide-rate sum for its own target, so its
    // fraction is 1.
    checkRate(result.productionTarget("M1"), 300.0 + 150.0);
    checkRate(result.productionTarget("M2"), 300.0 + 200.0);
    checkRate(result.productionTarget("M3"), 300.0 + 300.0);
    checkRate(result.productionTarget("M4"), 50.0);
    // The same for FIELD's LRAT limit: 2000 - 700 = 1300, plus M1's own liquid
    // rate, which equals its oil rate
    checkRate(result.productionLimit("M1", Mode::LRAT), 1300.0 + 150.0);
    // The same for FIELD's water injection target: 1000 less the water the four
    // groups inject, 700, plus the group's own
    checkRate(result.injectionTarget("M1", RcPhase::Water), 300.0 + 150.0);
    checkRate(result.injectionTarget("M2", RcPhase::Water), 300.0 + 200.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 300.0 + 300.0);
}

//! \brief A slave that is not coupled: its master groups get no constraints, and
//!   take no share of FIELD's targets
//!
//! S2 (M3 and M4) has not activated, and reports zeros.  Both its groups are
//! excluded from the distributions of production and of both injection phases,
//! including M3, which is available for higher-level control.  M1 and M2 share
//! all of FIELD's targets by 100 : 200.
BOOST_AUTO_TEST_CASE(UncoupledSlave)
{
    using Mode = Opm::Group::ProductionCMode;
    // Make S2 not coupled
    commonSetupForAllTests(flat_deck, {}, {}, /*uncoupled_slaves=*/{"S2"});
    BOOST_REQUIRE(!master().slaveIsCoupled(master().getSlaveIdx("S2")));
    const auto result = calculateConstraints();

    for (const std::string name : {"M3", "M4"}) {
        BOOST_CHECK_EQUAL(result.production.count(name), 0);
        BOOST_CHECK_EQUAL(result.injection.count({name, RcPhase::Water}), 0);
        BOOST_CHECK_EQUAL(result.injection.count({name, RcPhase::Gas}), 0);
        BOOST_CHECK_EQUAL(master().effectiveGCW(name), 0);
        BOOST_CHECK_EQUAL(master().effectiveInjectionGCW(name, Opm::Phase::WATER), 0);
        BOOST_CHECK_EQUAL(master().effectiveInjectionGCW(name, Opm::Phase::GAS), 0);
    }
    const int oil = phasePos(IndexTraits::oilPhaseIdx);
    checkRate(metric_rate(groupState().production_reduction_rates("FIELD")[oil]), 0.0);

    checkRate(result.productionTarget("M1"), 1000.0 * 100.0 / 300.0);
    checkRate(result.productionTarget("M2"), 1000.0 * 200.0 / 300.0);
    checkRate(result.productionLimit("M1", Mode::LRAT), 2000.0 * 100.0 / 300.0);
    checkRate(result.productionLimit("M2", Mode::LRAT), 2000.0 * 200.0 / 300.0);
    for (const auto phase : {RcPhase::Water, RcPhase::Gas}) {
        checkRate(result.injectionTarget("M1", phase), 1000.0 * 100.0 / 300.0);
        checkRate(result.injectionTarget("M2", phase), 1000.0 * 200.0 / 300.0);
    }
}

//! \brief A second calculation after an excluded group returns to group control
//!
//! M1 has no producers and no water injectors under group control in the first
//! calculation.  Then its slave reports group-controlled producers and water
//! injectors again.  The second calculation must not keep the exclusion of the
//! first one: it gives the targets of NoExclusion.
BOOST_AUTO_TEST_CASE(RecalculateAfterGroupReturns)
{
    // Make M1 have no producers and no water injectors under group control
    commonSetupForAllTests(flat_deck,
                           /*groups_without_group_controlled_producers=*/{"M1"},
                           /*groups_without_group_controlled_water_injectors=*/{"M1"});
    auto result = calculateConstraints();
    BOOST_REQUIRE_EQUAL(master().effectiveGCW("M1"), 0);
    BOOST_REQUIRE_EQUAL(master().effectiveInjectionGCW("M1", Opm::Phase::WATER), 0);
    checkRate(result.productionTarget("M2"), 800.0 * 200.0 / 500.0);

    slaveGroupReport("M1").num_group_controlled_producers = 1;
    slaveGroupReport("M1").num_group_controlled_water_injectors = 1;
    resendSlaveData();
    result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveGCW("M1"), 1);
    BOOST_CHECK_EQUAL(master().effectiveInjectionGCW("M1", Opm::Phase::WATER), 1);
    const int oil = phasePos(IndexTraits::oilPhaseIdx);
    checkRate(metric_rate(groupState().production_reduction_rates("FIELD")[oil]), 50.0);
    checkRate(result.productionTarget("M1"), 950.0 * 100.0 / 600.0);
    checkRate(result.productionTarget("M2"), 950.0 * 200.0 / 600.0);
    checkRate(result.productionTarget("M3"), 950.0 * 300.0 / 600.0);
    checkRate(result.injectionTarget("M1", RcPhase::Water), 950.0 * 100.0 / 600.0);
    checkRate(result.injectionTarget("M2", RcPhase::Water), 950.0 * 200.0 / 600.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 950.0 * 300.0 / 600.0);
}

//! \brief A second calculation with the same exclusion gives the same targets
//!
//! M1 has no producers and no water injectors under group control in both
//! calculations.  The second calculation excludes M1 again and redistributes
//! again, so M2 and M3 are sent their redistributed targets, not those before the
//! exclusion (316.67 and 475).
BOOST_AUTO_TEST_CASE(RecalculateWithSameExclusion)
{
    // Make M1 have no producers and no water injectors under group control
    commonSetupForAllTests(flat_deck,
                           /*groups_without_group_controlled_producers=*/{"M1"},
                           /*groups_without_group_controlled_water_injectors=*/{"M1"});
    calculateConstraints();
    resendSlaveData();
    const auto result = calculateConstraints();

    BOOST_CHECK_EQUAL(master().effectiveGCW("M1"), 0);
    BOOST_CHECK_EQUAL(master().effectiveInjectionGCW("M1", Opm::Phase::WATER), 0);
    checkRate(result.productionTarget("M2"), 800.0 * 200.0 / 500.0);
    checkRate(result.productionTarget("M3"), 800.0 * 300.0 / 500.0);
    checkRate(result.productionTarget("M1"), 950.0 * 100.0 / 600.0);
    checkRate(result.injectionTarget("M2", RcPhase::Water), 800.0 * 200.0 / 500.0);
    checkRate(result.injectionTarget("M3", RcPhase::Water), 800.0 * 300.0 / 500.0);
    checkRate(result.injectionTarget("M1", RcPhase::Water), 950.0 * 100.0 / 600.0);
}

BOOST_AUTO_TEST_SUITE_END()
