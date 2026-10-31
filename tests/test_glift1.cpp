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

#define BOOST_TEST_MODULE Glift1

#include "SimulatorFixture.hpp"

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/equil/EquilibrationHelpers.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/GasLiftStage2.hpp>
#include <opm/simulators/wells/BlackoilWellModelGasLift.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <optional>
#include <string>
#include <vector>

namespace Opm::Properties::TTag {
    struct TestGliftTypeTag {
        using InheritsFrom = std::tuple<TestTypeTag>;
    };
}

using SimulatorFixture = Opm::SimulatorFixture;
BOOST_GLOBAL_FIXTURE(SimulatorFixture);

BOOST_AUTO_TEST_CASE(G1)
{
    //using TypeTag = Opm::Properties::TTag::FlowProblemBlackoil;
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    //using EclProblem = Opm::EclProblem<TypeTag>;
    //using EclWellModel = typename EclProblem::EclWellModel;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using WellState = Opm::WellState<double, IndexTraits>;
    using StdWell = Opm::StandardWell<TypeTag>;
    using GasLiftSingleWell = Opm::GasLiftSingleWell<TypeTag>;
    using GasLiftGroupInfo = Opm::GasLiftGroupInfo<double, IndexTraits>;
    using GasLiftSingleWellGeneric = Opm::GasLiftSingleWellGeneric<double, IndexTraits>;
    using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
    const std::string filename = "GLIFT1.DATA";
    using GLiftSyncGroups = typename GasLiftSingleWellGeneric::GLiftSyncGroups;

    auto simulator = Opm::initSimulator<TypeTag>(filename.data(), "test_glift1", /*threads=*/2);

    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);  // 12 hours
    // Reset iteration context so code using problem().iterationContext() sees first iteration
    simulator->problem().resetIterationForNewTimestep();
    WellModel& well_model = simulator->problem().wellModel();

    // we tests 4 different setups given in the .DATA file
    // we only look at B-1H
    // 1) No rate limit -> well should be limited by alq
    // 2) No alq is needed. ORAT limit -> oil_rate = oil_limit and gas_rate is scaled
    // 3) Alq is needed. ORAT limit -> oil_rate = oil_limit and gas_rate is scaled
    // 4) Alq is needed. GRAT limit -> gas_rate = gas_limit and oil_rate is scaled
    for (int report_step_idx = 0; report_step_idx < 4; report_step_idx++ ) {
        well_model.beginReportStep(report_step_idx);
        well_model.beginTimeStep();
        auto logger_guard = well_model.groupStateHelper().pushLogger();
        auto& deferred_logger = well_model.groupStateHelper().deferredLogger();
        well_model.calculateExplicitQuantities();
        well_model.prepareTimeStep(deferred_logger);
        well_model.updateWellControls(deferred_logger);
        const Opm::WellInterface<TypeTag>* well_ptr = &well_model.getWell("B-1H");
        const StdWell *std_well = dynamic_cast<const StdWell *>(well_ptr);
        const auto& schedule = simulator->vanguard().schedule();
        auto wells_ecl = schedule.getWells(report_step_idx);
        std::optional<std::size_t> idx;
        int num_producers = 0;
        for(std::size_t i = 0; i < wells_ecl.size(); ++i) {
            const auto &well =  wells_ecl[i];
            if (well.isProducer() && well.name()=="B-1H") {
                idx = i;
                num_producers++;
            }
        }
        BOOST_CHECK_EQUAL( num_producers, 1);
        const auto &well = wells_ecl[*idx];
        BOOST_CHECK_EQUAL( well.name(), "B-1H");
        const auto& summary_state = simulator->vanguard().summaryState();
        WellState &well_state = well_model.wellState();
        const auto &group_state = well_model.groupState();
        GLiftEclWells ecl_well_map;
        Opm::BlackoilWellModelGasLift<TypeTag>::
            initGliftEclWellMap(well_model.localNonshutWells(), ecl_well_map);
        const auto& iterCtx = simulator->problem().iterationContext();
        const auto& comm = simulator->vanguard().grid().comm();
        GasLiftGroupInfo group_info {
            ecl_well_map,
            schedule,
            summary_state,
            simulator->episodeIndex(),
            iterCtx,
            deferred_logger,
            well_state,
            group_state,
            comm,
            /*glift_debug=*/false
        };
        GLiftSyncGroups sync_groups;
        auto std_well_unconst = *std_well;
        GasLiftSingleWell glift {std_well_unconst, *(simulator.get()), summary_state,
            deferred_logger, well_state, group_state, group_info, sync_groups,
            comm, /*glift_debug=*/false
        };
        group_info.initialize();
        const auto& controls = well.productionControls(summary_state);
        auto state = glift.runOptimize(iterCtx.iteration());
        auto pot = well_state[well.name()].well_potentials;

        const auto& glo = schedule.glo(report_step_idx);
        const auto increment = glo.gaslift_increment();
        const auto gl_well = glo.well(well.name());

        // no well limit -> alq limited
        if (report_step_idx == 0) {
            BOOST_CHECK(state->alqIsLimited());
            BOOST_CHECK(state->increase().has_value());
            BOOST_CHECK(!state->gasIsLimited());
            BOOST_CHECK(!state->oilIsLimited());
            BOOST_CHECK_CLOSE(state->alq(), *gl_well.max_rate(), 1e-8);
            BOOST_CHECK_CLOSE(state->oilRate(), pot[1], 1e-8);
            BOOST_CHECK_CLOSE(state->gasRate(), pot[2], 1e-8);
        }
        // ORAT limit no alq needed
        if (report_step_idx == 1) {
            BOOST_CHECK(!state->alqIsLimited());
            BOOST_CHECK(!state->gasIsLimited());
            BOOST_CHECK(state->oilIsLimited());
            // the oil rate is constraint by the oil target
            BOOST_CHECK_CLOSE(state->oilRate(), controls.oil_rate, 1e-8);
            // the gas rate is scaled by the oil target / oil potential
            BOOST_CHECK_CLOSE(state->gasRate(), pot[2] * controls.oil_rate / pot[1], 1e-6);
            BOOST_CHECK(!state->increase().has_value());
            BOOST_CHECK_CLOSE(state->alq(), 0.0, 1e-8);
        }
        // ORAT limit, alq needed
        if (report_step_idx == 2) {
            BOOST_CHECK(!state->alqIsLimited());
            BOOST_CHECK(!state->gasIsLimited());
            BOOST_CHECK(state->oilIsLimited());
            // the oil rate is constraint by the oil target
            BOOST_CHECK_CLOSE(state->oilRate(), controls.oil_rate, 1e-8);
            // the gas rate is scaled by the oil target / oil potential
            BOOST_CHECK_CLOSE(state->gasRate(), pot[2] * controls.oil_rate / pot[1], 1e-6);
            BOOST_CHECK(state->increase().has_value());
            // for this setup we needed one alq increment
            BOOST_CHECK_CLOSE(state->alq(), increment*1, 1e-8);
        }
        // GRAT limit, alq needed
        if (report_step_idx == 3) {
            BOOST_CHECK(!state->alqIsLimited());
            BOOST_CHECK(state->gasIsLimited());
            BOOST_CHECK(!state->oilIsLimited());
            // the gas rate is constraint by the gas target
            BOOST_CHECK_CLOSE(state->gasRate(), controls.gas_rate, 1e-8);
            // the oil rate is scaled by the gas target / gas potential
            BOOST_CHECK_CLOSE(state->oilRate(), pot[1] * controls.gas_rate / pot[2], 1e-6);
            BOOST_CHECK(state->increase().has_value());
            // for this setup we needed five alq increment
            BOOST_CHECK_CLOSE(state->alq(), increment*5, 1e-8);
        }
    }
}


BOOST_AUTO_TEST_CASE(G2)
{
    //using TypeTag = Opm::Properties::TTag::FlowProblemBlackoil;
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    //using EclProblem = Opm::EclProblem<TypeTag>;
    //using EclWellModel = typename EclProblem::EclWellModel;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using WellState = Opm::WellState<double, IndexTraits>;
    using GasLiftSingleWell = Opm::GasLiftSingleWell<TypeTag>;
    using GasLiftGroupInfo = Opm::GasLiftGroupInfo<double, IndexTraits>;
    using GasLiftSingleWellGeneric = Opm::GasLiftSingleWellGeneric<double, IndexTraits>;
    using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
    using GLiftOptWells = typename Opm::BlackoilWellModelGasLift<TypeTag>::GLiftOptWells;
    using GLiftProdWells = typename Opm::BlackoilWellModelGasLift<TypeTag>::GLiftProdWells;
    using GLiftWellStateMap = typename Opm::BlackoilWellModelGasLift<TypeTag>::GLiftWellStateMap;
    using GasLiftStage2 = Opm::GasLiftStage2<double, IndexTraits>;


    // we use the same file but start from report idx 4
    // where more wells and group controll are added
    const std::string filename = "GLIFT1.DATA";
    using GLiftSyncGroups = typename GasLiftSingleWellGeneric::GLiftSyncGroups;

    auto simulator = Opm::initSimulator<TypeTag>(filename.data(), "test_glift1", /*threads=*/2);

    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);  // 12 hours
    // Reset iteration context so code using problem().iterationContext() sees first iteration
    simulator->problem().resetIterationForNewTimestep();
    WellModel& well_model = simulator->problem().wellModel();

    const auto& comm = simulator->vanguard().grid().comm();
    const auto& summary_state = simulator->vanguard().summaryState();

    // we tests 4 different setups for stage 2 optimization
    // starting from report step 5 as given in the .DATA file
    // 1) No group rate limit -> group should be limited by max alq
    // 2) Group is limited by ORAT -> sufficient ALQ to production the target
    // 3) Group is limited by ORAT for PLAT-1 and GRAT for M5S
    // 4) No group limits but restriction on total gas
    for (int report_step_idx = 5; report_step_idx < 9; report_step_idx++ ) {
        GLiftOptWells glift_wells;
        GLiftProdWells prod_wells;
        GLiftWellStateMap state_map;
        simulator->setEpisodeIndex(report_step_idx);
        well_model.beginReportStep(report_step_idx);
        well_model.beginTimeStep();
        auto logger_guard = well_model.groupStateHelper().pushLogger();
        auto& deferred_logger = well_model.groupStateHelper().deferredLogger();
        well_model.calculateExplicitQuantities();
        const auto& schedule = simulator->vanguard().schedule();
        auto wells_ecl = schedule.getWells(report_step_idx);
        GLiftEclWells ecl_well_map;
        Opm::BlackoilWellModelGasLift<TypeTag>::
            initGliftEclWellMap(well_model.localNonshutWells(), ecl_well_map);
        const auto& iterCtx = simulator->problem().iterationContext();
        WellState& well_state = well_model.wellState();
        const auto& group_state = well_model.groupState();
        GasLiftGroupInfo group_info {
            ecl_well_map,
            schedule,
            summary_state,
            report_step_idx,
            iterCtx,
            deferred_logger,
            well_state,
            group_state,
            comm,
            /*glift_debug=*/false
        };
        group_info.initialize();
        GLiftSyncGroups sync_groups;
        for (const auto& well : well_model.localNonshutWells()) {
            if (group_info.hasWell(well->name())) {
                auto glift = std::make_unique<GasLiftSingleWell> (*well, *(simulator.get()), summary_state,
                    deferred_logger, well_state, group_state, group_info, sync_groups,
                    comm, /*glift_debug=*/false
                );
                auto state = glift->runOptimize(iterCtx.iteration());
                if (state) {
                    state_map.emplace(well->name(), std::move(state));
                    glift_wells.emplace(well->name(), std::move(glift));
                }
                prod_wells.insert({well->name(), well.get()});
            }
        }
        GasLiftStage2 glift2 {report_step_idx,
                            comm,
                            schedule,
                            summary_state,
                            deferred_logger,
                            well_state,
                            group_state,
                            prod_wells,
                            glift_wells,
                            group_info,
                            state_map,
                            /*glift_debug=*/false
        };
        glift2.runOptimize();

        // no group rate limits -> alq is limited by GLIFTOPT item 1
        if (report_step_idx == 5) {
            BOOST_CHECK_CLOSE(group_info.alqRate("PLAT-A")*86400, 450000, 1e-8);
        }
        // group oil rate limit -> give sufficient alq to produce group target
        if (report_step_idx == 6) {
            BOOST_CHECK_CLOSE(group_info.alqRate("PLAT-A")*86400, 150000, 1e-8);
            BOOST_CHECK(group_info.oilRate("PLAT-A") > *group_info.oilTarget("PLAT-A"));

            // also test wells
            std::vector<std::pair<std::string, double>> wells = {{"B-1H",25000}, {"B-2H",37500}, {"B-3H",25000},
                                                                 {"C-1H",37500}, {"C-2H",25000}};
            for (auto well: wells) {
                BOOST_CHECK_CLOSE(state_map[well.first]->alq()*86400, well.second, 1e-8);
            }
        }

        // same as above but with GRAT limit on sub-group (i.e need some more alq)
        if (report_step_idx == 7) {
            BOOST_CHECK_CLOSE(group_info.alqRate("PLAT-A")*86400, 187500, 1e-8);
            BOOST_CHECK(group_info.oilRate("PLAT-A") > *group_info.oilTarget("PLAT-A"));
        }

        // max alq + gas limit
        if (report_step_idx == 8) {
            BOOST_CHECK_CLOSE( (group_info.alqRate("PLAT-A") + group_info.gasRate("PLAT-A"))*86400, 800000, 1);
        }
        well_model.endTimeStep();
        well_model.endEpisode();
    }
}

namespace {

// A gas lift well with a synthetic inflow and VFP response. The well has no
// solution at the THP limit for ALQ values below min_flow_alq, as a well that
// needs lift gas to flow against its THP limit, and none above max_flow_alq.
// In between the lift gas lowers the bhp.
template<class IndexTraits>
class SyntheticGasLiftWell : public Opm::GasLiftSingleWellGeneric<double, IndexTraits>
{
    using Base = Opm::GasLiftSingleWellGeneric<double, IndexTraits>;

public:
    SyntheticGasLiftWell(const Opm::WellInterfaceGeneric<double, IndexTraits>& well,
                         Opm::DeferredLogger& deferred_logger,
                         Opm::WellState<double, IndexTraits>& well_state,
                         const Opm::GroupState<double>& group_state,
                         const Opm::SummaryState& summary_state,
                         Opm::GasLiftGroupInfo<double, IndexTraits>& group_info,
                         const Opm::Schedule& schedule,
                         const int report_step_idx,
                         typename Base::GLiftSyncGroups& sync_groups,
                         const Opm::Parallel::Communication& comm,
                         const double min_flow_alq,
                         const double max_flow_alq)
        : Base(deferred_logger, well_state, group_state, well.wellEcl(), summary_state,
               group_info, schedule, report_step_idx, sync_groups, comm, /*glift_debug=*/false)
        , well_(well)
        , min_flow_alq_(min_flow_alq)
        , max_flow_alq_(max_flow_alq)
    {
        this->optimize_ = true;
        this->min_alq_ = 0.0;
        this->max_alq_ = 20.0 * this->increment_;
        this->orig_alq_ = 0.0;
        this->alpha_w_ = 1.0;
        this->alpha_g_ = 0.0;
        this->max_iterations_ = 1000;
        this->oil_pos_ = IndexTraits::oilPhaseIdx;
        this->gas_pos_ = IndexTraits::gasPhaseIdx;
        this->water_pos_ = IndexTraits::waterPhaseIdx;
    }

    double increment() const { return this->increment_; }
    double maxAlq() const { return this->max_alq_; }
    void setMaxAlq(const double alq) { this->max_alq_ = alq; }
    void setOrigAlq(const double alq) { this->orig_alq_ = alq; }
    void limitAlq() { this->limitAlqToMax_(); }
    double initialOilRate() { return this->getInitialRatesWithLimit_().first->oil; }

    const Opm::WellInterfaceGeneric<double, IndexTraits>& getWell() const override
    {
        return well_;
    }

private:
    static constexpr double p_res = 200.0e5;

    bool checkThpControl_() const override { return true; }

    std::optional<double> computeBhpAtThpLimit_(double alq, double, bool) const override
    {
        if (alq < min_flow_alq_ - 1e-8 || alq > max_flow_alq_ + 1e-8) {
            return std::nullopt; // no solution at the THP limit
        }
        return p_res - 10.0e5 - 1.0e5 * (alq - min_flow_alq_) / this->increment_;
    }

    typename Base::RatesAndBhp computeWellRates_(double bhp, bool bhp_is_limited, bool) const override
    {
        const double oil = 1.0e-8 * (p_res - bhp);
        return {oil, 10.0 * oil, 0.0, bhp, bhp_is_limited};
    }

    const Opm::WellInterfaceGeneric<double, IndexTraits>& well_;
    double min_flow_alq_;
    double max_flow_alq_;
};

} // anonymous namespace

BOOST_AUTO_TEST_CASE(WellTestALQ)
{
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using GasLiftGroupInfo = Opm::GasLiftGroupInfo<double, IndexTraits>;
    using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
    using Synthetic = SyntheticGasLiftWell<IndexTraits>;

    auto simulator = Opm::initSimulator<TypeTag>("GLIFT1.DATA", "test_glift1_welltest", /*threads=*/2);
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);
    simulator->problem().resetIterationForNewTimestep();
    WellModel& well_model = simulator->problem().wellModel();

    well_model.beginReportStep(0);
    well_model.beginTimeStep();
    auto logger_guard = well_model.groupStateHelper().pushLogger();
    auto& deferred_logger = well_model.groupStateHelper().deferredLogger();
    well_model.calculateExplicitQuantities();
    well_model.prepareTimeStep(deferred_logger);
    well_model.updateWellControls(deferred_logger);

    const auto& well = well_model.getWell("B-1H");
    const auto& schedule = simulator->vanguard().schedule();
    const auto& summary_state = simulator->vanguard().summaryState();
    auto& well_state = well_model.wellState();
    const auto& group_state = well_model.groupState();
    GLiftEclWells ecl_well_map;
    Opm::BlackoilWellModelGasLift<TypeTag>::
        initGliftEclWellMap(well_model.localNonshutWells(), ecl_well_map);
    const auto& iterCtx = simulator->problem().iterationContext();
    const auto& comm = simulator->vanguard().grid().comm();
    GasLiftGroupInfo group_info {ecl_well_map, schedule, summary_state,
                                 simulator->episodeIndex(), iterCtx, deferred_logger,
                                 well_state, group_state, comm, /*glift_debug=*/false};
    group_info.initialize();

    const auto increment = schedule.glo(0).gaslift_increment();
    typename Synthetic::GLiftSyncGroups sync_groups;
    const double no_limit = 1.0e6;
    auto makeWell = [&](const double min_flow_alq, const double max_flow_alq) {
        return Synthetic {well, deferred_logger, well_state, group_state, summary_state,
                          group_info, schedule, simulator->episodeIndex(), sync_groups,
                          comm, min_flow_alq, max_flow_alq};
    };

    // The well flows at ALQ = 0: no ALQ increase is needed to start the search.
    {
        auto glift = makeWell(0.0, no_limit);
        const auto [alq, success] = glift.wellTestALQ();
        BOOST_CHECK(success);
        BOOST_CHECK_GE(alq, 0.0);
    }

    // The well has no solution at the THP limit until it gets three increments of
    // lift gas. The test must raise the ALQ until the well flows, not give up at 0.
    {
        const double min_flow_alq = 3.0 * increment;
        auto glift = makeWell(min_flow_alq, no_limit);
        const auto [alq, success] = glift.wellTestALQ();
        BOOST_CHECK(success);
        BOOST_CHECK_GE(alq, min_flow_alq - 1e-8);
        BOOST_CHECK_LE(alq, glift.maxAlq() + 1e-8);
    }

    // As above, but the well can only flow at that one ALQ, so a further increase
    // is not possible. The ALQ that made the well flow is still a valid test ALQ.
    {
        const double min_flow_alq = 3.0 * increment;
        auto glift = makeWell(min_flow_alq, min_flow_alq);
        const auto [alq, success] = glift.wellTestALQ();
        BOOST_CHECK(success);
        BOOST_CHECK_CLOSE(alq, min_flow_alq, 1e-6);
    }

    // The maximum ALQ is not reached by whole increments and the well only flows at
    // the maximum: the last increment is limited to the maximum, as in the optimization.
    {
        const double max_alq = 2.5 * increment;
        auto glift = makeWell(max_alq, no_limit);
        glift.setMaxAlq(max_alq);
        const auto [alq, success] = glift.wellTestALQ();
        BOOST_CHECK(success);
        BOOST_CHECK_CLOSE(alq, max_alq, 1e-6);
    }

    // The well cannot flow with any ALQ up to the maximum: the test fails.
    {
        auto glift = makeWell(100.0 * increment, no_limit);
        const auto [alq, success] = glift.wellTestALQ();
        BOOST_CHECK(!success);
    }
}

// The ALQ set by gas lift optimization is kept to the next report step, unless the
// well gets a new target from the schedule.
BOOST_AUTO_TEST_CASE(AlqCarryOver)
{
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;

    auto simulator = Opm::initSimulator<TypeTag>("GLIFT2.DATA", "test_glift1_alq", /*threads=*/2);
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);
    simulator->problem().resetIterationForNewTimestep();
    WellModel& well_model = simulator->problem().wellModel();

    const double day = 86400.0;
    const double alq = 50000.0 / day;
    auto alqOfWell = [&well_model]() {
        return well_model.wellState().well("B-1H").alq_state.get();
    };

    // Report step 0: set an ALQ as the optimization would do
    simulator->setEpisodeIndex(0);
    well_model.beginReportStep(0);
    BOOST_CHECK_SMALL(alqOfWell(), 1e-12);
    well_model.wellState().well("B-1H").alq_state.set(alq);
    well_model.commitWGState();

    // Report step 1: no event for the well, the ALQ is kept
    simulator->setEpisodeIndex(1);
    well_model.beginReportStep(1);
    BOOST_CHECK_CLOSE(alqOfWell(), alq, 1e-8);
    well_model.commitWGState();

    // Report step 2: WLIFTOPT lowers the maximum ALQ, which is no well event. The ALQ is
    // kept here and limited to the maximum when the optimization starts.
    simulator->setEpisodeIndex(2);
    well_model.beginReportStep(2);
    BOOST_CHECK_CLOSE(alqOfWell(), alq, 1e-8);
    well_model.commitWGState();

    // Report step 3: new target for the well, the ALQ is reset to the default
    simulator->setEpisodeIndex(3);
    well_model.beginReportStep(3);
    BOOST_CHECK_SMALL(alqOfWell(), 1e-12);
}

// An ALQ above the maximum ALQ is limited before the optimization starts, also in the group rates.
BOOST_AUTO_TEST_CASE(AlqLimitedToMax)
{
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using GasLiftGroupInfo = Opm::GasLiftGroupInfo<double, IndexTraits>;
    using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
    using Synthetic = SyntheticGasLiftWell<IndexTraits>;

    auto simulator = Opm::initSimulator<TypeTag>("GLIFT1.DATA", "test_glift1_maxalq", /*threads=*/2);
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);
    simulator->problem().resetIterationForNewTimestep();
    WellModel& well_model = simulator->problem().wellModel();

    // report step 5 has the group PLAT-A with a limit on the lift gas
    simulator->setEpisodeIndex(5);
    well_model.beginReportStep(5);
    well_model.beginTimeStep();
    auto logger_guard = well_model.groupStateHelper().pushLogger();
    auto& deferred_logger = well_model.groupStateHelper().deferredLogger();
    well_model.calculateExplicitQuantities();

    const auto& well = well_model.getWell("B-1H");
    const auto& schedule = simulator->vanguard().schedule();
    const auto& summary_state = simulator->vanguard().summaryState();
    auto& well_state = well_model.wellState();
    const auto& group_state = well_model.groupState();
    const auto& comm = simulator->vanguard().grid().comm();
    const double increment = schedule.glo(5).gaslift_increment();
    const double max_alq = 20.0 * increment; // the maximum ALQ of the synthetic well
    const double high_alq = 25.0 * increment;
    typename Synthetic::GLiftSyncGroups sync_groups;

    auto test = [&](const double start_alq, const double expected_alq, const bool init_groups = true) {
        well_state.well("B-1H").alq_state.set(start_alq);
        GLiftEclWells ecl_well_map;
        Opm::BlackoilWellModelGasLift<TypeTag>::
            initGliftEclWellMap(well_model.localNonshutWells(), ecl_well_map);
        const auto& iterCtx = simulator->problem().iterationContext();
        GasLiftGroupInfo group_info {ecl_well_map, schedule, summary_state, 5, iterCtx,
                                     deferred_logger, well_state, group_state, comm,
                                     /*glift_debug=*/false};
        // the group info is not initialized when only used for a well test
        if (init_groups) {
            group_info.initialize();
        }
        const double group_alq = init_groups ? group_info.alqRate("PLAT-A") : 0.0;
        Synthetic glift {well, deferred_logger, well_state, group_state, summary_state,
                         group_info, schedule, 5, sync_groups, comm, 0.0, 1.0e6};
        glift.setOrigAlq(start_alq);
        glift.limitAlq();
        BOOST_CHECK_CLOSE(well_state.well("B-1H").alq_state.get(), expected_alq, 1e-8);
        if (init_groups) {
            // the group ALQ follows the change of the well
            BOOST_CHECK_SMALL(group_info.alqRate("PLAT-A") - group_alq - (expected_alq - start_alq), 1e-10);
        }
    };

    // above the maximum: limited to the maximum
    test(high_alq, max_alq);
    // not above the maximum: unchanged
    test(10.0 * increment, 10.0 * increment);
    test(max_alq, max_alq);
    // The rates of the optimization are computed with the limited ALQ, the group oil rate
    // and the well potentials are moved from the rates of the well state to those rates.
    {
        const auto oil_pos = well_state.phaseUsageInfo().canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        well_state.well("B-1H").alq_state.set(high_alq);
        const double old_oil = well_state.well("B-1H").well_potentials[oil_pos];
        GLiftEclWells ecl_well_map;
        Opm::BlackoilWellModelGasLift<TypeTag>::
            initGliftEclWellMap(well_model.localNonshutWells(), ecl_well_map);
        const auto& iterCtx = simulator->problem().iterationContext();
        GasLiftGroupInfo group_info {ecl_well_map, schedule, summary_state, 5, iterCtx,
                                     deferred_logger, well_state, group_state, comm,
                                     /*glift_debug=*/false};
        group_info.initialize();
        const double group_oil = group_info.oilRate("PLAT-A");
        Synthetic glift {well, deferred_logger, well_state, group_state, summary_state,
                         group_info, schedule, 5, sync_groups, comm, 0.0, 1.0e6};
        glift.setOrigAlq(high_alq);
        glift.limitAlq();
        const double new_oil = glift.initialOilRate();
        BOOST_CHECK(std::abs(new_oil - old_oil) > 1e-12); // the test needs the rates to differ
        BOOST_CHECK_SMALL(group_info.oilRate("PLAT-A") - group_oil - (new_oil - old_oil), 1e-10);
        BOOST_CHECK_SMALL(well_state.well("B-1H").well_potentials[oil_pos] - new_oil, 1e-10);
        // doing it again changes nothing
        glift.initialOilRate();
        BOOST_CHECK_SMALL(group_info.oilRate("PLAT-A") - group_oil - (new_oil - old_oil), 1e-10);
    }
    // uninitialized group info (well test): the ALQ is limited without a group lookup
    test(high_alq, max_alq, /*init_groups=*/false);
}

// The group rates of a well that is limited by its own target are those of the limited rates.
// Reducing the ALQ to the maximum moves the group rates from the limited rates with the old ALQ
// to the limited rates with the new, and does not subtract the potentials, which are above the target.
// GLIFT3.DATA is GLIFT1.DATA with an ORAT limit on B-1H in report step 5.
BOOST_AUTO_TEST_CASE(AlqLimitedToMaxWithWellTarget)
{
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using GasLiftGroupInfo = Opm::GasLiftGroupInfo<double, IndexTraits>;
    using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
    using Synthetic = SyntheticGasLiftWell<IndexTraits>;

    auto simulator = Opm::initSimulator<TypeTag>("GLIFT3.DATA", "test_glift3_maxalq", /*threads=*/2);
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);
    simulator->problem().resetIterationForNewTimestep();
    WellModel& well_model = simulator->problem().wellModel();

    simulator->setEpisodeIndex(5);
    well_model.beginReportStep(5);
    well_model.beginTimeStep();
    auto logger_guard = well_model.groupStateHelper().pushLogger();
    auto& deferred_logger = well_model.groupStateHelper().deferredLogger();
    well_model.calculateExplicitQuantities();

    const auto& well = well_model.getWell("B-1H");
    const auto& schedule = simulator->vanguard().schedule();
    const auto& summary_state = simulator->vanguard().summaryState();
    auto& well_state = well_model.wellState();
    const auto& group_state = well_model.groupState();
    const auto& comm = simulator->vanguard().grid().comm();
    const double increment = schedule.glo(5).gaslift_increment();
    const double high_alq = 25.0 * increment; // above the maximum ALQ of the synthetic well
    typename Synthetic::GLiftSyncGroups sync_groups;

    const auto oil_pos = well_state.phaseUsageInfo().canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    const double target = well.wellEcl().productionControls(summary_state).oil_rate;
    well_state.well("B-1H").alq_state.set(high_alq);
    const double old_pot = well_state.well("B-1H").well_potentials[oil_pos];
    BOOST_REQUIRE_GT(old_pot, target);

    GLiftEclWells ecl_well_map;
    Opm::BlackoilWellModelGasLift<TypeTag>::
        initGliftEclWellMap(well_model.localNonshutWells(), ecl_well_map);
    const auto& iterCtx = simulator->problem().iterationContext();
    GasLiftGroupInfo group_info {ecl_well_map, schedule, summary_state, 5, iterCtx,
                                 deferred_logger, well_state, group_state, comm,
                                 /*glift_debug=*/false};
    group_info.initialize();
    const double group_oil = group_info.oilRate("PLAT-A");
    Synthetic glift {well, deferred_logger, well_state, group_state, summary_state,
                     group_info, schedule, 5, sync_groups, comm, 0.0, 1.0e6};
    glift.setOrigAlq(high_alq);
    glift.limitAlq();
    const double new_oil = glift.initialOilRate();
    const double new_pot = well_state.well("B-1H").well_potentials[oil_pos];
    // The test needs the potential to be above the target also with the limited ALQ, and to
    // differ from the old one
    BOOST_REQUIRE_GT(new_pot, target);
    BOOST_CHECK(std::abs(new_pot - old_pot) > 1e-12);
    // Limited by the target before and after: the group oil rate does not change
    BOOST_CHECK_SMALL(new_oil - target, 1e-10);
    BOOST_CHECK_SMALL(group_info.oilRate("PLAT-A") - group_oil, 1e-10);
    // doing it again changes nothing
    glift.initialOilRate();
    BOOST_CHECK_SMALL(group_info.oilRate("PLAT-A") - group_oil, 1e-10);
}

// When the ALQ is limited to the maximum and a group target limits the rate of the well, the
// reconciled rates are the baseline of both the increase and the decrease pass of the
// optimization: the second pass must not rebuild unlimited rates from the updated potentials,
// the cached group rates are those of the first.
// GLIFT4.DATA is GLIFT1.DATA with an ORAT target on PLAT-A in report step 5.  The oil rate of the
// group is above it, and the potential of B-1H is set below the rate it has with the maximum ALQ,
// so the group does not allow the well to increase its rate.
BOOST_AUTO_TEST_CASE(AlqLimitedToMaxWithGroupTarget)
{
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using GasLiftGroupInfo = Opm::GasLiftGroupInfo<double, IndexTraits>;
    using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
    using Synthetic = SyntheticGasLiftWell<IndexTraits>;

    auto simulator = Opm::initSimulator<TypeTag>("GLIFT4.DATA", "test_glift4_maxalq", /*threads=*/2);
    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);
    simulator->problem().resetIterationForNewTimestep();
    WellModel& well_model = simulator->problem().wellModel();

    simulator->setEpisodeIndex(5);
    well_model.beginReportStep(5);
    well_model.beginTimeStep();
    auto logger_guard = well_model.groupStateHelper().pushLogger();
    auto& deferred_logger = well_model.groupStateHelper().deferredLogger();
    well_model.calculateExplicitQuantities();

    const auto& well = well_model.getWell("B-1H");
    const auto& schedule = simulator->vanguard().schedule();
    const auto& summary_state = simulator->vanguard().summaryState();
    auto& well_state = well_model.wellState();
    const auto& group_state = well_model.groupState();
    const auto& comm = simulator->vanguard().grid().comm();
    const double increment = schedule.glo(5).gaslift_increment();
    const double high_alq = 25.0 * increment; // above the maximum ALQ of the synthetic well
    typename Synthetic::GLiftSyncGroups sync_groups;

    const auto oil_pos = well_state.phaseUsageInfo().canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    well_state.well("B-1H").alq_state.set(high_alq);
    well_state.well("B-1H").well_potentials[oil_pos] = 0.01;
    GLiftEclWells ecl_well_map;
    Opm::BlackoilWellModelGasLift<TypeTag>::
        initGliftEclWellMap(well_model.localNonshutWells(), ecl_well_map);
    const auto& iterCtx = simulator->problem().iterationContext();
    GasLiftGroupInfo group_info {ecl_well_map, schedule, summary_state, 5, iterCtx,
                                 deferred_logger, well_state, group_state, comm,
                                 /*glift_debug=*/false};
    group_info.initialize();
    const double target = *group_info.oilTarget("PLAT-A");
    const double initial_group_oil = group_info.oilRate("PLAT-A");
    BOOST_REQUIRE_GT(initial_group_oil, target);
    Synthetic glift {well, deferred_logger, well_state, group_state, summary_state,
                     group_info, schedule, 5, sync_groups, comm, 0.0, 1.0e6};
    glift.setOrigAlq(high_alq);
    glift.limitAlq();
    const double first = glift.initialOilRate();
    const double group_oil = group_info.oilRate("PLAT-A");
    const double new_pot = well_state.well("B-1H").well_potentials[oil_pos];
    // the group limits the rate of the well in the first pass to the rate it had, which the group
    // rate holds
    BOOST_REQUIRE_LT(first, new_pot - 1e-6);
    BOOST_CHECK_SMALL(first - 0.01, 1e-10);
    BOOST_CHECK_SMALL(group_oil - initial_group_oil, 1e-10);
    // the second pass starts from the same rates, and the group rates do not change
    const double second = glift.initialOilRate();
    BOOST_CHECK_SMALL(second - first, 1e-10);
    BOOST_CHECK_SMALL(group_info.oilRate("PLAT-A") - group_oil, 1e-10);
}
