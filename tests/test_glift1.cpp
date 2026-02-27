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
