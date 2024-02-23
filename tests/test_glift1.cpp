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
#include "config.h"

#define BOOST_TEST_MODULE Glift1

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <ebos/equil/equilibrationhelpers.hh>
#include <ebos/eclproblem.hh>
#include <ebos/ebos.hh>
#include <opm/models/utils/start.hh>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/flow/BlackoilModel.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/WellState.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <exception>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

namespace Opm::Properties {
    namespace TTag {
        struct TestGliftTypeTag {
            using InheritsFrom = std::tuple<EbosTypeTag>;
        };
    }
}

template <class TypeTag>
std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>>
initSimulator(const char *filename)
{
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    std::string filename_arg = "--ecl-deck-file-name=";
    filename_arg += filename;

    const char* argv[] = {
        "test_equil",
        filename_arg.c_str()
    };

    Opm::registerEclTimeSteppingParameters<TypeTag>();
    Opm::setupParameters_<TypeTag>(/*argc=*/sizeof(argv)/sizeof(argv[0]), argv, /*registerParams=*/true);

    Opm::EclGenericVanguard::readDeck(filename);

    return std::make_unique<Simulator>();
}


namespace {

struct GliftFixture {
    GliftFixture() {
    int argc = boost::unit_test::framework::master_test_suite().argc;
    char** argv = boost::unit_test::framework::master_test_suite().argv;
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif
        Opm::EclGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
        using TypeTag = Opm::Properties::TTag::FlowProblem;
        Opm::registerAllParameters_<TypeTag>();
    }
};

}

BOOST_GLOBAL_FIXTURE(GliftFixture);

BOOST_AUTO_TEST_CASE(G1)
{
    //using TypeTag = Opm::Properties::TTag::FlowProblem;
    using TypeTag = Opm::Properties::TTag::TestGliftTypeTag;
    //using EclProblem = Opm::EclProblem<TypeTag>;
    //using EclWellModel = typename EclProblem::EclWellModel;
    using WellModel = Opm::BlackoilWellModel<TypeTag>;
    using WellState = Opm::WellState;
    using StdWell = Opm::StandardWell<TypeTag>;
    using GasLiftSingleWell = Opm::GasLiftSingleWell<TypeTag>;
    using GasLiftGroupInfo = Opm::GasLiftGroupInfo;
    using GasLiftSingleWellGeneric = Opm::GasLiftSingleWellGeneric;
    using GLiftEclWells = typename GasLiftGroupInfo::GLiftEclWells;
    const std::string filename = "GLIFT1.DATA";
    using GLiftSyncGroups = typename GasLiftSingleWellGeneric::GLiftSyncGroups;

    auto simulator = initSimulator<TypeTag>(filename.data());

    simulator->model().applyInitialSolution();
    simulator->setEpisodeIndex(-1);
    simulator->setEpisodeLength(0.0);
    simulator->startNextEpisode(/*episodeStartTime=*/0.0, /*episodeLength=*/1e30);
    simulator->setTimeStepSize(43200);  // 12 hours
    simulator->model().newtonMethod().setIterationIndex(0);
    WellModel& well_model = simulator->problem().wellModel();
    int report_step_idx = 0;
    well_model.beginReportStep(report_step_idx);
    well_model.beginTimeStep();
    Opm::DeferredLogger deferred_logger;
    well_model.calculateExplicitQuantities(deferred_logger);
    well_model.prepareTimeStep(deferred_logger);
    well_model.updateWellControls(false, deferred_logger);
    well_model.initPrimaryVariablesEvaluation();
    Opm::WellInterface<TypeTag> *well_ptr = well_model.getWell("B-1H").get();
    StdWell *std_well = dynamic_cast<StdWell *>(well_ptr);

    const auto& schedule = simulator->vanguard().schedule();
    auto wells_ecl = schedule.getWells(report_step_idx);
    std::optional<std::size_t> idx;
    int num_producers = 0;
    for(std::size_t i = 0; i < wells_ecl.size(); ++i) {
        const auto &well =  wells_ecl[i];
        if (well.isProducer()) {
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
    well_model.initGliftEclWellMap(ecl_well_map);
    const int iteration_idx = simulator->model().newtonMethod().numIterations();
    const auto& comm = simulator->vanguard().grid().comm();
    GasLiftGroupInfo group_info {
        ecl_well_map,
        schedule,
        summary_state,
        simulator->episodeIndex(),
        iteration_idx,
        well_model.phaseUsage(),
        deferred_logger,
        well_state,
        group_state,
        comm,
        /*glift_debug=*/false
    };
    GLiftSyncGroups sync_groups;
    GasLiftSingleWell glift {*std_well, *(simulator.get()), summary_state,
        deferred_logger, well_state, group_state, group_info, sync_groups,
        comm, /*glift_debug=*/false
    };
    group_info.initialize();
    auto state = glift.runOptimize(iteration_idx);
    BOOST_CHECK_CLOSE(state->oilRate(), 0.01736111111111111, 1e-8);
    BOOST_CHECK_CLOSE(state->gasRate(), 1.6464, 1e-3);
    BOOST_CHECK(state->oilIsLimited());
    BOOST_CHECK(!state->gasIsLimited());
    BOOST_CHECK(!state->alqIsLimited());
    BOOST_CHECK_CLOSE(state->alq(), 0.0, 1e-8);
    BOOST_CHECK(!state->increase().has_value());
}

