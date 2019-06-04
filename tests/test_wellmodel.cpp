/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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

#define BOOST_TEST_MODULE WellModelTest

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/wells.h>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/grid/GridHelpers.hpp>
#include <opm/autodiff/FlowMainEbos.hpp>
#include <opm/autodiff/BlackoilModelEbos.hpp>
#include <opm/autodiff/createGlobalCellArray.hpp>
#include <opm/autodiff/GridInit.hpp>

#include <ebos/eclproblem.hh>
#include <ewoms/common/start.hh>

#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

using StandardWell = Opm::StandardWell<TTAG(EclFlowProblem)>;

struct SetupTest {

    using Grid = UnstructuredGrid;
    using GridInit = Opm::GridInit<Grid>;

    SetupTest ()
    {
        Opm::Parser parser;
        auto deck = parser.parseFile("TESTWELLMODEL.DATA");
        ecl_state.reset(new Opm::EclipseState(deck) );
        {
          const Opm::TableManager table ( deck );
          const Opm::Eclipse3DProperties eclipseProperties ( deck , table, ecl_state->getInputGrid());
          const Opm::Runspec runspec (deck);
          schedule.reset( new Opm::Schedule(deck, ecl_state->getInputGrid(), eclipseProperties, runspec));
        }

        // Create grid.
        const std::vector<double>& porv =
                            ecl_state->get3DProperties().getDoubleGridProperty("PORV").getData();

        std::unique_ptr<GridInit> grid_init(new GridInit(*ecl_state, porv));
        const Grid& grid = grid_init->grid();

        // Create material law manager.
        std::vector<int> compressed_to_cartesianIdx;
        Opm::createGlobalCellArray(grid, compressed_to_cartesianIdx);

        current_timestep = 0;

        // Create wells.
        wells_manager.reset(new Opm::WellsManager(*ecl_state,
                                                  *schedule,
                                                  summaryState,
                                                  current_timestep,
                                                  Opm::UgGridHelpers::numCells(grid),
                                                  Opm::UgGridHelpers::globalCell(grid),
                                                  Opm::UgGridHelpers::cartDims(grid),
                                                  Opm::UgGridHelpers::dimensions(grid),
                                                  Opm::UgGridHelpers::cell2Faces(grid),
                                                  Opm::UgGridHelpers::beginFaceCentroids(grid),
                                                  false,
                                                  std::unordered_set<std::string>() ) );

    };

    std::unique_ptr<const Opm::WellsManager> wells_manager;
    std::unique_ptr<const Opm::EclipseState> ecl_state;
    std::unique_ptr<const Opm::Schedule> schedule;
    Opm::SummaryState summaryState;
    int current_timestep;
};

struct GlobalFixture {
    GlobalFixture()
    {
        int argcDummy = 1;
        const char *tmp[] = {"test_wellmodel"};
        char **argvDummy = const_cast<char**>(tmp);

        // MPI setup.
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argcDummy, argvDummy);
#else
        Dune::MPIHelper::instance(argcDummy, argvDummy);
#endif

        Opm::FlowMainEbos<TTAG(EclFlowProblem)>::setupParameters_(argcDummy, argvDummy);
    }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_CASE(TestStandardWellInput) {
    const SetupTest setup_test;
    const Wells* wells = setup_test.wells_manager->c_wells();
    const auto& wells_ecl = setup_test.schedule->getWells2(setup_test.current_timestep);
    BOOST_CHECK_EQUAL( wells_ecl.size(), 2);
    const Opm::Well2& well = wells_ecl[1];
    const Opm::BlackoilModelParametersEbos<TTAG(EclFlowProblem) > param;

    // For the conversion between the surface volume rate and resrevoir voidage rate
    typedef Opm::BlackOilFluidSystem<double> FluidSystem;
    using RateConverterType = Opm::RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
    // Compute reservoir volumes for RESV controls.
    Opm::PhaseUsage phaseUsage;
    std::unique_ptr<RateConverterType> rateConverter;
    // Compute reservoir volumes for RESV controls.
    rateConverter.reset(new RateConverterType (phaseUsage,
                                     std::vector<int>(10, 0)));

    const int pvtIdx = 0;
    const int num_comp = wells->number_of_phases;

    BOOST_CHECK_THROW( StandardWell( well, -1, wells, param, *rateConverter, pvtIdx, num_comp), std::invalid_argument);
    BOOST_CHECK_THROW( StandardWell( well, 4, nullptr , param, *rateConverter, pvtIdx, num_comp), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE(TestBehavoir) {
    const SetupTest setup_test;
    const Wells* wells_struct = setup_test.wells_manager->c_wells();
    const auto& wells_ecl = setup_test.schedule->getWells2(setup_test.current_timestep);
    const int current_timestep = setup_test.current_timestep;
    std::vector<std::unique_ptr<const StandardWell> >  wells;

    {
        const int nw = wells_struct ? (wells_struct->number_of_wells) : 0;
        const Opm::BlackoilModelParametersEbos<TTAG(EclFlowProblem)> param;

        for (int w = 0; w < nw; ++w) {
            const std::string well_name(wells_struct->name[w]);

            size_t index_well = 0;
            for (; index_well < wells_ecl.size(); ++index_well) {
                if (well_name == wells_ecl[index_well].name()) {
                    break;
                }
            }
            // we should always be able to find the well in wells_ecl
            BOOST_CHECK(index_well !=  wells_ecl.size());
            // For the conversion between the surface volume rate and resrevoir voidage rate
            typedef Opm::BlackOilFluidSystem<double> FluidSystem;
            using RateConverterType = Opm::RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
            // Compute reservoir volumes for RESV controls.
            // TODO: not sure why for this class the initlizer list does not work
            // otherwise we should make a meaningful const PhaseUsage here.
            Opm::PhaseUsage phaseUsage;
            std::unique_ptr<RateConverterType> rateConverter;
            // Compute reservoir volumes for RESV controls.
            rateConverter.reset(new RateConverterType (phaseUsage,
                                             std::vector<int>(10, 0)));

            const int pvtIdx = 0;
            const int num_comp = wells_struct->number_of_phases;

            wells.emplace_back(new StandardWell(wells_ecl[index_well], current_timestep, wells_struct, param, *rateConverter, pvtIdx, num_comp) );
        }
    }

    // first well, it is a production well from the deck
    {
        const auto& well = wells[0];
        BOOST_CHECK_EQUAL(well->name(), "PROD1");
        BOOST_CHECK(well->wellType() == PRODUCER);
        BOOST_CHECK(well->numEq == 3);
        BOOST_CHECK(well->numWellEq == 4);
        const auto& wc = well->wellControls();
        const int ctrl_num = well_controls_get_num(wc);
        BOOST_CHECK(ctrl_num > 0);
        const auto& control = well_controls_get_current(wc);
        BOOST_CHECK(control >= 0);
        // GAS RATE CONTROL
        const auto& distr = well_controls_iget_distr(wc, control);
        BOOST_CHECK(distr[0] == 0.);
        BOOST_CHECK(distr[1] == 0.);
        BOOST_CHECK(distr[2] == 1.);
    }

    // second well, it is the injection well from the deck
    {
        const auto& well = wells[1];
        BOOST_CHECK_EQUAL(well->name(), "INJE1");
        BOOST_CHECK(well->wellType() == INJECTOR);
        BOOST_CHECK(well->numEq == 3);
        BOOST_CHECK(well->numWellEq == 4);
        const auto& wc = well->wellControls();
        const int ctrl_num = well_controls_get_num(wc);
        BOOST_CHECK(ctrl_num > 0);
        const auto& control = well_controls_get_current(wc);
        BOOST_CHECK(control >= 0);
        // WATER RATE CONTROL
        const auto& distr = well_controls_iget_distr(wc, control);
        BOOST_CHECK(distr[0] == 1.);
        BOOST_CHECK(distr[1] == 0.);
        BOOST_CHECK(distr[2] == 0.);
    }
}
