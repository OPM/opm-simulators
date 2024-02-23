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
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/grid/GridManager.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/common/utility/TimeService.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/grid/GridHelpers.hpp>
#include <opm/simulators/flow/FlowMain.hpp>
#include <opm/simulators/flow/BlackoilModel.hpp>

#include <ebos/eclproblem.hh>
#include <opm/models/utils/start.hh>

#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif



using StandardWell = Opm::StandardWell<Opm::Properties::TTag::FlowProblem>;

struct SetupTest {

    using Grid = UnstructuredGrid;

    SetupTest ()
    {
        Opm::Parser parser;
        auto deck = parser.parseFile("TESTWELLMODEL.DATA");
        ecl_state.reset(new Opm::EclipseState(deck) );
        {
          const Opm::TableManager table ( deck );
          const Opm::Runspec runspec (deck);
          python = std::make_shared<Opm::Python>();
          schedule.reset( new Opm::Schedule(deck, *ecl_state, python));
          summaryState.reset( new Opm::SummaryState(Opm::TimeService::from_time_t(schedule->getStartTime())));
        }
        current_timestep = 0;
    };

    std::unique_ptr<const Opm::EclipseState> ecl_state;
    std::shared_ptr<Opm::Python> python;
    std::unique_ptr<const Opm::Schedule> schedule;
    std::unique_ptr<Opm::SummaryState> summaryState;
    std::vector<std::vector<Opm::PerforationData>> well_perf_data;
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

        Opm::FlowMain<Opm::Properties::TTag::FlowProblem>::setupParameters_(argcDummy, argvDummy, Dune::MPIHelper::getCommunication());
    }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_CASE(TestStandardWellInput) {
    const SetupTest setup_test;
    const auto& wells_ecl = setup_test.schedule->getWells(setup_test.current_timestep);
    BOOST_CHECK_EQUAL( wells_ecl.size(), 2);
    const Opm::Well& well = wells_ecl[1];
    const Opm::BlackoilModelParameters<Opm::Properties::TTag::FlowProblem> param;

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

    Opm::PerforationData dummy;
    std::vector<Opm::PerforationData> pdata(well.getConnections().size(), dummy);
    for (auto c = 0*pdata.size(); c < pdata.size(); ++c) {
        pdata[c].ecl_index = c;
    }

    Opm::ParallelWellInfo pinfo{well.name()};

    BOOST_CHECK_THROW( StandardWell( well, pinfo, -1, param, *rateConverter, 0, 3, 3, 0, pdata), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE(TestBehavoir) {
    const SetupTest setup_test;
    const auto& wells_ecl = setup_test.schedule->getWells(setup_test.current_timestep);
    const int current_timestep = setup_test.current_timestep;
    std::vector<std::unique_ptr<const StandardWell> >  wells;

    {
        const int nw = wells_ecl.size();
        const Opm::BlackoilModelParameters<Opm::Properties::TTag::FlowProblem> param;

        for (int w = 0; w < nw; ++w) {
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
            Opm::PerforationData dummy;
            std::vector<Opm::PerforationData> pdata(wells_ecl[w].getConnections().size(), dummy);
            for (auto c = 0*pdata.size(); c < pdata.size(); ++c) {
                pdata[c].ecl_index = c;
            }

            Opm::ParallelWellInfo pinfo{wells_ecl[w].name()};
            wells.emplace_back(new StandardWell(wells_ecl[w], pinfo, current_timestep, param, *rateConverter, 0, 3, 3, w, pdata) );
        }
    }

    // first well, it is a production well from the deck
    {
        const auto& well = wells[0];
        BOOST_CHECK_EQUAL(well->name(), "PROD1");
        BOOST_CHECK(well->isProducer());
        BOOST_CHECK(StandardWell::Indices::numEq == 3);
        BOOST_CHECK(well->numStaticWellEq== 4);
    }

    // second well, it is the injection well from the deck
    {
        const auto& well = wells[1];
        BOOST_CHECK_EQUAL(well->name(), "INJE1");
        BOOST_CHECK(well->isInjector());
        BOOST_CHECK(StandardWell::Indices::numEq == 3);
        BOOST_CHECK(well->numStaticWellEq== 4);      
    }
}
