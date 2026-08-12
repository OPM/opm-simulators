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

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQParams.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/grid/GridManager.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/common/utility/TimeService.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/grid/GridHelpers.hpp>
#include <opm/simulators/flow/FlowMain.hpp>
#include <opm/simulators/flow/NonlinearSystemBlackOilReservoir.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <opm/models/utils/start.hh>

#include <opm/simulators/wells/StandardWell.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace Opm::Properties {

namespace TTag {

struct WellModelTestTypeTag
{ using InheritsFrom = std::tuple<FlowProblem>; };

}

// Disable convective mixing
template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::WellModelTestTypeTag>
{ static constexpr bool value = false; };

// Disable diffusion
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::WellModelTestTypeTag>
{ static constexpr bool value = false; };

}

using StandardWell = Opm::StandardWell<Opm::Properties::TTag::WellModelTestTypeTag>;

struct SetupTest {

    using Grid = UnstructuredGrid;

    SetupTest()
    {
        const auto deck = Opm::Parser{}.parseFile("TESTWELLMODEL.DATA");
        this->ecl_state = std::make_unique<const Opm::EclipseState>(deck);

        const Opm::TableManager table(deck);
        const Opm::Runspec runspec(deck);

        this->schedule = std::make_unique<const Opm::Schedule>
            (deck, *this->ecl_state, std::make_shared<Opm::Python>());

        this->summaryState = std::make_unique<Opm::SummaryState>
            (Opm::TimeService::from_time_t(schedule->getStartTime()),
             this->ecl_state->runspec().udqParams().undefinedValue());

        current_timestep = 0;
    };

    std::unique_ptr<const Opm::EclipseState> ecl_state;
    std::shared_ptr<Opm::Python> python;
    std::unique_ptr<const Opm::Schedule> schedule;
    std::unique_ptr<Opm::SummaryState> summaryState;
    std::vector<std::vector<Opm::PerforationData<double>>> well_perf_data;
    int current_timestep;
};

struct GlobalFixture {
    GlobalFixture()
    {
        int argcDummy = 1;
        const char *tmp[] = {"test_wellmodel", nullptr};
        char **argvDummy = const_cast<char**>(tmp);

        // MPI setup.
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argcDummy, argvDummy);
#else
        Dune::MPIHelper::instance(argcDummy, argvDummy);
#endif

        Opm::FlowMain<Opm::Properties::TTag::WellModelTestTypeTag>::setupParameters_(argcDummy, argvDummy, Dune::MPIHelper::getCommunication());
    }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_CASE(TestStandardWellInput) {
    const SetupTest setup_test;
    const auto& wells_ecl = setup_test.schedule->getWells(setup_test.current_timestep);
    BOOST_CHECK_EQUAL( wells_ecl.size(), 2);
    const Opm::Well& well = wells_ecl[1];
    const Opm::BlackoilModelParameters<double> param;

    // For the conversion between the surface volume rate and resrevoir voidage rate
    typedef Opm::BlackOilFluidSystem<double> FluidSystem;
    using RateConverterType = Opm::RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
    // Compute reservoir volumes for RESV controls.
    std::unique_ptr<RateConverterType> rateConverter;
    // Compute reservoir volumes for RESV controls.
    rateConverter.reset(new RateConverterType (std::vector<int>(10, 0)));

    Opm::PerforationData<double> dummy;
    std::vector<Opm::PerforationData<double>> pdata(well.getConnections().size(), dummy);
    for (auto c = 0*pdata.size(); c < pdata.size(); ++c) {
        pdata[c].ecl_index = c;
    }

    Opm::ParallelWellInfo<double> pinfo{well.name()};

    BOOST_CHECK_THROW( StandardWell( well, pinfo, -1, param, *rateConverter, 0, 3, 3, 0, pdata), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestBehavoir) {
    const SetupTest setup_test;
    const auto& wells_ecl = setup_test.schedule->getWells(setup_test.current_timestep);
    const int current_timestep = setup_test.current_timestep;
    std::vector<std::unique_ptr<const StandardWell> >  wells;

    {
        const int nw = wells_ecl.size();
        const Opm::BlackoilModelParameters<double> param;

        for (int w = 0; w < nw; ++w) {
            // For the conversion between the surface volume rate and resrevoir voidage rate
            using FluidSystem = Opm::BlackOilFluidSystem<double>;
            using RateConverterType = Opm::RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
            // Compute reservoir volumes for RESV controls.
            // TODO: not sure why for this class the initlizer list does not work
            // otherwise we should make a meaningful const PhaseUsage here.
            std::unique_ptr<RateConverterType> rateConverter;
            // Compute reservoir volumes for RESV controls.
            rateConverter.reset(new RateConverterType (std::vector<int>(10, 0)));
            Opm::PerforationData<double> dummy;
            std::vector<Opm::PerforationData<double>> pdata(wells_ecl[w].getConnections().size(), dummy);
            for (auto c = 0*pdata.size(); c < pdata.size(); ++c) {
                pdata[c].ecl_index = c;
            }

            Opm::ParallelWellInfo<double> pinfo{wells_ecl[w].name()};
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

BOOST_AUTO_TEST_CASE(TestPrimaryVariableScaling) {
    // The scaling must appear in the derivative only: the stored value and the
    // Evaluation's value stay physical. Set the parameters before the first
    // varScale() call - the scales are cached statically.
    // 2^16 rather than the recommended 2^23: SetDefault round-trips the value
    // through text at 6 significant digits, so it must be exactly representable
    // there (8388608 would arrive as 8388610). Command-line parsing is exact.
    Opm::Parameters::SetDefault<Opm::Parameters::WellBhpScaling<double>>(65536.0); // 2^16
    Opm::Parameters::SetDefault<Opm::Parameters::WellRateScaling<double>>(0.25);   // 2^-2

    const SetupTest setup_test;
    const auto& wells_ecl = setup_test.schedule->getWells(setup_test.current_timestep);
    const Opm::BlackoilModelParameters<double> param;

    using FluidSystem = Opm::BlackOilFluidSystem<double>;
    // Just enough initialisation for phase-usage queries to work
    // (same device as test_rftcontainer.cpp).
    FluidSystem::initBegin(/*numPvtRegions=*/1);
    using RateConverterType = Opm::RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int>>;
    RateConverterType rateConverter(std::vector<int>(10, 0));

    const auto& well_ecl = wells_ecl[0]; // PROD1
    Opm::PerforationData<double> dummy;
    std::vector<Opm::PerforationData<double>> pdata(well_ecl.getConnections().size(), dummy);
    for (auto c = 0*pdata.size(); c < pdata.size(); ++c) {
        pdata[c].ecl_index = c;
    }
    Opm::ParallelWellInfo<double> pinfo{well_ecl.name()};
    const StandardWell well(well_ecl, pinfo, setup_test.current_timestep,
                            param, rateConverter, 0, 3, 3, 0, pdata);

    constexpr double pv_bhp_scale = 65536.0;
    using PV = std::decay_t<decltype(well.primaryVariables())>;
    PV pv(well);
    pv.resize(well.numStaticWellEq);

    // Zero Newton update through the public interface triggers
    // setEvaluationsFromValues(); the bhp lands on its lower limit.
    typename PV::BVectorWell dwells(1);
    dwells[0].resize(well.numStaticWellEq);
    dwells[0] = 0.0;
    Opm::DeferredLogger logger;
    pv.updateNewton(dwells, /*stop_or_zero_rate_target=*/false,
                    /*dFLimit=*/0.2, /*dBHPLimit=*/0.1, logger);

    constexpr int numEq = StandardWell::Indices::numEq;

    // Values are physical, with and without scaling.
    BOOST_CHECK_EQUAL(pv.eval(PV::Bhp).value(), pv.value(PV::Bhp));
    BOOST_CHECK_EQUAL(pv.eval(PV::WQTotal).value(), pv.value(PV::WQTotal));

    // The derivatives carry the scale: d(x)/d(x/s) = s. These fail without the
    // scaling support (both were hard-coded 1.0).
    BOOST_CHECK_EQUAL(pv.eval(PV::Bhp).derivative(numEq + PV::Bhp), 65536.0);
    BOOST_CHECK_EQUAL(pv.eval(PV::WQTotal).derivative(numEq + PV::WQTotal), 0.25);

    // The dimensionless fractions stay unscaled.
    BOOST_CHECK_EQUAL(pv.eval(PV::WFrac).derivative(numEq + PV::WFrac), 1.0);
    BOOST_CHECK_EQUAL(pv.eval(PV::GFrac).derivative(numEq + PV::GFrac), 1.0);

    // Every slot: the stored value and the Evaluation agree, and only the own
    // column carries a derivative. This is the assertion a forgotten conversion
    // in setEvaluationsFromValues would break.
    for (int i = 0; i < well.numStaticWellEq; ++i) {
        BOOST_CHECK_EQUAL(pv.eval(i).value(), pv.value(i));
        for (int j = 0; j < well.numStaticWellEq; ++j) {
            if (i != j) {
                BOOST_CHECK_EQUAL(pv.eval(i).derivative(numEq + j), 0.0);
            }
        }
    }

    // value()/setValue() must be exact inverses, or the getPrimaryVars round
    // trip used by NLDD would rescale on every save/restore.
    const double bhp_phys = 300.0 * Opm::unit::barsa;
    pv.setValue(PV::Bhp, bhp_phys);
    BOOST_CHECK_EQUAL(pv.value(PV::Bhp), bhp_phys);

    // Newton step with the absolute bhp floor ACTIVE. The limit is physical, the
    // increment arrives scaled; a missing conversion on either side lands
    // somewhere other than exactly the floor. The unscaled test above never
    // reaches this branch because its increment is zero.
    constexpr double bhp_floor = 1.0 * Opm::unit::barsa - 1.0 * Opm::unit::Pascal;
    pv.setValue(PV::Bhp, 1.5 * Opm::unit::barsa);
    dwells[0][PV::Bhp] = (1.0 * Opm::unit::barsa) / pv_bhp_scale; // scaled: drives well below the floor
    pv.updateNewton(dwells, false, 0.2, 1.0, logger);
    BOOST_CHECK_EQUAL(pv.value(PV::Bhp), bhp_floor);
}
