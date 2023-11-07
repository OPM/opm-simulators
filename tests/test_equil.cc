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

#define BOOST_TEST_MODULE Equil

#include <ebos/equil/equilibrationhelpers.hh>
#include <ebos/eclproblem.hh>
#include <ebos/eclgenericvanguard.hh>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/GridManager.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/models/utils/start.hh>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>

#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <string.h>

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif


namespace Opm::Properties {

template<class TypeTag>
struct EnableTerminalOutput<TypeTag, TTag::EclBaseProblem> {
    static constexpr bool value = true;
};

namespace TTag {


struct TestEquilTypeTag {
    using InheritsFrom = std::tuple<FlowTimeSteppingParameters, FlowModelParameters, EclBaseProblem, BlackOilModel>;
};
struct TestEquilVapwatTypeTag {
    using InheritsFrom = std::tuple<FlowModelParameters, EclBaseProblem, BlackOilModel>;
};
}

template<class TypeTag>
struct EclWellModel<TypeTag, TTag::TestEquilTypeTag> {
    using type = BlackoilWellModel<TypeTag>;
};
template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::TestEquilTypeTag> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EclWellModel<TypeTag, TTag::TestEquilVapwatTypeTag> {
    using type = BlackoilWellModel<TypeTag>;
};
template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::TestEquilVapwatTypeTag> {
    static constexpr bool value = true;
};
} // namespace Opm::Properties

template <class TypeTag>
std::unique_ptr<Opm::GetPropType<TypeTag, Opm::Properties::Simulator>>
initSimulator(const char *filename)
{
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;

    std::string filenameArg = "--ecl-deck-file-name=";
    filenameArg += filename;

    const char* argv[] = {
        "test_equil",
        filenameArg.c_str()
    };

    Opm::setupParameters_<TypeTag>(/*argc=*/sizeof(argv)/sizeof(argv[0]), argv, /*registerParams=*/false);

    Opm::EclGenericVanguard::readDeck(filename);

    return std::make_unique<Simulator>();
}

template <class GridView>
static std::vector<std::pair<double,double>> cellVerticalExtent(const GridView& gridView)
{
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());

    int numElements = gridView.size(/*codim=*/0);
    std::vector<std::pair<double,double>> cellZMinMax(numElements);

    auto elemIt = gridView.template begin</*codim=*/0>();
    const auto& elemEndIt = gridView.template end</*codim=*/0>();
    for (; elemIt != elemEndIt; ++elemIt) {
        const auto& element = *elemIt;
        const unsigned int elemIdx = elemMapper.index(element);
        cellZMinMax[elemIdx] = Opm::EQUIL::Details::cellZMinMax(element);
    }
    return cellZMinMax;
}

template <class TypeTag>
static void initDefaultFluidSystem()
{
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;

    std::vector<std::pair<double, double> > Bo = {
        { 101353, 1. },
        { 6.21542e+07, 1 }
    };
    std::vector<std::pair<double, double> > muo = {
        { 101353, 1. },
        { 6.21542e+07, 1 }
    };

    std::vector<std::pair<double, double> > Bg = {
        { 101353, 1. },
        { 6.21542e+07, 1 }
    };
    std::vector<std::pair<double, double> > mug = {
        { 101353, 1. },
        { 6.21542e+07, 1 }
    };

    double rhoRefO = 700; // [kg/m3]
    double rhoRefG = 1000; // [kg/m3]
    double rhoRefW = 1000; // [kg/m3]

    FluidSystem::initBegin(/*numPvtRegions=*/1);
    FluidSystem::setEnableDissolvedGas(false);
    FluidSystem::setEnableVaporizedOil(false);
    FluidSystem::setReferenceDensities(rhoRefO, rhoRefW, rhoRefG, /*regionIdx=*/0);

    auto gasPvt = std::make_shared<Opm::GasPvtMultiplexer<double>>();
    gasPvt->setApproach(Opm::GasPvtApproach::DryGas);
    auto& dryGasPvt = gasPvt->getRealPvt<Opm::GasPvtApproach::DryGas>();
    dryGasPvt.setNumRegions(/*numPvtRegion=*/1);
    dryGasPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
    dryGasPvt.setGasFormationVolumeFactor(/*regionIdx=*/0, Bg);
    dryGasPvt.setGasViscosity(/*regionIdx=*/0, mug);

    auto oilPvt = std::make_shared<Opm::OilPvtMultiplexer<double>>();
    oilPvt->setApproach(Opm::OilPvtApproach::DeadOil);
    auto& deadOilPvt = oilPvt->getRealPvt<Opm::OilPvtApproach::DeadOil>();
    deadOilPvt.setNumRegions(/*numPvtRegion=*/1);
    deadOilPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
    deadOilPvt.setInverseOilFormationVolumeFactor(/*regionIdx=*/0, Bo);
    deadOilPvt.setOilViscosity(/*regionIdx=*/0, muo);

    auto waterPvt = std::make_shared<Opm::WaterPvtMultiplexer<double>>();
    waterPvt->setApproach(Opm::WaterPvtApproach::ConstantCompressibilityWater);
    auto& ccWaterPvt = waterPvt->getRealPvt<Opm::WaterPvtApproach::ConstantCompressibilityWater>();
    ccWaterPvt.setNumRegions(/*numPvtRegions=*/1);
    ccWaterPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
    ccWaterPvt.setViscosity(/*regionIdx=*/0, 1);
    ccWaterPvt.setCompressibility(/*regionIdx=*/0, 0);

    gasPvt->initEnd();
    oilPvt->initEnd();
    waterPvt->initEnd();

    FluidSystem::setGasPvt(std::move(gasPvt));
    FluidSystem::setOilPvt(std::move(oilPvt));
    FluidSystem::setWaterPvt(std::move(waterPvt));
    FluidSystem::initEnd();
}

static Opm::EquilRecord mkEquilRecord( double datd, double datp,
                                       double zwoc, double pcow_woc,
                                       double zgoc, double pcgo_goc )
{
    return Opm::EquilRecord( datd, datp, zwoc, pcow_woc, zgoc, pcgo_goc, true, true, 0, true);
}

template <typename Simulator>
double centerDepth(const Simulator& sim, const std::size_t cell)
{
    return Opm::UgGridHelpers::cellCenterDepth(sim.vanguard().grid(), cell);
}

namespace {

struct EquilFixture {
    EquilFixture() {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
#else
        Dune::MPIHelper::instance(argc, argv);
#endif
        Opm::EclGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());
        Opm::BlackoilModelParametersEbos<TypeTag>::registerParameters();
        Opm::AdaptiveTimeSteppingEbos<TypeTag>::registerParameters();
        Opm::Parameters::registerParam<TypeTag, bool>("EnableTerminalOutput",
                                                      "EnableTerminalOutput",
                                                      Opm::getPropValue<TypeTag, Opm::Properties::EnableTerminalOutput>(),
                                                      "Dummy added for the well model to compile.");
        Opm::registerAllParameters_<TypeTag>();
    }

    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    using Grid = Opm::GetPropType<TypeTag, Opm::Properties::Grid>;
    using GridView = Opm::GetPropType<TypeTag, Opm::Properties::GridView>;
    using ElementMapper = Opm::GetPropType<TypeTag, Opm::Properties::ElementMapper>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using Initializer = Opm::EQUIL::DeckDependent::InitialStateComputer<FluidSystem,
                                                                        Grid,
                                                                        GridView,
                                                                        ElementMapper,
                                                                        CartesianIndexMapper>;
};

}

BOOST_GLOBAL_FIXTURE(EquilFixture);

BOOST_AUTO_TEST_CASE(PhasePressure)
{
    const auto record = mkEquilRecord( 0, 1e5, 5, 0, 0, 0 );

    using TypeTag     = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;

    std::vector<double> x = {0.0,100.0};
    std::vector<double> y = {0.0,0.0};
    Opm::Tabulated1DFunction<double> trivialSaltVdTable{2,x,y};

    auto simulator = initSimulator<TypeTag>("equil_base.DATA");
    initDefaultFluidSystem<TypeTag>();

    const auto region = Opm::EQUIL::EquilReg {
        record,
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0
    };

    auto vspan = std::array<double, 2>{};
    {
        auto cells = std::vector<int>(simulator->vanguard().grid().size(0));
        std::iota(cells.begin(), cells.end(), 0);

        Opm::EQUIL::Details::verticalExtent(cells, cellVerticalExtent(simulator->vanguard().gridView()),
                                            simulator->vanguard().gridView().comm(), vspan);
    }

    const auto grav = 10.0;
    auto ptable = Opm::EQUIL::Details::PressureTable<
        FluidSystem, Opm::EQUIL::EquilReg
    >{ grav };

    ptable.equilibrate(region, vspan);

    const auto reltol = 1.0e-6;
    const auto first  = centerDepth(*simulator, 0);
    const auto last   = centerDepth(*simulator, simulator->vanguard().grid().size(0) - 1);

    BOOST_CHECK_CLOSE(ptable.water(first),  90e3  , reltol);
    BOOST_CHECK_CLOSE(ptable.water(last) , 180e3  , reltol);
    BOOST_CHECK_CLOSE(ptable.oil  (first), 103.5e3, reltol);
    BOOST_CHECK_CLOSE(ptable.oil  (last) , 166.5e3, reltol);
}

BOOST_AUTO_TEST_CASE(CellSubset)
{
    using PVal        = std::vector<double>;
    using PPress      = std::vector<PVal>;
    using TypeTag     = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;

    auto simulator = initSimulator<TypeTag>("equil_base.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();

    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());
    initDefaultFluidSystem<TypeTag>();

    const Opm::EquilRecord record[] = { mkEquilRecord( 0, 1e5, 2.5, -0.075e5, 0, 0 ),
                                  mkEquilRecord( 5, 1.35e5, 7.5, -0.225e5, 5, 0 ) };

    std::vector<double> x = {0.0,100.0};
    std::vector<double> y = {0.0,0.0};
    Opm::Tabulated1DFunction<double> trivialSaltVdTable{2, x, y};

    const Opm::EQUIL::EquilReg region[] =
    {
        Opm::EQUIL::EquilReg(record[0],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
        ,
        Opm::EQUIL::EquilReg(record[0],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
        ,
        Opm::EQUIL::EquilReg(record[1],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
        ,
        Opm::EQUIL::EquilReg(record[1],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
    };

    const int cdim[] = { 2, 1, 2 };
    int ncoarse = cdim[0];
    for (std::size_t d = 1; d < 3; ++d) { ncoarse *= cdim[d]; }

    std::vector< std::vector<int> > cells(ncoarse);
    for (int c = 0; c < simulator->vanguard().grid().size(0); ++c) {
        int ci = c;
        const int i = ci % grid.cartdims[0];  ci /= grid.cartdims[0];
        const int j = ci % grid.cartdims[1];
        const int k = ci / grid.cartdims[1];

        const int ic = (i / (grid.cartdims[0] / cdim[0]));
        const int jc = (j / (grid.cartdims[1] / cdim[1]));
        const int kc = (k / (grid.cartdims[2] / cdim[2]));
        const int ix = ic + cdim[0]*(jc + cdim[1]*kc);

        assert ((0 <= ix) && (ix < ncoarse));
        cells[ix].push_back(c);
    }

    auto vspan = std::array<double, 2>{};
    {
        auto vspancells = std::vector<int>(simulator->vanguard().grid().size(0));
        std::iota(vspancells.begin(), vspancells.end(), 0);

        Opm::EQUIL::Details::verticalExtent(vspancells, cellVerticalExtent(simulator->vanguard().gridView()),
                                            simulator->vanguard().gridView().comm(), vspan);
    }

    const auto grav = 10.0;
    auto ptable = Opm::EQUIL::Details::PressureTable<
        FluidSystem, Opm::EQUIL::EquilReg
    >{ grav };

    auto ppress = PPress(2, PVal(simulator->vanguard().grid().size(0), 0.0));
    for (auto r = cells.begin(), e = cells.end(); r != e; ++r) {
        const int rno = int(r - cells.begin());

        ptable.equilibrate(region[rno], vspan);

        PVal::size_type i = 0;
        for (auto c = r->begin(), ce = r->end(); c != ce; ++c, ++i) {
            const auto depth = centerDepth(*simulator, *c);

            ppress[0][*c] = ptable.water(depth);
            ppress[1][*c] = ptable.oil  (depth);
        }
    }

    const int first = 0, last = simulator->vanguard().grid().size(0) - 1;
    const double reltol = 1.0e-6;
    BOOST_CHECK_CLOSE(ppress[FluidSystem::waterPhaseIdx][first] , 105e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[FluidSystem::waterPhaseIdx][last ] , 195e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[FluidSystem::oilPhaseIdx][first] , 103.5e3 , reltol);
    BOOST_CHECK_CLOSE(ppress[FluidSystem::oilPhaseIdx][last ] , 166.5e3 , reltol);
}

BOOST_AUTO_TEST_CASE(RegMapping)
{
    const Opm::EquilRecord record[] = {
        mkEquilRecord( 0, 1e5, 2.5, -0.075e5, 0, 0 ),
        mkEquilRecord( 5, 1.35e5, 7.5, -0.225e5, 5, 0 ),
    };

    using PVal        = std::vector<double>;
    using PPress      = std::vector<PVal>;
    using TypeTag     = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;

    auto simulator = initSimulator<TypeTag>("equil_base.DATA");
    initDefaultFluidSystem<TypeTag>();

    std::vector<double> x = {0.0,100.0};
    std::vector<double> y = {0.0,0.0};
    Opm::Tabulated1DFunction<double> trivialSaltVdTable{2, x, y};

    const Opm::EQUIL::EquilReg region[] =
    {
        Opm::EQUIL::EquilReg(record[0],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
        ,
        Opm::EQUIL::EquilReg(record[0],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
        ,
        Opm::EQUIL::EquilReg(record[1],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
        ,
        Opm::EQUIL::EquilReg(record[1],
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Opm::EQUIL::Miscibility::NoMixing>(),
        trivialSaltVdTable,
        0)
        };

    auto vspan = std::array<double, 2>{};
    {
        auto cells = std::vector<int>(simulator->vanguard().grid().size(0));
        std::iota(cells.begin(), cells.end(), 0);

        Opm::EQUIL::Details::verticalExtent(cells, cellVerticalExtent(simulator->vanguard().gridView()),
                                            simulator->vanguard().gridView().comm(), vspan);
    }

    const auto grav = 10.0;
    auto ptable = Opm::EQUIL::Details::PressureTable<
        FluidSystem, Opm::EQUIL::EquilReg
    >{ grav };

    std::vector<int> eqlnum(simulator->vanguard().grid().size(0));
    // [ 0 1; 2 3]
    {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                eqlnum[i*10 + j] = 0;
            }
            for (int j = 5; j < 10; ++j) {
                eqlnum[i*10 + j] = 1;
            }
        }
        for (int i = 5; i < 10; ++i) {
            for (int j = 0; j < 5; ++j) {
                eqlnum[i*10 + j] = 2;
            }
            for (int j = 5; j < 10; ++j) {
                eqlnum[i*10 + j] = 3;
            }
        }
    }

    const Opm::RegionMapping<> eqlmap(eqlnum);

    auto ppress = PPress(2, PVal(simulator->vanguard().grid().size(0), 0.0));
    for (const auto& r : eqlmap.activeRegions()) {
        ptable.equilibrate(region[r], vspan);

        for (const auto& c : eqlmap.cells(r)) {
            const auto depth = centerDepth(*simulator, c);

            ppress[0][c] = ptable.water(depth);
            ppress[1][c] = ptable.oil  (depth);
        }
    }

    const int first = 0, last = simulator->vanguard().grid().size(0) - 1;
    const double reltol = 1.0e-6;
    BOOST_CHECK_CLOSE(ppress[FluidSystem::waterPhaseIdx][first] , 105e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[FluidSystem::waterPhaseIdx][last ] , 195e3   , reltol);
    BOOST_CHECK_CLOSE(ppress[FluidSystem::oilPhaseIdx][first] , 103.5e3 , reltol);
    BOOST_CHECK_CLOSE(ppress[FluidSystem::oilPhaseIdx][last ] , 166.5e3 , reltol);
}

BOOST_AUTO_TEST_CASE(DeckAllDead)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator = initSimulator<TypeTag>("equil_deadfluids.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());
    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                    eclipseState, 
                                    simulator->vanguard().grid(),
                                    simulator->vanguard().gridView(),
                                    simulator->vanguard().cartesianMapper(), 10.0);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-1;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first] , 1.496329839e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last ] , 1.504526940e7   , reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last] , 1.504526940e7   , reltol);
}

BOOST_AUTO_TEST_CASE(CapillaryInversion)
{
    // Test setup.
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;

    auto simulator = initSimulator<TypeTag>("equil_capillary.DATA");

    // Test the capillary inversion for oil-water.
    const int cell = 0;
    const double reltol = 1.0e-5;
    {
        const int phase = FluidSystem::waterPhaseIdx;
        const bool increasing = false;
        const std::vector<double> pc = { 10.0e5, 0.5e5, 0.4e5, 0.3e5, 0.2e5, 0.1e5, 0.099e5, 0.0e5, -10.0e5 };
        const std::vector<double> s = { 0.2, 0.2, 0.2, 0.466666666666, 0.733333333333, 1.0, 1.0, 1.0, 1.0 };
        BOOST_REQUIRE_EQUAL(pc.size(), s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Opm::EQUIL::satFromPc<FluidSystem>(*simulator->problem().materialLawManager(), phase, cell, pc[i], increasing);
            BOOST_CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }

    // Test the capillary inversion for gas-oil.
    {
        const int phase = FluidSystem::gasPhaseIdx;
        const bool increasing = true;
        const std::vector<double> pc = { 10.0e5, 0.6e5, 0.5e5, 0.4e5, 0.3e5, 0.2e5, 0.1e5, 0.0e5, -10.0e5 };
        const std::vector<double> s = { 0.8, 0.8, 0.8, 0.533333333333, 0.266666666666, 0.0, 0.0, 0.0, 0.0 };
        BOOST_REQUIRE_EQUAL(pc.size(), s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Opm::EQUIL::satFromPc<FluidSystem>(*simulator->problem().materialLawManager(), phase, cell, pc[i], increasing);
            BOOST_CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }

    // Test the capillary inversion for gas-water.
    {
        const int water = FluidSystem::waterPhaseIdx;
        const int gas = FluidSystem::gasPhaseIdx;
        const std::vector<double> pc = { 0.9e5, 0.8e5, 0.6e5, 0.4e5, 0.3e5 };
        const std::vector<double> s = { 0.2, 0.333333333333, 0.6, 0.866666666666, 1.0 };
        BOOST_REQUIRE_EQUAL(pc.size(), s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Opm::EQUIL::satFromSumOfPcs<FluidSystem>(*simulator->problem().materialLawManager(), water, gas, cell, pc[i]);
            BOOST_CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }
}

BOOST_AUTO_TEST_CASE(DeckWithCapillary)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator = initSimulator<TypeTag>("equil_capillary.DATA");
    auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState, 
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 10.0);

    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-4;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.469769063e7     , reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last ], 15452880.328284413, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx]  [last ], 15462880.328284413, reltol);

    const auto& sats = comp.saturation();
    std::vector<double> s[3];
    s[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.42190294373815257, 0.77800802072306474, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s[FluidSystem::oilPhaseIdx]   = { 0, 0, 0, 0.0073481611123183965, 0.79272270823081337, 0.8, 0.8, 0.8, 0.8, 0.57809705626184749, 0.22199197927693526, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s[FluidSystem::gasPhaseIdx]   = { 0.8, 0.8, 0.8, 0.79265183888768165, 0.0072772917691866562, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s[phase].size());
        for (size_t i = 0; i < s[phase].size(); ++i) {
            BOOST_CHECK_CLOSE(sats[phase][i], s[phase][i], reltol);
        }
    }
}

BOOST_AUTO_TEST_CASE(DeckWithCapillaryOverlap)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator = initSimulator<TypeTag>("equil_capillary_overlap.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState,
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-4;
    const double reltol_ecl = 100.0;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.48324e+07, reltol_ecl);  // eclipse
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.54801e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.49224e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.54901e+07, reltol_ecl);

    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first] , 14832467.14, reltol); // opm
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last ] , 15479883.47, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last ] , 15489883.47, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }

    std::vector<double> s_ecl[3]; // eclipse
    s_ecl[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22874042, 0.53397995, 0.78454906,  0.91542006, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_ecl[FluidSystem::oilPhaseIdx] =   { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.20039,     0.08458,    0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_ecl[FluidSystem::gasPhaseIdx] =   { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77125955, 0.46602005, 0.015063271, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0 };

    std::vector<double> s_opm[3]; // opm
    s_opm[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22892931226886132,  0.53406457830052489, 0.78457075254244724, 0.91539712466977541, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_opm[FluidSystem::oilPhaseIdx] = { 0,   0,   0,   0,   0,   0,   0,   0,                   0,                   0.20023624994125844,   0.084602875330224592, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_opm[FluidSystem::gasPhaseIdx] = { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77107068773113863, 0.46593542169947511, 0.015192997516294321, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0 };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
        }
    }
}

BOOST_AUTO_TEST_CASE(DeckWithLiveOil)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator = initSimulator<TypeTag>("equil_liveoil.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    // Initialize the fluid system
    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState,
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-4;
    const double reltol_ecl = 100.0;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.48324e+07, reltol_ecl);  // eclipse
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.54801e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.49224e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.54901e+07, reltol_ecl);

    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.483246714e7, reltol);  // opm
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.547991652e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.492246714e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.548991652e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    std::vector<double> s_ecl[3]; // eclipse
    s_ecl[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22898, 0.53422, 0.78470, 0.91531, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_ecl[FluidSystem::oilPhaseIdx] =   { 0,   0,   0,   0,   0,   0,   0,   0,       0,       0.20073, 0.08469, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_ecl[FluidSystem::gasPhaseIdx] =   { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77102, 0.46578, 0.01458, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector<double> s_opm[3]; // opm
    s_opm[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22916963446461344, 0.53430490523774521, 0.78471886612242092, 0.91528324362210933, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_opm[FluidSystem::oilPhaseIdx] = { 0,   0,   0,   0,   0,   0,   0,   0,            0,            0.20057438297017782,   0.084716756377890667, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_opm[FluidSystem::gasPhaseIdx] = { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77083036553538653, 0.46569509476225479, 0.014706750907401245,  0,       0, 0, 0, 0, 0, 0, 0, 0, 0 };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
        }
        std::cout << std::endl;
    }

    const auto& rs = comp.rs();
    const std::vector<double> rs_opm {74.61233568, 74.64905212, 74.68578656, 74.72253902, // opm
                                      74.75930951, 74.79609803, 74.83290459, 74.87519876,
                                      74.96925416, 75.09067512, 75.0,        75.0,
                                      75.0,        75.0,        75.0,        75.0,
                                      75.0,        75.0,        75.0,        75.0};
    const std::vector<double> rs_ecl {74.612228, 74.648956, 74.685707, 74.722473,  // eclipse
                                      74.759254, 74.796051, 74.832870, 74.875145,
                                      74.969231, 75.090706, 75.000000, 75.000000,
                                      75.000000, 75.000000, 75.000000, 75.000000,
                                      75.000000, 75.000000, 75.000000, 75.000000};
    for (size_t i = 0; i < rs_opm.size(); ++i) {
        //std::cout << std::setprecision(10) << rs[i] << '\n';
        BOOST_CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        BOOST_CHECK_CLOSE(rs[i], rs_ecl[i], reltol_ecl);
    }
}

BOOST_AUTO_TEST_CASE(DeckWithCO2STORE)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator1 = initSimulator<TypeTag>("equil_co2store_go.DATA");
    EquilFixture::Initializer comp_go(*simulator1->problem().materialLawManager(),
                                    simulator1->vanguard().eclState(),
                                    simulator1->vanguard().grid(),
                                    simulator1->vanguard().gridView(),
                                    simulator1->vanguard().cartesianMapper(), 9.80665);

    auto simulator2 = initSimulator<TypeTag>("equil_co2store_gw.DATA");
    EquilFixture::Initializer comp_gw(*simulator2->problem().materialLawManager(),
                                     simulator2->vanguard().eclState(),
                                     simulator2->vanguard().grid(),
                                     simulator2->vanguard().gridView(),
                                     simulator2->vanguard().cartesianMapper(), 9.80665);

    Opm::GridManager gm(simulator2->vanguard().eclState().getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    const double reltol = 1.0e-5;
    const auto& pressures_go = comp_go.press();
    BOOST_REQUIRE_EQUAL(pressures_go.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures_go[0].size()), grid.number_of_cells);

    const auto& pressures_gw = comp_gw.press();
    BOOST_REQUIRE_EQUAL(pressures_gw.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures_gw[0].size()), grid.number_of_cells);

    const auto& sats_go = comp_go.saturation();
    const auto& sats_gw = comp_gw.saturation();

    for (int i = 0; i < grid.number_of_cells; ++i) {
        BOOST_CHECK_CLOSE(pressures_go[FluidSystem::gasPhaseIdx][i],  pressures_gw[FluidSystem::gasPhaseIdx][i], reltol);
        BOOST_CHECK_CLOSE(pressures_go[FluidSystem::oilPhaseIdx][i], pressures_gw[FluidSystem::waterPhaseIdx][i], reltol);

        BOOST_CHECK_CLOSE(sats_go[FluidSystem::gasPhaseIdx][i], sats_gw[FluidSystem::gasPhaseIdx][i], reltol);
        BOOST_CHECK_CLOSE(sats_go[FluidSystem::oilPhaseIdx][i], sats_gw[FluidSystem::waterPhaseIdx][i], reltol);
    }
}


BOOST_AUTO_TEST_CASE(DeckWithWetGas)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator = initSimulator<TypeTag>("equil_wetgas.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState,
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-1;
    const double reltol_ecl = 100.0;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.48215e+07, reltol_ecl);  // eclipse
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.54801e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.49115e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.54901e+07, reltol_ecl);

    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.482150311e7, reltol);  // opm
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.547988347e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.491150311e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.548988347e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    std::vector<double> s_ecl[3]; // eclipse
    s_ecl[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.24285614, 0.53869015, 0.78454906,  0.91542006, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_ecl[FluidSystem::oilPhaseIdx] =   { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.18311,     0.08458,    0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_ecl[FluidSystem::gasPhaseIdx] =   { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75714386, 0.46130988, 0.032345835, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector<double> s_opm[3]; // opm
    s_opm[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.24310545, 0.5388, 0.78458,    0.91540, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_opm[FluidSystem::oilPhaseIdx] =   { 0,   0,   0,   0,   0,   0,   0,   0,          0,      0.18288667, 0.0846,  0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_opm[FluidSystem::gasPhaseIdx] =   { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75689455, 0.4612, 0.03253333, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0 };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], 100.*reltol);
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
        }
        std::cout << std::endl;
    }

    const auto& rv = comp.rv();
    const std::vector<double> rv_opm { // opm
        2.4884509e-4, 2.4910378e-4, 2.4936267e-4, 2.4962174e-4,
        2.4988100e-4, 2.5014044e-4, 2.5040008e-4, 2.5065990e-4,
        2.5091992e-4, 2.5118012e-4, 2.5223082e-4, 2.5105e-4,
        2.5105e-4,    2.5105e-4,    2.5105e-4,    2.5105e-4,
        2.5105e-4,    2.5105e-4,    2.5105e-4,    2.5105e-4};

    const std::vector<double> rv_ecl {  // eclipse
        0.24884584E-03,   0.24910446E-03,   0.24936325E-03,   0.24962222E-03,
        0.24988138E-03,   0.25014076E-03,   0.25040031E-03,   0.25066003E-03,
        0.25091995E-03,   0.25118008E-03,   0.25223137E-03,   0.25104999E-03,
        0.25104999E-03,   0.25104999E-03,   0.25104999E-03,   0.25104999E-03,
        0.25104999E-03,   0.25104999E-03,   0.25104999E-03,   0.25104999E-03};

    for (size_t i = 0; i < rv_opm.size(); ++i) {
        BOOST_CHECK_CLOSE(rv[i], rv_opm[i], reltol);
        BOOST_CHECK_CLOSE(rv[i], rv_ecl[i], reltol_ecl);
    }
}

BOOST_AUTO_TEST_CASE(DeckWithHumidWetGas)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    FluidSystem::setEnableVaporizedWater(true);
    auto simulator = initSimulator<TypeTag>("equil_humidwetgas.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState,
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    const double reltol = 1.0e-1;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.480599988e7, reltol);  
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.549297524e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.489599988e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.550297524e7, reltol);

    const auto& sats = comp.saturation();
    std::vector<double> s_opm[3]; 
    s_opm[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.35838026, 0.64069098, 0.9154626,    1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_opm[FluidSystem::oilPhaseIdx] =   { 0,   0,   0,   0,   0,   0,   0,   0,          0,      0.02738364, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_opm[FluidSystem::gasPhaseIdx] =   { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.64161973, 0.359309012, 0.057153701, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0 };
    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], 100.*reltol);
        }
        std::cout << std::endl;
    }

    const auto& rv = comp.rv();
    const std::vector<double> rv_opm { 
        0.00024837999651755729, 0.00024869285236692635, 0.00024900604366769004, 0.00024931957094322978,
        0.00024963343471801471, 0.00024994763551760586, 0.00025026217386865733, 0.00025057705029892072,
        0.00025089226533724643, 0.00025120780158539152, 0.00025105, 0.00025105,
        0.00025105, 0.00025105, 0.00025105, 0.00025105,
        0.00025105, 0.00025105, 0.00025105, 0.00025105};

    const auto& rvw = comp.rvw();
    const std::vector<double> rvw_opm {  
        0.00024837999651755729, 0.00024869285236692635, 0.00024900604366769004, 0.00024931957094322978,
        0.00024963343471801471, 0.00024994763551760586, 0.00025026217386865733, 0.00025057705029892072,
        0.00025089226533724643, 0.00025120780158539152, 0.00025236969680655122, 0.00025384953117447344,
        0.00025532939474124625, 0.00025680928750801825, 0.00025828920947593858, 0.00025976916064615645,
        0.00026124914101982041, 0.00026272915059807997, 0.0002642091893820838, 0.00026568925737298143};

    for (size_t i = 0; i < rv_opm.size(); ++i) {
        BOOST_CHECK_CLOSE(rv[i], rv_opm[i], reltol);
        BOOST_CHECK_CLOSE(rvw[i], rvw_opm[i], reltol);
    }
}

BOOST_AUTO_TEST_CASE(DeckWithRSVDAndRVVD)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator = initSimulator<TypeTag>("equil_rsvd_and_rvvd.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState,
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-4;
    const double reltol_ecl = 100.0;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.48350e+07, reltol_ecl);  // eclipse
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.54794e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.49250e+07, reltol_ecl);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.54894e+07, reltol_ecl);

    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.483499660e7, reltol);  // opm
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.547924516e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.492499660e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.548924516e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    std::vector<double> s_ecl[3]; // eclipse
    s_ecl[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22206347, 0.52871972, 0.78150368,  0.91819441,  1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_ecl[FluidSystem::oilPhaseIdx] =   { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.19656529,  0.081805572, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_ecl[FluidSystem::gasPhaseIdx] =   { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77793652, 0.47128031, 0.021931054, 0,           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    std::vector<double> s_opm[3]; // opm
    s_opm[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2223045711692897, 0.52882298575945874, 0.78152142505479982, 0.91816512259416283, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_opm[FluidSystem::oilPhaseIdx]   = { 0,   0,   0,   0,   0,   0,   0,   0,          0, 0.19637607881498206, 0.08183487740583717, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_opm[FluidSystem::gasPhaseIdx]   = { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7776954288307103, 0.47117701424054126, 0.02210249613021811,    0,          0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
            BOOST_CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
        }
        std::cout << std::endl;
    }

    const auto& rs = comp.rs();
    const std::vector<double> rs_opm { // opm
        74.62498302, 74.65959041, 74.69438035, 74.72935336,
        74.76450995, 74.79985061, 74.83537588, 74.87527065,
        74.96863769, 75.08891765, 52.5,        57.5,
        62.5,        67.5,        72.5,        76.45954841,
        76.70621045, 76.95287736, 77.19954913, 77.44622578};

    const std::vector<double> rs_ecl {  // eclipse
        74.625114, 74.659706, 74.694481, 74.729439,
        74.764580, 74.799904, 74.835419, 74.875252,
        74.968628, 75.088951, 52.500000, 57.500000,
        62.500000, 67.500000, 72.500000, 76.168388,
        76.349953, 76.531532, 76.713142, 76.894775,};

    const auto& rv = comp.rv();
    const std::vector<double> rv_opm { // opm
        2.50e-6, 7.50e-6,       1.25e-5,       1.75e-5,
        2.25e-5, 2.75e-5,       3.25e-5,       3.75e-5,
        4.25e-5, 2.51158386e-4, 2.52203372e-4, 5.75e-5,
        6.25e-5, 6.75e-5,       7.25e-5,       7.75e-5,
        8.25e-5, 8.75e-5,       9.25e-5,       9.75e-5};

    const std::vector<double> rv_ecl {  // eclipse
        0.24999999E-05, 0.74999998E-05, 0.12500000E-04, 0.17500000E-04,
        0.22500000E-04, 0.27500000E-04, 0.32500000E-04, 0.37500002E-04,
        0.42500000E-04, 0.25115837E-03, 0.25220393E-03, 0.57500001E-04,
        0.62500003E-04, 0.67499997E-04, 0.72499999E-04, 0.77500001E-04,
        0.82500002E-04, 0.87499997E-04, 0.92499999E-04, 0.97500000E-04};

    for (size_t i = 0; i < rv_opm.size(); ++i) {
        //std::cout << std::setprecision(10) << rs[i] << '\n';
        BOOST_CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        BOOST_CHECK_CLOSE(rs[i], rs_ecl[i], reltol_ecl);
        BOOST_CHECK_CLOSE(rv[i], rv_opm[i], reltol);
        BOOST_CHECK_CLOSE(rv[i], rv_ecl[i], reltol_ecl);
    }
}

BOOST_AUTO_TEST_CASE(DeckWithPBVDAndPDVD)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    auto simulator = initSimulator<TypeTag>("equil_pbvd_and_pdvd.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState,
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 5.0e-4;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 14821552.0, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  15479828.0, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 14911552.0, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  15489828.0, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }

    std::vector<double> s_opm[3]; // opm
    s_opm[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.24257337312592703, 0.53834824764362788, 0.7844998821510003, 0.9152832369551807, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_opm[FluidSystem::oilPhaseIdx] =   { 0,   0,   0,   0,   0,   0,   0,   0,          0, 0.18185970596719522, 0.084716763044819343, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_opm[FluidSystem::gasPhaseIdx] =   { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75742662687407303, 0.46165175235637212, 0.033640411881804465,    0,          0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
        }
    }

    const auto& rs = comp.rs();
    const std::vector<double> rs_opm { // opm
        74.55776480956456,
        74.6008507125663,
        74.6439680789467,
        74.68711693934459,
        74.73029732443825,
        74.77350926494491,
        74.81675279162118,
        74.86802321984302,
        74.96677993174352,
        75.09034523640406,
        75, 75, 75,75,75, 75, 75, 75, 75, 75 };

    const auto& rv = comp.rv();
    const std::vector<double> rv_opm {
        0.0002488465888573874,
        0.0002491051042753978,
        0.0002493638084736803,
        0.0002496227016360676,
        0.0002498817839466295,
        0.00025,
        0.00025,
        0.00025,
        0.00025,
        0.000251180039180951,
        0.0002522295187440788,
        0.0002275000000000001,
        0.0002125,
        0.0001975,
        0.0001825,
        0.0001675,
        0.0001525,
        0.0001375,
        0.0001225,
        0.0001075};

    for (size_t i = 0; i < rv_opm.size(); ++i) {
        BOOST_CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        BOOST_CHECK_CLOSE(rv[i], rv_opm[i], reltol);
    }
}

BOOST_AUTO_TEST_CASE(DeckWithRSVDAndRVVDAndRVWVD)
{
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    using FluidSystem = Opm::GetPropType<TypeTag, Opm::Properties::FluidSystem>;
    FluidSystem::setEnableVaporizedWater(true);
    auto simulator = initSimulator<TypeTag>("equil_rsvd_and_rvvd_and_rvwvd.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    EquilFixture::Initializer comp(*simulator->problem().materialLawManager(),
                                   eclipseState,
                                   simulator->vanguard().grid(),
                                   simulator->vanguard().gridView(),
                                   simulator->vanguard().cartesianMapper(), 9.80665);
    const auto& pressures = comp.press();
    BOOST_REQUIRE_EQUAL(pressures.size(), 3U);
    BOOST_REQUIRE_EQUAL(int(pressures[0].size()), grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    const double reltol = 1.0e-4;
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][first], 1.483359963e7, reltol);  
    BOOST_CHECK_CLOSE(pressures[FluidSystem::waterPhaseIdx][last],  1.549297524e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][first], 1.492359963e7, reltol);
    BOOST_CHECK_CLOSE(pressures[FluidSystem::oilPhaseIdx][last],  1.550297524e7, reltol);

    const auto& sats = comp.saturation();
    std::vector<double> s_opm[3]; 
    s_opm[FluidSystem::waterPhaseIdx] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.32527877578987319, 0.62976875867666171, 0.918795223850500588, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    s_opm[FluidSystem::oilPhaseIdx]   = { 0,   0,   0,   0,   0,   0,   0,   0,          0, 0.054786199472198836, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    s_opm[FluidSystem::gasPhaseIdx]   = { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.674721224210102681, 0.37023124132333829, 0.026418562022795279,    0,          0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for (int phase = 0; phase < 3; ++phase) {
        BOOST_REQUIRE_EQUAL(sats[phase].size(), s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            BOOST_CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
        }
        std::cout << std::endl;
    }

    const auto& rs = comp.rs();
    const std::vector<double> rs_opm { // opm
        74.617998198796087,74.652774471604374, 74.687905898686935, 74.723393674854691,
        74.759238999357947, 74.795443075905553, 74.832007112684167, 74.892422092838459, 
        74.986801564438935, 75.088917653469338, 52.5, 57.5, 
        62.5, 67.5, 72.5, 76.528193441026076,
        76.774856836636729, 77.021525099679991, 77.268198230347295, 77.514876228830232};

    const auto& rv = comp.rv();
    const std::vector<double> rv_opm { 
        2.5000000000000002e-06, 7.5000000000000002e-06, 1.2500000000000001e-05, 1.7500000000000002e-05,
        2.2500000000000001e-05, 2.7500000000000004e-05, 3.2500000000000004e-05, 3.7500000000000003e-05,
        4.2500000000000003e-05, 0.00025116322680309166, 5.2500000000000002e-05, 5.7500000000000002e-05,
        6.2500000000000001e-05, 6.7500000000000001e-05, 7.25e-05, 7.75e-05,
        8.25e-05, 8.7500000000000013e-05, 9.2500000000000012e-05, 9.7499999999999998e-05};

    const auto& rvw = comp.rvw();
    const std::vector<double> rvw_opm {  
        0.00024920798919277656, 0.00024941664682962629, 0.00024962743539212165, 0.00024984036204912818,
        0.00025005543399614773, 0.00025027265845543336, 0.000250492042676105, 0.00025071359393426718,
        0.00025093731953312241, 0.00025116322680309166, 0.00025236969680655122, 0.00025384953117447344,
        0.00025532939474124625, 0.00025680928750801825, 0.00025828920947593858, 0.00025976916064615645,
        0.00026124914101982041, 0.00026272915059807997, 0.0002642091893820838, 0.00026568925737298143};

    for (size_t i = 0; i < rv_opm.size(); ++i) {
        BOOST_CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        BOOST_CHECK_CLOSE(rv[i], rv_opm[i], reltol);
        BOOST_CHECK_CLOSE(rvw[i], rvw_opm[i], reltol);
    }
}

BOOST_AUTO_TEST_CASE(DeckWithSwatinit)
{
#if 0
    using TypeTag = Opm::Properties::TTag::TestEquilTypeTag;
    auto simulator = initSimulator<TypeTag>("equil_capillary_swatinit.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    // Create material law manager.
    std::vector<int> compressedToCartesianIdx
        = Opm::compressedToCartesian(grid.number_of_cells, grid.global_cell);
    MaterialLawManager materialLawManager = MaterialLawManager();
    materialLawManager.initFromDeck(deck, eclipseState, compressedToCartesianIdx);

    MaterialLawManager materialLawManagerScaled = MaterialLawManager();
    materialLawManagerScaled.initFromDeck(deck, eclipseState, compressedToCartesianIdx);

    // reference saturations
    const std::vector<double> s[3]{
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.42528761746004229, 0.77462669821009045, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0.014813991154779993, 0.78525420807446045, 0.8, 0.8, 0.8, 0.8, 0.57471238253995771, 0.22537330178990955, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.78518600884522005, 0.014745791925539575, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    // sw in cell 1-5 is forced to be 0.2 since swl=0.2
    // sw in cell 13 and 14 is forced to be swu=1 since P_oil - P_wat < 0.
    const std::vector<double> swatinit[3]{
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0.014813991154779993, 0.78525420807446045, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.78518600884522005, 0.014745791925539575, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    // Adjust oil pressure according to gas saturation and cap pressure
    typedef Opm::SimpleModularFluidState<double,
                                         /*numPhases=*/3,
                                         /*numComponents=*/3,
                                         FluidSystem,
                                         /*storePressure=*/false,
                                         /*storeTemperature=*/false,
                                         /*storeComposition=*/false,
                                         /*storeFugacity=*/false,
                                         /*storeSaturation=*/true,
                                         /*storeDensity=*/false,
                                         /*storeViscosity=*/false,
                                         /*storeEnthalpy=*/false> SatOnlyFluidState;

    SatOnlyFluidState fluidState;
    typedef MaterialLawManager::MaterialLaw MaterialLaw;

    // Initialize the fluid system
    FluidSystem::initFromDeck(deck, eclipseState);

    // reference pcs
    const int numCells = Opm::UgGridHelpers::numCells(grid);
    std::vector<double> pc_original(numCells * FluidSystem::numPhases);
    for (int c = 0; c < numCells; ++c) {
        std::vector<double> pc = {0,0,0};
        double sw = s[FluidSystem::waterPhaseIdx][c];
        double so = s[FluidSystem::oilPhaseIdx][c];
        double sg = s[FluidSystem::gasPhaseIdx][c];
        fluidState.setSaturation(FluidSystem::waterPhaseIdx, sw);
        fluidState.setSaturation(FluidSystem::oilPhaseIdx, so);
        fluidState.setSaturation(FluidSystem::gasPhaseIdx, sg);
        const auto& matParams = materialLawManager.materialLawParams(c);
        MaterialLaw::capillaryPressures(pc, matParams, fluidState);
        pc_original[3*c + 0] = pc[FluidSystem::oilPhaseIdx] - pc[FluidSystem::waterPhaseIdx];
        pc_original[3*c + 1] = 0.0;
        pc_original[3*c + 2] = pc[FluidSystem::oilPhaseIdx] + pc[FluidSystem::gasPhaseIdx];
    }

    std::vector<double> pc_scaled_truth = pc_original;

    // modify pcow for cell 1 - 12 (where sw is changed due to swatinit)
    // for the reference scaled pc.
    pc_scaled_truth[3*0 + 0] = 150031.3;
    pc_scaled_truth[3*1 + 0] = 136815.6;
    pc_scaled_truth[3*2 + 0] = 123612.7;
    pc_scaled_truth[3*3 + 0] = 110422.7;
    pc_scaled_truth[3*4 + 0] = 97245.4;
    pc_scaled_truth[3*5 + 0] = 84081;
    pc_scaled_truth[3*6 + 0] = 70929;
    pc_scaled_truth[3*7 + 0] = 57791;
    pc_scaled_truth[3*8 + 0] = 44665;
    pc_scaled_truth[3*9 + 0] = 31552;
    pc_scaled_truth[3*10 + 0] = 18451.5;
    pc_scaled_truth[3*11 + 0] =  5364.1;

    // compute the initial state
    // apply swatinit
    Opm::EQUIL::DeckDependent::InitialStateComputer<TypeTag> compScaled(materialLawManagerScaled, eclipseState, simulator->vanguard().gridView(), 9.81, true);
    // don't apply swatinit
    Opm::EQUIL::DeckDependent::InitialStateComputer<TypeTag> compUnscaled(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().gridView(), 9.81, false);

    // compute pc
    std::vector<double> pc_scaled(numCells * FluidSystem::numPhases);
    for (int c = 0; c < numCells; ++c) {
        std::vector<double> pc = {0,0,0};
        double sw = compScaled.saturation().data()[FluidSystem::waterPhaseIdx][c];
        double so = compScaled.saturation().data()[FluidSystem::oilPhaseIdx][c];
        double sg = compScaled.saturation().data()[FluidSystem::gasPhaseIdx][c];

        fluidState.setSaturation(FluidSystem::waterPhaseIdx, sw);
        fluidState.setSaturation(FluidSystem::oilPhaseIdx, so);
        fluidState.setSaturation(FluidSystem::gasPhaseIdx, sg);
        const auto& matParams = materialLawManagerScaled.materialLawParams(c);
        MaterialLaw::capillaryPressures(pc, matParams, fluidState);
        pc_scaled[3*c + 0] = pc[FluidSystem::oilPhaseIdx] - pc[FluidSystem::waterPhaseIdx];
        pc_scaled[3*c + 1] = 0.0;
        pc_scaled[3*c + 2] = pc[FluidSystem::oilPhaseIdx] + pc[FluidSystem::gasPhaseIdx];
    }
    std::vector<double> pc_unscaled(numCells * FluidSystem::numPhases);
    for (int c = 0; c < numCells; ++c) {
        std::vector<double> pc = {0,0,0};
        double sw = compUnscaled.saturation().data()[FluidSystem::waterPhaseIdx][c];
        double so = compUnscaled.saturation().data()[FluidSystem::oilPhaseIdx][c];
        double sg = compUnscaled.saturation().data()[FluidSystem::gasPhaseIdx][c];

        fluidState.setSaturation(FluidSystem::waterPhaseIdx, sw);
        fluidState.setSaturation(FluidSystem::oilPhaseIdx, so);
        fluidState.setSaturation(FluidSystem::gasPhaseIdx, sg);

        const auto& matParams = materialLawManager.materialLawParams(c);
        MaterialLaw::capillaryPressures(pc, matParams, fluidState);
        pc_unscaled[3*c + 0] = pc[FluidSystem::oilPhaseIdx] - pc[FluidSystem::waterPhaseIdx];
        pc_unscaled[3*c + 1] = 0.0;
        pc_unscaled[3*c + 2] = pc[FluidSystem::oilPhaseIdx] + pc[FluidSystem::gasPhaseIdx];
    }

    // test
    const double reltol = 1.0e-1;
    for (int phase = 0; phase < 3; ++phase) {
        for (size_t i = 0; i < 20; ++i) {
            BOOST_CHECK_CLOSE( pc_original[3*i + phase ], pc_unscaled[3*i + phase ], reltol);
            BOOST_CHECK_CLOSE( pc_scaled_truth[3*i + phase], pc_scaled[3*i + phase ], reltol);
        }
    }

    for (int phase = 0; phase < 3; ++phase) {
        for (size_t i = 0; i < 20; ++i) {
            BOOST_CHECK_CLOSE(compUnscaled.saturation()[phase][i], s[phase][i], reltol);
            BOOST_CHECK_CLOSE(compScaled.saturation()[phase][i], swatinit[phase][i], reltol);
        }
    }
#endif
}
