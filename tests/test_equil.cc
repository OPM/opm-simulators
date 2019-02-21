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

#include <ebos/equil/equilibrationhelpers.hh>
#include <ebos/eclproblem.hh>
#include <ewoms/common/start.hh>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/GridManager.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <string.h>

#define CHECK(value, expected)             \
    {                                      \
        if ((value) != (expected))         \
            std::abort();                  \
    }

#define CHECK_CLOSE(value, expected, reltol)                            \
    {                                                                   \
        if (std::fabs((expected) - (value)) > 1e-14 &&                  \
            std::fabs(((expected) - (value))/((expected) + (value))) > reltol) \
            std::abort();                                               \
    }

#define REQUIRE(cond)                      \
    {                                      \
        if (!(cond))                       \
            std::abort();                  \
    }

BEGIN_PROPERTIES

NEW_TYPE_TAG(TestEquilTypeTag, INHERITS_FROM(BlackOilModel, EclBaseProblem));

END_PROPERTIES

template <class TypeTag>
std::unique_ptr<typename GET_PROP_TYPE(TypeTag, Simulator)>
initSimulator(const char *filename)
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    std::string filenameArg = "--ecl-deck-file-name=";
    filenameArg += filename;

    const char* argv[] = {
        "test_equil",
        filenameArg.c_str()
    };

    Ewoms::setupParameters_<TypeTag>(/*argc=*/sizeof(argv)/sizeof(argv[0]), argv, /*registerParams=*/false);

    return std::unique_ptr<Simulator>(new Simulator);
}

template <class TypeTag>
static void initDefaultFluidSystem()
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

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
    gasPvt->setApproach(Opm::GasPvtMultiplexer<double>::DryGasPvt);
    auto& dryGasPvt = gasPvt->getRealPvt<Opm::GasPvtMultiplexer<double>::DryGasPvt>();
    dryGasPvt.setNumRegions(/*numPvtRegion=*/1);
    dryGasPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
    dryGasPvt.setGasFormationVolumeFactor(/*regionIdx=*/0, Bg);
    dryGasPvt.setGasViscosity(/*regionIdx=*/0, mug);

    auto oilPvt = std::make_shared<Opm::OilPvtMultiplexer<double>>();
    oilPvt->setApproach(Opm::OilPvtMultiplexer<double>::DeadOilPvt);
    auto& deadOilPvt = oilPvt->getRealPvt<Opm::OilPvtMultiplexer<double>::DeadOilPvt>();
    deadOilPvt.setNumRegions(/*numPvtRegion=*/1);
    deadOilPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
    deadOilPvt.setInverseOilFormationVolumeFactor(/*regionIdx=*/0, Bo);
    deadOilPvt.setOilViscosity(/*regionIdx=*/0, muo);

    auto waterPvt = std::make_shared<Opm::WaterPvtMultiplexer<double>>();
    waterPvt->setApproach(Opm::WaterPvtMultiplexer<double>::ConstantCompressibilityWaterPvt);
    auto& ccWaterPvt = waterPvt->getRealPvt<Opm::WaterPvtMultiplexer<double>::ConstantCompressibilityWaterPvt>();
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
    using namespace Opm;

    DeckItem dd( "datdep", double() );
    dd.push_back( datd  );
    Opm::Dimension dd_dim( "dddim", 1 );
    dd.push_backDimension( dd_dim, dd_dim );

    DeckItem dp( "datps", double() );
    dp.push_back( datp );
    Opm::Dimension dp_dim( "dpdim", 1 );
    dp.push_backDimension( dp_dim, dp_dim );

    DeckItem zw( "zwoc", double() );
    zw.push_back( zwoc );
    Opm::Dimension zw_dim( "zwdim", 1 );
    zw.push_backDimension( zw_dim, zw_dim );

    DeckItem pcow( "pcow", double() );
    pcow.push_back( pcow_woc );
    Opm::Dimension pcow_dim( "pcowdim", 1 );
    pcow.push_backDimension( pcow_dim, pcow_dim );

    DeckItem zg( "zgoc", double() );
    zg.push_back( zgoc );
    Opm::Dimension zg_dim( "zgdim", 1 );
    zg.push_backDimension( zg_dim, zg_dim );

    DeckItem pcgo( "pcgo", double() );
    pcgo.push_back( pcgo_goc );
    Opm::Dimension pcgo_dim( "pcgodim", 1 );
    pcgo.push_backDimension( pcgo_dim, pcgo_dim );

    DeckItem i1( "i1", int() );
    DeckItem i2( "i2", int() );
    DeckItem i3( "i3", int() );
    i1.push_back( 0 );
    i2.push_back( 0 );
    i3.push_back( 0 );

    DeckRecord rec;
    rec.addItem( std::move( dd ) );
    rec.addItem( std::move( dp ) );
    rec.addItem( std::move( zw ) );
    rec.addItem( std::move( pcow ) );
    rec.addItem( std::move( zg ) );
    rec.addItem( std::move( pcgo ) );
    rec.addItem( std::move( i1 ) );
    rec.addItem( std::move( i2 ) );
    rec.addItem( std::move( i3 ) );

    return EquilRecord( rec );
}

void test_PhasePressure();
void test_PhasePressure()
{
    typedef std::vector<double> PVal;
    typedef std::vector<PVal> PPress;

    auto record = mkEquilRecord( 0, 1e5, 5, 0, 0, 0 );

    typedef TTAG(TestEquilTypeTag) TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    auto simulator = initSimulator<TypeTag>("equil_base.DATA");
    initDefaultFluidSystem<TypeTag>();

    Ewoms::EQUIL::EquilReg
        region(record,
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0);

    std::vector<int> cells(simulator->vanguard().grid().size(0));
    std::iota(cells.begin(), cells.end(), 0);

    const double grav   = 10;
    const PPress ppress = Ewoms::EQUIL::phasePressures<FluidSystem>(simulator->vanguard().grid(), region, cells, grav);

    const int first = 0, last = simulator->vanguard().grid().size(0) - 1;
    const double reltol = 1.0e-8;
    CHECK_CLOSE(ppress[0][first] ,  90e3   , reltol);
    CHECK_CLOSE(ppress[0][last ] , 180e3   , reltol);
    CHECK_CLOSE(ppress[1][first] , 103.5e3 , reltol);
    CHECK_CLOSE(ppress[1][last ] , 166.5e3 , reltol);
}

void test_CellSubset();
void test_CellSubset()
{
    typedef std::vector<double> PVal;
    typedef std::vector<PVal> PPress;

    typedef TTAG(TestEquilTypeTag) TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    auto simulator = initSimulator<TypeTag>("equil_base.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());
    initDefaultFluidSystem<TypeTag>();

    Opm::EquilRecord record[] = { mkEquilRecord( 0, 1e5, 2.5, -0.075e5, 0, 0 ),
                                  mkEquilRecord( 5, 1.35e5, 7.5, -0.225e5, 5, 0 ) };

    Ewoms::EQUIL::EquilReg region[] =
    {
        Ewoms::EQUIL::EquilReg(record[0],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0)
        ,
        Ewoms::EQUIL::EquilReg(record[0],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0)
        ,
        Ewoms::EQUIL::EquilReg(record[1],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0)
        ,
        Ewoms::EQUIL::EquilReg(record[1],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
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

    PPress ppress(2, PVal(simulator->vanguard().grid().size(0), 0));
    for (std::vector< std::vector<int> >::const_iterator
             r = cells.begin(), e = cells.end();
         r != e; ++r)
    {
        const int    rno  = int(r - cells.begin());
        const double grav = 10;
        const PPress p    =
            Ewoms::EQUIL::phasePressures<FluidSystem>(simulator->vanguard().grid(), region[rno], *r, grav);

        PVal::size_type i = 0;
        for (std::vector<int>::const_iterator
                 c = r->begin(), ce = r->end();
             c != ce; ++c, ++i)
        {
            assert (i < p[0].size());

            ppress[0][*c] = p[0][i];
            ppress[1][*c] = p[1][i];
        }
    }

    const int first = 0, last = simulator->vanguard().grid().size(0) - 1;
    const double reltol = 1.0e-8;
    CHECK_CLOSE(ppress[0][first] , 105e3   , reltol);
    CHECK_CLOSE(ppress[0][last ] , 195e3   , reltol);
    CHECK_CLOSE(ppress[1][first] , 103.5e3 , reltol);
    CHECK_CLOSE(ppress[1][last ] , 166.5e3 , reltol);
}

void test_RegMapping();
void test_RegMapping()
{
    typedef std::vector<double> PVal;
    typedef std::vector<PVal> PPress;

    Opm::EquilRecord record[] = { mkEquilRecord( 0, 1e5, 2.5, -0.075e5, 0, 0 ),
                                  mkEquilRecord( 5, 1.35e5, 7.5, -0.225e5, 5, 0 ) };

    typedef TTAG(TestEquilTypeTag) TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    auto simulator = initSimulator<TypeTag>("equil_base.DATA");
    initDefaultFluidSystem<TypeTag>();

    Ewoms::EQUIL::EquilReg region[] =
    {
        Ewoms::EQUIL::EquilReg(record[0],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0)
        ,
        Ewoms::EQUIL::EquilReg(record[0],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0)
        ,
        Ewoms::EQUIL::EquilReg(record[1],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0)
        ,
        Ewoms::EQUIL::EquilReg(record[1],
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        std::make_shared<Ewoms::EQUIL::Miscibility::NoMixing>(),
        0)
        };

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

    Ewoms::RegionMapping<> eqlmap(eqlnum);

    PPress ppress(2, PVal(simulator->vanguard().grid().size(0), 0));
    for (const auto& r : eqlmap.activeRegions()) {
        const auto& rng = eqlmap.cells(r);

        const int    rno  = r;
        const double grav = 10;
        const PPress p    =
            Ewoms::EQUIL::phasePressures<FluidSystem>(simulator->vanguard().grid(), region[rno], rng, grav);

        PVal::size_type i = 0;
        for (const auto& c : rng) {
            assert (i < p[0].size());

            ppress[0][c] = p[0][i];
            ppress[1][c] = p[1][i];

            ++i;
        }
    }

    const int first = 0, last = simulator->vanguard().grid().size(0) - 1;
    const double reltol = 1.0e-8;
    CHECK_CLOSE(ppress[0][first] , 105e3   , reltol);
    CHECK_CLOSE(ppress[0][last ] , 195e3   , reltol);
    CHECK_CLOSE(ppress[1][first] , 103.5e3 , reltol);
    CHECK_CLOSE(ppress[1][last ] , 166.5e3 , reltol);
}

void test_DeckAllDead();
void test_DeckAllDead()
{
    typedef TTAG(TestEquilTypeTag) TypeTag;
    auto simulator = initSimulator<TypeTag>("equil_deadfluids.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> comp(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 10.0);
    const auto& pressures = comp.press();
    REQUIRE(pressures.size() == 3);
    REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-3;
    CHECK_CLOSE(pressures[0][first] , 1.496329839e7   , reltol);
    CHECK_CLOSE(pressures[0][last ] , 1.504526940e7   , reltol);
    CHECK_CLOSE(pressures[1][last] , 1.504526940e7   , reltol);
}

void test_CapillaryInversion();
void test_CapillaryInversion()
{
    // Test setup.
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager MaterialLawManager;
    auto simulator = initSimulator<TypeTag>("equil_capillary.DATA");

    // Test the capillary inversion for oil-water.
    const int cell = 0;
    const double reltol = 1.0e-7;
    {
        const int phase = 0;
        const bool increasing = false;
        const std::vector<double> pc = { 10.0e5, 0.5e5, 0.4e5, 0.3e5, 0.2e5, 0.1e5, 0.099e5, 0.0e5, -10.0e5 };
        const std::vector<double> s = { 0.2, 0.2, 0.2, 0.466666666666, 0.733333333333, 1.0, 1.0, 1.0, 1.0 };
        REQUIRE(pc.size() == s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Ewoms::EQUIL::satFromPc<FluidSystem, MaterialLaw, MaterialLawManager>(*simulator->problem().materialLawManager(), phase, cell, pc[i], increasing);
            CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }

    // Test the capillary inversion for gas-oil.
    {
        const int phase = 2;
        const bool increasing = true;
        const std::vector<double> pc = { 10.0e5, 0.6e5, 0.5e5, 0.4e5, 0.3e5, 0.2e5, 0.1e5, 0.0e5, -10.0e5 };
        const std::vector<double> s = { 0.8, 0.8, 0.8, 0.533333333333, 0.266666666666, 0.0, 0.0, 0.0, 0.0 };
        REQUIRE(pc.size() == s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Ewoms::EQUIL::satFromPc<FluidSystem, MaterialLaw, MaterialLawManager>(*simulator->problem().materialLawManager(), phase, cell, pc[i], increasing);
            CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }

    // Test the capillary inversion for gas-water.
    {
        const int water = 0;
        const int gas = 2;
        const std::vector<double> pc = { 0.9e5, 0.8e5, 0.6e5, 0.4e5, 0.3e5 };
        const std::vector<double> s = { 0.2, 0.333333333333, 0.6, 0.866666666666, 1.0 };
        REQUIRE(pc.size() == s.size());
        for (size_t i = 0; i < pc.size(); ++i) {
            const double s_computed = Ewoms::EQUIL::satFromSumOfPcs<FluidSystem, MaterialLaw, MaterialLawManager>(*simulator->problem().materialLawManager(), water, gas, cell, pc[i]);
            CHECK_CLOSE(s_computed, s[i], reltol);
        }
    }
}

void test_DeckWithCapillary();
void test_DeckWithCapillary()
{
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
    auto simulator = initSimulator<TypeTag>("equil_capillary.DATA");
    auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> comp(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 10.0);

    const auto& pressures = comp.press();
    REQUIRE(pressures.size() == 3);
    REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    CHECK_CLOSE(pressures[0][first] , 1.469769063e7   , reltol);
    CHECK_CLOSE(pressures[0][last ] , 15452880.328284413   , reltol);
    CHECK_CLOSE(pressures[1][last] , 15462880.328284413   , reltol);

    const auto& sats = comp.saturation();
    const std::vector<double> s[3]{
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.42190294373815257, 0.77800802072306474, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0, 0, 0, 0.0073481611123183965, 0.79272270823081337, 0.8, 0.8, 0.8, 0.8, 0.57809705626184749, 0.22199197927693526, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.79265183888768165, 0.0072772917691866562, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        REQUIRE(sats[phase].size() == s[phase].size());
        for (size_t i = 0; i < s[phase].size(); ++i) {
            CHECK_CLOSE(sats[phase][i], s[phase][i], reltol);
        }
    }
}

void test_DeckWithCapillaryOverlap();
void test_DeckWithCapillaryOverlap()
{
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
    auto simulator = initSimulator<TypeTag>("equil_capillary_overlap.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> comp(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 9.80665);
    const auto& pressures = comp.press();
    REQUIRE(pressures.size() == 3);
    REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    const double reltol_ecl = 1.0;
    CHECK_CLOSE(pressures[0][first], 1.48324e+07, reltol_ecl);  // eclipse
    CHECK_CLOSE(pressures[0][last],  1.54801e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][first], 1.49224e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][last],  1.54901e+07, reltol_ecl);

    CHECK_CLOSE(pressures[0][first] , 14832467.14, reltol); // opm
    CHECK_CLOSE(pressures[0][last ] , 15479883.47, reltol);
    CHECK_CLOSE(pressures[1][last ] , 15489883.47, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }

    const std::vector<double> s_ecl[3]{// eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22874042, 0.53397995, 0.78454906,  0.91542006, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.20039,     0.08458,    0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77125955, 0.46602005, 0.015063271, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    const std::vector<double> s_opm[3]{ // opm
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22892931226886132,  0.53406457830052489, 0.78457075254244724, 0.91539712466977541, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,                   0,                   0.20023624994125844,   0.084602875330224592, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77107068773113863, 0.46593542169947511, 0.015192997516294321, 0,      0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
            CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
        }
    }
}

void test_DeckWithLiveOil();
void test_DeckWithLiveOil()
{
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
    auto simulator = initSimulator<TypeTag>("equil_liveoil.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    // Initialize the fluid system
    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> comp(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 9.80665);
    const auto& pressures = comp.press();
    REQUIRE(pressures.size() == 3);
    REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    const double reltol_ecl = 1.0;
    CHECK_CLOSE(pressures[0][first], 1.48324e+07, reltol_ecl);  // eclipse
    CHECK_CLOSE(pressures[0][last],  1.54801e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][first], 1.49224e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][last],  1.54901e+07, reltol_ecl);

    CHECK_CLOSE(pressures[0][first], 1.483246714e7, reltol);  // opm
    CHECK_CLOSE(pressures[0][last],  1.547991652e7, reltol);
    CHECK_CLOSE(pressures[1][first], 1.492246714e7, reltol);
    CHECK_CLOSE(pressures[1][last],  1.548991652e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s_ecl[3]{ // eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22898, 0.53422, 0.78470, 0.91531, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,       0,       0.20073, 0.08469, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77102, 0.46578, 0.01458, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    const std::vector<double> s_opm[3]{ // opm
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22916963446461344, 0.53430490523774521, 0.78471886612242092, 0.91528324362210933, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,            0,            0.20057438297017782,   0.084716756377890667, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77083036553538653, 0.46569509476225479, 0.014706750907401245,  0,       0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
            CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
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
        CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        CHECK_CLOSE(rs[i], rs_ecl[i], reltol_ecl);
    }
}

void test_DeckWithLiveGas();
void test_DeckWithLiveGas()
{
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
    auto simulator = initSimulator<TypeTag>("equil_livegas.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> comp(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 9.80665);
    const auto& pressures = comp.press();
    REQUIRE(pressures.size() == 3);
    REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-3;
    const double reltol_ecl = 1.0;
    CHECK_CLOSE(pressures[0][first], 1.48215e+07, reltol_ecl);  // eclipse
    CHECK_CLOSE(pressures[0][last],  1.54801e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][first], 1.49115e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][last],  1.54901e+07, reltol_ecl);

    CHECK_CLOSE(pressures[0][first], 1.482150311e7, reltol);  // opm
    CHECK_CLOSE(pressures[0][last],  1.547988347e7, reltol);
    CHECK_CLOSE(pressures[1][first], 1.491150311e7, reltol);
    CHECK_CLOSE(pressures[1][last],  1.548988347e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s_ecl[3]{ // eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.24285614, 0.53869015, 0.78454906,  0.91542006, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.18311,     0.08458,    0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75714386, 0.46130988, 0.032345835, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    const std::vector<double> s_opm[3]{ // opm
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.24310545, 0.5388, 0.78458,    0.91540, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,      0.18288667, 0.0846,  0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.75689455, 0.4612, 0.03253333, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };
    for (int phase = 0; phase < 3; ++phase) {
        REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            CHECK_CLOSE(sats[phase][i], s_opm[phase][i], 100.*reltol);
            CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
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
        CHECK_CLOSE(rv[i], rv_opm[i], reltol);
        CHECK_CLOSE(rv[i], rv_ecl[i], reltol_ecl);
    }
}

void test_DeckWithRSVDAndRVVD();
void test_DeckWithRSVDAndRVVD()
{
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
    auto simulator = initSimulator<TypeTag>("equil_rsvd_and_rvvd.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> comp(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 9.80665);
    const auto& pressures = comp.press();
    REQUIRE(pressures.size() == 3);
    REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    const double reltol_ecl = 1.0;
    CHECK_CLOSE(pressures[0][first], 1.48350e+07, reltol_ecl);  // eclipse
    CHECK_CLOSE(pressures[0][last],  1.54794e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][first], 1.49250e+07, reltol_ecl);
    CHECK_CLOSE(pressures[1][last],  1.54894e+07, reltol_ecl);

    CHECK_CLOSE(pressures[0][first], 1.483499660e7, reltol);  // opm
    CHECK_CLOSE(pressures[0][last],  1.547924516e7, reltol);
    CHECK_CLOSE(pressures[1][first], 1.492499660e7, reltol);
    CHECK_CLOSE(pressures[1][last],  1.548924516e7, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    const std::vector<double> s_ecl[3]{ // eclipse
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22206347, 0.52871972, 0.78150368,  0.91819441,  1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0,          0.19656529,  0.081805572, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77793652, 0.47128031, 0.021931054, 0,           0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    const std::vector<double> s_opm[3]{ // opm
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.22231423543119974, 0.52882640735211706, 0.78152142505479982, 0.91816512259416283, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0, 0.19636279642563928, 0.08183487740583717, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.77768576456880023, 0.47117359264788294, 0.022115778519560897,    0,          0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    for (int phase = 0; phase < 3; ++phase) {
        REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
            CHECK_CLOSE(sats[phase][i], s_ecl[phase][i], reltol_ecl);
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
        CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        CHECK_CLOSE(rs[i], rs_ecl[i], reltol_ecl);
        CHECK_CLOSE(rv[i], rv_opm[i], reltol);
        CHECK_CLOSE(rv[i], rv_ecl[i], reltol_ecl);
    }
}


void test_DeckWithPBVDAndPDVD();
void test_DeckWithPBVDAndPDVD()
{
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
    auto simulator = initSimulator<TypeTag>("equil_pbvd_and_pdvd.DATA");
    const auto& eclipseState = simulator->vanguard().eclState();
    Opm::GridManager gm(eclipseState.getInputGrid());
    const UnstructuredGrid& grid = *(gm.c_grid());

    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> comp(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 9.80665);
    const auto& pressures = comp.press();
    REQUIRE(pressures.size() == 3);
    REQUIRE(int(pressures[0].size()) == grid.number_of_cells);

    const int first = 0, last = grid.number_of_cells - 1;
    // The relative tolerance is too loose to be very useful,
    // but the answer we are checking is the result of an ODE
    // solver, and it is unclear if we should check it against
    // the true answer or something else.
    const double reltol = 1.0e-6;
    CHECK_CLOSE(pressures[0][first], 14821552, reltol);
    CHECK_CLOSE(pressures[0][last],  15479828, reltol);
    CHECK_CLOSE(pressures[1][first], 14911552, reltol);
    CHECK_CLOSE(pressures[1][last],  15489828, reltol);

    const auto& sats = comp.saturation();
    // std::cout << "Saturations:\n";
    // for (const auto& sat : sats) {
    //     for (const double s : sat) {
    //         std::cout << s << ' ';
    //     }
    //     std::cout << std::endl;
    // }

    const std::vector<double> s_opm[3]{ // opm
        { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2426402656423233, 0.5383705390740118, 0.7844998821510003, 0.9152832369551807, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 0,   0,   0,   0,   0,   0,   0,   0,          0, 0.1817779931230221, 0.08471676304481934, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7573597343576767, 0.4616294609259882, 0.03372212472597758,    0,          0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    for (int phase = 0; phase < 3; ++phase) {
        REQUIRE(sats[phase].size() == s_opm[phase].size());
        for (size_t i = 0; i < s_opm[phase].size(); ++i) {
            //std::cout << std::setprecision(10) << sats[phase][i] << '\n';
            CHECK_CLOSE(sats[phase][i], s_opm[phase][i], reltol);
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
        CHECK_CLOSE(rs[i], rs_opm[i], reltol);
        CHECK_CLOSE(rv[i], rv_opm[i], reltol);
    }
}

void test_DeckWithSwatinit();
void test_DeckWithSwatinit()
{
#if 0
    typedef typename TTAG(TestEquilTypeTag) TypeTag;
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
        double sw = s[0][c];
        double so = s[1][c];
        double sg = s[2][c];
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
    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> compScaled(materialLawManagerScaled, eclipseState, simulator->vanguard().grid(), 9.81, true);
    // don't apply swatinit
    Ewoms::EQUIL::DeckDependent::InitialStateComputer<TypeTag> compUnscaled(*simulator->problem().materialLawManager(), eclipseState, simulator->vanguard().grid(), 9.81, false);

    // compute pc
    std::vector<double> pc_scaled(numCells * FluidSystem::numPhases);
    for (int c = 0; c < numCells; ++c) {
        std::vector<double> pc = {0,0,0};
        double sw = compScaled.saturation().data()[0][c];
        double so = compScaled.saturation().data()[1][c];
        double sg = compScaled.saturation().data()[2][c];

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
        double sw = compUnscaled.saturation().data()[0][c];
        double so = compUnscaled.saturation().data()[1][c];
        double sg = compUnscaled.saturation().data()[2][c];

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
    const double reltol = 1.0e-3;
    for (int phase = 0; phase < 3; ++phase) {
        for (size_t i = 0; i < 20; ++i) {
            CHECK_CLOSE( pc_original[3*i + phase ], pc_unscaled[3*i + phase ], reltol);
            CHECK_CLOSE( pc_scaled_truth[3*i + phase], pc_scaled[3*i + phase ], reltol);
        }
    }

    for (int phase = 0; phase < 3; ++phase) {
        for (size_t i = 0; i < 20; ++i) {
            CHECK_CLOSE(compUnscaled.saturation()[phase][i], s[phase][i], reltol);
            CHECK_CLOSE(compScaled.saturation()[phase][i], swatinit[phase][i], reltol);
        }
    }
#endif
}

int main(int argc, char** argv)
{
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    typedef TTAG(TestEquilTypeTag) TypeTag;
    Ewoms::registerAllParameters_<TypeTag>();

    test_PhasePressure();
    test_CellSubset();
    test_RegMapping();
    test_DeckAllDead();
    test_CapillaryInversion();
    test_DeckWithCapillary();
    test_DeckWithCapillaryOverlap();
    test_DeckWithLiveOil();
    test_DeckWithLiveGas();
    test_DeckWithRSVDAndRVVD();
    test_DeckWithPBVDAndPDVD();
    //test_DeckWithSwatinit();

    return 0;
}
