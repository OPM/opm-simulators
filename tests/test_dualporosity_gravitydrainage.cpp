/*
  Copyright 2026 TNO

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

#define BOOST_TEST_MODULE TestDualPorosityGravityDrainage
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/models/blackoil/blackoildualporositygravitydrainage.hh>

#include <opm/simulators/flow/Transmissibility.hpp>

#include <array>
#include <string>

using namespace Opm;

bool
init_unit_test_func()
{
    return true;
}

namespace {

// End points used by the fraction tests: the mobile water fraction moves
// between swco = 0.2 and 1 - scrhy = 0.7 while the level rises, and between
// swcr = 0.25 and 1 - scohy = 0.85 while it falls.
constexpr DualPorosityFractions::WaterFractionEndPoints waterEndPoints{0.2, 0.25, 0.15, 0.3};

} // anonymous namespace

BOOST_AUTO_TEST_CASE(InitialWaterFractionLocatesTheContact)
{
    using DualPorosityFractions::initialWaterFraction;

    // (swi - swco) / (1 - scohy - swco), clamped to [0, 1].
    BOOST_CHECK_CLOSE(initialWaterFraction(0.525, waterEndPoints), 0.5, 1e-12);
    BOOST_CHECK_EQUAL(initialWaterFraction(0.2, waterEndPoints), 0.0);
    BOOST_CHECK_EQUAL(initialWaterFraction(0.1, waterEndPoints), 0.0);
    BOOST_CHECK_EQUAL(initialWaterFraction(0.9, waterEndPoints), 1.0);
}

BOOST_AUTO_TEST_CASE(WaterFractionIsContinuousAtTheInitialSaturation)
{
    using DualPorosityFractions::waterFraction;
    using Eval = DenseAd::Evaluation<double, 1>;

    const double swi = 0.525;
    const double xwi = 0.5;

    // Both branches give the initial fraction at the initial saturation.
    BOOST_CHECK_CLOSE(waterFraction(swi, swi, xwi, waterEndPoints), xwi, 1e-12);
    BOOST_CHECK_CLOSE(waterFraction(swi - 1e-12, swi, xwi, waterEndPoints), xwi, 1e-9);

    // Rising level: (sw - xwi (scrhy - scohy) - swco) / (1 - scrhy - swco),
    // slope 1 / 0.5 = 2.
    const Eval rising = waterFraction(Eval::createVariable(0.6, 0), swi, xwi, waterEndPoints);
    BOOST_CHECK_CLOSE(rising.value(), (0.6 - 0.5 * 0.15 - 0.2) / 0.5, 1e-12);
    BOOST_CHECK_CLOSE(rising.derivative(0), 2.0, 1e-12);

    // Falling level: (sw - xwi (swcr - swco) - swco) / (1 - scohy - swcr),
    // slope 1 / 0.6.
    const Eval falling = waterFraction(Eval::createVariable(0.4, 0), swi, xwi, waterEndPoints);
    BOOST_CHECK_CLOSE(falling.value(), (0.4 - 0.5 * 0.05 - 0.2) / 0.6, 1e-12);
    BOOST_CHECK_CLOSE(falling.derivative(0), 1.0 / 0.6, 1e-12);

    // Clamped to [0, 1], with no slope where clamped.
    const Eval full = waterFraction(Eval::createVariable(1.0, 0), swi, xwi, waterEndPoints);
    BOOST_CHECK_EQUAL(full.value(), 1.0);
    BOOST_CHECK_EQUAL(full.derivative(0), 0.0);
    BOOST_CHECK_EQUAL(waterFraction(0.0, swi, xwi, waterEndPoints), 0.0);
}

BOOST_AUTO_TEST_CASE(GasFractionMirrorsTheWaterFraction)
{
    using namespace DualPorosityFractions;

    // The same numbers in the gas roles give the same fractions.
    const GasFractionEndPoints gasEndPoints{0.2, 0.25, 0.15, 0.3};
    BOOST_CHECK_CLOSE(initialGasFraction(0.525, gasEndPoints), 0.5, 1e-12);
    BOOST_CHECK_CLOSE(gasFraction(0.6, 0.525, 0.5, gasEndPoints),
                      waterFraction(0.6, 0.525, 0.5, waterEndPoints), 1e-12);
    BOOST_CHECK_CLOSE(gasFraction(0.4, 0.525, 0.5, gasEndPoints),
                      waterFraction(0.4, 0.525, 0.5, waterEndPoints), 1e-12);
}

BOOST_AUTO_TEST_CASE(GravityHeadsAndTheirPhaseShares)
{
    using namespace DualPorosityFractions;

    // g dz (rho_o - rho_g) (xg_a - xg_b) = 9.81 * 10 * 600 * 0.4
    const double gasOil = gasOilGravityHead(9.81, 10.0, 800.0, 200.0, 0.6, 0.2);
    BOOST_CHECK_CLOSE(gasOil, 23544.0, 1e-10);

    // g dz (rho_w - rho_o) (xw_a - xw_b) = 9.81 * 10 * 200 * 0.2
    const double waterOil = waterOilGravityHead(9.81, 10.0, 1000.0, 800.0, 0.3, 0.1);
    BOOST_CHECK_CLOSE(waterOil, 3924.0, 1e-10);

    // A head changes sign with the order of the two cells, and vanishes for a
    // zero block height.
    BOOST_CHECK_CLOSE(gasOilGravityHead(9.81, 10.0, 800.0, 200.0, 0.2, 0.6), -gasOil, 1e-10);
    BOOST_CHECK_EQUAL(waterOilGravityHead(9.81, 0.0, 1000.0, 800.0, 0.3, 0.1), 0.0);

    // Each pair head is shared half and half; oil takes part in both pairs.
    const auto shares = phaseGravityHeads(gasOil, waterOil);
    BOOST_CHECK_CLOSE(shares.gas, 11772.0, 1e-10);
    BOOST_CHECK_CLOSE(shares.water, 1962.0, 1e-10);
    BOOST_CHECK_CLOSE(shares.oil, -13734.0, 1e-10);
}

namespace {

Deck gravityDrainageDeck(const std::string& runspecKeywords)
{
    // 2x1x4: matrix cells cart 0-3 (k = 1, 2), fracture cells cart 4-7 (k = 3, 4),
    // co-located twins.
    return Parser{}.parseString(std::string(R"(RUNSPEC
OIL
WATER
DIMENS
 2 1 4 /
)") + runspecKeywords + R"(GRID
DX
 8*100. /
DY
 8*100. /
DZ
 8*10. /
TOPS
 2*2000. 2*2010. 2*2000. 2*2010. /
PORO
 4*0.20 4*0.01 /
PERMX
 4*1.0 4*1000.0 /
PERMY
 4*1.0 4*1000.0 /
PERMZ
 4*0.5 4*500.0 /
SIGMA
 0.12 /
DZMTRX
 5.0 /
END)");
}

template <class Check>
void withTransmissibility(const Deck& deck, Check&& check)
{
    using Grid = Dune::CpGrid;
    using GridView = Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using Trans = Transmissibility<Grid, GridView, ElementMapper, CartesianIndexMapper, double>;

    Grid grid;
    EclipseState eclState(deck);
    grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);
    const auto& gridView = grid.leafGridView();
    CartesianIndexMapper cartMapper(grid);
    auto centroids = [](int) { return std::array<double, Dune::CpGrid::dimensionworld>{}; };

    Trans trans(eclState, gridView, cartMapper, grid, centroids, false, false, false);
    trans.update(true);

    // All cells are active, so element indices coincide with Cartesian ones.
    check(trans);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(TwinPairsAreRegisteredForGravityDrainage)
{
    withTransmissibility(gravityDrainageDeck("DUALPORO\nGRAVDR\n"), [](const auto& trans)
    {
        for (unsigned matrix = 0; matrix < 4; ++matrix) {
            const unsigned fracture = matrix + 4;
            const auto fromMatrix = trans.dualPorosityGravityDrainagePair(matrix, fracture);
            BOOST_CHECK(fromMatrix.isTwinPair);
            BOOST_CHECK(fromMatrix.firstIsMatrix);

            const auto fromFracture = trans.dualPorosityGravityDrainagePair(fracture, matrix);
            BOOST_CHECK(fromFracture.isTwinPair);
            BOOST_CHECK(!fromFracture.firstIsMatrix);
        }

        // Neighbours within a continuum are not twin pairs.
        BOOST_CHECK(!trans.dualPorosityGravityDrainagePair(0, 1).isTwinPair);
        BOOST_CHECK(!trans.dualPorosityGravityDrainagePair(4, 5).isTwinPair);
        BOOST_CHECK(!trans.dualPorosityGravityDrainagePair(0, 2).isTwinPair);
    });
}

BOOST_AUTO_TEST_CASE(NoTwinPairsWithoutGravityDrainage)
{
    withTransmissibility(gravityDrainageDeck("DUALPORO\n"), [](const auto& trans)
    {
        for (unsigned matrix = 0; matrix < 4; ++matrix) {
            BOOST_CHECK(!trans.dualPorosityGravityDrainagePair(matrix, matrix + 4).isTwinPair);
        }
    });
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
