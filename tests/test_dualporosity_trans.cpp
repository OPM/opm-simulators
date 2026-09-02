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

#define BOOST_TEST_MODULE DualPorosityTransmissibility
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/material/densead/Evaluation.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/flow/Transmissibility.hpp>

#include <algorithm>
#include <array>
#include <map>
#include <string>

using namespace Opm;

bool
init_unit_test_func()
{
    return true;
}

namespace {

template <class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
class TestTransmissibility : public Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>
{
    using ParentType = Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>;
public:
    TestTransmissibility(const EclipseState& eclState,
                         const GridView& gridView,
                         const CartesianIndexMapper& cartMapper,
                         const Grid& grid,
                         std::function<std::array<double,Dune::CpGrid::dimensionworld>(int)> centroids,
                         bool enableEnergy,
                         bool enableDiffusivity,
                         bool enableDispersivity)
        : ParentType(eclState,gridView,cartMapper,grid,centroids,
                     enableEnergy,enableDiffusivity,enableDispersivity) {}

    auto getTransmissibilitymap() { return this->trans_; }
};

Deck dualContinuumDeck(bool nodppm, bool dualperm)
{
    // 2x1x4: matrix cells cart 0-3 (k=0,1), fracture cells cart 4-7 (k=2,3),
    // co-located twins, sigma one field value.
    const std::string deckData = std::string(R"(RUNSPEC
OIL
WATER
DIMENS
 2 1 4 /
)") + (dualperm ? "DUALPERM\n" : "DUALPORO\n")
    + (nodppm ? "NODPPM\n" : "") + R"(GRID
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
)" "END";
    return Parser{}.parseString(deckData);
}

Deck singlePorosityDeck(double permMd, double permzMd, double poro, double topDepth)
{
    // The 2x1x2 single-porosity analogue of one half of the dual-continuum
    // grid: identical cell boxes, one permeability, one porosity.  The
    // analogue must sit at the SAME depth as the half it mirrors: the test
    // harness supplies placeholder centroids, which makes the computed
    // transmissibilities depth-sensitive, so only same-depth grids compare
    // exactly.
    const std::string deckData = std::string(R"(RUNSPEC
OIL
WATER
DIMENS
 2 1 2 /
GRID
DX
 4*100. /
DY
 4*100. /
DZ
 4*10. /
TOPS
 2*)") + std::to_string(topDepth) + " 2*" + std::to_string(topDepth + 10.0) + R"( /
PORO
 4*)" + std::to_string(poro) + R"( /
PERMX
 4*)" + std::to_string(permMd) + R"( /
PERMY
 4*)" + std::to_string(permMd) + R"( /
PERMZ
 4*)" + std::to_string(permzMd) + R"( /
END)";
    return Parser{}.parseString(deckData);
}


struct TransResult {
    std::map<std::size_t, double> byId;
    double trans(std::size_t c1, std::size_t c2) const {
        auto it = byId.find(details::isId(c1, c2));
        return it == byId.end() ? 0.0 : it->second;
    }
};

TransResult computeTransFromDeck(const Deck& deck)
{
    using Grid = Dune::CpGrid;
    using GridView = Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using Transmissibility = TestTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,double>;

    Grid grid;
    EclipseState eclState(deck);
    grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);
    const auto& gridView = grid.leafGridView();
    CartesianIndexMapper cartMapper = Dune::CartesianIndexMapper<Grid>(grid);
    auto centroids = [](int) { return std::array<double,Dune::CpGrid::dimensionworld>{}; };

    Transmissibility eclTransmissibility(eclState, gridView, cartMapper, grid,
                                         centroids, false, false, false);
    eclTransmissibility.update(true);

    // All cells are active, so compressed indices coincide with cartesian ones.
    TransResult result;
    for (const auto& t : eclTransmissibility.getTransmissibilitymap()) {
        result.byId.emplace(t.first, t.second);
    }
    return result;
}

TransResult computeTrans(bool nodppm, bool dualperm = false)
{
    return computeTransFromDeck(dualContinuumDeck(nodppm, dualperm));
}

// The coupling transmissibility of every twin pair:
// matrix PERMX (1 mD, SI) * bulk volume (1e5 m3) * sigma (0.12 1/m2).
constexpr double expectedCoupling = 9.869232667160130e-16 * 1.0e5 * 0.12;

// No non-zero entry may cross the halves except an exact twin pair.
void checkOnlySameHalfOrTwin(const TransResult& res)
{
    for (const auto& t : res.byId) {
        if (t.second != 0.0) {
            const auto elements = details::isIdReverse(t.first);
            const auto g1 = static_cast<std::size_t>(std::min(elements.first, elements.second));
            const auto g2 = static_cast<std::size_t>(std::max(elements.first, elements.second));
            const bool sameHalf = (g1 < 4) == (g2 < 4);
            const bool twinPair = (g2 - g1 == 4);
            BOOST_CHECK(sameHalf || twinPair);
        }
    }
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(DualPorosityTransPolicy)
{
    const auto res = computeTrans(/*nodppm=*/true);

    // Matrix-matrix connections carry NO flow (single permeability).
    for (const auto& [c1, c2] : {std::pair<std::size_t,std::size_t>{0,1}, {2,3}, {0,2}, {1,3}}) {
        BOOST_CHECK_EQUAL(res.trans(c1, c2), 0.0);
    }

    // Fracture-fracture connections flow through their real geometry.
    BOOST_CHECK_GT(res.trans(4, 5), 0.0);
    BOOST_CHECK_GT(res.trans(6, 7), 0.0);
    BOOST_CHECK_GT(res.trans(4, 6), 0.0);
    BOOST_CHECK_GT(res.trans(5, 7), 0.0);

    // Every twin pair is coupled with exactly the sigma-NNC transmissibility.
    for (std::size_t g = 0; g < 4; ++g) {
        BOOST_CHECK_CLOSE(res.trans(g, g + 4), expectedCoupling, 1e-4);
    }

    checkOnlySameHalfOrTwin(res);
}

BOOST_AUTO_TEST_CASE(DualPorosityPermScaling)
{
    // Without NODPPM the fracture permeability is scaled by the fracture
    // porosity (0.01): fracture-fracture transmissibilities scale with it,
    // the coupling (computed from MATRIX permeability) does not.
    const auto plain  = computeTrans(/*nodppm=*/true);
    const auto scaled = computeTrans(/*nodppm=*/false);

    BOOST_CHECK_CLOSE(scaled.trans(4, 5), 0.01 * plain.trans(4, 5), 1e-6);
    BOOST_CHECK_CLOSE(scaled.trans(4, 6), 0.01 * plain.trans(4, 6), 1e-6);

    for (std::size_t g = 0; g < 4; ++g) {
        BOOST_CHECK_CLOSE(scaled.trans(g, g + 4), plain.trans(g, g + 4), 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(DualPermeabilityTransPolicy)
{
    const auto res = computeTrans(/*nodppm=*/true, /*dualperm=*/true);

    // Matrix-to-matrix transmissibilities are calculated from the matrix
    // cells' permeability with the standard formulas: every matrix-matrix
    // value must equal the corresponding connection of a single-porosity
    // grid with the same cell boxes, permeability and porosity.
    // Matrix half occupies depths 2000-2020 in the doubled grid.
    const auto spMatrix = computeTransFromDeck(singlePorosityDeck(1.0, 0.5, 0.20, 2000.0));
    for (const auto& [c1, c2] : {std::pair<std::size_t,std::size_t>{0,1}, {2,3}, {0,2}, {1,3}}) {
        BOOST_CHECK_GT(res.trans(c1, c2), 0.0);
        BOOST_CHECK_CLOSE(res.trans(c1, c2), spMatrix.trans(c1, c2), 1e-6);
    }

    // The fracture half likewise matches its own single-porosity analogue.
    // Fracture half occupies depths 2020-2040 (natural stacking).
    const auto spFracture = computeTransFromDeck(singlePorosityDeck(1000.0, 500.0, 0.01, 2020.0));
    for (const auto& [c1, c2] : {std::pair<std::size_t,std::size_t>{0,1}, {2,3}, {0,2}, {1,3}}) {
        BOOST_CHECK_CLOSE(res.trans(c1 + 4, c2 + 4), spFracture.trans(c1, c2), 1e-6);
    }

    // Twin coupling is unchanged: still exactly the sigma-NNC value.
    for (std::size_t g = 0; g < 4; ++g) {
        BOOST_CHECK_CLOSE(res.trans(g, g + 4), expectedCoupling, 1e-4);
    }

    checkOnlySameHalfOrTwin(res);
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
