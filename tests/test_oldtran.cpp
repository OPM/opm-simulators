/*
  Copyright 2026 Equinor ASA

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

#define BOOST_TEST_MODULE OLDTRANTest
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/flow/Transmissibility.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <string>

using namespace Opm;

const int dimWorld = 3;

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
    std::string errorstring;
    int errorcode;
};

void MPI_err_handler(MPI_Comm*, int* err_code, ...)
{
    std::vector<char> err_string(MPI_MAX_ERROR_STRING);
    int err_length;
    MPI_Error_string(*err_code, err_string.data(), &err_length);
    std::string s(err_string.data(), err_length);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *err_code);
}
#endif

bool
init_unit_test_func()
{
    return true;
}

// Class extending EclTransmissibility, such that we can access the protected member trans_ to check its contents
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
class TestTransmissibility : public Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>
{
    using ParentType = Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>;
    public:
        TestTransmissibility(const EclipseState& eclState,
                             const GridView& gridView,
                             const CartesianIndexMapper& cartMapper,
                             const Grid& grid,
                             std::function<std::array<double,dimWorld>(int)> centroids,
                             bool enableEnergy,
                             bool enableDiffusivity,
                             bool enableDispersivity)
            : ParentType(eclState,gridView,cartMapper,grid,centroids,
                         enableEnergy,enableDiffusivity,enableDispersivity) {}
        auto getTransmissibilitymap() {
            return this->trans_;
        }
};

namespace {

using Grid = Dune::CpGrid;
using GridView = Grid::LeafGridView;
using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
using Trans = TestTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,double>;

// 3x1x2 grid of 100x100x10 m cells with the tops of the first layer given. The cells of the
// columns with the tops 1000, 1030 and 1015 are
//   column 1: 1000-1010 and 1010-1020
//   column 2: 1030-1040 and 1040-1050
//   column 3: 1015-1025 and 1025-1035
// so the columns 1 and 2 do not overlap, and the first layer of column 2 only overlaps the
// second layer of column 3.
std::string deckString(const std::string& tops, const std::string& tranKeyword,
                       const std::string& runspecKeywords = "")
{
    return "RUNSPEC\n"
           + runspecKeywords +
           "DIMENS\n"
           " 3 1 2 /\n"
           "GRID\n"
           + tranKeyword +
           "DX\n"
           " 6*100 /\n"
           "DY\n"
           " 6*100 /\n"
           "DZ\n"
           " 6*10 /\n"
           "TOPS\n"
           " " + tops + " /\n"
           "PORO\n"
           " 6*0.25 /\n"
           "PERMX\n"
           " 6*1000 /\n"
           "PERMY\n"
           " 6*1000 /\n"
           "PERMZ\n"
           " 6*100 /\n"
           "END\n";
}

struct Setup
{
    explicit Setup(const std::string& deckData)
        : deck(Parser{}.parseString(deckData))
        , eclState(deck)
    {
        grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);
        cartMapper = std::make_unique<CartesianIndexMapper>(grid);
        const auto& gridView = grid.leafGridView();
        ElementMapper elementMapper(gridView, Dune::mcmgElementLayout());
        std::vector<std::array<double,Dune::CpGrid::dimensionworld>> cellCenters(gridView.size(0));
        for (const auto& element : elements(gridView)) {
            const auto center = element.geometry().center();
            cellCenters[elementMapper.index(element)] = {center[0], center[1], center[2]};
        }
        auto centroids = [cellCenters](int cell) { return cellCenters[cell]; };
        transmissibility = std::make_unique<Trans>(eclState, grid.leafGridView(), *cartMapper, grid,
                                                   centroids, false, false, false);
        transmissibility->update(true);
        map = transmissibility->getTransmissibilitymap();
    }

    // Transmissibility between the cells, or zero if the grid does not connect them
    double trans(const int cell1, const int cell2) const
    {
        const auto it = map.find(details::isId(cell1, cell2));
        return it == map.end() ? 0.0 : it->second;
    }

    Deck deck;
    EclipseState eclState;
    Grid grid;
    std::unique_ptr<CartesianIndexMapper> cartMapper;
    std::unique_ptr<Trans> transmissibility;
    decltype(std::declval<Trans>().getTransmissibilitymap()) map;
};

// cell index = i + 3*k, as all cells are active and the grid has one row
constexpr int c11 = 0, c21 = 1, c31 = 2, c12 = 3, c22 = 4, c32 = 5;

const double perm = 1000 * Opm::prefix::milli * Opm::unit::darcy;
const double permZ = 100 * Opm::prefix::milli * Opm::unit::darcy;
// The full face between lateral neighbours is 100 m x 10 m and the cells are 100 m long
const double lateralTrans = perm * (100.0 * 10.0) / 100.0;
const double verticalTrans = permZ * (100.0 * 100.0) / 10.0;

} // anonymous namespace

BOOST_AUTO_TEST_CASE(OldTranStepped)
{
    Setup setup(deckString("1000 1030 1015", ""));
    BOOST_REQUIRE(setup.eclState.getInputGrid().useOldTran());

    // neighbours in a layer have the full face, also if the cells do not overlap
    BOOST_CHECK_CLOSE(setup.trans(c11, c21), lateralTrans, 1e-6);
    BOOST_CHECK_CLOSE(setup.trans(c21, c31), lateralTrans, 1e-6);
    BOOST_CHECK_CLOSE(setup.trans(c12, c22), lateralTrans, 1e-6);
    BOOST_CHECK_CLOSE(setup.trans(c22, c32), lateralTrans, 1e-6);

    // the cell in the first layer of column 2 overlaps the one in the second layer of
    // column 3, but they are not connected
    BOOST_CHECK_EQUAL(setup.trans(c21, c32), 0.0);

    // vertical connections are not changed
    BOOST_CHECK_CLOSE(setup.trans(c11, c12), verticalTrans, 1e-6);
    BOOST_CHECK_CLOSE(setup.trans(c21, c22), verticalTrans, 1e-6);
    BOOST_CHECK_CLOSE(setup.trans(c31, c32), verticalTrans, 1e-6);
}

BOOST_AUTO_TEST_CASE(NewTranStepped)
{
    Setup setup(deckString("1000 1030 1015", "NEWTRAN\n"));
    BOOST_REQUIRE(!setup.eclState.getInputGrid().useOldTran());

    // the connections follow the overlap of the cells, as for corner-point grids
    BOOST_CHECK_EQUAL(setup.trans(c11, c21), 0.0);
    BOOST_CHECK_EQUAL(setup.trans(c21, c31), 0.0);
    BOOST_CHECK_GT(setup.trans(c21, c32), 0.0);
    BOOST_CHECK_CLOSE(setup.trans(c11, c12), verticalTrans, 1e-6);
}

// OLDTRAN is not implemented for the thermal, diffusive and dispersive transmissibilities,
// NEWTRAN is used for all of them in these runs
BOOST_AUTO_TEST_CASE(NewTranWithThermalDeck)
{
    Setup setup(deckString("1000 1030 1015", "", "THERMAL\n"));

    BOOST_CHECK_EQUAL(setup.trans(c11, c21), 0.0);
    BOOST_CHECK_EQUAL(setup.trans(c21, c31), 0.0);
    BOOST_CHECK_GT(setup.trans(c21, c32), 0.0);
}

BOOST_AUTO_TEST_CASE(OldTranPartialOverlap)
{
    // The cells of column 2 are 5 m below those of column 1, so they overlap over 5 m
    Setup stepped(deckString("1000 1005 1000", ""));
    BOOST_REQUIRE(stepped.eclState.getInputGrid().useOldTran());

    // OLDTRAN uses the full face, not the part where the cells overlap
    BOOST_CHECK_CLOSE(stepped.trans(c11, c21), lateralTrans, 1e-6);
    BOOST_CHECK_CLOSE(stepped.trans(c12, c22), lateralTrans, 1e-6);
    // and the overlap of the first layer of column 2 with the second layer of column 1
    // is not a connection
    BOOST_CHECK_EQUAL(stepped.trans(c12, c21), 0.0);

    Setup newtran(deckString("1000 1005 1000", "NEWTRAN\n"));
    BOOST_CHECK_LT(newtran.trans(c11, c21), lateralTrans);
    BOOST_CHECK_GT(newtran.trans(c11, c21), 0.0);
    BOOST_CHECK_GT(newtran.trans(c12, c21), 0.0);
}

BOOST_AUTO_TEST_CASE(OldTranFlat)
{
    // Flat cells give the same transmissibilities with OLDTRAN and NEWTRAN
    Setup newtran(deckString("1000 1000 1000", "NEWTRAN\n"));
    Setup oldtran(deckString("1000 1000 1000", ""));
    BOOST_CHECK(!newtran.eclState.getInputGrid().useOldTran());
    BOOST_CHECK(oldtran.eclState.getInputGrid().useOldTran());

    for (int cell1 = 0; cell1 < 6; ++cell1) {
        for (int cell2 = cell1 + 1; cell2 < 6; ++cell2) {
            // the +1 avoids a relative comparison of zeros
            BOOST_CHECK_CLOSE(oldtran.trans(cell1, cell2) + 1.0, newtran.trans(cell1, cell2) + 1.0, 1e-9);
        }
    }
    BOOST_CHECK_CLOSE(oldtran.trans(c11, c21), lateralTrans, 1e-6);
    BOOST_CHECK_CLOSE(newtran.trans(c11, c21), lateralTrans, 1e-6);
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
