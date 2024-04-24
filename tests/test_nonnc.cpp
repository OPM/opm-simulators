/*
  Copyright 2024 Equinor ASA

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

#define BOOST_TEST_MODULE NONNCTest
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/flow/Transmissibility.hpp>

using namespace Opm;

const int dimWorld = 3;
const int cartDims[3] = {8,15,3};

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

BOOST_AUTO_TEST_CASE(NoNNC)
{
    using Grid = Dune::CpGrid;
    using GridView = Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using Transmissibility = TestTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,double>;

    Parser parser;

    // Deck with 2 NNCs
    const auto deck = parser.parseString( R"(RUNSPEC
                                            NONNC
                                            DIMENS
                                             8 15 3 /
                                            GRID
                                            DX
                                              360*10./
                                            DY
                                               360*10./
                                            DZ
                                               360*1./
                                            TOPS
                                               360*100./

                                            PORO
                                               0. 0.25 0. 357*0.25/
                                            PERMX
                                                360*1000./
                                            PERMY
                                                360*1000./
                                            PERMZ
                                                360*10./

                                            NNC
                                            -- I1 J1 K1  I2 J2 K2 Trans
                                                1  1  1   2  2  2  1000.0 / --- connection between 0 and 129
                                                1  1  1   3  3  3  1000.0 / --- connection between 0 and 258
                                            /

                                            END)");
    Grid grid;
    EclipseState eclState(deck);

    grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);
    const auto& gridView = grid.leafGridView();

    CartesianIndexMapper cartMapper =  Dune::CartesianIndexMapper<Grid>(grid);

    auto centroids = [](int) { return std::array<double,Dune::CpGrid::dimensionworld>{}; };
    Transmissibility eclTransmissibility(eclState,
                                         gridView,
                                         cartMapper,
                                         grid,
                                         centroids,
                                         false,false,false);
    // Call update, true indicates that update is called on all processes
    eclTransmissibility.update(true);

    auto transmissibilityMap = eclTransmissibility.getTransmissibilitymap();

    // Check that the transmissibilities of the NNCs that were added manually are
    // not contained in the transmissibility array or 0.0
    BOOST_CHECK(transmissibilityMap.count(details::isId(0,129)) == 0 ||
                eclTransmissibility.transmissibility(0,129) == 0.0);
    BOOST_CHECK(transmissibilityMap.count(details::isId(0,258)) == 0 ||
                eclTransmissibility.transmissibility(0,258) == 0.0);

    // If there is a non-zero transmissibility in the map, ensure that it is form a neighboring connection
    for (auto&& trans : transmissibilityMap) {
        if (trans.second != 0.0) {
            const auto& id = trans.first;
            const auto& elements = details::isIdReverse(id);
            int gc1 = std::min(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));
            int gc2 = std::max(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));
            BOOST_CHECK(gc2 - gc1 == 1 ||
                        gc2 - gc1 == cartDims[0] ||
                        gc2 - gc1 == cartDims[0]*cartDims[1] ||
                        gc2 - gc1 == 0);
        }
    }
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
