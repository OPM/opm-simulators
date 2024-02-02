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

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>

#include <ebos/ecltransmissibility.hh>

using namespace Opm;

const int dimWorld = 3;
const int cartDims[3] = {8,15,3};

// Class extending EclTransmissibility, such that we can access the protected member trans_ to check its contents
template<class Grid, class GridView, class ElementMapper, class CartesianIndexMapper, class Scalar>
class Transmissibility : public EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>
{
    using ParentType = EclTransmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>;
    public:
        Transmissibility(const EclipseState& eclState,
                        const GridView& gridView,
                        const CartesianIndexMapper& cartMapper,
                        const Grid& grid,
                        std::function<std::array<double,dimWorld>(int)> centroids,
                        bool enableEnergy,
                        bool enableDiffusivity,
                        bool enableDispersivity) : ParentType(eclState,gridView,cartMapper,grid,centroids,enableEnergy,enableDiffusivity,enableDispersivity) {}
        auto getTransmissibilitymap() {
            return this->trans_;
        }
};

int main(int argc, char** argv )
{
    Dune::MPIHelper::instance(argc, argv);

    using Grid = Dune::CpGrid;
    using GridView = Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using Transmissibility = Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,double>;

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
    EclipseGrid eclGrid(deck);
    EclipseState eclState(deck);

    grid.processEclipseFormat(&eclGrid, &eclState, false, false, false);
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

    // Check if the transmissibilities of the NNCs that were added manually are either not contained in the transmissibility array (because they might be on a different process) or 0.0
    if (transmissibilityMap.count(details::isId(0,129)) > 0)
        assert(eclTransmissibility.transmissibility(0,129) == 0.0);
    if (transmissibilityMap.count(details::isId(0,258)) > 0)
        assert(eclTransmissibility.transmissibility(0,258) == 0.0);

    // If there is a non-zero transmissibility in the map, ensure that it is form a neighboring connection
    for (auto&& trans : transmissibilityMap) {
        if (trans.second != 0.0) {
            const auto& id = trans.first;
            const auto& elements = details::isIdReverse(id);
            int gc1 = std::min(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));
            int gc2 = std::max(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));
            assert(gc2 - gc1 == 1 || gc2 - gc1 == cartDims[0] || gc2 - gc1 == cartDims[0]*cartDims[1] || gc2 - gc1 == 0);
        }
    }
    return 0;
}
