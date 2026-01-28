// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
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
#include <config.h>

#define BOOST_TEST_MODULE TestFacePropertiesTPSA
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/flow/FacePropertiesTPSA.hpp>

#include <string>
#include <vector>

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

bool init_unit_test_func()
{
    return true;
}

static Opm::Deck createDeck()
{
    // Deck to test
    std::string deck_string = R"(
        RUNSPEC
        DIMENS
            2 2 2 /

        GRID
        DX
            8*50.0 /
        DY
            8*50.0 /
        DZ
            4*30.0 4*50.0 /
        TOPS
            3*1000.0 1*1020 /

        PORO
            8*0.3 /
        SMODULUS
            8*1.0 /

        END
    )";
    Opm::Parser parser;
    return parser.parseString(deck_string);
}

BOOST_AUTO_TEST_CASE(SimpleGridWithNNC)
{
    // using declarations
    using Grid = Dune::CpGrid;
    using GridView = Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using FacePropertiesTPSA = Opm::FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, double>;
    using DimVector = Dune::FieldVector<double, GridView::dimensionworld>;

    // Test simple 2x2x2 grid with one column shifted downwards
    Opm::Deck deck = createDeck();
    Opm::EclipseState eclState(deck);

    Grid grid;
    grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);
    const auto& gridView = grid.leafGridView();

    CartesianIndexMapper cartMapper =  Dune::CartesianIndexMapper<Grid>(grid);
    auto centroids = [&eclState, &cartMapper](int index)
        { return eclState.getInputGrid().getCellCenter(cartMapper.cartesianIndex(index)); };

    // Init. FacePropertiesTPSA and calculate all properties
    FacePropertiesTPSA faceProps(eclState,
                             gridView,
                             cartMapper,
                             grid,
                             centroids);
    faceProps.update();

    // Check face properties
    // Natural neighbor checks between cells 0 and 4
    const unsigned elem1 = 0;
    const unsigned elem2 = 4;
    const double normDist = 40.0;
    const double weightAvg_04 = 0.375;
    const double weightAvg_40 = 0.625;  // = 1 - weightAvg_04
    const double weightProd = 3.75e-16;
    const DimVector normal_04 = {0.0, 0.0, 1.0};
    const DimVector normal_40 = {0.0, 0.0, -1.0};  // = -1 * normal_04
    BOOST_CHECK_CLOSE(faceProps.weightAverage(elem1, elem2), weightAvg_04, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.weightAverage(elem2, elem1), weightAvg_40, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.weightProduct(elem1, elem2), weightProd, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.normalDistance(elem1, elem2), normDist, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.normalDistance(elem2, elem1), normDist, 1.0e-8);
    BOOST_CHECK_EQUAL(faceProps.cellFaceNormal(elem1, elem2), normal_04);
    BOOST_CHECK_EQUAL(faceProps.cellFaceNormal(elem2, elem1), normal_40);

    // NNC checks between cell 3 and 5
    const unsigned elemNNC1 = 3;
    const unsigned elemNNC2 = 5;
    const double normDistNNC = 50.0;
    const double weightAvgNNC = 0.5;
    const double weightProdNNC = 6.25e-16;
    const DimVector normalNNC_35 = {0.0, -1.0, 0.0};
    const DimVector normalNNC_53 = {0.0, 1.0, 0.0};
    BOOST_CHECK_CLOSE(faceProps.weightAverage(elemNNC1, elemNNC2), weightAvgNNC, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.weightAverage(elemNNC2, elemNNC1), weightAvgNNC, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.weightProduct(elemNNC1, elemNNC2), weightProdNNC, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.weightProduct(elemNNC2, elemNNC1), weightProdNNC, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.normalDistance(elemNNC1, elemNNC2), normDistNNC, 1.0e-8);
    BOOST_CHECK_CLOSE(faceProps.normalDistance(elemNNC2, elemNNC1), normDistNNC, 1.0e-8);
    BOOST_CHECK_EQUAL(faceProps.cellFaceNormal(elemNNC1, elemNNC2), normalNNC_35);
    BOOST_CHECK_EQUAL(faceProps.cellFaceNormal(elemNNC2, elemNNC1), normalNNC_53);

    // Boundary interfaces in cell 7
    const unsigned elemBnd = 7;
    const double normDistBnd = 25.0;
    const double weightAvgBnd = 1.0;
    const double weightProdBnd = 0.0;
    const std::vector<DimVector> normalBnd = { {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0},  // x-dir
                                               {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0},  // y-dir
                                               {0.0, 0.0, 1.0} };  // z-dir
    for (std::size_t bndIdx = 0; bndIdx < 5; ++bndIdx) {
        BOOST_CHECK_CLOSE(faceProps.weightAverageBoundary(elemBnd, bndIdx), weightAvgBnd, 1.0e-8);
        BOOST_CHECK_CLOSE(faceProps.weightProductBoundary(elemBnd, bndIdx), weightProdBnd, 1.0e-8);
        BOOST_CHECK_CLOSE(faceProps.normalDistanceBoundary(elemBnd, bndIdx), normDistBnd, 1.0e-8);
        BOOST_CHECK_EQUAL(faceProps.cellFaceNormalBoundary(elemBnd, bndIdx), normalBnd[bndIdx]);
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
