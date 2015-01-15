/*
  Copyright 2014 Andreas Lauser

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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE TransmissibilityMultipliers

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>

#include <opm/core/grid/GridManager.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#if HAVE_DUNE_CORNERPOINT
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <dune/grid/CpGrid.hpp>
#endif

#include <opm/core/grid/GridHelpers.hpp>
#include <boost/test/unit_test.hpp>

#include <string>

#include <stdlib.h>

// as surprising as it seems, this is a minimal deck required to get to the point where
// the transmissibilities are calculated. The problem here is that the OPM property
// objects mix fluid and rock properties, so that all properties need to be defined even
// if they are not of interest :/
std::string deckPreMult =
    "RUNSPEC\n"
    "TABDIMS\n"
    "/\n"
    "OIL\n"
    "GAS\n"
    "WATER\n"
    "METRIC\n"
    "DIMENS\n"
    "2 2 2/\n"
    "GRID\n"
    "DXV\n"
    "1.0 2.0 /\n"
    "DYV\n"
    "3.0 4.0 /\n"
    "DZV\n"
    "5.0 6.0/\n"
    "TOPS\n"
    "4*100 /\n";
std::string deckPostMult =
    "PROPS\n"
    "DENSITY\n"
    "100 200 300 /\n"
    "PVTW\n"
    " 100 1 1e-6 1.0 0 /\n"
    "PVDG\n"
    "1 1 1e-2\n"
    "100 0.25 2e-2 /\n"
    "PVTO\n"
    "1e-3 1.0 1.0 1.0\n"
    "     100.0 1.0 1.0\n"
    "/\n"
    "1.0 10.0 1.1 0.9\n"
    "    100.0 1.1 0.9\n"
    "/\n"
    "/\n"
    "SWOF\n"
    "0.0 0.0 1.0 0.0\n"
    "1.0 1.0 0.0 1.0/\n"
    "SGOF\n"
    "0.0 0.0 1.0 0.0\n"
    "1.0 1.0 0.0 1.0/\n"
    "PORO\n"
    "8*0.3 /\n"
    "PERMX\n"
    "8*1 /\n"
    "SCHEDULE\n"
    "TSTEP\n"
    "1.0 2.0 3.0 4.0 /\n"
    "/\n";

std::string origDeckString = deckPreMult + deckPostMult;
std::string multDeckString =
    deckPreMult +
    "MULTX\n" +
    "1 2 3 4 5 6 7 8 /\n" +
    "MULTY\n" +
    "1 2 3 4 5 6 7 8 /\n" +
    "MULTZ\n" +
    "1 2 3 4 5 6 7 8 /\n" +
    deckPostMult;

std::string multMinusDeckString =
    deckPreMult +
    "MULTX-\n" +
    "1 2 3 4 5 6 7 8 /\n" +
    "MULTY-\n" +
    "1 2 3 4 5 6 7 8 /\n" +
    "MULTZ-\n" +
    "1 2 3 4 5 6 7 8 /\n" +
    deckPostMult;

// the NTG values get harmonically averaged for the transmissibilites. If the value
// is the same on both sides, the averaging buils down to a no-op, though...
std::string ntgDeckString =
    deckPreMult +
    "NTG\n" +
    "8*0.5 /\n" +
    deckPostMult;


/// \brief The check of the transmissibility values
///
/// Written along the UgGridHelpers interface to make it usable for both
/// UnstructuredGrid and CpGrid.
/// \tparam Grid The type of the grid.
/// \param grid The grid where the transmissibility values live
/// \param origGeology The originaly geology without multipliers.
/// \param multGeology The geology using posititve multipliers.
/// \param multMinusGeology The geology uing negative multipliers
/// \param ntgGeology
template<class G>
void checkTransmissibilityValues(const G&                  grid,
                                 const Opm::DerivedGeology& origGeology,
                                 const Opm::DerivedGeology& multGeology,
                                 const Opm::DerivedGeology& multMinusGeology,
                                 const Opm::DerivedGeology& ntgGeology);


BOOST_AUTO_TEST_CASE(TransmissibilityMultipliersLegacyGridInterface)
{
    Opm::parameter::ParameterGroup param;
    Opm::ParserPtr parser(new Opm::Parser() );

    /////
    // create a DerivedGeology object without any multipliers involved
    Opm::DeckConstPtr origDeck = parser->parseString(origDeckString);
    Opm::EclipseStateConstPtr origEclipseState(new Opm::EclipseState(origDeck));

    auto origGridManager = std::make_shared<Opm::GridManager>(origEclipseState->getEclipseGrid());
    auto origProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(origDeck, origEclipseState, *(origGridManager->c_grid()));

    Opm::DerivedGeology origGeology(*(origGridManager->c_grid()), *origProps, origEclipseState, false);
    /////

    /////
    // create a DerivedGeology object _with_ transmissibility multipliers involved
    Opm::DeckConstPtr multDeck = parser->parseString(multDeckString);
    Opm::EclipseStateConstPtr multEclipseState(new Opm::EclipseState(multDeck));

    auto multGridManager = std::make_shared<Opm::GridManager>(multEclipseState->getEclipseGrid());
    auto multProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(multDeck, multEclipseState, *(multGridManager->c_grid()));

    Opm::DerivedGeology multGeology(*(multGridManager->c_grid()), *multProps, multEclipseState, false);
    /////

    /////
    // create a DerivedGeology object _with_ transmissibility multipliers involved for
    // the negative faces
    Opm::DeckConstPtr multMinusDeck = parser->parseString(multMinusDeckString);
    Opm::EclipseStateConstPtr multMinusEclipseState(new Opm::EclipseState(multMinusDeck));

    auto multMinusGridManager = std::make_shared<Opm::GridManager>(multMinusEclipseState->getEclipseGrid());
    auto multMinusProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(multMinusDeck, multMinusEclipseState, *(multMinusGridManager->c_grid()));

    Opm::DerivedGeology multMinusGeology(*(multMinusGridManager->c_grid()), *multMinusProps, multMinusEclipseState, false);
    /////

    /////
    // create a DerivedGeology object with the NTG keyword involved
    Opm::DeckConstPtr ntgDeck = parser->parseString(ntgDeckString);
    Opm::EclipseStateConstPtr ntgEclipseState(new Opm::EclipseState(ntgDeck));

    auto ntgGridManager = std::make_shared<Opm::GridManager>(ntgEclipseState->getEclipseGrid());
    auto ntgProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(ntgDeck, ntgEclipseState, *(ntgGridManager->c_grid()));

    Opm::DerivedGeology ntgGeology(*(ntgGridManager->c_grid()), *ntgProps, ntgEclipseState, false);
    /////

    checkTransmissibilityValues(*(origGridManager->c_grid()), origGeology, multGeology,
                                multMinusGeology, ntgGeology);
}

template<class G>
void checkTransmissibilityValues(const G&                  grid,
                                 const Opm::DerivedGeology& origGeology,
                                 const Opm::DerivedGeology& multGeology,
                                 const Opm::DerivedGeology& multMinusGeology,
                                 const Opm::DerivedGeology& ntgGeology)
{
    // compare the transmissibilities (note that for this we assume that the multipliers
    // do not change the grid topology)
    int numFaces=Opm::UgGridHelpers::numFaces(grid);
    typename Opm::UgGridHelpers::FaceCellTraits<G>::Type face_cells=
        Opm::UgGridHelpers::faceCells(grid);
    for (int faceIdx = 0; faceIdx < numFaces; ++ faceIdx) {
        // in DUNE-speak, a face here is more like an intersection which is not specific
        // to a codim-0 entity (i.e., cell)

        // get the cell indices of the compressed grid for the face's interior and
        // exterior cell
        int insideCellIdx  = face_cells(faceIdx, 0);
        int outsideCellIdx = face_cells(faceIdx, 1);

        if (insideCellIdx < 0 || outsideCellIdx < 0) {
            // do not consider cells at the domain boundary: Their would only be used for
            // Dirichlet-like boundary conditions which have not been implemented so
            // far...
            continue;
        }

        // translate these to canonical indices (i.e., the logically Cartesian ones used by the deck)
        const int* global_cell=Opm::UgGridHelpers::globalCell(grid);

        if (global_cell) {
            insideCellIdx = global_cell[insideCellIdx];
            outsideCellIdx = global_cell[outsideCellIdx];
        }

        double origTrans = origGeology.transmissibility()[faceIdx];
        double multTrans = multGeology.transmissibility()[faceIdx];
        double multMinusTrans = multMinusGeology.transmissibility()[faceIdx];
        double ntgTrans = ntgGeology.transmissibility()[faceIdx];
        BOOST_CHECK_CLOSE(origTrans*(insideCellIdx + 1), multTrans, 1e-6);
        BOOST_CHECK_CLOSE(origTrans*(outsideCellIdx + 1), multMinusTrans, 1e-6);

        const int* cartdims = Opm::UgGridHelpers::cartDims(grid);
        int insideCellKIdx = insideCellIdx/(cartdims[0]*cartdims[1]);
        int outsideCellKIdx = outsideCellIdx/(cartdims[0]*cartdims[1]);
        if (insideCellKIdx == outsideCellKIdx)
            // NTG only reduces the permebility of the X-Y plane
            BOOST_CHECK_CLOSE(origTrans*0.5, ntgTrans, 1e-6);
    }
}

#if HAVE_DUNE_CORNERPOINT
BOOST_AUTO_TEST_CASE(TransmissibilityMultipliersCpGrid)
{
    int argc = 1;
    char **argv;
    argv = new (char*);
    argv[0] = strdup("footest");

    Dune::MPIHelper::instance(argc, argv);

    Opm::parameter::ParameterGroup param;
    Opm::ParserPtr parser(new Opm::Parser() );

    /////
    // create a DerivedGeology object without any multipliers involved
    Opm::DeckConstPtr origDeck = parser->parseString(origDeckString);
    Opm::EclipseStateConstPtr origEclipseState(new Opm::EclipseState(origDeck));

    auto origGrid = std::make_shared<Dune::CpGrid>();
    origGrid->processEclipseFormat(origEclipseState->getEclipseGrid(), 0.0, false);

    auto origProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(origDeck,
                                                                    origEclipseState,
                                                                    *origGrid);

    Opm::DerivedGeology origGeology(*origGrid, *origProps, origEclipseState, false);
    /////

    /////
    // create a DerivedGeology object _with_ transmissibility multipliers involved
    Opm::DeckConstPtr multDeck = parser->parseString(multDeckString);
    Opm::EclipseStateConstPtr multEclipseState(new Opm::EclipseState(multDeck));

    auto multGrid = std::make_shared<Dune::CpGrid>();
    multGrid->processEclipseFormat(multEclipseState->getEclipseGrid(), 0.0, false);

    auto multProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(multDeck, multEclipseState, *multGrid);

    Opm::DerivedGeology multGeology(*multGrid, *multProps, multEclipseState, false);
    /////

    /////
    // create a DerivedGeology object _with_ transmissibility multipliers involved for
    // the negative faces
    Opm::DeckConstPtr multMinusDeck = parser->parseString(multMinusDeckString);
    Opm::EclipseStateConstPtr multMinusEclipseState(new Opm::EclipseState(multMinusDeck));

    auto multMinusGrid = std::make_shared<Dune::CpGrid>();
    multMinusGrid->processEclipseFormat(multMinusEclipseState->getEclipseGrid(), 0.0, false);

    auto multMinusProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(multMinusDeck, multMinusEclipseState, *multMinusGrid);

    Opm::DerivedGeology multMinusGeology(*multMinusGrid, *multMinusProps, multMinusEclipseState,
                                         false);
    /////

    /////
    // create a DerivedGeology object with the NTG keyword involved
    Opm::DeckConstPtr ntgDeck = parser->parseString(ntgDeckString);
    Opm::EclipseStateConstPtr ntgEclipseState(new Opm::EclipseState(ntgDeck));

    auto ntgGrid = std::make_shared<Dune::CpGrid>();
    ntgGrid->processEclipseFormat(ntgEclipseState->getEclipseGrid(), 0.0, false);

    auto ntgProps = std::make_shared<Opm::BlackoilPropsAdFromDeck>(ntgDeck, ntgEclipseState, *ntgGrid);

    Opm::DerivedGeology ntgGeology(*ntgGrid, *ntgProps, ntgEclipseState, false);
    /////

    return checkTransmissibilityValues(*origGrid, origGeology, multGeology,
                                       multMinusGeology, ntgGeology);
}

#endif
