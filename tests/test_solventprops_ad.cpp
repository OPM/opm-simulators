/*
  Copyright 2016 IRIS AS

  This file is part of the Open Porous Media Project (OPM).

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

#define BOOST_TEST_MODULE SolventPropertiesTest

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/autodiff/SolventPropsAdFromDeck.hpp>

#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <fstream>
#include <iostream>


const std::string deckData = "\n\
RUNSPEC \n\
\n\
SOLVENT \n\
\n\
MISCIBLE\n\
1  3 /\n\
\n\
DIMENS \n\
1 1 1 \n\
/\n\
TABDIMS\n\
/\n\
GRID \n\
\n\
DXV \n\
1 \n\
/\n\
DYV \n\
1 \n\
/\n\
DZV \n\
1 \n\
/\n";

const std::string solventData = "\n\
SDENSITY \n\
0.1 / \n\
PVDS \n\
1 1 0.1 / \n\
SSFN \n\
0.0 0.0 0.0 \n\
1.0 1.0 1.0 \n\
/ \n\
MISC \n\
0.0 0.0 \n\
1.0 1.0 \n\
/ \n\
SOF2 \n\
0 0  \n\
0.88 1 / \n";

BOOST_AUTO_TEST_CASE(Construction)
{
        Opm::ParseContext parseContext;
        Opm::ParserPtr parser(new Opm::Parser());
        Opm::DeckPtr deck =  parser->parseString(deckData + solventData, parseContext);
        Opm::EclipseStateConstPtr eclState;
        eclState.reset(new Opm::EclipseState(deck , parseContext));
        const int global_ind = 0;
        Opm::SolventPropsAdFromDeck solventprops(deck, eclState, 1, &global_ind);
}

BOOST_AUTO_TEST_CASE(SolventData)
{
        Opm::ParseContext parseContext;
        Opm::ParserPtr parser(new Opm::Parser());
        Opm::DeckPtr deck =  parser->parseString(deckData + solventData, parseContext);
        Opm::EclipseStateConstPtr eclState;
        eclState.reset(new Opm::EclipseState(deck , parseContext));
        const int global_ind = 0;
        Opm::SolventPropsAdFromDeck solventprops(deck, eclState, 1, &global_ind);

        const Opm::SolventPropsAdFromDeck::Cells cells(1, 0);
        typedef Opm::SolventPropsAdFromDeck::V V;
        V rho = solventprops.solventSurfaceDensity(cells);
        BOOST_REQUIRE_EQUAL(rho.size(), cells.size());
        BOOST_CHECK_EQUAL(rho[0], 0.1);
}

const std::string pmiscData = "\n\
PMISC\n\
100 0.0 \n\
200 0.0 \n\
500 1.0 \n\
1000 1.0 /\n\
\n";

BOOST_AUTO_TEST_CASE(PMISC)
{
        Opm::ParseContext parseContext;
        Opm::ParserPtr parser(new Opm::Parser());
        Opm::DeckPtr deck =  parser->parseString(deckData + solventData + pmiscData, parseContext);
        Opm::EclipseStateConstPtr eclState;
        eclState.reset(new Opm::EclipseState(deck , parseContext));
        const Opm::SolventPropsAdFromDeck::Cells cells(3, 0);
        typedef Opm::SolventPropsAdFromDeck::V V;
        const int* global_ind = new int[3] {0 , 1 , 2};
        Opm::SolventPropsAdFromDeck solventprops(deck, eclState, 3, global_ind);
        V po(3);
        po << 150,250,550;
        po = po * Opm::unit::barsa;
        BOOST_REQUIRE_EQUAL(po.size(), cells.size());
        V pmisc = solventprops.pressureMiscibilityFunction(Opm::SolventPropsAdFromDeck::ADB::constant(po), cells).value();
        BOOST_REQUIRE_EQUAL(pmisc.size(), cells.size());
        BOOST_CHECK_EQUAL(pmisc[0], 0.0);
        const double tol = 1e-12;
        const double value = (250.0 - 200.0) / (500.0 - 200.0); // linear interpolation
        BOOST_CHECK_CLOSE(pmisc[1], value, tol);
        BOOST_CHECK_EQUAL(pmisc[2], 1.0);
}

const std::string tlpmixpaData = "\n\
TLPMIXPA\n\
100 0.0 \n\
200 0.0 \n\
500 1.0 \n\
1000 1.0 /\n\
\n";


BOOST_AUTO_TEST_CASE(TLPMIXPA)
{
        Opm::ParseContext parseContext;
        Opm::ParserPtr parser(new Opm::Parser());
        Opm::DeckPtr deck =  parser->parseString(deckData + solventData + tlpmixpaData, parseContext);
        Opm::EclipseStateConstPtr eclState;
        eclState.reset(new Opm::EclipseState(deck , parseContext));
        const Opm::SolventPropsAdFromDeck::Cells cells(3, 0);
        typedef Opm::SolventPropsAdFromDeck::V V;
        const int* global_ind = new int[3] {0 , 1 , 2};
        Opm::SolventPropsAdFromDeck solventprops(deck, eclState, 3, global_ind);
        V po(3);
        po << 150,250,550;
        po = po * Opm::unit::barsa;
        BOOST_REQUIRE_EQUAL(po.size(), cells.size());
        V tlpmixpa = solventprops.pressureMixingParameter(Opm::SolventPropsAdFromDeck::ADB::constant(po), cells).value();
        BOOST_REQUIRE_EQUAL(tlpmixpa.size(), cells.size());
        BOOST_CHECK_EQUAL(tlpmixpa[0], 0.0);
        const double tol = 1e-12;
        const double value = (250.0 - 200.0) / (500.0 - 200.0); // linear interpolation
        BOOST_CHECK_CLOSE(tlpmixpa[1], value, tol);
        BOOST_CHECK_EQUAL(tlpmixpa[2], 1.0);
}

BOOST_AUTO_TEST_CASE(TLPMIXPA_NOT_SPECIFIED)
{
    Opm::ParseContext parseContext;
    Opm::ParserPtr parser(new Opm::Parser());
    // no pmisc data and default tlpmixdata i.e it should throw
    Opm::DeckPtr deck =  parser->parseString(deckData + solventData, parseContext);
    Opm::EclipseStateConstPtr eclState;
    eclState.reset(new Opm::EclipseState(deck , parseContext));
    const Opm::SolventPropsAdFromDeck::Cells cells(3, 0);
    const int* global_ind = new int[3] {0 , 1 , 2};
    Opm::SolventPropsAdFromDeck solventprops(deck, eclState, 3, global_ind);
    typedef Opm::SolventPropsAdFromDeck::V V;
    V po(3);
    po << 150,250,550;
    po = po * Opm::unit::barsa;
    BOOST_REQUIRE_EQUAL(po.size(), cells.size());
    V tlpmixpa = solventprops.pressureMixingParameter(Opm::SolventPropsAdFromDeck::ADB::constant(po), cells).value();
    BOOST_REQUIRE_EQUAL(tlpmixpa.size(), cells.size());
    // if not specified tlpmixpa is 1.0 for all cells.
    BOOST_CHECK_EQUAL(tlpmixpa[0], 1.0);
    BOOST_CHECK_EQUAL(tlpmixpa[1], 1.0);
    BOOST_CHECK_EQUAL(tlpmixpa[2], 1.0);

}

const std::string tlpmixpaDataDefault = "\n\
TLPMIXPA\n\
/\n\
\n";

BOOST_AUTO_TEST_CASE(TLPMIXPA_DEFAULT)
{
    Opm::ParseContext parseContext;
    Opm::ParserPtr parser(new Opm::Parser());
    Opm::DeckPtr deck =  parser->parseString(deckData + solventData + pmiscData + tlpmixpaDataDefault, parseContext);
    Opm::EclipseStateConstPtr eclState;
    eclState.reset(new Opm::EclipseState(deck , parseContext));
    const Opm::SolventPropsAdFromDeck::Cells cells(3, 0);
    typedef Opm::SolventPropsAdFromDeck::V V;
    const int* global_ind = new int[3] {0 , 1 , 2};
    Opm::SolventPropsAdFromDeck solventprops(deck, eclState, 3, global_ind);
    V po(3);
    po << 150,250,550;
    po = po * Opm::unit::barsa;
    BOOST_REQUIRE_EQUAL(po.size(), cells.size());
    V tlpmixpa = solventprops.pressureMixingParameter(Opm::SolventPropsAdFromDeck::ADB::constant(po), cells).value();
    BOOST_REQUIRE_EQUAL(tlpmixpa.size(), cells.size());
    BOOST_CHECK_EQUAL(tlpmixpa[0], 0.0);
    const double tol = 1e-12;
    const double value = (250.0 - 200.0) / (500.0 - 200.0); // linear interpolation
    BOOST_CHECK_CLOSE(tlpmixpa[1], value, tol);
    BOOST_CHECK_EQUAL(tlpmixpa[2], 1.0);
}

BOOST_AUTO_TEST_CASE(TLPMIXPA_DEFAULT_NOPMISC)
{
    Opm::ParseContext parseContext;
    Opm::ParserPtr parser(new Opm::Parser());
    // no pmisc data and default tlpmixdata i.e it should throw
    Opm::DeckPtr deck =  parser->parseString(deckData + solventData + tlpmixpaDataDefault, parseContext);
    Opm::EclipseStateConstPtr eclState;
    eclState.reset(new Opm::EclipseState(deck , parseContext));
    const Opm::SolventPropsAdFromDeck::Cells cells(3, 0);
    const int* global_ind = new int[3] {0 , 1 , 2};
    BOOST_CHECK_THROW(Opm::SolventPropsAdFromDeck solventprops(deck, eclState, 3, global_ind), std::invalid_argument);
}

