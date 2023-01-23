/*
  Copyright 2015 Statoil ASA.

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

#define BOOST_TEST_MODULE NORNE_PVT_TESTS

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <sstream>
#include <iostream>

#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>


using namespace Opm;

/*
  This file contains a regression test of the LiveOilPvt class for the
  Norne data. The test is created by calling the PvtLiveOil::mu( )
  function for a list of P and Rs values and storing the results in
  the XX_expected vectors. No actual validation of the data has been
  done, so the test will only serve to further cement possible bugs.

  The columns of input data and expected data is organized in two
  columns, that is mainly because two (P,Rs) values were selected from
  each of the 'inner tables' in the PVTO keyword, and carries no
  further semantic meaning.
*/

void verify_norne_oil_pvt_region1(const Opm::EclipseState& eclState, const Opm::Schedule& schedule) {
    Opm::LiveOilPvt<double> oilPvt;
    oilPvt.initFromState(eclState, schedule);

    std::vector<double> rs = {33, 33,
                              43, 43,
                              53, 53,
                              61, 61,
                              70, 70,
                              80, 80,
                              100, 100 ,
                              100};


    std::vector<double> P = {114, 148,
                             134, 168,
                             154, 188,
                             174, 208,
                             194, 228,
                             214, 248,
                             234, 268,
                             270 };


    std::vector<double> mu_expected = {0.00106736588,   0.00113961037,
                                       0.00093801366,   0.00099871729,
                                       0.00083529743,   0.00088728769,
                                       0.00077986989,   0.00082627508,
                                       0.00072883113,   0.00076988665,
                                       0.00068250424,   0.00072040786,
                                       0.00062347677,   0.00064963306,
                                       0.00065122911};


    std::vector<double> b_expected = {0.88421444595,   0.88893909117,
                                      0.86493342861,   0.86978957420,
                                      0.84676402016,   0.85171762998,
                                      0.83354279748,   0.83851861429,
                                      0.81904041272,   0.82404719615,
                                      0.80341044483,   0.80845950744,
                                      0.77131381726,   0.77661604334,
                                      0.77691738473};

    {
        std::vector<int> tableIndex(P.size() , 0);

        // convert the pressures to SI units (bar to Pascal)
        for (auto& value : P)
            value *= Metric::Pressure;

        // convert the gas dissolution factors to SI units
        for (auto& value : rs)
            value *=  Metric::GasDissolutionFactor;

        for (unsigned i = 0; i < P.size(); ++i) {
            double mu;
            double b;
            double RsSat = oilPvt.saturatedGasDissolutionFactor(/*tableIndex=*/0, /*T=*/273.15, P[i]);
            if (rs[i] >= RsSat) {
                mu = oilPvt.saturatedViscosity(/*tableIndex=*/0, /*T=*/273.15, P[i]);
                b = oilPvt.saturatedInverseFormationVolumeFactor(/*tableIndex=*/0, /*T=*/273.15, P[i]);
            }
            else {
                mu = oilPvt.viscosity(/*tableIndex=*/0, /*T=*/273.15, P[i], rs[i]);
                b = oilPvt.inverseFormationVolumeFactor(/*tableIndex=*/0, /*T=*/273.15, P[i], rs[i]);
            }

            BOOST_CHECK_CLOSE( mu , mu_expected[i], 1e-5 );
            BOOST_CHECK_CLOSE( b , b_expected[i], 1e-5 );
        }
    }
}


void verify_norne_oil_pvt_region2(const Opm::EclipseState& eclState, const Opm::Schedule& schedule) {
    Opm::LiveOilPvt<double> oilPvt;
    oilPvt.initFromState(eclState, schedule);

    std::vector<double> rs = {21 , 21,
                              30 , 30,
                              38 , 38,
                              48 , 48,
                              55 , 55,
                              65 , 65,
                              75 , 75,
                              85 , 85,
                              95 , 95,
                              105 , 105,
                              115 , 115,
                              125 , 125,
                              135 , 135,
                              145 , 145,
                              155 , 155,
                              165 , 165,
                              175 , 175,
                              185 , 185,
                              195 , 195,
                              205 , 205,
                              215 , 215,
                              225 , 225,
                              234 , 234,
                              240 , 240,
                              252 , 252,
                              262 , 262,
                              272 , 272,
                              280 , 280,
                              410, 410, 410};


    std::vector<double> P = {70,  110,
                             95,  145,
                             115, 165,
                             135, 185,
                             155, 205,
                             195, 245,
                             215, 265,
                             235, 285,
                             255, 305,
                             275, 325,
                             293, 343,
                             310, 360,
                             326, 376,
                             342, 392,
                             357, 407,
                             371, 420,
                             385, 435,
                             399, 450,
                             420, 480,
                             437, 487,
                             449, 499,
                             460, 510,
                             471, 521,
                             482, 532,
                             503, 553,
                             650, 680, 710};


    std::vector<double> mu_expected = {0.00120767750,   0.00129077352,
                                       0.00111063039,   0.00119627038,
                                       0.00103118116,   0.00110633521,
                                       0.00094413471,   0.00100998373,
                                       0.00090320931,   0.00096374536,
                                       0.00086714481,   0.00092142974,
                                       0.00081811098,   0.00086735227,
                                       0.00077704364,   0.00082229010,
                                       0.00070975205,   0.00076029164,
                                       0.00065679329,   0.00071124175,
                                       0.00061496175,   0.00067213642,
                                       0.00058000381,   0.00064115346,
                                       0.00055124739,   0.00061633274,
                                       0.00052840888,   0.00059781928,
                                       0.00050926184,   0.00058323394,
                                       0.00049295739,   0.00056996321,
                                       0.00048026810,   0.00056474486,
                                       0.00047088998,   0.00056427878,
                                       0.00047649659,   0.00060774836,
                                       0.00048006188,   0.00059909192,
                                       0.00026623648,   0.00060915386,
                                       0.00025670489,   0.00062157315,
                                       0.00024760210,   0.00064290735,
                                       0.00023889979,   0.00067946283,
                                       0.00022330662,   0.00077837223,
                                       0.01142273040,  -0.00351292519,  -0.00129867195};

    std::vector<double> b_expected = {0.90699449462,   0.91120449633,
                                      0.89040695696,   0.89551008140,
                                      0.87548859167,   0.88062965205,
                                      0.85697013389,   0.86224235632,
                                      0.84533618728,   0.85061301709,
                                      0.83069819286,   0.83585867335,
                                      0.81473536808,   0.81994107210,
                                      0.79955491390,   0.80479144821,
                                      0.78507711370,   0.79032915313,
                                      0.77073097762,   0.77596189361,
                                      0.75627401890,   0.76141290296,
                                      0.74161331648,   0.74678198081,
                                      0.72686889575,   0.73206734035,
                                      0.71214353439,   0.71737175926,
                                      0.69733207231,   0.70259007745,
                                      0.68243272267,   0.68761475238,
                                      0.66755004999,   0.67286761567,
                                      0.65268405426,   0.65813834713,
                                      0.63858753316,   0.64504008462,
                                      0.62408347496,   0.62949038145,
                                      0.61223874629,   0.61449268543,
                                      0.60422344638,   0.59939995459,
                                      0.59620814647,   0.58594855211,
                                      0.58819284656,   0.57739165219,
                                      0.57289091037,   0.56019050084,
                                      0.55474601877,   0.55809201119,   0.54526832277};

    // convert the pressures to SI units (bar to Pascal)
    for (auto& value : P)
        value *= Metric::Pressure;

    // convert the gas dissolution factors to SI units
    for (auto& value : rs)
        value *=  Metric::GasDissolutionFactor;

    for (unsigned i = 0; i < P.size(); ++i) {
        double mu;
        double b;
        double RsSat = oilPvt.saturatedGasDissolutionFactor(/*tableIndex=*/1, /*T=*/273.15, P[i]);
        if (rs[i] >= RsSat) {
            mu = oilPvt.saturatedViscosity(/*tableIndex=*/1, /*T=*/273.15, P[i]);
            b = oilPvt.saturatedInverseFormationVolumeFactor(/*tableIndex=*/1, /*T=*/273.15, P[i]);
        }
        else {
            mu = oilPvt.viscosity(/*tableIndex=*/1, /*T=*/273.15, P[i], rs[i]);
            b = oilPvt.inverseFormationVolumeFactor(/*tableIndex=*/1, /*T=*/273.15, P[i], rs[i]);
        }

        BOOST_CHECK_CLOSE( mu , mu_expected[i], 1e-5 );
        BOOST_CHECK_CLOSE( b , b_expected[i], 1e-5 );
    }
}

BOOST_AUTO_TEST_CASE( Test_Norne_PVT) {
    Opm::ParseContext parseContext({{ ParseContext::PARSE_RANDOM_SLASH , InputErrorAction::IGNORE }});
    Opm::ErrorGuard errorGuard;
    Opm::Parser parser;
    auto python = std::make_shared<Opm::Python>();

    auto deck = parser.parseFile("norne_pvt.data", parseContext, errorGuard);

    Opm::EclipseState eclState(deck);
    Opm::Schedule schedule(deck, eclState, python);

    verify_norne_oil_pvt_region1(eclState, schedule);
    verify_norne_oil_pvt_region2(eclState, schedule);
}
