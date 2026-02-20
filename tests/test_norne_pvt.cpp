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

#include <algorithm>
#include <fstream>
#include <iostream>

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


// The TEST_OR_PRINT macro has been added to simplify updating this
// regression test.

#define TEST_OR_PRINT BOOST_CHECK_CLOSE
// #define TEST_OR_PRINT(val, expected, tolerance) printVal(val);

inline void printVal(const double val)
{
    std::cout.precision(16);
    std::cout << val << '\n';
}


void
verify_norne_oil_pvt_region1(const Opm::EclipseState& eclState,
                             const Opm::Schedule& schedule,
                             std::istream& expectedData)
{
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

    std::vector<double> mu_expected(P.size());
    std::vector<double> b_expected(P.size());
    for (std::size_t ii = 0; ii < P.size(); ++ii) {
        expectedData >> mu_expected[ii] >> b_expected[ii];
    }

    {
        // convert the pressures to SI units (bar to Pascal)
        std::ranges::transform(P, P.begin(),
                               [](const auto value) { return value * Metric::Pressure; });

        // convert the gas dissolution factors to SI units
        std::ranges::transform(rs, rs.begin(),
                               [](const auto value) { return value * Metric::GasDissolutionFactor; });

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

            TEST_OR_PRINT( mu , mu_expected[i], 1e-5 );
            TEST_OR_PRINT( b , b_expected[i], 1e-5 );
        }
    }
}


void
verify_norne_oil_pvt_region2(const Opm::EclipseState& eclState,
                             const Opm::Schedule& schedule,
                             std::istream& expectedData)
{
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

    std::vector<double> mu_expected(P.size());
    std::vector<double> b_expected(P.size());
    for (std::size_t ii = 0; ii < P.size(); ++ii) {
        expectedData >> mu_expected[ii] >> b_expected[ii];
    }


    // convert the pressures to SI units (bar to Pascal)
    std::ranges::transform(P, P.begin(),
                           [](const auto value) { return value * Metric::Pressure; });

    // convert the gas dissolution factors to SI units
    std::ranges::transform(rs, rs.begin(),
                           [](const auto value) { return value * Metric::GasDissolutionFactor; });

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

        TEST_OR_PRINT( mu , mu_expected[i], 1e-5 );
        TEST_OR_PRINT( b , b_expected[i], 1e-5 );
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

    // The expected data are stored in the order they are printed (if
    // setting TEST_OR_PRINT to printing), i.e. mu and b for each
    // sample point in order.
    std::ifstream expectedData("norne_pvt_expected.txt");
    BOOST_REQUIRE(expectedData);

    verify_norne_oil_pvt_region1(eclState, schedule, expectedData);
    verify_norne_oil_pvt_region2(eclState, schedule, expectedData);
}
