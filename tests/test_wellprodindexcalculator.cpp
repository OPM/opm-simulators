/*
  Copyright 2020 Equinor.

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

#define BOOST_TEST_MODULE TestWellProdIndexCalculator

#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/WellProdIndexCalculator.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace {
    double liquid_PI_unit()
    {
        return Opm::UnitSystem::newMETRIC().to_si(Opm::UnitSystem::measure::liquid_productivity_index, 1.0);
    }

    double cp_rm3_per_db()
    {
        return Opm::UnitSystem::newMETRIC().to_si(Opm::UnitSystem::measure::transmissibility, 1.0);
    }

    std::string drainRadDefaulted()
    {
        return { R"(
  'P' 'G' 10 10 2005 'LIQ' /
)" };
    }

    std::string explicitDrainRad()
    {
        // rd = exp(2)
        return { R"(
  'P' 'G' 10 10 2005 'LIQ' 7.38905609893065 /
)" };
    }

    std::string noSkinFactor_SameCF()
    {
        // r0 = exp(1)
        return { R"(
  'P' 0 0 1 3 OPEN 1 100 2.0 4* 2.718281828459045 /
)" };
    }

    std::string noSkinFactor_DifferentCF()
    {
        // r0 = exp(1)
        return { R"(
  'P' 0 0 1 1 OPEN 1  50 2.0 4* 2.718281828459045 /
  'P' 0 0 2 2 OPEN 1 100 2.0 4* 2.718281828459045 /
  'P' 0 0 3 3 OPEN 1 200 2.0 4* 2.718281828459045 /
)" };
    }

    std::string skin2_SameCF()
    {
        // r0 = exp(1), Skin = 2
        return { R"(
  'P' 0 0 1 3 OPEN 1 100 2.0 1* 2.0 2* 2.718281828459045 /
)" };
    }

    std::string skin421_DifferentCF()
    {
        // r0 = exp(1), Skin = 4, 2, 1
        return { R"(
  'P' 0 0 1 1 OPEN 1  50 2.0 1* 4.0 2* 2.718281828459045 /
  'P' 0 0 2 2 OPEN 1 100 2.0 1* 2.0 2* 2.718281828459045 /
  'P' 0 0 3 3 OPEN 1 200 2.0 1* 1.0 2* 2.718281828459045 /
)" };
    }

    Opm::Well createWell(const std::string& welspecs,
                         const std::string& compdat)
    {
        const auto deck = Opm::Parser{}.parseString(R"(RUNSPEC
DIMENS
  10 10 3 /
START
 8 OCT 2020 /
GRID
DXV
  10*100.0 /
DYV
  10*100.0 /
DZV
  3*10.0 /
DEPTHZ
  121*2000.0 /
PERMX
  300*100.0 /
PERMY
  300*100.0 /
PERMZ
  300*10.0 /
PORO
  300*0.3 /
SCHEDULE
WELSPECS
)" + welspecs + R"(
/

COMPDAT
)" + compdat + R"(
/

TSTEP
  10
/
END
)");

        const auto es    = Opm::EclipseState{ deck };
        const auto sched = Opm::Schedule{ deck, es };

        return sched.getWell("P", 0);
    }
}

BOOST_AUTO_TEST_SUITE(ConnectionLevel)

BOOST_AUTO_TEST_CASE(allDefaulted_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), noSkinFactor_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 2.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 4.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(allDefaulted_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), noSkinFactor_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 2.0), expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 1.0), expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 0.5), expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_Skin2_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), skin2_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 2.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 4.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_skin421_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), skin421_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 2.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 1.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 0.5), 1.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), noSkinFactor_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 0.5 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 2.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), noSkinFactor_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 0.25 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.0  * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 4.0  * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_Skin2_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), skin2_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 0.75 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.5  * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 3.0  * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_skin421_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), skin421_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 100*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), (5.0 / 6.0) * 0.5 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.5         * 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), (8.0 / 3.0) * 2.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END() // ConnectionLevel

// ===========================================================================

BOOST_AUTO_TEST_SUITE(AllConnections)

BOOST_AUTO_TEST_CASE(allDefaulted_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), noSkinFactor_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        1.0*expectCF, 2.0*expectCF, 4.0*expectCF
    };

    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_CASE(allDefaulted_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), noSkinFactor_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 2.0, 1.0, 0.5 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        expectCF, expectCF, expectCF
    };

    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_Skin2_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), skin2_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        1.0*expectCF, 2.0*expectCF, 4.0*expectCF
    };

    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_skin421_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), skin421_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 2.0, 1.0, 0.5 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        expectCF, expectCF, expectCF
    };

    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), noSkinFactor_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        0.5*expectCF, 1.0*expectCF, 2.0*expectCF
    };

    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), noSkinFactor_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        0.25*expectCF, 1.0*expectCF, 4.0*expectCF
    };

    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_Skin2_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), skin2_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        0.75*expectCF, 1.5*expectCF, 3.0*expectCF
    };

    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_skin421_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), skin421_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = std::vector<double> {
        (5.0 / 6.0) * 0.5 * expectCF,
        1.5         * 1.0 * expectCF,
        (8.0 / 3.0) * 2.0 * expectCF,
    };
    
    const auto connPI = connectionProdIndStandard(wpiCalc, connMobility);
    BOOST_CHECK_CLOSE(connPI[0], expectPI[0], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[1], expectPI[1], 1.0e-10);
    BOOST_CHECK_CLOSE(connPI[2], expectPI[2], 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END() // AllConnections

// ===========================================================================

BOOST_AUTO_TEST_SUITE(WellLevel)

BOOST_AUTO_TEST_CASE(allDefaulted_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), noSkinFactor_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = (1.0 + 2.0 + 4.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(allDefaulted_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), noSkinFactor_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 2.0, 1.0, 0.5 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = (1.0 + 1.0 + 1.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_Skin2_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), skin2_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = (1.0 + 2.0 + 4.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_skin421_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(drainRadDefaulted(), skin421_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 2.0, 1.0, 0.5 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = (1.0 + 1.0 + 1.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), noSkinFactor_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = (0.5 + 1.0 + 2.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), noSkinFactor_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = (0.25 + 1.0 + 4.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_Skin2_SameCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), skin2_SameCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = (0.75 + 1.5 + 3.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_skin421_DifferentCF)
{
    const auto wpiCalc = Opm::WellProdIndexCalculator {
        createWell(explicitDrainRad(), skin421_DifferentCF())
    };

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto connMobility = std::vector<double> { 1.0, 2.0, 4.0 };

    const auto expectCF = 100*cp_rm3_per_db();
    const auto expectPI = ((5.0 / 6.0) * 0.5 +
                           1.5         * 1.0 +
                           (8.0 / 3.0) * 2.0)*expectCF;

    BOOST_CHECK_CLOSE(wellProdIndStandard(wpiCalc, connMobility), expectPI, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END() // WellLevel

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Re_Init_Connection_Level)

BOOST_AUTO_TEST_CASE(allDefaulted_SameCF)
{
    auto well = createWell(drainRadDefaulted(), noSkinFactor_SameCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 2.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 4.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(allDefaulted_DifferentCF)
{
    auto well = createWell(drainRadDefaulted(), noSkinFactor_DifferentCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 2.0), expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 1.0), expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 0.5), expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_Skin2_SameCF)
{
    auto well = createWell(drainRadDefaulted(), skin2_SameCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 2.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 4.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(defaultedDRad_skin421_DifferentCF)
{
    auto well = createWell(drainRadDefaulted(), skin421_DifferentCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 2.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 1.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 0.5), 1.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_SameCF)
{
    auto well = createWell(explicitDrainRad(), noSkinFactor_SameCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 0.5 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 2.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_DifferentCF)
{
    auto well = createWell(explicitDrainRad(), noSkinFactor_DifferentCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 0.25 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.0  * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 4.0  * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_Skin2_SameCF)
{
    auto well = createWell(explicitDrainRad(), skin2_SameCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), 0.75 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.5  * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), 3.0  * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_CASE(logarithmic_skin421_DifferentCF)
{
    auto well = createWell(explicitDrainRad(), skin421_DifferentCF());
    auto wpiCalc = Opm::WellProdIndexCalculator { well };

    well.updateWellProductivityIndex( );
    const auto scalingFactor = well.convertDeckPI(2.0) / (1.0*liquid_PI_unit());

    BOOST_CHECK_CLOSE(scalingFactor, 2.0, 1.0e-10);

    std::vector<bool> scalingApplicable;
    well.applyWellProdIndexScaling(scalingFactor, scalingApplicable);

    wpiCalc.reInit(well);

    BOOST_REQUIRE_EQUAL(wpiCalc.numConnections(), std::size_t{3});

    const auto expectCF = 200*cp_rm3_per_db();
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(0, 1.0), (5.0 / 6.0) * 0.5 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(1, 2.0), 1.5         * 1.0 * expectCF, 1.0e-10);
    BOOST_CHECK_CLOSE(wpiCalc.connectionProdIndStandard(2, 4.0), (8.0 / 3.0) * 2.0 * expectCF, 1.0e-10);
}

BOOST_AUTO_TEST_SUITE_END() // Re_Init_Connection_Level
