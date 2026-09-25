/*
  Copyright 2026 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM. If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#define BOOST_TEST_MODULE ExtboContainer
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/ExtboContainer.hpp>

#include <opm/common/OpmLog/CounterLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <map>
#include <memory>
#include <string>

namespace {
using Container = Opm::ExtboContainer<double>;

// Oil carries 100 kg of solvent per 950 kg (2/19) and gas 1 kg per 2.3 kg
// (10/23).
Container::PhaseFractionInput twoPhaseCell()
{
    return {
        .oilSaturation = 0.5,
        .gasSaturation = 0.3,
        .rs = 100.0,
        .rv = 1.0e-3,
        .xVolume = 0.5,
        .yVolume = 0.5,
        .oilDensity = 800.0,
        .gasDensity = 1.0,
        .solventDensity = 2.0,
    };
}

std::map<std::string, int> requests(const int oil, const int gas)
{
    return {{"SOLVMFO", oil}, {"SOLVMFG", gas}};
}
} // anonymous namespace

BOOST_AUTO_TEST_CASE(PhaseSolventMassFractions)
{
    const auto [oil, gas] = Container::phaseSolventMassFractions(twoPhaseCell());
    BOOST_CHECK_CLOSE(oil, 2.0 / 19.0, 1.0e-12);
    BOOST_CHECK_CLOSE(gas, 10.0 / 23.0, 1.0e-12);
}

BOOST_AUTO_TEST_CASE(AbsentPhaseHasZeroFraction)
{
    auto noOil = twoPhaseCell();
    noOil.oilSaturation = 0.0;
    const auto [oil1, gas1] = Container::phaseSolventMassFractions(noOil);
    BOOST_CHECK_EQUAL(oil1, 0.0);
    BOOST_CHECK_CLOSE(gas1, 10.0 / 23.0, 1.0e-12);

    auto noGas = twoPhaseCell();
    noGas.gasSaturation = 0.0;
    const auto [oil2, gas2] = Container::phaseSolventMassFractions(noGas);
    BOOST_CHECK_CLOSE(oil2, 2.0 / 19.0, 1.0e-12);
    BOOST_CHECK_EQUAL(gas2, 0.0);
}

BOOST_AUTO_TEST_CASE(LimitingCompositions)
{
    // Oil without dissolved gas holds no solvent; dry gas of pure solvent is
    // all solvent.
    auto cell = twoPhaseCell();
    cell.rs = 0.0;
    cell.rv = 0.0;
    cell.yVolume = 1.0;
    const auto [oil, gas] = Container::phaseSolventMassFractions(cell);
    BOOST_CHECK_EQUAL(oil, 0.0);
    BOOST_CHECK_CLOSE(gas, 1.0, 1.0e-12);
}

BOOST_AUTO_TEST_CASE(RequestedFractionsAreWrittenAsOpmExtended)
{
    Container container;
    auto keywords = requests(1, 1);
    container.allocate(2, keywords, /*extendedOutput=*/true, /*log=*/false);
    BOOST_CHECK_EQUAL(keywords.at("SOLVMFO"), 0);
    BOOST_CHECK_EQUAL(keywords.at("SOLVMFG"), 0);
    BOOST_REQUIRE(container.phaseMassFractionsRequested());

    container.assignPhaseMassFractions(0, 0.1, 0.6);
    container.assignPhaseMassFractions(1, 0.2, 0.7);
    Opm::data::Solution sol;
    container.outputRestart(sol);
    BOOST_CHECK(!container.allocated());

    BOOST_REQUIRE(sol.has("SOLVMFO"));
    BOOST_REQUIRE(sol.has("SOLVMFG"));
    for (const auto* name : {"SOLVMFO", "SOLVMFG"}) {
        BOOST_CHECK(sol.at(name).target == Opm::data::TargetType::RESTART_OPM_EXTENDED);
        BOOST_CHECK(sol.at(name).dim == Opm::UnitSystem::measure::identity);
    }

    const auto& oil = sol.data<double>("SOLVMFO");
    BOOST_REQUIRE_EQUAL(oil.size(), 2);
    BOOST_CHECK_EQUAL(oil[0], 0.1);
    BOOST_CHECK_EQUAL(oil[1], 0.2);
    const auto& gas = sol.data<double>("SOLVMFG");
    BOOST_REQUIRE_EQUAL(gas.size(), 2);
    BOOST_CHECK_EQUAL(gas[0], 0.6);
    BOOST_CHECK_EQUAL(gas[1], 0.7);
}

BOOST_AUTO_TEST_CASE(OnlyRequestedFractionIsWritten)
{
    Container container;
    auto keywords = requests(0, 1);
    container.allocate(2, keywords, /*extendedOutput=*/true, /*log=*/false);
    Opm::data::Solution sol;
    container.outputRestart(sol);
    BOOST_CHECK(!sol.has("SOLVMFO"));
    BOOST_CHECK(sol.has("SOLVMFG"));
}

BOOST_AUTO_TEST_CASE(FractionsNeedExtendedOutput)
{
    // Without the extended restart file, or with NORST > 0, the request is
    // consumed but nothing is computed or written.
    Container container;
    auto keywords = requests(1, 1);
    container.allocate(2, keywords, /*extendedOutput=*/false, /*log=*/false);
    BOOST_CHECK_EQUAL(keywords.at("SOLVMFO"), 0);
    BOOST_CHECK_EQUAL(keywords.at("SOLVMFG"), 0);
    BOOST_CHECK(!container.phaseMassFractionsRequested());

    Opm::data::Solution sol;
    container.outputRestart(sol);
    BOOST_CHECK(!sol.has("SOLVMFO"));
    BOOST_CHECK(!sol.has("SOLVMFG"));

    // NORST may change between report steps.
    keywords = requests(1, 1);
    container.allocate(2, keywords, /*extendedOutput=*/true, /*log=*/false);
    BOOST_CHECK(container.phaseMassFractionsRequested());
    keywords = requests(1, 1);
    container.allocate(2, keywords, /*extendedOutput=*/false, /*log=*/false);
    BOOST_CHECK(!container.phaseMassFractionsRequested());
}

BOOST_AUTO_TEST_CASE(UnwritableRequestWarnsOnce)
{
    auto counter = std::make_shared<Opm::CounterLog>();
    Opm::OpmLog::addBackend("COUNTER", counter);

    Container container;
    for (int step = 0; step < 3; ++step) {
        auto keywords = requests(1, 1);
        container.allocate(2, keywords, /*extendedOutput=*/false, /*log=*/true);
    }
    BOOST_CHECK_EQUAL(counter->numMessages(Opm::Log::MessageType::Warning), 1);

    // No warning without a request, when the arrays can be written, or when
    // not logging.
    Container other;
    auto none = requests(0, 0);
    other.allocate(2, none, /*extendedOutput=*/false, /*log=*/true);
    auto writable = requests(1, 1);
    other.allocate(2, writable, /*extendedOutput=*/true, /*log=*/true);
    Container silent;
    auto keywords = requests(1, 1);
    silent.allocate(2, keywords, /*extendedOutput=*/false, /*log=*/false);
    BOOST_CHECK_EQUAL(counter->numMessages(Opm::Log::MessageType::Warning), 1);

    Opm::OpmLog::removeBackend("COUNTER");
}
