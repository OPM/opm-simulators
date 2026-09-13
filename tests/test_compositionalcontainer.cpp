/*
  Copyright 2026 SINTEF Digital

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

#define BOOST_TEST_MODULE CompositionalContainer
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/CompositionalContainer.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>

#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>
#include <opm/output/data/Solution.hpp>

#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <vector>

namespace {
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<double, 3, false>;
using Container = Opm::CompositionalContainer<FluidSystem>;
using RestartOutput = Container::RestartOutput;

// The three-component CO2/methane/decane fluid of the saturation-pressure
// solver's own tests, so their reference values can be reused here.
struct Fixture
{
    Fixture()
    {
        using CompParam = FluidSystem::ComponentParam;
        FluidSystem::init();
        FluidSystem::addComponent(CompParam{"CO2", 44.0, 304.128, 73.773e5, 0.09412, 0.22394});
        FluidSystem::addComponent(CompParam{"C1", 16.04, 190.564, 45.992e5, 0.09863, 0.01142});
        FluidSystem::addComponent(CompParam{"C10", 142.28, 617.7, 21.03e5, 0.60980, 0.4884});
    }
};
} // anonymous namespace

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(SaturationPressureRequiresRestartOutput)
{
    Container container;
    std::map<std::string, int> keywords{{"PSAT", 1}};
    container.allocate(2, keywords, RestartOutput::Disabled);
    BOOST_CHECK(!container.saturationPressureAllocated());
    BOOST_CHECK(!container.saturationPressureRequested());
    BOOST_CHECK_EQUAL(keywords.at("PSAT"), 1);

    container.allocate(2, keywords, RestartOutput::Enabled);
    BOOST_REQUIRE(container.saturationPressureAllocated());
    BOOST_CHECK(container.saturationPressureRequested());
    BOOST_CHECK_EQUAL(keywords.at("PSAT"), 0);
    container.assignSaturationPressure(0, 100.0e5);

    // A summary-only pass must discard the old buffer, even before export.
    keywords["PSAT"] = 1;
    container.allocate(2, keywords, RestartOutput::Disabled);
    BOOST_CHECK(!container.saturationPressureAllocated());
    BOOST_CHECK(!container.saturationPressureRequested());
    BOOST_CHECK_EQUAL(keywords.at("PSAT"), 1);
    Opm::data::Solution substep;
    std::vector<double> oilSaturation;
    container.outputRestart(substep, oilSaturation);
    BOOST_CHECK(!substep.has("PSAT"));

    // Re-enable the request on a later restart snapshot.
    keywords["PSAT"] = 1;
    container.allocate(2, keywords, RestartOutput::Enabled);
    BOOST_REQUIRE(container.saturationPressureAllocated());
    BOOST_CHECK(container.saturationPressureRequested());
    container.assignSaturationPressure(0, 150.0e5);
    container.assignSaturationPressure(1, 200.0e5);
    Opm::data::Solution restart;
    container.outputRestart(restart, oilSaturation);
    BOOST_REQUIRE(restart.has("PSAT"));
    const auto& pressure = restart.data<double>("PSAT");
    BOOST_REQUIRE_EQUAL(pressure.size(), 2);
    BOOST_CHECK_EQUAL(pressure[0], 150.0e5);
    BOOST_CHECK_EQUAL(pressure[1], 200.0e5);
    BOOST_CHECK(restart.at("PSAT").dim == Opm::UnitSystem::measure::pressure);
    BOOST_CHECK(!container.saturationPressureAllocated());
}

BOOST_AUTO_TEST_CASE(DisabledSaturationPressureClearsPreviousRequest)
{
    Container container;
    std::map<std::string, int> keywords{{"PSAT", 1}};
    container.allocate(2, keywords, RestartOutput::Enabled);
    BOOST_REQUIRE(container.saturationPressureAllocated());
    keywords["PSAT"] = 0;
    container.allocate(2, keywords, RestartOutput::Enabled);
    BOOST_CHECK(!container.saturationPressureAllocated());
    BOOST_CHECK(!container.saturationPressureRequested());
}

BOOST_AUTO_TEST_CASE(SaturationPressureRequestIsIndependentOfLocalBufferSize)
{
    Container container;
    std::map<std::string, int> keywords{{"PSAT", 1}};
    container.allocate(0, keywords, RestartOutput::Enabled);

    BOOST_CHECK(container.saturationPressureRequested());
    BOOST_CHECK(!container.saturationPressureAllocated());
    BOOST_CHECK_EQUAL(keywords.at("PSAT"), 0);
}

BOOST_AUTO_TEST_CASE(CellSaturationPressureSelectsThePhaseAndTheBranch)
{
    // 100 degC with Peng-Robinson, as in the solver's tests.  Phase presence
    // comes from the flash liquid fraction L, not from the saturations.
    constexpr double temperature = 373.15;
    constexpr auto eos = Opm::CompositionalConfig::EOSType::PR;
    using CompVec = std::array<double, 3>;

    // Both hydrocarbon phases present: the cell is at its saturation pressure
    // and the solver is not consulted.
    {
        const auto psat = Container::cellSaturationPressure(
            0.3, 75.0e5, CompVec{0.0, 0.5, 0.5}, temperature, eos);
        BOOST_REQUIRE(psat.has_value());
        BOOST_CHECK_CLOSE(*psat, 75.0e5, 1.0e-12);
    }

    // Liquid only (L == 1): the bubble point of the total composition. The
    // reference bubble pressure is 160.5601 bar.
    {
        const auto psat = Container::cellSaturationPressure(
            1.0, 150.0e5, CompVec{0.0, 0.5, 0.5}, temperature, eos);
        BOOST_REQUIRE(psat.has_value());
        BOOST_CHECK_CLOSE(*psat / 1.0e5, 160.56010, 1.0e-3);
    }

    // Vapour only (L == 0): the retrograde dew point of that liquid's
    // equilibrium vapour, which recovers the same pressure.  The vapour is
    // known to single precision, hence the tolerance.
    {
        const auto psat = Container::cellSaturationPressure(
            0.0, 150.0e5, CompVec{0.0, 0.987784, 0.012216}, temperature, eos);
        BOOST_REQUIRE(psat.has_value());
        BOOST_CHECK_CLOSE(*psat / 1.0e5, 160.56010, 1.0e-2);
    }

    // Pure methane far above its critical temperature has no saturation
    // pressure, and the cell pressure must not stand in for one.
    {
        const auto psat = Container::cellSaturationPressure(
            0.0, 150.0e5, CompVec{0.0, 1.0, 0.0}, temperature, eos);
        BOOST_CHECK(!psat.has_value());
    }

    // Rachford-Rice values near the exact single-phase labels still denote
    // two phases, even when round-off puts them slightly outside (0, 1).
    for (const double liquidFraction : {1.0 + 2.0e-16, -1.0e-16, 1.0 - 1.0e-12, 1.0e-12}) {
        const auto psat = Container::cellSaturationPressure(
            liquidFraction, 75.0e5, CompVec{0.0, 0.5, 0.5}, temperature, eos);
        BOOST_REQUIRE_MESSAGE(psat.has_value(), "L = " << liquidFraction);
        BOOST_CHECK_CLOSE(*psat, 75.0e5, 1.0e-12);
    }
}

BOOST_AUTO_TEST_CASE(RoundOffInTheGasSaturationDoesNotHideTheBubblePoint)
{
    // Sg = max(1 - So - Sw, 0) leaves a positive round-off residual at Sw = 0.3.
    // Classifying this liquid-only cell by saturation would return the cell
    // pressure instead of the bubble point.
    constexpr double temperature = 373.15;
    constexpr auto eos = Opm::CompositionalConfig::EOSType::PR;
    const std::array<double, 3> liquid{0.0, 0.5, 0.5};

    constexpr double waterSaturation = 0.3;
    const double oilSaturation = 1.0 - waterSaturation;
    const double gasSaturation = std::max(1.0 - oilSaturation - waterSaturation, 0.0);
    BOOST_REQUIRE_GT(gasSaturation, 0.0);
    BOOST_REQUIRE_LT(gasSaturation, 1.0e-15);

    // The flash's liquid-only label must select the bubble point.
    const auto psat = Container::cellSaturationPressure(
        1.0, 200.0e5, liquid, temperature, eos);
    BOOST_REQUIRE(psat.has_value());
    BOOST_CHECK_CLOSE(*psat / 1.0e5, 160.56010, 1.0e-3);
}
