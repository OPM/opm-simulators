/*
  Copyright 2026, SINTEF Digital

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

#include "config.h"

#define BOOST_TEST_MODULE CompositionalEquil
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/equil/InitStateEquilComp.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <opm/material/constraintsolvers/SaturationPressure.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace {

using Scalar = double;
using FluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, 3, false>;
using InitialStateComputer = Opm::EQUIL::Comp::InitialStateComputer<FluidSystem>;
using SatP = Opm::SaturationPressure<Scalar, FluidSystem>;
using CompVec = std::array<Scalar, 3>;

constexpr Scalar barsa = 1.0e5;
constexpr Scalar gravity = 9.80665;

// A 1x1x20 vertical column from 2000 m to 2100 m in 5 m cells, filled with a
// CO2/methane/decane mixture that grades from methane-rich at the top to
// decane-rich at the bottom.  This is the geometry of the deck the
// implementation was verified against, so the reference-simulator numbers
// quoted below apply to it directly.
std::string deckString(const std::string& equil,
                       const std::string& runspecExtra = "EQLDIMS\n/\n",
                       const std::string& regions = "",
                       const std::string& zmfvd =
                           "ZMFVD\n"
                           " 2000   0 0.7 0.3\n"
                           " 2100   0 0.3 0.7  /\n")
{
    return
        "RUNSPEC\n"
        "METRIC\n"
        "TABDIMS\n/\n"
        "OIL\nGAS\n"
        "DIMENS\n1 1 20 /\n"
        "COMPS\n3 /\n"
        "START\n  1 'JAN' 2016  /\n"
        + runspecExtra +
        "GRID\n"
        "DX\n20*5 /\n"
        "DY\n20*100 /\n"
        "DZ\n20*5 /\n"
        "TOPS\n1*2000.0 /\n"
        "PORO\n20*0.3 /\n"
        "PERMX\n20*2000 /\n"
        "PERMY\n20*2000 /\n"
        "PERMZ\n20*2000 /\n"
        "PROPS\n"
        "CNAMES\nCO2\nMETHANE\nDECANE\n/\n"
        "ROCK\n68.9476 0 /\n"
        + zmfvd +
        "RTEMP\n100\n/\n"
        "EOS\nPR /\n"
        "BIC\n0\n0\n0\n/\n"
        "ACF\n0.22394\n0.01142\n0.4884\n/\n"
        "PCRIT\n73.773\n45.992\n21.03\n/\n"
        "TCRIT\n304.128\n190.564\n617.7\n/\n"
        "MW\n44.00\n16.04\n142.28\n/\n"
        "VCRIT\n0.09412\n0.09863\n0.60980\n/\n"
        "STCOND\n15.0 1.0 /\n"
        + regions +
        "SOLUTION\n"
        + equil +
        "END\n";
}

struct EquilFixture
{
    explicit EquilFixture(const std::string& deck_string)
        : deck(Opm::Parser{}.parseString(deck_string))
        , eclState(deck)
        , schedule(deck, eclState, std::make_shared<Opm::Python>())
    {
        FluidSystem::initFromState(eclState, schedule);
        for (std::size_t c = 0; c < depths.size(); ++c) {
            depths[c] = 2002.5 + 5.0 * static_cast<Scalar>(c);
        }
    }

    InitialStateComputer compute(const std::vector<int>& eqlnum) const
    {
        return InitialStateComputer(eclState,
                                    eclState.compositionalConfig().eosType(0),
                                    {depths.begin(), depths.end()},
                                    eqlnum,
                                    Opm::Parallel::Communication{},
                                    gravity,
                                    /*numSamplePoints=*/100);
    }

    Opm::Deck deck;
    Opm::EclipseState eclState;
    Opm::Schedule schedule;
    std::array<Scalar, 20> depths{};
};

// The mixture composition the ZMFVD table prescribes at a depth.
CompVec tableComposition(const Scalar depth)
{
    const Scalar t = (depth - 2000.0) / 100.0;
    return {0.0, 0.7 - 0.4 * t, 0.3 + 0.4 * t};
}

// The density a pressure difference between two neighbouring cells implies
// through the hydrostatic relation dp = rho g dz.
Scalar impliedDensity(const Scalar pAbove, const Scalar pBelow)
{
    return (pBelow - pAbove) / (gravity * 5.0);
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(ContinuousLiquidColumn)
{
    // EQUIL item 10 is 1 and the gas-oil contact sits at the top of the
    // column, i.e. no free gas: a continuous liquid whose composition follows
    // ZMFVD, integrated hydrostatically from the datum.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2000 0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    for (std::size_t c = 0; c < states.size(); ++c) {
        const auto& fs = states[c];
        BOOST_CHECK_CLOSE(fs.saturation(FluidSystem::oilPhaseIdx), 1.0, 1e-10);
        BOOST_CHECK_SMALL(fs.saturation(FluidSystem::gasPhaseIdx), 1e-10);
        BOOST_CHECK_CLOSE(fs.temperature(FluidSystem::oilPhaseIdx), 373.15, 1e-10);

        // The total composition is the ZMFVD interpolant at the cell centre.
        const CompVec z = tableComposition(fix.depths[c]);
        for (int comp = 0; comp < 3; ++comp) {
            BOOST_CHECK_SMALL(std::abs(Opm::getValue(fs.moleFraction(comp)) - z[comp]),
                              1e-10);
        }
    }

    // The pressure increases with depth at a liquid-like gradient.
    for (std::size_t c = 0; c + 1 < states.size(); ++c) {
        const Scalar rho =
            impliedDensity(Opm::getValue(states[c].pressure(FluidSystem::oilPhaseIdx)),
                           Opm::getValue(states[c + 1].pressure(FluidSystem::oilPhaseIdx)));
        BOOST_CHECK_GT(rho, 400.0);
        BOOST_CHECK_LT(rho, 800.0);
    }

    // The reference simulator initializes this column to 149.666 barsa in the
    // top cell and 154.694 barsa in the bottom one.
    BOOST_CHECK_SMALL(std::abs(Opm::getValue(states.front().pressure(FluidSystem::oilPhaseIdx))
                               - 149.666 * barsa), 0.05 * barsa);
    BOOST_CHECK_SMALL(std::abs(Opm::getValue(states.back().pressure(FluidSystem::oilPhaseIdx))
                               - 154.694 * barsa), 0.05 * barsa);
}

BOOST_AUTO_TEST_CASE(GasCapAboveContact)
{
    // EQUIL item 10 is 3: ZMFVD is the liquid composition and the gas-oil
    // contact at 2050 m lies inside the column.  The contact pressure is the
    // saturation pressure of the contact liquid (item 11 is defaulted and the
    // datum pressure is more than one atmosphere away), and the gas above the
    // contact holds the equilibrium vapour of that liquid.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2050 0 3* 3 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();

    // Independently computed saturation point of the contact liquid; the
    // reference simulator puts it at 160.5601 barsa.
    const CompVec liquid = tableComposition(2050.0);
    Scalar psat = 0.0;
    CompVec vapor{};
    BOOST_REQUIRE(SatP::bubblePressure(liquid, 373.15,
                                       fix.eclState.compositionalConfig().eosType(0),
                                       psat, vapor));
    BOOST_CHECK_SMALL(std::abs(psat - 160.5601 * barsa), 0.05 * barsa);

    for (std::size_t c = 0; c < states.size(); ++c) {
        const auto& fs = states[c];
        const bool inGas = fix.depths[c] < 2050.0;
        BOOST_CHECK_CLOSE(fs.saturation(inGas ? FluidSystem::gasPhaseIdx
                                              : FluidSystem::oilPhaseIdx), 1.0, 1e-10);

        // Every gas cell carries the contact vapour; the liquid cells follow
        // the table.
        const CompVec z = inGas ? vapor : tableComposition(fix.depths[c]);
        for (int comp = 0; comp < 3; ++comp) {
            BOOST_CHECK_SMALL(std::abs(Opm::getValue(fs.moleFraction(comp)) - z[comp]),
                              1e-6);
        }
    }
    // The vapour is far richer in the light component than the table value at
    // the contact depth.
    BOOST_CHECK_GT(vapor[1], 0.9);

    // The pressure passes through the saturation pressure at the contact,
    // with a gas-like gradient above it and a liquid-like one below.
    const Scalar pAbove = Opm::getValue(states[9].pressure(FluidSystem::gasPhaseIdx));
    const Scalar pBelow = Opm::getValue(states[10].pressure(FluidSystem::oilPhaseIdx));
    BOOST_CHECK_LT(pAbove, psat);
    BOOST_CHECK_GT(pBelow, psat);

    for (std::size_t c = 0; c + 1 < states.size(); ++c) {
        const Scalar rho =
            impliedDensity(Opm::getValue(states[c].pressure(FluidSystem::oilPhaseIdx)),
                           Opm::getValue(states[c + 1].pressure(FluidSystem::oilPhaseIdx)));
        if (c + 1 <= 9) {         // both cells in the gas cap
            BOOST_CHECK_GT(rho, 30.0);
            BOOST_CHECK_LT(rho, 300.0);
        }
        else if (c >= 10) {       // both cells in the liquid leg
            BOOST_CHECK_GT(rho, 350.0);
            BOOST_CHECK_LT(rho, 900.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(GasCapKeepingDatumPressure)
{
    // As GasCapAboveContact, but EQUIL item 11 is 1: the datum pressure is
    // kept at the contact instead of being replaced by the saturation
    // pressure of the contact liquid.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2050 0 3* 3 1 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();

    const Scalar pContact = 150.0 * barsa;
    BOOST_CHECK_LT(Opm::getValue(states[9].pressure(FluidSystem::gasPhaseIdx)), pContact);
    BOOST_CHECK_GT(Opm::getValue(states[10].pressure(FluidSystem::oilPhaseIdx)), pContact);

    // Still a gas cap over a liquid leg.
    BOOST_CHECK_CLOSE(states[0].saturation(FluidSystem::gasPhaseIdx), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(states[19].saturation(FluidSystem::oilPhaseIdx), 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(TwoIndependentRegions)
{
    // Two equilibration regions splitting the column in half, each with its
    // own datum pressure.  The regions must be integrated independently: the
    // 50 bar difference between the records shows up as a pressure jump at
    // the region boundary that a single hydrostatic column could never have.
    const EquilFixture fix(deckString(
        "EQUIL\n"
        " 2010 150 2300 0 2000 0 /\n"
        " 2060 200 2300 0 2050 0 /\n",
        "EQLDIMS\n2 /\n"
        "REGDIMS\n2 1 0 0 /\n",
        "REGIONS\n"
        "EQLNUM\n10*1 10*2 /\n",
        "ZMFVD\n"
        " 2000   0 0.7 0.3\n"
        " 2100   0 0.3 0.7  /\n"
        " 2000   0 0.7 0.3\n"
        " 2100   0 0.3 0.7  /\n"));

    std::vector<int> eqlnum(20, 0);
    std::fill(eqlnum.begin() + 10, eqlnum.end(), 1);
    const auto states = fix.compute(eqlnum).fluidStates();

    // Both regions hold single-phase liquid rising in pressure with depth.
    for (std::size_t c = 0; c + 1 < states.size(); ++c) {
        if (c == 9) {
            continue;
        }
        BOOST_CHECK_LT(Opm::getValue(states[c].pressure(FluidSystem::oilPhaseIdx)),
                       Opm::getValue(states[c + 1].pressure(FluidSystem::oilPhaseIdx)));
    }
    const Scalar jump = Opm::getValue(states[10].pressure(FluidSystem::oilPhaseIdx))
                      - Opm::getValue(states[9].pressure(FluidSystem::oilPhaseIdx));
    BOOST_CHECK_GT(jump, 40.0 * barsa);
}

BOOST_AUTO_TEST_CASE(MissingZmfvdIsAnError)
{
    // The composition versus depth is the one piece of input the
    // equilibration cannot invent.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2000 0 /\n",
                                      "EQLDIMS\n/\n", "", ""));
    BOOST_CHECK_THROW(fix.compute(std::vector<int>(20, 0)),
                      std::runtime_error);
}

namespace {

struct MpiFixture
{
    MpiFixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};

} // Anonymous namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);
