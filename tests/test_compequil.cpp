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
// The water zone is only reachable when the water phase is compiled in.
using WaterFluidSystem = Opm::GenericOilGasWaterFluidSystem<Scalar, 3, true>;
using InitialStateComputer = Opm::EQUIL::Comp::InitialStateComputer<FluidSystem>;
using SatP = Opm::SaturationPressure<Scalar, FluidSystem>;
using CompVec = std::array<Scalar, 3>;

constexpr Scalar barsa = 1.0e5;
constexpr Scalar gravity = 9.80665;

// A 1x1x20 vertical column from 2000 m to 2100 m in 5 m cells, filled with a
// CO2/methane/decane mixture that grades from methane-rich at the top to
// decane-rich at the bottom. This is the geometry used by the numeric
// expectations below.
std::string deckString(const std::string& equil,
                       const std::string& runspecExtra = "EQLDIMS\n/\n",
                       const std::string& regions = "",
                       const std::string& zmfvd =
                           "ZMFVD\n"
                           " 2000   0 0.7 0.3\n"
                           " 2100   0 0.3 0.7  /\n",
                       const std::string& rtemp = "RTEMP\n100\n/\n",
                       const std::string& phases = "OIL\nGAS\n",
                       const std::string& propsExtra = "")
{
    return
        "RUNSPEC\n"
        "METRIC\n"
        "TABDIMS\n/\n"
        + phases +
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
        + zmfvd
        + rtemp +
        "EOS\nPR /\n"
        "BIC\n0\n0\n0\n/\n"
        "ACF\n0.22394\n0.01142\n0.4884\n/\n"
        "PCRIT\n73.773\n45.992\n21.03\n/\n"
        "TCRIT\n304.128\n190.564\n617.7\n/\n"
        "MW\n44.00\n16.04\n142.28\n/\n"
        "VCRIT\n0.09412\n0.09863\n0.60980\n/\n"
        "STCOND\n15.0 1.0 /\n"
        + propsExtra
        + regions +
        "SOLUTION\n"
        + equil +
        "END\n";
}

template <class FS>
struct BasicEquilFixture
{
    using Computer = Opm::EQUIL::Comp::InitialStateComputer<FS>;

    explicit BasicEquilFixture(const std::string& deck_string)
        : deck(Opm::Parser{}.parseString(deck_string))
        , eclState(deck)
        , schedule(deck, eclState, std::make_shared<Opm::Python>())
    {
        FS::initFromState(eclState, schedule);
        for (std::size_t c = 0; c < depths.size(); ++c) {
            depths[c] = 2002.5 + 5.0 * static_cast<Scalar>(c);
        }
    }

    /// \param connateWater, maxWater  Per-cell water saturation endpoints, as
    ///        the simulator reads them from the scaled saturation functions.
    Computer compute(const std::vector<int>& eqlnum,
                     const std::vector<Scalar>& connateWater = {},
                     const std::vector<Scalar>& maxWater = {}) const
    {
        return Computer(eclState,
                        eclState.compositionalConfig().eosType(0),
                        {depths.begin(), depths.end()},
                        eqlnum,
                        Opm::Parallel::Communication{},
                        gravity,
                        /*numSamplePoints=*/100,
                        connateWater,
                        maxWater);
    }

    Opm::Deck deck;
    Opm::EclipseState eclState;
    Opm::Schedule schedule;
    std::array<Scalar, 20> depths{};
};

using EquilFixture = BasicEquilFixture<FluidSystem>;
using WaterEquilFixture = BasicEquilFixture<WaterFluidSystem>;

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

// The connate water saturation of the saturation function below, i.e. the
// water the hydrocarbon column retains above the water-oil contact.
constexpr Scalar connateSw = 0.01;

// The three-phase counterpart of deckString().
std::string waterDeckString(const std::string& equil,
                            const std::string& composition =
                                "ZMFVD\n"
                                " 2000   0 0.7 0.3\n"
                                " 2100   0 0.3 0.7  /\n")
{
    return deckString(equil, "EQLDIMS\n/\n", "", composition,
                      "RTEMP\n100\n/\n",
                      "OIL\nGAS\nWATER\n",
                      "SWFN\n 0.01 0.0 0.0\n 1.00 1.0 0.0 /\n"
                      "SGFN\n 0.00 0.0 0.0\n 0.99 1.0 0.0 /\n"
                      "SOF3\n 0.00 0.0 0.0\n 0.99 1.0 1.0 /\n");
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Type1LiquidRootPressureIntegration)
{
    // Type 1 takes ZMFVD as the total composition. Since the datum lies below
    // the gas-oil contact, the initializer uses the liquid EOS root to integrate
    // pressure from the datum and marks oil as the nominal phase. The downstream
    // flash may still split the mixture into oil and gas.
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

    // The expected pressures are 149.666 barsa in the top cell and
    // 154.694 barsa in the bottom one.
    BOOST_CHECK_SMALL(std::abs(Opm::getValue(states.front().pressure(FluidSystem::oilPhaseIdx))
                               - 149.666 * barsa), 0.05 * barsa);
    BOOST_CHECK_SMALL(std::abs(Opm::getValue(states.back().pressure(FluidSystem::oilPhaseIdx))
                               - 154.694 * barsa), 0.05 * barsa);
}

BOOST_AUTO_TEST_CASE(Type1VapourRootPressureIntegration)
{
    // The contact is below the whole column, so type 1 must use the vapour
    // root. At 10 bar and 100 C this mixture has distinct gas and liquid roots
    // with densities around 40 and 494 kg/m^3, respectively. The pressure
    // gradient therefore detects an incorrect choice of EOS root.
    const EquilFixture fix(deckString("EQUIL\n 2012.5 10 2300 0 2200 0 /\n",
                                      "EQLDIMS\n/\n", "",
                                      "ZMFVD\n 2000 0 0.5 0.5 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    // The datum coincides with the third cell centre.
    BOOST_CHECK_CLOSE(states[2].pressure(FluidSystem::gasPhaseIdx), 10.0 * barsa, 1e-10);

    const CompVec z{0.0, 0.5, 0.5};
    for (const auto& fs : states) {
        BOOST_CHECK_CLOSE(fs.temperature(FluidSystem::gasPhaseIdx), 373.15, 1e-10);
        for (int comp = 0; comp < 3; ++comp) {
            BOOST_CHECK_SMALL(fs.moleFraction(comp) - z[comp], 1e-10);
        }
    }

    for (std::size_t c = 0; c + 1 < states.size(); ++c) {
        const Scalar rho = impliedDensity(states[c].pressure(FluidSystem::gasPhaseIdx),
                                         states[c + 1].pressure(FluidSystem::gasPhaseIdx));
        BOOST_CHECK_GT(rho, 30.0);
        BOOST_CHECK_LT(rho, 60.0);
    }
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

    // Independently compute the saturation point of the contact liquid and
    // check the expected value of 160.5601 barsa.
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
    // As GasCapAboveContact, but EQUIL item 11 is 1: the reference depth is
    // reset from 2010 m to the contact while the numeric 150 bar input pressure
    // is retained there. This intentionally need not be an equilibrium
    // saturation pressure.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2050 0 3* 3 1 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();

    const Scalar pContact = 150.0 * barsa;
    BOOST_CHECK_LT(Opm::getValue(states[9].pressure(FluidSystem::gasPhaseIdx)), pContact);
    BOOST_CHECK_GT(Opm::getValue(states[10].pressure(FluidSystem::oilPhaseIdx)), pContact);

    // Still a gas cap over a liquid leg.
    BOOST_CHECK_CLOSE(states[0].saturation(FluidSystem::gasPhaseIdx), 1.0, 1e-10);
    BOOST_CHECK_CLOSE(states[19].saturation(FluidSystem::oilPhaseIdx), 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(CompositionTableIsInheritedAcrossRegions)
{
    // A region without its own composition record reuses the nearest preceding
    // record, following TableContainer semantics.
    const EquilFixture fix(deckString(
        "EQUIL\n 2010 150 2300 0 2000 0 /\n 2010 150 2300 0 2000 0 /\n",
        "EQLDIMS\n2 /\n", "REGIONS\nEQLNUM\n10*1 10*2 /\n",
        "ZMFVD\n"
        " 2000   0 0.7 0.3\n"
        " 2100   0 0.3 0.7  /\n"
        "/\n"));

    std::vector<int> eqlnum(20, 0);
    std::fill(eqlnum.begin() + 10, eqlnum.end(), 1);
    const auto states = fix.compute(eqlnum).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    // Both regions equilibrate off the same table, so the column is the one a
    // single region would have produced.
    for (std::size_t c = 0; c < states.size(); ++c) {
        BOOST_TEST_CONTEXT("Cell " << c) {
            const CompVec z = tableComposition(fix.depths[c]);
            for (int comp = 0; comp < 3; ++comp) {
                BOOST_CHECK_SMALL(states[c].moleFraction(comp) - z[comp], Scalar{1.0e-10});
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(BothCompositionKeywordsUseEachRegionsOwnRecord)
{
    // A deck may give ZMFVD for one region and COMPVD for another as long as
    // every region states which one it reads.
    const EquilFixture fix(deckString(
        "EQUIL\n 2010 150 2300 0 2000 0 /\n 2010 150 2300 0 2000 0 /\n",
        "EQLDIMS\n2 /\n", "REGIONS\nEQLNUM\n10*1 10*2 /\n",
        "ZMFVD\n 2000 0 0.7 0.3\n 2100 0 0.3 0.7 /\n/\n"
        "COMPVD\n/\n"
        " 2000 0 0.5 0.5 1 150.0\n"
        " 2100 0 0.5 0.5 1 150.0 /\n"));

    std::vector<int> eqlnum(20, 0);
    std::fill(eqlnum.begin() + 10, eqlnum.end(), 1);
    const auto states = fix.compute(eqlnum).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    for (std::size_t c = 0; c < states.size(); ++c) {
        BOOST_TEST_CONTEXT("Cell " << c) {
            const CompVec z = (c < 10) ? tableComposition(fix.depths[c])
                                       : CompVec{0.0, 0.5, 0.5};
            for (int comp = 0; comp < 3; ++comp) {
                BOOST_CHECK_SMALL(states[c].moleFraction(comp) - z[comp], Scalar{1.0e-10});
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(BothCompositionKeywordsNeedARecordForEveryRegion)
{
    // Neither keyword can be inherited once a deck uses both: there is no rule
    // saying which one a region that states nothing would follow.
    const EquilFixture fix(deckString(
        "EQUIL\n 2010 150 2300 0 2000 0 /\n"
        " 2010 150 2300 0 2000 0 /\n"
        " 2010 150 2300 0 2000 0 /\n",
        "EQLDIMS\n3 /\n", "REGIONS\nEQLNUM\n7*1 7*2 6*3 /\n",
        "ZMFVD\n 2000 0 0.7 0.3\n 2100 0 0.3 0.7 /\n/\n/\n"
        "COMPVD\n/\n"
        " 2000 0 0.5 0.5 1 150.0\n"
        " 2100 0 0.5 0.5 1 150.0 /\n/\n"));

    std::vector<int> eqlnum(20, 0);
    std::fill(eqlnum.begin() + 7, eqlnum.begin() + 14, 1);
    std::fill(eqlnum.begin() + 14, eqlnum.end(), 2);
    BOOST_CHECK_EXCEPTION(fix.compute(eqlnum), std::runtime_error,
                          [](const std::runtime_error& error) {
                              const std::string message = error.what();
                              return message.find("Region 3") != std::string::npos
                                  && message.find("record of its own") != std::string::npos;
                          });
}

BOOST_AUTO_TEST_CASE(CompvdTwoZoneContactOutsideTheCells)
{
    // The contact can lie below every cell with the datum beyond it. The
    // column carrying the datum is then integrated to the contact before the
    // other one starts, so the interval between the contact and the cells
    // keeps the density of its own phase.
    const auto bottomCellPressure = [](const std::string& goc) {
        const EquilFixture fix(deckString(
            "EQUIL\n 2300 300 2400 0 " + goc + " 0 /\n", "EQLDIMS\n/\n", "",
            "COMPVD\n"
            " 2000   0 0.95 0.05  0  150.0\n"
            " 2100   0 0.95 0.05  0  150.0\n"
            " 2250   0 0.60 0.40  1  150.0\n"
            " 2400   0 0.40 0.60  1  150.0 /\n"));
        const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
        // Every cell lies above either contact, so the whole column is gas.
        BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});
        return states.back().pressure(FluidSystem::gasPhaseIdx);
    };

    // Moving the contact down turns part of the path from the datum from
    // liquid into gas, which is lighter, so the cells gain pressure. Pinning
    // the contact onto the cells would make the two runs identical instead.
    const Scalar shallowContact = bottomCellPressure("2100");
    const Scalar deepContact = bottomCellPressure("2250");
    BOOST_CHECK_GT(deepContact - shallowContact, 3.0 * barsa);
}

BOOST_AUTO_TEST_CASE(NonzeroContactCapillaryPressureIsAnError)
{
    for (const int initType : std::array{1, 3}) {
        for (const int capillaryPressure : std::array{-10, 10}) {
            BOOST_TEST_CONTEXT("Type " << initType << ", PC_GOC " << capillaryPressure) {
                const EquilFixture fix(deckString(
                    "EQUIL\n 2010 150 2300 0 2050 " + std::to_string(capillaryPressure)
                    + " 3* " + std::to_string(initType) + " /\n"));
                BOOST_CHECK_EXCEPTION(fix.compute(std::vector<int>(20, 0)),
                                      std::runtime_error,
                                      [](const std::runtime_error& error) {
                                          const std::string message = error.what();
                                          return message.find("EQUIL item 6") != std::string::npos
                                              && message.find("region 1") != std::string::npos;
                                      });
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(CompositionalEquilAccuracy)
{
    for (const int initType : std::array{1, 3}) {
        for (const int accuracy : std::array{-20, -1, 0, 1, 20}) {
            BOOST_TEST_CONTEXT("Type " << initType << ", accuracy " << accuracy) {
                const EquilFixture fix(deckString(
                    "EQUIL\n 2010 150 2300 0 2050 0 2* " + std::to_string(accuracy)
                    + " " + std::to_string(initType) + " /\n"));
                if (accuracy == 0) {
                    BOOST_CHECK_NO_THROW(fix.compute(std::vector<int>(20, 0)));
                }
                else {
                    BOOST_CHECK_EXCEPTION(fix.compute(std::vector<int>(20, 0)),
                                          std::runtime_error,
                                          [](const std::runtime_error& error) {
                                              const std::string message = error.what();
                                              return message.find("EQUIL item 9") != std::string::npos
                                                  && message.find("region 1") != std::string::npos;
                                          });
                }
            }
        }
    }
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

BOOST_AUTO_TEST_CASE(InvalidEqlnumFailsOnAllRanks)
{
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2000 0 /\n"));
    std::vector<int> eqlnum(20, 0);
    const Opm::Parallel::Communication comm;
    if (comm.rank() == 0) {
        eqlnum.back() = 1;
    }

    BOOST_CHECK_THROW(fix.compute(eqlnum), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(MismatchedEqlnumSizeFailsOnAllRanks)
{
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2000 0 /\n"));
    const Opm::Parallel::Communication comm;
    for (const std::size_t size : std::array<std::size_t, 3>{0, 19, 21}) {
        BOOST_TEST_CONTEXT("EQLNUM size " << size) {
            std::vector<int> eqlnum(fix.depths.size(), 0);
            // Only rank 0 has inconsistent input, but every rank must throw
            // before entering the region setup collectives or indexing cells.
            if (comm.rank() == 0) {
                eqlnum.resize(size);
            }
            BOOST_CHECK_EXCEPTION(fix.compute(eqlnum), std::runtime_error,
                                  [&](const std::runtime_error& error) {
                                      const std::string expected = "EQLNUM contains "
                                          + std::to_string(size) + " entries for 20 cell depths";
                                      return comm.rank() != 0
                                          || std::string(error.what()).find(expected)
                                              != std::string::npos;
                                  });
        }
    }
}

BOOST_AUTO_TEST_CASE(ConstantTemperatureFromRtempvd)
{
    // A depth-independent reservoir temperature is a single-row RTEMPVD
    // table, which is as valid a way to state it as RTEMP. The column should
    // reproduce the RTEMP case in Type1LiquidRootPressureIntegration.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2000 0 /\n",
                                      "EQLDIMS\n/\n", "",
                                      "ZMFVD\n"
                                      " 2000   0 0.7 0.3\n"
                                      " 2100   0 0.3 0.7  /\n",
                                      "RTEMPVD\n 2000 100 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    for (const auto& fs : states) {
        BOOST_CHECK_CLOSE(fs.temperature(FluidSystem::oilPhaseIdx), 373.15, 1e-10);
    }
    BOOST_CHECK_SMALL(std::abs(Opm::getValue(states.front().pressure(FluidSystem::oilPhaseIdx))
                               - 149.666 * barsa), 0.05 * barsa);
    BOOST_CHECK_SMALL(std::abs(Opm::getValue(states.back().pressure(FluidSystem::oilPhaseIdx))
                               - 154.694 * barsa), 0.05 * barsa);
}

BOOST_AUTO_TEST_CASE(GradedTemperatureFromRtempvd)
{
    // A two-row RTEMPVD table is interpolated at the cell centre, so the
    // column carries the geothermal gradient the table states.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2000 0 /\n",
                                      "EQLDIMS\n/\n", "",
                                      "ZMFVD\n"
                                      " 2000   0 0.7 0.3\n"
                                      " 2100   0 0.3 0.7  /\n",
                                      "RTEMPVD\n 2000 100\n 2100 120 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    for (std::size_t c = 0; c < states.size(); ++c) {
        const Scalar expected = 373.15 + 0.2 * (fix.depths[c] - 2000.0);
        BOOST_CHECK_CLOSE(states[c].temperature(FluidSystem::oilPhaseIdx), expected, 1e-10);
    }
}

BOOST_AUTO_TEST_CASE(DepthTablesUseConstantEndpoints)
{
    // Extending the tables with constant endpoint values must leave both the
    // cell states and the hydrostatic pressure integration unchanged. The
    // datum and the outer cells lie outside the original table ranges.
    for (const auto* keyword : std::array{"RTEMPVD", "TEMPVD"}) {
        BOOST_TEST_CONTEXT(keyword) {
            const EquilFixture narrow(deckString(
                "EQUIL\n 2010 150 2300 0 2000 0 /\n", "EQLDIMS\n/\n", "",
                "ZMFVD\n 2040 0 0.7 0.3\n 2060 0 0.3 0.7 /\n",
                std::string(keyword) + "\n 2040 100\n 2060 120 /\n"));
            const auto states = narrow.compute(std::vector<int>(20, 0)).fluidStates();

            const EquilFixture padded(deckString(
                "EQUIL\n 2010 150 2300 0 2000 0 /\n", "EQLDIMS\n/\n", "",
                "ZMFVD\n 2000 0 0.7 0.3\n 2040 0 0.7 0.3\n"
                " 2060 0 0.3 0.7\n 2100 0 0.3 0.7 /\n",
                std::string(keyword) + "\n 2000 100\n 2040 100\n"
                                       " 2060 120\n 2100 120 /\n"));
            const auto expected = padded.compute(std::vector<int>(20, 0)).fluidStates();

            for (std::size_t c = 0; c < states.size(); ++c) {
                const Scalar t = std::clamp((narrow.depths[c] - 2040.0) / 20.0, 0.0, 1.0);
                const CompVec z{0.0, 0.7 - 0.4 * t, 0.3 + 0.4 * t};
                for (int comp = 0; comp < 3; ++comp) {
                    BOOST_CHECK_SMALL(states[c].moleFraction(comp) - z[comp], 1e-10);
                }
                BOOST_CHECK_CLOSE(states[c].temperature(FluidSystem::oilPhaseIdx),
                                  373.15 + 20.0 * t, 1e-10);
                BOOST_CHECK_CLOSE(states[c].pressure(FluidSystem::oilPhaseIdx),
                                  expected[c].pressure(FluidSystem::oilPhaseIdx), 1e-10);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(GasOilContactOutsideDepthTables)
{
    // The bubble point and equilibrium vapour must use the endpoint liquid
    // composition and temperature when the contact is outside the tables.
    for (const Scalar contact : std::array{2025.0, 2075.0}) {
        BOOST_TEST_CONTEXT("Contact depth " << contact) {
            EquilFixture fix(deckString(
                "EQUIL\n 2010 150 2300 0 " + std::to_string(contact) + " 0 3* 3 /\n",
                "EQLDIMS\n/\n", "",
                "ZMFVD\n 2040 0 0.7 0.3\n 2060 0 0.3 0.7 /\n",
                "RTEMPVD\n 2040 100\n 2060 120 /\n"));
            const bool aboveTable = contact < 2040.0;
            const std::size_t contactCell = aboveTable ? 4 : 15;
            fix.depths[contactCell] = contact;
            const CompVec liquid = aboveTable ? CompVec{0.0, 0.7, 0.3}
                                             : CompVec{0.0, 0.3, 0.7};
            const Scalar temperature = aboveTable ? 373.15 : 393.15;
            Scalar psat{};
            CompVec vapor{};
            BOOST_REQUIRE(SatP::bubblePressure(liquid, temperature,
                                               fix.eclState.compositionalConfig().eosType(0),
                                               psat, vapor));

            const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
            BOOST_CHECK_CLOSE(states[contactCell].pressure(FluidSystem::oilPhaseIdx),
                              psat, 1e-8);
            BOOST_CHECK_CLOSE(states[contactCell].temperature(FluidSystem::oilPhaseIdx),
                              temperature, 1e-10);
            for (int comp = 0; comp < 3; ++comp) {
                BOOST_CHECK_SMALL(states[contactCell].moleFraction(comp) - liquid[comp], 1e-10);
                BOOST_CHECK_SMALL(states.front().moleFraction(comp) - vapor[comp], 1e-8);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(MissingZmfvdIsAnError)
{
    // Compositional equilibration requires composition-versus-depth input.
    const EquilFixture fix(deckString("EQUIL\n 2010 150 2300 0 2000 0 /\n",
                                      "EQLDIMS\n/\n", "", ""));
    BOOST_CHECK_THROW(fix.compute(std::vector<int>(20, 0)),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(CompvdSinglePhaseMatchesZmfvd)
{
    // The COMPVD liquid flag selects the same EOS root that this ZMFVD case
    // obtains from the datum and contact depths. Its saturation-pressure column
    // is unused, so the initialized states must agree exactly.
    const std::string equil = "EQUIL\n 2010 150 2300 0 2000 0 /\n";
    const EquilFixture zmfvd(deckString(equil));
    const EquilFixture compvd(deckString(equil, "EQLDIMS\n/\n", "",
                                         "COMPVD\n"
                                         " 2000   0 0.7 0.3  1  150.0\n"
                                         " 2100   0 0.3 0.7  1  150.0 /\n"));

    const auto expected = zmfvd.compute(std::vector<int>(20, 0)).fluidStates();
    const auto states = compvd.compute(std::vector<int>(20, 0)).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), expected.size());

    for (std::size_t c = 0; c < states.size(); ++c) {
        BOOST_TEST_CONTEXT("Cell " << c) {
            BOOST_CHECK_CLOSE(states[c].pressure(FluidSystem::oilPhaseIdx),
                              expected[c].pressure(FluidSystem::oilPhaseIdx), 1e-10);
            BOOST_CHECK_CLOSE(states[c].temperature(FluidSystem::oilPhaseIdx),
                              expected[c].temperature(FluidSystem::oilPhaseIdx), 1e-10);
            for (int comp = 0; comp < 3; ++comp) {
                BOOST_CHECK_SMALL(states[c].moleFraction(comp)
                                  - expected[c].moleFraction(comp), 1e-12);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(CompvdStatedPhaseSelectsEosRoot)
{
    // The contact lies above the whole column, which on its own selects the
    // liquid root. COMPVD names the phase its composition belongs to, and that
    // naming decides instead. At 10 bar and 100 C the two roots of this mixture
    // differ by an order of magnitude in density, so the hydrostatic gradient
    // identifies which one was used.
    const auto columnDensity = [](const std::string& phaseFlag) {
        const EquilFixture fix(deckString("EQUIL\n 2012.5 10 2300 0 2000 0 /\n",
                                          "EQLDIMS\n/\n", "",
                                          "COMPVD\n 2000 0 0.5 0.5 " + phaseFlag + " 10.0 /\n"));
        const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
        return impliedDensity(states[9].pressure(FluidSystem::oilPhaseIdx),
                              states[10].pressure(FluidSystem::oilPhaseIdx));
    };

    const Scalar vapour = columnDensity("0");
    BOOST_CHECK_GT(vapour, 30.0);
    BOOST_CHECK_LT(vapour, 60.0);

    const Scalar liquid = columnDensity("1");
    BOOST_CHECK_GT(liquid, 400.0);
    BOOST_CHECK_LT(liquid, 800.0);
}

BOOST_AUTO_TEST_CASE(CompvdTwoZoneGasOverLiquid)
{
    // COMPVD naming both phases describes a gas zone over a liquid one meeting
    // at the gas-oil contact. Each zone takes the composition of its own rows.
    // Here the liquid column is anchored at the datum and the gas column at the
    // contact, making pressure continuous while the gradient changes there.
    const EquilFixture fix(deckString(
        "EQUIL\n 2062.5 200 2300 0 2050 0 /\n", "EQLDIMS\n/\n", "",
        "COMPVD\n"
        " 2000   0 0.95 0.05  0  150.0\n"
        " 2049   0 0.95 0.05  0  150.0\n"
        " 2051   0 0.60 0.40  1  150.0\n"
        " 2100   0 0.40 0.60  1  150.0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    for (std::size_t c = 0; c < states.size(); ++c) {
        BOOST_TEST_CONTEXT("Cell " << c) {
            const Scalar depth = fix.depths[c];
            const bool inGas = depth < 2050.0;
            BOOST_CHECK_CLOSE(states[c].saturation(inGas ? FluidSystem::gasPhaseIdx
                                                         : FluidSystem::oilPhaseIdx),
                              1.0, 1e-10);

            // The gas zone takes the vapour rows, which are constant here; the
            // liquid zone interpolates its own rows.
            const Scalar t = (depth - 2051.0) / (2100.0 - 2051.0);
            const CompVec z = inGas ? CompVec{0.0, 0.95, 0.05}
                                    : CompVec{0.0, 0.60 - 0.20 * t, 0.40 + 0.20 * t};
            for (int comp = 0; comp < 3; ++comp) {
                BOOST_CHECK_SMALL(states[c].moleFraction(comp) - z[comp], 1e-10);
            }
        }
    }

    // The datum is the centre of cell 12, so its pressure must come back exactly.
    BOOST_CHECK_CLOSE(states[12].pressure(FluidSystem::oilPhaseIdx), 200.0 * barsa, 1e-8);

    // Pressure rises monotonically through the contact, and the gas zone is far
    // lighter than the liquid one below it.
    for (std::size_t c = 0; c + 1 < states.size(); ++c) {
        BOOST_CHECK_LT(states[c].pressure(FluidSystem::oilPhaseIdx),
                       states[c + 1].pressure(FluidSystem::oilPhaseIdx));
    }
    const Scalar rhoGas = impliedDensity(states[2].pressure(FluidSystem::gasPhaseIdx),
                                         states[3].pressure(FluidSystem::gasPhaseIdx));
    const Scalar rhoLiquid = impliedDensity(states[15].pressure(FluidSystem::oilPhaseIdx),
                                            states[16].pressure(FluidSystem::oilPhaseIdx));
    BOOST_CHECK_LT(rhoGas, 0.5 * rhoLiquid);
}

BOOST_AUTO_TEST_CASE(CompvdTwoZoneDatumInGasZone)
{
    // The datum lies in the gas cap, so the gas column is the one anchored at it
    // and the liquid column picks its pressure up at the contact. Each zone is
    // integrated with the root its own COMPVD rows name, so the liquid zone keeps
    // a liquid gradient even though the datum sits above it.
    const EquilFixture fix(deckString(
        "EQUIL\n 2012.5 150 2300 0 2050 0 /\n", "EQLDIMS\n/\n", "",
        "COMPVD\n"
        " 2000   0 0.95 0.05  0  150.0\n"
        " 2049   0 0.95 0.05  0  150.0\n"
        " 2051   0 0.60 0.40  1  150.0\n"
        " 2100   0 0.40 0.60  1  150.0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();

    // The datum is the centre of cell 2, in the gas zone.
    BOOST_CHECK_CLOSE(states[2].pressure(FluidSystem::gasPhaseIdx), 150.0 * barsa, 1e-8);

    // Each zone keeps its own gradient: the liquid one stays liquid-like even
    // though the datum is above it.
    const Scalar rhoGas = impliedDensity(states[2].pressure(FluidSystem::gasPhaseIdx),
                                         states[3].pressure(FluidSystem::gasPhaseIdx));
    const Scalar rhoLiquid = impliedDensity(states[15].pressure(FluidSystem::oilPhaseIdx),
                                            states[16].pressure(FluidSystem::oilPhaseIdx));
    BOOST_CHECK_GT(rhoLiquid, 350.0);
    BOOST_CHECK_LT(rhoGas, 0.5 * rhoLiquid);
}

BOOST_AUTO_TEST_CASE(CompvdTwoZoneClampsCompositionOutsideRows)
{
    // The vapour rows cover 2020 m to 2040 m only. Gas cells outside that band
    // take the endpoint composition rather than a continued slope, as the
    // ZMFVD tables do.
    const EquilFixture fix(deckString(
        "EQUIL\n 2062.5 200 2300 0 2050 0 /\n", "EQLDIMS\n/\n", "",
        "COMPVD\n"
        " 2020   0 0.99 0.01  0  150.0\n"
        " 2040   0 0.90 0.10  0  150.0\n"
        " 2051   0 0.60 0.40  1  150.0\n"
        " 2100   0 0.40 0.60  1  150.0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();

    for (std::size_t c = 0; c < 10; ++c) {
        BOOST_TEST_CONTEXT("Cell " << c << " at " << fix.depths[c] << " m") {
            const Scalar t = std::clamp((fix.depths[c] - 2020.0) / 20.0, 0.0, 1.0);
            BOOST_CHECK_SMALL(states[c].moleFraction(1) - (0.99 - 0.09 * t), 1e-10);
            BOOST_CHECK_SMALL(states[c].moleFraction(2) - (0.01 + 0.09 * t), 1e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(WaterZoneBelowContact)
{
    // The water-oil contact at 2050 m lies inside the column. Above it the
    // hydrocarbon leaves room for the connate water of the saturation function;
    // below it the pore space holds water alone, following its own hydrostatic
    // column rather than the lighter hydrocarbon gradient.
    const WaterEquilFixture fix(waterDeckString("EQUIL\n 2010 150 2050 0 2000 0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0),
                                    std::vector<Scalar>(20, connateSw),
                                    std::vector<Scalar>(20, 1.0)).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});

    for (std::size_t c = 0; c < states.size(); ++c) {
        BOOST_TEST_CONTEXT("Cell " << c) {
            const auto& fs = states[c];
            if (fix.depths[c] > 2050.0) {
                BOOST_CHECK_CLOSE(fs.saturation(WaterFluidSystem::waterPhaseIdx), 1.0, 1e-10);
                BOOST_CHECK_SMALL(fs.saturation(WaterFluidSystem::oilPhaseIdx), 1e-10);
            }
            else {
                BOOST_CHECK_CLOSE(fs.saturation(WaterFluidSystem::waterPhaseIdx),
                                  connateSw, 1e-10);
                BOOST_CHECK_CLOSE(fs.saturation(WaterFluidSystem::oilPhaseIdx),
                                  1.0 - connateSw, 1e-10);
            }
        }
    }

    // The water column is the heavier of the two.
    const Scalar rhoHc = impliedDensity(states[2].pressure(WaterFluidSystem::oilPhaseIdx),
                                        states[3].pressure(WaterFluidSystem::oilPhaseIdx));
    const Scalar rhoWater = impliedDensity(states[15].pressure(WaterFluidSystem::waterPhaseIdx),
                                           states[16].pressure(WaterFluidSystem::waterPhaseIdx));
    BOOST_CHECK_GT(rhoWater, rhoHc);
}

BOOST_AUTO_TEST_CASE(CoincidentContactsKeepTheGasRoot)
{
    // Contacts that coincide leave no liquid: every hydrocarbon cell sits above
    // the gas-oil contact. The datum states the water pressure, so the
    // hydrocarbon is anchored on the contact itself, and the root has to follow
    // the gas just above it rather than the water below.
    const auto hydrocarbonDensity = [](const std::string& composition) {
        const WaterEquilFixture fix(waterDeckString(
            "EQUIL\n 2060 10 2050 0 2050 0 /\n", composition));
        const auto states = fix.compute(std::vector<int>(20, 0),
                                        std::vector<Scalar>(20, connateSw),
                                        std::vector<Scalar>(20, 1.0)).fluidStates();
        BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});
        return impliedDensity(states[4].pressure(WaterFluidSystem::oilPhaseIdx),
                              states[5].pressure(WaterFluidSystem::oilPhaseIdx));
    };

    // The root COMPVD states outright is the one the inferred root has to match.
    const Scalar stated = hydrocarbonDensity("COMPVD\n 2000 0 0.5 0.5 0 10.0 /\n");
    const Scalar inferred = hydrocarbonDensity("ZMFVD\n 2000 0 0.5 0.5 /\n");
    BOOST_CHECK_LT(stated, 60.0);
    BOOST_CHECK_CLOSE(inferred, stated, 1e-8);
}

BOOST_AUTO_TEST_CASE(DatumOnCoincidentContactsKeepsTheGasRoot)
{
    // The datum may sit exactly on the contacts rather than below them. That
    // depth is no more part of the hydrocarbon than one in the water, so the
    // root still has to follow the gas above it.
    const auto hydrocarbonDensity = [](const std::string& composition) {
        const WaterEquilFixture fix(waterDeckString(
            "EQUIL\n 2050 10 2050 0 2050 0 /\n", composition));
        const auto states = fix.compute(std::vector<int>(20, 0),
                                        std::vector<Scalar>(20, connateSw),
                                        std::vector<Scalar>(20, 1.0)).fluidStates();
        BOOST_REQUIRE_EQUAL(states.size(), std::size_t{20});
        return impliedDensity(states[4].pressure(WaterFluidSystem::oilPhaseIdx),
                              states[5].pressure(WaterFluidSystem::oilPhaseIdx));
    };

    const Scalar stated = hydrocarbonDensity("COMPVD\n 2000 0 0.5 0.5 0 10.0 /\n");
    const Scalar inferred = hydrocarbonDensity("ZMFVD\n 2000 0 0.5 0.5 /\n");
    BOOST_CHECK_LT(stated, 60.0);
    BOOST_CHECK_CLOSE(inferred, stated, 1e-8);
}

BOOST_AUTO_TEST_CASE(WaterEndpointsAreReadPerCell)
{
    // The endpoints come from the scaled saturation functions, so they may
    // differ from cell to cell. Each cell has to take its own, not the value
    // of some other cell.
    std::vector<Scalar> connate(20), maximum(20);
    for (std::size_t c = 0; c < connate.size(); ++c) {
        connate[c] = 0.01 + 0.002 * static_cast<Scalar>(c);
        maximum[c] = 1.00 - 0.003 * static_cast<Scalar>(c);
    }

    const WaterEquilFixture fix(waterDeckString("EQUIL\n 2010 150 2050 0 2000 0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0), connate, maximum).fluidStates();
    BOOST_REQUIRE_EQUAL(states.size(), connate.size());

    for (std::size_t c = 0; c < states.size(); ++c) {
        BOOST_TEST_CONTEXT("Cell " << c) {
            const bool belowContact = fix.depths[c] > 2050.0;
            BOOST_CHECK_CLOSE(states[c].saturation(WaterFluidSystem::waterPhaseIdx),
                              belowContact ? maximum[c] : connate[c], 1e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(WaterEndpointsMustCoverEveryCell)
{
    // A non-empty endpoint vector is read for every cell, so one of the wrong
    // length would be indexed past its end.
    const WaterEquilFixture fix(waterDeckString("EQUIL\n 2010 150 2050 0 2000 0 /\n"));
    const std::vector<int> eqlnum(20, 0);
    const std::vector<Scalar> full(20, connateSw);
    const std::vector<Scalar> tooShort(19, connateSw);

    BOOST_CHECK_NO_THROW(fix.compute(eqlnum, full, std::vector<Scalar>(20, 1.0)));

    // The endpoints share a validation block with EQLNUM, so the message has
    // to name the input that is actually wrong.
    const auto namesTheEndpoint = [](const std::string& endpoint) {
        return [endpoint](const std::runtime_error& error) {
            const std::string message = error.what();
            return (message.find(endpoint) != std::string::npos)
                && (message.find("EQLNUM") == std::string::npos);
        };
    };
    BOOST_CHECK_EXCEPTION(fix.compute(eqlnum, tooShort, std::vector<Scalar>(20, 1.0)),
                          std::runtime_error, namesTheEndpoint("connate water"));
    BOOST_CHECK_EXCEPTION(fix.compute(eqlnum, full, std::vector<Scalar>(19, 1.0)),
                          std::runtime_error, namesTheEndpoint("maximum water"));
}

BOOST_AUTO_TEST_CASE(DatumBelowWaterOilContact)
{
    // A datum below the water-oil contact states the pressure of the water, not
    // of the hydrocarbon. The water column is anchored at the datum and hands
    // the hydrocarbon its pressure at the contact, so the datum pressure must
    // appear in the water phase of the cell containing the datum (cell 12,
    // whose centre is at 2062.5 m).
    const WaterEquilFixture fix(waterDeckString("EQUIL\n 2062.5 200 2050 0 2000 0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0),
                                    std::vector<Scalar>(20, connateSw),
                                    std::vector<Scalar>(20, 1.0)).fluidStates();

    BOOST_CHECK_CLOSE(states[12].pressure(WaterFluidSystem::waterPhaseIdx),
                      200.0 * barsa, 1e-8);

    // The hydrocarbon above the contact is lighter, so its pressure has fallen
    // below the datum value, and it keeps only the connate water.
    BOOST_CHECK_LT(states[0].pressure(WaterFluidSystem::oilPhaseIdx), 200.0 * barsa);
    BOOST_CHECK_CLOSE(states[0].saturation(WaterFluidSystem::waterPhaseIdx),
                      connateSw, 1e-10);
}

BOOST_AUTO_TEST_CASE(WaterOilCapillaryPressureIsPreserved)
{
    // Retain both pressure columns even in the water zone. A separate reference
    // pressure supplies the single-pressure flash until it supports capillarity.
    // Exercise both directions of anchoring the columns at the contact.
    for (const Scalar datum : std::array{2010.0, 2062.5}) {
        BOOST_TEST_CONTEXT("Datum depth " << datum) {
            const WaterEquilFixture fix(waterDeckString(
                "EQUIL\n " + std::to_string(datum) + " 150 2052.5 10 2000 0 /\n"));
            const auto initial = fix.compute(std::vector<int>(20, 0),
                                              std::vector<Scalar>(20, connateSw),
                                              std::vector<Scalar>(20, 0.8));
            const auto& states = initial.fluidStates();
            const auto& referencePressures = initial.referencePressures();
            BOOST_REQUIRE_EQUAL(referencePressures.size(), states.size());

            const auto& contact = states[10];
            BOOST_CHECK_CLOSE(contact.pressure(WaterFluidSystem::oilPhaseIdx)
                                  - contact.pressure(WaterFluidSystem::waterPhaseIdx),
                              10.0 * barsa, 1e-6);

            // The hydrocarbon keeps its own, lighter gradient below the contact.
            const Scalar rhoHc = impliedDensity(states[15].pressure(WaterFluidSystem::oilPhaseIdx),
                                                states[16].pressure(WaterFluidSystem::oilPhaseIdx));
            const Scalar rhoWater = impliedDensity(states[15].pressure(WaterFluidSystem::waterPhaseIdx),
                                                   states[16].pressure(WaterFluidSystem::waterPhaseIdx));
            BOOST_CHECK_GT(rhoWater, rhoHc);
            for (std::size_t c = 0; c < states.size(); ++c) {
                const auto phase = fix.depths[c] > 2052.5 ? WaterFluidSystem::waterPhaseIdx
                                                         : WaterFluidSystem::oilPhaseIdx;
                BOOST_CHECK_EQUAL(referencePressures[c], states[c].pressure(phase));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(WaterContactBelowCells)
{
    // With the datum at the last cell centre, querying a deeper contact used to
    // evaluate a zero-length interval. Adding a cell at that contact must not
    // change the pressures of the original cells, for any integration path.
    const std::array decks{
        waterDeckString("EQUIL\n 2097.5 150 2200 0 2000 0 /\n"),
        waterDeckString("EQUIL\n 2097.5 150 2200 0 2097.5 0 3* 3 /\n"),
        waterDeckString("EQUIL\n 2097.5 150 2200 0 2050 0 /\n",
                        "COMPVD\n 2000 0 1 0 0 100\n"
                        " 2050 0 0.3 0.7 1 100\n"
                        " 2100 0 0.3 0.7 1 100 /\n")};
    for (std::size_t path = 0; path < decks.size(); ++path) {
        BOOST_TEST_CONTEXT("Integration path " << path) {
            WaterEquilFixture fix(decks[path]);
            const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
            fix.depths.back() = 2200.0;
            const auto expected = fix.compute(std::vector<int>(20, 0)).fluidStates();

            for (std::size_t c = 0; c + 1 < states.size(); ++c) {
                for (unsigned phase = 0; phase < WaterFluidSystem::numPhases; ++phase) {
                    BOOST_TEST_CONTEXT("Cell " << c << ", phase " << phase) {
                        BOOST_CHECK(std::isfinite(states[c].pressure(phase)));
                        BOOST_CHECK_CLOSE(states[c].pressure(phase),
                                          expected[c].pressure(phase), 1e-8);
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(WaterContactAboveCells)
{
    // The water-anchored path must also integrate all the way to an off-grid
    // contact before handing its pressure to the hydrocarbon column.
    WaterEquilFixture fix(waterDeckString("EQUIL\n 2002.5 150 1900 0 1800 0 /\n"));
    const auto states = fix.compute(std::vector<int>(20, 0)).fluidStates();
    fix.depths.front() = 1900.0;
    const auto expected = fix.compute(std::vector<int>(20, 0)).fluidStates();

    for (std::size_t c = 1; c < states.size(); ++c) {
        for (unsigned phase = 0; phase < WaterFluidSystem::numPhases; ++phase) {
            BOOST_TEST_CONTEXT("Cell " << c << ", phase " << phase) {
                BOOST_CHECK(std::isfinite(states[c].pressure(phase)));
                BOOST_CHECK_CLOSE(states[c].pressure(phase),
                                  expected[c].pressure(phase), 1e-8);
            }
        }
    }
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
