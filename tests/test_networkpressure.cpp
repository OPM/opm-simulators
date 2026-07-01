/*
  Copyright 2026 Equinor ASA

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

#define BOOST_TEST_MODULE NetworkPressureTests

#include <opm/simulators/wells/BlackoilWellModelNetworkPressureComputation.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

// The pressure computation uses a Dune parallel Communication, which requires
// MPIHelper::instance() to have been called first when built with MPI. Use the
// shared global fixture (a no-op in a serial build), which also installs an MPI
// error handler that prints a diagnostic string on failure.
#include "MpiFixture.hpp"

BOOST_GLOBAL_FIXTURE(MPIFixture);

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/Network/Branch.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Network/Node.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <map>
#include <sstream>
#include <limits>
#include <string>
#include <vector>

using namespace Opm;
using namespace Opm::unit;

// -------------------------------------------------------
// Strings with VFP tables, to be used in the test below.
// -------------------------------------------------------

// Extracted from opm-tests file network/include/vfp_gi_flowl_from_plata_to_m5s.inc
std::string vfp_string_gi3 = R"(
VFPINJ
      3        285.0         GAS /
-- gas rates Sm3/d
    5000   35227   85606  135985  186364  236742  287121  337500
  387879  488636  589394  690151  790909  891667  992424 1193939
 1395454 1596969 1798485 2000000 /

-- Tubing head pressure [bar]
  50.00 100.00 150.00 200.00 250.00 300.00 350.00 400.00
 450.00 500.00 /

  1   68.834  67.861  64.174  57.211  45.128  16.789   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  2  135.406 134.907 133.320 130.594 126.600 121.173 114.056
 104.827  92.718  49.403   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  3  201.928 201.480 200.430 198.702 196.232 192.981 188.910
 183.957 178.030 162.798 141.512 109.806  42.709   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  4  267.925 267.490 266.645 265.275 263.368 260.885 257.813
 254.126 249.813 239.150 225.493 208.238 186.338 157.525
 115.285   0.000   0.000   0.000   0.000   0.000 /

  5  333.410 333.013 332.258 331.080 329.442 327.317 324.706
 321.608 317.986 309.166 298.133 284.667 268.488 249.147
 225.851 159.496   0.000   0.000   0.000   0.000 /

  6  398.562 398.178 397.499 396.424 394.952 393.045 390.702
 387.925 384.699 376.891 367.189 355.515 341.742 325.666
 307.029 260.283 193.390   0.000   0.000   0.000 /

  7  463.483 463.125 462.485 461.499 460.130 458.363 456.200
 453.640 450.670 443.490 434.632 424.034 411.618 397.282
 380.898 341.205 289.966 220.258  86.011   0.000 /

  8  528.251 527.918 527.304 526.370 525.077 523.413 521.378
 518.971 516.194 509.474 501.192 491.323 479.803 466.581
 451.566 415.765 371.106 315.080 241.288 113.915 /

  9  592.917 592.584 592.008 591.112 589.870 588.296 586.363
 584.072 581.422 575.035 567.189 557.858 546.990 534.549
 520.494 487.240 446.434 396.770 335.726 256.968 /

 10  657.480 657.173 656.610 655.752 654.562 653.038 651.182
 648.981 646.446 640.328 632.814 623.893 613.525 601.685
 588.334 556.910 518.715 472.942 418.222 351.918 /
)";

// Extracted from opm-tests file network/include/vfp_wi_flowl_from_m5s_to_m5n.inc
// (but changed table ID to 3)
std::string vfp_string_wi3 = R"(
VFPINJ
      3        235.0         WAT /
-- water rates Sm3/d
      50     100     200     300     500     750    1000    1500
    2000    3000    5000    7500   10000   12000 /

-- Tubing head pressure [bar]
  15.00  20.00  25.00  50.00  75.00 100.00 150.00 200.00
 250.00 300.00 350.00 /

1   16.785    16.644    16.360    16.073    15.488    14.729    13.928    12.162
   10.108     4.805     0.000     0.000     0.000     0.000  /
2   21.785    21.644    21.360    21.073    20.488    19.729    18.928    17.162
   15.108     9.805     0.000     0.000     0.000     0.000  /
3   26.785    26.644    26.360    26.073    25.488    24.729    23.928    22.162
   20.108    14.805     0.000     0.000     0.000     0.000  /
4   51.785    51.644    51.360    51.073    50.488    49.729    48.928    47.162
   45.108    39.805    22.114     0.000     0.000     0.000  /
5   76.785    76.644    76.360    76.073    75.488    74.729    73.928    72.162
   70.108    64.805    47.114     4.911     0.000     0.000  /
6  101.785   101.644   101.360   101.073   100.488    99.729    98.928    97.162
   95.108    89.805    72.114    29.911     0.000     0.000  /
7  151.785   151.644   151.360   151.073   150.488   149.729   148.928   147.162
  145.108   139.805   122.114    79.911     5.769     0.000  /
8  201.785   201.644   201.360   201.073   200.488   199.729   198.928   197.162
  195.108   189.805   172.114   129.911    55.769     0.000  /
9  251.785   251.644   251.360   251.073   250.488   249.729   248.928   247.162
  245.108   239.805   222.114   179.911   105.769    16.533  /
10  301.785   301.644   301.360   301.073   300.488   299.729   298.928   297.162
  295.108   289.805   272.114   229.911   155.769    66.533  /
11  351.785   351.644   351.360   351.073   350.488   349.729   348.928   347.162
  345.108   339.805   322.114   279.911   205.769   116.533  /
)";

// Production table, hand-written to get a simple case.
std::string vfp_string_prod = R"(
VFPPROD
-- Table  Datum Depth  Rate Type  WFR Type   GFR Type   THP Type   ALQ Type    UNITS     TAB Type
     3     250.00      LIQ        WCT         GOR         THP        GRAT      METRIC   BHP      /
-- FLO: LIQ rates
       20.0       100.0    1000.0     2000.0 /
-- thp values
      10.00      30.00 /
-- WFR: WCT values
      0.000      0.5      1.0 /
-- GFR: GOR values
       100.0 /
-- ALQ: GRAT values
        0.0 /

  1  1  1  1    12.0   15.0   20.0   30.0 /
  1  2  1  1    13.0   16.0   21.0   31.0 /
  1  3  1  1    14.0   17.0   22.0   32.0 /
  2  1  1  1    32.0   35.0   40.0   50.0 /
  2  2  1  1    33.0   36.0   41.0   51.0 /
  2  3  1  1    34.0   37.0   42.0   52.0 /
/
)";

// Mock objects to represent minimally sufficient data from a single point in time.
struct MockWellModel
{
    using Scalar = double;
    struct IndexTraits
    {
        static constexpr std::size_t waterPhaseIdx = 0;
        static constexpr std::size_t oilPhaseIdx = 1;
        static constexpr std::size_t gasPhaseIdx = 2;
    };

    struct MockSchedule
    {
        struct MockGroup
        {
            const std::vector<std::string>& wells() const
            {
                static const std::vector<std::string> wells {"B-1H", "B-2H", "B-3H",
                                                             "G-3H", "G-4H",
                                                             "C-1H", "C-2H",
                                                             "F-1H", "F-2H"};
                return wells;
            }
            bool hasSatelliteProduction() const { return false; }
        };
        struct MockGSatProdValue { std::vector<double> rate {0.0}; };
        struct MockGSatProd
        {
            template <typename SummaryState>
            MockGSatProdValue get(const std::string&, const SummaryState&) const { return {}; }
        };
        struct MockScheduleStep { MockGSatProd gsatprod() const { return {}; } };
        const MockGroup& getGroup(const std::string&, int) const
        {
            static const MockGroup group {};
            return group;
        }
        const Well& getWell(const std::string&, int) const { static const Well well {}; return well; }
        MockScheduleStep operator[](int) const { return {}; }
    };

    struct MockGroupStateHelper
    {
        struct MockWellState
        {
            struct MockALQState { Scalar get() const { return 0.0; } };
            struct MockWellRates { MockALQState alq_state {}; };
            bool isOpen(const std::string&) const { return true; }
            Scalar getGlobalEfficiencyScalingFactor(const std::string&) const { return 1.0; }
            std::optional<std::size_t> index(const std::string&) const { return {}; }
            bool wellIsOwned(std::size_t, const std::string&) const { return true; }
            MockWellRates well(const std::string&) const { return {}; }
        };

        struct MockGroupState
        {
            bool has_production_rates(const std::string) const { return true; }
            std::vector<double> network_leaf_node_injection_rates(const std::string) const
            {
                // Phase order water, oil, gas.
                return {convert::from(500.0, cubic(meter) / day),
                        0.0,
                        convert::from(5000.0, cubic(meter) / day)};
            }
            std::vector<double> network_leaf_node_production_rates(const std::string) const
            {
                // Phase order water, oil, gas.
                return {convert::from(500.0, cubic(meter) / day),
                        convert::from(500.0, cubic(meter) / day),
                        convert::from(5000.0, cubic(meter) / day)};
            }
            Scalar well_group_thp(const std::string&) const { return convert::from(100.0, bars); }
        };

        MockGroupState groupState() const { return {}; }
        MockWellState wellState() const { return {}; }
    };

    MockGroupStateHelper groupStateHelper() const { return MockGroupStateHelper {}; }
    MockSchedule schedule() const { return {}; }
    int summaryState() const { return 0; }
};

enum class NetworkScenario
{
    GasInjection,
    WaterInjection,
    Production
};

std::string inputString(NetworkScenario scenario)
{
    switch (scenario) {
    case NetworkScenario::GasInjection:
        return vfp_string_gi3;
    case NetworkScenario::WaterInjection:
        return vfp_string_wi3;
    case NetworkScenario::Production:
        return vfp_string_prod;
    }
    throw std::runtime_error("Unknown scenario.");
}

double terminalPressure(NetworkScenario scenario)
{
    switch (scenario) {
    case NetworkScenario::GasInjection:
        return convert::from(350.0, bars);
    case NetworkScenario::WaterInjection:
        return convert::from(150.0, bars);
    case NetworkScenario::Production:
        return convert::from(20.0, bars);
    }
    throw std::runtime_error("Unknown scenario.");
}

struct NetworkSetup
{
    NetworkSetup(NetworkScenario scenario)
        : deck{Parser{}.parseString(inputString(scenario))}
    {
        // Set up VFP property objects.
        if (scenario == NetworkScenario::Production) {
            vfp_prod_table = VFPProdTable{deck["VFPPROD"].front(), false, UnitSystem {}};
            vfp_prod_props.addTable(vfp_prod_table);
        } else {
            vfp_inj_table = VFPInjTable{deck["VFPINJ"].front(), UnitSystem {}};
            vfp_inj_props.addTable(vfp_inj_table);
        }

        // Build a simple network.
        network.add_branch(Network::Branch{"M5S", "PLAT-A", 3, 0.0});
        network.add_branch(Network::Branch{"G1", "M5S", 9999, 0.0});
        Network::Node node{"PLAT-A"};
        node.terminal_pressure(terminalPressure(scenario));
        network.update_node(node);
    }

    Deck deck;
    VFPInjTable vfp_inj_table;
    VFPInjProperties<double> vfp_inj_props;
    VFPProdTable vfp_prod_table;
    VFPProdProperties<double> vfp_prod_props;
    MockWellModel well_model;
    Network::ExtNetwork network;
};

BOOST_AUTO_TEST_SUITE(NetworkPressureComputationTests)

BOOST_AUTO_TEST_CASE(gas_injection_pressure_computation)
{
    auto s = NetworkSetup{NetworkScenario::GasInjection};

    // Sanity check to ensure the table is read correctly.
    const auto gasrate = convert::from(5000.0, cubic(meter) / day);
    const auto thp = convert::from(100.0, bars);
    const auto expected_bhp = convert::from(135.406, bars);
    BOOST_CHECK_CLOSE(s.vfp_inj_props.bhp(3, 0.0, 0.0, gasrate, thp), expected_bhp, 1e-7);
    using Comm = Dune::Communication<int>;

    // Test using mock setup.
    NetworkPressureComputation<MockWellModel, VFPInjProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_inj_props, UnitSystem {}, 0, Comm{});
    const auto [pressures, branch_data] = comp.run();
    BOOST_REQUIRE(pressures.find("G1") != pressures.end());
    const auto expected_pressure = convert::from(463.483, bars);
    BOOST_CHECK_CLOSE(pressures.at("G1"), expected_pressure, 1e-7);
}

BOOST_AUTO_TEST_CASE(water_injection_pressure_computation)
{
    auto s = NetworkSetup{NetworkScenario::WaterInjection};

    // Sanity check to ensure the table is read correctly.
    const auto waterrate = convert::from(50.0, cubic(meter) / day);
    const auto thp = convert::from(20.0, bars);
    const auto expected_bhp = convert::from(21.785, bars);
    BOOST_CHECK_CLOSE(s.vfp_inj_props.bhp(3, waterrate, 0.0, 0.0, thp), expected_bhp, 1e-7);

    // Test using mock setup.
    using Comm = Dune::Communication<int>;
    NetworkPressureComputation<MockWellModel, VFPInjProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_inj_props, UnitSystem {}, 0, Comm{});
    const auto [pressures, branch_data] = comp.run();
    BOOST_REQUIRE(pressures.find("G1") != pressures.end());
    const auto expected_pressure = convert::from(150.488, bars);
    BOOST_CHECK_CLOSE(pressures.at("G1"), expected_pressure, 1e-7);
}

BOOST_AUTO_TEST_CASE(production_pressure_computation)
{
    auto s = NetworkSetup{NetworkScenario::Production};

    // Sanity check to ensure the table is read correctly.
    const auto oilrate = -convert::from(1000.0, cubic(meter) / day);
    const auto waterrate = oilrate;
    const auto thp = convert::from(30.0, bars);
    const auto expected_bhp = convert::from(51.0, bars);
    BOOST_CHECK_CLOSE(s.vfp_prod_props.bhp(3, waterrate, oilrate, 0.0, thp, 0.0, 0.5, 0.0, false), expected_bhp, 1e-7);

    // Test using mock setup.
    using Comm = Dune::Communication<int>;
    NetworkPressureComputation<MockWellModel, VFPProdProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_prod_props, UnitSystem {}, 0, Comm{});
    const auto [pressures, branch_data] = comp.run();
    BOOST_REQUIRE(pressures.find("G1") != pressures.end());
    const auto expected_pressure = convert::from(31.0, bars);
    BOOST_CHECK_CLOSE(pressures.at("G1"), expected_pressure, 1e-7);
}

BOOST_AUTO_TEST_SUITE_END() // NetworkPressureComputationTests
