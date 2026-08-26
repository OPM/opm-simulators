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
#include <opm/simulators/wells/NetworkNodePressureUpdater.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <boost/test/unit_test.hpp>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>

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
#include <optional>
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
            // Leaf rates in Sm3/day, phase order water, oil, gas. Tests may override
            // these before running the computation. Injection leaf rates are per network
            // phase: the gas network sees only the gas rate, the water network only water.
            static inline std::vector<double> injection_rates_sm3_day {500.0, 0.0, 5000.0};
            static inline std::vector<double> production_rates_sm3_day {500.0, 500.0, 5000.0};

            static std::vector<double> toSI(const std::vector<double>& r)
            {
                std::vector<double> out(r.size());
                std::ranges::transform(r, out.begin(), [](double v) { return convert::from(v, cubic(meter) / day); });
                return out;
            }

            bool has_production_rates(const std::string) const { return true; }
            bool has_network_leaf_node_production_rates(const std::string) const { return true; }
            // Every phase's injection rate at the leaf; the VFP table picks the one
            // its FLO type names, so it is not masked per network here.
            bool has_network_leaf_node_injection_rates(const std::string) const { return true; }
            std::vector<double> network_leaf_node_injection_rates(const std::string) const
            {
                return toSI(injection_rates_sm3_day);
            }
            std::vector<double> network_leaf_node_production_rates(const std::string) const
            {
                return toSI(production_rates_sm3_day);
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
    NetworkSetup(NetworkScenario scenario, std::optional<double> terminal_pressure_override = std::nullopt)
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
        node.terminal_pressure(terminal_pressure_override.value_or(terminalPressure(scenario)));
        network.update_node(node);
        // Restore the default leaf rates so tests do not leak state into each other.
        MockWellModel::MockGroupStateHelper::MockGroupState::injection_rates_sm3_day = {500.0, 0.0, 5000.0};
        MockWellModel::MockGroupStateHelper::MockGroupState::production_rates_sm3_day = {500.0, 500.0, 5000.0};
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

    // NetworkPressureComputation stores const references to comm and unit system, hence
    // we need to make sure that their lifetime is longer than the constructor lasts
    auto comm = Comm{};
    auto unit_system = UnitSystem {};
    // Test using mock setup.
    NetworkPressureComputation<MockWellModel, VFPInjProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_inj_props, unit_system, 0, comm);
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
    // NetworkPressureComputation stores const references to comm and unit system, hence
    // we need to make sure that their lifetime is longer than the constructor lasts
    auto comm = Comm{};
    auto unit_system = UnitSystem {};
    NetworkPressureComputation<MockWellModel, VFPInjProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_inj_props, unit_system, 0, comm);
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
    using Comm = Dune::Communication<int>;    // NetworkPressureComputation stores const references to comm and unit system, hence
    // we need to make sure that their lifetime is longer than the constructor lasts
    auto comm = Comm{};
    auto unit_system = UnitSystem {};
    NetworkPressureComputation<MockWellModel, VFPProdProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_prod_props, unit_system, 0, comm);
    const auto [pressures, branch_data] = comp.run();
    BOOST_REQUIRE(pressures.find("G1") != pressures.end());
    const auto expected_pressure = convert::from(31.0, bars);
    BOOST_CHECK_CLOSE(pressures.at("G1"), expected_pressure, 1e-7);
}

// The tables below use zero-filled cells for (rate, THP) combinations the flow line
// cannot deliver, and their axes do not cover every state the wells may be in during
// network iterations. A network branch lookup must never extrapolate into that region
// (it gives negative pressures) nor accept a zero-filled cell as a node pressure.

BOOST_AUTO_TEST_CASE(gas_injection_rate_beyond_flow_axis)
{
    auto s = NetworkSetup{NetworkScenario::GasInjection};
    // 2.5e6 Sm3/d is beyond the last flow-axis point (2.0e6); a linear extrapolation of the
    // THP=350 row (..., 86.011, 0.000) gives a pressure of about -213 bar.
    MockWellModel::MockGroupStateHelper::MockGroupState::injection_rates_sm3_day = {0.0, 0.0, 2.5e6};

    using Comm = Dune::Communication<int>;
    auto comm = Comm{};
    auto unit_system = UnitSystem {};
    NetworkPressureComputation<MockWellModel, VFPInjProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_inj_props, unit_system, 0, comm);
    const auto [pressures, branch_data] = comp.run();
    BOOST_REQUIRE(pressures.find("G1") != pressures.end());
    // Clamped to the axis end the table gives 0.0 -> no solution: the node is flagged and the
    // placeholder pressure is the upstream (terminal) pressure, never a negative value.
    BOOST_CHECK(pressures.at("M5S") >= unit::atm);
    BOOST_CHECK(pressures.at("G1") >= unit::atm);
    BOOST_CHECK(comp.invalidNodes().count("M5S") == 1);
    BOOST_CHECK(comp.invalidNodes().count("G1") == 1);
}

BOOST_AUTO_TEST_CASE(gas_injection_zero_cell_region)
{
    // At THP=100 bar the table is zero for rates >= 589394 Sm3/d.
    auto s = NetworkSetup{NetworkScenario::GasInjection, convert::from(100.0, bars)};
    MockWellModel::MockGroupStateHelper::MockGroupState::injection_rates_sm3_day = {0.0, 0.0, 6.0e5};

    using Comm = Dune::Communication<int>;
    auto comm = Comm{};
    auto unit_system = UnitSystem {};
    NetworkPressureComputation<MockWellModel, VFPInjProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_inj_props, unit_system, 0, comm);
    const auto [pressures, branch_data] = comp.run();
    BOOST_REQUIRE(pressures.find("G1") != pressures.end());
    BOOST_CHECK(pressures.at("G1") >= unit::atm);
    BOOST_CHECK(comp.invalidNodes().count("M5S") == 1);
    BOOST_CHECK(comp.invalidNodes().count("G1") == 1);
}

BOOST_AUTO_TEST_CASE(gas_injection_thp_below_axis)
{
    // Terminal pressure 20 bar is below the first THP-axis point (50 bar). Extrapolating the
    // first interval gives 68.834 - 0.6*(135.406 - 68.834) = 28.9 bar; clamping to the axis
    // gives the THP=50 row value 68.834 bar.
    auto s = NetworkSetup{NetworkScenario::GasInjection, convert::from(20.0, bars)};

    using Comm = Dune::Communication<int>;
    auto comm = Comm{};
    auto unit_system = UnitSystem {};
    NetworkPressureComputation<MockWellModel, VFPInjProperties<double>, Comm> comp(
        s.well_model, s.network, s.vfp_inj_props, unit_system, 0, comm);
    const auto [pressures, branch_data] = comp.run();
    BOOST_REQUIRE(pressures.find("G1") != pressures.end());
    BOOST_CHECK_CLOSE(pressures.at("G1"), convert::from(68.834, bars), 1e-7);
    BOOST_CHECK(comp.invalidNodes().empty());
}

BOOST_AUTO_TEST_SUITE_END() // NetworkPressureComputationTests

BOOST_AUTO_TEST_SUITE(NodePressureUpdaterTests)

namespace {
    // A synthetic well/network response modelled on the GNETINJE gas case: the wells are on
    // group control (constant rate) while the applied leaf pressure is above their THP,
    // then on THP control with a very steep rate response, then rate-limited. The network
    // pressure falls linearly with the leaf rate.
    struct StiffResponse
    {
        double q_group = 800e3;      // group-controlled rate  [sm3/d]
        double q_max = 2000e3;       // WCONINJE rate limit
        double thp_switch = 210.0;   // wells go to THP control below this leaf pressure [bar]
        double dq_dthp = 60e3;       // THP-mode rate response [sm3/d/bar] (measured ~58e3 on GNETINJE_GAS-01)
        double p0 = 500.0;           // network pressure at zero rate [bar]
        double dp_dq = -0.5e-3;      // network response [bar per sm3/d]

        double rate(double p_leaf) const
        {
            if (p_leaf >= thp_switch) {
                return q_group;
            }
            return std::clamp(q_group + dq_dthp * (p_leaf - thp_switch), 0.0, q_max);
        }
        double computed(double p_leaf) const
        {
            return p0 + dp_dq * rate(p_leaf);
        }
        // Fixed point: p = p0 + dp_dq * rate(p); on the THP branch this is linear.
        double fixed_point() const
        {
            const double a = dp_dq * dq_dthp;
            return (p0 + dp_dq * (q_group - dq_dthp * thp_switch)) / (1.0 - a);
        }
    };

    int iterationsToConverge(bool secant, const StiffResponse& r, double tol_bar, int max_it)
    {
        using namespace Opm::unit;
        NodePressureUpdater<double> updater;
        double applied = convert::from(r.p0, bars);
        for (int it = 1; it <= max_it; ++it) {
            const double computed = convert::from(r.computed(convert::to(applied, bars)), bars);
            if (std::abs(computed - applied) < convert::from(tol_bar, bars)) {
                return it;
            }
            applied = secant
                ? updater.next(applied, computed, /*valid=*/true, 0.1, convert::from(5.0, bars))
                : NodePressureUpdater<double>::damped(applied, computed - applied, 0.1, convert::from(5.0, bars));
        }
        return max_it + 1;
    }
}

BOOST_AUTO_TEST_CASE(stiff_response_converges_with_bracketing)
{
    const StiffResponse r{};
    // Loop gain |dp_dq * dq_dthp| = 30 on the THP branch, so the damped update
    // (0.1, capped at 5 bar/step) is unstable there (needs 0.1 < 2/31) and never settles.
    BOOST_CHECK_GT(iterationsToConverge(false, r, 0.5, 100), 100);
    // The bracketing update does, and to the analytic fixed point.
    const int its = iterationsToConverge(true, r, 0.5, 100);
    BOOST_CHECK_LE(its, 25);

    NodePressureUpdater<double> updater;
    double applied = convert::from(r.p0, bars);
    for (int it = 0; it < its; ++it) {
        const double computed = convert::from(r.computed(convert::to(applied, bars)), bars);
        applied = updater.next(applied, computed, true, 0.1, convert::from(5.0, bars));
    }
    BOOST_CHECK_CLOSE(convert::to(applied, bars), r.fixed_point(), 1.0);
}

BOOST_AUTO_TEST_CASE(invalid_evaluation_moves_pressure_down)
{
    NodePressureUpdater<double> updater;
    const double applied = convert::from(300.0, bars);
    // No VFP solution: treated as "pressure too high"; the update must move down and
    // not freeze the node.
    const double next = updater.next(applied, applied, /*valid=*/false, 0.1, convert::from(5.0, bars));
    BOOST_CHECK_LT(next, applied);
}

BOOST_AUTO_TEST_CASE(stale_bracket_end_is_dropped)
{
    // Seen on GNETINJE_GAS-01: a residual > 0 recorded while the wells were momentarily
    // stopped (Q=0) leaves a stale lower bracket end; all later evaluations at or above it
    // have r < 0, and the bracket must not collapse onto the stale end and stall.
    using namespace Opm::unit;
    NodePressureUpdater<double> updater;
    const auto bar = [](double v) { return convert::from(v, bars); };
    // stale lower end: applied 207.17, computed 499.4 (wells stopped)
    updater.next(bar(207.17), bar(499.4), true, 0.1, bar(5.0));
    // then the real response: computed 94.78 whenever applied >= 207.17
    double applied = bar(247.27);
    double prev = applied;
    for (int it = 0; it < 30; ++it) {
        applied = updater.next(applied, bar(94.78), true, 0.1, bar(5.0));
        BOOST_CHECK(applied <= prev);          // must keep moving down (or stay converged) ...
        prev = applied;
    }
    BOOST_CHECK(applied < bar(207.0));         // ... and get past the stale end
}

BOOST_AUTO_TEST_CASE(plateau_floor_limits_first_step)
{
    // Wells on group control at THP 207: the leaf rate, and hence the computed pressure
    // (94.8), is independent of the applied pressure down to 207 bar; below it the wells
    // become THP-limited. The first (unbracketed) step from 499 must land just below the
    // kink, not at the computed pressure deep in the dead zone.
    using namespace Opm::unit;
    NodePressureUpdater<double> updater;
    const auto bar = [](double v) { return convert::from(v, bars); };
    // No history: first step is damped/capped, but the plateau rule takes it to the kink.
    const double next = updater.next(bar(499.4), bar(94.8), true, 0.1, bar(100.0), bar(207.0));
    BOOST_CHECK_CLOSE(convert::to(next, bars), 206.5, 1e-6);
}

BOOST_AUTO_TEST_CASE(unbracketed_steps_are_capped)
{
    using namespace Opm::unit;
    NodePressureUpdater<double> updater;
    const auto bar = [](double v) { return convert::from(v, bars); };
    // Two points on a plateau (r' = -1) predict the fixed point at 94.8; without a
    // bracket the step is limited to 25% of the pressure.
    updater.next(bar(499.4), bar(94.8), true, 0.1, bar(100.0));
    const double next = updater.next(bar(459.4), bar(94.8), true, 0.1, bar(100.0));
    BOOST_CHECK_CLOSE(convert::to(next, bars), 0.75 * 459.4, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END() // NodePressureUpdaterTests
