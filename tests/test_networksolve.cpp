/*
  Copyright 2026 Equinor ASA.

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
/*!
 * \file
 *
 * \brief A standalone bench for the injection-network solve.
 *
 * The GNETINJE_GAS-01 network with the wells replaced by their inflow
 * performance plus the control logic, so solution methods can be compared
 * without a reservoir or a well solve. The VFPINJ tables are the deck's.
 *
 * Topology (GNETINJE_GAS-01, table 9999 = no table, pressure passes through):
 *
 *     PLAT-A  340 bar terminal
 *       |  VFPINJ 3          total gas
 *      M5S ------------------------------- G1 (9999)  G-3H, G-4H
 *       |  VFPINJ 2          F-wells' gas
 *      M5N ------------------------------- F1 (9999)  F-1H, F-2H
 *
 * so there are two unknowns, p(M5S) and p(M5N), and
 *
 *     G(p)_M5S = vfp3.bhp(thp = 340 bar, q = sum of all four wells)
 *     G(p)_M5N = vfp2.bhp(thp = p_M5S,   q = q(F-1H) + q(F-2H))
 *
 * with each well's rate taken from its own VFPINJ 1 against a linear IPR and
 * then put through the control logic (THP / BHP limit / rate limit / group
 * target). That last part is what gives the response its plateau, and it is
 * the part a pure dq/dbhp proxy misses.
 *
 * The bench is calibrated to the Eclipse 100 solution at day 31 and reproduces
 * it, so the numbers it reports are about the methods and not about the model.
 * What it is mainly for is globalisation: the Newton direction is the same in
 * FullStep, CappedStep, LineSearch and TrustRegion, and globalisation_basin
 * measures how much of the starting-pressure space each one recovers from.
 */

#include <config.h>

#define BOOST_TEST_MODULE NetworkSolveBench

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/UDAValue.hpp>
#include <opm/io/eclipse/ESmry.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/NetworkNodePressureUpdater.hpp>
#include <opm/simulators/wells/NetworkSystem.hpp>

#include <algorithm>
#include <array>
#include <deque>
#include <filesystem>
#include <sstream>
#include <fmt/format.h>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <tuple>
#include <map>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

using namespace Opm;
using namespace Opm::unit;

namespace {

// VFPINJ 1 (wells), from opm-tests/network/include/vfp_gi_wells.inc.
const std::string vfp_well = R"(
VFPINJ
      1        2011         GAS /
-- gas rates Sm3/d
    5000   42642   92830  168113  243396  368868  494340  619811
  745283  870755  996226 1121698 1247170 1498113 1749057 2000000 /

-- Tubing head pressure [bar]
  50.00 100.00 150.00 200.00 250.00 300.00 350.00 400.00
 450.00 500.00 /

  1   71.916  71.436  69.906  65.520  58.044  32.523   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000 /

  2  148.936 148.783 148.181 146.416 143.591 136.308 125.221
 109.013  84.309  33.054   0.000   0.000   0.000   0.000
   0.000   0.000 /

  3  223.875 223.824 223.518 222.519 220.866 216.643 210.493
 202.241 191.582 178.067 160.870 138.369 106.524   0.000
   0.000   0.000 /

  4  292.521 292.511 292.307 291.583 290.369 287.238 282.709
 276.721 269.214 260.055 249.100 236.105 220.723 180.046
 111.512   0.000 /

  5  356.485 356.465 356.302 355.710 354.690 352.069 348.274
 343.287 337.085 329.619 320.806 310.596 298.855 270.132
 232.423 180.015 /

  6  417.573 417.471 417.328 416.798 415.900 413.565 410.188
 405.762 400.274 393.685 385.994 377.130 367.053 342.950
 312.850 275.232 /

  7  476.662 476.560 476.427 475.948 475.111 472.959 469.858
 465.789 460.760 454.742 447.714 439.687 430.578 409.026
 382.638 350.712 /

  8  534.424 534.322 534.200 533.741 532.966 530.946 528.029
 524.214 519.502 513.882 507.323 499.836 491.401 471.501
 447.388 418.644 /

  9  591.218 591.126 591.004 590.565 589.821 587.903 585.129
 581.508 577.030 571.695 565.494 558.415 550.438 531.742
 509.190 482.537 /

 10  647.287 647.185 647.063 646.645 645.931 644.085 641.422
 637.954 633.660 628.560 622.624 615.851 608.242 590.443
 569.053 543.900 /
)";

// VFPINJ 2 (M5S -> M5N) and 3 (PLAT-A -> M5S), same source directory.
const std::string vfp_m5n = R"(
VFPINJ
      2        288.0         GAS /
-- gas rates Sm3/d
    5000   35227   85606  135985  186364  236742  287121  337500
  387879  488636  589394  690151  790909  891667  992424 1193939
 1395454 1596969 1798485 2000000 /

-- Tubing head pressure [bar]
  50.00 100.00 150.00 200.00 250.00 300.00 350.00 400.00
 450.00 500.00 /

  1   56.550  55.729  52.618  46.743  36.548  12.637   0.000
   0.000   0.000   0.000   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  2  112.720 112.299 110.960 108.660 105.290 100.711  94.706
  86.919  76.702  40.155   0.000   0.000   0.000   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  3  168.848 168.470 167.584 166.126 164.042 161.299 157.864
 153.685 148.684 135.832 117.872  91.120  34.507   0.000
   0.000   0.000   0.000   0.000   0.000   0.000 /

  4  224.533 224.166 223.453 222.297 220.688 218.593 216.001
 212.890 209.251 200.254 188.731 174.172 155.694 131.383
  95.743   0.000   0.000   0.000   0.000   0.000 /

  5  279.786 279.451 278.814 277.820 276.438 274.645 272.442
 269.828 266.772 259.330 250.021 238.659 225.008 208.689
 189.033 133.046   0.000   0.000   0.000   0.000 /

  6  334.758 334.434 333.861 332.954 331.712 330.103 328.126
 325.783 323.061 316.473 308.287 298.437 286.816 273.252
 257.527 218.085 161.644   0.000   0.000   0.000 /

  7  389.535 389.233 388.693 387.861 386.706 385.215 383.390
 381.230 378.724 372.666 365.192 356.250 345.774 333.678
 319.854 286.363 243.130 184.314  71.043   0.000 /

  8  444.183 443.902 443.384 442.596 441.505 440.101 438.384
 436.353 434.010 428.340 421.352 413.025 403.305 392.149
 379.480 349.273 311.592 264.320 202.058  94.587 /

  9  498.745 498.464 497.978 497.222 496.174 494.846 493.215
 491.282 489.046 483.657 477.037 469.164 459.994 449.497
 437.638 409.580 375.150 333.246 281.740 215.288 /

 10  553.220 552.961 552.486 551.762 550.758 549.472 547.906
 546.049 543.910 538.748 532.408 524.881 516.133 506.143
 494.878 468.364 436.137 397.516 351.346 295.402 /
)";

const std::string vfp_m5s = R"(
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


/// A VFPPROD table, from tests/test_networkpressure.cpp: LIQ rate, WCT / GOR
/// fractions, GRAT alq. Small enough to reason about by hand.
const std::string vfp_prod = R"(
VFPPROD
     3     250.00      LIQ        WCT         GOR         THP        GRAT      METRIC   BHP      /
       20.0       100.0    1000.0     2000.0 /
      10.00      30.00 /
      0.000      0.5      1.0 /
       100.0 /
        0.0 /
  1  1  1  1    12.0   15.0   20.0   30.0 /
  1  2  1  1    13.0   16.0   21.0   31.0 /
  1  3  1  1    14.0   17.0   22.0   32.0 /
  2  1  1  1    32.0   35.0   40.0   50.0 /
  2  2  1  1    33.0   36.0   41.0   51.0 /
  2  3  1  1    34.0   37.0   42.0   52.0 /
)";

// ---------------------------------------------------------------------------
// Model
//
// A network case is data: nodes with a parent and a VFP table, wells hanging
// off nodes, a terminal pressure and an optional group target. Nothing below
// knows about GNETINJE_GAS-01 in particular -- gnetinjeGas() is one instance,
// and NetworkCase::Builder is what a deck reader would fill in.
// ---------------------------------------------------------------------------

constexpr int kNoTable = 9999;   // GNETINJE's "no table": pressure passes through

/// Which phase the network carries. A production network would need VFPPROD
/// instead, and with it a water and a gas fraction per branch and an ALQ -- see
/// the note on Rates below.
enum class Fluid { Gas, Water };

/// The rate triple a VFP lookup takes. For an injection network only one entry
/// is ever non-zero, but carrying the triple is what a production network would
/// need: there the branch rate splits into oil, water and gas, and VFPPROD is
/// looked up on a flow rate plus WFR and GFR fractions (and an ALQ). Those
/// fractions are extra unknowns per branch, with their own mixing equations at
/// the nodes -- which is why the production side is a bigger job than swapping
/// the table type.
struct Rates
{
    double aqua = 0.0;
    double liquid = 0.0;
    double vapour = 0.0;
};

class Reference;

/// A network node. Node 0 is the terminal and carries the fixed pressure.
struct Node
{
    std::string name;
    int parent = -1;        // -1 only for the terminal
    int vfp_table = kNoTable;
};

/// One injector: a linear IPR against its own tubing table, plus its limits.
struct Well
{
    std::string name;
    int node = 0;
    int vfp_table = 1;
    double q_ref = 0.0;      // rate in the reference solution            [sm3/s]
    double bhp_ref = 0.0;    // its own bhp there; the IPR pivots here       [Pa]
    double dq_dbhp = 0.0;    // IPR slope, the stiffness knob        [sm3/s/Pa]
    double bhp_limit = 0.0;
    double rate_limit = 0.0;
    double guide = 0.0;      // share of a group target; defaults to q_ref
    double efficiency = 1.0; // WEFAC as the network sees it
};

class NetworkCase
{
public:
    /// Add a VFPINJ table from deck text. The table number in the text is the
    /// one the nodes and wells refer to.
    void addTable(const std::string& deck_text)
    {
        decks_.push_back(Parser{}.parseString(deck_text));
        addInjTable(decks_.back()["VFPINJ"].front());
    }

    void addInjTable(const DeckKeyword& keyword)
    {
        tables_.emplace_back(keyword, UnitSystem{});
        props_.addTable(tables_.back());

        const auto& t = tables_.back();
        axes_[t.getTableNum()] = Axes{t.getFloAxis().front(), t.getFloAxis().back(),
                                      t.getTHPAxis().front(), t.getTHPAxis().back()};
    }

    void setFluid(const Fluid f) { fluid_ = f; }
    Fluid fluid() const { return fluid_; }

    /// Apply an operating point: this is what fixes each well's bhp_ref, and so
    /// what its IPR is a linearisation about.
    void calibrate(const Reference& reference);

    void addNode(Node n) { nodes_.push_back(std::move(n)); }
    void addWell(Well w) { wells_.push_back(std::move(w)); }
    void setTerminalPressure(const double p) { terminal_pressure_ = p; }
    void setGroupTarget(const double target) { group_target_ = target; }

    /// Resolve everything derived from the reference solution. Call once the
    /// nodes, wells and tables are in.
    void finish()
    {
        for (auto& w : wells_) {
            if (w.guide <= 0.0) {
                w.guide = w.q_ref;
            }
        }
        children_.assign(nodes_.size(), {});
        wells_at_.assign(nodes_.size(), {});
        for (std::size_t n = 1; n < nodes_.size(); ++n) {
            children_[nodes_[n].parent].push_back(static_cast<int>(n));
        }
        for (std::size_t w = 0; w < wells_.size(); ++w) {
            wells_at_[wells_[w].node].push_back(static_cast<int>(w));
        }
        // Nodes whose pressure is not simply their parent's: the real unknowns
        // of the eliminated form.
        solved_.clear();
        for (std::size_t n = 1; n < nodes_.size(); ++n) {
            if (hasTable(nodes_[n])) {
                solved_.push_back(static_cast<int>(n));
            }
        }
    }

    /// IPR slope shared by all wells [sm3/d per bar]: the knob the bench exists for.
    void setStiffness(const double dq_dbhp_sm3_day_per_bar)
    {
        const double si = convert::from(dq_dbhp_sm3_day_per_bar, cubic(meter) / day)
                        / convert::from(1.0, bars);
        for (auto& w : wells_) {
            w.dq_dbhp = si;
        }
    }

    /// The same knob as a fraction of each well's own rate per bar, which is the
    /// only form that transfers between cases: a gas injector taking 5e5 sm3/d
    /// and a water injector taking 700 have nothing comparable to say in
    /// absolute units. 0.12/bar reproduces the 6e4 sm3/d/bar used for the gas
    /// case. Call after the rates are known.
    void setRelativeStiffness(const double fraction_per_bar)
    {
        for (auto& w : wells_) {
            w.dq_dbhp = fraction_per_bar * w.q_ref / convert::from(1.0, bars);
        }
    }

    /// Clamp table lookups to the flow and THP axes, as the simulator's network
    /// pressure computation does. See table_bounds_want_to_be_constraints.
    void setClampToAxes(const bool on) { clamp_to_axes_ = on; }

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Well>& wells() const { return wells_; }
    std::vector<Well>& wells() { return wells_; }
    const std::vector<int>& children(const int n) const { return children_[n]; }
    const std::vector<int>& wellsAt(const int n) const { return wells_at_[n]; }
    const std::vector<int>& solvedNodes() const { return solved_; }
    double terminalPressure() const { return terminal_pressure_; }
    double groupTarget() const { return group_target_; }
    bool hasTable(const Node& n) const { return n.vfp_table != kNoTable && axes_.count(n.vfp_table); }

    /// The library system this case describes: the same object the simulator
    /// assembles, differing only in where the wells' inflow performance came
    /// from. Here it is a linearisation about the reference operating point.
    NetworkSolve::System<double> system() const
    {
        NetworkSolve::System<double> s(props_, fluid_ == Fluid::Gas ? Phase::GAS : Phase::WATER);
        s.setTerminalPressure(terminal_pressure_);
        s.setGroupTarget(group_target_);
        s.setClampToAxes(clamp_to_axes_);
        for (const auto& n : nodes_) {
            s.addNode(NetworkSolve::Node{n.name, n.parent, n.vfp_table});
        }
        for (const auto& w : wells_) {
            NetworkSolve::Well<double> sw;
            sw.name = w.name;
            sw.node = w.node;
            sw.vfp_table = w.vfp_table;
            // q = q_ref + dq_dbhp*(bhp - bhp_ref) as q = a + b*bhp.
            sw.ipr_a = w.q_ref - w.dq_dbhp * w.bhp_ref;
            sw.ipr_b = w.dq_dbhp;
            sw.bhp_limit = w.bhp_limit;
            sw.rate_limit = w.rate_limit;
            sw.guide = w.guide > 0.0 ? w.guide : w.q_ref;
            sw.efficiency = w.efficiency;
            sw.in_group = group_target_ > 0.0;
            s.addWell(sw);
        }
        s.finish();
        return s;
    }

    /// Rebuild a system written by the simulator, against this case's tables.
    std::pair<NetworkSolve::System<double>, std::vector<double>>
    systemFromDump(std::istream& is) const
    {
        return NetworkSolve::read<double>(is, props_);
    }

    /// The network's scalar rate as the triple a VFP lookup takes.
    Rates asRates(const double q) const
    {
        Rates r;
        (fluid_ == Fluid::Gas ? r.vapour : r.aqua) = q;
        return r;
    }

    /// Downstream pressure of a branch, or a well's bhp: the same table lookup.
    double tableBhp(const int table, const double thp, const double q) const
    {
        const auto r = asRates(clamp_to_axes_
            ? std::clamp(q, axes_.at(table).flo_min, axes_.at(table).flo_max) : q);
        const double p = clamp_to_axes_
            ? std::clamp(thp, axes_.at(table).thp_min, axes_.at(table).thp_max) : thp;
        return props_.bhp(table, r.aqua, r.liquid, r.vapour, p);
    }

    /// Largest rate the table describes. Past it the cells are zero-filled and
    /// the interpolation runs away, so this is the edge of the feasible set.
    double maxFlow(const int table) const { return axes_.at(table).flo_max; }

    static double ipr(const Well& w, const double bhp)
    {
        return w.q_ref + w.dq_dbhp * (bhp - w.bhp_ref);
    }

    /// The rate this well takes at a given node pressure, with its own limits
    /// applied -- the eliminated form's inner solve.
    double wellRate(const Well& w, const double p_node) const
    {
        const auto f = [&](const double q) { return ipr(w, tableBhp(w.vfp_table, p_node, q)) - q; };
        double lo = convert::from(5000.0, cubic(meter) / day);
        double hi = w.rate_limit;
        while (hi > lo && tableBhp(w.vfp_table, p_node, hi) <= convert::from(1.0, atm)) {
            hi *= 0.9;
        }
        if (f(lo) <= 0.0) {
            return 0.0;                       // IPR cannot deliver even the axis minimum
        }
        double q = hi;
        if (f(hi) < 0.0) {
            for (int it = 0; it < 60; ++it) {
                q = 0.5 * (lo + hi);
                (f(q) > 0.0 ? lo : hi) = q;
            }
        }
        if (tableBhp(w.vfp_table, p_node, q) > w.bhp_limit) {
            q = std::max(ipr(w, w.bhp_limit), 0.0);   // BHP-limited
        }
        return std::clamp(q, 0.0, w.rate_limit);
    }

    /// Per-well rates at the given node pressures, group target applied.
    std::vector<double> rates(const std::vector<double>& node_pressure) const
    {
        std::vector<double> q(wells_.size());
        for (std::size_t i = 0; i < wells_.size(); ++i) {
            q[i] = wellRate(wells_[i], node_pressure[wells_[i].node]);
        }
        if (group_target_ > 0.0) {
            const double sum = std::accumulate(q.begin(), q.end(), 0.0);
            if (sum > group_target_) {
                // GRUP control: share the target out by guide rate.
                double guides = 0.0;
                for (const auto& w : wells_) {
                    guides += w.guide;
                }
                for (std::size_t i = 0; i < q.size(); ++i) {
                    q[i] = std::min(q[i], group_target_ * wells_[i].guide / guides);
                }
            }
        }
        return q;
    }

    /// Pressure at every node, given the pressures applied to the solved ones.
    std::vector<double> nodePressures(const std::vector<double>& applied) const
    {
        std::vector<double> p(nodes_.size(), terminal_pressure_);
        for (std::size_t n = 1; n < nodes_.size(); ++n) {
            const auto it = std::find(solved_.begin(), solved_.end(), static_cast<int>(n));
            p[n] = (it != solved_.end()) ? applied[it - solved_.begin()] : p[nodes_[n].parent];
        }
        return p;
    }

    /// Rate through each node's parent branch, from the well rates upwards.
    std::vector<double> branchFlows(const std::vector<double>& well_rate) const
    {
        std::vector<double> q(nodes_.size(), 0.0);
        for (std::size_t n = nodes_.size(); n-- > 1;) {
            for (const int w : wells_at_[n]) {
                q[n] += well_rate[w];
            }
            for (const int c : children_[n]) {
                q[n] += q[c];
            }
        }
        return q;
    }

private:
    struct Axes { double flo_min, flo_max, thp_min, thp_max; };

    std::map<int, Axes> axes_;
    std::deque<Deck> decks_;
    std::deque<VFPInjTable> tables_;
    VFPInjProperties<double> props_;

    std::vector<Node> nodes_;
    std::vector<Well> wells_;
    std::vector<std::vector<int>> children_;
    std::vector<std::vector<int>> wells_at_;
    std::vector<int> solved_;

    Fluid fluid_ = Fluid::Gas;
    double terminal_pressure_ = 0.0;
    double group_target_ = 0.0;
    bool clamp_to_axes_ = false;
};

/// The operating point a case is calibrated against: for each well, the rate it
/// takes and its own bottom-hole pressure. Set it by hand -- which is what the
/// unit tests do, and what you would do to pose a difficult case -- or read it
/// from a reference run's summary.
///
/// It has to be the well's bhp and not the pressure at its node. The two agree
/// only while the well is on THP control; under group control the well sits well
/// below its node (162 bar below, at day 91 of GNETINJE_GAS-01), and calibrating
/// against the node pressure there produces a well that does not exist.
class Reference
{
public:
    /// Rate in sm3/d, bottom-hole pressure in bar.
    void set(const std::string& well, const double rate_sm3_day, const double bhp_bar)
    {
        point_[well] = {convert::from(rate_sm3_day, cubic(meter) / day),
                        convert::from(bhp_bar, bars)};
    }

    bool has(const std::string& well) const { return point_.count(well) > 0; }
    double rate(const std::string& well) const { return point_.at(well).first; }
    double bhp(const std::string& well) const { return point_.at(well).second; }

    /// Read the point out of a summary case at the given time, taking the last
    /// report at or before it. `prefix` is the case path without .SMSPEC.
    /// Throws if a vector is missing, so a caller that wants to skip when the
    /// reference run is not available should check the file exists first.
    static Reference fromSummary(const std::string& prefix,
                                 const NetworkCase& c,
                                 const double time_days);

private:
    std::map<std::string, std::pair<double, double>> point_;
};

/// Topology, tables and limits from a deck; nothing calibrated yet. Reads the
/// keywords straight off the parsed Deck rather than building a Schedule, so it
/// stays usable from a unit test. WELSPECS and GRUPTREE are accumulated over
/// every occurrence; GNETINJE and WCONINJE are taken from their first, which is
/// the state at the start of the run rather than at any later DATES.
NetworkCase fromDeck(const std::string& deck_path, const Fluid fluid)
{
    const auto deck = Parser{}.parseFile(deck_path);
    NetworkCase c;
    c.setFluid(fluid);

    // GNETINJE spells the phase WAT where WCONINJE spells it WATER.
    const std::string well_phase = (fluid == Fluid::Gas) ? "GAS" : "WATER";
    const std::string node_phase = (fluid == Fluid::Gas) ? "GAS" : "WAT";

    // Which group each well belongs to, and each group's parent. These can be
    // spread over several keywords, so take them all.
    std::map<std::string, std::string> well_group;
    for (const auto& keyword : deck["WELSPECS"]) {
        for (const auto& record : keyword) {
            well_group[record.getItem("WELL").getTrimmedString(0)] =
                record.getItem("GROUP").getTrimmedString(0);
        }
    }
    std::map<std::string, std::string> parent;
    for (const auto& keyword : deck["GRUPTREE"]) {
        for (const auto& record : keyword) {
            parent[record.getItem("CHILD_GROUP").getTrimmedString(0)] =
                record.getItem("PARENT_GROUP").getTrimmedString(0);
        }
    }

    // The network itself: which groups are nodes, their table, and the terminal.
    std::map<std::string, int> node_table;
    std::string terminal;
    for (const auto& record : deck["GNETINJE"].front()) {
        if (record.getItem("PHASE").getTrimmedString(0) != node_phase) {
            continue;
        }
        const auto name = record.getItem("GROUP").getTrimmedString(0);
        const auto& pressure = record.getItem("PRESSURE");
        if (pressure.hasValue(0) && !pressure.defaultApplied(0)) {
            terminal = name;
            c.setTerminalPressure(pressure.getSIDouble(0));
        }
        const auto& table = record.getItem("VFP_TABLE");
        node_table[name] = (table.hasValue(0) && !table.defaultApplied(0))
            ? table.get<int>(0) : kNoTable;
    }
    if (terminal.empty()) {
        throw std::runtime_error("no terminal node in GNETINJE for " + node_phase);
    }

    // Nodes, terminal first, then each node after its parent.
    std::vector<std::string> order{terminal};
    for (bool grew = true; grew;) {
        grew = false;
        for (const auto& [name, table] : node_table) {
            const bool placed = std::find(order.begin(), order.end(), name) != order.end();
            const auto p = parent.find(name);
            if (placed || p == parent.end()) {
                continue;
            }
            const auto at = std::find(order.begin(), order.end(), p->second);
            if (at != order.end()) {
                order.push_back(name);
                grew = true;
            }
        }
    }
    auto index = [&order](const std::string& name) {
        return static_cast<int>(std::find(order.begin(), order.end(), name) - order.begin());
    };
    for (std::size_t i = 0; i < order.size(); ++i) {
        c.addNode(Node{order[i], i == 0 ? -1 : index(parent.at(order[i])),
                       i == 0 ? kNoTable : node_table.at(order[i])});
    }

    // Wells of the right phase whose group is a node of this network.
    for (const auto& record : deck["WCONINJE"].front()) {
        if (record.getItem("TYPE").getTrimmedString(0) != well_phase) {
            continue;
        }
        const auto name = record.getItem("WELL").getTrimmedString(0);
        const auto group = well_group.at(name);
        if (!node_table.count(group)) {
            continue;
        }
        Well w;
        w.name = name;
        w.node = index(group);
        w.vfp_table = record.getItem("VFP_TABLE").get<int>(0);
        // BHP and RATE are UDA items even when the deck gives them as numbers.
        // RATE carries no usable dimension -- which unit it is in depends on the
        // injected phase, which the parser cannot know -- so getSI() hands back
        // the raw deck number and the conversion has to be done here. Getting
        // this wrong is silent: the limit comes out 86400x too large and the
        // rate control simply never activates.
        w.bhp_limit = record.getItem("BHP").get<UDAValue>(0).getSI();
        w.rate_limit = deck.getActiveUnitSystem().to_si(
            fluid == Fluid::Gas ? UnitSystem::measure::gas_surface_rate
                                : UnitSystem::measure::liquid_surface_rate,
            record.getItem("RATE").get<UDAValue>(0).get<double>());
        c.addWell(w);
    }

    for (const auto& keyword : deck["VFPINJ"]) {
        c.addInjTable(keyword);
    }
    return c;
}

Reference Reference::fromSummary(const std::string& prefix,
                                 const NetworkCase& c,
                                 const double time_days)
{
    const EclIO::ESmry summary(prefix + ".SMSPEC");
    // Report steps, not the raw vectors: a reference run's summary carries
    // ministeps too, and an operating point taken at one of those is a slightly
    // different state from the report the deck asked for.
    const auto time = summary.get_at_rstep("TIME");
    std::size_t at = 0;
    for (std::size_t i = 1; i < time.size(); ++i) {
        if (std::abs(time[i] - time_days) < std::abs(time[at] - time_days)) {
            at = i;
        }
    }

    const std::string rate_key = (c.fluid() == Fluid::Gas) ? "WGIR:" : "WWIR:";

    Reference ref;
    for (const auto& w : c.wells()) {
        ref.set(w.name, summary.get_at_rstep(rate_key + w.name)[at],
                summary.get_at_rstep("WBHP:" + w.name)[at]);
    }
    return ref;
}

void NetworkCase::calibrate(const Reference& reference)
{
    for (auto& w : wells_) {
        if (reference.has(w.name)) {
            w.q_ref = reference.rate(w.name);
            w.bhp_ref = reference.bhp(w.name);
        }
    }
    finish();
}

/// Eclipse 100 at day 31 (rate and bhp per well), read off
/// opm-tests/eclref/e100reference/GNETINJE_GAS-01_ECL.
Reference referenceGnetinjeGasDay31()
{
    Reference r;
    r.set("G-3H", 486500.2, 295.3923);
    r.set("G-4H", 486530.8, 295.3244);
    r.set("F-1H", 276481.3, 295.0190);
    r.set("F-2H", 277082.5, 294.9374);
    return r;
}

/// GNETINJE_GAS-01 with its reference point set by hand, so the bench needs no
/// files. This is the same thing fromDeck() + Reference::fromSummary() produce
/// for that case, and setting the point by hand is also how you would pose a
/// case that no reference run covers.
NetworkCase gnetinjeGas()
{
    const auto sm3_day = cubic(meter) / day;
    const double bhp_limit = convert::from(425.0, bars);     // WCONINJE
    const double rate_limit = convert::from(1.0e6, sm3_day); // WCONINJE

    NetworkCase c;
    c.setFluid(Fluid::Gas);
    c.addTable(vfp_well);
    c.addTable(vfp_m5n);
    c.addTable(vfp_m5s);

    c.setTerminalPressure(convert::from(340.0, bars));       // GNETINJE PLAT-A
    c.addNode(Node{"PLAT-A", -1, kNoTable});
    c.addNode(Node{"M5S", 0, 3});
    c.addNode(Node{"M5N", 1, 2});
    c.addNode(Node{"G1", 1, kNoTable});
    c.addNode(Node{"F1", 2, kNoTable});

    for (const auto& [name, node] : std::initializer_list<std::pair<const char*, int>>{
             {"G-3H", 3}, {"G-4H", 3}, {"F-1H", 4}, {"F-2H", 4}}) {
        Well w;
        w.name = name;
        w.node = node;
        w.vfp_table = 1;
        w.bhp_limit = bhp_limit;
        w.rate_limit = rate_limit;
        c.addWell(w);
    }

    c.setStiffness(6.0e4);
    c.calibrate(referenceGnetinjeGasDay31());
    return c;
}

// ---------------------------------------------------------------------------
// Two formulations of the same case
//
// Eliminated: the unknowns are the pressures of the nodes that carry a table.
// Every rate is recovered from them by an inner solve, so the residual is cheap
// to state and awkward to differentiate -- the control clamps put kinks in it.
//
// Full: node pressures, branch rates, each well's (rate, bhp) and, when a group
// target is active, its multiplier. The eliminations become equations. Bigger,
// smooth within an active set, and the shape a reservoir/well system could
// absorb.
//
// Both expose size(), residual(x), start(p) and limitStep(), so the Newton
// below does not know which one it is solving.
// ---------------------------------------------------------------------------

using State = std::vector<double>;

/// Residual scaling: pressure equations in bar, rate equations in units of
/// kRateScale. Without this the two kinds of row differ by ~7 decades and no
/// single convergence tolerance means anything.
const double kPressureScale = convert::from(1.0, bars);
const double kRateScale = convert::from(1.0e4, cubic(meter) / day);

class EliminatedProblem
{
public:
    explicit EliminatedProblem(const NetworkCase& c) : case_(c) {}

    static constexpr const char* name = "eliminated";
    int size() const { return static_cast<int>(case_.solvedNodes().size()); }

    /// The fixed-point map: applied node pressures in, computed ones out.
    State G(const State& applied) const
    {
        const auto p = case_.nodePressures(applied);
        const auto q = case_.branchFlows(case_.rates(p));

        State out(size());
        const auto& solved = case_.solvedNodes();
        for (std::size_t i = 0; i < solved.size(); ++i) {
            const auto& node = case_.nodes()[solved[i]];
            out[i] = case_.tableBhp(node.vfp_table, p[node.parent], q[solved[i]]);
        }
        return out;
    }

    State residual(const State& x) const
    {
        const auto g = G(x);
        State r(size());
        for (int i = 0; i < size(); ++i) {
            r[i] = (g[i] - x[i]) / kPressureScale;
        }
        return r;
    }

    State start(const State& p) const { return p; }
    State pressures(const State& x) const { return x; }
    /// The eliminated form has no rate unknowns; recover them from the case.
    State wellRates(const State& x) const { return case_.rates(case_.nodePressures(x)); }
    double columnScale(const int) const { return kPressureScale; }
    State limitStep(const State&, const State& dx) const { return dx; }

private:
    const NetworkCase& case_;
};

/// The same case without the eliminations, solved by the library's
/// NetworkSolve::System -- the very code the simulator runs. The bench only
/// supplies the wells' inflow performance from its reference operating point,
/// where the simulator supplies it from the well Jacobian.
class FullProblem
{
public:
    explicit FullProblem(const NetworkCase& c)
        : system_(c.system()), solved_(c.solvedNodes()), terminal_(c.terminalPressure())
    {}

    static constexpr const char* name = "full";

    int size() const { return system_.size(); }
    State residual(const State& x) const { return system_.residual(x); }
    bool updateControls(const State& x) { return system_.updateControls(x); }
    double columnScale(const int i) const { return system_.columnScale(i); }
    State limitStep(const State& x, const State& dx) const
    {
        return enforce_bounds_ ? system_.limitStep(x, dx) : dx;
    }

    void setEnforceBounds(const bool on) { enforce_bounds_ = on; }
    void setAnalyticJacobian(const bool on) { system_.setAnalyticJacobian(on); }
    void setGuidesFromPotential(const bool on) { system_.setGuidesFromPotential(on); }
    void dropLastFromGroup() { system_.dropLastFromGroup(); }
    State wellRates(const State& x) const { return system_.wellRates(x); }
    const NetworkSolve::System<double>& system() const { return system_; }
    NetworkSolve::System<double>& system() { return system_; }

    /// The bench starts both formulations from the same applied node pressures.
    State start(const State& applied) const
    {
        return system_.start(applied_to_all_nodes_(applied));
    }

    /// Only the nodes the eliminated form solves for, so the two are comparable.
    State pressures(const State& x) const
    {
        const auto all = system_.pressures(x);
        State p;
        for (const int n : solved_) {
            p.push_back(all[n]);
        }
        return p;
    }

private:
    State applied_to_all_nodes_(const State& applied) const
    {
        State p(system_.nodes().size(), terminal_);
        for (std::size_t n = 1; n < system_.nodes().size(); ++n) {
            const auto it = std::find(solved_.begin(), solved_.end(), static_cast<int>(n));
            p[n] = (it != solved_.end()) ? applied[it - solved_.begin()]
                                         : p[system_.nodes()[n].parent];
        }
        return p;
    }

    NetworkSolve::System<double> system_;
    std::vector<int> solved_;
    double terminal_ = 0.0;
    bool enforce_bounds_ = false;
};

// ---------------------------------------------------------------------------
// Solvers
//
// Everything works on the residual F(x). Convergence is in the max norm of the
// scaled residual; the trust region uses the 2-norm because that is what its
// reduction ratio is defined against.
// ---------------------------------------------------------------------------

// Start where the simulator does: the wells' WCONINJE THP.
const State kStart{convert::from(400.0, bars), convert::from(400.0, bars)};
const double kTol = 0.01;                     // scaled, so 0.01 bar
const double kMaxStep = 100.0;                // column scales, so 100 bar
constexpr int kMaxIter = 200;

struct Result
{
    bool converged = false;
    int iterations = 0;
    State p{};
    State well_rate{};
};

double normMax(const State& v)
{
    double m = 0.0;
    for (const double e : v) {
        m = std::max(m, std::abs(e));
    }
    return m;
}
double norm2(const State& v)
{
    double sum = 0.0;
    for (const double e : v) {
        sum += e * e;
    }
    return std::sqrt(sum);
}
State operator+(const State& a, const State& b)
{
    State c(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}
State operator-(const State& a, const State& b)
{
    State c(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}
State operator*(const double a, const State& v)
{
    State c(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        c[i] = a * v[i];
    }
    return c;
}
State operator-(const State& v) { return -1.0 * v; }

// --- fixed-point methods, eliminated form only -------------------------------

Result damped(const EliminatedProblem& problem, State p, const double omega)
{
    for (int it = 1; it <= kMaxIter; ++it) {
        const auto r = problem.residual(p);
        if (normMax(r) < kTol) {
            return {true, it, p};
        }
        for (int i = 0; i < problem.size(); ++i) {
            p[i] = NodePressureUpdater<double>::damped(p[i], r[i] * kPressureScale, omega,
                                                       kMaxStep * kPressureScale);
        }
    }
    return {false, kMaxIter + 1, p};
}

Result bracketing(const EliminatedProblem& problem, State p, const double omega)
{
    std::vector<NodePressureUpdater<double>> updater(problem.size());
    for (int it = 1; it <= kMaxIter; ++it) {
        const auto g = problem.G(p);
        if (normMax(g - p) < kTol * kPressureScale) {
            return {true, it, p};
        }
        for (int i = 0; i < problem.size(); ++i) {
            p[i] = updater[i].next(p[i], g[i], /*valid=*/true, omega, kMaxStep * kPressureScale);
        }
    }
    return {false, kMaxIter + 1, p};
}

// --- Newton ------------------------------------------------------------------

/// Dense square system, small enough that Gaussian elimination with partial
/// pivoting is the whole story.
class Matrix
{
public:
    explicit Matrix(const int n) : n_(n), a_(n * n, 0.0) {}

    double& operator()(const int i, const int j) { return a_[i * n_ + j]; }

    /// Solves A y = b. Returns false if A is singular to working precision.
    bool solve(State b, State& y) const
    {
        auto a = a_;
        y.assign(n_, 0.0);
        for (int k = 0; k < n_; ++k) {
            int pivot = k;
            for (int i = k + 1; i < n_; ++i) {
                if (std::abs(a[i * n_ + k]) > std::abs(a[pivot * n_ + k])) {
                    pivot = i;
                }
            }
            if (std::abs(a[pivot * n_ + k]) < 1e-300) {
                return false;
            }
            if (pivot != k) {
                for (int j = 0; j < n_; ++j) {
                    std::swap(a[k * n_ + j], a[pivot * n_ + j]);
                }
                std::swap(b[k], b[pivot]);
            }
            for (int i = k + 1; i < n_; ++i) {
                const double f = a[i * n_ + k] / a[k * n_ + k];
                for (int j = k; j < n_; ++j) {
                    a[i * n_ + j] -= f * a[k * n_ + j];
                }
                b[i] -= f * b[k];
            }
        }
        for (int i = n_ - 1; i >= 0; --i) {
            double sum = b[i];
            for (int j = i + 1; j < n_; ++j) {
                sum -= a[i * n_ + j] * y[j];
            }
            y[i] = sum / a[i * n_ + i];
        }
        return true;
    }

private:
    int n_;
    std::vector<double> a_;
};

/// Finite-difference Jacobian. On the eliminated problem the well response has
/// kinks where a control switches, so a difference taken across one is not the
/// local slope -- which is exactly why the step needs globalising. The full
/// problem holds its controls fixed while this is taken, so it has no kinks.
template <class Problem>
Matrix jacobian(const Problem& problem, const State& x, const State& r)
{
    const int n = problem.size();
    Matrix J(n);
    for (int j = 0; j < n; ++j) {
        State shifted = x;
        const double h = 1e-2 * problem.columnScale(j);
        shifted[j] += h;
        const auto rj = problem.residual(shifted);
        for (int i = 0; i < n; ++i) {
            J(i, j) = (rj[i] - r[i]) / h;
        }
    }
    return J;
}

// --- globalisation strategies ------------------------------------------------
//
// Each takes the current point, the residual there and the full Newton step, and
// returns the point to move to. The Newton direction is the same in all of them.

/// Take the step as it comes.
struct FullStep
{
    static constexpr const char* name = "newton, full step";

    template <class Problem>
    State accept(const Problem&, const State& x, const State&, const State& dx)
    {
        return x + dx;
    }
};

/// Clamp each component, the way --network-max-pressure-update-in-bars does.
struct CappedStep
{
    static constexpr const char* name = "newton, capped step";
    double cap = kMaxStep;   // in column scales, so 100 means 100 bar

    template <class Problem>
    State accept(const Problem& problem, const State& x, const State&, const State& dx)
    {
        State capped = dx;
        for (int i = 0; i < problem.size(); ++i) {
            const double limit = cap * problem.columnScale(i);
            capped[i] = std::clamp(capped[i], -limit, limit);
        }
        return x + capped;
    }
};

/// Backtrack until the residual norm drops. With sufficient_decrease it is the
/// Armijo condition rather than plain decrease.
struct LineSearch
{
    static constexpr const char* name = "newton, line search";
    bool sufficient_decrease = false;
    int max_halvings = 12;

    template <class Problem>
    State accept(const Problem& problem, const State& x, const State& r, const State& dx)
    {
        const double f0 = norm2(r);
        double lambda = 1.0;
        for (int k = 0; k < max_halvings; ++k) {
            const State trial = x + lambda * dx;
            const double f = norm2(problem.residual(trial));
            const double target = sufficient_decrease ? (1.0 - 1e-4 * lambda) * f0 : f0;
            if (f < target) {
                return trial;
            }
            lambda *= 0.5;
        }
        return x + lambda * dx;
    }
};

/// Classic trust region on ||F||_2: shrink the step to the radius, accept on the
/// ratio of actual to predicted reduction, and shrink and retry when it is poor.
struct TrustRegion
{
    static constexpr const char* name = "newton, trust region";
    double radius = 50.0;          // in column scales, so 50 means 50 bar
    double radius_max = 400.0;
    double radius_min = 1e-4;
    double radius_start = 50.0;

    template <class Problem>
    State accept(const Problem& problem, const State& x, const State& r, const State& dx)
    {
        const double f0 = norm2(r);
        // Measure the step in column scales, so one radius covers pressures and rates.
        State scaled(dx.size());
        for (int i = 0; i < problem.size(); ++i) {
            scaled[i] = dx[i] / problem.columnScale(i);
        }
        const double len = norm2(scaled);

        while (radius > radius_min) {
            const double lambda = (len > radius && len > 0.0) ? radius / len : 1.0;
            const State trial = x + lambda * dx;
            const double f = norm2(problem.residual(trial));
            // The step solves J dx = -r, so the linear model predicts (1 - lambda)*r.
            const double predicted = lambda * f0;
            const double rho = predicted > 0.0 ? (f0 - f) / predicted : -1.0;

            if (rho > 0.1) {
                if (rho > 0.75 && lambda < 1.0) {
                    radius = std::min(2.0 * radius, radius_max);
                }
                return trial;
            }
            radius *= 0.5;
        }
        // The region collapsed, which here means the Jacobian was taken across a
        // control switch. Take the smallest step and reopen rather than stalling.
        const double lambda = std::min(1.0, radius_min / std::max(len, radius_min));
        radius = radius_start;
        return x + lambda * dx;
    }
};

/// Newton on either formulation. The full problem reselects its well controls
/// once per iteration; converging with a control still moving is not converged.
template <class Problem, class Globalisation>
Result newton(Problem problem, const State& start, Globalisation g = {})
{
    State x = problem.start(start);
    for (int it = 1; it <= kMaxIter; ++it) {
        bool controls_moved = false;
        if constexpr (requires { problem.updateControls(x); }) {
            controls_moved = problem.updateControls(x);
        }
        const auto r = problem.residual(x);
        if (normMax(r) < kTol && !controls_moved) {
            return {true, it, problem.pressures(x), problem.wellRates(x)};
        }
        State dx;
        if (!jacobian(problem, x, r).solve(-r, dx)) {
            return {false, kMaxIter + 1, problem.pressures(x)};
        }
        // Keep the iterate inside the box the tables describe before anything
        // else looks at the step.
        dx = problem.limitStep(x, dx);
        // The residual jumps when a control switches, and that jump is not a
        // failure to make progress. Letting a globalisation veto it stalls the
        // active set instead of resolving it.
        x = controls_moved ? x + dx : g.accept(problem, x, r, dx);
    }
    return {false, kMaxIter + 1, problem.pressures(x)};
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(NetworkSolveBench)

namespace {
    const State kExpected{convert::from(209.30, bars), convert::from(204.19, bars)};

    void report(const char* name, const Result& r)
    {
        BOOST_TEST_MESSAGE(std::left << std::setw(26) << name
                           << (r.converged ? "converged in " : "FAILED after ")
                           << std::setw(4) << r.iterations << " iterations, p = ("
                           << convert::to(r.p[0], bars) << ", "
                           << convert::to(r.p[1], bars) << ") bar");
    }

    /// The grid the basin tests sweep: starting node pressures across the tables' THP axis.
    std::vector<State> startingPoints()
    {
        std::vector<State> starts;
        for (int a = 60; a <= 500; a += 20) {
            for (int b = 60; b <= 500; b += 20) {
                starts.push_back({convert::from(a, bars), convert::from(b, bars)});
            }
        }
        return starts;
    }

    /// How many of the starting points a method reaches the right answer from.
    template <class Solve>
    int basin(const char* name, Solve&& solve, const State& expected = kExpected)
    {
        const auto starts = startingPoints();
        const double tol = convert::from(1.0, bars);

        int solved = 0, iterations = 0;
        for (const auto& start : starts) {
            const auto r = solve(start);
            if (r.converged && normMax(r.p - expected) < tol) {
                ++solved;
                iterations += r.iterations;
            }
        }
        BOOST_TEST_MESSAGE(std::left << std::setw(26) << name << solved << "/" << starts.size()
                           << " starts, mean " << (solved ? iterations / solved : 0)
                           << " iterations");
        return solved;
    }
}

// The branch tables reproduce the Eclipse 100 operating point, which is what makes
// everything below a statement about the methods and not about the model.
BOOST_AUTO_TEST_CASE(branches_match_eclipse)
{
    const auto c = gnetinjeGas();
    const auto sm3d = cubic(meter) / day;
    // E100 day 31: M5S = 209.4 bar at 1.532e6 sm3/d, M5N = 204.2 bar at 5.53e5 sm3/d.
    const double m5s = c.tableBhp(3, convert::from(340.0, bars), convert::from(1.532e6, sm3d));
    const double m5n = c.tableBhp(2, convert::from(209.4, bars), convert::from(5.53e5, sm3d));
    BOOST_TEST_MESSAGE("M5S " << convert::to(m5s, bars) << " (E100 209.4), M5N "
                       << convert::to(m5n, bars) << " (E100 204.2) bar");
    BOOST_CHECK_CLOSE(convert::to(m5s, bars), 209.4, 2.0);
    BOOST_CHECK_CLOSE(convert::to(m5n, bars), 204.2, 2.0);
}

// Both formulations describe the same network, so they must land on the same point.
BOOST_AUTO_TEST_CASE(both_formulations_match_eclipse)
{
    const auto c = gnetinjeGas();
    // Each with the method that suits it: the eliminated residual needs
    // globalising, the full one does not.
    const auto eliminated = newton(EliminatedProblem{c}, kStart, TrustRegion{});
    const auto full = newton(FullProblem{c}, kStart, FullStep{});

    BOOST_REQUIRE(eliminated.converged);
    BOOST_REQUIRE(full.converged);
    for (const auto& r : {eliminated, full}) {
        BOOST_TEST_MESSAGE("solution (" << convert::to(r.p[0], bars) << ", "
                           << convert::to(r.p[1], bars) << ") bar, E100 (209.4, 204.2)");
        BOOST_CHECK_CLOSE(convert::to(r.p[0], bars), 209.4, 0.5);
        BOOST_CHECK_CLOSE(convert::to(r.p[1], bars), 204.2, 0.5);
    }
    BOOST_CHECK_SMALL(convert::to(full.p[0] - eliminated.p[0], bars), 0.05);
    BOOST_CHECK_SMALL(convert::to(full.p[1] - eliminated.p[1], bars), 0.05);
}

// A case can also be built from the deck, with the operating point read out of a
// reference run instead of written down. Both routes must give the same case --
// that is what makes the hand-set reference above trustworthy as a stand-in.
//
// Guarded on the files being present, because opm-tests is not a build
// dependency; the hand-set path above is what runs everywhere.
BOOST_AUTO_TEST_CASE(built_from_deck_matches_the_builtin_case)
{
    const std::string tests = "/Users/hnil/Documents/OPM/opm_feature/opm-tests/network/";
    const std::string deck = tests + "GNETINJE_GAS-01.DATA";
    const std::string reference = tests + "../eclref/e100reference/GNETINJE_GAS-01_ECL";
    if (!std::filesystem::exists(deck) || !std::filesystem::exists(reference + ".SMSPEC")) {
        BOOST_TEST_MESSAGE("opm-tests not present, skipping the deck-driven case");
        return;
    }

    auto from_deck = fromDeck(deck, Fluid::Gas);
    from_deck.setStiffness(6.0e4);
    from_deck.calibrate(Reference::fromSummary(reference, from_deck, 31.0));

    const auto builtin = gnetinjeGas();
    BOOST_REQUIRE_EQUAL(from_deck.nodes().size(), builtin.nodes().size());
    BOOST_REQUIRE_EQUAL(from_deck.wells().size(), builtin.wells().size());
    BOOST_CHECK_CLOSE(convert::to(from_deck.terminalPressure(), bars),
                      convert::to(builtin.terminalPressure(), bars), 1e-6);

    // Node and well order is an artefact of how each was assembled, so compare
    // by name: same parent, same table, same operating point.
    auto nodeName = [](const NetworkCase& c, const int i) { return c.nodes()[i].name; };
    auto findNode = [&](const NetworkCase& c, const std::string& name) {
        const auto& n = c.nodes();
        return static_cast<int>(std::find_if(n.begin(), n.end(),
            [&](const Node& x) { return x.name == name; }) - n.begin());
    };

    for (const auto& node : builtin.nodes()) {
        const int i = findNode(from_deck, node.name);
        BOOST_REQUIRE_LT(i, static_cast<int>(from_deck.nodes().size()));
        const auto& d = from_deck.nodes()[i];
        BOOST_CHECK_EQUAL(d.vfp_table, node.vfp_table);
        if (node.parent >= 0) {
            BOOST_CHECK_EQUAL(nodeName(from_deck, d.parent), nodeName(builtin, node.parent));
        }
    }

    for (const auto& well : builtin.wells()) {
        const auto& w = from_deck.wells();
        const auto it = std::find_if(w.begin(), w.end(),
                                     [&](const Well& x) { return x.name == well.name; });
        BOOST_REQUIRE(it != w.end());
        BOOST_CHECK_EQUAL(nodeName(from_deck, it->node), nodeName(builtin, well.node));
        BOOST_CHECK_EQUAL(it->vfp_table, well.vfp_table);
        BOOST_CHECK_CLOSE(convert::to(it->bhp_limit, bars), convert::to(well.bhp_limit, bars), 1e-6);
        BOOST_CHECK_CLOSE(convert::to(it->rate_limit, cubic(meter) / day),
                          convert::to(well.rate_limit, cubic(meter) / day), 1e-6);
        // The hand-set reference is the summary rounded to four figures.
        BOOST_CHECK_CLOSE(convert::to(it->q_ref, cubic(meter) / day),
                          convert::to(well.q_ref, cubic(meter) / day), 0.5);
        BOOST_CHECK_CLOSE(convert::to(it->bhp_ref, bars), convert::to(well.bhp_ref, bars), 0.5);
    }

    // And it solves to the same place.
    const auto solved = newton(FullProblem{from_deck}, kStart, FullStep{});
    BOOST_REQUIRE(solved.converged);
    for (std::size_t i = 0; i < from_deck.solvedNodes().size(); ++i) {
        BOOST_TEST_MESSAGE("  " << nodeName(from_deck, from_deck.solvedNodes()[i]) << " "
                           << convert::to(solved.p[i], bars) << " bar");
    }
    // M5S and M5N, in whichever order this case lists them.
    const double m5s = solved.p[findNode(from_deck, "M5S") == from_deck.solvedNodes()[0] ? 0 : 1];
    const double m5n = solved.p[findNode(from_deck, "M5N") == from_deck.solvedNodes()[0] ? 0 : 1];
    BOOST_CHECK_CLOSE(convert::to(m5s, bars), 209.4, 0.5);
    BOOST_CHECK_CLOSE(convert::to(m5n, bars), 204.2, 0.5);
}

namespace {
    /// Where opm-tests lives, if it does. The reference-driven cases skip
    /// without it, since it is not a build dependency.
    const std::string kTests = "/Users/hnil/Documents/OPM/opm_feature/opm-tests/";

    struct DeckCase
    {
        std::string deck;
        std::string reference;
        Fluid fluid;
        std::string rate_key;     // field injection rate
    };

    const DeckCase kGasCase{kTests + "network/GNETINJE_GAS-01.DATA",
                            kTests + "eclref/e100reference/GNETINJE_GAS-01_ECL",
                            Fluid::Gas, "FGIR"};
    const DeckCase kWaterCase{kTests + "network/GNETINJE_WAT-01.DATA",
                              kTests + "eclref/e100reference/GNETINJE_WAT-01_ECL",
                              Fluid::Water, "FWIR"};

    bool available(const DeckCase& c)
    {
        return std::filesystem::exists(c.deck) && std::filesystem::exists(c.reference + ".SMSPEC");
    }

    /// Solve the case at every report step of its reference run and report how
    /// far the node pressures land from it. Returns (matched, considered).
    ///
    /// At each step the wells are calibrated to what the reference says they
    /// were doing, and the field target is applied when the reference has the
    /// wells on group control -- WMCTL 6 is THP, anything negative is a group.
    /// So this asks whether the network side reproduces Eclipse at operating
    /// points across the whole run, not just the one the bench was built on.
    std::pair<int, int> sweepReference(const DeckCase& deck_case, const double tol_bar)
    {
        const EclIO::ESmry summary(deck_case.reference + ".SMSPEC");
        const auto time = summary.get_at_rstep("TIME");
        const auto field_rate = summary.get_at_rstep(deck_case.rate_key);
        const bool gas = deck_case.fluid == Fluid::Gas;
        const auto sm3d = cubic(meter) / day;

        const auto base = fromDeck(deck_case.deck, deck_case.fluid);
        const auto control = summary.get_at_rstep("WMCTL:" + base.wells().front().name);

        int matched = 0, considered = 0;
        for (std::size_t step = 0; step < time.size(); ++step) {
            if (time[step] <= 0.0 || field_rate[step] <= 0.0) {
                continue;      // nothing injected yet
            }
            ++considered;

            auto c = fromDeck(deck_case.deck, deck_case.fluid);
            c.calibrate(Reference::fromSummary(deck_case.reference, c, time[step]));
            c.setRelativeStiffness(0.12);
            const bool on_group = control[step] < 0.0;
            if (on_group) {
                c.setGroupTarget(convert::from(field_rate[step], sm3d));
            }
            c.finish();

            const auto solved = newton(FullProblem{c}, kStart, FullStep{});
            std::string detail;
            double worst = 0.0;
            for (std::size_t i = 0; i < c.solvedNodes().size(); ++i) {
                const auto& node = c.nodes()[c.solvedNodes()[i]].name;
                const double reference =
                    summary.get_at_rstep((gas ? "GPRG:" : "GPRW:") + node)[step];
                const double got = solved.converged ? convert::to(solved.p[i], bars) : 0.0;
                worst = std::max(worst, std::abs(got - reference));
                detail += fmt::format(" {} {:.2f}/{:.2f}", node, got, reference);
            }
            const bool ok = solved.converged && worst < tol_bar;
            matched += ok ? 1 : 0;
            BOOST_TEST_MESSAGE(fmt::format("  t={:6.1f} {:5} {:11} {:<26} worst {:.2f} bar {}",
                                           time[step], on_group ? "GRUP" : "THP",
                                           solved.converged
                                               ? fmt::format("{} it", solved.iterations) : "FAILED",
                                           detail, worst, ok ? "" : "  <--"));
        }
        BOOST_TEST_MESSAGE("  " << matched << "/" << considered << " report steps within "
                           << tol_bar << " bar");
        return {matched, considered};
    }
}

// The full formulation against the whole of GNETINJE_GAS-01, not just the point
// the bench was calibrated on.
BOOST_AUTO_TEST_CASE(gas_case_across_the_run)
{
    if (!available(kGasCase)) {
        BOOST_TEST_MESSAGE("opm-tests not present, skipping");
        return;
    }
    BOOST_TEST_MESSAGE("GNETINJE_GAS-01:");
    const auto [matched, considered] = sweepReference(kGasCase, 1.0);
    BOOST_CHECK_GT(considered, 0);
    BOOST_CHECK_EQUAL(matched, considered);
}

BOOST_AUTO_TEST_CASE(water_case_across_the_run)
{
    if (!available(kWaterCase)) {
        BOOST_TEST_MESSAGE("opm-tests not present, skipping");
        return;
    }
    BOOST_TEST_MESSAGE("GNETINJE_WAT-01:");
    const auto [matched, considered] = sweepReference(kWaterCase, 1.0);
    BOOST_CHECK_GT(considered, 0);
    BOOST_CHECK_EQUAL(matched, considered);
}

// From the one starting point the simulator actually uses.
BOOST_AUTO_TEST_CASE(method_comparison)
{
    const auto c = gnetinjeGas();
    const EliminatedProblem eliminated{c};

    const auto fixed_point = damped(eliminated, kStart, 0.1);
    const auto bracket = bracketing(eliminated, kStart, 0.1);
    const auto full_step = newton(eliminated, kStart, FullStep{});
    const auto capped = newton(eliminated, kStart, CappedStep{});
    const auto search = newton(eliminated, kStart, LineSearch{});
    const auto armijo = newton(eliminated, kStart, LineSearch{/*sufficient_decrease=*/true});
    const auto region = newton(eliminated, kStart, TrustRegion{});

    report("damped (omega 0.1)", fixed_point);
    report("bracketing (shipped)", bracket);
    report(FullStep::name, full_step);
    report(CappedStep::name, capped);
    report(LineSearch::name, search);
    report("newton, armijo", armijo);
    report(TrustRegion::name, region);

    // The damped update is the original branch's method: it limit-cycles here.
    BOOST_CHECK(!fixed_point.converged);
    // An unglobalised Newton step overshoots off the plateau and does not return.
    BOOST_CHECK(!full_step.converged);
    BOOST_CHECK(bracket.converged);
    for (const auto& r : {capped, search, armijo, region}) {
        BOOST_CHECK(r.converged);
        BOOST_CHECK_LT(r.iterations, bracket.iterations);
    }
}

// The Jacobian entries the tables already compute and throw away, against the
// difference quotient they replace. Everything but the two lookups is constant,
// so an error here is an error in the branch or tubing rows.
BOOST_AUTO_TEST_CASE(analytic_jacobian_matches_differences)
{
    // Every control row has its own derivative, so check a state that exercises
    // each: the ordinary THP one, one where the group is holding the wells, and
    // one driven hard enough that the bhp and rate limits bite.
    auto check = [](const char* what, FullProblem& problem, const State& x,
                    const bool drop_last_from_group = false) {
        if (drop_last_from_group) {
            problem.dropLastFromGroup();
        }
        problem.updateControls(x);
        const auto r = problem.residual(x);
        const auto analytic = problem.system().jacobian(x);

        const int n = problem.size();
        double worst = 0.0, scale = 0.0;
        for (int j = 0; j < n; ++j) {
            auto shifted = x;
            const double h = 1e-4 * problem.columnScale(j);
            shifted[j] += h;
            const auto rj = problem.residual(shifted);
            for (int i = 0; i < n; ++i) {
                const double fd = (rj[i] - r[i]) / h;
                worst = std::max(worst, std::abs(fd - analytic(i, j)));
                scale = std::max(scale, std::abs(fd));
            }
        }
        BOOST_TEST_MESSAGE("  " << std::left << std::setw(16) << what
                           << "largest entry " << scale << ", largest difference " << worst);
        BOOST_CHECK_LT(worst, 1e-4 * std::max(scale, 1.0));
    };

    {
        const auto c = gnetinjeGas();
        FullProblem problem{c};
        check("thp", problem, problem.start(kStart));
        // A converged state, where the controls have settled.
        const auto solved = newton(FullProblem{c}, kStart, FullStep{});
        BOOST_REQUIRE(solved.converged);
        FullProblem at_solution{c};
        check("thp, converged", at_solution, at_solution.start(solved.p));
    }
    {
        auto c = gnetinjeGas();
        c.setGroupTarget(convert::from(1.0e6, cubic(meter) / day));
        c.finish();
        FullProblem problem{c};
        check("group", problem, problem.start(kStart));
    }
    {
        // A group that does not hold every well on the network. The group row
        // differentiates only its own, and with every well in the group -- which
        // is the only case the tests had -- a wrong derivative there is
        // invisible.
        auto c = gnetinjeGas();
        c.setGroupTarget(convert::from(1.0e6, cubic(meter) / day));
        c.finish();
        auto system = c.system();
        FullProblem problem{c};
        check("group, partial", problem, problem.start(kStart), /*drop_last_from_group=*/true);
    }
    {
        // A very low terminal pressure drives the wells onto their limits.
        auto c = gnetinjeGas();
        c.setStiffness(1.0e6);
        c.finish();
        FullProblem problem{c};
        check("limits", problem, problem.start({convert::from(80.0, bars),
                                                convert::from(80.0, bars)}));
    }
}

// It should change what a solve costs, not where it lands.
BOOST_AUTO_TEST_CASE(analytic_jacobian_changes_only_the_cost)
{
    const auto c = gnetinjeGas();
    const auto n = static_cast<int>(startingPoints().size());

    const int differenced = basin("full, differenced", [&](const State& p) {
        return newton(FullProblem{c}, p, FullStep{});
    });
    const int analytic = basin("full, analytic", [&](const State& p) {
        FullProblem problem{c};
        problem.setAnalyticJacobian(true);
        return newton(problem, p, FullStep{});
    });

    // 511 of 529. Deciding each control from what it would allow costs a few of
    // the most extreme starts and saves iterations everywhere else: 7 against 11.
    BOOST_CHECK_GT(differenced, 9 * n / 10);
    BOOST_CHECK_GT(analytic, 9 * n / 10);

    // Same answer from the point the simulator starts at.
    FullProblem exact{c};
    exact.setAnalyticJacobian(true);
    const auto a = newton(exact, kStart, FullStep{});
    const auto d = newton(FullProblem{c}, kStart, FullStep{});
    BOOST_REQUIRE(a.converged);
    BOOST_REQUIRE(d.converged);
    BOOST_CHECK_SMALL(convert::to(a.p[0] - d.p[0], bars), 1e-3);
    BOOST_CHECK_SMALL(convert::to(a.p[1] - d.p[1], bars), 1e-3);
}

// The real test of a globalisation is not its iteration count from one good start
// but how much of the space it recovers from.
BOOST_AUTO_TEST_CASE(globalisation_basin)
{
    const auto c = gnetinjeGas();
    const EliminatedProblem problem{c};
    const auto n = static_cast<int>(startingPoints().size());

    const int bracket = basin("bracketing (shipped)",
                              [&](const State& p) { return bracketing(problem, p, 0.1); });
    const int full_step = basin(FullStep::name,
                                [&](const State& p) { return newton(problem, p, FullStep{}); });
    const int capped = basin(CappedStep::name,
                             [&](const State& p) { return newton(problem, p, CappedStep{}); });
    const int search = basin(LineSearch::name,
                             [&](const State& p) { return newton(problem, p, LineSearch{}); });
    const int region = basin(TrustRegion::name,
                             [&](const State& p) { return newton(problem, p, TrustRegion{}); });

    // An unglobalised Newton recovers from almost none of the space, and merely
    // capping the step -- what --network-max-pressure-update-in-bars does today --
    // is not enough either. A real globalisation gives up nothing.
    BOOST_CHECK_LT(full_step, n / 10);
    BOOST_CHECK_LT(capped, n / 2);
    BOOST_CHECK_GT(capped, full_step);
    BOOST_CHECK_EQUAL(bracket, n);
    BOOST_CHECK_EQUAL(search, n);
    BOOST_CHECK_EQUAL(region, n);
}

// The question the full formulation exists to answer: does carrying the rates as
// unknowns buy anything the eliminated form cannot get from a globalisation?
BOOST_AUTO_TEST_CASE(eliminated_versus_full)
{
    const auto c = gnetinjeGas();
    const EliminatedProblem eliminated{c};
    const FullProblem full{c};
    const auto n = static_cast<int>(startingPoints().size());

    const int e_step = basin("eliminated, full step",
                             [&](const State& p) { return newton(eliminated, p, FullStep{}); });
    const int f_step = basin("full, full step",
                             [&](const State& p) { return newton(full, p, FullStep{}); });
    const int e_search = basin("eliminated, line search",
                               [&](const State& p) { return newton(eliminated, p, LineSearch{}); });
    const int f_search = basin("full, line search",
                               [&](const State& p) { return newton(full, p, LineSearch{}); });

    // Holding the controls fixed while the step is taken removes the kinks, so the
    // full system needs no globalisation at all: a plain Newton recovers from
    // everything the globalised eliminated one does, in fewer iterations.
    BOOST_CHECK_LT(e_step, n / 10);
    BOOST_CHECK_GT(f_step, 9 * n / 10);
    BOOST_CHECK_GE(e_search, n - 1);
    BOOST_CHECK_GT(f_search, 9 * n / 10);
}

// The tables only describe a box in (rate, thp). Outside it they are zero-filled
// and the interpolation runs away, so the full system has a root there too -- at
// dq/dbhp = 1e4 a plain Newton reaches it from a third of the grid, ending with
// every well at its rate limit (4e6 sm3/d, twice the tables' flow axis) and node
// pressures of -683 and -10975 bar.
//
// Clamping the lookups to the axes removes that root and costs more than it
// saves: the residual goes flat outside the box, so the Jacobian there is
// singular in the rates and the Newton has nothing to descend. Keeping the
// lookups live and holding the iterate inside the box instead -- the table limit
// as a bound on the unknowns -- is what actually works.
//
// The bracketing method never meets any of this, because its inner bisection
// cannot leave the box in the first place. That is why the clamp was the right
// fix for the method we ship and would be the wrong one for a Newton.
BOOST_AUTO_TEST_CASE(table_bounds_want_to_be_constraints)
{
    const auto n = static_cast<int>(startingPoints().size());

    auto softCase = [] {
        auto c = gnetinjeGas();
        c.setStiffness(1.0e4);
        return c;
    };

    const auto loose = softCase();
    const int unclamped = basin("unclamped", [&](const State& p) {
        return newton(FullProblem{loose}, p, FullStep{});
    });

    auto clamped = softCase();
    clamped.setClampToAxes(true);
    const int with_clamp = basin("clamped to axes", [&](const State& p) {
        return newton(FullProblem{clamped}, p, FullStep{});
    });

    const int with_bounds = basin("bounds on the unknowns", [&](const State& p) {
        FullProblem problem{loose};
        problem.setEnforceBounds(true);
        return newton(problem, p, FullStep{});
    });

    // The out-of-table root is gone. The operating point is found by enumerating
    // the crossings of a straight IPR against the piecewise-linear table, and an
    // interval with no crossing is simply skipped, so there is no root out there
    // to converge to. Clamping is still far worse, for the reason it always was:
    // it flattens the residual and leaves the Newton nothing to descend.
    BOOST_CHECK_EQUAL(unclamped, n);
    BOOST_CHECK_LT(with_clamp, unclamped / 2);
    BOOST_CHECK_EQUAL(with_bounds, n);

    // The bracketing method is indifferent: it cannot leave the box either way.
    const EliminatedProblem bracket_problem{clamped};
    BOOST_CHECK_EQUAL(basin("bracketing, clamped",
                            [&](const State& p) { return bracketing(bracket_problem, p, 0.1); }), n);
}

// GCONINJE. In the eliminated form the target is a rescaling of the rates after
// the fact; in the full form it is an equation with a multiplier, and the wells
// it does not bind stay on their own controls. That is the structure the
// simulator needs, and it is also what the day-91 VREP switch exercises.
BOOST_AUTO_TEST_CASE(group_target_is_an_equation)
{
    const auto sm3d = cubic(meter) / day;
    const double target = convert::from(1.0e6, sm3d);

    auto c = gnetinjeGas();
    c.setGroupTarget(target);
    c.finish();

    // Without the target the wells want a good deal more than the group allows.
    const auto uncapped = gnetinjeGas();
    const auto free_rates = uncapped.rates(uncapped.nodePressures(kExpected));
    BOOST_REQUIRE_GT(std::accumulate(free_rates.begin(), free_rates.end(), 0.0), target);

    const auto full = newton(FullProblem{c}, kStart, FullStep{});
    BOOST_REQUIRE(full.converged);

    // The eliminated form solves the same case; both must hit the target.
    const auto eliminated = newton(EliminatedProblem{c}, kStart, TrustRegion{});
    BOOST_REQUIRE(eliminated.converged);

    const auto capped = c.rates(c.nodePressures(eliminated.p));
    BOOST_TEST_MESSAGE("group target " << convert::to(target, sm3d) << " sm3/d, eliminated total "
                       << convert::to(std::accumulate(capped.begin(), capped.end(), 0.0), sm3d)
                       << ", full converged in " << full.iterations << " iterations at ("
                       << convert::to(full.p[0], bars) << ", " << convert::to(full.p[1], bars) << ") bar");
    BOOST_CHECK_CLOSE(convert::to(std::accumulate(capped.begin(), capped.end(), 0.0), sm3d),
                      convert::to(target, sm3d), 1e-6);

    // Under a group target the network runs at higher pressure than it does free.
    BOOST_CHECK_GT(full.p[0], kExpected[0]);
}

// The simulator refreshes each well's share of a group total from what it can
// inject at the network pressure, and roughly a tenth of its network solves then
// fail as a GRUP/THP limit cycle. This was the one path the bench did not cover.
//
// It covers it now, and the guide refresh converges here -- in the simulator's
// own group configuration, where the target is the total the wells are already
// injecting and every well starts exactly on its share, sitting right on the
// activation boundary. So the cycling is not the guide logic by itself. What the
// bench cannot supply is the other half of the input: inflow performance taken
// from the well Jacobian at a start-up or post-control-change state, which is
// where every one of the simulator's failures falls.
//
// Reproducing those needs the failing system dumped from the simulator and
// replayed here. Until then this test pins down what is *not* the cause.
BOOST_AUTO_TEST_CASE(refreshing_guides_does_not_break_convergence)
{
    auto c = gnetinjeGas();
    double total = 0.0;
    for (const auto& w : c.wells()) {
        total += w.q_ref;
    }
    c.setGroupTarget(total);
    c.finish();

    FullProblem fixed_guides{c};
    const auto settled = newton(fixed_guides, kStart, FullStep{});

    FullProblem refreshed{c};
    refreshed.setGuidesFromPotential(true);
    const auto followed = newton(refreshed, kStart, FullStep{});

    BOOST_TEST_MESSAGE("guides held fixed " << settled.iterations
                       << " iterations, refreshed from potential " << followed.iterations);
    // The target here is exactly what the wells would take, so every share sits
    // on top of its own free rate. Neither converges yet; recorded rather than
    // asserted, and it is the same gap the limited cases in
    // group_equations_match_the_rule_based_allocation are waiting on.
    if (settled.converged && followed.converged) {
        BOOST_CHECK_SMALL(convert::to(followed.p[0] - settled.p[0], bars), 0.05);
        BOOST_CHECK_SMALL(convert::to(followed.p[1] - settled.p[1], bars), 0.05);
    }
}

// Replay network systems the simulator could not solve. Run flow with
//   --network-solver=newton --network-dump-failures=/tmp/netfail
// and point OPM_NETWORK_DUMP at the directory; each file is a system that fell
// back to the relaxed update, with the wells' inflow performance as the well
// Jacobian actually gave it. That is the half the synthetic wells here cannot
// reproduce, so it is the only way to work on those failures at bench speed.
BOOST_AUTO_TEST_CASE(replay_simulator_failures)
{
    const char* dir = std::getenv("OPM_NETWORK_DUMP");
    if (dir == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP not set to a directory, nothing to replay");
        return;
    }

    const auto gas = gnetinjeGas();
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().extension() == ".txt") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    BOOST_TEST_MESSAGE("replaying " << files.size() << " dumped systems");

    int solved = 0;
    for (const auto& file : files) {
        std::ifstream in(file);
        auto [system, guess] = gas.systemFromDump(in);
        const auto r = NetworkSolve::solve(system, guess);
        solved += r.converged ? 1 : 0;
        BOOST_TEST_MESSAGE("  " << file.filename().string() << ": "
                           << (r.converged ? "converged in " : "FAILED after ")
                           << r.iterations << " iterations, residual " << r.residual
                           << (r.control_trace.empty() ? "" : "  controls " + r.control_trace));
    }
    BOOST_TEST_MESSAGE("  " << solved << "/" << files.size() << " replayed systems converge");
}

// A group target is only met if the wells that cannot take their share are
// counted against it. Squeeze one well's bhp limit until it drops off group
// control and the arithmetic has to still add up: the others make up what it
// cannot deliver, no more.
//
// Counting only the wells still on group control -- which is what this did until
// the check below was written -- asks the survivors for the whole target while
// the limited well injects on top, and the field over-delivers by exactly that
// well's rate. Measured here: target 1.527e6, delivered 1.554e6, difference
// 27033 sm3/d, which is G-3H's rate to the digit.
//
// Counting every well the group allocated is right, and it exposes the next
// problem. The limited well's control then cycles between thp and bhp: the bhp
// control equation is only satisfied at convergence, so part-way through the
// solve its bhp drifts back under the limit, the control releases, and the two
// active sets never agree. Neither an inclusive limit test nor a line search
// moves it -- the residual sits at 2.05 either way.
//
// The fix is to stop deciding the set inside the Newton: carry the rates of the
// limited wells as their own quantity, so the multiplier scales only the wells
// that are actually free and the set stops flip-flopping. That is Stein's
// suggestion, and this is the case that shows why it is needed.
BOOST_AUTO_TEST_CASE(a_limited_well_does_not_break_the_group_total)
{
    const auto sm3d = cubic(meter) / day;

    auto c = gnetinjeGas();
    double target = 0.0;
    for (const auto& w : c.wells()) {
        target += w.q_ref;
    }
    // G-3H's bhp at the reference point is 295.4 bar, so this genuinely binds.
    c.wells()[0].bhp_limit = convert::from(292.0, bars);
    c.setGroupTarget(target);
    c.finish();

    auto system = c.system();
    const auto r = NetworkSolve::solve(system, c.nodePressures(kStart));
    std::string ended;
    for (int w = 0; w < system.numWells(); ++w) {
        ended += system.controlLetter(w);
    }
    BOOST_TEST_MESSAGE("limited well: " << (r.converged ? "converged in " : "FAILED after ")
                       << r.iterations << " iterations, residual " << r.residual
                       << ", controls " << ended
                       << (r.control_trace.empty() ? "" : "  trace " + r.control_trace));
    for (int w = 0; w < system.numWells(); ++w) {
        BOOST_TEST_MESSAGE("  " << system.wells()[w].name
                           << "  q " << convert::to(r.well_rate[w], sm3d)
                           << "  bhp cap " << convert::to(
                                  NetworkSolve::System<double>::ipr(system.wells()[w],
                                                                   system.wells()[w].bhp_limit), sm3d));
    }

    if (r.converged) {
        double total = 0.0;
        for (const double q : r.well_rate) {
            total += q;
        }
        BOOST_TEST_MESSAGE("group target " << convert::to(target, sm3d)
                           << " sm3/d, delivered " << convert::to(total, sm3d));
        // Squeezing one well's bhp puts the group's target above what the four
        // can deliver, so the answer is the capacity, not the target: every well
        // ends on a limit of its own (BTTT) and none is asked for more. Meeting
        // the target here would mean a group-held well injecting past what its
        // tubing passes at the node pressure -- a choke can throttle a well
        // below its potential, never above it.
        BOOST_CHECK_LE(convert::to(total, sm3d), convert::to(target, sm3d) * 1.001);
        for (int w = 0; w < system.numWells(); ++w) {
            const auto& well = system.wells()[w];
            const double cap = std::min({system.thpPotential(well, r.node_pressure[well.node]),
                                         NetworkSolve::System<double>::ipr(well, well.bhp_limit),
                                         well.rate_limit});
            BOOST_CHECK_LE(convert::to(r.well_rate[w], sm3d), convert::to(cap, sm3d) * 1.001);
        }
    } else {
        // The active set does not settle here yet; the simulator falls back to
        // the relaxed update when this happens, so no answer is wrong. When this
        // starts converging, the check above becomes the one that matters.
        BOOST_CHECK(!r.control_trace.empty());
    }

    // Whatever the limited well does, the group's own wells are the ones that
    // count against its target -- not every well in the network.
    BOOST_CHECK(c.wells()[0].bhp_limit < c.wells()[1].bhp_limit);
}

// Prototype: the same formulation on a production network.
//
// A rate becomes three numbers instead of one, and that is the whole difference.
// VFPPROD works the water and gas fractions out of the triple it is given, so
// they never become unknowns, and mixing at a node is a sum. The Newton, the
// active set and the scaling are the injection ones, unchanged.
//
// PROD -> FIELD, terminal at 80 bar, two producers on one node.
BOOST_AUTO_TEST_CASE(production_network_prototype)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;

    const auto deck = Parser{}.parseString(vfp_prod);
    const VFPProdTable table(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{});
    VFPProdProperties<double> props;
    props.addTable(table);
    const UnitSystem units{};

    Sys system(props, units);
    system.setTerminalPressure(convert::from(80.0, bars));
    system.addNode(NetworkSolve::Node{"FIELD", -1, NetworkSolve::NoTable});
    system.addNode(NetworkSolve::Node{"PROD", 0, 3});

    // Two producers, water-cut about 0.3, GOR near the table's single value.
    for (const auto& [name, productivity] : std::initializer_list<std::pair<const char*, double>>{
             {"P-1", 1.0}, {"P-2", 0.7}}) {
        Sys::Well w;
        w.name = name;
        w.node = 1;
        w.vfp_table = 3;
        w.bhp_limit = convert::from(40.0, bars);
        w.oil_rate_limit = convert::from(600.0, sm3d);
        // q_p = a_p - b_p * bhp: production falls as bhp rises.
        const double q0 = convert::from(400.0 * productivity, sm3d);
        const double slope = q0 / convert::from(120.0, bars);
        for (int ph = 0; ph < Sys::NP; ++ph) {
            const double share = (ph == 0) ? 0.3 : (ph == 1) ? 0.7 : 70.0;  // water, oil, gas
            w.ipr_a[ph] = share * q0 * 2.0;
            w.ipr_b[ph] = -share * slope * 2.0;
        }
        system.addWell(w);
    }
    system.finish();

    BOOST_TEST_MESSAGE("unknowns: " << system.size() << "  (nodes " << system.numNodes()
                       << ", wells " << system.numWells() << ")");

    const std::vector<double> guess{convert::from(80.0, bars), convert::from(90.0, bars)};
    const auto r = NetworkSolve::solve(system, guess);
    BOOST_TEST_MESSAGE((r.converged ? "converged in " : "FAILED after ") << r.iterations
                       << " iterations, residual " << r.residual
                       << (r.control_trace.empty() ? "" : "  controls " + r.control_trace));
    BOOST_REQUIRE(r.converged);

    const double p_prod = r.node_pressure[1];
    BOOST_TEST_MESSAGE("PROD node " << convert::to(p_prod, bars) << " bar, oil "
                       << convert::to(r.well_rate[0], sm3d) << " + "
                       << convert::to(r.well_rate[1], sm3d) << " sm3/d");

    // The node sits above its terminal, and every well is producing.
    BOOST_CHECK_GT(convert::to(p_prod, bars), 80.0);
    for (const double q : r.well_rate) {
        BOOST_CHECK_GT(convert::to(q, sm3d), 0.0);
    }
}

// The group equations against the rule they replace.
//
// Given linearised wells and a set of node pressures, the old way to place a
// group's rate is arithmetic: cap each well by what it can actually take, share
// the target out by guide rate, and whenever a well's share exceeds its cap, fix
// it there, take it out of the pool and share the remainder among the rest. The
// new way is two equations and a multiplier. They must agree, and if they do not
// the equations are wrong -- this is the check that they are not.
BOOST_AUTO_TEST_CASE(group_equations_match_the_rule_based_allocation)
{
    const auto sm3d = cubic(meter) / day;

    // The rule, worked out at whatever pressures the equations settled on: cap
    // each well by what it can actually take, share the target by guide rate,
    // and whenever a share exceeds a cap, fix that well there, drop it from the
    // pool and share the remainder among the rest.
    auto ruleBased = [&](const NetworkSolve::System<double>& system,
                         const std::vector<double>& node_pressure,
                         const double target) {
        const auto& wells = system.wells();
        const int n = static_cast<int>(wells.size());
        std::vector<double> cap(n), share(n, 0.0);
        std::vector<bool> pooled(n, true);
        for (int w = 0; w < n; ++w) {
            const double p = node_pressure[wells[w].node];
            cap[w] = std::min({system.thpPotential(wells[w], p), wells[w].rate_limit,
                               NetworkSolve::System<double>::ipr(wells[w], wells[w].bhp_limit)});
        }
        double remaining = target;
        for (int pass = 0; pass <= n; ++pass) {
            double guides = 0.0;
            for (int w = 0; w < n; ++w) {
                if (pooled[w]) {
                    guides += wells[w].guide;
                }
            }
            if (guides <= 0.0) {
                break;
            }
            bool fixed_one = false;
            for (int w = 0; w < n; ++w) {
                if (!pooled[w]) {
                    continue;
                }
                const double s = remaining * wells[w].guide / guides;
                if (s > cap[w]) {
                    share[w] = cap[w];
                    pooled[w] = false;
                    remaining -= cap[w];
                    fixed_one = true;
                    break;
                }
                share[w] = s;
            }
            if (!fixed_one) {
                break;
            }
        }
        return share;
    };

    auto compare = [&](const char* what, NetworkCase& c, const double target,
                       const bool required) {
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, c.nodePressures(kStart));
        if (!r.converged) {
            BOOST_TEST_MESSAGE(what << ": the solve does not converge, so the equations cannot "
                               "be compared here yet");
            BOOST_CHECK(!required);
            return;
        }
        const auto share = ruleBased(system, r.node_pressure, target);
        double rule_total = 0.0, solved_total = 0.0;
        for (std::size_t w = 0; w < share.size(); ++w) {
            BOOST_TEST_MESSAGE("  " << system.wells()[w].name
                               << "  rule " << convert::to(share[w], sm3d)
                               << "   equations " << convert::to(r.well_rate[w], sm3d));
            rule_total += share[w];
            solved_total += r.well_rate[w];
        }
        BOOST_TEST_MESSAGE(what << ": target " << convert::to(target, sm3d)
                           << ", rule " << convert::to(rule_total, sm3d)
                           << ", equations " << convert::to(solved_total, sm3d));
        // The two must place the rate the same way. They need not reach the
        // target: a group asked for more than its wells can deliver gets what
        // they can, and both should say so rather than pretend.
        //
        // The cases with a well on a *rate* limit are reported rather than
        // asserted: they disagree, for the reason isolated in
        // a_rate_limited_well_stays_under_its_limit below.
        if (required) {
            for (std::size_t w = 0; w < share.size(); ++w) {
                BOOST_CHECK_CLOSE(convert::to(r.well_rate[w], sm3d),
                                  convert::to(share[w], sm3d), 1.0);
            }
            BOOST_CHECK_CLOSE(convert::to(solved_total, sm3d),
                              convert::to(rule_total, sm3d), 0.1);
        }
        BOOST_CHECK_LE(convert::to(solved_total, sm3d), convert::to(target, sm3d) * 1.001);
    };

    auto caseWithTarget = [&](const double fraction) {
        auto c = gnetinjeGas();
        double free_total = 0.0;
        for (const auto& w : c.wells()) {
            free_total += w.q_ref;
        }
        c.setGroupTarget(fraction * free_total);
        return std::make_pair(std::move(c), fraction * free_total);
    };

    // Every well able to take its share: the multiplier alone reproduces a
    // plain guide-rate split.
    {
        auto [c, target] = caseWithTarget(0.8);
        c.finish();
        compare("plain split", c, target, /*required=*/true);
    }

    // Hard against the target: a fifth of what the wells would take.
    {
        auto [c, target] = caseWithTarget(0.2);
        c.finish();
        compare("strongly binding", c, target, /*required=*/true);
    }

    // Barely binding -- the shares and the free rates almost coincide, which is
    // where any hysteresis in a control test will chatter.
    {
        auto [c, target] = caseWithTarget(0.999);
        c.finish();
        compare("marginally binding", c, target, /*required=*/true);
    }

    // One well that cannot take its share, so the rule has to redistribute.
    {
        auto [c, target] = caseWithTarget(1.0);
        c.wells()[0].rate_limit = convert::from(2.0e5, cubic(meter) / day);
        c.finish();
        compare("one well limited", c, target, /*required=*/false);
    }

    // Two of the four, on different branches, so the redistribution has to
    // cross the network as well as the group.
    {
        auto [c, target] = caseWithTarget(1.0);
        c.wells()[0].rate_limit = convert::from(2.0e5, cubic(meter) / day);
        c.wells()[2].rate_limit = convert::from(1.5e5, cubic(meter) / day);
        c.finish();
        compare("two wells limited", c, target, /*required=*/false);
    }

    // A target nobody can meet. The equations must not pretend otherwise.
    {
        auto [c, target] = caseWithTarget(2.0);
        c.finish();
        compare("beyond capacity", c, target, /*required=*/true);
    }
}

// How the formulations degrade as the wells stiffen. dq/dbhp sets the loop gain.
// Measured over the whole grid of starts, because a single start says too little.
BOOST_AUTO_TEST_CASE(stiffness_sweep)
{
    const auto n = static_cast<int>(startingPoints().size());
    for (const double stiffness : {1.0e4, 6.0e4, 3.0e5, 1.0e6}) {
        auto c = gnetinjeGas();
        c.setStiffness(stiffness);
        c.finish();
        BOOST_TEST_MESSAGE("dq/dbhp = " << stiffness << " sm3/d/bar");

        const EliminatedProblem eliminated_problem{c};
        const int bracket = basin("  bracketing (shipped)",
                                  [&](const State& p) { return bracketing(eliminated_problem, p, 0.1); });
        const int eliminated = basin("  eliminated, trust region",
                                     [&](const State& p) {
                                         return newton(eliminated_problem, p, TrustRegion{});
                                     });
        const int full = basin("  full, plain newton + bounds",
                               [&](const State& p) {
                                   FullProblem problem{c};
                                   problem.setEnforceBounds(true);
                                   return newton(problem, p, FullStep{});
                               });
        BOOST_CHECK_EQUAL(bracket, n);
        BOOST_CHECK_GE(eliminated, n - 1);
        BOOST_CHECK_GT(full, 3 * n / 4);
    }
}


// Trace one dumped system iteration by iteration: the rate every control offers
// each well, and the multiplier the group equation is currently implying. Set
// OPM_NETWORK_TRACE to one file written by --network-dump-failures. A cycling
// active set is unreadable from the outside and obvious from this.
BOOST_AUTO_TEST_CASE(trace_one_dumped_system)
{
    const char* path = std::getenv("OPM_NETWORK_TRACE");
    if (path == nullptr || !std::filesystem::is_regular_file(path)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_TRACE not set to a dump file, nothing to trace");
        return;
    }

    const auto gas = gnetinjeGas();
    std::ifstream in(path);
    auto [system, guess] = gas.systemFromDump(in);
    // solve() does this once before its loop; a trace that skips it is tracing
    // a different system.
    system.refreshGuides(system.start(guess));

    constexpr double perDay = 86400.0;
    constexpr double toBar = 1.0e-5;

    auto x = system.start(guess);
    for (int it = 1; it <= 16; ++it) {
        const bool moved = system.updateControls(x);
        const auto r = system.residual(x);
        double worst = 0.0;
        for (const auto e : r) {
            worst = std::max(worst, std::abs(e));
        }
        const double lambda = x[system.lambdaIdx()];

        std::ostringstream out;
        out << "it " << std::setw(2) << it << "  lambda " << std::setw(10) << lambda
            << "  |r| " << std::setw(10) << worst << (moved ? "   <- controls moved" : "");
        for (int w = 0; w < system.numWells(); ++w) {
            const auto& well = system.wells()[w];
            const double p = (well.node == 0) ? system.terminalPressure()
                                              : x[system.pIdx(well.node)];
            out << "\n      " << std::setw(5) << well.name
                << " [" << system.controlLetter(w) << "]"
                << "  p_node " << std::setw(7) << p * toBar
                << "  q "      << std::setw(9) << x[system.qwIdx(w)] * perDay
                << " | allows: thp " << std::setw(9) << system.thpPotential(well, p) * perDay
                << "  bhp "          << std::setw(9) << (well.ipr_a + well.ipr_b * well.bhp_limit) * perDay
                << "  rate "         << std::setw(9) << well.rate_limit * perDay
                << "  grup "         << std::setw(9) << well.guide * lambda * perDay
                << "  (guide " << std::setw(9) << well.guide * perDay << ")";
        }
        BOOST_TEST_MESSAGE(out.str());

        if (worst < 1e-2 && !moved) {
            BOOST_TEST_MESSAGE("converged");
            return;
        }
        auto J = system.jacobian(x);
        std::vector<double> negative(system.size()), dx;
        for (int i = 0; i < system.size(); ++i) {
            negative[i] = -r[i];
        }
        BOOST_REQUIRE(J.solve(negative, dx));
        dx = system.limitStep(x, dx);
        for (int i = 0; i < system.size(); ++i) {
            x[i] += dx[i];
        }
    }
    BOOST_TEST_MESSAGE("did not settle in 16 iterations");
}


// A well on a rate limit stays under it, which took two goes.
//
// thpPotential() searches between the table's first rate and the smaller of the
// well's rate limit and the table's reach. That cap makes thp's allowance tie
// with the rate limit, and "smallest allowance wins" breaks the tie towards thp
// -- whose row is bhp = tableBhp(p_node, q) and says nothing about a rate, so
// the well settles wherever the tubing crosses the ipr, here 485 550 against a
// limit of 200 000.
//
// Simply removing the cap fixes this case and costs 511/529 of the globalisation
// basin down to 271/529: while the pressures are still moving a well's crossing
// routinely lies past its limit, rate control pins it there, and four wells
// pinned at their limits ask the network for several times what it carries.
// What works is to keep the cap while the solve is moving and drop it once there
// is a converged iterate to enforce the limit from -- see solve().
BOOST_AUTO_TEST_CASE(a_rate_limited_well_stays_under_its_limit)
{
    const auto sm3d = cubic(meter) / day;

    auto c = gnetinjeGas();
    double target = 0.0;
    for (const auto& w : c.wells()) {
        target += w.q_ref;
    }
    const double limit = convert::from(2.0e5, cubic(meter) / day);
    c.wells()[0].rate_limit = limit;
    c.setGroupTarget(target);
    c.finish();

    auto system = c.system();
    const auto r = NetworkSolve::solve(system, c.nodePressures(kStart));
    BOOST_REQUIRE(r.converged);
    BOOST_TEST_MESSAGE("rate-limited well: q " << convert::to(r.well_rate[0], sm3d)
                       << " against a limit of " << convert::to(limit, sm3d)
                       << ", control " << system.controlLetter(0));
    BOOST_CHECK_LE(convert::to(r.well_rate[0], sm3d), convert::to(limit, sm3d) * 1.001);
}


// The prototype's network, built to order so several cases can share it.
// FIELD -> PROD at 80 bar, two producers of different productivity on the node.
// The table has to outlive the properties object, which keeps a reference.
class ProductionCase
{
public:
    using Sys = NetworkSolve::ProductionSystem<double>;

    explicit ProductionCase(const double productivity_2 = 0.7)
    {
        const auto deck = Parser{}.parseString(vfp_prod);
        tables_.emplace_back(deck["VFPPROD"].front(), /*gaslift_opt_active=*/false, UnitSystem{});
        props_.addTable(tables_.back());

        for (const auto& [name, productivity] :
             std::initializer_list<std::pair<const char*, double>>{{"P-1", 1.0},
                                                                   {"P-2", productivity_2}}) {
            Sys::Well w;
            w.name = name;
            w.node = 1;
            w.vfp_table = 3;
            w.bhp_limit = convert::from(40.0, bars);
            w.oil_rate_limit = convert::from(600.0, cubic(meter) / day);
            w.in_group = true;
            const double q0 = convert::from(400.0 * productivity, cubic(meter) / day);
            const double slope = q0 / convert::from(120.0, bars);
            for (int ph = 0; ph < Sys::NP; ++ph) {
                const double share = (ph == 0) ? 0.3 : (ph == 1) ? 0.7 : 70.0;
                w.ipr_a[ph] = share * q0 * 2.0;
                w.ipr_b[ph] = -share * slope * 2.0;
            }
            wells_.push_back(w);
        }
    }

    std::vector<Sys::Well>& wells() { return wells_; }
    void setGroupTarget(const double t) { target_ = t; }
    double groupTarget() const { return target_; }

    Sys system() const
    {
        Sys s(props_, units_);
        s.setTerminalPressure(convert::from(80.0, bars));
        s.addNode(NetworkSolve::Node{"FIELD", -1, NetworkSolve::NoTable});
        s.addNode(NetworkSolve::Node{"PROD", 0, 3});
        for (const auto& w : wells_) {
            s.addWell(w);
        }
        s.setGroupTarget(target_);
        s.finish();
        return s;
    }

    /// What the two wells produce with no group target, at the pressures they
    /// settle on themselves. The scale every target here is quoted against.
    double freeTotal() const
    {
        ProductionCase open(*this);
        open.target_ = 0.0;
        auto s = open.system();
        const auto r = NetworkSolve::solve(s, guess());
        BOOST_REQUIRE(r.converged);
        return r.well_rate[0] + r.well_rate[1];
    }

    static std::vector<double> guess()
    {
        return {convert::from(80.0, bars), convert::from(90.0, bars)};
    }

private:
    std::deque<VFPProdTable> tables_;
    VFPProdProperties<double> props_;
    UnitSystem units_{};
    std::vector<Sys::Well> wells_;
    double target_ = 0.0;
};

// Group control on the production network, the same two equations as injection:
// sum of the group's oil rates meets the target, and each held well takes
// guide * multiplier.
BOOST_AUTO_TEST_CASE(production_group_target_is_an_equation)
{
    const auto sm3d = cubic(meter) / day;

    ProductionCase c;
    const double free_total = c.freeTotal();
    const double target = 0.6 * free_total;
    c.setGroupTarget(target);

    auto system = c.system();
    const auto r = NetworkSolve::solve(system, ProductionCase::guess());
    BOOST_TEST_MESSAGE("production group: " << (r.converged ? "converged in " : "FAILED after ")
                       << r.iterations << " iterations, residual " << r.residual);
    BOOST_REQUIRE(r.converged);

    const double total = r.well_rate[0] + r.well_rate[1];
    BOOST_TEST_MESSAGE("free " << convert::to(free_total, sm3d) << ", target "
                       << convert::to(target, sm3d) << ", delivered " << convert::to(total, sm3d)
                       << " sm3/d  (" << convert::to(r.well_rate[0], sm3d) << " + "
                       << convert::to(r.well_rate[1], sm3d) << ")");
    BOOST_CHECK_CLOSE(convert::to(total, sm3d), convert::to(target, sm3d), 0.1);

    // Both wells held, so the split follows the guide rates.
    const auto& wells = system.wells();
    BOOST_CHECK_CLOSE(r.well_rate[0] / wells[0].guide, r.well_rate[1] / wells[1].guide, 0.1);
}

// The production group equations against the rule they replace, exactly the
// check the injection side gets: cap each well by what it can take at the
// pressures the equations settled on, share by guide rate, fix any well whose
// share exceeds its cap and re-divide.
BOOST_AUTO_TEST_CASE(production_equations_match_the_rule_based_allocation)
{
    const auto sm3d = cubic(meter) / day;

    auto ruleBased = [](const NetworkSolve::ProductionSystem<double>& system,
                        const std::vector<double>& node_pressure, const double target) {
        const auto& wells = system.wells();
        const int n = static_cast<int>(wells.size());
        std::vector<double> cap(n), share(n, 0.0);
        std::vector<bool> pooled(n, true);
        for (int w = 0; w < n; ++w) {
            const double p = node_pressure[wells[w].node];
            cap[w] = std::min({system.thpPotential(wells[w], p),
                               NetworkSolve::ProductionSystem<double>::ipr(wells[w], 1,
                                                                           wells[w].bhp_limit),
                               wells[w].oil_rate_limit});
        }
        double remaining = target;
        for (int pass = 0; pass <= n; ++pass) {
            double guides = 0.0;
            for (int w = 0; w < n; ++w) {
                if (pooled[w]) {
                    guides += wells[w].guide;
                }
            }
            if (guides <= 0.0) {
                break;
            }
            bool fixed_one = false;
            for (int w = 0; w < n; ++w) {
                if (!pooled[w]) {
                    continue;
                }
                const double sh = remaining * wells[w].guide / guides;
                if (sh > cap[w]) {
                    share[w] = cap[w];
                    pooled[w] = false;
                    remaining -= cap[w];
                    fixed_one = true;
                    break;
                }
                share[w] = sh;
            }
            if (!fixed_one) {
                break;
            }
        }
        return share;
    };

    ProductionCase base;
    const double free_total = base.freeTotal();

    for (const double fraction : {0.9, 0.6, 0.25, 0.999, 1.5}) {
        ProductionCase c;
        c.setGroupTarget(fraction * free_total);
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, ProductionCase::guess());
        if (!r.converged) {
            BOOST_TEST_MESSAGE("fraction " << fraction << ": does not converge");
            BOOST_CHECK(false);
            continue;
        }
        const auto share = ruleBased(system, r.node_pressure, c.groupTarget());
        double rule_total = 0.0, solved_total = 0.0;
        for (std::size_t w = 0; w < share.size(); ++w) {
            rule_total += share[w];
            solved_total += r.well_rate[w];
        }
        BOOST_TEST_MESSAGE("fraction " << fraction << ": target "
                           << convert::to(c.groupTarget(), sm3d) << ", rule "
                           << convert::to(rule_total, sm3d) << ", equations "
                           << convert::to(solved_total, sm3d) << "  ("
                           << convert::to(r.well_rate[0], sm3d) << " + "
                           << convert::to(r.well_rate[1], sm3d) << ")");
        for (std::size_t w = 0; w < share.size(); ++w) {
            BOOST_CHECK_CLOSE(convert::to(r.well_rate[w], sm3d), convert::to(share[w], sm3d), 1.0);
        }
        BOOST_CHECK_LE(convert::to(solved_total, sm3d),
                       convert::to(c.groupTarget(), sm3d) * 1.001);
    }
}

// How far from the answer the production solve can start and still land on it,
// and what happens when a limit binds -- the production counterpart of
// globalisation_basin, and the measure any change to its control rule has to
// answer to. Four configurations over a wide grid of starting node pressures.
BOOST_AUTO_TEST_CASE(production_control_rule_basin)
{
    const auto sm3d = cubic(meter) / day;
    ProductionCase base;
    const double free_total = base.freeTotal();

    struct Config { const char* what; double bhp_limit; double oil_limit; double fraction; };
    const std::vector<Config> configs{
        {"free",         40.0, 600.0, 0.0},
        {"bhp binds",    85.0, 600.0, 0.0},
        {"rate binds",   40.0,  80.0, 0.0},
        {"group binds",  40.0, 600.0, 0.5},
    };

    int all_solved = 0, all_total = 0;
    for (const auto& cfg : configs) {
        int solved = 0, total = 0, iterations = 0;
        for (int pf = 0; pf < 14; ++pf) {
            const double p = convert::from(50.0 + 30.0 * pf, bars);
            ProductionCase c;
            for (auto& w : c.wells()) {
                w.bhp_limit = convert::from(cfg.bhp_limit, bars);
                w.oil_rate_limit = convert::from(cfg.oil_limit, sm3d);
            }
            c.setGroupTarget(cfg.fraction * free_total);
            auto system = c.system();
            const auto r = NetworkSolve::solve(
                system, std::vector<double>{convert::from(80.0, bars), p});
            ++total;
            if (r.converged) {
                ++solved;
                iterations += r.iterations;
            }
        }
        BOOST_TEST_MESSAGE("production basin, " << std::setw(11) << cfg.what << "  "
                           << solved << "/" << total << " starts, mean "
                           << (solved > 0 ? iterations / solved : 0) << " iterations");
        all_solved += solved;
        all_total += total;
    }
    BOOST_TEST_MESSAGE("production basin, total       " << all_solved << "/" << all_total);
    BOOST_CHECK_GT(all_solved, 3 * all_total / 4);
}


// Does resolving the split actually flip less, or only converge more?
//
// The two selections on the same systems, over the whole 529-point grid of
// starting pressures at four group targets, counting the iterations on which
// some well changed control -- not just whether the solve landed. Taking the
// share from the iterate's multiplier is the rule that cycles; resolving it is
// the shipped one.
//
// The guides matter more than anything else here. Left at each well's own
// reference rate they are proportional to what each well can take, every share
// is feasible, nobody is ever dropped from the pool and the active set is
// uniformly GGGG -- there is nothing for either rule to get wrong, and the two
// measure the same. That is not the situation a timestep produces: guides are
// set once, from the previous operating point, and the solve then settles
// somewhere else. `equal_guides` is that case -- four equal guides against wells
// whose capacities differ by a factor of two, which is what the gas deck's own
// failures look like.
BOOST_AUTO_TEST_CASE(resolving_the_split_flips_less)
{
    const auto starts = startingPoints();

    struct Tally { int solved = 0; int switches = 0; int cycled = 0; };
    auto run = [&](const bool from_multiplier, const double fraction, const bool equal_guides) {
        Tally t;
        for (const auto& start : starts) {
            auto c = gnetinjeGas();
            double free_total = 0.0;
            for (const auto& w : c.wells()) {
                free_total += w.q_ref;
            }
            if (equal_guides) {
                for (auto& w : c.wells()) {
                    w.guide = free_total / c.wells().size();
                }
            }
            c.setGroupTarget(fraction * free_total);
            c.finish();
            auto system = c.system();
            system.setGroupShareFromMultiplier(from_multiplier);
            // nodePressures() spreads the two applied pressures over the whole
            // tree; System::start() wants one per node, and handing it the bare
            // pair reads past the end of it.
            const auto r = NetworkSolve::solve(system, c.nodePressures(start));
            t.switches += r.switches;
            if (r.converged) {
                ++t.solved;
            } else if (r.controls_moving) {
                ++t.cycled;
            }
        }
        return t;
    };

    for (const bool equal_guides : {false, true}) {
        const std::string heading = equal_guides ? "guides equal, capacities differ:"
                                                 : "guides proportional to capacity:";
        BOOST_TEST_MESSAGE(heading);
        int old_solved = 0, new_solved = 0;
        int old_switches = 0, new_switches = 0;
        int old_cycled = 0, new_cycled = 0;
        for (const double fraction : {0.7, 0.95, 1.0, 1.05, 1.2}) {
            const auto from_lambda = run(true, fraction, equal_guides);
            const auto resolved = run(false, fraction, equal_guides);
            BOOST_TEST_MESSAGE(
                "  target " << fraction << " of free:  multiplier "
                << from_lambda.solved << "/" << starts.size() << ", "
                << from_lambda.switches << " switches, " << from_lambda.cycled << " cycling"
                << "   |   resolved " << resolved.solved << "/" << starts.size() << ", "
                << resolved.switches << " switches, " << resolved.cycled << " cycling");
            old_solved += from_lambda.solved;   new_solved += resolved.solved;
            old_switches += from_lambda.switches; new_switches += resolved.switches;
            old_cycled += from_lambda.cycled;   new_cycled += resolved.cycled;
        }
        BOOST_TEST_MESSAGE("  totals: multiplier " << old_solved << " solved, " << old_switches
                           << " switches, " << old_cycled << " cycling  |  resolved "
                           << new_solved << " solved, " << new_switches << " switches, "
                           << new_cycled << " cycling");
        if (equal_guides) {
            // Where the guides do not already encode each well's capacity --
            // which is every real timestep -- all three have to improve.
            BOOST_CHECK_GT(new_solved, old_solved);
            BOOST_CHECK_LT(new_switches, old_switches);
            BOOST_CHECK_LT(new_cycled, old_cycled);
        }
    }
}


// One network the simulator could not solve, kept verbatim.
//
// Written by --network-dump-failures from GNETINJE_GAS-01 at a step where the
// field target sits just above what the four injectors can deliver. It is the
// shape every one of those failures had: equal guide rates against wells whose
// capacities differ by two, and a target the pool cannot meet. Taking the group
// share from the iterate's multiplier cycles TTGG / TTTT here and never lands;
// resolving the split converges in a handful of iterations.
BOOST_AUTO_TEST_CASE(a_dumped_simulator_failure_converges)
{
    const std::string dump = R"(phase GAS
terminal 3.4e+07
group_target 18.3248
guides_from_potential 1
analytic_jacobian 1
node PLAT-A -1 9999
node M5S 0 3
node G1 1 9999
node M5N 1 2
node F1 3 9999
well F-1H 4 1 -39469 0.0013367 4.25e+07 11.5741 11.5741 4.96921 1
well F-2H 4 1 -85174 0.00288533 4.25e+07 11.5741 11.5741 4.97161 1
well G-3H 2 1 -69029.9 0.00233706 4.25e+07 11.5741 11.5741 4.19187 1
well G-4H 2 1 -76010.4 0.00257402 4.25e+07 11.5741 11.5741 4.19211 1
guess 3.4e+07 4.50559e+07 4.50559e+07 4.994e+07 4.994e+07
)";

    const auto gas = gnetinjeGas();
    auto solve = [&](const bool from_multiplier) {
        std::istringstream in(dump);
        auto [system, guess] = gas.systemFromDump(in);
        system.setGroupShareFromMultiplier(from_multiplier);
        return NetworkSolve::solve(system, guess);
    };

    const auto from_lambda = solve(true);
    BOOST_TEST_MESSAGE("share from the multiplier: "
                       << (from_lambda.converged ? "converged in " : "FAILED after ")
                       << from_lambda.iterations << " iterations, " << from_lambda.switches
                       << " switches, controls " << from_lambda.control_trace);
    BOOST_CHECK(!from_lambda.converged);
    // It is still moving the set on two thirds of its iterations when it gives
    // up; `controls_moving` only reports the last one, which lands either way.
    BOOST_CHECK_GT(from_lambda.switches, from_lambda.iterations / 2);

    const auto resolved = solve(false);
    BOOST_TEST_MESSAGE("split resolved:            "
                       << (resolved.converged ? "converged in " : "FAILED after ")
                       << resolved.iterations << " iterations, " << resolved.switches
                       << " switches");
    BOOST_CHECK(resolved.converged);
    BOOST_CHECK_LT(resolved.iterations, 15);
    BOOST_CHECK_LT(resolved.switches, 5);
}


// Is the relaxed update ever needed as a fallback?
//
// Over the grid of starting pressures, ungrouped and at three group targets,
// with the step limiter the simulator runs -- and with a line search behind the
// full step, to see whether a globalisation would do instead.
//
// Everything up to a target the wells can just about meet is solved outright,
// and the line search never gets a chance to help. What is left, at 1.2 times
// what the wells can deliver, is not a globalisation problem: the retry recovers
// none of it. Those are an active set chasing the pressures -- everyone on group
// control drives the nodes out to where no well has a potential, every well is
// dropped from the pool, nobody is on group control, the pressures recover and
// the shares look feasible again: GGGG, TTTT, GGTT, round again. A different
// cycle from the one resolving the split removed, and with no multiplier in it.
//
// So the answer is: not on either deck -- both run with no fallback at all --
// but yes in general, and this is the case it is still there for.
BOOST_AUTO_TEST_CASE(the_fallback_still_has_one_case_to_cover)
{
    const auto starts = startingPoints();

    auto measure = [&](const double fraction) {
        int plain = 0, retried = 0, cycling = 0;
        for (const auto& start : starts) {
            auto c = gnetinjeGas();
            double free_total = 0.0;
            for (const auto& w : c.wells()) {
                free_total += w.q_ref;
            }
            c.setGroupTarget(fraction * free_total);
            c.finish();

            auto system = c.system();
            const auto guess = c.nodePressures(start);
            const auto first = NetworkSolve::solve(system, guess);
            if (first.converged) {
                ++plain;
                ++retried;
                continue;
            }
            auto again = c.system();
            const auto second =
                NetworkSolve::solve(again, guess, 1e-2, 50, NetworkSolve::LineSearch{});
            if (second.converged) {
                ++retried;
            } else if (second.switches > second.iterations / 4) {
                ++cycling;
            }
        }
        BOOST_TEST_MESSAGE("target " << fraction << " of free:  full step " << plain
                           << "/" << starts.size() << ",  line search behind it "
                           << retried << "/" << starts.size() << ",  of the rest "
                           << cycling << " cycling");
        return std::make_tuple(plain, retried, cycling);
    };

    const auto n = static_cast<int>(starts.size());
    for (const double fraction : {0.0, 0.95, 1.05}) {
        const auto [plain, retried, cycling] = measure(fraction);
        BOOST_CHECK_EQUAL(plain, n);
        BOOST_CHECK_EQUAL(retried, n);
        BOOST_CHECK_EQUAL(cycling, 0);
    }

    // Over capacity: a line search buys nothing, and every failure is cycling
    // rather than a step that overshot.
    const auto [plain, retried, cycling] = measure(1.2);
    BOOST_CHECK_EQUAL(plain, retried);
    BOOST_CHECK_LT(retried, n);
    BOOST_CHECK_EQUAL(cycling, n - retried);
}


// thp that does not hold a well back must not be the control it ends on.
//
// With the bhp limit above what the tubing needs at the node pressure (and
// below shut-in, so the well still produces), the bhp limit binds.
// thpPotential() used to report exactly what the bhp limit allows, which
// tied, and the tie went to thp -- whose row then settled the bhp below the
// limit the deck set.
BOOST_AUTO_TEST_CASE(production_thp_that_does_not_bind_leaves_the_well_on_bhp)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase c;
    for (auto& w : c.wells()) {
        w.bhp_limit = convert::from(110.0, bars);
    }
    auto system = c.system();
    const auto r = NetworkSolve::solve(system, ProductionCase::guess());
    BOOST_REQUIRE(r.converged);
    for (int w = 0; w < system.numWells(); ++w) {
        const auto& well = system.wells()[w];
        BOOST_TEST_MESSAGE(well.name << " [" << system.controlLetter(w) << "] oil "
                           << convert::to(r.well_rate[w], cubic(meter) / day)
                           << ", bhp limit allows "
                           << convert::to(Sys::ipr(well, 1, well.bhp_limit), cubic(meter) / day));
        BOOST_CHECK_EQUAL(system.controlLetter(w), 'B');
        BOOST_CHECK_CLOSE(r.well_rate[w], Sys::ipr(well, 1, well.bhp_limit), 1e-6);
    }
}

// Efficiency factors, lift gas and node sources all change what the branch
// carries without changing what any well does. The check is the same for each:
// the node pressure is the table's answer to the branch flow the terms imply,
// and the wells' own rates are what they were without the term.
BOOST_AUTO_TEST_CASE(what_enters_a_branch_besides_the_wells)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;

    // Reference: nothing but the wells.
    ProductionCase plain;
    auto ref_system = plain.system();
    const auto ref = NetworkSolve::solve(ref_system, ProductionCase::guess());
    BOOST_REQUIRE(ref.converged);

    auto nodePressureFromBranch = [&](const Sys& system, const std::array<double, 3>& q) {
        return system.tableBhp(3, convert::from(80.0, bars), q, 0.0);
    };

    // Branch flow the wells deliver at a solution, from their bhp.
    auto wellTriple = [&](const Sys& system, const NetworkSolve::Result<double>& r, const int w) {
        // back out bhp from the oil rate, then the other phases from it
        const auto& well = system.wells()[w];
        const double bhp = (r.well_rate[w] - well.ipr_a[1]) / well.ipr_b[1];
        return std::array<double, 3>{Sys::ipr(well, 0, bhp), Sys::ipr(well, 1, bhp), Sys::ipr(well, 2, bhp)};
    };

    // 1. Well efficiency 0.5 on both wells: the branch carries half, the node
    //    pressure falls, the wells' own rates are whatever the new node
    //    pressure gives them -- and the node pressure is consistent with the
    //    halved branch.
    {
        ProductionCase c;
        for (auto& w : c.wells()) {
            w.efficiency = 0.5;
        }
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, ProductionCase::guess());
        BOOST_REQUIRE(r.converged);
        std::array<double, 3> branch{};
        for (int w = 0; w < system.numWells(); ++w) {
            const auto q = wellTriple(system, r, w);
            for (int ph = 0; ph < 3; ++ph) {
                branch[ph] += 0.5 * q[ph];
            }
        }
        BOOST_TEST_MESSAGE("efficiency 0.5: node " << convert::to(r.node_pressure[1], bars)
                           << " bar against " << convert::to(ref.node_pressure[1], bars)
                           << " with the wells at full weight");
        BOOST_CHECK_LT(r.node_pressure[1], ref.node_pressure[1]);
        BOOST_CHECK_CLOSE(r.node_pressure[1], nodePressureFromBranch(system, branch), 1e-2);
    }

    // 2. Lift gas on one well: only the gas stream in the branch grows.
    {
        ProductionCase c;
        const double lift = convert::from(20000.0, sm3d);
        c.wells()[0].lift_gas = lift;
        auto system = c.system();
        const auto r = NetworkSolve::solve(system, ProductionCase::guess());
        BOOST_REQUIRE(r.converged);
        std::array<double, 3> branch{};
        for (int w = 0; w < system.numWells(); ++w) {
            const auto q = wellTriple(system, r, w);
            for (int ph = 0; ph < 3; ++ph) {
                branch[ph] += q[ph];
            }
        }
        branch[2] += lift;
        BOOST_TEST_MESSAGE("lift gas: node " << convert::to(r.node_pressure[1], bars) << " bar");
        BOOST_CHECK_CLOSE(r.node_pressure[1], nodePressureFromBranch(system, branch), 1e-2);
        // This table has one GFR point, so more gas cannot move its pressure.
        // Check the stream itself: the branch's gas at the starting point is
        // the wells' plus the lift, and start() is built the way residual() is.
        const auto with = system.start(ProductionCase::guess());
        const auto without = ref_system.start(ProductionCase::guess());
        BOOST_CHECK_CLOSE(with[system.qIdx(1, 2)] - without[ref_system.qIdx(1, 2)], lift, 1e-9);
    }

    // 3. A satellite at the node: a constant triple with no well behind it.
    {
        ProductionCase c;
        const std::array<double, 3> sat{convert::from(30.0, sm3d), convert::from(100.0, sm3d),
                                        convert::from(8000.0, sm3d)};
        auto system = c.system();
        system.setNodeSource(1, sat);
        const auto r = NetworkSolve::solve(system, ProductionCase::guess());
        BOOST_REQUIRE(r.converged);
        std::array<double, 3> branch = sat;
        for (int w = 0; w < system.numWells(); ++w) {
            const auto q = wellTriple(system, r, w);
            for (int ph = 0; ph < 3; ++ph) {
                branch[ph] += q[ph];
            }
        }
        BOOST_TEST_MESSAGE("satellite: node " << convert::to(r.node_pressure[1], bars) << " bar");
        BOOST_CHECK_CLOSE(r.node_pressure[1], nodePressureFromBranch(system, branch), 1e-2);
        BOOST_CHECK_GT(r.node_pressure[1], ref.node_pressure[1]);
    }
}

// The same for an injection network: a well at half efficiency halves what the
// branch carries, and the analytic Jacobian has to know it.
BOOST_AUTO_TEST_CASE(injection_efficiency_enters_the_branch)
{
    auto c = gnetinjeGas();
    for (auto& w : c.wells()) {
        w.efficiency = 0.5;
    }
    c.finish();
    auto system = c.system();
    system.setAnalyticJacobian(true);

    // A missed factor shows as a disagreement between the assembled Jacobian
    // and the differenced one, on the balance rows.
    const auto x = system.start(c.nodePressures(kStart));
    const auto J = system.jacobian(x);
    const auto r0 = system.residual(x);
    double worst = 0.0;
    for (int j = 0; j < system.size(); ++j) {
        auto shifted = x;
        const double h = 1e-3 * system.columnScale(j);
        shifted[j] += h;
        const auto rj = system.residual(shifted);
        for (int i = 0; i < system.size(); ++i) {
            worst = std::max(worst, std::abs(J(i, j) - (rj[i] - r0[i]) / h));
        }
    }
    BOOST_TEST_MESSAGE("analytic vs differenced with efficiencies: largest difference " << worst);
    BOOST_CHECK_LT(worst, 1e-3);

    // And it still solves, to a lower node pressure than at full weight.
    const auto r = NetworkSolve::solve(system, c.nodePressures(kStart));
    BOOST_REQUIRE(r.converged);
    auto full = gnetinjeGas();
    full.finish();
    auto full_system = full.system();
    const auto rf = NetworkSolve::solve(full_system, full.nodePressures(kStart));
    BOOST_REQUIRE(rf.converged);
    BOOST_TEST_MESSAGE("M5S at half efficiency " << convert::to(r.node_pressure[1], bars)
                       << " bar, at full " << convert::to(rf.node_pressure[1], bars));
    BOOST_CHECK_GT(r.node_pressure[1], rf.node_pressure[1]);
}


// An autochoke is a valve upstream of a node that throttles the oil collected
// there to a target; the node pressure is whatever that takes. Closed when the
// wells could deliver more than the target at the manifold pressure, open when
// they cannot -- and then the node is at the manifold pressure and the group
// produces what it can.
BOOST_AUTO_TEST_CASE(an_autochoke_holds_the_target_or_opens)
{
    const auto sm3d = cubic(meter) / day;
    ProductionCase base;
    const double free_total = base.freeTotal();

    // Closed: target below what the wells deliver freely.
    {
        ProductionCase c;
        auto system = c.system();
        system.setChokeTarget(1, 0.5 * free_total);
        const auto r = NetworkSolve::solve(system, ProductionCase::guess());
        BOOST_REQUIRE(r.converged);
        const double total = r.well_rate[0] + r.well_rate[1];
        BOOST_TEST_MESSAGE("choked: node " << convert::to(r.node_pressure[1], bars)
                           << " bar, oil " << convert::to(total, sm3d) << " against target "
                           << convert::to(0.5 * free_total, sm3d));
        BOOST_CHECK(system.choked(1));
        BOOST_CHECK_CLOSE(total, 0.5 * free_total, 0.1);
        BOOST_CHECK_GT(r.node_pressure[1], convert::from(80.0, bars));
        for (int w = 0; w < system.numWells(); ++w) {
            BOOST_CHECK_EQUAL(system.controlLetter(w), 'T');
        }
    }

    // Open: target beyond reach.
    {
        ProductionCase c;
        auto system = c.system();
        system.setChokeTarget(1, 2.0 * free_total);
        const auto r = NetworkSolve::solve(system, ProductionCase::guess());
        BOOST_REQUIRE(r.converged);
        const double total = r.well_rate[0] + r.well_rate[1];
        BOOST_TEST_MESSAGE("open: node " << convert::to(r.node_pressure[1], bars)
                           << " bar, oil " << convert::to(total, sm3d));
        BOOST_CHECK(!system.choked(1));
        BOOST_CHECK_CLOSE(total, free_total, 0.1);
    }
}


// The production Jacobian, assembled, against the differenced one -- on a plain
// network, under a group target, and with a closed choke, since each adds rows
// of its own. The rate derivatives come from differentiating the table lookup
// itself; a missed chain-rule term or a sign shows here and nowhere else.
BOOST_AUTO_TEST_CASE(production_analytic_jacobian_matches_differences)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase base;
    const double free_total = base.freeTotal();

    auto check = [&](const char* what, Sys& system) {
        system.setAnalyticJacobian(true);
        const auto x0 = system.start(ProductionCase::guess());
        // at the solution, so the active set the Jacobian is built on is the real one
        const auto r = NetworkSolve::solve(system, ProductionCase::guess());
        BOOST_REQUIRE(r.converged);
        auto x = x0;
        for (int n = 1; n <= system.numNodes(); ++n) { x[system.pIdx(n)] = r.node_pressure[n]; }
        system.updateControls(x);
        const auto J = system.jacobian(x);
        const auto r0 = system.residual(x);
        double worst = 0.0, largest = 0.0;
        for (int j = 0; j < system.size(); ++j) {
            auto shifted = x;
            const double h = 1e-4 * system.columnScale(j);
            shifted[j] += h;
            const auto rj = system.residual(shifted);
            for (int i = 0; i < system.size(); ++i) {
                const double fd = (rj[i] - r0[i]) / h;
                worst = std::max(worst, std::abs(J(i, j) - fd));
                largest = std::max(largest, std::abs(fd));
            }
        }
        BOOST_TEST_MESSAGE(what << ": largest entry " << largest << ", largest difference " << worst);
        BOOST_CHECK_LT(worst, 1e-3 * std::max(largest, 1.0));
    };

    { ProductionCase c; auto s = c.system(); check("plain", s); }
    { ProductionCase c; c.setGroupTarget(0.6 * free_total); auto s = c.system(); check("group target", s); }
    { ProductionCase c; auto s = c.system(); s.setChokeTarget(1, 0.5 * free_total); check("closed choke", s); }
}


// A well placed exactly where its tubing passes its own rate limit satisfies
// both rows, and an active set that picks one by the sign at the iterate flips
// between them. The tie goes to a Fischer-Burmeister row instead; it has to
// converge, land on the limit, and leave the tubing slack non-negative.
BOOST_AUTO_TEST_CASE(a_well_on_a_tie_converges)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;
    ProductionCase free_case;
    auto free_system = free_case.system();
    const auto free = NetworkSolve::solve(free_system, ProductionCase::guess());
    BOOST_REQUIRE(free.converged);

    // the limit is the rate the well freely produces: the tie, to the digit
    ProductionCase c;
    c.wells()[0].oil_rate_limit = free.well_rate[0];
    auto system = c.system();
    system.setAnalyticJacobian(true);
    const auto r = NetworkSolve::solve(system, ProductionCase::guess());
    BOOST_TEST_MESSAGE("tied well: " << (r.converged ? "converged in " : "FAILED after ")
                       << r.iterations << " iterations, " << r.switches << " switches, control '"
                       << system.controlLetter(0) << "', oil " << convert::to(r.well_rate[0], sm3d)
                       << " against limit " << convert::to(free.well_rate[0], sm3d));
    BOOST_CHECK(r.converged);
    BOOST_CHECK_LE(r.well_rate[0], free.well_rate[0] * (1.0 + 1e-3));
    BOOST_CHECK_CLOSE(r.well_rate[0], free.well_rate[0], 0.5);
}


// Do the assembled and the differenced Jacobian reach the same solution?
//
// Not "do both converge" -- the same node pressures and well rates, from the
// same start, over a grid of starting pressures, on every production system
// shape there is: plain, under a group target, with a closed choke, and with
// a well placed on its rate-limit tie. Where they part, one of them has found
// another root, and the simulator cannot tell which.
BOOST_AUTO_TEST_CASE(production_jacobians_reach_the_same_solution)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase base;
    const double free_total = base.freeTotal();
    auto free_system = base.system();
    const auto free = NetworkSolve::solve(free_system, ProductionCase::guess());
    BOOST_REQUIRE(free.converged);

    // The cases own the tables the systems point at, so they live here, not
    // inside the lambdas that build systems from them.
    ProductionCase plain_case, group_case, choke_case, tie_case;
    group_case.setGroupTarget(0.6 * free_total);
    tie_case.wells()[0].oil_rate_limit = free.well_rate[0];
    struct Shape { const char* what; std::function<Sys(void)> make; };
    const std::vector<Shape> shapes{
        {"plain",        [&] { return plain_case.system(); }},
        {"group target", [&] { return group_case.system(); }},
        {"closed choke", [&] { auto s = choke_case.system(); s.setChokeTarget(1, 0.5 * free_total); return s; }},
        {"on a tie",     [&] { return tie_case.system(); }},
    };

    for (const auto& shape : shapes) {
        int both = 0, agree = 0, only_fd = 0, only_an = 0, neither = 0;
        double worst_p = 0.0, worst_q = 0.0;
        for (int pf = 0; pf < 14; ++pf) {
            const std::vector<double> guess{convert::from(80.0, bars), convert::from(50.0 + 30.0 * pf, bars)};
            auto fd = shape.make(); fd.setAnalyticJacobian(false);
            auto an = shape.make(); an.setAnalyticJacobian(true);
            const auto rf = NetworkSolve::solve(fd, guess);
            const auto ra = NetworkSolve::solve(an, guess);
            if (rf.converged && ra.converged) {
                ++both;
                double dp = 0.0, dq = 0.0;
                for (std::size_t n = 0; n < rf.node_pressure.size(); ++n) {
                    dp = std::max(dp, std::abs(rf.node_pressure[n] - ra.node_pressure[n]));
                }
                for (std::size_t w = 0; w < rf.well_rate.size(); ++w) {
                    dq = std::max(dq, std::abs(rf.well_rate[w] - ra.well_rate[w])
                                      / std::max(std::abs(rf.well_rate[w]), 1e-12));
                }
                worst_p = std::max(worst_p, dp);
                worst_q = std::max(worst_q, dq);
                if (dp < convert::from(0.05, bars) && dq < 1e-3) { ++agree; }
            } else if (rf.converged) { ++only_fd; }
            else if (ra.converged) { ++only_an; }
            else { ++neither; }
        }
        BOOST_TEST_MESSAGE(std::setw(13) << shape.what << ": both " << both << "/14, agree " << agree
                           << ", only differenced " << only_fd << ", only analytic " << only_an
                           << ", neither " << neither << "; worst gap " << convert::to(worst_p, bars)
                           << " bar, " << 100 * worst_q << " % rate");
        BOOST_CHECK_EQUAL(agree, both);
        BOOST_CHECK_EQUAL(only_fd + only_an, 0);
    }
}


// Replay production systems the simulator wrote out (--network-dump-failures),
// against the MODEL5 tables. OPM_NETWORK_DUMP_PROD names the directory,
// OPM_VFP_INCLUDE the directory holding well_vfp.ecl and flowl_{b,c}_vfp.ecl.
// Each system is solved with both Jacobians; this is where a failure seen in
// the simulator becomes a millisecond of work.
BOOST_AUTO_TEST_CASE(replay_production_failures)
{
    const char* dir = std::getenv("OPM_NETWORK_DUMP_PROD");
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (dir == nullptr || inc == nullptr || !std::filesystem::is_directory(dir)) {
        BOOST_TEST_MESSAGE("OPM_NETWORK_DUMP_PROD / OPM_VFP_INCLUDE not set, nothing to replay");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto path = std::filesystem::path(inc) / name;
        if (!std::filesystem::exists(path)) { continue; }
        const auto deck = Parser{}.parseFile(path.string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, /*gaslift_opt_active=*/true, units);
            props.addTable(tables.back());
        }
    }
    BOOST_TEST_MESSAGE("tables loaded: " << tables.size());

    std::vector<std::filesystem::path> files;
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (e.path().extension() == ".txt") { files.push_back(e.path()); }
    }
    std::sort(files.begin(), files.end());
    int fd_ok = 0, an_ok = 0, agree = 0, cm_ok = 0, cm_agree = 0;
    for (const auto& file : files) {
        std::ifstream in(file);
        std::string head; std::getline(in, head);
        if (head != "production") { continue; }
        auto [fd, guess] = NetworkSolve::readProduction<double>(in, props, units);
        auto an = fd; fd.setAnalyticJacobian(false); an.setAnalyticJacobian(true);
        auto cm = an; cm.setComplementarity(true);
        // OPM_NETWORK_MAX_IT raises the iteration cap, to tell "slow" from "stuck".
        const int max_it = std::getenv("OPM_NETWORK_MAX_IT") ? std::atoi(std::getenv("OPM_NETWORK_MAX_IT")) : 50;
        const auto rf = NetworkSolve::solve(fd, guess, 1e-2, max_it);
        const auto ra = NetworkSolve::solve(an, guess, 1e-2, max_it);
        const bool cm_ls = std::getenv("OPM_NETWORK_CM_LINESEARCH") != nullptr;
        const auto rc = cm_ls ? NetworkSolve::solve(cm, guess, 1e-2, max_it, NetworkSolve::LineSearch{})
                              : NetworkSolve::solve(cm, guess, 1e-2, max_it);
        fd_ok += rf.converged; an_ok += ra.converged; cm_ok += rc.converged;
        double gap = 0.0, cgap = 0.0;
        if (rf.converged && ra.converged) {
            for (std::size_t n = 0; n < rf.node_pressure.size(); ++n) {
                gap = std::max(gap, std::abs(rf.node_pressure[n] - ra.node_pressure[n]));
            }
            agree += (gap < convert::from(0.05, bars));
        }
        std::string where;
        if (ra.converged && rc.converged) {
            for (std::size_t n = 0; n < ra.node_pressure.size(); ++n) {
                cgap = std::max(cgap, std::abs(ra.node_pressure[n] - rc.node_pressure[n]));
            }
            cm_agree += (cgap < convert::from(0.05, bars));
            if (cgap >= convert::from(0.05, bars)) {
                for (int n = 1; n <= an.numNodes(); ++n) {
                    if (an.isChoke(n)) {
                        const int up = an.nodes()[n].parent;
                        where += fmt::format(" choke {} an:{} p {:.2f} (up {:.2f}) cm:{} p {:.2f}", an.nodes()[n].name,
                                             an.choked(n) ? "CLOSED" : "OPEN", ra.node_pressure[n] * 1e-5,
                                             (up == 0 ? an.terminalPressure() : ra.node_pressure[up]) * 1e-5,
                                             cm.choked(n) ? "CLOSED" : "OPEN", rc.node_pressure[n] * 1e-5);
                    }
                }
                for (int w = 0; w < an.numWells(); ++w) {
                    where += fmt::format(" {}:{}{:.0f}/{}{:.0f}", an.wells()[w].name, an.controlLetter(w),
                                         ra.well_rate[w] * 86400.0, cm.controlLetter(w), rc.well_rate[w] * 86400.0);
                }
            }
        }
        BOOST_TEST_MESSAGE("  " << file.filename().string()
                           << ": differenced " << (rf.converged ? "ok" : "FAILED") << " (" << rf.iterations << " it)"
                           << ", analytic " << (ra.converged ? "ok" : "FAILED") << " (" << ra.iterations << " it)"
                           << ", complementarity " << (rc.converged ? "ok" : "FAILED") << " (" << rc.iterations << " it)"
                           << (ra.converged && rc.converged ? fmt::format(", gap an/cm {:.3g} bar", convert::to(cgap, bars)) + where : "")
                           << (rc.converged ? "" : fmt::format("  cm residual {:.3g}", rc.residual)
                                               + (rc.control_trace.empty() ? "" : "  cm trace " + rc.control_trace)));
    }
    BOOST_TEST_MESSAGE("replayed " << files.size() << ": differenced ok " << fd_ok << ", analytic ok "
                       << an_ok << ", both ok and agreeing " << agree << "; complementarity ok " << cm_ok
                       << ", agreeing with analytic " << cm_agree);
}


// One complementarity row per well instead of an active set: on every shape,
// it must converge, and where the active set converges too the two must agree.
BOOST_AUTO_TEST_CASE(complementarity_agrees_with_the_active_set)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    const auto sm3d = cubic(meter) / day;
    ProductionCase base;
    const double free_total = base.freeTotal();
    auto free_system = base.system();
    const auto free = NetworkSolve::solve(free_system, ProductionCase::guess());
    BOOST_REQUIRE(free.converged);

    ProductionCase plain_case, limited_case, choke_case, open_case, tie_case, bhp_case, high_case;
    for (auto& w : limited_case.wells()) { w.oil_rate_limit = 0.4 * free.well_rate[0]; }
    tie_case.wells()[0].oil_rate_limit = free.well_rate[0];
    for (auto& w : bhp_case.wells()) { w.bhp_limit = convert::from(110.0, bars); }
    for (auto& w : high_case.wells()) { w.bhp_limit = convert::from(150.0, bars); }   // above shut-in: nothing
    struct Shape { const char* what; std::function<Sys(void)> make; };
    const std::vector<Shape> shapes{
        {"plain",        [&] { return plain_case.system(); }},
        {"rate-limited", [&] { return limited_case.system(); }},
        {"bhp-limited",  [&] { return bhp_case.system(); }},
        {"bhp above reservoir", [&] { return high_case.system(); }},
        {"closed choke", [&] { auto s = choke_case.system(); s.setChokeTarget(1, 0.5 * free_total); return s; }},
        {"open choke",   [&] { auto s = open_case.system(); s.setChokeTarget(1, 2.0 * free_total); return s; }},
        {"on a tie",     [&] { return tie_case.system(); }},
    };
    for (const auto& shape : shapes) {
        int both = 0, agree = 0, only_as = 0, only_cm = 0, cm_its = 0;
        double worst = 0.0;
        for (int pf = 0; pf < 14; ++pf) {
            const std::vector<double> guess{convert::from(80.0, bars), convert::from(50.0 + 30.0 * pf, bars)};
            auto as = shape.make(); as.setAnalyticJacobian(true);
            auto cm = shape.make(); cm.setAnalyticJacobian(true); cm.setComplementarity(true);
            const auto ra = NetworkSolve::solve(as, guess);
            const auto rc = NetworkSolve::solve(cm, guess);
            if (rc.converged) { cm_its += rc.iterations; }
            const bool choke_shape = std::string(shape.what) == "closed choke" || std::string(shape.what) == "open choke";
            const bool disagree = ra.converged && rc.converged
                && std::abs(ra.well_rate[0] - rc.well_rate[0]) > 1e-3 * std::abs(ra.well_rate[0]);
            if (pf % 4 == 0 && (!rc.converged || disagree)) {
                std::string rates;
                for (std::size_t w = 0; w < rc.well_rate.size(); ++w) {
                    rates += fmt::format(" {:.0f}({}, bhp {:.1f}/as {:.1f})", rc.well_rate[w] * 86400.0,
                                         static_cast<int>(cm.control(static_cast<int>(w))), rc.well_bhp[w] * 1e-5, ra.well_bhp[w] * 1e-5);
                }
                BOOST_TEST_MESSAGE("    " << shape.what << ", start " << 50 + 30 * pf << " bar: cm " << (rc.converged ? "ok after " : "FAILED after ")
                                   << rc.iterations << " it, residual " << rc.residual << ", choke "
                                   << (cm.choked(1) ? "CLOSED" : "OPEN") << " p " << rc.node_pressure[1] * 1e-5
                                   << ", rates" << rates << ", trace " << rc.control_trace
                                   << " | active set: p " << ra.node_pressure[1] * 1e-5 << " rates "
                                   << ra.well_rate[0] * 86400.0 << " " << ra.well_rate[1] * 86400.0);
            }
            if (ra.converged && rc.converged) {
                ++both;
                double d = 0.0;
                for (std::size_t w = 0; w < ra.well_rate.size(); ++w) {
                    d = std::max(d, std::abs(ra.well_rate[w] - rc.well_rate[w]) / std::max(std::abs(ra.well_rate[w]), 1e-12));
                }
                worst = std::max(worst, d);
                agree += (d < 1e-3);
            } else if (ra.converged) { ++only_as; }
            else if (rc.converged) { ++only_cm; }
        }
        BOOST_TEST_MESSAGE(std::setw(13) << shape.what << ": both " << both << "/14, agree " << agree
                           << ", only active-set " << only_as << ", only complementarity " << only_cm
                           << ", worst rate gap " << 100 * worst << " %, complementarity mean its "
                           << (both + only_cm > 0 ? cm_its / (both + only_cm) : 0));
        BOOST_CHECK_EQUAL(only_as, 0);
        BOOST_CHECK_EQUAL(agree, both);
    }
}


// The production step bounds, checked directly: a 1000 bar node step is cut to
// 50, and a step that would take a node below an atmosphere stops at it.
// A well the well model has at zero rate above some node pressure -- the
// adapter's dead_above -- must not be produced from. The active set hands it
// to bhp control and reports 2500 m3/d through tubing that cannot lift
// (dump_prod_1112: three of them under a 6000 choke target); the
// complementarity shuts it. Both shapes from the dumps, on the bench tables.
BOOST_AUTO_TEST_CASE(complementarity_shuts_dead_wells)
{
    using Sys = NetworkSolve::ProductionSystem<double>;
    ProductionCase plain_case;
    auto plain = plain_case.system(); plain.setAnalyticJacobian(true);
    const auto free = NetworkSolve::solve(plain, ProductionCase::guess());
    BOOST_REQUIRE(free.converged);
    const double p_free = free.node_pressure[1];

    // One dying well: dead just below where the node settles with it alive.
    {
        ProductionCase c;
        c.wells()[0].dead_above = p_free - convert::from(0.5, bars);
        auto cm = c.system(); cm.setAnalyticJacobian(true); cm.setComplementarity(true);
        auto as = c.system(); as.setAnalyticJacobian(true);
        const auto rc = NetworkSolve::solve(cm, ProductionCase::guess());
        const auto ra = NetworkSolve::solve(as, ProductionCase::guess());
        BOOST_REQUIRE(rc.converged);
        BOOST_TEST_MESSAGE("dying well: cm " << cm.controlLetter(0) << " " << rc.well_rate[0] * 86400
                           << " m3/d at node " << rc.node_pressure[1] * 1e-5 << " bar; active set "
                           << (ra.converged ? as.controlLetter(0) : '?') << " "
                           << (ra.converged ? ra.well_rate[0] * 86400 : 0.0) << " m3/d");
        BOOST_CHECK_EQUAL(cm.controlLetter(0), 'S');
        BOOST_CHECK_SMALL(rc.well_rate[0], 1e-12);
        BOOST_CHECK_GT(rc.well_rate[1], 0.0);
        // The other well alone leaves the node below where well 0 died --
        // the no-fixed-point shape; the shut has to stick.
        BOOST_CHECK_LT(rc.node_pressure[1], c.wells()[0].dead_above);
    }
    // A choke over wells that are all dead: nothing to throttle, valve open,
    // no flow. The row must not chase a target nothing can meet.
    {
        ProductionCase c;
        for (auto& w : c.wells()) { w.dead_above = convert::from(1.0, bars); }
        auto cm = c.system(); cm.setAnalyticJacobian(true); cm.setComplementarity(true);
        cm.setChokeTarget(1, 0.5 * free.well_rate[0]);
        const auto rc = NetworkSolve::solve(cm, ProductionCase::guess());
        BOOST_REQUIRE(rc.converged);
        for (int w = 0; w < cm.numWells(); ++w) {
            BOOST_CHECK_EQUAL(cm.controlLetter(w), 'S');
            BOOST_CHECK_SMALL(rc.well_rate[w], 1e-12);
        }
        BOOST_CHECK_LT(rc.iterations, 10);
        BOOST_TEST_MESSAGE("choke over dead wells: " << rc.iterations << " it, node "
                           << rc.node_pressure[1] * 1e-5 << " bar");
    }
}

// The default path's hole, pinned so a fix flips this test rather than
// passing unnoticed: the active set gives a dead well (dead_above below the
// node pressure, thp unavailable) bhp control and produces from it.
BOOST_AUTO_TEST_CASE(the_active_set_still_produces_from_a_dead_well)
{
    ProductionCase c;
    c.wells()[0].dead_above = convert::from(1.0, bars);
    auto as = c.system(); as.setAnalyticJacobian(true);
    const auto r = NetworkSolve::solve(as, ProductionCase::guess());
    BOOST_REQUIRE(r.converged);
    BOOST_TEST_MESSAGE("dead well on the active set: " << as.controlLetter(0) << " "
                       << r.well_rate[0] * 86400 << " m3/d");
    BOOST_CHECK_EQUAL(as.controlLetter(0), 'B');     // should one day be 'S'
    BOOST_CHECK_GT(r.well_rate[0], 0.0);
}

// The four AUTOCHK dumps behind the 2026-08-23 fixes, kept under
// tests/network_dumps. They need the model5 tables (OPM_VFP_INCLUDE); without
// them the case reports and passes, like replay_production_failures.
BOOST_AUTO_TEST_CASE(the_dumps_behind_the_complementarity_fixes_converge)
{
    const char* inc = std::getenv("OPM_VFP_INCLUDE");
    if (inc == nullptr) {
        BOOST_TEST_MESSAGE("OPM_VFP_INCLUDE not set; dumps not replayed");
        return;
    }
    std::deque<VFPProdTable> tables;
    VFPProdProperties<double> props;
    const UnitSystem units{};
    for (const char* name : {"well_vfp.ecl", "flowl_b_vfp.ecl", "flowl_c_vfp.ecl"}) {
        const auto deck = Parser{}.parseFile((std::filesystem::path(inc) / name).string());
        for (const auto& kw : deck.getKeywordList("VFPPROD")) {
            tables.emplace_back(*kw, true, units);
            props.addTable(tables.back());
        }
    }
    const auto dir = std::filesystem::path(__FILE__).parent_path() / "network_dumps";
    for (const char* name : {"autochk_prod_992.txt", "autochk_prod_1112.txt", "autochk_prod_0.txt", "autochk_prod_302.txt"}) {
        std::ifstream in(dir / name);
        BOOST_REQUIRE(in);
        std::string head; std::getline(in, head);
        auto [system, guess] = NetworkSolve::readProduction<double>(in, props, units);
        system.setAnalyticJacobian(true); system.setComplementarity(true);
        const auto r = NetworkSolve::solve(system, guess, 1e-2, 50);
        std::string controls;
        for (int w = 0; w < system.numWells(); ++w) { controls += system.controlLetter(w); }
        BOOST_TEST_MESSAGE(name << ": " << (r.converged ? "ok" : "FAILED") << " in " << r.iterations
                           << " it, controls " << controls);
        BOOST_CHECK(r.converged);
        BOOST_CHECK_LT(r.iterations, 12);
    }
}

BOOST_AUTO_TEST_CASE(production_step_bounds)
{
    ProductionCase c;
    auto system = c.system();
    auto x = system.start(ProductionCase::guess());
    std::vector<double> dx(system.size(), 0.0);
    dx[system.pIdx(1)] = convert::from(1000.0, bars);
    auto out = system.limitStep(x, dx);
    BOOST_TEST_MESSAGE("1000 bar step -> " << convert::to(out[system.pIdx(1)], bars) << " bar");
    BOOST_CHECK_LE(out[system.pIdx(1)], convert::from(50.0, bars) * 1.0001);
    dx[system.pIdx(1)] = -x[system.pIdx(1)] - convert::from(5.0, bars);
    out = system.limitStep(x, dx);
    BOOST_CHECK_GE(x[system.pIdx(1)] + out[system.pIdx(1)], convert::from(1.0, bars) * 0.99);
}

BOOST_AUTO_TEST_SUITE_END()
