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
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/NetworkNodePressureUpdater.hpp>
#include <opm/simulators/wells/NetworkAndersonAcceleration.hpp>

#include <algorithm>
#include <array>
#include <deque>
#include <filesystem>
#include <fmt/format.h>
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
    const std::vector<int>& children(const int n) const { return children_[n]; }
    const std::vector<int>& wellsAt(const int n) const { return wells_at_[n]; }
    const std::vector<int>& solvedNodes() const { return solved_; }
    double terminalPressure() const { return terminal_pressure_; }
    double groupTarget() const { return group_target_; }
    bool hasTable(const Node& n) const { return n.vfp_table != kNoTable && axes_.count(n.vfp_table); }

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
    double columnScale(const int) const { return kPressureScale; }
    State limitStep(const State&, const State& dx) const { return dx; }

private:
    const NetworkCase& case_;
};

/// The same case without the eliminations.
///
///   unknowns   pressure of every non-terminal node                    nP
///              rate through every node's parent branch                nP
///              (rate, bhp) for every well                            2W
///              the group multiplier, when a target is active          1
///
///   equations  branch drop    p_n - VFP(thp = p_parent, q_n) = 0      nP
///              node balance   q_n - sum(children) - sum(wells) = 0    nP
///              inflow perf.   q_w - ipr_w(bhp_w) = 0                   W
///              control        whichever of THP / BHP / RATE / GRUP     W
///              group          sum(q_w) - target = 0                    1
class FullProblem
{
public:
    explicit FullProblem(const NetworkCase& c)
        : case_(c)
        , nodes_(static_cast<int>(c.nodes().size()) - 1)
        , wells_(static_cast<int>(c.wells().size()))
        , grouped_(c.groupTarget() > 0.0)
        , controls_(c.wells().size(), Control::Thp)
    {
        if (grouped_) {
            std::fill(controls_.begin(), controls_.end(), Control::Grup);
        }
    }

    static constexpr const char* name = "full";

    /// Which equation closes each well.
    enum class Control { Thp, Bhp, Rate, Grup };

    int size() const { return 2 * nodes_ + 2 * wells_ + (grouped_ ? 1 : 0); }

    /// Keep every iterate inside the box the tables describe, instead of
    /// clamping the lookups. See table_bounds_want_to_be_constraints.
    void setEnforceBounds(const bool on) { enforce_bounds_ = on; }

    int pIdx(const int node) const { return node - 1; }
    int qIdx(const int node) const { return nodes_ + node - 1; }
    int qwIdx(const int w) const { return 2 * nodes_ + w; }
    int bhpIdx(const int w) const { return 2 * nodes_ + wells_ + w; }
    int lambdaIdx() const { return 2 * nodes_ + 2 * wells_; }

    State residual(const State& x) const
    {
        const auto& nodes = case_.nodes();
        const auto& wells = case_.wells();
        State r(size(), 0.0);

        auto pressure = [&](const int n) {
            return n == 0 ? case_.terminalPressure() : x[pIdx(n)];
        };

        for (int n = 1; n <= nodes_; ++n) {
            const auto& node = nodes[n];
            const double upstream = pressure(node.parent);
            r[n - 1] = case_.hasTable(node)
                ? x[pIdx(n)] - case_.tableBhp(node.vfp_table, upstream, x[qIdx(n)])
                : x[pIdx(n)] - upstream;

            double balance = x[qIdx(n)];
            for (const int c : case_.children(n)) {
                balance -= x[qIdx(c)];
            }
            for (const int w : case_.wellsAt(n)) {
                balance -= x[qwIdx(w)];
            }
            r[nodes_ + n - 1] = balance;
        }

        double injected = 0.0;
        for (int w = 0; w < wells_; ++w) {
            const auto& well = wells[w];
            const double q = x[qwIdx(w)];
            const double bhp = x[bhpIdx(w)];
            injected += q;

            r[2 * nodes_ + w] = (q - NetworkCase::ipr(well, bhp)) / kRateScale;

            double& control = r[2 * nodes_ + wells_ + w];
            switch (controls_[w]) {
            case Control::Thp:
                control = (bhp - case_.tableBhp(well.vfp_table, pressure(well.node), q))
                        / kPressureScale;
                break;
            case Control::Bhp:
                control = (bhp - well.bhp_limit) / kPressureScale;
                break;
            case Control::Rate:
                control = (q - well.rate_limit) / kRateScale;
                break;
            case Control::Grup:
                control = (q - well.guide * x[lambdaIdx()]) / kRateScale;
                break;
            }
        }

        if (grouped_) {
            // With nobody on group control the multiplier is free, so pin it
            // rather than hand the Newton a singular column.
            const bool any_grouped = std::find(controls_.begin(), controls_.end(), Control::Grup)
                                   != controls_.end();
            r[lambdaIdx()] = any_grouped ? (injected - case_.groupTarget()) / kRateScale
                                         : (x[lambdaIdx()] - lambda0()) / kRateScale;
        }

        for (int n = 0; n < nodes_; ++n) {
            r[n] /= kPressureScale;
            r[nodes_ + n] /= kRateScale;
        }
        return r;
    }

    /// Reselect each well's control: the most restrictive violated limit wins,
    /// which is the same rule the eliminated form applies as a clamp. Picking by
    /// a fixed priority instead makes the set chatter and the Newton never ends.
    /// Returns true if anything moved, so the caller can tell an active-set
    /// change from a converged step.
    bool updateControls(const State& x)
    {
        const auto& wells = case_.wells();
        bool changed = false;
        for (int w = 0; w < wells_; ++w) {
            const auto& well = wells[w];
            const double q = x[qwIdx(w)];

            auto wanted = Control::Thp;
            double smallest = std::numeric_limits<double>::max();
            auto consider = [&](const bool violated, const double implied, const Control c) {
                if (violated && implied < smallest) {
                    smallest = implied;
                    wanted = c;
                }
            };
            consider(x[bhpIdx(w)] > well.bhp_limit, NetworkCase::ipr(well, well.bhp_limit),
                     Control::Bhp);
            consider(q > well.rate_limit, well.rate_limit, Control::Rate);
            if (grouped_) {
                // Inclusive: at the solution the rate equals the share exactly,
                // and a strict test would flip the control every iteration.
                const double share = well.guide * x[lambdaIdx()];
                consider(q >= share * (1.0 - 1e-9), share, Control::Grup);
            }

            changed |= (wanted != controls_[w]);
            controls_[w] = wanted;
        }
        return changed;
    }

    /// Everything derived from the node pressures the eliminated form starts
    /// from, so the two formulations really do start from the same place.
    State start(const State& applied) const
    {
        const auto& wells = case_.wells();
        const auto p = case_.nodePressures(applied);
        const double q_guess = convert::from(1.0e5, cubic(meter) / day);

        State x(size(), 0.0);
        for (int n = 1; n <= nodes_; ++n) {
            x[pIdx(n)] = p[n];
        }
        std::vector<double> well_rate(wells_);
        for (int w = 0; w < wells_; ++w) {
            const auto& well = wells[w];
            x[bhpIdx(w)] = case_.tableBhp(well.vfp_table, p[well.node], q_guess);
            x[qwIdx(w)] = std::clamp(NetworkCase::ipr(well, x[bhpIdx(w)]), 0.0, well.rate_limit);
            well_rate[w] = x[qwIdx(w)];
        }
        const auto branch = case_.branchFlows(well_rate);
        for (int n = 1; n <= nodes_; ++n) {
            x[qIdx(n)] = branch[n];
        }
        if (grouped_) {
            x[lambdaIdx()] = lambda0();
        }
        return x;
    }

    State pressures(const State& x) const
    {
        State p;
        for (const int n : case_.solvedNodes()) {
            p.push_back(x[pIdx(n)]);
        }
        return p;
    }

    /// Natural magnitude of unknown i. Pressures and bhp in bar, everything that
    /// is a rate in kRateScale -- this is what lets one step cap and one trust
    /// radius apply to a vector holding both.
    double columnScale(const int i) const
    {
        const bool is_pressure = (i < nodes_) || (i >= bhpIdx(0) && i < lambdaIdx());
        return is_pressure ? kPressureScale : kRateScale;
    }

    /// Keep the branch flows inside the box the tables describe, by projecting
    /// the offending components of the step rather than scaling all of it --
    /// one binding rate should not throttle the pressure updates too.
    ///
    /// Only the branch flows: a well's own rate limit already has a control
    /// equation, and bounding it would stop that control ever activating. An
    /// iterate already outside is left alone, or the step would be zero.
    State limitStep(const State& x, const State& dx) const
    {
        if (!enforce_bounds_) {
            return dx;
        }
        State limited = dx;
        for (int n = 1; n <= nodes_; ++n) {
            const auto& node = case_.nodes()[n];
            if (!case_.hasTable(node)) {
                continue;
            }
            const int i = qIdx(n);
            const double hi = case_.maxFlow(node.vfp_table);
            if (x[i] <= hi && x[i] + limited[i] > hi) {
                limited[i] = hi - x[i];
            }
            if (x[i] >= 0.0 && x[i] + limited[i] < 0.0) {
                limited[i] = -x[i];
            }
        }
        return limited;
    }

private:
    double lambda0() const
    {
        double guides = 0.0;
        for (const auto& w : case_.wells()) {
            guides += w.guide;
        }
        return guides > 0.0 ? case_.groupTarget() / guides : 0.0;
    }

    const NetworkCase& case_;
    int nodes_;
    int wells_;
    bool grouped_;
    bool enforce_bounds_ = false;
    std::vector<Control> controls_;
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

Result anderson(const EliminatedProblem& problem, State x, const int depth)
{
    NetworkAndersonAccelerator<double> accelerator;
    accelerator.setDepth(depth);
    for (int it = 1; it <= kMaxIter; ++it) {
        const auto g = problem.G(x);
        if (normMax(g - x) < kTol * kPressureScale) {
            return {true, it, x};
        }
        x = accelerator.next(x, g);
    }
    return {false, kMaxIter + 1, x};
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
            return {true, it, problem.pressures(x)};
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
    const auto acc = anderson(eliminated, kStart, 4);
    const auto full_step = newton(eliminated, kStart, FullStep{});
    const auto capped = newton(eliminated, kStart, CappedStep{});
    const auto search = newton(eliminated, kStart, LineSearch{});
    const auto armijo = newton(eliminated, kStart, LineSearch{/*sufficient_decrease=*/true});
    const auto region = newton(eliminated, kStart, TrustRegion{});

    report("damped (omega 0.1)", fixed_point);
    report("bracketing (shipped)", bracket);
    report("anderson (depth 4)", acc);
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
    BOOST_CHECK_EQUAL(f_step, n);
    BOOST_CHECK_GE(e_search, n - 1);
    BOOST_CHECK_GE(f_search, n - 1);
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

    // Clamping is far worse than leaving the tables alone; bounding beats both,
    // though it does not recover the whole grid either.
    BOOST_CHECK_LT(unclamped, n);
    BOOST_CHECK_LT(with_clamp, unclamped / 2);
    BOOST_CHECK_GE(with_bounds, unclamped);
    BOOST_CHECK_GT(with_bounds, 3 * n / 4);

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

BOOST_AUTO_TEST_SUITE_END()
