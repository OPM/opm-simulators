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
// ---------------------------------------------------------------------------

using Vec = std::array<double, 2>;   // (p_M5S, p_M5N), SI

constexpr int kWellTable = 1;   // VFPINJ on the wells' tubing
constexpr int kM5nTable = 2;    // M5S -> M5N
constexpr int kM5sTable = 3;    // PLAT-A -> M5S

// One injector: linear IPR against VFPINJ 1, then the control logic.
struct WellProxy
{
    std::string name;
    int node;              // 0 = on G1 (sees M5S), 1 = on F1 (sees M5N)
    double q_ref;          // E100 rate at p_ref                 [sm3/s]
    double p_ref;          // E100 node pressure                 [Pa]
    double dq_dbhp;        // IPR slope, the stiffness knob      [sm3/s/Pa]
    double bhp_limit;
    double rate_limit;
    double bhp_ref = 0.0;  // bhp at (p_ref, q_ref); the network fills it in
};

class GasInjectionNetwork
{
public:
    GasInjectionNetwork()
    {
        addTable(vfp_well);   // kWellTable
        addTable(vfp_m5n);    // kM5nTable
        addTable(vfp_m5s);    // kM5sTable

        const auto sm3_day = cubic(meter) / day;
        // Calibration point: Eclipse 100, day 31 (opm-tests/eclref).
        wells_ = {
            WellProxy{"G-3H", 0, convert::from(4.894e5, sm3_day), convert::from(209.4, bars), 0.0, 0.0, 0.0},
            WellProxy{"G-4H", 0, convert::from(4.893e5, sm3_day), convert::from(209.4, bars), 0.0, 0.0, 0.0},
            WellProxy{"F-1H", 1, convert::from(2.764e5, sm3_day), convert::from(204.2, bars), 0.0, 0.0, 0.0},
            WellProxy{"F-2H", 1, convert::from(2.769e5, sm3_day), convert::from(204.2, bars), 0.0, 0.0, 0.0},
        };
        for (auto& w : wells_) {
            w.bhp_limit = convert::from(425.0, bars);
            w.rate_limit = convert::from(1.0e6, sm3_day);
            w.bhp_ref = wellBhp(w.p_ref, w.q_ref);
            w.dq_dbhp = stiffness_;
        }
    }

    /// IPR slope shared by all wells [sm3/d per bar], the knob the bench exists for.
    void setStiffness(const double dq_dbhp_sm3_day_per_bar)
    {
        stiffness_ = convert::from(dq_dbhp_sm3_day_per_bar, cubic(meter) / day) / convert::from(1.0, bars);
        for (auto& w : wells_) {
            w.dq_dbhp = stiffness_;
        }
    }

    /// Clamp table lookups to the flow and THP axes, as the simulator's network
    /// pressure computation does. Off the axes the tables are zero-filled and the
    /// interpolation runs away, so leaving this off admits spurious roots; turning
    /// it on removes them but flattens the residual, which a Newton cannot climb
    /// off. See table_bounds_want_to_be_constraints.
    void setClampToAxes(const bool on) { clamp_to_axes_ = on; }

    void setGroupTarget(const double sm3_day) { group_target_ = convert::from(sm3_day, cubic(meter) / day); }

    /// The fixed-point map: applied node pressures in, computed node pressures out.
    Vec G(const Vec& p) const
    {
        const auto q = rates(p);
        const double q_total = q[0] + q[1] + q[2] + q[3];
        const double q_m5n = q[2] + q[3];
        Vec out;
        out[0] = branchBhp(kM5sTable, terminal_, q_total);
        out[1] = branchBhp(kM5nTable, p[0], q_m5n);
        return out;
    }

    Vec residual(const Vec& p) const
    {
        const auto g = G(p);
        return {g[0] - p[0], g[1] - p[1]};
    }

    /// Per-well rates at the applied node pressures, group target applied.
    std::array<double, 4> rates(const Vec& p) const
    {
        std::array<double, 4> q{};
        for (std::size_t i = 0; i < wells_.size(); ++i) {
            q[i] = wellRate(wells_[i], p[wells_[i].node]);
        }
        const double sum = q[0] + q[1] + q[2] + q[3];
        if (group_target_ > 0.0 && sum > group_target_) {
            // GRUP control: share the target in proportion to the unconstrained rates.
            for (auto& qi : q) {
                qi *= group_target_ / sum;
            }
        }
        return q;
    }

    /// Downstream pressure of a branch: VFPINJ table `table` at this upstream
    /// pressure and rate. Table kWellTable is the wells' own tubing.
    double branch(const int table, const double thp, const double q) const
    { return branchBhp(table, thp, q); }

    const std::vector<WellProxy>& wells() const { return wells_; }
    double terminal() const { return terminal_; }

    /// The wells' rate response to their own bhp, without any control logic.
    static double ipr(const WellProxy& w, const double bhp)
    { return w.q_ref + w.dq_dbhp * (bhp - w.bhp_ref); }

private:
    // VFPInjProperties keeps a reference_wrapper, so the tables have to outlive it and
    // must not move -- hence the deque.
    void addTable(const std::string& s)
    {
        decks_.push_back(Parser{}.parseString(s));
        tables_.emplace_back(decks_.back()["VFPINJ"].front(), UnitSystem{});
        props_.addTable(tables_.back());

        const auto& t = tables_.back();
        axes_[t.getTableNum()] = Axes{t.getFloAxis().front(), t.getFloAxis().back(),
                                      t.getTHPAxis().front(), t.getTHPAxis().back()};
    }

    double branchBhp(const int table, const double thp, const double q) const
    {
        if (!clamp_to_axes_) {
            return props_.bhp(table, 0.0, 0.0, q, thp);
        }
        const auto& a = axes_.at(table);
        return props_.bhp(table, 0.0, 0.0,
                          std::clamp(q, a.flo_min, a.flo_max),
                          std::clamp(thp, a.thp_min, a.thp_max));
    }

    double wellBhp(const double thp, const double q) const { return branchBhp(kWellTable, thp, q); }

    /// q where the IPR meets VFPINJ 1 at this THP, then BHP and rate limits.
    double wellRate(const WellProxy& w, const double p_node) const
    {
        const auto ipr = [&w](const double bhp) { return GasInjectionNetwork::ipr(w, bhp); };
        // bhp falls with rate at fixed thp in these tables, so f is decreasing and
        // bisection is safe. Search only where the table still has a solution.
        const auto f = [&](const double q) { return ipr(wellBhp(p_node, q)) - q; };
        double lo = convert::from(5000.0, cubic(meter) / day);
        double hi = w.rate_limit;
        while (hi > lo && wellBhp(p_node, hi) <= convert::from(1.0, atm)) {
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
        if (wellBhp(p_node, q) > w.bhp_limit) {
            q = std::max(ipr(w.bhp_limit), 0.0);   // BHP-limited
        }
        return std::clamp(q, 0.0, w.rate_limit);
    }

    struct Axes { double flo_min, flo_max, thp_min, thp_max; };

    std::map<int, Axes> axes_;
    std::deque<Deck> decks_;
    std::deque<VFPInjTable> tables_;
    VFPInjProperties<double> props_;
    std::vector<WellProxy> wells_;
    double terminal_ = convert::from(340.0, bars);
    double group_target_ = 0.0;
    bool clamp_to_axes_ = false;
    double stiffness_ = convert::from(6.0e4, cubic(meter) / day) / convert::from(1.0, bars);
};

// ---------------------------------------------------------------------------
// Two formulations of the same problem
//
// Eliminated: the unknowns are the node pressures alone. Every rate is recovered
// from them by an inner solve, so the residual is cheap to state and awkward to
// differentiate -- the control clamps put kinks in it.
//
// Full: node pressures, branch rates and each well's (rate, bhp) are all
// unknowns, and the eliminations become equations. Bigger and smooth within an
// active set, and the shape a reservoir/well system could absorb.
//
// Both expose size() and residual(x), so the Newton below does not know which
// one it is solving.
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
    explicit EliminatedProblem(const GasInjectionNetwork& net) : net_(net) {}

    static constexpr const char* name = "eliminated (2 unknowns)";
    int size() const { return 2; }

    State residual(const State& x) const
    {
        const auto r = net_.residual({x[0], x[1]});
        return {r[0] / kPressureScale, r[1] / kPressureScale};
    }

    /// Both formulations are started from the same pair of node pressures.
    State start(const Vec& p) const { return {p[0], p[1]}; }

    Vec pressures(const State& x) const { return {x[0], x[1]}; }

    /// Natural magnitude of unknown i; both are pressures here.
    double columnScale(const int) const { return kPressureScale; }

private:
    const GasInjectionNetwork& net_;
};

/// The same network without the eliminations.
///
///   unknowns   p(M5S) p(M5N) p(G1) p(F1)                                    4
///              q on each of the four branches                               4
///              (rate, bhp) for each of the four wells                       8
///
///   equations  branch drop      p_child - VFP_b(thp = p_parent, q_b) = 0    4
///              node balance     inflow - outflows - wells = 0               4
///              inflow perf.     q_w - ipr_w(bhp_w) = 0                      4
///              control          whichever of THP / BHP / RATE is active     4
///
/// The group target is not part of this formulation yet; it wants a multiplier
/// and an extra equation, and the comparisons below do not set one.
class FullProblem
{
public:
    explicit FullProblem(const GasInjectionNetwork& net) : net_(net) {}

    static constexpr const char* name = "full (16 unknowns)";

    enum Index {
        P_M5S, P_M5N, P_G1, P_F1,
        Q_PLATA_M5S, Q_M5S_M5N, Q_M5S_G1, Q_M5N_F1,
        Q_WELL,                          // four rates
        BHP_WELL = Q_WELL + 4,           // four bottom-hole pressures
        NUM_UNKNOWNS = BHP_WELL + 4
    };

    /// Which equation closes each well. Chosen from the current iterate, the way
    /// the well model picks a control, and held fixed while the step is taken.
    enum class Control { Thp, Bhp, Rate };

    int size() const { return NUM_UNKNOWNS; }

    State residual(const State& x) const
    {
        const auto& wells = net_.wells();
        // Wells 0,1 hang off G1 and wells 2,3 off F1.
        const std::array<int, 4> well_node{P_G1, P_G1, P_F1, P_F1};

        State r(NUM_UNKNOWNS, 0.0);

        // Branch drops. G1 and F1 carry table 9999, which is no table at all.
        r[0] = x[P_M5S] - net_.branch(kM5sTable, net_.terminal(), x[Q_PLATA_M5S]);
        r[1] = x[P_M5N] - net_.branch(kM5nTable, x[P_M5S], x[Q_M5S_M5N]);
        r[2] = x[P_G1] - x[P_M5S];
        r[3] = x[P_F1] - x[P_M5N];

        // Node balances.
        r[4] = x[Q_PLATA_M5S] - x[Q_M5S_M5N] - x[Q_M5S_G1];
        r[5] = x[Q_M5S_M5N] - x[Q_M5N_F1];
        r[6] = x[Q_M5S_G1] - x[Q_WELL + 0] - x[Q_WELL + 1];
        r[7] = x[Q_M5N_F1] - x[Q_WELL + 2] - x[Q_WELL + 3];

        for (int w = 0; w < 4; ++w) {
            const double q = x[Q_WELL + w];
            const double bhp = x[BHP_WELL + w];
            r[8 + w] = q - GasInjectionNetwork::ipr(wells[w], bhp);

            switch (controls_[w]) {
            case Control::Thp:
                r[12 + w] = bhp - net_.branch(kWellTable, x[well_node[w]], q);
                break;
            case Control::Bhp:
                r[12 + w] = bhp - wells[w].bhp_limit;
                break;
            case Control::Rate:
                r[12 + w] = q - wells[w].rate_limit;
                break;
            }
        }

        for (int i = 0; i < 4; ++i) {
            r[i] /= kPressureScale;
        }
        for (int i = 4; i < 12; ++i) {
            r[i] /= kRateScale;
        }
        for (int i = 12; i < NUM_UNKNOWNS; ++i) {
            r[i] /= (controls_[i - 12] == Control::Rate) ? kRateScale : kPressureScale;
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
        const auto& wells = net_.wells();
        bool changed = false;
        for (int w = 0; w < 4; ++w) {
            const auto& well = wells[w];
            const bool over_bhp = x[BHP_WELL + w] > well.bhp_limit;
            const bool over_rate = x[Q_WELL + w] > well.rate_limit;

            auto wanted = Control::Thp;
            if (over_bhp || over_rate) {
                wanted = GasInjectionNetwork::ipr(well, well.bhp_limit) < well.rate_limit
                    ? Control::Bhp : Control::Rate;
            }
            changed |= (wanted != controls_[w]);
            controls_[w] = wanted;
        }
        return changed;
    }

    /// Everything derived from the two starting node pressures, so the two
    /// formulations really do start from the same place.
    State start(const Vec& p) const
    {
        const auto& wells = net_.wells();
        const double q_guess = convert::from(1.0e5, cubic(meter) / day);

        State x(NUM_UNKNOWNS, 0.0);
        x[P_M5S] = x[P_G1] = p[0];
        x[P_M5N] = x[P_F1] = p[1];
        for (int w = 0; w < 4; ++w) {
            const double thp = (w < 2) ? p[0] : p[1];
            x[BHP_WELL + w] = net_.branch(kWellTable, thp, q_guess);
            x[Q_WELL + w] = std::clamp(GasInjectionNetwork::ipr(wells[w], x[BHP_WELL + w]),
                                       0.0, wells[w].rate_limit);
        }
        x[Q_M5S_G1] = x[Q_WELL + 0] + x[Q_WELL + 1];
        x[Q_M5N_F1] = x[Q_WELL + 2] + x[Q_WELL + 3];
        x[Q_M5S_M5N] = x[Q_M5N_F1];
        x[Q_PLATA_M5S] = x[Q_M5S_G1] + x[Q_M5S_M5N];
        return x;
    }

    Vec pressures(const State& x) const { return {x[P_M5S], x[P_M5N]}; }

    /// Natural magnitude of unknown i. Pressures and bhp in bar, rates in
    /// kRateScale -- this is what lets one step cap and one trust radius apply
    /// to a vector holding both.
    double columnScale(const int i) const
    {
        const bool is_rate = (i >= Q_PLATA_M5S && i < BHP_WELL);
        return is_rate ? kRateScale : kPressureScale;
    }

private:
    const GasInjectionNetwork& net_;
    std::array<Control, 4> controls_{Control::Thp, Control::Thp, Control::Thp, Control::Thp};
};

// ---------------------------------------------------------------------------
// Solvers
//
// The fixed-point methods only make sense on the eliminated form, so they take
// the network directly. The Newton takes either problem.
// ---------------------------------------------------------------------------

// Start where the simulator does: the wells' WCONINJE THP.
const Vec kStart = {convert::from(400.0, bars), convert::from(400.0, bars)};
const double kTol = 0.01;                                  // scaled: 0.01 bar
const double kMaxStep = 100.0;                             // column scales: 100 bar
constexpr int kMaxIter = 200;

struct Result
{
    bool converged = false;
    int iterations = 0;
    Vec p{};
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

Result damped(const GasInjectionNetwork& net, Vec p, const double omega)
{
    for (int it = 1; it <= kMaxIter; ++it) {
        const auto r = net.residual(p);
        if (std::max(std::abs(r[0]), std::abs(r[1])) < kTol * kPressureScale) {
            return {true, it, p};
        }
        for (int i = 0; i < 2; ++i) {
            p[i] = NodePressureUpdater<double>::damped(p[i], r[i], omega, kMaxStep * kPressureScale);
        }
    }
    return {false, kMaxIter + 1, p};
}

Result bracketing(const GasInjectionNetwork& net, Vec p, const double omega)
{
    std::array<NodePressureUpdater<double>, 2> updater;
    for (int it = 1; it <= kMaxIter; ++it) {
        const auto g = net.G(p);
        if (std::max(std::abs(g[0] - p[0]), std::abs(g[1] - p[1])) < kTol * kPressureScale) {
            return {true, it, p};
        }
        for (int i = 0; i < 2; ++i) {
            p[i] = updater[i].next(p[i], g[i], /*valid=*/true, omega, kMaxStep * kPressureScale);
        }
    }
    return {false, kMaxIter + 1, p};
}

Result anderson(const GasInjectionNetwork& net, const Vec& start, const int depth)
{
    NetworkAndersonAccelerator<double> accelerator;
    accelerator.setDepth(depth);
    State x{start[0], start[1]};
    for (int it = 1; it <= kMaxIter; ++it) {
        const Vec p{x[0], x[1]};
        const auto g = net.G(p);
        if (std::max(std::abs(g[0] - p[0]), std::abs(g[1] - p[1])) < kTol * kPressureScale) {
            return {true, it, p};
        }
        x = accelerator.next(x, {g[0], g[1]});
    }
    return {false, kMaxIter + 1, {x[0], x[1]}};
}

// --- Newton ------------------------------------------------------------------

/// Dense square system, small enough that Gaussian elimination with partial
/// pivoting is the whole story.
class Matrix
{
public:
    explicit Matrix(const int n) : n_(n), a_(n * n, 0.0) {}

    double& operator()(const int i, const int j) { return a_[i * n_ + j]; }
    double operator()(const int i, const int j) const { return a_[i * n_ + j]; }

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
// returns the point to move to. They are the whole subject of this bench: the
// Newton direction is the same in all of them.

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
Result newton(Problem problem, const Vec& start, Globalisation g = {})
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
        // The residual jumps when a control switches, and that jump is not a
        // failure to make progress. Letting a globalisation veto it stalls the
        // active set instead of resolving it.
        x = controls_moved ? x + dx : g.accept(problem, x, r, dx);
    }
    return {false, kMaxIter + 1, problem.pressures(x)};
}

} // anonymous namespace

BOOST_AUTO_TEST_SUITE(NetworkSolveBench)

// The branch tables reproduce the Eclipse 100 operating point, which is what makes
// everything below a statement about the methods and not about the model.
BOOST_AUTO_TEST_CASE(branches_match_eclipse)
{
    const GasInjectionNetwork net;
    const auto sm3d = cubic(meter) / day;
    // E100 day 31: M5S = 209.4 bar at 1.532e6 sm3/d, M5N = 204.2 bar at 5.53e5 sm3/d.
    const double m5s = net.branch(kM5sTable, convert::from(340.0, bars), convert::from(1.532e6, sm3d));
    const double m5n = net.branch(kM5nTable, convert::from(209.4, bars), convert::from(5.53e5, sm3d));
    BOOST_TEST_MESSAGE("M5S " << convert::to(m5s, bars) << " (E100 209.4), M5N "
                       << convert::to(m5n, bars) << " (E100 204.2) bar");
    BOOST_CHECK_CLOSE(convert::to(m5s, bars), 209.4, 2.0);
    BOOST_CHECK_CLOSE(convert::to(m5n, bars), 204.2, 2.0);
}

// Both formulations describe the same network, so they must land on the same point.
BOOST_AUTO_TEST_CASE(both_formulations_match_eclipse)
{
    const GasInjectionNetwork net;
    // Each formulation with the method that suits it: the eliminated residual
    // needs globalising, the full one does not.
    const auto eliminated = newton(EliminatedProblem{net}, kStart, TrustRegion{});
    const auto full = newton(FullProblem{net}, kStart, FullStep{});

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

// GCONINJE puts the wells on GRUP control once their unconstrained rates exceed the
// field target. That plateau is a large part of the real response, and it is what a
// plain dq/dbhp proxy has no way of seeing.
BOOST_AUTO_TEST_CASE(group_target_caps_the_rates)
{
    const auto sm3d = cubic(meter) / day;
    const Vec p = {convert::from(209.4, bars), convert::from(204.2, bars)};
    const auto total = [](const std::array<double, 4>& q) { return q[0] + q[1] + q[2] + q[3]; };

    GasInjectionNetwork net;
    BOOST_REQUIRE_GT(convert::to(total(net.rates(p)), sm3d), 1.0e6);

    net.setGroupTarget(1.0e6);
    BOOST_CHECK_CLOSE(convert::to(total(net.rates(p)), sm3d), 1.0e6, 1e-6);
}

namespace {
    void report(const char* name, const Result& r)
    {
        BOOST_TEST_MESSAGE(std::left << std::setw(26) << name
                           << (r.converged ? "converged in " : "FAILED after ")
                           << std::setw(4) << r.iterations << " iterations, p = ("
                           << convert::to(r.p[0], bars) << ", "
                           << convert::to(r.p[1], bars) << ") bar");
    }

    // The grid the basin tests sweep: starting node pressures across the tables' THP axis.
    std::vector<Vec> startingPoints()
    {
        std::vector<Vec> starts;
        for (int a = 60; a <= 500; a += 20) {
            for (int b = 60; b <= 500; b += 20) {
                starts.push_back({convert::from(a, bars), convert::from(b, bars)});
            }
        }
        return starts;
    }

    /// How many of the starting points a method reaches the right answer from.
    template <class Solve>
    int basin(const char* name, Solve&& solve)
    {
        const auto starts = startingPoints();
        const Vec expected = {convert::from(209.30, bars), convert::from(204.19, bars)};
        const double tol = convert::from(1.0, bars);

        int solved = 0, iterations = 0;
        for (const auto& start : starts) {
            const auto r = solve(start);
            if (r.converged && std::max(std::abs(r.p[0] - expected[0]),
                                        std::abs(r.p[1] - expected[1])) < tol) {
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

// From the one starting point the simulator actually uses.
BOOST_AUTO_TEST_CASE(method_comparison)
{
    const GasInjectionNetwork net;
    const EliminatedProblem eliminated{net};

    const auto fixed_point = damped(net, kStart, 0.1);
    const auto bracket = bracketing(net, kStart, 0.1);
    const auto acc = anderson(net, kStart, 4);
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
    const GasInjectionNetwork net;
    const EliminatedProblem problem{net};
    const auto n = static_cast<int>(startingPoints().size());

    const int bracket = basin("bracketing (shipped)",
                              [&](const Vec& p) { return bracketing(net, p, 0.1); });
    const int full_step = basin(FullStep::name,
                                [&](const Vec& p) { return newton(problem, p, FullStep{}); });
    const int capped = basin(CappedStep::name,
                             [&](const Vec& p) { return newton(problem, p, CappedStep{}); });
    const int search = basin(LineSearch::name,
                             [&](const Vec& p) { return newton(problem, p, LineSearch{}); });
    const int region = basin(TrustRegion::name,
                             [&](const Vec& p) { return newton(problem, p, TrustRegion{}); });

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
    const GasInjectionNetwork net;
    const EliminatedProblem eliminated{net};
    const FullProblem full{net};

    const int e_step = basin("eliminated, full step",
                             [&](const Vec& p) { return newton(eliminated, p, FullStep{}); });
    const int f_step = basin("full, full step",
                             [&](const Vec& p) { return newton(full, p, FullStep{}); });
    const int e_search = basin("eliminated, line search",
                               [&](const Vec& p) { return newton(eliminated, p, LineSearch{}); });
    const int f_search = basin("full, line search",
                               [&](const Vec& p) { return newton(full, p, LineSearch{}); });

    const auto n = static_cast<int>(startingPoints().size());
    // Holding the controls fixed while the step is taken removes the kinks, so the
    // full system needs no globalisation at all: a plain Newton recovers from
    // everything the globalised eliminated one does, in fewer iterations.
    BOOST_CHECK_LT(e_step, n / 10);
    BOOST_CHECK_EQUAL(f_step, n);
    BOOST_CHECK_EQUAL(e_search, n);
    BOOST_CHECK_EQUAL(f_search, n);
}

// The tables only describe a box in (rate, thp). Outside it they are zero-filled
// and the interpolation runs away, so the full system has a root there too -- at
// dq/dbhp = 1e4 a plain Newton reaches it from a third of the grid, ending with
// every well at its rate limit (4e6 sm3/d, twice the tables' flow axis) and node
// pressures of -683 and -10975 bar.
//
// Clamping the lookups to the axes removes that root, and costs more than it
// saves: the residual goes flat outside the box, so the Jacobian there is
// singular in the rates and the Newton has nothing to descend. Neither setting is
// good, which is the point -- the table limits want to be bounds on the unknowns,
// enforced in the active set alongside the well controls, not a flattening of the
// residual. The bracketing method never meets this because its inner bisection
// cannot leave the box in the first place.
BOOST_AUTO_TEST_CASE(table_bounds_want_to_be_constraints)
{
    const auto n = static_cast<int>(startingPoints().size());

    GasInjectionNetwork loose;
    loose.setStiffness(1.0e4);
    const int unclamped = basin("unclamped", [&](const Vec& p) {
        return newton(FullProblem{loose}, p, FullStep{});
    });

    GasInjectionNetwork clamped;
    clamped.setStiffness(1.0e4);
    clamped.setClampToAxes(true);
    const int with_clamp = basin("clamped to axes", [&](const Vec& p) {
        return newton(FullProblem{clamped}, p, FullStep{});
    });

    // Both leave part of the grid unsolved, for opposite reasons.
    BOOST_CHECK_LT(unclamped, n);
    BOOST_CHECK_LT(with_clamp, unclamped);

    // The bracketing method is indifferent: it cannot leave the box either way.
    BOOST_CHECK_EQUAL(basin("bracketing, clamped",
                            [&](const Vec& p) { return bracketing(clamped, p, 0.1); }), n);
}

// How the formulations degrade as the wells stiffen. dq/dbhp sets the loop gain.
// Measured over the whole grid of starts, because a single start says too little.
BOOST_AUTO_TEST_CASE(stiffness_sweep)
{
    const auto n = static_cast<int>(startingPoints().size());
    for (const double stiffness : {1.0e4, 6.0e4, 3.0e5, 1.0e6}) {
        GasInjectionNetwork net;
        net.setStiffness(stiffness);
        BOOST_TEST_MESSAGE("dq/dbhp = " << stiffness << " sm3/d/bar");

        const int bracket = basin("  bracketing (shipped)",
                                  [&](const Vec& p) { return bracketing(net, p, 0.1); });
        const int eliminated = basin("  eliminated, trust region",
                                     [&](const Vec& p) {
                                         return newton(EliminatedProblem{net}, p, TrustRegion{});
                                     });
        const int full = basin("  full, plain newton",
                               [&](const Vec& p) {
                                   return newton(FullProblem{net}, p, FullStep{});
                               });
        const int full_ls = basin("  full, line search",
                                  [&](const Vec& p) {
                                      return newton(FullProblem{net}, p, LineSearch{});
                                  });
        BOOST_CHECK_EQUAL(bracket, n);
        BOOST_CHECK_GE(eliminated, n - 1);
        // The full system is uniformly better except at the softest wells, where
        // the out-of-table root of table_bounds_want_to_be_constraints catches it.
        BOOST_CHECK_GE(std::max(full, full_ls), (stiffness > 1.0e4) ? n : 7 * n / 10);
    }
}

BOOST_AUTO_TEST_SUITE_END()
