/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE WellDensitySegmentedTest

#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/core/wells.h>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/Units.hpp>

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <memory>

using namespace Opm;

BOOST_AUTO_TEST_CASE(TestPressureDeltas)
{
    // Simple water injector.
    const int np = 3;
    const int nperf = 10;
    const double ref_depth = 0.0;
    const double comp_frac_w[np] = { 1.0, 0.0, 0.0 };
    const double comp_frac_o[np] = { 0.0, 1.0, 0.0 };
    const int cells[nperf/2] = { 0, 1, 2, 3, 4 };
    const double WI[nperf/2] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    std::shared_ptr<Wells> wells(create_wells(np, 2, nperf), destroy_wells);
    BOOST_REQUIRE(wells);
    int ok = add_well(INJECTOR, ref_depth, nperf/2, comp_frac_w, cells, WI, "INJ", wells.get());
    BOOST_REQUIRE(ok);
    ok = add_well(PRODUCER, ref_depth, nperf/2, comp_frac_o, cells, WI, "PROD", wells.get());
    BOOST_REQUIRE(ok);
    std::vector<double> rates = { 1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0 };
    WellStateFullyImplicitBlackoil wellstate;
    wellstate.perfPhaseRates() = rates;
    PhaseUsage pu;
    pu.num_phases = 3;
    pu.phase_used[0] = true;
    pu.phase_used[1] = true;
    pu.phase_used[2] = true;
    pu.phase_pos[0] = 0;
    pu.phase_pos[1] = 1;
    pu.phase_pos[2] = 2;
    const std::vector<double> b_perf = { 2.0, 3.0, 100,
                                         2.1, 3.3, 110,
                                         2.2, 3.6, 120,
                                         2.3, 4.0, 130,
                                         2.4, 4.5, 140,
                                         2.0, 3.0, 100,
                                         2.1, 3.3, 110,
                                         2.2, 3.6, 120,
                                         2.3, 4.0, 130,
                                         2.4, 4.5, 140 };
    const std::vector<double> rsmax_perf = { 50, 50, 50, 50, 50, 50, 50, 50, 50, 50 };
    const std::vector<double> rvmax_perf = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
    const std::vector<double> z_perf = { 10, 30, 50, 70, 90, 10, 30, 50, 70, 90 };
    const std::vector<double> surf_dens = { 1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0,
                                          1000.0, 800.0, 10.0};
    const double gravity = Opm::unit::gravity;

    std::vector<double> cd =
            WellDensitySegmented::computeConnectionDensities(
                    *wells, wellstate, pu,
                    b_perf, rsmax_perf, rvmax_perf, surf_dens);
    std::vector<double> dp =
            WellDensitySegmented::computeConnectionPressureDelta(
                    *wells, z_perf, cd, gravity);

    const std::vector<double> answer = { 20e3*gravity, 62e3*gravity, 106e3*gravity, 152e3*gravity, 200e3*gravity, 
                                         20e3*gravity, 62e3*gravity, 106e3*gravity, 152e3*gravity, 200e3*gravity };
    BOOST_REQUIRE_EQUAL(dp.size(), answer.size());
    // for (auto p : dp) { std::cout << p << std::endl; }
    for (size_t i = 0; i < dp.size(); ++i) {
        BOOST_CHECK_CLOSE(dp[i], answer[i], 1e-8);
    }
}
