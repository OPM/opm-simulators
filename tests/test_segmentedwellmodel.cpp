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

#define BOOST_TEST_MODULE SegmentedWellModelTest

#include <opm/autodiff/SegmentedWellModel.hpp>
#include <opm/core/wells.h>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/Units.hpp>

#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace Opm;

BOOST_AUTO_TEST_CASE(TestPressureDeltas)
{
    // Simple water injector.
    const int np = 3;
    const int nperf = 4;
    const double ref_depth = 1000.0;
    const double comp_frac[np] = { 1.0, 0.0, 0.0 };
    const int cells[nperf] = { 0, 1, 2, 3 };
    const double WI[nperf] = { 1.0, 1.0, 1.0, 1.0 };
    const char *name = "INJECTOR";
    Wells* wells = create_wells(np, 1, nperf);
    BOOST_REQUIRE(wells != NULL);
    const int ok = add_well(INJECTOR, ref_depth, nperf, comp_frac, cells, WI, name, wells);
    BOOST_REQUIRE(ok);
    std::vector<double> rates = { 1.0, 0.0, 0.0,
                                  1.0, 0.0, 0.0,
                                  1.5, 0.0, 0.0,
                                  1.5, 0.0, 0.0 };
    WellState wellstate;
    wellstate.perfRates() = rates;
    PhaseUsage pu;
    pu.num_phases = 3;
    pu.phase_used[0] = true;
    pu.phase_used[1] = true;
    pu.phase_used[2] = true;
    pu.phase_pos[0] = 0;
    pu.phase_pos[1] = 1;
    pu.phase_pos[2] = 2;
    const std::vector<double> b_perf = { 2.0, 3.0, 100.0,
                                         2.0, 3.5, 110.0,
                                         2.0, 4.0, 120.0,
                                         2.0, 4.5, 130.0 };
    const std::vector<double> rsmax_perf;
    const std::vector<double> rvmax_perf;
    const std::vector<double> z_perf = { 1100.0, 1200.0, 1300.0, 1400.0 };
    const std::vector<double> surf_dens = { 1000.0, 800.0, 10.0 };
    const double gravity = Opm::unit::gravity;
    const std::vector<double> dp = SegmentedWellModel::computeConnectionPressureDelta(*wells, wellstate, pu, b_perf, rsmax_perf, rvmax_perf, z_perf, surf_dens, gravity);
    const std::vector<double> answer = { 2e5*gravity, 4e5*gravity, 6e5*gravity, 8e5*gravity };
    BOOST_REQUIRE_EQUAL(dp.size(), answer.size());
    for (size_t i = 0; i < dp.size(); ++i) {
        BOOST_CHECK_CLOSE(dp[i], answer[i], 1e-8);
    }
}
