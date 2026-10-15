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

#define BOOST_TEST_MODULE CompWellZeroRateTargets
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <flowexperimental/comp/flow_comp.hpp>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/common/utility/FileSystem.hpp>
#include <opm/io/eclipse/ESmry.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

namespace {

// A row of single-connection wells. INJ injects gas at a fixed bhp, which
// raises the pressure around the other wells, all of which have a zero rate
// target. The water around PROD_WI is immobile. A zero target holds the
// surface rate at zero; wellbore storage still lets the connections flow
// while the pressure rises.
constexpr auto deck = R"(RUNSPEC
METRIC
TABDIMS
/
OIL
GAS
WATER
DIMENS
30 1 1 /
COMPS
3 /
FULLIMP
WELLDIMS
8 1 1 8 /
START
1 'JAN' 2016 /
UNIFOUT

GRID
DX
30*10 /
DY
30*1 /
DZ
30*1 /
TOPS
30*0 /
PERMX
30*100 /
PERMY
30*100 /
PERMZ
30*100 /
PORO
30*0.1 /

PROPS
CNAMES
CO2
METHANE
DECANE
/
ROCK
68.9476 0 /
EOS
PR /
BIC
0
0
0
/
ACF
0.22394 0.01142 0.4884 /
PCRIT
73.773 45.992 21.03 /
TCRIT
304.128 190.564 617.7 /
MW
44.00 16.04 142.28 /
VCRIT
0.09412 0.09863 0.60980 /
SGFN
0.0  0.00 0.0
0.1  0.00 0.0
0.2  0.05 0.0
0.4  0.10 0.0
0.6  0.40 0.0
0.8  0.60 0.0
0.99 0.70 0.0
/
SWFN
0.01 0.00 0.0
0.20 0.00 0.0
0.40 0.01 0.0
0.60 0.20 0.0
0.80 0.80 0.0
1.00 1.00 0.0
/
SOF3
0.0  0.00 0.0
0.1  0.00 0.0
0.2  0.05 0.05
0.4  0.10 0.10
0.6  0.40 0.40
0.8  0.60 0.60
0.99 0.70 0.70
/
STCOND
15.0 1.0 /
PVTW
75.0 1.03 4.0E-5 0.3 0.0 /
DENSITY
800.0 1000.0 1.0 /

SOLUTION
PRESSURE
30*75 /
SGAS
30*0.5 /
-- krw is zero up to SW = 0.2, so the water around PROD_WI cannot flow
SWAT
27*0.5 3*0.2 /
TEMPI
30*150 /
ZMF
30*0.1
30*0.3
30*0.6 /

SUMMARY
WOPR
/
WWPR
/
WGPR
/
WOIR
/
WWIR
/
WGIR
/
WBHP
/
WMCTL
/

SCHEDULE
WELSPECS
INJ     FIELD  1 1 1* GAS /
INJ_G   FIELD  5 1 1* GAS /
INJ_W   FIELD  9 1 1* WATER /
PROD_O  FIELD 13 1 1* OIL /
PROD_G  FIELD 17 1 1* GAS /
PROD_L  FIELD 21 1 1* OIL /
PROD_W  FIELD 25 1 1* WATER /
PROD_WI FIELD 30 1 1* WATER /
/
COMPDAT
INJ      1 1 1 1 OPEN 2* 0.0151 /
INJ_G    5 1 1 1 OPEN 2* 0.0151 /
INJ_W    9 1 1 1 OPEN 2* 0.0151 /
PROD_O  13 1 1 1 OPEN 2* 0.0151 /
PROD_G  17 1 1 1 OPEN 2* 0.0151 /
PROD_L  21 1 1 1 OPEN 2* 0.0151 /
PROD_W  25 1 1 1 OPEN 2* 0.0151 /
PROD_WI 30 1 1 1 OPEN 2* 0.0151 /
/
WELLSTRE
ISTR 1.0 0.0 0.0 /
/
WINJGAS
INJ   STREAM ISTR /
INJ_G STREAM ISTR /
/
WCONINJE
INJ   GAS   OPEN BHP  100000 1* 100 /
INJ_G GAS   OPEN RATE 0      1* 150 /
INJ_W WATER OPEN RATE 0      1* 150 /
/
WCONPROD
PROD_O  OPEN ORAT 0 4* 50 /
PROD_G  OPEN GRAT 2* 0 2* 50 /
PROD_L  OPEN LRAT 3* 0 1* 50 /
PROD_W  OPEN WRAT 1* 0 3* 50 /
PROD_WI OPEN WRAT 1* 0 3* 50 /
/
TSTEP
0.01 0.1 1 5 10 /
END
)";

struct MpiFixture
{
    MpiFixture()
    {
        auto& suite = boost::unit_test::framework::master_test_suite();
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(suite.argc, suite.argv);
#else
        Dune::MPIHelper::instance(suite.argc, suite.argv);
#endif
        Opm::FlowGenericVanguard::setCommunication(
            std::make_unique<Opm::Parallel::Communication>());
    }
};

} // namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);

BOOST_AUTO_TEST_CASE(ZeroSurfaceRatesUnderZeroTargets)
{
    // a directory of its own, so that concurrent runs leave each other's files alone
    const auto dir = std::filesystem::temp_directory_path()
        / Opm::unique_path("compwell_zero_rate_targets-%%%%%%%%");
    std::filesystem::create_directories(dir);
    const auto deck_file = dir / "ZERO_RATE_TARGETS.DATA";
    std::ofstream(deck_file) << deck;

    std::vector<std::string> args = {"test_compwell_zero_rate_targets", deck_file.string(),
                                     "--output-dir=" + dir.string()};
    std::vector<char*> argv;
    std::ranges::transform(args, std::back_inserter(argv),
                           [](std::string& arg) { return arg.data(); });
    argv.push_back(nullptr);
    Opm::Parameters::reset();
    const int status = Opm::Main(static_cast<int>(args.size()), argv.data(), /*ownMPI=*/false)
                           .runStatic<Opm::Properties::TTag::FlowCompProblem<3, true>>();
    BOOST_REQUIRE_EQUAL(status, EXIT_SUCCESS);

    const Opm::EclIO::ESmry smry(dir / "ZERO_RATE_TARGETS.SMSPEC");
    smry.loadData();

    // each well with a zero rate target and the WMCTL value of its control
    const std::vector<std::pair<std::string, float>> wells = {
        {"INJ_G", 3.0f}, {"INJ_W", 2.0f}, {"PROD_O", 1.0f}, {"PROD_G", 3.0f},
        {"PROD_L", 4.0f}, {"PROD_W", 2.0f}, {"PROD_WI", 2.0f},
    };
    for (const auto& [well, control] : wells) {
        for (const std::string rate : {"WOPR", "WWPR", "WGPR", "WOIR", "WWIR", "WGIR"}) {
            const auto key = rate + ":" + well;
            const auto magnitudes = smry.get(key) | std::views::transform([](const float value) {
                return std::abs(value);
            });
            BOOST_TEST_CONTEXT(key) {
                BOOST_CHECK_SMALL(std::ranges::max(magnitudes), 1.0e-6f);
            }
        }
        const auto [lowest, highest] = std::ranges::minmax(smry.get("WMCTL:" + well));
        BOOST_TEST_CONTEXT("WMCTL:" << well) {
            BOOST_CHECK_EQUAL(lowest, control);
            BOOST_CHECK_EQUAL(highest, control);
        }
        // INJ has raised the pressure around the well from 75 bar.
        BOOST_TEST_CONTEXT("WBHP:" << well) {
            BOOST_CHECK_GT(smry.get("WBHP:" + well).back(), 85.0f);
        }
    }

    std::filesystem::remove_all(dir);
}
