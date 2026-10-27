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

#define BOOST_TEST_MODULE CompWellPvtRegion
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
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <vector>

namespace {

// A row of cells with a gas injector at one end, a water injector on a rate
// target in the middle and a producer with a water rate limit at the other
// end. With two PVT regions, every cell lies in region 2, which holds the
// water of the one-region deck, and region 1 holds a different water that no
// cell uses. Both decks thus describe the same reservoir.
std::string deck(const bool two_regions)
{
    const std::string tabdims = two_regions ? "1* 2 /" : "/";
    const std::string water = two_regions ? R"(
ROCK
68.9476 0 /
68.9476 0 /
PVTW
75.0 1.01 2.5E-5 0.9 0.0 /
75.0 1.03 4.0E-5 0.3 0.0 /
DENSITY
800.0 1150.0 1.0 /
800.0 1000.0 1.0 /

REGIONS
PVTNUM
30*2 /
)" : R"(
ROCK
68.9476 0 /
PVTW
75.0 1.03 4.0E-5 0.3 0.0 /
DENSITY
800.0 1000.0 1.0 /
)";

    return R"(RUNSPEC
METRIC
TABDIMS
)" + tabdims + R"(
OIL
GAS
WATER
DIMENS
30 1 1 /
COMPS
3 /
FULLIMP
WELLDIMS
3 1 1 3 /
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
)" + water + R"(
SOLUTION
PRESSURE
30*75 /
SGAS
30*0.5 /
SWAT
30*0.5 /
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
WWIR
/
WGIR
/
WBHP
/

SCHEDULE
WELSPECS
INJ   FIELD  1 1 1* GAS /
INJ_W FIELD 15 1 1* WATER /
PROD  FIELD 30 1 1* GAS /
/
COMPDAT
INJ    1 1 1 1 OPEN 2* 0.0151 /
INJ_W 15 1 1 1 OPEN 2* 0.0151 /
PROD  30 1 1 1 OPEN 2* 0.0151 /
/
WELLSTRE
ISTR 1.0 0.0 0.0 /
/
WINJGAS
INJ STREAM ISTR /
/
WCONINJE
INJ   GAS   OPEN BHP  100000 1* 150 /
INJ_W WATER OPEN RATE 0.5    1* 150 /
/
WCONPROD
PROD OPEN BHP 1* 0.1 100000 2* 50 /
/
TSTEP
0.01 0.02 0.04 0.08 0.15 0.3 0.4 9*1.0 /
END
)";
}

const std::vector<std::string> wells = {"INJ", "INJ_W", "PROD"};
const std::vector<std::string> vectors = {"WBHP", "WOPR", "WWPR", "WGPR", "WWIR", "WGIR"};

// Runs the deck and returns the summary vectors of every well, keyed by
// VECTOR:WELL.
std::map<std::string, std::vector<float>>
run(const std::filesystem::path& dir, const std::string& name, const std::string& deck_text)
{
    const auto deck_file = dir / (name + ".DATA");
    std::ofstream(deck_file) << deck_text;

    std::vector<std::string> args = {"test_compwell_pvt_region", deck_file.string(),
                                     "--output-dir=" + dir.string()};
    std::vector<char*> argv;
    std::ranges::transform(args, std::back_inserter(argv),
                           [](std::string& arg) { return arg.data(); });
    argv.push_back(nullptr);
    Opm::Parameters::reset();
    const int status = Opm::Main(static_cast<int>(args.size()), argv.data(), /*ownMPI=*/false)
                           .runStatic<Opm::Properties::TTag::FlowCompProblem<3, true>>();
    BOOST_REQUIRE_EQUAL(status, EXIT_SUCCESS);

    const Opm::EclIO::ESmry smry(dir / (name + ".SMSPEC"));
    smry.loadData();
    std::map<std::string, std::vector<float>> result;
    for (const auto& well : wells) {
        for (const auto& vector : vectors) {
            const auto key = vector + ":" + well;
            result[key] = smry.get(key);
        }
    }
    return result;
}

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

// The wells take their water properties from the PVT region of their
// connection cell, so the unused region of the two-region deck changes nothing.
BOOST_AUTO_TEST_CASE(WaterOfTheConnectionCellPvtRegion)
{
    // a directory of its own, so that concurrent runs leave each other's files alone
    const auto dir = std::filesystem::temp_directory_path()
        / Opm::unique_path("compwell_pvt_region-%%%%%%%%");
    std::filesystem::create_directories(dir);

    const auto one_region = run(dir, "ONE_REGION", deck(false));
    const auto two_regions = run(dir, "TWO_REGIONS", deck(true));

    // water flows through both the injector and the producer
    BOOST_CHECK_GT(std::ranges::max(one_region.at("WWIR:INJ_W")), 0.4f);
    BOOST_CHECK_GT(std::ranges::max(one_region.at("WWPR:PROD")), 0.05f);

    for (const auto& [key, expected] : one_region) {
        const auto& actual = two_regions.at(key);
        BOOST_TEST_CONTEXT(key) {
            BOOST_REQUIRE_EQUAL(actual.size(), expected.size());
            for (std::size_t step = 0; step < expected.size(); ++step) {
                BOOST_TEST_CONTEXT("step " << step) {
                    BOOST_CHECK_LE(std::abs(actual[step] - expected[step]),
                                   1.0e-5f * std::max(1.0f, std::abs(expected[step])));
                }
            }
        }
    }

    std::filesystem::remove_all(dir);
}
