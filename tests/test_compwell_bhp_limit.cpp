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

#define BOOST_TEST_MODULE CompWellBhpLimit
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

// A 30-cell row with a gas injector at 150 bar and a producer at 50 bar.
// After four report steps, the appended schedule changes one BHP limit so
// that the corresponding well initially cannot flow.
std::string deck(const std::string& schedule)
{
    return R"(RUNSPEC
METRIC
TABDIMS
/
OIL
GAS
DIMENS
30 1 1 /
COMPS
3 /
FULLIMP
WELLDIMS
2 1 1 2 /
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
SGOF
0.0 0.00 1.00 0.0
0.2 0.04 0.64 0.0
0.4 0.16 0.36 0.0
0.6 0.36 0.16 0.0
0.8 0.64 0.04 0.0
1.0 1.00 0.00 0.0
/
STCOND
15.0 1.0 /

SOLUTION
PRESSURE
30*75 /
SGAS
30*1 /
TEMPI
30*150 /
ZMF
30*0.1
30*0.3
30*0.6
/

SUMMARY
WBHP
/
WOPR
/
WGPR
/
WOIR
/
WGIR
/
BPR
1 1 1 /
30 1 1 /
/

SCHEDULE
WELSPECS
INJ  FIELD  1 1 1* GAS /
PROD FIELD 30 1 1* GAS /
/
COMPDAT
INJ   1 1 1 1 OPEN 2* 0.0151 /
PROD 30 1 1 1 OPEN 2* 0.0151 /
/
WELLSTRE
ISTR 1.0 0.0 0.0 /
/
WINJGAS
INJ STREAM ISTR /
/
WCONINJE
INJ GAS OPEN BHP 100000 1* 150 /
/
WCONPROD
PROD OPEN BHP 2* 100000 2* 50 /
/
TSTEP
0.01 0.1 0.4 1 /
)" + schedule + "END\n";
}

// Index of the first report step after the appended schedule takes effect.
constexpr std::size_t first_step = 4;

const std::vector<std::string> rates = {"WOPR", "WGPR", "WOIR", "WGIR"};

// Run the deck and return time, well BHPs and rates, and connection-cell
// pressures under their summary keys.
std::map<std::string, std::vector<float>>
run(const std::filesystem::path& dir, const std::string& name, const std::string& deck_text)
{
    const auto deck_file = dir / (name + ".DATA");
    std::ofstream(deck_file) << deck_text;

    std::vector<std::string> args = {"test_compwell_bhp_limit", deck_file.string(),
                                     "--output-dir=" + dir.string()};
    std::vector<char*> argv;
    std::ranges::transform(args, std::back_inserter(argv),
                           [](std::string& arg) { return arg.data(); });
    argv.push_back(nullptr);
    Opm::Parameters::reset();
    // The deck has no water phase.
    const int status = Opm::Main(static_cast<int>(args.size()), argv.data(), /*ownMPI=*/false)
                           .runStatic<Opm::Properties::TTag::FlowCompProblem<3, false>>();
    BOOST_REQUIRE_EQUAL(status, EXIT_SUCCESS);

    const Opm::EclIO::ESmry smry(dir / (name + ".SMSPEC"));
    smry.loadData();
    std::map<std::string, std::vector<float>> result;
    for (const std::string key : {"TIME", "BPR:1,1,1", "BPR:30,1,1"}) {
        result[key] = smry.get(key);
    }
    for (const std::string well : {"INJ", "PROD"}) {
        result["WBHP:" + well] = smry.get("WBHP:" + well);
        for (const auto& rate : rates) {
            result[rate + ":" + well] = smry.get(rate + ":" + well);
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

// Raise the producer's minimum BHP above its cell pressure. It should stop
// rather than inject, then reopen when the injector raises the pressure enough.
BOOST_AUTO_TEST_CASE(ProducerStopsAndReopensAtItsBhpLimit)
{
    // Give each run a separate output directory.
    const auto dir = std::filesystem::temp_directory_path()
        / Opm::unique_path("compwell_bhp_limit_producer-%%%%%%%%");
    std::filesystem::create_directories(dir);

    const auto smry = run(dir, "PRODUCER", deck(R"(WCONPROD
PROD OPEN BHP 2* 100000 2* 120 /
/
TSTEP
8*0.25 4*1 /
)"));
    const auto& time = smry.at("TIME");
    const auto& bhp = smry.at("WBHP:PROD");
    const auto& cell_pressure = smry.at("BPR:30,1,1");

    // Find the first report step at which the producer reaches its BHP limit again.
    auto reopened = first_step;
    while (reopened < bhp.size() && std::abs(bhp[reopened] - 120.0f) > 1.0e-3f) {
        ++reopened;
    }
    BOOST_CHECK_GT(reopened, first_step);
    BOOST_CHECK_LT(reopened, bhp.size());

    for (auto step = first_step; step < reopened; ++step) {
        BOOST_TEST_CONTEXT("stopped producer at day " << time[step]) {
            for (const auto& rate : rates) {
                BOOST_CHECK_SMALL(smry.at(rate + ":PROD")[step], 1.0e-6f);
            }
            BOOST_CHECK_CLOSE(bhp[step], cell_pressure[step], 0.05f);
            BOOST_CHECK_GT(bhp[step], bhp[step - 1]);
        }
    }
    for (auto step = reopened; step < bhp.size(); ++step) {
        BOOST_TEST_CONTEXT("reopened producer at day " << time[step]) {
            BOOST_CHECK_GT(smry.at("WOPR:PROD")[step], 0.1f);
            BOOST_CHECK_GT(smry.at("WGPR:PROD")[step], 1.0f);
            BOOST_CHECK_SMALL(smry.at("WOIR:PROD")[step], 1.0e-6f);
            BOOST_CHECK_SMALL(smry.at("WGIR:PROD")[step], 1.0e-6f);
        }
    }

    std::filesystem::remove_all(dir);
}

// Lower the injector's maximum BHP below its cell pressure. It should stop
// rather than produce. Its BHP falls with the cell pressure as the producer
// depletes the row, but the injector remains stopped.
BOOST_AUTO_TEST_CASE(InjectorStopsAtItsBhpLimit)
{
    // Give each run a separate output directory.
    const auto dir = std::filesystem::temp_directory_path()
        / Opm::unique_path("compwell_bhp_limit_injector-%%%%%%%%");
    std::filesystem::create_directories(dir);

    const auto smry = run(dir, "INJECTOR", deck(R"(WCONINJE
INJ GAS OPEN BHP 100000 1* 40 /
/
TSTEP
8*0.25 4*1 /
)"));
    const auto& time = smry.at("TIME");
    const auto& bhp = smry.at("WBHP:INJ");
    const auto& cell_pressure = smry.at("BPR:1,1,1");

    BOOST_CHECK_GT(time.size(), first_step);
    for (auto step = first_step; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("stopped injector at day " << time[step]) {
            for (const auto& rate : rates) {
                BOOST_CHECK_SMALL(smry.at(rate + ":INJ")[step], 1.0e-6f);
            }
            BOOST_CHECK_CLOSE(bhp[step], cell_pressure[step], 0.05f);
            BOOST_CHECK_LT(bhp[step], bhp[step - 1]);
            BOOST_CHECK_GT(smry.at("WGPR:PROD")[step], 1.0f);
        }
    }

    std::filesystem::remove_all(dir);
}
