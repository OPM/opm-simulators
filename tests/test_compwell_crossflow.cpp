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

#define BOOST_TEST_MODULE CompWellCrossflow
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

// Two layers without flow between them, at 100 and 60 bar, and a producer
// at 70 bar connected to both. With crossflow, the producer injects into the
// lower-pressure layer through its connection there.
std::string deck(const std::string& crossflow)
{
    return R"(RUNSPEC
METRIC
TABDIMS
/
OIL
GAS
DIMENS
1 1 2 /
COMPS
3 /
FULLIMP
WELLDIMS
1 2 1 1 /
START
1 'JAN' 2016 /
UNIFOUT

GRID
DX
2*1000 /
DY
2*1 /
DZ
2*1 /
TOPS
0 /
PERMX
2*100 /
PERMY
2*100 /
PERMZ
2*0 /
PORO
2*0.1 /

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
100 60 /
SGAS
2*1 /
TEMPI
2*150 /
ZMF
2*0.1
2*0.3
2*0.6
/

SUMMARY
WGPR
/
BPR
1 1 1 /
1 1 2 /
/

SCHEDULE
WELSPECS
PROD FIELD 1 1 1* GAS 3* )" + crossflow + R"( /
/
COMPDAT
PROD 1 1 1 2 OPEN 2* 0.0151 /
/
WCONPROD
PROD OPEN BHP 2* 100000 2* 70 /
/
TSTEP
5*0.1 /
END
)";
}

// Run the deck and return the producer's gas rate and the layer pressures
// under their summary keys.
std::map<std::string, std::vector<float>>
run(const std::filesystem::path& dir, const std::string& name, const std::string& deck_text)
{
    const auto deck_file = dir / (name + ".DATA");
    std::ofstream(deck_file) << deck_text;

    std::vector<std::string> args = {"test_compwell_crossflow", deck_file.string(),
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
    for (const std::string key : {"TIME", "WGPR:PROD", "BPR:1,1,1", "BPR:1,1,2"}) {
        result[key] = smry.get(key);
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

// Without crossflow, the producer's connection in the lower-pressure layer is
// closed, and that layer keeps its pressure, while the producer still produces
// from the other layer. With crossflow, the layer takes up fluid.
BOOST_AUTO_TEST_CASE(ProducerConnectionDoesNotInjectWithoutCrossflow)
{
    // Give each run a separate output directory.
    const auto dir = std::filesystem::temp_directory_path()
        / Opm::unique_path("compwell_crossflow-%%%%%%%%");
    std::filesystem::create_directories(dir);

    const auto without = run(dir, "WITHOUT_CROSSFLOW", deck("NO"));
    const auto& time = without.at("TIME");
    BOOST_CHECK_GT(time.size(), 0u);
    for (std::size_t step = 0; step < time.size(); ++step) {
        BOOST_TEST_CONTEXT("without crossflow at day " << time[step]) {
            BOOST_CHECK_SMALL(without.at("BPR:1,1,2")[step] - 60.0f, 1.0e-3f);
            BOOST_CHECK_LT(without.at("BPR:1,1,1")[step], 100.0f);
            BOOST_CHECK_GT(without.at("WGPR:PROD")[step], 1.0f);
        }
    }

    const auto with = run(dir, "WITH_CROSSFLOW", deck("YES"));
    BOOST_CHECK_GT(with.at("BPR:1,1,2").back(), 61.0f);

    std::filesystem::remove_all(dir);
}
