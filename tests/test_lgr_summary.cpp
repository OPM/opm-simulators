// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \brief Test the summary of a serial run with an LGR
 *
 * LGR_SUMMARY.DATA refines the centre cell (2,2,1) of a 3x3x1 grid into 3x3x1
 * and produces 100 STB/day from the centre cell of the LGR.  The run writes
 * its summary file, which must hold the field, well, connection and block
 * values, the LGR ones included.
 */
#include <config.h>

#define BOOST_TEST_MODULE LgrSummary
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/io/eclipse/EclFile.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>
#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/TTagFlowProblemTPFA.hpp>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <cstddef>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace Opm::Properties::TTag {

struct TestTypeTag
{ using InheritsFrom = std::tuple<FlowProblemTPFA>; };

}

namespace Opm::Properties {

// Disable convective mixing
template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::TestTypeTag>
{ static constexpr bool value = false; };

// Disable diffusion
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::TestTypeTag>
{ static constexpr bool value = false; };

}

namespace Opm {

// Runs a deck with the black-oil simulator.
class LgrSummaryMain : public Main
{
public:
    LgrSummaryMain(int argc, char** argv)
        // ownMPI=false: MPI is initialised once, in main() below.
        : Main{argc, argv, /*ownMPI=*/false}
    {}

    int run()
    {
        int exitCode = EXIT_SUCCESS;
        if (this->initialize_<Properties::TTag::FlowEarlyBird>(exitCode, /*keep_keywords=*/false)) {
            this->setupVanguard();
            FlowMain<Properties::TTag::TestTypeTag> flowMain(this->argc_, this->argv_,
                                                             this->outputCout_, this->outputFiles_);
            exitCode = flowMain.execute();
        }
        return exitCode;
    }
};

} // namespace Opm

namespace {

// Last value of each summary vector, keyed by (keyword, well, cell, LGR).
using Key = std::tuple<std::string, std::string, int, std::string>;

std::map<Key, float> lastValues(const std::filesystem::path& smspecPath,
                                const std::filesystem::path& unsmryPath)
{
    Opm::EclIO::EclFile smspec(smspecPath.string());
    const auto& keywords = smspec.get<std::string>("KEYWORDS");
    const auto& names = smspec.get<std::string>("WGNAMES");
    const auto& nums = smspec.get<int>("NUMS");
    const auto lgrs = smspec.hasKey("LGRS") ? smspec.get<std::string>("LGRS") : std::vector<std::string>{};

    Opm::EclIO::EclFile unsmry(unsmryPath.string());
    std::size_t last = 0;
    const auto arrays = unsmry.getList();
    for (std::size_t i = 0; i < arrays.size(); ++i) {
        if (std::get<0>(arrays[i]) == "PARAMS") {
            last = i;
        }
    }
    const auto& params = unsmry.get<float>(last);

    std::map<Key, float> values;
    for (std::size_t i = 0; i < keywords.size(); ++i) {
        const auto lgr = (i < lgrs.size()) ? lgrs[i] : std::string{};
        values[{keywords[i], names[i], nums[i], lgr}] = params[i];
    }
    return values;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(SerialRunWithAnLgr)
{
    const auto outputDir = std::filesystem::temp_directory_path() / "test_lgr_summary";
    std::filesystem::remove_all(outputDir);

    // CARFIN is parsed at low strictness, as in the other runs with local grids.
    const std::string outputArg = "--output-dir=" + outputDir.string();
    const char* argv[] = {"test_lgr_summary", "LGR_SUMMARY.DATA", "--parsing-strictness=low",
                          outputArg.c_str(), nullptr};
    Opm::LgrSummaryMain main(4, const_cast<char**>(argv));
    BOOST_REQUIRE_EQUAL(main.run(), EXIT_SUCCESS);

    auto values = lastValues(outputDir / "LGR_SUMMARY.SMSPEC", outputDir / "LGR_SUMMARY.UNSMRY");
    const auto value = [&values](const std::string& keyword, const std::string& name,
                                 const int num, const std::string& lgr) -> double
    {
        const auto it = values.find({keyword, name, num, lgr});
        BOOST_REQUIRE_MESSAGE(it != values.end(), "no summary vector " + keyword);
        return it->second;
    };
    const std::string none = ":+:+:+:+";

    // Field and well vectors: the producer holds its oil rate of 100 STB/day for 20 days.
    BOOST_CHECK_CLOSE(value("TIME", none, 0, ""), 20.0, 1.0e-4);
    BOOST_CHECK_CLOSE(value("FOPR", none, 0, ""), 100.0, 1.0e-4);
    BOOST_CHECK_CLOSE(value("FOPT", none, 0, ""), 2000.0, 1.0e-4);
    BOOST_CHECK_CLOSE(value("WOPR", "PROD", 0, ""), 100.0, 1.0e-4);
    BOOST_CHECK_GT(value("FPR", none, 0, ""), 0.0);
    BOOST_CHECK_GT(value("WBHP", "PROD", 0, ""), 0.0);

    // FOE needs the initial in-place volumes.
    BOOST_CHECK_GT(value("FOE", none, 0, ""), 0.0);

    // The LGR vectors of the well and its single connection equal the well's own.
    BOOST_CHECK_CLOSE(value("LWOPR", "PROD", 0, "LGR1"), value("WOPR", "PROD", 0, ""), 1.0e-4);
    BOOST_CHECK_CLOSE(value("LWBHP", "PROD", 0, "LGR1"), value("WBHP", "PROD", 0, ""), 1.0e-4);
    BOOST_CHECK_CLOSE(value("LCOPR", "PROD", 5, "LGR1"), value("WOPR", "PROD", 0, ""), 1.0e-4);

    // Block vectors: an unrefined cell, and the LGR cell (2,2,1) of the producer.
    BOOST_CHECK_GT(value("BPR", none, 1, ""), 0.0);
    BOOST_CHECK_GT(value("LBPR", none, 5, "LGR1"), 0.0);

    std::filesystem::remove_all(outputDir);
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    // MPI setup.
    int argcDummy = 1;
    const char *tmp[] = {"test_lgr_summary", nullptr};
    char **argvDummy = const_cast<char**>(tmp);
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argcDummy, argvDummy);
#else
    Dune::MPIHelper::instance(argcDummy, argvDummy);
#endif

    Opm::FlowGenericVanguard::setCommunication(std::make_unique<Opm::Parallel::Communication>());

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
