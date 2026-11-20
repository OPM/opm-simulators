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

#include <config.h>

#define BOOST_TEST_MODULE TestTomlConfigExamples

#ifndef HAVE_MPI
#define HAVE_MPI 0
#endif

#include <boost/test/unit_test.hpp>

#include <opm/common/utility/FileSystem.hpp>

#include <opm/simulators/flow/HybridNewtonConfig.hpp>
#include <opm/simulators/flow/TomlConfigExamples.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/TPSALinearSolverParameters.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/SetupPartitioningParams.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

class TempDir
{
public:
    TempDir()
        : dir_ { std::filesystem::temp_directory_path() /
                 Opm::unique_path("toml-examples-%%%%%%") }
    {
        std::filesystem::create_directories(dir_);
    }

    ~TempDir()
    {
        std::filesystem::remove_all(dir_);
    }

    const std::filesystem::path& path() const { return dir_; }

    std::string write(const std::string& name, const std::string& content) const
    {
        const auto file = dir_ / name;
        std::ofstream { file } << content;
        return file.generic_string();
    }

private:
    std::filesystem::path dir_;
};

bool parsesAsNumber(const std::string& s, double& value)
{
    if (s.empty()) {
        return false;
    }
    char* end = nullptr;
    value = std::strtod(s.c_str(), &end);
    return end == s.c_str() + s.size();
}

/// Check that two trees hold the same keys and the same values, comparing
/// numbers by value rather than by representation.
void checkSameTree(const Opm::PropertyTree& expected,
                   const Opm::PropertyTree& actual,
                   const std::string& path = "")
{
    // TOML writes sub tables after the plain values, so only compare the
    // sets of keys.
    auto keys = expected.get_child_keys();
    auto actualKeys = actual.get_child_keys();
    std::ranges::sort(keys);
    std::ranges::sort(actualKeys);
    BOOST_CHECK_MESSAGE(keys == actualKeys, "Different keys at '" << path << "'");

    if (keys.empty()) {
        const auto e = expected.get<std::string>("");
        const auto a = actual.get<std::string>("");
        double ev{}, av{};
        if (parsesAsNumber(e, ev) && parsesAsNumber(a, av)) {
            BOOST_CHECK_MESSAGE(std::abs(ev - av) <= 1.0e-14 * std::abs(ev),
                                "Value " << a << " differs from " << e
                                << " at '" << path << "'");
        }
        else {
            BOOST_CHECK_MESSAGE(e == a, "Value '" << a << "' differs from '" << e
                                << "' at '" << path << "'");
        }
        return;
    }

    for (const auto& key : keys) {
        if (actual.get_child_optional(key).has_value()) {
            checkSameTree(expected.get_child(key), actual.get_child(key),
                          path.empty() ? key : path + '.' + key);
        }
    }
}

/// Check that each key/value pair and each table header is documented by a
/// comment.
void checkAllDocumented(const std::string& toml)
{
    std::istringstream lines { toml };
    std::string previous;
    for (std::string line; std::getline(lines, line); previous = line) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        BOOST_CHECK_MESSAGE(!previous.empty() && previous.front() == '#',
                            "Undocumented entry: " << line);
    }
}

void checkLinearSolverExample(const Opm::PropertyTree& prm, const std::string& name)
{
    BOOST_TEST_MESSAGE("Linear solver preset " << name);

    const auto toml = Opm::linearSolverTomlExample(prm, "--linear-solver", name);
    BOOST_TEST_MESSAGE(toml);

    checkAllDocumented(toml);
    checkSameTree(prm, Opm::PropertyTree::fromTomlString(toml));
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(LinearSolverPresets)
{
    for (const auto* preset : { "cprw", "cpr", "cpr_trueimpes", "cpr_quasiimpes",
                                "cpr_trueimpesanalytic", "system_cpr", "system_cprw",
                                "amg", "ilu0", "dilu", "umfpack", "mixed-ilu0",
                                "mixed-dilu", "mixed-cprw" })
    {
        Opm::FlowLinearSolverParameters p;
        p.linsolver_ = preset;
        checkLinearSolverExample(Opm::setupPropertyTree(p, false, false), preset);
    }
}

BOOST_AUTO_TEST_CASE(TpsaLinearSolverPresets)
{
    for (const auto* preset : { "ilu0", "dilu", "amg" }) {
        Opm::TpsaLinearSolverParameters p;
        p.reset();
        p.linsolver_ = preset;
        checkLinearSolverExample(Opm::setupPropertyTree(p, false, false, true), preset);
    }
}

BOOST_AUTO_TEST_CASE(LinearSolverFromTomlFile)
{
    // A configuration file written from the example must give the same
    // configuration when used with --linear-solver.
    Opm::FlowLinearSolverParameters p;
    p.linsolver_ = "cprw";
    const auto prm = Opm::setupPropertyTree(p, false, false);

    const auto dir = TempDir{};
    const auto file = dir.write("cprw.toml",
                                Opm::linearSolverTomlExample(prm, "--linear-solver", "cprw"));

    p.linsolver_ = file;
    checkSameTree(prm, Opm::setupPropertyTree(p, false, false));
}

#if HAVE_MPI
BOOST_AUTO_TEST_CASE(ZoltanParams)
{
    for (const auto* conf : { "graph", "hypergraph", "scotch" }) {
        const auto params = Opm::setupZoltanParams(conf, 0.35);
        const auto toml = Opm::zoltanParamsTomlExample(params, conf, 1.1);
        BOOST_TEST_MESSAGE(toml);

        checkAllDocumented(toml);

        const auto dir = TempDir{};
        const auto file = dir.write("zoltan-params.toml", toml);
        BOOST_CHECK(Opm::setupZoltanParams(file) == params);
    }
}

BOOST_AUTO_TEST_CASE(MetisParams)
{
    const auto dir = TempDir{};
    {
        const auto params = Opm::setupMetisParams("default");
        const auto toml = Opm::metisParamsTomlExample(params, "default");
        BOOST_TEST_MESSAGE(toml);

        const auto file = dir.write("metis-default.toml", toml);
        BOOST_CHECK(Opm::setupMetisParams(file).empty());
    }

    {
        const auto params = std::map<std::string, std::string> {
            { "METIS_OPTION_PTYPE", "METIS_PTYPE_KWAY" },
            { "METIS_OPTION_NCUTS", "2" },
        };
        const auto toml = Opm::metisParamsTomlExample(params, "test");
        BOOST_TEST_MESSAGE(toml);

        checkAllDocumented(toml);
        const auto file = dir.write("metis.toml", toml);
        BOOST_CHECK(Opm::setupMetisParams(file) == params);
    }
}
#endif // HAVE_MPI

BOOST_AUTO_TEST_CASE(HybridNewton)
{
    const auto toml = Opm::hybridNewtonTomlExample();
    BOOST_TEST_MESSAGE(toml);

    auto pt = Opm::PropertyTree::fromTomlString(toml);
    BOOST_REQUIRE(pt.get_child_optional("model_1").has_value());

    const auto dir = TempDir{};
    pt.put("model_1.cell_indices_file", dir.write("cells.txt", "# cells\n0\n1\n\n2\n"));

    const auto config = Opm::HybridNewtonConfig { pt.get_child("model_1") };

    BOOST_CHECK_EQUAL(config.n_cells, std::size_t{3});
    BOOST_REQUIRE_EQUAL(config.apply_times.size(), std::size_t{2});
    BOOST_CHECK_EQUAL(config.apply_times[0], 0.0);
    BOOST_CHECK_EQUAL(config.apply_times[1], 8640000.0);

    BOOST_REQUIRE_EQUAL(config.input_features.size(), std::size_t{4});
    BOOST_CHECK_EQUAL(config.input_features[0].first, "PRESSURE");
    BOOST_CHECK(config.input_features[0].second.scaler.type == Opm::Scaler::Type::Standard);
    BOOST_CHECK_EQUAL(config.input_features[1].first, "SGAS");
    BOOST_CHECK(config.input_features[1].second.scaler.type == Opm::Scaler::Type::MinMax);
    BOOST_CHECK_EQUAL(config.input_features[2].first, "PERMX");
    BOOST_CHECK(config.input_features[2].second.transform.type == Opm::Transform::Type::Log10);
    BOOST_CHECK_EQUAL(config.input_features[3].first, "TIMESTEP");
    BOOST_CHECK(config.input_features[3].second.scaler.type == Opm::Scaler::Type::None);

    BOOST_REQUIRE_EQUAL(config.output_features.size(), std::size_t{2});
    BOOST_CHECK_EQUAL(config.output_features[0].first, "DELTA_PRESSURE");
    BOOST_CHECK(config.output_features[0].second.is_delta);
    BOOST_CHECK_EQUAL(config.output_features[0].second.actual_name, "PRESSURE");
}

BOOST_AUTO_TEST_CASE(HybridNewtonApplyTimes)
{
    // Single value, array and the object form used in JSON files.
    const auto dir = TempDir{};
    const auto cells = dir.write("cells.txt", "0\n");

    for (const auto& [applyTimes, expected] :
             std::vector<std::pair<std::string, std::vector<double>>> {
                 { "apply_times = 86400.0", { 86400.0 } },
                 { "apply_times = [0.0, 86400]", { 0.0, 86400.0 } },
                 { "apply_times = { \"0\" = 1.0, \"1\" = 2.0 }", { 1.0, 2.0 } },
             })
    {
        const auto pt = Opm::PropertyTree::fromTomlString(
            "model_path = \"m.model\"\n" + applyTimes + "\n"
            "cell_indices_file = \"" + cells + "\"\n");

        const auto config = Opm::HybridNewtonConfig { pt };
        BOOST_CHECK_EQUAL_COLLECTIONS(config.apply_times.begin(), config.apply_times.end(),
                                      expected.begin(), expected.end());
    }
}

BOOST_AUTO_TEST_CASE(WriteExamples)
{
    const auto dir = TempDir{};

    BOOST_CHECK_THROW(Opm::writeTomlExamples("no-such-example", dir.path()),
                      std::invalid_argument);
    BOOST_CHECK_THROW(Opm::writeTomlExamples(" , ", dir.path()),
                      std::invalid_argument);

    const auto written = Opm::writeTomlExamples("hybrid-newton", dir.path());
    BOOST_REQUIRE_EQUAL(written.size(), std::size_t{1});
    BOOST_CHECK(written[0] == dir.path() / "hybrid-newton.toml");
    BOOST_CHECK(std::filesystem::exists(written[0]));

    // Existing files must not be overwritten.
    BOOST_CHECK_THROW(Opm::writeTomlExamples("hybrid-newton", dir.path()),
                      std::invalid_argument);
}
