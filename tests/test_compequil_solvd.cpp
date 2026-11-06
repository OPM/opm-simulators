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

#define BOOST_TEST_MODULE CompEquilSolvd
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include <flowexperimental/comp/flow_comp.hpp>

#include <opm/io/eclipse/ERst.hpp>

#include "RunDeckHelper.hpp"

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

// Two columns of ten cells, each its own equilibration region with water
// under a liquid column. The layers thicken with depth, and column 2 lies
// 10 m below column 1. The bottom cell of column 1 uses saturation table 2,
// whose connate water is 0.15 against 0.01 for table 1.
std::string deck(const std::string& rptsol)
{
    return R"(RUNSPEC
METRIC
TABDIMS
2 1 /
EQLDIMS
2 /
OIL
GAS
WATER
DIMENS
2 1 10 /
COMPS
3 /
FULLIMP
START
1 'JAN' 2016 /
UNIFOUT

GRID
DX
20*10 /
DY
20*10 /
DZ
2*2 2*2 2*4 2*4 2*6 2*6 2*8 2*8 2*10 2*10 /
BOX
1 2 1 1 1 1 /
TOPS
2000 2010 /
ENDBOX
PERMX
20*100 /
PERMY
20*100 /
PERMZ
20*100 /
PORO
20*0.2 /

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
0.5  0.30 0.0
0.99 0.70 0.0
/
0.0  0.00 0.0
0.1  0.00 0.0
0.5  0.30 0.0
0.85 0.70 0.0
/
SWFN
0.01 0.00 0.0
0.40 0.05 0.0
0.70 0.40 0.0
1.00 1.00 0.0
/
0.15 0.00 0.0
0.40 0.05 0.0
0.70 0.40 0.0
1.00 1.00 0.0
/
SOF3
0.0  0.00 0.0
0.2  0.00 0.0
0.5  0.20 0.2
0.99 0.70 0.7
/
0.0  0.00 0.0
0.2  0.00 0.0
0.5  0.20 0.2
0.85 0.70 0.7
/
STCOND
15.0 1.0 /
ROCK
68.9476 0 /
PVTW
200.0 1.01 4.0E-5 0.3 0.0 /
DENSITY
800.0 1000.0 1.0 /
ZMFVD
2000 0.1 0.3 0.6
2100 0.1 0.3 0.6 /
2000 0.1 0.3 0.6
2100 0.1 0.3 0.6 /
RTEMP
100 /

REGIONS
SATNUM
18*1 2 1 /
EQLNUM
1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 /

SOLUTION
EQUIL
2020 200 2041.5 0 2000 0 /
2030 205 2050.0 0 2010 0 /
RPTSOL
)" + rptsol + R"( /

SUMMARY
FPR

SCHEDULE
TSTEP
1 /
END
)";
}

constexpr std::array<double, 10> thickness{2, 2, 4, 4, 6, 6, 8, 8, 10, 10};
constexpr std::array<double, 2> columnTop{2000.0, 2010.0};
constexpr std::array<double, 2> waterOilContact{2041.5, 2050.0};
constexpr std::array<double, 2> connateWater{0.15, 0.01};

using TypeTag = Opm::Properties::TTag::FlowCompProblem<3, true>;

// The PRT file and the initial cell pressures in bar of a run, read on rank 0.
struct RunOutput
{
    std::string prt;
    std::vector<float> pressure;
};

// Runs the deck through the simulator's main entry point on every rank. Rank 0
// writes the deck to a directory of its own and shares its path, so that the
// ranks write to one output directory, and reads back what the run wrote.
RunOutput run(const std::string& name, const std::string& deck_text)
{
    const auto& comm = Opm::FlowGenericVanguard::comm();
    std::optional<Opm::TemporaryDirectory> dir;
    std::string path;
    if (comm.rank() == 0) {
        dir.emplace("compequil_solvd");
        path = dir->path().string();
        std::ofstream(dir->path() / (name + ".DATA")) << deck_text;
    }
    int size = static_cast<int>(path.size());
    comm.broadcast(&size, 1, 0);
    path.resize(size);
    comm.broadcast(path.data(), size, 0);

    const std::filesystem::path out(path);
    std::vector<std::string> args = {boost::unit_test::framework::master_test_suite().argv[0],
                                     (out / (name + ".DATA")).string(),
                                     "--output-dir=" + path};
    std::vector<char*> argv;
    for (auto& arg : args) {
        argv.push_back(arg.data());
    }
    argv.push_back(nullptr);
    Opm::Parameters::reset();
    const int status = Opm::Main(static_cast<int>(args.size()), argv.data(), /*ownMPI=*/false)
                           .runStatic<TypeTag>();
    // Every rank is done with the directory, which goes with rank 0's copy,
    // and none waits for another once a check below fails.
    comm.barrier();
    BOOST_REQUIRE_EQUAL(status, EXIT_SUCCESS);

    RunOutput output;
    if (comm.rank() == 0) {
        std::ifstream prt(out / (name + ".PRT"));
        output.prt.assign(std::istreambuf_iterator<char>(prt), std::istreambuf_iterator<char>());
        Opm::EclIO::ERst rst(out / (name + ".UNRST"));
        BOOST_REQUIRE(rst.hasReportStepNumber(0));
        rst.loadReportStepNumber(0);
        output.pressure = rst.getRestartData<float>("PRESSURE", 0);
    }
    return output;
}

// Columns of a report's main table.
enum Column : std::size_t { Depth, Poil, Pwat, Pgas, Deno, Denw, Deng,
                            Soil, Swat, Sgas, Psat, StatedPsat };

// One region's report as printed: the numbers of its heading and its main table.
struct PrintedReport
{
    int region{};
    int satTable{};
    std::vector<std::vector<double>> rows;
};

std::vector<PrintedReport> printedReports(const std::string& prt)
{
    constexpr std::string_view heading = "Composition versus depth in equilibration region";
    std::vector<PrintedReport> reports;
    std::istringstream lines(prt);
    std::string line;
    while (std::getline(lines, line)) {
        const auto pos = line.find(heading);
        if (pos != std::string::npos) {
            auto& report = reports.emplace_back();
            BOOST_REQUIRE_EQUAL(std::sscanf(line.c_str() + pos + heading.size(),
                                            " %d (saturation table %d)",
                                            &report.region, &report.satTable), 2);
            continue;
        }
        if (reports.empty() || reports.back().rows.size() == 50) {
            continue;
        }
        std::istringstream fields(line);
        std::vector<double> row;
        double value{};
        while (fields >> value) {
            row.push_back(value);
        }
        // The heading's units and rules hold no numbers.
        if (row.size() == 12 && fields.eof()) {
            reports.back().rows.push_back(std::move(row));
        }
    }
    return reports;
}

// A column of the main table interpolated linearly to a depth between two rows.
double atDepth(const std::vector<std::vector<double>>& rows,
               const Column column,
               const double depth)
{
    std::size_t i = 1;
    while (i + 1 < rows.size() && rows[i][Depth] < depth) {
        ++i;
    }
    const auto& upper = rows[i - 1];
    const auto& lower = rows[i];
    const double t = (depth - upper[Depth]) / (lower[Depth] - upper[Depth]);
    return (1.0 - t) * upper[column] + t * lower[column];
}

} // namespace

using RunDeckFixture = Opm::RunDeckFixture;
BOOST_GLOBAL_FIXTURE(RunDeckFixture);

// With RPTSOL SOLVD each region is reported from the top face of its cells to
// the bottom face, under the saturation table of its last cell, whose connate
// water fills the rows above the water-oil contact. The report's pressures
// pass through those the cells start from: the water column's below the
// contact and the liquid column's above it.
BOOST_AUTO_TEST_CASE(ReportsEachRegionFromItsCells)
{
    const auto output = run("SOLVD", deck("RESTART=2 SOLVD"));
    if (Opm::FlowGenericVanguard::comm().rank() != 0) {
        return;
    }

    const auto reports = printedReports(output.prt);
    BOOST_REQUIRE_EQUAL(reports.size(), std::size_t{2});
    BOOST_REQUIRE_EQUAL(output.pressure.size(), std::size_t{20});
    for (std::size_t r = 0; r < reports.size(); ++r) {
        BOOST_TEST_CONTEXT("Region " << r + 1) {
            const auto& report = reports[r];
            BOOST_CHECK_EQUAL(report.region, static_cast<int>(r) + 1);
            BOOST_CHECK_EQUAL(report.satTable, (r == 0) ? 2 : 1);
            BOOST_REQUIRE_EQUAL(report.rows.size(), std::size_t{50});
            BOOST_CHECK_EQUAL(report.rows.front()[Depth], columnTop[r]);
            BOOST_CHECK_EQUAL(report.rows.back()[Depth], columnTop[r] + 60.0);

            for (const auto& row : report.rows) {
                const bool water = row[Depth] > waterOilContact[r];
                BOOST_CHECK_EQUAL(row[Swat], water ? 1.0 : connateWater[r]);
            }

            double top = columnTop[r];
            for (std::size_t k = 0; k < thickness.size(); ++k) {
                const double depth = top + thickness[k] / 2;
                top += thickness[k];
                const auto column = (depth > waterOilContact[r]) ? Pwat : Poil;
                BOOST_TEST_CONTEXT("Cell at " << depth << " m") {
                    // The report prints pressures to 0.001 bar.
                    BOOST_CHECK_SMALL(atDepth(report.rows, column, depth)
                                      - output.pressure[r + 2 * k], 0.002);
                }
            }
        }
    }
}

// Without the SOLVD mnemonic there is no report.
BOOST_AUTO_TEST_CASE(NoReportUnlessRequested)
{
    const auto output = run("NO_SOLVD", deck("RESTART=2 PRESSURE"));
    if (Opm::FlowGenericVanguard::comm().rank() != 0) {
        return;
    }
    BOOST_REQUIRE(!output.prt.empty());
    BOOST_CHECK(output.prt.find("Composition versus depth") == std::string::npos);
}
