/*
  Copyright 2026 SINTEF Digital

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

// Checks the results of COMP_NUMERICAL_AQUIFER.DATA, run in the current directory.

#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/ESmry.hpp>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace
{

const std::string caseName {"COMP_NUMERICAL_AQUIFER"};

// Cells 5 and 6 of the row are the aquifer, given 230 bar in AQUNUM.
constexpr std::size_t firstAquiferCell = 4;
constexpr std::size_t numCells = 6;
constexpr float aquiferPressure = 230.0f;

bool
checkInitialState()
{
    Opm::EclIO::ERst restart(caseName + ".UNRST");
    if (!restart.hasReportStepNumber(0)) {
        std::cerr << "No initial restart step\n";
        return false;
    }

    restart.loadReportStepNumber(0);
    const auto& swat = restart.getRestartData<float>("SWAT", 0);
    const auto& pressure = restart.getRestartData<float>("PRESSURE", 0);
    if (swat.size() != numCells || pressure.size() != numCells) {
        std::cerr << "Expected " << numCells << " cells in the initial restart step\n";
        return false;
    }

    bool ok = true;
    for (std::size_t cell = 0; cell < numCells; ++cell) {
        const bool aquifer = cell >= firstAquiferCell;
        // Above the water-oil contact the reservoir holds connate water only.
        const float expectedSwat = aquifer ? 1.0f : 0.2f;
        if (!(std::abs(swat[cell] - expectedSwat) <= 1.0e-5f)) {
            std::cerr << "Cell " << cell + 1 << " initial SWAT: " << swat[cell] << " (expected "
                      << expectedSwat << ")\n";
            ok = false;
        }
        if (aquifer && !(std::abs(pressure[cell] - aquiferPressure) <= 1.0e-3f)) {
            std::cerr << "Cell " << cell + 1 << " initial pressure: " << pressure[cell]
                      << " (expected " << aquiferPressure << ")\n";
            ok = false;
        }
        if (!aquifer && !(pressure[cell] < aquiferPressure - 20.0f)) {
            std::cerr << "Cell " << cell + 1 << " initial pressure: " << pressure[cell]
                      << " (expected the equilibrium pressure of about 200)\n";
            ok = false;
        }
    }
    return ok;
}

bool
checkAquiferVectors()
{
    Opm::EclIO::ESmry summary(caseName + ".SMSPEC");
    const std::vector<std::string> keys {
        "TIME", "ANQR:1", "ANQT:1", "ANQP:1", "BPR:1,1,1", "BSWAT:1,1,1"};
    for (const auto& key : keys) {
        if (!summary.hasKey(key)) {
            std::cerr << "Missing summary vector " << key << '\n';
            return false;
        }
    }
    summary.loadData(keys);

    const auto& time = summary.get("TIME");
    const auto& rate = summary.get("ANQR:1");
    const auto& total = summary.get("ANQT:1");
    const auto& aquPressure = summary.get("ANQP:1");
    const auto& resPressure = summary.get("BPR:1,1,1");
    const auto& swat = summary.get("BSWAT:1,1,1");
    // The values at the start show no flow yet.
    std::size_t first = 0;
    while (first < time.size() && time[first] == 0.0f) {
        ++first;
    }
    if (rate.size() < first + 2) {
        std::cerr << "Expected several report values, got " << rate.size() - first << '\n';
        return false;
    }

    bool ok = true;
    for (std::size_t step = first; step < rate.size(); ++step) {
        // Water flows from the aquifer into the reservoir as long as the
        // aquifer pressure is the higher one.
        if (!(rate[step] > 0.0f) || !(aquPressure[step] > resPressure[step])
            || !(aquPressure[step] < aquiferPressure)) {
            std::cerr << "Report " << step << ": ANQR " << rate[step] << ", ANQP "
                      << aquPressure[step] << ", BPR:1,1,1 " << resPressure[step] << '\n';
            ok = false;
        }
        if (step > first && !(rate[step] < rate[step - 1] && total[step] > total[step - 1])) {
            std::cerr << "Report " << step << ": ANQR " << rate[step] << " after " << rate[step - 1]
                      << ", ANQT " << total[step] << " after " << total[step - 1] << '\n';
            ok = false;
        }
    }

    // The influx reaches the connected cell.
    if (!(swat.back() > 0.2f)) {
        std::cerr << "BSWAT:1,1,1 at the end: " << swat.back() << '\n';
        ok = false;
    }

    if (ok) {
        std::cout << "ANQR:1 from " << rate[first] << " to " << rate.back() << ", ANQT:1 "
                  << total.back() << ", ANQP:1 " << aquPressure.back() << '\n';
    }
    return ok;
}

} // namespace

int
main()
{
    try {
        const bool initial = checkInitialState();
        const bool vectors = checkAquiferVectors();
        return (initial && vectors) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
