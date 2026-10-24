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

#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/ESmry.hpp>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace
{

enum class SaturationCheck { WaterOnly, AllPhases, Dry };

struct ProbeCase
{
    std::vector<float> expectedSwat;
    SaturationCheck check = SaturationCheck::WaterOnly;
};

// RTEMP is 150 C in every cell. These values also apply when SWAT is
// inferred from SOIL and SGAS, or the supplied pair is rounded from float32.
const std::vector<float> rtempSwat {
    0.1626093f, 0.0392756f, 0.4371724f, 0.7565116f, 0.1626093f,
    0.0343808f, 0.0278085f, 0.2000000f, 0.2000000f};

// The EXPLICIT_SWAT probes use flash SWAT normalization. Their cells 6-8
// distinguish hydrocarbon input as gas or oil; SSHIFT changes molar volumes.
const std::map<std::string_view, ProbeCase> probeCases {
    {"EXPLICIT_SWAT_PROBES_SGAS",
     {{0.1626093f,
      0.0392756f,
      0.4371724f,
      0.7565116f,
      0.1111248f,
      0.2801100f,
      0.2697011f,
      0.8957410f,
      0.2000000f}}},
    {"EXPLICIT_SWAT_PROBES_SOIL",
     {{0.1626093f,
      0.0392756f,
      0.4371724f,
      0.7565116f,
      0.1111248f,
      0.0343808f,
      0.0278085f,
      0.2000000f,
      0.2000000f}}},
    {"EXPLICIT_SWAT_PROBES_SSHIFT",
     {{0.1614406f,
      0.0389520f,
      0.4350555f,
      0.7549224f,
      0.1074639f,
      0.2814685f,
      0.2704936f,
      0.9109980f,
      0.2000000f}}},
    {"EXPLICIT_SWAT_PROBES_XMF_YMF",
     {{0.1626093f,
      0.0392756f,
      0.4371724f,
      0.7565116f,
      0.1111248f,
      0.0612444f,
      0.0504184f,
      0.3269900f,
      0.2000000f}}},
    {"EXPLICIT_SWAT_NEAR_DEW_POINT", {{0.2000029f, 0.1626093f}}},
    {"EXPLICIT_SWAT_ZERO_COMPONENT", {{0.2000000f, 0.1830071f}}},
    {"COMP_EXPLICIT_SOIL_RTEMP", {rtempSwat, SaturationCheck::AllPhases}},
    {"COMP_EXPLICIT_SOIL_SGAS", {rtempSwat, SaturationCheck::AllPhases}},
    {"COMP_EXPLICIT_SATURATIONS_SCALED",
     {{0.2f, 0.05f, 0.5f, 0.8f, 0.2f, 0.2f, 0.2f, 0.2f, 0.2f}, SaturationCheck::AllPhases}},
    {"COMP_EXPLICIT_DRY_SOIL", {std::vector<float>(9), SaturationCheck::Dry}},
};

int
checkProbeCase(std::string_view caseName, const ProbeCase& probe)
{
    try {
        const auto& expected = probe.expectedSwat;
        Opm::EclIO::ERst restart(std::string {caseName} + ".UNRST");
        if (!restart.hasReportStepNumber(0)) {
            std::cerr << caseName << ": no initial restart step\n";
            return 1;
        }

        restart.loadReportStepNumber(0);
        std::vector<float> sgas;
        std::vector<float> soil;
        if (probe.check != SaturationCheck::WaterOnly) {
            sgas = restart.getRestartData<float>("SGAS", 0);
            soil = restart.getRestartData<float>("SOIL", 0);
            if (sgas.size() != expected.size() || soil.size() != expected.size()) {
                std::cerr << caseName << ": unexpected hydrocarbon saturation count\n";
                return 1;
            }
        }
        if (probe.check == SaturationCheck::Dry) {
            for (std::size_t cell = 0; cell < sgas.size(); ++cell) {
                if (!(sgas[cell] >= 0.0f && soil[cell] >= 0.0f
                      && std::abs(sgas[cell] + soil[cell] - 1.0f) <= 1.0e-5f)) {
                    std::cerr << caseName << " cell " << cell + 1
                              << ": dry phase saturations do not fill the pore volume\n";
                    return 1;
                }
            }
            std::cout << caseName << ": dry phase saturations fill the pore volume\n";
            return 0;
        }
        const auto& swat = restart.getRestartData<float>("SWAT", 0);
        if (swat.size() != expected.size()) {
            std::cerr << caseName << ": expected " << expected.size() << " SWAT values, got "
                      << swat.size() << '\n';
            return 1;
        }

        bool ok = true;
        for (std::size_t cell = 0; cell < swat.size(); ++cell) {
            if (!(std::abs(swat[cell] - expected[cell]) <= 1.0e-5f
                  && (probe.check == SaturationCheck::WaterOnly
                      || (sgas[cell] >= 0.0f && soil[cell] >= 0.0f
                          && std::abs(swat[cell] + sgas[cell] + soil[cell] - 1.0f) <= 1.0e-5f)))) {
                std::cerr << caseName << " cell " << cell + 1 << " initial SWAT: " << swat[cell]
                          << " (expected " << expected[cell] << ")\n";
                ok = false;
            }
        }

        if (ok) {
            std::cout << caseName << ": initial SWAT as expected in all cells\n";
        }
        return ok ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

} // namespace

int
main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: test_comp_explicit_swat_init <default|normalized|case>\n";
        return 1;
    }

    const std::string_view mode {argv[1]};
    if (const auto probe = probeCases.find(mode); probe != probeCases.end()) {
        return checkProbeCase(probe->first, probe->second);
    }

    if (mode != "default" && mode != "normalized") {
        std::cerr << "Unknown initialization mode: " << mode << '\n';
        return 1;
    }

    try {
        constexpr auto key = "BSWAT:15,1,1";
        Opm::EclIO::ESmry summary("1D_COMP_NO_WELLS_DUMMY_WATER.SMSPEC");
        if (!summary.hasKey(key)) {
            std::cerr << "Missing summary vector " << key << '\n';
            return 1;
        }

        summary.loadData({key});
        const auto& swat = summary.get(key);
        if (swat.size() != 100) {
            std::cerr << "Expected 100 report values for " << key << ", got " << swat.size()
                      << '\n';
            return 1;
        }

        // The middle cell barely changes over the first 0.01-day step. Its
        // first reported SWAT is about 0.163 with the scaling and 0.200
        // without it. The bounds leave room for small numerical changes but
        // cannot both pass if the option is ignored.
        const float lower = mode == "normalized" ? 0.15f : 0.19f;
        const float upper = mode == "normalized" ? 0.18f : 0.21f;
        const float first = swat.front();
        if (!std::isfinite(first) || first < lower || first > upper) {
            std::cerr << mode << " " << key << " first report: " << first << " (expected between "
                      << lower << " and " << upper << ")\n";
            return 1;
        }

        std::cout << mode << " " << key << " first report: " << first << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
