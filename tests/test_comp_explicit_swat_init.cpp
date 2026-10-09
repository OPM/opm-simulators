#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/ESmry.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <map>
#include <string>
#include <string_view>

namespace
{

// Initial SWAT per cell with --normalize-explicit-swat=true. Cells 6-8 depend
// on whether SGAS or SOIL is given, and SSHIFT moves every two-phase cell.
const std::map<std::string_view, std::array<float, 9>> probeCases {
    {"EXPLICIT_SWAT_PROBES_SGAS",
     {0.1626093f,
      0.0392756f,
      0.4371724f,
      0.7565116f,
      0.1111248f,
      0.2801100f,
      0.2697011f,
      0.8957410f,
      0.2000000f}},
    {"EXPLICIT_SWAT_PROBES_SOIL",
     {0.1626093f,
      0.0392756f,
      0.4371724f,
      0.7565116f,
      0.1111248f,
      0.0343808f,
      0.0278085f,
      0.2000000f,
      0.2000000f}},
    {"EXPLICIT_SWAT_PROBES_SSHIFT",
     {0.1614406f,
      0.0389520f,
      0.4350555f,
      0.7549224f,
      0.1074639f,
      0.2814685f,
      0.2704936f,
      0.9109980f,
      0.2000000f}},
    {"EXPLICIT_SWAT_PROBES_XMF_YMF",
     {0.1626093f,
      0.0392756f,
      0.4371724f,
      0.7565116f,
      0.1111248f,
      0.0612444f,
      0.0504184f,
      0.3269900f,
      0.2000000f}},
};

int
checkProbeCase(std::string_view caseName, const std::array<float, 9>& expected)
{
    try {
        Opm::EclIO::ERst restart(std::string {caseName} + ".UNRST");
        if (!restart.hasReportStepNumber(0)) {
            std::cerr << caseName << ": no initial restart step\n";
            return 1;
        }

        restart.loadReportStepNumber(0);
        const auto& swat = restart.getRestartData<float>("SWAT", 0);
        if (swat.size() != expected.size()) {
            std::cerr << caseName << ": expected " << expected.size() << " SWAT values, got "
                      << swat.size() << '\n';
            return 1;
        }

        bool ok = true;
        for (std::size_t cell = 0; cell < swat.size(); ++cell) {
            if (!(std::abs(swat[cell] - expected[cell]) <= 1.0e-5f)) {
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
        std::cerr << "Usage: test_comp_explicit_swat_init <default|normalized|probe case>\n";
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
