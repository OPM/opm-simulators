#include <opm/io/eclipse/ESmry.hpp>

#include <cmath>
#include <exception>
#include <iostream>
#include <string_view>

int
main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: test_comp_explicit_swat_init <default|normalized>\n";
        return 1;
    }

    const std::string_view mode {argv[1]};
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
