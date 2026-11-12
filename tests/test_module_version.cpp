#include <opm/simulators/utils/moduleVersion.hpp>

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Expected the build configuration\n";
        return 1;
    }

    const auto name = Opm::moduleVersionName();
    const auto hash = Opm::moduleVersionHash();
    const auto version = Opm::moduleVersion();
    const auto timestamp = Opm::compileTimestamp();
    const auto buildType = Opm::moduleBuildType();

    if (name.empty() || version != name + " (" + hash + ")") {
        std::cerr << "Inconsistent module version information\n";
        return 1;
    }
    if ("config=" + buildType != argv[1]) {
        std::cerr << "Unexpected build configuration: " << buildType << '\n';
        return 1;
    }

    std::cout << version << '\n' << timestamp << '\n';
    return 0;
}
