// SPDX-License-Identifier: GPL-3.0-or-later
#include <opm/simulators/utils/moduleVersion.hpp>
namespace Opm
{
std::string
moduleVersionName()
{
#ifdef RESLAB_CPU_REFERENCE
    return "2026.04-cpuref";
#else
    return "2026.04-vulkan";
#endif
}
std::string
moduleVersionHash()
{
    return "opm-2026.04-reslab-vulkan";
}
std::string
moduleVersion()
{
    return "2026.04";
}
std::string
compileTimestamp()
{
    return __DATE__ " " __TIME__;
}
} // namespace Opm
