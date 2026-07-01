/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
  Copyright 2014 STATOIL ASA.

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
#include <opm/simulators/flow/Banners.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <fmt/format.h>

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <cstdlib> // getenv
#else
#include <sys/utsname.h>
#include <unistd.h>
#endif

namespace {

unsigned long long getTotalSystemMemory()
{
#if defined(_WIN32)
    // Windows has no sysconf(_SC_PHYS_PAGES); query the total physical memory.
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status)) {
        return static_cast<unsigned long long>(status.ullTotalPhys);
    }
    return 0;
#else
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#endif
}

}

namespace Opm {

void printPRTHeader(const int nprocs, const int nthreads,
                    const std::string& parameters,
                    std::string_view moduleVersion,
                    std::string_view compileTimestamp)
{
    const double megabyte = 1024 * 1024;
    unsigned num_cpu = std::thread::hardware_concurrency();
#if defined(_WIN32)
    const char* user = getenv("USERNAME");  // no getlogin() on Windows
#else
    const char* user = getlogin();
#endif
    std::time_t now = std::time(0);
    struct std::tm  tstruct;
    char      tmstr[80];
    tstruct = *std::localtime(&now);
    std::strftime(tmstr, sizeof(tmstr), "%d-%m-%Y at %X", &tstruct);
    const double mem_size = getTotalSystemMemory() / megabyte;
    std::ostringstream ss;
    ss << "\n\n\n";
    ss << " ########  #          ######   #           #\n";
    ss << " #         #         #      #   #         # \n";
    ss << " #####     #         #      #    #   #   #  \n";
    ss << " #         #         #      #     # # # #   \n";
    ss << " #         #######    ######       #   #    \n\n";
    ss << "Flow is a simulator for fully implicit three-phase black-oil flow,";
    ss << " and is part of OPM.\nFor more information visit: https://opm-project.org \n\n";
    ss << "Flow Version     =  " << moduleVersion << "\n";
#if defined(_WIN32)
    // Windows has no uname()/utsname; query the machine name via the Win32 API.
    char computer_name[MAX_COMPUTERNAME_LENGTH + 1] = {};
    DWORD computer_name_len = sizeof(computer_name);
    const char* nodename =
        GetComputerNameA(computer_name, &computer_name_len) ? computer_name : "unknown";
    ss << "Machine name     =  " << nodename << " (Number of logical cores: " << num_cpu;
    ss << ", Memory size: " << std::fixed << std::setprecision(2) << mem_size << " MB) \n";
    ss << "Build time       =  " << compileTimestamp << "\n";
#else
    struct utsname arch;
    if (uname(&arch) == 0) {
       ss << "Machine name     =  " << arch.nodename << " (Number of logical cores: " << num_cpu;
       ss << ", Memory size: " << std::fixed << std::setprecision (2) << mem_size << " MB) \n";
       ss << "Operating system =  " << arch.sysname << " " << arch.machine << " (Kernel: " << arch.release;
       ss << ", " << arch.version << " )\n";
       ss << "Build time       =  " << compileTimestamp << "\n";
    }
#endif
    if (user) {
       ss << "User             =  " << user << std::endl;
    }
    ss << "Simulation started on " << tmstr << " hrs\n";
    ss << "Using "<< nprocs << " MPI processes with "<< nthreads <<" OMP threads on each \n";
    ss << "Parameters used by Flow:\n" << parameters;

    OpmLog::note(ss.str());
}

void printFlowBanner(int nprocs, int nthreads, std::string_view moduleVersionName)
{
    const int lineLen = 70;
    std::string banner = "This is flow ";
    banner += moduleVersionName;
    const int bannerPreLen = (lineLen - 2 - banner.size())/2;
    const int bannerPostLen = bannerPreLen + (lineLen - 2 - banner.size())%2;
    std::cout << "**********************************************************************\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*" << std::string(bannerPreLen, ' ') << banner << std::string(bannerPostLen, ' ') << "*\n";
    std::cout << "*                                                                    *\n";
    std::cout << "* Flow is a simulator for fully implicit three-phase black-oil flow, *\n";
    std::cout << "*             including solvent and polymer capabilities.            *\n";
    std::cout << "*          For more information, see https://opm-project.org         *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "**********************************************************************\n\n";

    std::cout << "Using "<< nprocs << " MPI processes with "<< nthreads <<" OMP threads on each \n\n";
}

void printFlowTrailer(int nprocs,
                      int nthreads,
                      const double total_setup_time,
                      const double deck_read_time,
                      const SimulatorReport& report)
{
    std::ostringstream ss;
    ss << "\n\n================    End of simulation     ===============\n\n";
    ss << fmt::format("Number of MPI processes: {:9}\n", nprocs);
    ss << fmt::format("Threads per MPI process: {:9}\n", nthreads);
    ss << fmt::format("Setup time:                 {:9.2f} s\n", total_setup_time);
    ss << fmt::format("  Deck input:               {:9.2f} s\n", deck_read_time);
    report.reportFullyImplicit(ss);
    OpmLog::info(ss.str());
}

} // namespace Opm
