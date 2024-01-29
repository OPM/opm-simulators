/*
  Copyright 2020 Equinor.

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

#include <opm/simulators/utils/HDF5Serializer.hpp>

#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <boost/date_time.hpp>

#include <fmt/format.h>

#include <array>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Need one parameter, the .OPMRST file to inspect\n";
        return 1;
    }

    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    Opm::Parallel::Communication comm{MPI_COMM_SELF};
#else
    Opm::Parallel::Communication comm{};
#endif

    Opm::HDF5Serializer ser(argv[1], Opm::HDF5File::OpenMode::READ, comm);

    std::tuple<std::array<std::string,5>,int> header;
    try {
        ser.read(header, "/", "simulator_info", Opm::HDF5File::DataSetMode::ROOT_ONLY);
    } catch(...) {
        std::cerr << "Error reading data from file, is it really a .OPMRST file?\n";
        return 2;
    }

    const auto& [strings, procs] = header;
    std::cout << "Info for " << argv[1] <<":\n";
    std::cout << fmt::format("\tSimulator name: {}\n"
                             "\tSimulator version: {}\n"
                             "\tCompile time stamp: {}\n"
                             "\tCase name: {}\n"
                             "\tNumber of processes: {}\n",
                             strings[0], strings[1], strings[2], strings[3], procs);

    std::cout << fmt::format("\tLast report step: {}\n", ser.lastReportStep());
    const std::vector<int> reportSteps = ser.reportSteps();
    for (int step : reportSteps) {
        Opm::SimulatorTimer timer;
        try {
            ser.read(timer, fmt::format("/report_step/{}", step),
                     "simulator_timer", Opm::HDF5File::DataSetMode::ROOT_ONLY);
        } catch (...) {
            std::cerr << "*** Failed to read timer info for level " << step << std::endl;
        }
        std::cout << "\t\tReport step id " << step << ": Time "
                  << timer.currentDateTime() << std::endl;
    }

    std::cout << "====== Parameter values ====\n"
              << strings[4];

    return 0;
}
