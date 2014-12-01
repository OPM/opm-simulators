/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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



#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/core/tof/AnisotropicEikonal.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <boost/filesystem.hpp>
#include <memory>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <iterator>


namespace
{
    void warnIfUnusedParams(const Opm::parameter::ParameterGroup& param)
    {
        if (param.anyUnused()) {
            std::cout << "--------------------   Warning: unused parameters:   --------------------\n";
            param.displayUsage();
            std::cout << "-------------------------------------------------------------------------" << std::endl;
        }
    }
} // anon namespace



// ----------------- Main program -----------------
int
main(int argc, char** argv)
try
{
    using namespace Opm;

    parameter::ParameterGroup param(argc, argv, false);

    // Read grid.
    GridManager grid_manager(param.get<std::string>("grid_filename"));
    const UnstructuredGrid& grid = *grid_manager.c_grid();

    // Read metric tensor.
    std::vector<double> metric;
    {
        std::ifstream metric_stream(param.get<std::string>("metric_filename").c_str());
        std::istream_iterator<double> beg(metric_stream);
        std::istream_iterator<double> end;
        metric.assign(beg, end);
        if (int(metric.size()) != grid.number_of_cells*grid.dimensions*grid.dimensions) {
            OPM_THROW(std::runtime_error, "Size of metric field differs from (dim^2 * number of cells).");
        }
    }

    // Read starting cells.
    std::vector<int> startcells;
    {
        std::ifstream start_stream(param.get<std::string>("startcells_filename").c_str());
        std::istream_iterator<int> beg(start_stream);
        std::istream_iterator<int> end;
        startcells.assign(beg, end);
    }

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    std::string output_dir;
    if (output) {
        output_dir =
            param.getDefault("output_dir", std::string("output"));
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        param.writeParam(output_dir + "/eikonal.param");
    }

    // Issue a warning if any parameters were unused.
    warnIfUnusedParams(param);

    // Solve eikonal equation.
    Opm::time::StopWatch timer;
    timer.start();
    std::vector<double> solution;
    AnisotropicEikonal2d ae(grid);
    ae.solve(metric.data(), startcells, solution);
    timer.stop();
    double tt = timer.secsSinceStart();
    std::cout << "Eikonal solver took: " << tt << " seconds." << std::endl;

    // Output.
    if (output) {
        std::string filename = output_dir + "/solution.txt";
        std::ofstream stream(filename.c_str());
        stream.precision(16);
        std::copy(solution.begin(), solution.end(), std::ostream_iterator<double>(stream, "\n"));
    }
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
