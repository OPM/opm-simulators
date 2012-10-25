/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/pressure/FlowBCManager.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/newwells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/initState.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/transport/reorder/TransportModelTracerTof.hpp>
#include <opm/core/transport/reorder/TransportModelTracerTofDiscGal.hpp>

#include <boost/scoped_ptr.hpp>
#include <boost/filesystem.hpp>

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
{
    using namespace Opm;

    parameter::ParameterGroup param(argc, argv, false);

    // Read grid.
    GridManager grid_manager(param.get<std::string>("grid_filename"));
    const UnstructuredGrid& grid = *grid_manager.c_grid();

    // Read porosity, compute pore volume.
    std::vector<double> porevol;
    {
        std::ifstream poro_stream(param.get<std::string>("poro_filename").c_str());
        std::istream_iterator<double> beg(poro_stream);
        std::istream_iterator<double> end;
        porevol.assign(beg, end); // Now contains poro.
        if (int(porevol.size()) != grid.number_of_cells) {
            THROW("Size of porosity field differs from number of cells.");
        }
        for (int i = 0; i < grid.number_of_cells; ++i) {
            porevol[i] *= grid.cell_volumes[i];
        }
    }

    // Read flux.
    std::vector<double> flux;
    {
        std::ifstream flux_stream(param.get<std::string>("flux_filename").c_str());
        std::istream_iterator<double> beg(flux_stream);
        std::istream_iterator<double> end;
        flux.assign(beg, end);
        if (int(flux.size()) != grid.number_of_faces) {
            THROW("Size of flux field differs from number of faces.");
        }
    }

    // Read source terms.
    std::vector<double> src;
    {
        std::ifstream src_stream(param.get<std::string>("src_filename").c_str());
        std::istream_iterator<double> beg(src_stream);
        std::istream_iterator<double> end;
        src.assign(beg, end);
        if (int(src.size()) != grid.number_of_cells) {
            THROW("Size of source term field differs from number of cells.");
        }
    }

    // Choice of tof solver.
    bool use_dg = param.getDefault("use_dg", false);
    int dg_degree = -1;
    bool use_cvi = false;
    if (use_dg) {
        dg_degree = param.getDefault("dg_degree", 0);
        use_cvi = param.getDefault("use_cvi", false);
    }

    // Write parameters used for later reference.
    bool output = param.getDefault("output", true);
    std::ofstream epoch_os;
    std::string output_dir;
    if (output) {
        output_dir =
            param.getDefault("output_dir", std::string("output"));
        boost::filesystem::path fpath(output_dir);
        try {
            create_directories(fpath);
        }
        catch (...) {
            THROW("Creating directories failed: " << fpath);
        }
        param.writeParam(output_dir + "/simulation.param");
    }

    // Issue a warning if any parameters were unused.
    warnIfUnusedParams(param);

    // Solve time-of-flight.
    Opm::time::StopWatch transport_timer;
    transport_timer.start();
    std::vector<double> tof;
    if (use_dg) {
        Opm::TransportModelTracerTofDiscGal tofsolver(grid, use_cvi);
        tofsolver.solveTof(&flux[0], &porevol[0], &src[0], dg_degree, tof);
    } else {
        Opm::TransportModelTracerTof tofsolver(grid);
        tofsolver.solveTof(&flux[0], &porevol[0], &src[0], tof);
    }
    transport_timer.stop();
    double tt = transport_timer.secsSinceStart();
    std::cout << "Transport solver took: " << tt << " seconds." << std::endl;

    // Output.
    if (output) {
        std::string tof_filename = output_dir + "/tof.txt";
        std::ofstream tof_stream(tof_filename.c_str());
        std::copy(tof.begin(), tof.end(), std::ostream_iterator<double>(tof_stream, "\n"));
    }
}
