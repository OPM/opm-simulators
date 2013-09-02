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
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/wells.h>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/SparseTable.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/props/IncompPropertiesFromDeck.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/initState.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/tof/TofReorder.hpp>
#include <opm/core/tof/TofDiscGalReorder.hpp>

#include <memory>
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
try
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
            OPM_THROW(std::runtime_error, "Size of porosity field differs from number of cells.");
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
            OPM_THROW(std::runtime_error, "Size of flux field differs from number of faces.");
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
            OPM_THROW(std::runtime_error, "Size of source term field differs from number of cells.");
        }
    }

    const bool compute_tracer = param.getDefault("compute_tracer", false);
    Opm::SparseTable<int> tracerheads;
    if (compute_tracer) {
        std::ifstream tr_stream(param.get<std::string>("tracerheads_filename").c_str());
        int num_rows;
        tr_stream >> num_rows;
        for (int row = 0; row < num_rows; ++row) {
            int row_size;
            tr_stream >> row_size;
            std::vector<int> rowdata(row_size);
            for (int elem = 0; elem < row_size; ++elem) {
                tr_stream >> rowdata[elem];
            }
            tracerheads.appendRow(rowdata.begin(), rowdata.end());
        }
    }

    // Choice of tof solver.
    bool use_dg = param.getDefault("use_dg", false);
    bool use_multidim_upwind = false;
    // Need to initialize dg solver here, since it uses parameters now.
    std::unique_ptr<Opm::TofDiscGalReorder> dg_solver;
    if (use_dg) {
        dg_solver.reset(new Opm::TofDiscGalReorder(grid, param));
    } else {
        use_multidim_upwind = param.getDefault("use_multidim_upwind", false);
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
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        param.writeParam(output_dir + "/simulation.param");
    }

    // Issue a warning if any parameters were unused.
    warnIfUnusedParams(param);

    // Solve time-of-flight.
    Opm::time::StopWatch transport_timer;
    transport_timer.start();
    std::vector<double> tof;
    std::vector<double> tracer;
    if (use_dg) {
        if (compute_tracer) {
            dg_solver->solveTofTracer(&flux[0], &porevol[0], &src[0], tracerheads, tof, tracer);
        } else {
            dg_solver->solveTof(&flux[0], &porevol[0], &src[0], tof);
        }
    } else {
        Opm::TofReorder tofsolver(grid, use_multidim_upwind);
        if (compute_tracer) {
            tofsolver.solveTofTracer(&flux[0], &porevol[0], &src[0], tracerheads, tof, tracer);
        } else {
            tofsolver.solveTof(&flux[0], &porevol[0], &src[0], tof);
        }
    }
    transport_timer.stop();
    double tt = transport_timer.secsSinceStart();
    std::cout << "Transport solver took: " << tt << " seconds." << std::endl;

    // Output.
    if (output) {
        std::string tof_filename = output_dir + "/tof.txt";
        std::ofstream tof_stream(tof_filename.c_str());
        tof_stream.precision(16);
        std::copy(tof.begin(), tof.end(), std::ostream_iterator<double>(tof_stream, "\n"));
        if (compute_tracer) {
            std::string tracer_filename = output_dir + "/tracer.txt";
            std::ofstream tracer_stream(tracer_filename.c_str());
            tracer_stream.precision(16);
            const int nt = tracer.size()/grid.number_of_cells;
            for (int i = 0; i < nt*grid.number_of_cells; ++i) {
                tracer_stream << tracer[i] << (((i + 1) % nt == 0) ? '\n' : ' ');
            }
        }
    }
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
