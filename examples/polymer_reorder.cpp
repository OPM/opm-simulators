/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include "config.h"

#include "Utilities.hpp"

#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/utility/cart_grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/cpgpreprocess/cgridinterface.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/SimpleFluid2p.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/transport/CSRMatrixUmfpackSolver.hpp>

#include <opm/core/transport/reorder/twophasetransport.hpp>

#include <boost/filesystem/convenience.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <tr1/array>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <vector>
#include <numeric>






class ReservoirState {
public:
    ReservoirState(const UnstructuredGrid* g, const int num_phases = 2)
        : press_ (g->number_of_cells, 0.0),
          fpress_(g->number_of_faces, 0.0),
          flux_  (g->number_of_faces, 0.0),
          sat_   (num_phases * g->number_of_cells, 0.0)
    {
	for (int cell = 0; cell < g->number_of_cells; ++cell) {
	    sat_[num_phases*cell + num_phases - 1] = 1.0;
	}
    }

    int numPhases() const { return sat_.size()/press_.size(); }

    ::std::vector<double>& pressure    () { return press_ ; }
    ::std::vector<double>& facepressure() { return fpress_; }
    ::std::vector<double>& faceflux    () { return flux_  ; }
    ::std::vector<double>& saturation  () { return sat_   ; }

    const ::std::vector<double>& pressure    () const { return press_ ; }
    const ::std::vector<double>& facepressure() const { return fpress_; }
    const ::std::vector<double>& faceflux    () const { return flux_  ; }
    const ::std::vector<double>& saturation  () const { return sat_   ; }

private:
    ::std::vector<double> press_ ;
    ::std::vector<double> fpress_;
    ::std::vector<double> flux_  ;
    ::std::vector<double> sat_   ;
};





template <class State>
void outputState(const UnstructuredGrid* grid,
		 const State& state,
		 const int step,
		 const std::string& output_dir)
{
    std::ostringstream vtkfilename;
    vtkfilename << output_dir << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
	THROW("Failed to open " << vtkfilename.str());
    }
    Opm::writeVtkDataGeneralGrid(grid, state.pressure(), state.saturation(), vtkfile);
}






// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    std::cout << "\n================    Test program for incompressible two-phase flow     ===============\n\n";
    Opm::parameter::ParameterGroup param(argc, argv, false);
    std::cout << "---------------    Reading parameters     ---------------" << std::endl;

    // Reading various control parameters.
    const int num_psteps = param.getDefault("num_psteps", 1);
    const double stepsize_days = param.getDefault("stepsize_days", 1.0);
    const double stepsize = Opm::unit::convert::from(stepsize_days, Opm::unit::day);
    const bool output = param.getDefault("output", true);
    std::string output_dir;
    if (output) {
	output_dir = param.getDefault("output_dir", std::string("output"));
	// Ensure that output dir exists
	boost::filesystem::path fpath(output_dir);
	create_directories(fpath);
    }

    // If we have a "deck_filename", grid and props will be read from that.
    bool use_deck = param.has("deck_filename");
    boost::scoped_ptr<Opm::Grid> grid;
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
    if (use_deck) {
	std::string deck_filename = param.get<std::string>("deck_filename");
	Opm::EclipseGridParser deck(deck_filename);
	// Grid init
	grid.reset(new Opm::Grid(deck));
	// Rock and fluid init
	const int* gc = grid->c_grid()->global_cell;
	std::vector<int> global_cell(gc, gc + grid->c_grid()->number_of_cells);
	props.reset(new Opm::IncompPropertiesFromDeck(deck, global_cell));
    } else {
	// Grid init.
	const int nx = param.getDefault("nx", 100);
	const int ny = param.getDefault("ny", 100);
	const int nz = param.getDefault("nz", 1);
	grid.reset(new Opm::Grid(nx, ny, nz));
	// Rock and fluid init.
	props.reset(new Opm::IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
    }

    // Extra rock init.
    std::vector<double> porevol;
    compute_porevolume(grid->c_grid(), *props, porevol);
    double tot_porevol = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Solvers init.
    Opm::PressureSolver psolver(grid->c_grid(), *props);

    // State-related and source-related variables init.
    std::vector<double> totmob;
    ReservoirState state(grid->c_grid(), props->numPhases());
    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(grid->c_grid()->number_of_cells);
    double flow_per_sec = 0.1*tot_porevol/Opm::unit::day;
    std::vector<double> src   (grid->c_grid()->number_of_cells, 0.0);
    src[0]                         =  flow_per_sec;
    src[grid->c_grid()->number_of_cells - 1] = -flow_per_sec;
    std::vector<double> reorder_src = src;

    // Control init.
    double current_time = 0.0;
    double total_time = stepsize*num_psteps;

    // Warn if any parameters are unused.
    if (param.anyUnused()) {
	std::cout << "--------------------   Unused parameters:   --------------------\n";
	param.displayUsage();
	std::cout << "----------------------------------------------------------------" << std::endl;
    }

    // Write parameters used for later reference.
    if (output) {
	param.writeParam(output_dir + "/spu_2p.param");
    }

    // Main simulation loop.
    std::cout << "\n\n================    Starting main simulation loop     ===============" << std::endl;
    for (int pstep = 0; pstep < num_psteps; ++pstep) {
        std::cout << "\n\n---------------    Simulation step number " << pstep
                  << "    ---------------"
                  << "\n      Current time (days)     " << Opm::unit::convert::to(current_time, Opm::unit::day)
                  << "\n      Current stepsize (days) " << Opm::unit::convert::to(stepsize, Opm::unit::day)
                  << "\n      Total time (days)       " << Opm::unit::convert::to(total_time, Opm::unit::day)
                  << "\n" << std::endl;

	if (output) {
	    outputState(grid->c_grid(), state, pstep, output_dir);
	}

	compute_totmob(*props, state.saturation(), totmob);
	psolver.solve(grid->c_grid(), totmob, src, state);

	Opm::toWaterSat(state.saturation(), reorder_sat);
	// We must treat reorder_src here,
	// if we are to handle anything but simple water
	// injection, since it is expected to be
	// equal to total outflow (if negative)
	// and water inflow (if positive).
	// Also, for anything but noflow boundaries,
	// boundary flows must be accumulated into
	// source term following the same convention.
	twophasetransport(&porevol[0],
			  &reorder_src[0],
			  stepsize,
			  const_cast<UnstructuredGrid*>(grid->c_grid()),
			  props.get(),
			  &state.faceflux()[0],
			  &reorder_sat[0]);
	Opm::toBothSat(reorder_sat, state.saturation());

	current_time += stepsize;
    }

    if (output) {
	outputState(grid->c_grid(), state, num_psteps, output_dir);
    }
}
