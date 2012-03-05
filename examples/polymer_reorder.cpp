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

#include <opm/core/pressure/IncompTpfa.hpp>

#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/fluid/IncompPropertiesBasic.hpp>
#include <opm/core/fluid/IncompPropertiesFromDeck.hpp>

#include <opm/core/linalg/LinearSolverUmfpack.hpp>
// #include <opm/core/linalg/LinearSolverIstl.hpp>

#include <opm/polymer/TransportModelPolymer.hpp>
#include <opm/polymer/PolymerProperties.hpp>

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



class AdHocProps : public Opm::IncompPropertiesBasic
{
public:
    AdHocProps(const Opm::parameter::ParameterGroup& param, int dim, int num_cells)
	: Opm::IncompPropertiesBasic(param, dim, num_cells)
    {
	ASSERT(numPhases() == 2);
	sw_.resize(3);
	sw_[0] = 0.2;
	sw_[1] = 0.7;
	sw_[2] = 1.0;
	krw_.resize(3);
	krw_[0] = 0.0;
	krw_[1] = 0.7;
	krw_[2] = 1.0;
	so_.resize(2);
	so_[0] = 0.3;
	so_[1] = 0.8;
	kro_.resize(2);
	kro_[0] = 0.0;
	kro_[1] = 1.0;
    }

    virtual void relperm(const int n,
			 const double* s,
			 const int* /*cells*/,
			 double* kr,
			 double* dkrds) const
    {
	// ASSERT(dkrds == 0);
	// We assume two phases flow
	for (int i = 0; i < n; ++i) {
	    kr[2*i] = krw(s[2*i]);
	    kr[2*i+1] = kro(s[2*i+1]);
	    if (dkrds != 0) {
		dkrds[2*i] = krw_dsw(s[2*i]);
		dkrds[2*i+3] = kro_dso(s[2*i+1]);
		dkrds[2*i+1] = -dkrds[2*i+3];
		dkrds[2*i+2] = -dkrds[2*i];
	    }
	}
    }


private:
    double krw(double s) const
    {
	return Opm::linearInterpolation(sw_, krw_, s);
    }

    double krw_dsw(double s) const
    {
	return Opm::linearInterpolationDerivative(sw_, krw_, s);
    }


    double kro(double s) const
    {
	return Opm::linearInterpolation(so_, kro_, s);
    }

    double kro_dso(double s) const
    {
	return Opm::linearInterpolationDerivative(so_, kro_, s);
    }

    std::vector<double> sw_;
    std::vector<double> krw_;
    std::vector<double> so_;
    std::vector<double> kro_;
};


class ReservoirState
{
public:
    ReservoirState(const UnstructuredGrid* g, const int num_phases = 2, const double init_sat = 0.0)
        : press_ (g->number_of_cells, 0.0),
          fpress_(g->number_of_faces, 0.0),
          flux_  (g->number_of_faces, 0.0),
          sat_   (num_phases * g->number_of_cells, 0.0),
	  concentration_(g->number_of_cells, 0.0),
	  cmax_(g->number_of_cells, 0.0)
    {
	for (int cell = 0; cell < g->number_of_cells; ++cell) {
	    sat_[num_phases*cell] = init_sat;
	    sat_[num_phases*cell + num_phases - 1] = 1.0 - init_sat;
	}
    }

    int numPhases() const { return sat_.size()/press_.size(); }

    std::vector<double>& pressure    () { return press_ ; }
    std::vector<double>& facepressure() { return fpress_; }
    std::vector<double>& faceflux    () { return flux_  ; }
    std::vector<double>& saturation  () { return sat_   ; }
    std::vector<double>& concentration() { return concentration_; }
    std::vector<double>& cmax() { return cmax_; }

    const std::vector<double>& pressure    () const { return press_ ; }
    const std::vector<double>& facepressure() const { return fpress_; }
    const std::vector<double>& faceflux    () const { return flux_  ; }
    const std::vector<double>& saturation  () const { return sat_   ; }
    const std::vector<double>& concentration() const { return concentration_; }
    const std::vector<double>& cmax() const { return cmax_; }

private:
    std::vector<double> press_ ;
    std::vector<double> fpress_;
    std::vector<double> flux_  ;
    std::vector<double> sat_   ;
    std::vector<double> concentration_;
    std::vector<double> cmax_;
};



class PolymerInflow
{
public:
    PolymerInflow(const double starttime,
		  const double endtime,
		  const double amount)
	: stime_(starttime), etime_(endtime), amount_(amount)
    {
    }
    double operator()(double time)
    {
	if (time >= stime_ && time < etime_) {
	    return amount_;
	} else {
	    return 0.0;
	}
    }
private:
    double stime_;
    double etime_;
    double amount_;
};




template <class State>
void outputState(const UnstructuredGrid* grid,
		 const State& state,
		 const int step,
		 const std::string& output_dir)
{
    // Write data in VTK format.
    std::ostringstream vtkfilename;
    vtkfilename << output_dir << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
    std::ofstream vtkfile(vtkfilename.str().c_str());
    if (!vtkfile) {
	THROW("Failed to open " << vtkfilename.str());
    }
    Opm::DataMap dm;
    dm["saturation"] = &state.saturation();
    dm["pressure"] = &state.pressure();
    dm["concentration"] = &state.concentration();
    Opm::writeVtkData(grid, dm, vtkfile);

    // Write data (not grid) in Matlab format
    for (Opm::DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
	std::ostringstream fname;
	fname << output_dir << "/" << it->first << "-" << std::setw(3) << std::setfill('0') << step << ".dat";
	std::ofstream file(fname.str().c_str());
	if (!file) {
	    THROW("Failed to open " << fname.str());
	}
	const std::vector<double>& d = *(it->second);
	std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
    }
}






// ----------------- Main program -----------------
int
main(int argc, char** argv)
{
    std::cout << "\n================    Test program for incompressible two-phase flow with polymer    ===============\n\n";
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
    boost::scoped_ptr<Opm::GridManager> grid;
    boost::scoped_ptr<Opm::IncompPropertiesInterface> props;
    Opm::PolymerProperties polydata;
    if (use_deck) {
	std::string deck_filename = param.get<std::string>("deck_filename");
	Opm::EclipseGridParser deck(deck_filename);
	polydata.readFromDeck(deck);
	// Grid init
	// grid.reset(new Opm::GridManager(deck));
	const int nx = param.getDefault("nx", 100);
	const int ny = param.getDefault("ny", 100);
	const int nz = param.getDefault("nz", 1);
	const double dx = param.getDefault("dx", 1.0);
	const double dy = param.getDefault("dy", 1.0);
	const double dz = param.getDefault("dz", 1.0);
	grid.reset(new Opm::GridManager(nx, ny, nz, dx, dy, dz));
	// Rock and fluid init
	const int* gc = grid->c_grid()->global_cell;
	std::vector<int> global_cell(gc, gc + grid->c_grid()->number_of_cells);
	props.reset(new Opm::IncompPropertiesFromDeck(deck, global_cell));
	// props.reset(new AdHocProps(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
    } else {
	// Grid init.
	const int nx = param.getDefault("nx", 100);
	const int ny = param.getDefault("ny", 100);
	const int nz = param.getDefault("nz", 1);
	const double dx = param.getDefault("dx", 1.0);
	const double dy = param.getDefault("dy", 1.0);
	const double dz = param.getDefault("dz", 1.0);
	grid.reset(new Opm::GridManager(nx, ny, nz, dx, dy, dz));
	// Rock and fluid init.
	// props.reset(new Opm::IncompPropertiesBasic(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
	props.reset(new AdHocProps(param, grid->c_grid()->dimensions, grid->c_grid()->number_of_cells));
	// Setting polydata defaults to mimic a simple example case.

	double c_max = param.getDefault("c_max_limit", 5.0);
	double mix_param = param.getDefault("mix_param", 1.0);
	double rock_density = param.getDefault("rock_density", 1000.0);
	double dead_pore_vol = param.getDefault("dead_pore_vol", 0.15);
	std::vector<double> c_vals_visc(2, -1e100);
	c_vals_visc[0] = 0.0;
	c_vals_visc[1] = 7.0;
	std::vector<double> visc_mult_vals(2, -1e100);
	visc_mult_vals[0] = 1.0;
	// polydata.visc_mult_vals[1] = param.getDefault("c_max_viscmult", 30.0);
	visc_mult_vals[1] = 20.0;
	std::vector<double> c_vals_ads(3, -1e100);
	c_vals_ads[0] = 0.0;
	c_vals_ads[1] = 2.0;
	c_vals_ads[2] = 8.0;
	std::vector<double> ads_vals(3, -1e100);
	ads_vals[0] = 0.0;
	// polydata.ads_vals[1] = param.getDefault("c_max_ads", 0.0025);
	ads_vals[1] = 0.0015;
	ads_vals[2] = 0.0025;
	polydata.set(c_max, mix_param, rock_density, dead_pore_vol, c_vals_visc, visc_mult_vals, c_vals_ads, ads_vals);
    }



    double poly_start = param.getDefault("poly_start_days", 300.0)*Opm::unit::day;
    double poly_end = param.getDefault("poly_end_days", 800.0)*Opm::unit::day;
    double poly_amount = param.getDefault("poly_amount", 5.0);
    PolymerInflow poly_inflow(poly_start, poly_end, poly_amount);

    // Extra rock init.
    std::vector<double> porevol;
    computePorevolume(*grid->c_grid(), *props, porevol);
    double tot_porevol = std::accumulate(porevol.begin(), porevol.end(), 0.0);

    // Gravity init.
    double gravity[3] = { 0.0 };
    double g = param.getDefault("gravity", 0.0);
    bool use_gravity = g != 0.0;
    if (use_gravity) {
	gravity[grid->c_grid()->dimensions - 1] = g;
	if (props->density()[0] == props->density()[1]) {
	    std::cout << "**** Warning: nonzero gravity, but zero density difference." << std::endl;
	}
    }

    // Solvers init.
    Opm::LinearSolverUmfpack linsolver;
    // Opm::LinearSolverIstl linsolver(param);
    const double *grav = use_gravity ? &gravity[0] : 0;
    Opm::IncompTpfa psolver(*grid->c_grid(), props->permeability(), grav, linsolver);

    Opm::TransportModelPolymer::SingleCellMethod method;
    std::string method_string = param.getDefault("single_cell_method", std::string("Bracketing"));
    if (method_string == "Bracketing") {
	method = Opm::TransportModelPolymer::Bracketing;
    } else if (method_string == "Newton") {
	method = Opm::TransportModelPolymer::Newton;
    } else {
	THROW("Unknown method: " << method_string);
    }
    const double nltol = param.getDefault("nl_tolerance", 1e-9);
    const int maxit = param.getDefault("nl_maxiter", 30);
    Opm::TransportModelPolymer tmodel(*grid->c_grid(), props->porosity(), &porevol[0], *props, polydata,
				      method, nltol, maxit);

    // State-related and source-related variables init.
    int num_cells = grid->c_grid()->number_of_cells;
    std::vector<double> totmob;
    std::vector<double> omega; // Empty dummy unless/until we include gravity here.
    double init_sat = param.getDefault("init_sat", 0.0);
    ReservoirState state(grid->c_grid(), props->numPhases(), init_sat);
    // We need a separate reorder_sat, because the reorder
    // code expects a scalar sw, not both sw and so.
    std::vector<double> reorder_sat(num_cells);
    double flow_per_sec = 0.1*tot_porevol/Opm::unit::day;
    if (param.has("injection_rate_per_day")) {
	flow_per_sec = param.get<double>("injection_rate_per_day")/Opm::unit::day;
    }
    std::vector<double> src(num_cells, 0.0);
    src[0]             =  flow_per_sec;
    src[num_cells - 1] = -flow_per_sec;
    std::vector<double> reorder_src = src;

    // Control init.
    double current_time = 0.0;
    double total_time = stepsize*num_psteps;

    // The allcells vector is used in calls to computeTotalMobility()
    // and computeTotalMobilityOmega().
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
	allcells[cell] = cell;
    }

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
    Opm::time::StopWatch pressure_timer;
    double ptime = 0.0;
    Opm::time::StopWatch transport_timer;
    double ttime = 0.0;
    Opm::time::StopWatch total_timer;
    total_timer.start();
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

	if (use_gravity) {
	    computeTotalMobilityOmega(*props, allcells, state.saturation(), totmob, omega);
	} else {
	    computeTotalMobility(*props, allcells, state.saturation(), totmob);
	}
	pressure_timer.start();
	psolver.solve(totmob, omega, src, state.pressure(), state.faceflux());
	pressure_timer.stop();
	double pt = pressure_timer.secsSinceStart();
	std::cout << "Pressure solver took:  " << pt << " seconds." << std::endl;
	ptime += pt;

	const double inflowc0 = poly_inflow(current_time + 1e-5*stepsize);
	const double inflowc1 = poly_inflow(current_time + (1.0 - 1e-5)*stepsize);
	if (inflowc0 != inflowc1) {
	    std::cout << "**** Warning: polymer inflow rate changes during timestep. Using rate near start of step.";
	}
	const double inflow_c = inflowc0;
	Opm::toWaterSat(state.saturation(), reorder_sat);
	// We must treat reorder_src here,
	// if we are to handle anything but simple water
	// injection, since it is expected to be
	// equal to total outflow (if negative)
	// and water inflow (if positive).
	// Also, for anything but noflow boundaries,
	// boundary flows must be accumulated into
	// source term following the same convention.
	transport_timer.start();
	tmodel.solve(&state.faceflux()[0], &reorder_src[0], stepsize, inflow_c,
		     &reorder_sat[0], &state.concentration()[0], &state.cmax()[0]);
	transport_timer.stop();
	double tt = transport_timer.secsSinceStart();
	std::cout << "Transport solver took: " << tt << " seconds." << std::endl;
	ttime += tt;
	Opm::toBothSat(reorder_sat, state.saturation());

	current_time += stepsize;
    }
    total_timer.stop();

    std::cout << "\n\n================    End of simulation     ===============\n"
	      << "Total time taken: " << total_timer.secsSinceStart()
	      << "\n  Pressure time:  " << ptime
	      << "\n  Transport time: " << ttime << std::endl;

    if (output) {
	outputState(grid->c_grid(), state, num_psteps, output_dir);
    }
}
