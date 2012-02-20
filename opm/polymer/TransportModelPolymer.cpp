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


#include <opm/polymer/TransportModelPolymer.hpp>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/RootFinders.hpp>
#include <cmath>

#define DEBUG

static double norm(double* res);

namespace Opm
{


    TransportModelPolymer::TransportModelPolymer(const UnstructuredGrid& grid,
						 const double* porosity,
						 const double* porevolume,
						 const IncompPropertiesInterface& props,
						 const PolymerData& polyprops,
						 const int method,
						 const double tol,
						 const int maxit)
	: grid_(grid),
	  porosity_(porosity),
	  porevolume_(porevolume),
	  props_(props),
	  polyprops_(polyprops),
	  tol_(tol),
	  maxit_(maxit),
	  darcyflux_(0),
	  source_(0),
	  dt_(0.0),
	  inflow_c_(0.0),
	  saturation_(0),
	  concentration_(0),
	  cmax_(0),
	  fractionalflow_(grid.number_of_cells, -1.0),
	  mc_(grid.number_of_cells, -1.0),
	  method_(method)
    {
	if (props.numPhases() != 2) {
	    THROW("Property object must have 2 phases");
	}
	visc_ = props.viscosity();

	// Set up smin_ and smax_
	int num_cells = props.numCells();
	smin_.resize(props.numPhases()*num_cells);
	smax_.resize(props.numPhases()*num_cells);
	std::vector<int> cells(num_cells);
	for (int i = 0; i < num_cells; ++i) {
	    cells[i] = i;
	}
	props.satRange(props.numCells(), &cells[0], &smin_[0], &smax_[0]);
    }



    void TransportModelPolymer::solve(const double* darcyflux,
				      const double* source,
				      const double dt,
				      const double inflow_c,
				      double* saturation,
				      double* concentration,
				      double* cmax)
    {
	darcyflux_ = darcyflux;
	source_ = source;
	dt_ = dt;
	inflow_c_ = inflow_c;
	saturation_ = saturation;
	concentration_ = concentration;
	cmax_ = cmax;
	reorderAndTransport(grid_, darcyflux);
    }




    // Residual for saturation equation, single-cell implicit Euler transport
    //
    //     r(s) = s - s0 + dt/pv*( influx + outflux*f(s) )
    //
    // where influx is water influx, outflux is total outflux.
    // Influxes are negative, outfluxes positive.
    struct TransportModelPolymer::ResidualS
    {
	const TransportModelPolymer& tm_;
	const int cell_;
	const double s0_;
	const double influx_;  // sum_j min(v_ij, 0)*f(s_j)
	const double outflux_; // sum_j max(v_ij, 0)
	const double dtpv_;    // dt/pv(i)
	const double c_;
	explicit ResidualS(const TransportModelPolymer& tmodel,
			   const int cell,
			   const double s0,
			   const double influx,
			   const double outflux,
			   const double dtpv,
			   const double c)
	    : tm_(tmodel),
	      cell_(cell),
	      s0_(s0),
	      influx_(influx),
	      outflux_(outflux),
	      dtpv_(dtpv),
	      c_(c)
	{
	}
	double operator()(double s) const
	{
	    return s - s0_ +  dtpv_*(outflux_*tm_.fracFlow(s, c_, cell_) + influx_);
	}
    };



    // Residual for concentration equation, single-cell implicit Euler transport
    //
    //  \TODO doc me
    // where ...
    // Influxes are negative, outfluxes positive.
    struct TransportModelPolymer::ResidualC
    {
	int cell;
	double s0;
	double c0;
	double cmax0;
	double influx;  // sum_j min(v_ij, 0)*f(s_j)
	double influx_polymer;  // sum_j min(v_ij, 0)*f(s_j)*mc(c_j)
	double outflux; // sum_j max(v_ij, 0)
	double porosity;
	double dtpv;    // dt/pv(i)
	mutable double s; // Mutable in order to change it with every operator() call to be the last computed s value.
	const TransportModelPolymer& tm;
	explicit ResidualC(const TransportModelPolymer& tmodel, int cell_index)
	    : tm(tmodel)
	{
	    cell    = cell_index;
	    s0      = tm.saturation_[cell];
	    c0      = tm.concentration_[cell];
	    cmax0   = tm.cmax_[cell];
            double dflux       = -tm.source_[cell];
            bool src_is_inflow = dflux < 0.0;
	    influx  =  src_is_inflow ? dflux : 0.0;
	    influx_polymer = src_is_inflow ? dflux*tm.computeMc(tm.inflow_c_) : 0.0;
	    outflux = !src_is_inflow ? dflux : 0.0;
	    dtpv    = tm.dt_/tm.porevolume_[cell];
	    porosity = tm.porosity_[cell];
	    s = -1e100;

	    for (int i = tm.grid_.cell_facepos[cell]; i < tm.grid_.cell_facepos[cell+1]; ++i) {
		int f = tm.grid_.cell_faces[i];
		double flux;
		int other;
		// Compute cell flux
		if (cell == tm.grid_.face_cells[2*f]) {
		    flux  = tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f+1];
		} else {
		    flux  =-tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f];
		}
		// Add flux to influx or outflux, if interior.
		if (other != -1) {
		    if (flux < 0.0) {
			influx  += flux*tm.fractionalflow_[other];
			influx_polymer += flux*tm.fractionalflow_[other]*tm.mc_[other];
		    } else {
			outflux += flux;
		    }
		}
	    }
	}
	double operator()(double c) const
	{
	    ResidualS res_s(tm, cell, s0, influx, outflux, dtpv, c);
	    int iters_used;
	    // Solve for s first.
	    s = modifiedRegulaFalsi(res_s, tm.smin_[2*cell], tm.smax_[2*cell], tm.maxit_, tm.tol_, iters_used);
	    double ff = tm.fracFlow(s, c, cell);
	    double mc = tm.computeMc(c);
	    double dps = tm.polyprops_.dps;
	    double rhor = tm.polyprops_.rhor;
	    double ads0 = tm.polyprops_.adsorbtion(std::max(c0, cmax0));
	    double ads = tm.polyprops_.adsorbtion(std::max(c, cmax0));
	    double res = (s - dps)*c - (s0 - dps)*c0
		+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
		+ dtpv*(outflux*ff*mc + influx_polymer);
#ifdef EXTRA_DEBUG_OUTPUT
	    std::cout << "c = " << c << "    s = " << s << "    c-residual = " << res << std::endl;
#endif
	    return res;
	}

	double lastSaturation() const
	{
	    return s;
	}
    };


    struct TransportModelPolymer::Residual
    {
	int cell;
	double s0;
	double c0;
	double cmax0;
	double influx;  // sum_j min(v_ij, 0)*f(s_j)
	double influx_polymer;  // sum_j min(v_ij, 0)*f(s_j)*mc(c_j)
	double outflux; // sum_j max(v_ij, 0)
	double porosity;
	double dtpv;    // dt/pv(i)
	const TransportModelPolymer& tm;

	Residual(const TransportModelPolymer& tmodel, int cell_index)
	    : tm(tmodel)
	{
	    cell    = cell_index;
	    s0      = tm.saturation_[cell];
	    c0      = tm.concentration_[cell];
	    cmax0   = tm.cmax_[cell];
            double dflux       = -tm.source_[cell];
            bool src_is_inflow = dflux < 0.0;
	    influx  =  src_is_inflow ? dflux : 0.0;
	    influx_polymer = src_is_inflow ? dflux*tm.computeMc(tm.inflow_c_) : 0.0;
	    outflux = !src_is_inflow ? dflux : 0.0;
	    dtpv    = tm.dt_/tm.porevolume_[cell];
	    porosity = tm.porosity_[cell];

	    for (int i = tm.grid_.cell_facepos[cell]; i < tm.grid_.cell_facepos[cell+1]; ++i) {
		int f = tm.grid_.cell_faces[i];
		double flux;
		int other;
		// Compute cell flux
		if (cell == tm.grid_.face_cells[2*f]) {
		    flux  = tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f+1];
		} else {
		    flux  =-tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f];
		}
		// Add flux to influx or outflux, if interior.
		if (other != -1) {
		    if (flux < 0.0) {
			influx  += flux*tm.fractionalflow_[other];
			influx_polymer += flux*tm.fractionalflow_[other]*tm.mc_[other];
		    } else {
			outflux += flux;
		    }
		}
	    }
	}

	void computeResidual(const double* x, double* res) const
	{
	    double s = x[0];
	    double c = x[1];
	    double ff = tm.fracFlow(s, c, cell);
	    double mc = tm.computeMc(c);
	    double dps = tm.polyprops_.dps;
	    double rhor = tm.polyprops_.rhor;
	    double ads0 = tm.polyprops_.adsorbtion(std::max(c0, cmax0));
	    double ads = tm.polyprops_.adsorbtion(std::max(c, cmax0));
	    res[0] = s - s0 +  dtpv*(outflux*tm.fracFlow(s, c, cell) + influx);
	    res[1] = (s - dps)*c - (s0 - dps)*c0
		+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
		+ dtpv*(outflux*ff*mc + influx_polymer);
	}

	void computeGradient(const double* x, double* res, double* gradient, bool s_or_c) const
	// If s_or_c == true, compute the gradient of s-residual, if s_or_c == false, compute gradient of c-residual
	{
	    double s = x[0];
	    double c = x[1];
	    double ff_ds_dc[2];
	    double ff = tm.fracFlowWithDer(s, c, cell, ff_ds_dc);
	    double mc_dc;
	    double mc = tm.computeMcWithDer(c, &mc_dc);
	    double dps = tm.polyprops_.dps;
	    double rhor = tm.polyprops_.rhor;
	    double ads0 = tm.polyprops_.adsorbtion(std::max(c0, cmax0));
	    double ads;
	    double ads_dc;
	    if (c < cmax0) {
		ads = tm.polyprops_.adsorbtion(cmax0);
		ads_dc = 0;
	    } else {
		ads = tm.polyprops_.adsorbtionWithDer(c, &ads_dc);
	    }
	    res[0] = s - s0 +  dtpv*(outflux*tm.fracFlow(s, c, cell) + influx);
	    res[1] = (s - dps)*c - (s0 - dps)*c0
		+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
		+ dtpv*(outflux*ff*mc + influx_polymer);
	    if (s_or_c == true) {
		gradient[0] = 1 + dtpv*outflux*ff_ds_dc[0];
		gradient[1] = dtpv*outflux*ff_ds_dc[1];
	    } else if (s_or_c == false) {
		gradient[0] = c + dtpv*outflux*(ff_ds_dc[0])*mc;
		gradient[1] = s - dps + rhor*((1.0 - porosity)/porosity)*(ads_dc - ads0)
		    + dtpv*outflux*(ff_ds_dc[1]*mc + ff*mc_dc);
	    }
	}
    };



    struct TransportModelPolymer::ResidualSDir
    {
	int cell;
	double s0;
	double c0;
	double cmax0;
	double influx;  // sum_j min(v_ij, 0)*f(s_j)
	double influx_polymer;  // sum_j min(v_ij, 0)*f(s_j)*mc(c_j)
	double outflux; // sum_j max(v_ij, 0)
	double porosity;
	double dtpv;    // dt/pv(i)
	double direction[2];
	double x[2];
	const TransportModelPolymer& tm;

	ResidualSDir(const TransportModelPolymer& tmodel, int cell_index)
	    : tm(tmodel)
	{
	    cell    = cell_index;
	    s0      = tm.saturation_[cell];
	    c0      = tm.concentration_[cell];
	    cmax0   = tm.cmax_[cell];
            double dflux       = -tm.source_[cell];
            bool src_is_inflow = dflux < 0.0;
	    influx  =  src_is_inflow ? dflux : 0.0;
	    influx_polymer = src_is_inflow ? dflux*tm.computeMc(tm.inflow_c_) : 0.0;
	    outflux = !src_is_inflow ? dflux : 0.0;
	    dtpv    = tm.dt_/tm.porevolume_[cell];
	    porosity = tm.porosity_[cell];

	    for (int i = tm.grid_.cell_facepos[cell]; i < tm.grid_.cell_facepos[cell+1]; ++i) {
		int f = tm.grid_.cell_faces[i];
		double flux;
		int other;
		// Compute cell flux
		if (cell == tm.grid_.face_cells[2*f]) {
		    flux  = tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f+1];
		} else {
		    flux  =-tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f];
		}
		// Add flux to influx or outflux, if interior.
		if (other != -1) {
		    if (flux < 0.0) {
			influx  += flux*tm.fractionalflow_[other];
			influx_polymer += flux*tm.fractionalflow_[other]*tm.mc_[other];
		    } else {
			outflux += flux;
		    }
		}
	    }
	}

	void setup(const double* x_arg, const double* direction_arg) {
	    x[0] = x_arg[0];
	    x[1] = x_arg[1];
	    direction[0] = direction_arg[0];
	    direction[1] = direction_arg[1];
	}

	double operator()(double t) const
	{
	    double s = x[0] + t*direction[0];
	    double c = x[1] + t*direction[1];
	    return s - s0 +  dtpv*(outflux*tm.fracFlow(s, c, cell) + influx);
	}
    };

    struct TransportModelPolymer::ResidualCDir
    {
	int cell;
	double s0;
	double c0;
	double cmax0;
	double influx;  // sum_j min(v_ij, 0)*f(s_j)
	double influx_polymer;  // sum_j min(v_ij, 0)*f(s_j)*mc(c_j)
	double outflux; // sum_j max(v_ij, 0)
	double porosity;
	double dtpv;    // dt/pv(i)
	double direction[2];
	double x[2];
	const TransportModelPolymer& tm;

	ResidualCDir(const TransportModelPolymer& tmodel, int cell_index)
	    : tm(tmodel)
	{
	    cell    = cell_index;
	    s0      = tm.saturation_[cell];
	    c0      = tm.concentration_[cell];
	    cmax0   = tm.cmax_[cell];
            double dflux       = -tm.source_[cell];
            bool src_is_inflow = dflux < 0.0;
	    influx  =  src_is_inflow ? dflux : 0.0;
	    influx_polymer = src_is_inflow ? dflux*tm.computeMc(tm.inflow_c_) : 0.0;
	    outflux = !src_is_inflow ? dflux : 0.0;
	    dtpv    = tm.dt_/tm.porevolume_[cell];
	    porosity = tm.porosity_[cell];

	    for (int i = tm.grid_.cell_facepos[cell]; i < tm.grid_.cell_facepos[cell+1]; ++i) {
		int f = tm.grid_.cell_faces[i];
		double flux;
		int other;
		// Compute cell flux
		if (cell == tm.grid_.face_cells[2*f]) {
		    flux  = tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f+1];
		} else {
		    flux  =-tm.darcyflux_[f];
		    other = tm.grid_.face_cells[2*f];
		}
		// Add flux to influx or outflux, if interior.
		if (other != -1) {
		    if (flux < 0.0) {
			influx  += flux*tm.fractionalflow_[other];
			influx_polymer += flux*tm.fractionalflow_[other]*tm.mc_[other];
		    } else {
			outflux += flux;
		    }
		}
	    }
	}

	void setup(const double* x_arg, const double* direction_arg) {
	    x[0] = x_arg[0];
	    x[1] = x_arg[1];
	    direction[0] = direction_arg[0];
	    direction[1] = direction_arg[1];
	}

	double operator()(double t) const
	{
	    double s = x[0] + t*direction[0];
	    double c = x[1] + t*direction[1];
	    double ff = tm.fracFlow(s, c, cell);
	    double mc = tm.computeMc(c);
	    double dps = tm.polyprops_.dps;
	    double rhor = tm.polyprops_.rhor;
	    double ads0 = tm.polyprops_.adsorbtion(std::max(c0, cmax0));
	    double ads = tm.polyprops_.adsorbtion(std::max(c, cmax0));
	    return (s - dps)*c - (s0 - dps)*c0
		+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
		+ dtpv*(outflux*ff*mc + influx_polymer);
	}
    };


    void TransportModelPolymer::solveSingleCell(const int cell)
    {
	if (method_ == 1) {
	    solveSingleCellBracketing(cell);
	} else if (method_ == 2) {
	    solveSingleCellSplitting(cell);
	} else {
	    THROW("Method is " << method_ << "");
	}
    }

    void TransportModelPolymer::solveSingleCellBracketing(int cell)
    {
	    ResidualC res(*this, cell);
	    const double a = 0.0;
	    const double b = polyprops_.c_max_limit;
	    int iters_used;
	    concentration_[cell] = modifiedRegulaFalsi(res, a, b, maxit_, tol_, iters_used);
	    cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
	    saturation_[cell] = res.lastSaturation();
	    fractionalflow_[cell] = fracFlow(saturation_[cell], concentration_[cell], cell);
	    mc_[cell] = computeMc(concentration_[cell]);
    }

    void TransportModelPolymer::solveSingleCellSplitting(int cell)
    {
	const int max_iters_falsi = 20;
	const double tol = 1e-9;
	int iters_used_falsi = 0;
	const int max_iters_split = 20;
	int iters_used_split = 0;

	Residual residual(*this, cell);
	ResidualSDir residual_s_dir(*this, cell);
	ResidualCDir residual_c_dir(*this, cell);
	double x[2] = {saturation_[cell], concentration_[cell]};
	double res[2];
	residual.computeResidual(x, res);

	if (norm(res) < tol) {
	    return;
#ifdef DEBUG
	    std::cout << "short escape" << std::endl;
#endif
	}

	bool use_zero_search = true;
	bool res_s_done;
	double x_min[2] = {0.0, 0.0};
	double x_max[2] = {1.0, polyprops_.c_max_limit};
	double t_min;
	double t_max;
	double t;
	double res_s_t_min;
	double res_s_t_max;
	double res_c_t_min;
	double res_c_t_max;
	double direction[2];
	double gradient[2];

#ifdef DEBUG
	std::cout << "Initial step" << std::endl;
#endif

	if (std::abs(res[0]) < std::abs(res[1])) {
	    // solve for s-residual in a 45 degree diagonal direction (down right)
	    direction[0] = 1;
	    direction[1] = -1;
	    residual_s_dir.setup(x, direction);
	    if (res[0] < 0) {
		t_min = 0.;
		res_s_t_min = res[0];
		t_max = std::min((x_max[0]-x[0])/direction[0], (x_min[1]-x[1])/direction[1]);
		res_s_t_max = residual_s_dir(t_max);
		if (res_s_t_max*res_s_t_min >= 0) {
		    use_zero_search = false;
		    t = t_max;
		}
	    } else {
		t_max = 0.;
		res_s_t_max = res[0];
		t_min = -std::min((x[0]-x_min[0])/direction[0], (x[1]-x_max[1])/direction[1]);
		res_s_t_min = residual_s_dir(t_min);
		if (res_s_t_max*res_s_t_min >= 0) {
		    use_zero_search = false;
		    t = t_min;
		}
	    }
	    if (use_zero_search) {
		t = modifiedRegulaFalsi(residual_s_dir, t_min, t_max, max_iters_falsi, tol, iters_used_falsi);
	    }
	    x[0] +=  t*direction[0];
	    x[1] +=  t*direction[1];
	    res_s_done = true;
	    residual.computeGradient(x, res, gradient, true);
	} else {
	    // solve for c-residual in 45 degree diagonal direction (up-right)
	    direction[0] = 1.;
	    direction[1] = 1.;
	    residual_c_dir.setup(x, direction);
	    if (res[1] < 0) {
		t_min = 0.;
		res_c_t_min = res[1];
		t_max = std::min((x_max[0]-x[0])/direction[0], (x_max[1]-x[1])/direction[1]);
		res_c_t_max = residual_c_dir(t_max);
		if (res_c_t_max*res_c_t_min >= 0) {
		    use_zero_search = false;
		    t = t_max;
		}
	    } else {
		t_max = 0;
		res_c_t_max = res[1];
		t_min = -std::min((x[0]-x_min[0])/direction[0], (x[1]-x_min[1])/direction[1]);
		res_c_t_min = residual_c_dir(t_min);
		if (res_c_t_max*res_c_t_min >= 0) {
		    use_zero_search = false;
		    t = t_min;
		}
	    }
	    if (use_zero_search) {
		t = modifiedRegulaFalsi(residual_c_dir, t_min, t_max, max_iters_falsi, tol, iters_used_falsi);
	    }
	    x[0] +=  t*direction[0];
	    x[1] +=  t*direction[1];
	    res_s_done = false;
	    residual.computeGradient(x, res, gradient, false);
	}

#ifdef DEBUG
	std::cout << "s: " << x[0] << std::endl;
	std::cout << "c: " << x[1] << std::endl;
	std::cout << "res[0]: " << res[0] << std::endl;
	std::cout << "res[1]: " << res[1] << std::endl;
	std::cout << "res_s_done" << res_s_done << std::endl;
	std::cout << "gradient[0]: " << gradient[0] << std::endl;
	std::cout << "gradient[1]: " << gradient[1] << std::endl;
#endif


	while ((norm(res) > tol) && (iters_used_split < max_iters_split)) {
	    use_zero_search = true;
	    if (res_s_done) { // solve for c-residual
		direction[0] = -gradient[1];
		direction[1] = gradient[0];
		residual_c_dir.setup(x, direction);
		if (res[1] < 0) {
		    t_min = 0.;
		    res_c_t_min = res[1];
		    t_max = std::min((x_max[0]-x[0])/direction[0], (x_max[1]-x[1])/direction[1]);
		    residual_c_dir(t_max);
		    if (res_c_t_max*res_c_t_min >= 0) {
			use_zero_search = false;
			t = t_max;
		    }
		} else {
		    t_max = 0;
		    res_c_t_max = res[1];
		    t_min = -std::min((x[0]-x_min[0])/direction[0], (x[1]-x_min[1])/direction[1]);
		    res_c_t_min = residual_c_dir(t_min);
		    if (res_c_t_max*res_c_t_min >= 0) {
			use_zero_search = false;
			t = t_min;
		    }
		}
		if (use_zero_search) {
		    t = modifiedRegulaFalsi(residual_c_dir, t_min, t_max, max_iters_falsi, tol, iters_used_falsi);
		}
		x[0] +=  t*direction[0];
		x[1] +=  t*direction[1];
		res_s_done = false;
		residual.computeGradient(x, res, gradient, false);
	    } else { // solve for s residual
		use_zero_search = true;
		direction[0] = gradient[1];
		direction[1] = -gradient[0];
		residual_s_dir.setup(x, direction);
		if (res[0] < 0) {
		    t_min = 0.;
		    res_s_t_min = res[0];
		    t_max = std::min((x_max[0]-x[0])/direction[0], (x_min[1]-x[1])/direction[1]);
		    res_s_t_max = residual_s_dir(t_max);
		    if (res_s_t_max*res_s_t_min >= 0) {
			use_zero_search = false;
			t = t_max;
		    }
		} else {
		    t_max = 0.;
		    res_s_t_max = res[0];
		    t_min = -std::min((x[0]-x_min[0])/direction[0], (x[1]-x_max[1])/direction[1]);
		    res_s_t_min = residual_s_dir(t_min);
		    if (res_s_t_max*res_s_t_min >= 0) {
			use_zero_search = false;
			t = t_min;
		    }
		}
		if (use_zero_search) {
		    t = modifiedRegulaFalsi(residual_s_dir, t_min, t_max, max_iters_falsi, tol, iters_used_falsi);
		}
		x[0] +=  t*direction[0];
		x[1] +=  t*direction[1];
		res_s_done = true;
		residual.computeGradient(x, res, gradient, true);
	    }
#ifdef DEBUG
	    std::cout << "s: " << x[0] << std::endl;
	    std::cout << "c: " << x[1] << std::endl;
	    std::cout << "res[0]: " << res[0] << std::endl;
	    std::cout << "res[1]: " << res[1] << std::endl;
	    std::cout << "res_s_done" << res_s_done << std::endl;
	    std::cout << "gradient[0]: " << gradient[0] << std::endl;
	    std::cout << "gradient[1]: " << gradient[1] << std::endl;
#endif

	    iters_used_split += 1;
	}

	if ((iters_used_split >=  max_iters_split) && (norm(res) >= tol)) {
	    solveSingleCellBracketing(cell);
	    std::cout << "splitting did not work" << std::endl;
	} else {
	    concentration_[cell] = x[1];
	    cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
	    saturation_[cell] = x[0];
	    fractionalflow_[cell] = fracFlow(saturation_[cell], concentration_[cell], cell);
	    mc_[cell] = computeMc(concentration_[cell]);
	}
    }

    void TransportModelPolymer::solveMultiCell(const int num_cells, const int* /*cells*/)
    {
	THROW("TransportModelPolymer::solveMultiCell() not yet implemented, "
	      "got a component of size " << num_cells);
    }




    double TransportModelPolymer::fracFlow(double s, double c, int cell) const
    {
	double c_max_limit = polyprops_.c_max_limit;
	double cbar = c/c_max_limit;
	double mu_w = visc_[0];
	double mu_m = polyprops_.viscMult(c)*mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.c_max_limit)*mu_w;
	double omega = polyprops_.omega;
	double mu_m_omega = std::pow(mu_m, omega);
	double mu_w_e   = mu_m_omega*std::pow(mu_w, 1.0 - omega);
	double mu_p_eff = mu_m_omega*std::pow(mu_p, 1.0 - omega);
	double inv_mu_w_eff = (1.0 - cbar)/mu_w_e + cbar/mu_p_eff;
	double inv_visc_eff[2] = { inv_mu_w_eff, 1.0/visc_[1] };
	double sat[2] = { s, 1.0 - s };
	double mob[2];
	props_.relperm(1, sat, &cell, mob, 0);
	mob[0] *= inv_visc_eff[0];
	mob[1] *= inv_visc_eff[1];
	return mob[0]/(mob[0] + mob[1]);
    }

    double TransportModelPolymer::fracFlowWithDer(double s, double c, int cell, double* der) const
    {
	// We should check the dimension of der
	double c_max_limit = polyprops_.c_max_limit;
	double cbar = c/c_max_limit;
	double mu_w = visc_[0];
	double mu_m_dc; // derivative of mu_m with respect to c
	double mu_m = polyprops_.viscMultWithDer(c, &mu_m_dc)*mu_w;
	mu_m_dc *= mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.c_max_limit)*mu_w;
	double omega = polyprops_.omega;
	double mu_m_omega = std::pow(mu_m, omega);
	double mu_m_omega_minus1 = std::pow(mu_m, omega - 1.0);
	double mu_w_omega = std::pow(mu_w, 1.0 - omega);
	double mu_w_e   = mu_m_omega*mu_w_omega;
	double mu_w_e_dc = omega*mu_m_dc*mu_m_omega_minus1*mu_w_omega;
	double mu_p_omega = std::pow(mu_p, 1.0 - omega);
	double mu_p_eff = mu_m_omega*mu_p_omega;
	double mu_p_eff_dc = omega*mu_m_dc*mu_m_omega_minus1*mu_p_omega;
	double mu_w_eff = 1./((1 - cbar)/mu_w_e + cbar/mu_p_eff);
	double inv_mu_w_eff_dc = -mu_w_e_dc/(mu_w_e*mu_w_e)*(1. - cbar) - mu_p_eff_dc/(mu_p_eff*mu_p_eff)*cbar + (1./mu_p_eff - 1./mu_w_e);
	double mu_w_eff_dc = -mu_w_eff*mu_w_eff*inv_mu_w_eff_dc;
	double visc_eff[2] = { mu_w_eff, visc_[1] };
	double sat[2] = { s, 1.0 - s };
	double mob[2];
	double mob_ds[2];
	double mob_dc[2];
	double perm[2];
	double perm_ds[4];
	props_.relperm(1, sat, &cell, perm, perm_ds);
	mob[0] = perm[0]/visc_eff[0];
	mob[1] = perm[1]/visc_eff[1];
	mob_ds[0] = perm_ds[0]/mu_w_eff;
	mob_ds[1] = perm_ds[1]/mu_w_eff;
	mob_dc[0] = - perm[0]*mu_w_eff_dc/(mu_w_eff*mu_w_eff);
	mob_dc[1] = - perm[1]*mu_p_eff_dc/(mu_p_eff*mu_p_eff);
	der[0] = (mob_ds[0]*mob[1] - mob_ds[1]*mob[0])/((mob[0] + mob[1])*(mob[0] + mob[1]));
	der[1] = (mob_dc[0]*mob[1] - mob_dc[1]*mob[0])/((mob[0] + mob[1])*(mob[0] + mob[1]));
 	return mob[0]/(mob[0] + mob[1]);
    }


    double TransportModelPolymer::computeMc(double c) const
    {
	double c_max_limit = polyprops_.c_max_limit;
	double cbar = c/c_max_limit;
	double mu_w = visc_[0];
	double mu_m = polyprops_.viscMult(c)*mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.c_max_limit)*mu_w;
	double omega = polyprops_.omega;
	double mu_m_omega = std::pow(mu_m, omega);
	double mu_w_e   = mu_m_omega*std::pow(mu_w, 1.0 - omega);
	double mu_p_eff = mu_m_omega*std::pow(mu_p, 1.0 - omega);
	double inv_mu_w_eff = (1.0 - cbar)/mu_w_e + cbar/mu_p_eff;
	return c/(inv_mu_w_eff*mu_p_eff);
    }

    double TransportModelPolymer::computeMcWithDer(double c, double* der) const
    {
	double c_max_limit = polyprops_.c_max_limit;
	double cbar = c/c_max_limit;
	double mu_w = visc_[0];
	double mu_m_dc; // derivative of mu_m with respect to c
	double mu_m = polyprops_.viscMultWithDer(c, &mu_m_dc)*mu_w;
	mu_m_dc *= mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.c_max_limit)*mu_w;
	double omega = polyprops_.omega;
	double mu_m_omega = std::pow(mu_m, omega);
	double mu_m_omega_minus1 = std::pow(mu_m, omega-1);
	double mu_w_omega = std::pow(mu_w, 1.0 - omega);
	double mu_w_e   = mu_m_omega*mu_w_omega;
	double mu_w_e_dc = omega*mu_m_dc*mu_m_omega_minus1*mu_w_omega;
	double mu_p_omega = std::pow(mu_p, 1.0 - omega);
	double mu_p_eff = mu_m_omega*mu_p_omega;
	double mu_p_eff_dc = omega*mu_m_dc*mu_m_omega_minus1*mu_p_omega;
	double mu_w_eff = 1./((1 - cbar)/mu_w_e + cbar/mu_p_eff);
	double inv_mu_w_eff_dc = -mu_w_e_dc/(mu_w_e*mu_w_e)*(1. - cbar) - mu_p_eff_dc/(mu_p_eff*mu_p_eff)*cbar + (1./mu_p_eff - 1./mu_w_e);
	double mu_w_eff_dc = -mu_w_eff*mu_w_eff*inv_mu_w_eff_dc;
	*der = mu_w_eff/mu_p_eff + c*mu_w_eff_dc/mu_p_eff - c*mu_p_eff_dc*mu_w_eff/(mu_p_eff*mu_p_eff);
	return c*mu_w_eff/mu_p_eff;
    }

} // namespace Opm


static double norm(double* res) {
    double absres0 = std::abs(res[0]);
    double absres1 = std::abs(res[1]);
    if (absres0 <= absres1) {
	return absres1;
    }
    else {
	return absres0;
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
