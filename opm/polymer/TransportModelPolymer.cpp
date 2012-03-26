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


class Opm::TransportModelPolymer::ResidualEquation
{
public:
    int cell;
    double s0;
    double c0;
    double cmax0;
    double influx;  // sum_j min(v_ij, 0)*f(s_j)
    double influx_polymer;  // sum_j min(v_ij, 0)*f(s_j)*mc(c_j)
    double outflux; // sum_j max(v_ij, 0)
    double porosity;
    double dtpv;    // dt/pv(i)
    double dps;
    double res_factor;
    double c_max_ads;
    double rhor;
    double ads0;
    int gradient_method;
    const TransportModelPolymer& tm;

    ResidualEquation(const TransportModelPolymer& tmodel, int cell_index);
    void computeResidual(const double* x, double* res) const;
    void computeResidual(const double* x, double* res, double& mc, double& ff) const;
    double computeResidualS(const double* x) const;
    double computeResidualC(const double* x) const;
    void computeGradientResS(const double* x, double* res, double* gradient) const;
    void computeGradientResC(const double* x, double* res, double* gradient) const;
    void computeJacobiRes(const double* x, double* res_s_ds_dc, double* res_c_ds_dc) const;

};


namespace
{
    double norm(double* res)
    {
	return std::max(std::abs(res[0]), std::abs(res[1]));
    }

    bool solveNewtonStep(const double* , const Opm::TransportModelPolymer::ResidualEquation&,
			 const double*, double*);


    // Define a piecewise linear curve along which we will look for zero of the "s" or "r" residual.
    // The curve starts at "x", goes along the direction "direction" until it hits the boundary of the box of
    // admissible values for "s" and "x" (which is given by "[x_min[0], x_max[0]]x[x_min[1], x_max[1]]").
    // Then it joins in a straight line the point "end_point".
    class CurveInSCPlane{
    public:
	CurveInSCPlane();
	void setup(const double* x, const double* direction,
		   const double* end_point, const double* x_min,
		   const double* x_max, double& t_max_out,
		   double& t_out_out);
	void computeXOfT(double*, const double) const;

    private:
	double direction_[2];
	double end_point_[2];
	double x_max_[2];
	double x_min_[2];
	double t_out_;
	double t_max_; // t_max = t_out + 1
	double x_out_[2];
	double x_[2];
    };


    // Compute the "s" residual along the curve "curve" for a given residual equation "res_eq".
    // The operator() is sent to a root solver.
    class ResSOnCurve
    {
    public:
	ResSOnCurve(const Opm::TransportModelPolymer::ResidualEquation& res_eq);
	double operator()(const double t) const;
	CurveInSCPlane curve;
    private:
	Opm::TransportModelPolymer::ResidualEquation res_eq_;
    };

    // Compute the "c" residual along the curve "curve" for a given residual equation "res_eq".
    // The operator() is sent to a root solver.
    class ResCOnCurve
    {
    public:
	ResCOnCurve(const Opm::TransportModelPolymer::ResidualEquation& res_eq);
	double operator()(const double t) const;
	CurveInSCPlane curve;
    private:
	Opm::TransportModelPolymer::ResidualEquation res_eq_;
    };

}


namespace Opm
{
    TransportModelPolymer::TransportModelPolymer(const UnstructuredGrid& grid,
						 const double* porosity,
						 const double* porevolume,
						 const IncompPropertiesInterface& props,
						 const PolymerProperties& polyprops,
						 const SingleCellMethod method,
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
	const double cmax0_;
	const double influx_;  // sum_j min(v_ij, 0)*f(s_j)
	const double outflux_; // sum_j max(v_ij, 0)
	const double dtpv_;    // dt/pv(i)
	const double c_;
	explicit ResidualS(const TransportModelPolymer& tmodel,
			   const int cell,
			   const double s0,
			   const double cmax0,
			   const double influx,
			   const double outflux,
			   const double dtpv,
			   const double c)
	    : tm_(tmodel),
	      cell_(cell),
	      s0_(s0),
	      cmax0_(cmax0),
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

	void computeBothResiduals(const double s_arg, const double c_arg, double& res_s, double& res_c, double& mc, double& ff) const
	{
	    double dps = tm.polyprops_.deadPoreVol();
	    ff = tm.fracFlow(s_arg, c_arg, cell);
	    mc = tm.computeMc(c_arg);
	    double rhor = tm.polyprops_.rockDensity();
	    double ads0 = tm.polyprops_.adsorbtion(c0, cmax0);
	    double ads = tm.polyprops_.adsorbtion(c_arg, cmax0);
	    res_s =  s_arg - s0 +  dtpv*(outflux*ff + influx);
	    res_c = s_arg*(1 - dps)*c_arg - (s0 - dps)*c0
		+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
		+ dtpv*(outflux*ff*mc + influx_polymer);

	}

	double operator()(double c) const
	{
	    double dps = tm.polyprops_.deadPoreVol();
	    ResidualS res_s(tm, cell, s0, cmax0, influx, outflux, dtpv, c);
	    int iters_used;
	    // Solve for s first.
	    s = modifiedRegulaFalsi(res_s, std::max(tm.smin_[2*cell], dps), tm.smax_[2*cell],
				    tm.maxit_, tm.tol_, iters_used);
	    double ff = tm.fracFlow(s, c, cell);
	    double mc = tm.computeMc(c);
	    double rhor = tm.polyprops_.rockDensity();
	    double ads0 = tm.polyprops_.adsorbtion(c0, cmax0);
	    double ads = tm.polyprops_.adsorbtion(c, cmax0);
	    double res = (1 - dps)*s*c - (1 - dps)*s0*c0
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


    // ResidualEquation gathers parameters to construct the residual, computes its
    // value and the values of its derivatives.

    TransportModelPolymer::ResidualEquation::ResidualEquation(const TransportModelPolymer& tmodel, int cell_index)
	: tm(tmodel)
    {
	gradient_method = 1;
	cell    = cell_index;
	s0      = tm.saturation_[cell];
	c0      = tm.concentration_[cell];
	cmax0   = tm.cmax_[cell];
	dps = tm.polyprops_.deadPoreVol();
	rhor = tm.polyprops_.rockDensity();
	ads0 = tm.polyprops_.adsorbtion(c0, cmax0);
	res_factor = tm.polyprops_.resFactor();
	c_max_ads = tm.polyprops_.cMaxAds();
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

    void TransportModelPolymer::ResidualEquation::computeResidual(const double* x, double* res) const
    {
	double s = x[0];
	double c = x[1];
	double ff = tm.fracFlow(s, c, cell);
	double mc = tm.computeMc(c);
	double dps = tm.polyprops_.deadPoreVol();
	double rhor = tm.polyprops_.rockDensity();
	double ads0 = tm.polyprops_.adsorbtion(c0, cmax0);
	double ads = tm.polyprops_.adsorbtion(c, cmax0);
	res[0] = s - s0 +  dtpv*(outflux*ff + influx);
	res[1] = (1 - dps)*s*c - (1 - dps)*s0*c0
	    + rhor*((1.0 - porosity)/porosity)*(ads - ads0)
	    + dtpv*(outflux*ff*mc + influx_polymer);
    }

    void TransportModelPolymer::ResidualEquation::computeResidual(const double* x, double* res, double& mc, double& ff) const
    {
	double s = x[0];
	double c = x[1];
	ff = tm.fracFlow(s, c, cell);
	mc = tm.computeMc(c);
	double dps = tm.polyprops_.deadPoreVol();
	double rhor = tm.polyprops_.rockDensity();
	double ads0 = tm.polyprops_.adsorbtion(c0, cmax0);
	double ads = tm.polyprops_.adsorbtion(c, cmax0);
	res[0] = s - s0 +  dtpv*(outflux*ff + influx);
	res[1] = (1 - dps)*s*c - (1 - dps)*s0*c0
	    + rhor*((1.0 - porosity)/porosity)*(ads - ads0)
	    + dtpv*(outflux*ff*mc + influx_polymer);
    }


    double TransportModelPolymer::ResidualEquation::computeResidualS(const double* x) const
    {
	double s = x[0];
	double c = x[1];
	double ff = tm.fracFlow(s, c, cell);
	return s - s0 +  dtpv*(outflux*ff + influx);
    }

    double TransportModelPolymer::ResidualEquation::computeResidualC(const double* x) const
    {
	double s = x[0];
	double c = x[1];
	double ff = tm.fracFlow(s, c, cell);
	double mc = tm.computeMc(c);
	double ads = tm.polyprops_.adsorbtion(c, cmax0);
	return (1 - dps)*s*c - (1 - dps)*s0*c0
	    + rhor*((1.0 - porosity)/porosity)*(ads - ads0)
	    + dtpv*(outflux*ff*mc + influx_polymer);
    }

    void TransportModelPolymer::ResidualEquation::computeGradientResS(const double* x, double* res, double* gradient) const
    // If gradient_method == 1, use finite difference
    // If gradient_method == 2, use analytic expresions
    {
	if (gradient_method == 1) {
	    double epsi = 1e-8;
	    double res_epsi[2];
	    double x_epsi[2];
	    computeResidual(x, res);
	    x_epsi[0] = x[0] + epsi;
	    x_epsi[1] = x[1];
	    computeResidual(x_epsi, res_epsi);
	    gradient[0] =  (res_epsi[0] - res[0])/epsi;
	    x_epsi[0] = x[0];
	    x_epsi[1] = x[1] + epsi;
	    computeResidual(x_epsi, res_epsi);
	    gradient[1] = (res_epsi[0] - res[0])/epsi;
	} else if (gradient_method == 2)  {
	    double s = x[0];
	    double c = x[1];
	    double ff_ds_dc[2];
	    double ff = tm.fracFlowWithDer(s, c, cell, ff_ds_dc);
	    double mc_dc;
	    double mc = tm.computeMcWithDer(c, &mc_dc);
	    double ads_dc;
            double ads = tm.polyprops_.adsorbtionWithDer(c, cmax0, &ads_dc);
	    res[0] = s - s0 +  dtpv*(outflux*ff + influx);
	    res[1] = (1 - dps)*s*c - (1 - dps)*s0*c0
		+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
		+ dtpv*(outflux*ff*mc + influx_polymer);
	    gradient[0] = 1 + dtpv*outflux*ff_ds_dc[0];
	    gradient[1] = dtpv*outflux*ff_ds_dc[1];
	}
    }

    void TransportModelPolymer::ResidualEquation::computeGradientResC(const double* x, double* res, double* gradient) const
    // If gradient_method == 1, use finite difference
    // If gradient_method == 2, use analytic expresions
    {
	if (gradient_method == 1) {
	    double epsi = 1e-8;
	    double res_epsi[2];
	    double x_epsi[2];
	    computeResidual(x, res);
	    x_epsi[0] = x[0] + epsi;
	    x_epsi[1] = x[1];
	    computeResidual(x_epsi, res_epsi);
	    gradient[0] =  (res_epsi[1] - res[1])/epsi;
	    x_epsi[0] = x[0];
	    x_epsi[1] = x[1] + epsi;
	    computeResidual(x_epsi, res_epsi);
	    gradient[1] = (res_epsi[1] - res[1])/epsi;
	} else if (gradient_method == 2)  {
	    double s = x[0];
	    double c = x[1];
	    double ff_ds_dc[2];
	    double ff = tm.fracFlowWithDer(s, c, cell, ff_ds_dc);
	    double mc_dc;
	    double mc = tm.computeMcWithDer(c, &mc_dc);
	    double dps = tm.polyprops_.deadPoreVol();
	    double rhor = tm.polyprops_.rockDensity();
	    double ads0 = tm.polyprops_.adsorbtion(c0, cmax0);
	    double ads_dc;
            double ads = tm.polyprops_.adsorbtionWithDer(c, cmax0, &ads_dc);
	    res[0] = s - s0 +  dtpv*(outflux*ff + influx);
	    res[1] = (1 - dps)*s*c - (1 - dps)*s0*c0
		+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
		+ dtpv*(outflux*ff*mc + influx_polymer);
	    gradient[0] = (1 - dps)*c + dtpv*outflux*(ff_ds_dc[0])*mc;
	    gradient[1] = (1 - dps)*s + rhor*((1.0 - porosity)/porosity)*ads_dc
		+ dtpv*outflux*(ff_ds_dc[1]*mc + ff*mc_dc);
	}
    }

    // Compute the Jacobian of the residual equations.
    void TransportModelPolymer::ResidualEquation::computeJacobiRes(const double* x, double* res_s_ds_dc, double* res_c_ds_dc) const
    {
	if (gradient_method == 1) {
	    double epsi = 1e-8;
	    double res[2];
	    double res_epsi[2];
	    double x_epsi[2];
	    computeResidual(x, res);
	    x_epsi[0] = x[0] + epsi;
	    x_epsi[1] = x[1];
	    computeResidual(x_epsi, res_epsi);
	    res_s_ds_dc[0] =  (res_epsi[0] - res[0])/epsi;
	    x_epsi[0] = x[0];
	    x_epsi[1] = x[1] + epsi;
	    computeResidual(x_epsi, res_epsi);
	    res_s_ds_dc[1] = (res_epsi[0] - res[0])/epsi;
	    x_epsi[0] = x[0] + epsi;
	    x_epsi[1] = x[1];
	    computeResidual(x_epsi, res_epsi);
	    res_c_ds_dc[0] =  (res_epsi[1] - res[1])/epsi;
	    x_epsi[0] = x[0];
	    x_epsi[1] = x[1] + epsi;
	    computeResidual(x_epsi, res_epsi);
	    res_c_ds_dc[1] = (res_epsi[1] - res[1])/epsi;
	} else if (gradient_method == 2)  {
	    double s = x[0];
	    double c = x[1];
	    double ff_ds_dc[2];
	    double ff = tm.fracFlowWithDer(s, c, cell, ff_ds_dc);
	    double mc_dc;
	    double mc = tm.computeMcWithDer(c, &mc_dc);
	    double ads_dc;
            tm.polyprops_.adsorbtionWithDer(c, cmax0, &ads_dc);
	    res_s_ds_dc[0] = 1 + dtpv*outflux*ff_ds_dc[0];
	    res_s_ds_dc[1] = dtpv*outflux*ff_ds_dc[1];
	    res_c_ds_dc[0] = (1 - dps)*c + dtpv*outflux*(ff_ds_dc[0])*mc;
	    res_c_ds_dc[1] = (1 - dps)*s + rhor*((1.0 - porosity)/porosity)*ads_dc
		+ dtpv*outflux*(ff_ds_dc[1]*mc + ff*mc_dc);
	}
    }


    void TransportModelPolymer::solveSingleCell(const int cell)
    {
	switch (method_) {
	case Bracketing:
	    solveSingleCellBracketing(cell);
	    break;
	case Newton:
	    solveSingleCellNewton(cell);
	    break;
	default:
	    THROW("Unknown method " << method_);
	}
    }


    void TransportModelPolymer::solveSingleCellBracketing(int cell)
    {
	ResidualC res(*this, cell);
	const double a = 0.0;
	const double b = polyprops_.cMax();
	int iters_used;

	// Check if current state is an acceptable solution.
	double res_sc[2];
	double mc, ff;
	res.computeBothResiduals(saturation_[cell], concentration_[cell], res_sc[0], res_sc[1], mc, ff);
	if (norm(res_sc) < tol_) {
	    fractionalflow_[cell] = ff;
	    mc_[cell] = mc;
	    return;
	}

	concentration_[cell] = modifiedRegulaFalsi(res, a, b, maxit_, tol_, iters_used);
	cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
	saturation_[cell] = res.lastSaturation();
	fractionalflow_[cell] = fracFlow(saturation_[cell], concentration_[cell], cell);
	mc_[cell] = computeMc(concentration_[cell]);
    }



    // Newton method, where we first try a Newton step. Then, if it does not work well, we look for
    // the zero of either the residual in s or the residual in c along a specified piecewise linear
    // curve. In these cases, we can use a robust 1d solver.
    void TransportModelPolymer::solveSingleCellNewton(int cell)
    {
	// the tolerance for 1d solver is set as a function of the residual, because if we are far
	// from the solution we do not need a very accurate 1d solver (recall that the 1d solver
	// solves for only one of the two residuals)
	// The tolerance falsi_tol is improved by
	// (reduc_factor_falsi_tol * "previous residual") at each step
	double falsi_tol;
	const double reduc_factor_falsi_tol = 1e-2;
	int iters_used_falsi = 0;
	const int max_iters_split = maxit_;
	int iters_used_split = 0;

	// Check if current state is an acceptable solution.
	ResidualEquation res_eq(*this, cell);
	double x[2] = {saturation_[cell], concentration_[cell]};
	double res[2];
	double mc;
	double ff;
	res_eq.computeResidual(x, res, mc, ff);
	if (norm(res) <= tol_) {
	    cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
 	    fractionalflow_[cell] = ff;
	    mc_[cell] = mc;
	    return;
	}

	falsi_tol =  std::max(reduc_factor_falsi_tol*norm(res), tol_);
	double x_min[2] = { std::max(polyprops_.deadPoreVol(), smin_[2*cell]), 0.0 };
	double x_max[2] = { 1.0, polyprops_.cMax() };
	double t;
	double t_max;
	double t_out;
	double direction[2];
	double end_point[2];
	double gradient[2];
	bool unsuccessfull_newton_step = true;
	double x_new[2];
	double res_new[2];
	ResSOnCurve res_s_on_curve(res_eq);
	ResCOnCurve res_c_on_curve(res_eq);
	bool if_res_s;
	int counter_drop_newton = 0;
	bool not_so_successfull_newton_step = false;


 	while ((norm(res) > tol_) && (iters_used_split < max_iters_split)) {
	    // We first try a Newton step
	    if (counter_drop_newton == 0 && solveNewtonStep(x, res_eq, res, x_new)) {
		res_eq.computeResidual(x_new, res_new, mc, ff);
		unsuccessfull_newton_step = false;
		not_so_successfull_newton_step = false;
		if (norm(res_new) > norm(res) || x_new[0] < x_min[0] || x_new[1] < x_min[1] || x_new[0] > x_max[0] || x_new[1] > x_max[1]) {
		    unsuccessfull_newton_step = true;
		} else  {
		    x[0] = x_new[0];
		    x[1] = x_new[1];
		    if (norm(res_new) > 1e-1*norm(res) && norm(res_new) < 1e1*tol_) {
			// We are close to the solution and Newton does not perform well.
			// Then, we drop Newton for a given number of iterations.
			not_so_successfull_newton_step = true;
			counter_drop_newton = 4;
		    }
		    res[0] = res_new[0];
		    res[1] = res_new[1];
		    iters_used_split += 1;
		}
	    } else {
		unsuccessfull_newton_step = true;
	    }

	    if (not_so_successfull_newton_step || unsuccessfull_newton_step) {
		// Newton was not satisfactory. We start 1d solvers.
		if (not_so_successfull_newton_step) {
		    counter_drop_newton -= 1;
		}
		// General comment on the zero curves of the s and c residuals:
		// Typically res_s(s,c)=0 defines an increasing curve in the s-c plane while
		// res_c(s,c)=0 defines a decreasing curve. However, we do not assume that in the algorithm.
		// We know that res_s(x_top_left)<0, res_s(x_bottom_right)>0
		// and res_c(x_bottom_left)<0, res_c(x_top_right)>0
		// Here, "top", "bottom", ... refer to the corner of the admissible box of (s,c) values.
		// We use these results to construct a 1d curve for which we are sure that res_s or res_c change sign
		// and which can therefore be used by a 1d solver.

		// We start with the zero curve of the s and r residual we are closest to.
		if (std::abs(res[0]) < std::abs(res[1])) {
		    falsi_tol =  std::max(reduc_factor_falsi_tol*std::abs(res[0]), tol_);
		    if (res[0] < -falsi_tol) {
			direction[0] = x_max[0] - x[0];
			direction[1] = x_min[1] - x[1];
			if_res_s = true;
		    } else if (res[0] > falsi_tol) {
			direction[0] = x_min[0] - x[0];
			direction[1] = x_max[1] - x[1];
			if_res_s = true;
		    } else {
			res_eq.computeGradientResS(x, res, gradient);
			direction[0] = -gradient[1];
			direction[1] = gradient[0];
			if_res_s = false;
		    }
		} else {
		    falsi_tol =  std::max(reduc_factor_falsi_tol*std::abs(res[1]), tol_);
		    if (res[1] < -falsi_tol) {
			direction[0] = x_max[0] - x[0];
			direction[1] = x_max[1] - x[1];
			if_res_s = false;
		    } else if (res[1] > falsi_tol) {
			direction[0] = x_min[0] - x[0];
			direction[1] = x_min[1] - x[1];
			if_res_s = false;
		    } else {
			res_eq.computeGradientResC(x, res, gradient);
			direction[0] = -gradient[1];
			direction[1] = gradient[0];
			if_res_s = true;
		    }
		}
		if (if_res_s) {
		    if (res[0] < 0) {
			end_point[0] = x_max[0];
			end_point[1] = x_min[1];
			res_s_on_curve.curve.setup(x, direction, end_point, x_min, x_max, t_max, t_out);
			if (res_s_on_curve(t_out) >= 0) {
			    t_max = t_out;
			}
		    } else {
			end_point[0] = x_min[0];
			end_point[1] = x_max[1];
			res_s_on_curve.curve.setup(x, direction, end_point, x_min, x_max, t_max, t_out);
			if (res_s_on_curve(t_out) <= 0) {
			    t_max = t_out;
			}
		    }
		    // Note: In some experiments modifiedRegularFalsi does not yield a result under the given tolerance.
		    t = modifiedRegulaFalsi(res_s_on_curve, 0., t_max, maxit_, falsi_tol, iters_used_falsi);
		    res_s_on_curve.curve.computeXOfT(x, t);
		} else {
		    if (res[1] < 0) {
			end_point[0] = x_max[0];
			end_point[1] = x_max[1];
			res_c_on_curve.curve.setup(x, direction, end_point, x_min, x_max, t_max, t_out);
			if (res_c_on_curve(t_out) >= 0) {
			    t_max = t_out;
			}
		    } else {
			end_point[0] = x_min[0];
			end_point[1] = x_min[1];
			res_c_on_curve.curve.setup(x, direction, end_point, x_min, x_max, t_max, t_out);
			if (res_c_on_curve(t_out) <= 0) {
			    t_max = t_out;
			}
		    }
		    t = modifiedRegulaFalsi(res_c_on_curve, 0., t_max, maxit_, falsi_tol, iters_used_falsi);
		    res_c_on_curve.curve.computeXOfT(x, t);

		}

		res_eq.computeResidual(x, res, mc, ff);
		iters_used_split += 1;
	    }
	}


	if ((iters_used_split >=  max_iters_split) && (norm(res) > tol_)) {
	    MESSAGE("Newton for single cell did not work in cell number " << cell);
	    solveSingleCellBracketing(cell);
	} else {
	    concentration_[cell] = x[1];
	    cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
	    saturation_[cell] = x[0];
	    fractionalflow_[cell] = ff;
	    mc_[cell] = mc;
	}
    }

    void TransportModelPolymer::solveMultiCell(const int num_cells, const int* cells)
    {
	double max_s_change = 0.0;
	double max_c_change = 0.0;
	int num_iters = 0;
	// Must store state variables before we start.
	std::vector<double> s0(num_cells);
	std::vector<double> c0(num_cells);
	std::vector<double> cmax0(num_cells);
	// Must set initial fractional flows etc. before we start.
	for (int i = 0; i < num_cells; ++i) {
	    const int cell = cells[i];
	    fractionalflow_[cell] = fracFlow(saturation_[cell], concentration_[cell], cell);
	    mc_[cell] = computeMc(concentration_[cell]);
	    s0[i] = saturation_[cell];
	    c0[i] = concentration_[cell];
	    cmax0[i] = cmax_[i];
	}
	do {
	    int max_s_change_cell = -1;
	    int max_c_change_cell = -1;
	    max_s_change = 0.0;
	    max_c_change = 0.0;
	    for (int i = 0; i < num_cells; ++i) {
		const int cell = cells[i];
		const double old_s = saturation_[cell];
		const double old_c = concentration_[cell];
		saturation_[cell] = s0[i];
		concentration_[cell] = c0[i];
		cmax_[cell] = cmax0[i];
		solveSingleCell(cell);
		// std::cout << "cell = " << cell << "    delta s = " << saturation_[cell] - old_s << std::endl;
		if (max_s_change < std::fabs(saturation_[cell] - old_s)) {
		    max_s_change_cell = cell;
		}
		if (max_c_change < std::fabs(concentration_[cell] - old_c)) {
		    max_c_change_cell = cell;
		}
		max_s_change = std::max(max_s_change, std::fabs(saturation_[cell] - old_s));
		max_c_change = std::max(max_c_change, std::fabs(concentration_[cell] - old_c));
	    }
	    // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change
	    // 	      << "    in cell " << max_change_cell << std::endl;
	} while (((max_s_change > tol_) || (max_c_change > tol_)) && ++num_iters < maxit_);
	if (max_s_change > tol_) {
	    THROW("In solveMultiCell(), we did not converge after "
		  << num_iters << " iterations. Delta s = " << max_s_change);
	}
	if (max_c_change > tol_) {
	    THROW("In solveMultiCell(), we did not converge after "
		  << num_iters << " iterations. Delta c = " << max_c_change);
	}
	std::cout << "Solved " << num_cells << " cell multicell problem in "
		  << num_iters << " iterations." << std::endl;
    }




    double TransportModelPolymer::fracFlow(double s, double c, int cell) const
    {
	double c_max_limit = polyprops_.cMax();
	double cbar = c/c_max_limit;
	double c_ads = polyprops_.adsorbtion(c, cmax_[cell]);
	double c_max_ads = polyprops_.cMaxAds();
        double res_factor = polyprops_.resFactor();
        double res_k = 1 + (res_factor - 1)*c_ads/c_max_ads;
	double mu_w = visc_[0];
	double mu_m = polyprops_.viscMult(c)*mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.cMax())*mu_w;
	double omega = polyprops_.mixParam();
	double mu_m_omega = std::pow(mu_m, omega);
	double mu_w_e   = mu_m_omega*std::pow(mu_w, 1.0 - omega);
	double mu_p_eff = mu_m_omega*std::pow(mu_p, 1.0 - omega);
	double inv_mu_w_eff = (1.0 - cbar)/mu_w_e + cbar/mu_p_eff;
	double inv_visc_eff[2] = { inv_mu_w_eff, 1.0/visc_[1] };
	double sat[2] = { s, 1.0 - s };
	double mob[2];
	props_.relperm(1, sat, &cell, mob, 0);
	mob[0] *= inv_visc_eff[0]/res_k;
	mob[1] *= inv_visc_eff[1];
	return mob[0]/(mob[0] + mob[1]);
    }

    double TransportModelPolymer::fracFlowWithDer(double s, double c, int cell, double* der) const
    {
	double c_max_limit = polyprops_.cMax();
	double cbar = c/c_max_limit;
        double c_ads_dc;
        double c_max = cmax_[cell];
        double c_ads = polyprops_.adsorbtionWithDer(c, c_max, &c_ads_dc);
	double c_max_ads = polyprops_.cMaxAds();
        double res_factor = polyprops_.resFactor();
        double res_k = 1 + (res_factor - 1)*c_ads/c_max_ads;
        double res_k_dc = (res_factor - 1)*c_ads_dc/c_max_ads;
	double mu_w = visc_[0];
	double mu_m_dc; // derivative of mu_m with respect to c
	double mu_m = polyprops_.viscMultWithDer(c, &mu_m_dc)*mu_w;
	mu_m_dc *= mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.cMax())*mu_w;
	double omega = polyprops_.mixParam();
 	double mu_w_e   = std::pow(mu_m, omega)*std::pow(mu_w, 1 - omega);
	double mu_w_e_dc = omega*mu_m_dc*std::pow(mu_m, omega - 1)*std::pow(mu_w, 1 - omega);
	double mu_p_eff = std::pow(mu_m, omega)*std::pow(mu_p, 1 - omega);
	double mu_p_eff_dc = omega*mu_m_dc*std::pow(mu_m, omega - 1)*std::pow(mu_p, 1 - omega);
	double mu_w_eff = 1./((1 - cbar)/mu_w_e + cbar/mu_p_eff);
	double mu_w_eff_dc = -1./c_max_limit*mu_w_eff*mu_w_eff*(1./mu_p_eff - 1./mu_w_e)
	    + (1-cbar)*(mu_w_eff*mu_w_eff/(mu_w_e*mu_w_e))*mu_w_e_dc
	    + cbar*(mu_w_eff*mu_w_eff/(mu_p_eff*mu_p_eff))*mu_p_eff_dc;
	double visc_eff[2] = { mu_w_eff, visc_[1] };
	double sat[2] = { s, 1.0 - s };
	double mob[2];
	double mob_ds[2];
	double mob_dc[2];
	double perm[2];
	double perm_ds[4];
	props_.relperm(1, sat, &cell, perm, perm_ds);
	mob[0] = perm[0]/visc_eff[0]/res_k;
	mob[1] = perm[1]/visc_eff[1];
	mob_ds[0] = perm_ds[0]/visc_eff[0]/res_k;
	mob_ds[1] = perm_ds[1]/visc_eff[1];
	mob_dc[0] = - perm[0]*(mu_w_eff_dc/(mu_w_eff*mu_w_eff*res_k) + res_k_dc/(res_k*res_k*mu_w_eff));
	mob_dc[1] = 0.;
	der[0] = (mob_ds[0]*mob[1] - mob_ds[1]*mob[0])/((mob[0] + mob[1])*(mob[0] + mob[1]));
	der[1] = (mob_dc[0]*mob[1] - mob_dc[1]*mob[0])/((mob[0] + mob[1])*(mob[0] + mob[1]));
 	return mob[0]/(mob[0] + mob[1]);
    }

    double TransportModelPolymer::computeMc(double c) const
    {
	double c_max_limit = polyprops_.cMax();
	double cbar = c/c_max_limit;
	double mu_w = visc_[0];
	double mu_m = polyprops_.viscMult(c)*mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.cMax())*mu_w;
	double omega = polyprops_.mixParam();
	double mu_m_omega = std::pow(mu_m, omega);
	double mu_w_e   = mu_m_omega*std::pow(mu_w, 1.0 - omega);
	double mu_p_eff = mu_m_omega*std::pow(mu_p, 1.0 - omega);
	double inv_mu_w_eff = (1.0 - cbar)/mu_w_e + cbar/mu_p_eff;
	return c/(inv_mu_w_eff*mu_p_eff);
    }

    double TransportModelPolymer::computeMcWithDer(double c, double* der) const
    {
	double c_max_limit = polyprops_.cMax();
	double cbar = c/c_max_limit;
	double mu_w = visc_[0];
	double mu_m_dc; // derivative of mu_m with respect to c
	double mu_m = polyprops_.viscMultWithDer(c, &mu_m_dc)*mu_w;
	mu_m_dc *= mu_w;
	double mu_p = polyprops_.viscMult(polyprops_.cMax())*mu_w;
	double omega = polyprops_.mixParam();
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


namespace
{

    CurveInSCPlane::CurveInSCPlane()
	{
	}

    // Setup the curve (see comment above).
    // The curve is parametrized by t in [0, t_max], t_out is equal to t when the curve hits the bounding
    // rectangle. x_out=(s_out, c_out) denotes the values of s and c at that point.
    void CurveInSCPlane::setup(const double* x, const double* direction,
		   const double* end_point, const double* x_min,
		   const double* x_max, double& t_max_out,
		   double& t_out_out)
    {
	x_[0] = x[0];
	x_[1] = x[1];
	x_max_[0] = x_max[0];
	x_max_[1] = x_max[1];
	x_min_[0] = x_min[0];
	x_min_[1] = x_min[1];
	direction_[0] = direction[0];
	direction_[1] = direction[1];
	end_point_[0] = end_point[0];
	end_point_[1] = end_point[1];
	if ((end_point_[0]-x_[0])*direction_[0] + (end_point_[1]-x_[1])*direction_[1] < 0) {
	    direction_[0] *= -1.0;
	    direction_[1] *= -1.0;
	}
	if ((std::abs(direction_[0]) + std::abs(direction_[0])) == 0) {
	    direction_[0] = end_point_[0]-x_[0];
	    direction_[1] = end_point_[1]-x_[1];
	}
	bool t0_exists = true;
	double t0 = 0; // dummy default value (so that compiler does not complain).
	if (direction_[0] > 0) {
	    t0 = (x_max_[0] - x_[0])/direction_[0];
	} else if (direction_[0] < 0) {
	    t0 = (x_min_[0] - x_[0])/direction_[0];
	} else {
	    t0_exists = false;
	}
	bool t1_exists = true;
	double t1 = 0; // dummy default value.
	if (direction_[1] > 0) {
	    t1 = (x_max_[1] - x_[1])/direction_[1];
	} else if (direction[1] < 0) {
	    t1 = (x_min_[1] - x_[1])/direction_[1];
	} else {
	    t1_exists = false;
	}
	if (t0_exists) {
	    if (t1_exists) {
		t_out_ = std::min(t0, t1);
	    } else {
		t_out_ = t0;
	    }
	} else if (t1_exists) {
	    t_out_ = t1;
	} else {
	    THROW("Direction illegal: is a zero vector.");
	}
	x_out_[0] = x_[0] + t_out_*direction_[0];
	x_out_[1] = x_[1] + t_out_*direction_[1];
	t_max_ = t_out_ + 1;
	t_max_out = t_max_;
	t_out_out = t_out_;
    }


    // Compute x=(s,c) for a given t (t is the parameter for the piecewise linear curve)
    void CurveInSCPlane::computeXOfT(double* x_of_t, const double t) const {
	if (t <= t_out_) {
	    x_of_t[0] = x_[0] + t*direction_[0];
		x_of_t[1] = x_[1] + t*direction_[1];
	} else {
	    x_of_t[0] = 1/(t_max_-t_out_)*((t_max_ - t)*x_out_[0] + end_point_[0]*(t - t_out_));
	    x_of_t[1] = 1/(t_max_-t_out_)*((t_max_ - t)*x_out_[1] + end_point_[1]*(t - t_out_));
	}
    }


    ResSOnCurve::ResSOnCurve(const Opm::TransportModelPolymer::ResidualEquation& res_eq)
	: res_eq_(res_eq)
    {
    }

    double ResSOnCurve::operator()(const double t) const
    {
	double x_of_t[2];
	curve.computeXOfT(x_of_t, t);
	return res_eq_.computeResidualS(x_of_t);
    }

    ResCOnCurve::ResCOnCurve(const Opm::TransportModelPolymer::ResidualEquation& res_eq)
	: res_eq_(res_eq)
    {
    }

    double ResCOnCurve::operator()(const double t) const
    {
	double x_of_t[2];
	curve.computeXOfT(x_of_t, t);
	return res_eq_.computeResidualC(x_of_t);
    }

    bool solveNewtonStep(const double* x, const  Opm::TransportModelPolymer::ResidualEquation& res_eq,
			 const double* res, double* x_new) {

    	double res_s_ds_dc[2];
    	double res_c_ds_dc[2];

	res_eq.computeJacobiRes(x, res_s_ds_dc, res_c_ds_dc);

    	double det = res_s_ds_dc[0]*res_c_ds_dc[1] - res_c_ds_dc[0]*res_s_ds_dc[1];
    	if (std::abs(det) < 1e-8) {
    	    return false;
    	} else {
    	    x_new[0] = x[0] - (res[0]*res_c_ds_dc[1] - res[1]*res_s_ds_dc[1])/det;
    	    x_new[1] = x[1] - (res[1]*res_s_ds_dc[0] - res[0]*res_c_ds_dc[0])/det;
	    return true;
	}
    }

} // Anonymous namespace


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
