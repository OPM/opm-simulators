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

namespace Opm
{


    TransportModelPolymer::TransportModelPolymer(const UnstructuredGrid& grid,
						 const double* porosity,
						 const double* porevolume,
						 const IncompPropertiesInterface& props,
						 const PolymerData& polyprops)
	: grid_(grid),
	  porosity_(porosity),
	  porevolume_(porevolume),
	  props_(props),
	  polyprops_(polyprops),
	  darcyflux_(0),
	  source_(0),
	  dt_(0.0),
	  inflow_c_(0.0),
	  saturation_(0),
	  concentration_(0),
	  cmax_(0),
	  fractionalflow_(grid.number_of_cells, -1.0),
	  mc_(grid.number_of_cells, -1.0)
    {
	if (props.numPhases() != 2) {
	    THROW("Property object must have 2 phases");
	}
	visc_ = props.viscosity();
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
	    const double a = 0.2;   // TODO: Make this a proper s_min value.
	    const double b = 1.0;
	    const int maxit = 20;
	    const double tol = 1e-9;
	    int iters_used;
	    // Solve for s first.
	    s = modifiedRegulaFalsi(res_s, a, b, maxit, tol, iters_used);
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



    void TransportModelPolymer::solveSingleCell(const int cell)
    {
	ResidualC res(*this, cell);
	const double a = 0.0;
	const double b = polyprops_.c_max_limit;
	const int maxit = 20;
	const double tol = 1e-9;
	int iters_used;
	concentration_[cell] = modifiedRegulaFalsi(res, a, b, maxit, tol, iters_used);
	cmax_[cell] = std::max(cmax_[cell], concentration_[cell]);
	saturation_[cell] = res.lastSaturation();
	fractionalflow_[cell] = fracFlow(saturation_[cell], concentration_[cell], cell);
	mc_[cell] = computeMc(concentration_[cell]);
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




} // namespace Opm



/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
