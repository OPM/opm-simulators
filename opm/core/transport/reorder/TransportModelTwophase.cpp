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

#include <opm/core/transport/reorder/TransportModelTwophase.hpp>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/utility/RootFinders.hpp>

#include <fstream>
#include <iterator>

namespace Opm
{


    TransportModelTwophase::TransportModelTwophase(const UnstructuredGrid& grid,
						   const double* porevolume,
						   const Opm::IncompPropertiesInterface& props)
	: grid_(grid),
	  porevolume_(porevolume),
	  props_(props),
	  darcyflux_(0),
	  source_(0),
	  dt_(0.0),
	  saturation_(0),
	  fractionalflow_(grid.number_of_cells, -1.0)
	  // ia_(grid.number_of_cells + 1, -1),
	  // ja_(grid.number_of_faces, -1)
    {
	if (props.numPhases() != 2) {
	    THROW("Property object must have 2 phases");
	}
	visc_ = props.viscosity();
    }

    void TransportModelTwophase::solve(const double* darcyflux,
				       const double* source,
				       const double dt,
				       double* saturation)
    {
	darcyflux_ = darcyflux;
	source_ = source;
	dt_ = dt;
	saturation_ = saturation;

	// std::vector<int> seq(grid_.number_of_cells);
	// std::vector<int> comp(grid_.number_of_cells + 1);
	// int ncomp;
	// compute_sequence_graph(&grid_, darcyflux_,
	// 		       &seq[0], &comp[0], &ncomp,
	// 		       &ia_[0], &ja_[0]);


	reorderAndTransport(grid_, darcyflux);
    }

    // Residual function r(s) for a single-cell implicit Euler transport
    //
    //     r(s) = s - s0 + dt/pv*( influx + outflux*f(s) )
    //
    // where influx is water influx, outflux is total outflux.
    // Influxes are negative, outfluxes positive.
    struct TransportModelTwophase::Residual
    {
	int cell;
	double s0;
	double influx;  // sum_j min(v_ij, 0)*f(s_j)
	double outflux; // sum_j max(v_ij, 0)
	double dtpv;    // dt/pv(i)
	const TransportModelTwophase& tm;
	explicit Residual(const TransportModelTwophase& tmodel, int cell_index)
	    : tm(tmodel)
	{
	    cell    = cell_index;
	    s0      = tm.saturation_[cell];
            double dflux       = -tm.source_[cell];
            bool src_is_inflow = dflux < 0.0;
	    influx  =  src_is_inflow ? dflux : 0.0;
	    outflux = !src_is_inflow ? dflux : 0.0;
	    dtpv    = tm.dt_/tm.porevolume_[cell];

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
		    } else {
			outflux += flux;
		    }
		}
	    }
	}
	double operator()(double s) const
	{
	    return s - s0 +  dtpv*(outflux*tm.fracFlow(s, cell) + influx);
	}
    };


    void TransportModelTwophase::solveSingleCell(const int cell)
    {
	Residual res(*this, cell);
	const double tol = 1e-9;
	// const double r0 = res(saturation_[cell]);
	// if (std::fabs(r0) < tol) {
	//     return;
	// }
	const double a = 0.0;
	const double b = 1.0;
	const int maxit = 20;
	int iters_used;
	// saturation_[cell] = modifiedRegulaFalsi(res, a, b, saturation_[cell], maxit, tol, iters_used);
	saturation_[cell] = modifiedRegulaFalsi(res, a, b, maxit, tol, iters_used);
	fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
    }

    void TransportModelTwophase::solveMultiCell(const int num_cells, const int* cells)
    {
	// std::ofstream os("dump");
	// std::copy(cells, cells + num_cells, std::ostream_iterator<double>(os, "\n"));

	// Experiment: try a breath-first search to build a more suitable ordering.
	// Verdict: failed to improve #iterations.
	// {
	//     std::vector<int> pos(grid_.number_of_cells, -1);
	//     for (int i = 0; i < num_cells; ++i) {
	// 	const int cell = cells[i];
	// 	pos[cell] = i;
	//     }
	//     std::vector<int> done_pos(num_cells, 0);
	//     std::vector<int> upstream_pos;
	//     std::vector<int> new_pos;
	//     upstream_pos.push_back(0);
	//     done_pos[0] = 1;
	//     int current = 0;
	//     while (new_pos.size() < num_cells) {
	// 	const int i = upstream_pos[current++];
	// 	new_pos.push_back(i);
	// 	const int cell = cells[i];
	// 	for (int j = ia_[cell]; j < ia_[cell+1]; ++j) {
	// 	    const int opos = pos[ja_[j]];
	// 	    if (!done_pos[opos]) {
	// 		upstream_pos.push_back(opos);
	// 		done_pos[opos] = 1;
	// 	    }
	// 	}
	//     }
	//     std::reverse(new_pos.begin(), new_pos.end());
	//     std::copy(new_pos.begin(), new_pos.end(), const_cast<int*>(cells));
	// }

	// Experiment: try a random ordering.
	// Verdict: amazingly, reduced #iterations by approx. 25%!
	// int* c = const_cast<int*>(cells);
	// std::random_shuffle(c, c + num_cells);

	// Experiment: implement a metric measuring badness of ordering
	//             as average distance in (cyclic) ordering from
	//             upstream neighbours.
	// Verdict: does not seem to predict #iterations very well, if at all.
	// std::vector<int> pos(grid_.number_of_cells, -1);
	// for (int i = 0; i < num_cells; ++i) {
	//     const int cell = cells[i];
	//     pos[cell] = i;
	// }
	// double diffsum = 0;
	// for (int i = 0; i < num_cells; ++i) {
	//     const int cell = cells[i];
	//     int num_upstream = 0;
	//     int loc_diffsum = 0;
	//     for (int j = ia_[cell]; j < ia_[cell+1]; ++j) {
	// 	const int opos = pos[ja_[j]];
	// 	if (opos == -1) {
	// 	    std::cout << "Hmmm" << std::endl;
	// 	    continue;
	// 	}
	// 	++num_upstream;
	// 	const int diff = (i - opos + num_cells) % num_cells;
	// 	loc_diffsum += diff;
	//     }
	//     diffsum += double(loc_diffsum)/double(num_upstream);
	// }
	// std::cout << "Average distance from upstream neighbours: " << diffsum/double(num_cells)
	// 	  << std::endl;

	double max_s_change = 0.0;
	const double tol = 1e-9;
	const int max_iters = 300;
	int num_iters = 0;
	// Must store s0 before we start.
	std::vector<double> s0(num_cells);
	// Must set initial fractional flows before we start.
	for (int i = 0; i < num_cells; ++i) {
	    const int cell = cells[i];
	    fractionalflow_[cell] = fracFlow(saturation_[cell], cell);
	    s0[i] = saturation_[cell];
	}
	do {
	    int max_change_cell = -1;
	    max_s_change = 0.0;
	    for (int i = 0; i < num_cells; ++i) {
		const int cell = cells[i];
		const double old_s = saturation_[cell];
		saturation_[cell] = s0[i];
		solveSingleCell(cell);
		// std::cout << "cell = " << cell << "    delta s = " << saturation_[cell] - old_s << std::endl;
		if (max_s_change < std::fabs(saturation_[cell] - old_s)) {
		    max_change_cell = cell;
		}
		max_s_change = std::max(max_s_change, std::fabs(saturation_[cell] - old_s));
	    }
	    // std::cout << "Iter = " << num_iters << "    max_s_change = " << max_s_change
	    // 	      << "    in cell " << max_change_cell << std::endl;
	} while (max_s_change > tol && ++num_iters < max_iters);
	if (max_s_change > tol) {
	    THROW("In solveMultiCell(), we did not converge after "
	    	  << num_iters << " iterations. Delta s = " << max_s_change);
	}
	std::cout << "Solved " << num_cells << " cell multicell problem in "
		  << num_iters << " iterations." << std::endl;
    }

    double TransportModelTwophase::fracFlow(double s, int cell) const
    {
	double sat[2] = { s, 1.0 - s };
	double mob[2];
	props_.relperm(1, sat, &cell, mob, 0);
	mob[0] /= visc_[0];
	mob[1] /= visc_[1];
	return mob[0]/(mob[0] + mob[1]);
    }




} // namespace Opm



/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
