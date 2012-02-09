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
#include <opm/core/transport/reorder/nlsolvers.h>

struct Opm::Parameters
{
    double s0;
    double influx;  /* sum_j min(v_ij, 0)*f(s_j) */
    double outflux; /* sum_j max(v_ij, 0)        */
    double dtpv;    /* dt/pv(i)                  */
    int cell;
    const Opm::IncompPropertiesInterface* props;
};

static double residual(double s, void *data);
static double fluxfun_props(double s, int cell, const Opm::IncompPropertiesInterface* props);

namespace Opm
{


    TransportModelTwophase::TransportModelTwophase(const UnstructuredGrid* grid,
						   const Opm::IncompPropertiesInterface* props,
						   const double* darcyflux,
						   const double* porevolume,
						   const double* source,
						   const double dt,
						   double* saturation)
	: grid_(grid),
	  props_(props),
	  darcyflux_(darcyflux),
	  porevolume_(porevolume),
	  source_(source),
	  dt_(dt),
	  saturation_(saturation),
	  fractionalflow_(grid->number_of_cells, 0.0)
    {
    }


    void TransportModelTwophase::solveSingleCell(int cell)
    {
	NonlinearSolverCtrl ctrl;
	ctrl.maxiterations = 20;
	ctrl.nltolerance   = 1e-9;
	ctrl.method        = NonlinearSolverCtrl::REGULAFALSI;
	ctrl.initialguess  = 0.5;
	ctrl.min_bracket   = 0.0;
	ctrl.max_bracket   = 1.0;

	Parameters prm = getParameters(cell);
	saturation_[cell] = find_zero(residual, &prm, &ctrl);
	// double ff1 = fluxfun_props(d->saturation[cell], cell, d->props);
	// double ff2 = fluxfun(d->saturation[cell], -999);
	// printf("New = %f   old = %f\n", ff1, ff2);
	fractionalflow_[cell] = fluxfun_props(saturation_[cell], cell, props_);
    }


    Opm::Parameters TransportModelTwophase::getParameters(int cell)
    {
	int i;
	Opm::Parameters p;
	double flux;
	int f, other;

	p.s0      = saturation_[cell];
	p.influx  = source_[cell] >  0 ? -source_[cell] : 0.0;
	p.outflux = source_[cell] <= 0 ? -source_[cell] : 0.0;
	p.dtpv    = dt_/porevolume_[cell];
	p.cell    = cell;
	p.props   = props_;

	for (i=grid_->cell_facepos[cell]; i<grid_->cell_facepos[cell+1]; ++i) {
	    f = grid_->cell_faces[i];

	    /* Compute cell flux*/
	    if (cell == grid_->face_cells[2*f]) {
		flux  = darcyflux_[f];
		other = grid_->face_cells[2*f+1];
	    }
	    else {
		flux  =-darcyflux_[f];
		other = grid_->face_cells[2*f];
	    }

	    if (other != -1) {
		if (flux < 0.0) {
		    p.influx  += flux*fractionalflow_[other];
		}
		else {
		    p.outflux += flux;
		}
	    }
	}
	return p;
    }


} // namespace Opm

/* ====================== Internals =================================*/


/* Residual function r(s) for a single-cell bvp */
/*
 *     r(s) = s - s0 + dt/pv*(influx - outflux*f(s) )
 */
/* influx is water influx, outflux is total outflux */
static double
residual(double s, void *data)
{
    Opm::Parameters *p = (Opm::Parameters*) data;
    return s - p->s0 +  p->dtpv*(p->outflux*fluxfun_props(s, p->cell, p->props) + p->influx);
}


static double fluxfun_props(double s, int cell, const Opm::IncompPropertiesInterface* props)
{
    const double* visc = props->viscosity();
    double sat[2] = { s, 1.0 - s };
    double mob[2];
    props->relperm(1, sat, &cell, mob, 0);
    mob[0] /= visc[0];
    mob[1] /= visc[1];
    return mob[0]/(mob[0] + mob[1]);
}


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
