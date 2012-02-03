/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */
/* Copyright 2012 (c) SINTEF */


#include <opm/polymer/polymermodel.hpp>
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/nlsolvers.h>
#include <opm/core/fluid/IncompPropertiesInterface.hpp>

#include <cstdlib>
#include <cstdio>
#include <cmath>


/* Parameters used in solution of single-cell boundary-value problem */


struct Parameters
{
    double c;
    double s0;
    // double c0;
    double dtpv;    /* dt/pv(i)                  */
    double influx;  /* sum_j min(v_ij, 0)*f(s_j) */
    // double influx_polymer;
    double outflux; /* sum_j max(v_ij, 0)        */
    // double dps;
    // double rhor;
    // double phi;
    int cell;
    const Opm::IncompPropertiesInterface* props;
    const PolymerData* polydata;
};


static struct Parameters get_parameters(struct PolymerSolverData *d, int cell);
static double residual_s(double s, void *data);
static double residual_c(double c, void *data);
static double fluxfun_props(double s,
			    double c,
			    int cell,
			    const Opm::IncompPropertiesInterface* props, 
			    const PolymerData* polydata);

void
destroy_solverdata(struct PolymerSolverData *d)
{
    if (d!=NULL)
    {
        free(d->fractionalflow);
    }
    free(d);
}

struct PolymerSolverData *
init_solverdata(struct UnstructuredGrid *grid,
		const Opm::IncompPropertiesInterface* props,
		const PolymerData* polydata,
		const double *darcyflux,
                const double *porevolume,
		const double *source,
                const double dt,
		double *saturation,
		double *concentration)
{
    int i;
    struct PolymerSolverData *d = (struct PolymerSolverData*) malloc(sizeof *d);

    if(d!=NULL)
    {
        d->grid       = grid;
	d->props      = props;
	d->polydata   = polydata;
        d->darcyflux  = darcyflux;
        d->porevolume = porevolume;
        d->source     = source;
        d->dt         = dt;

        d->saturation     = saturation;
	d->concentration  = concentration;
        d->fractionalflow = (double*) malloc(grid->number_of_cells *
					     sizeof *d->fractionalflow);
        if (d->fractionalflow == NULL)
        {
            destroy_solverdata(d);
            d = NULL;
        }
        for(i=0; i<grid->number_of_cells; ++i)
        {
            d->fractionalflow[i] = 0.0;
        }
    }
    return d;
}

/* Solver for single-cell bvp calls root-finder in nlsolvers.c */
void polymer_solvecell(void *data, struct NonlinearSolverCtrl *ctrl, int cell)
{
    struct PolymerSolverData   *d   = (struct PolymerSolverData*) data;
    struct Parameters   prm = get_parameters(d, cell);

    d->saturation[cell] = find_zero(residual_s, &prm, ctrl);
    // double ff1 = fluxfun_props(d->saturation[cell], cell, d->props);
    // double ff2 = fluxfun(d->saturation[cell], -999);
    // printf("New = %f   old = %f\n", ff1, ff2);
    d->fractionalflow[cell] = fluxfun_props(d->saturation[cell], d->concentration[cell], cell, d->props, d->polydata);
}


/* ====================== Internals =================================*/


/* Residual function r(s) for a single-cell bvp */
/*
 *     r(s) = s - s0 + dt/pv*(influx - outflux*f(s) )
 */
/* influx is water influx, outflux is total outflux */
static double
residual_s(double s, void *data)
{
    struct Parameters *p = (struct Parameters*) data;
    double c = p->c;
    return s - p->s0 +  p->dtpv*(p->outflux*fluxfun_props(s, c, p->cell, p->props, p->polydata) + p->influx);
}

static double
residual_c(double c, void *data)
{
    (void) c;
    (void) data;
    // struct Parameters *p = (struct Parameters*) data;
    //    return s - p->s0 +  p->dtpv*(p->outflux*fluxfun_props(s, p->cell, p->props) + p->influx);
    return 0.0;
}

static struct Parameters
get_parameters(struct PolymerSolverData *d, int cell)
{
    int i;
    struct UnstructuredGrid *g  = d->grid;
    struct Parameters        p;
    double flux;
    int f, other;

    p.c       = d->concentration[cell];
    p.s0      = d->saturation[cell];
    p.influx  = d->source[cell] >  0 ? -d->source[cell] : 0.0;
    p.outflux = d->source[cell] <= 0 ? -d->source[cell] : 0.0;
    p.dtpv    = d->dt/d->porevolume[cell];
    p.cell    = cell;
    p.props   = d->props;

    d->saturation[cell] = 0;
    for (i=g->cell_facepos[cell]; i<g->cell_facepos[cell+1]; ++i) {
        f = g->cell_faces[i];

        /* Compute cell flux*/
        if (cell == g->face_cells[2*f]) {
            flux  = d->darcyflux[f];
            other = g->face_cells[2*f+1];
        }
        else {
            flux  =-d->darcyflux[f];
            other = g->face_cells[2*f];
        }

        if (other != -1) {
            if (flux < 0.0) {
                p.influx  += flux*d->fractionalflow[other];
            }
            else {
                p.outflux += flux;
            }
        }
    }
    p.polydata = d->polydata;
    return p;
}

static double fluxfun_props(double s, double c, int cell,
			    const Opm::IncompPropertiesInterface* props,
			    const PolymerData* pd)
{
    const double* visc = props->viscosity();
    double c_max_limit = pd->c_max_limit;
    double cbar = c/c_max_limit;
    double mu_w = visc[0];
    double mu_m = pd->viscMult(c)*mu_w;
    double mu_p = pd->viscMult(pd->c_max_limit)*mu_w;
    double omega = pd->omega;
    double mu_m_omega = std::pow(mu_m, omega);
    double mu_w_e   = mu_m_omega*std::pow(mu_w, 1.0 - omega);
    double mu_p_eff = mu_m_omega*std::pow(mu_p, 1.0 - omega);
    double inv_mu_w_eff = (1.0 - cbar)/mu_w_e + cbar/mu_p_eff;
    double inv_visc_eff[2] = { inv_mu_w_eff, 1.0/visc[1] };
    double sat[2] = { s, 1.0 - s };
    double mob[2];
    props->relperm(1, sat, &cell, mob, NULL);
    mob[0] *= inv_visc_eff[0];
    mob[1] *= inv_visc_eff[1];
    return mob[0]/(mob[0] + mob[1]);
}


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
