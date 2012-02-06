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


struct ParametersSRes
{
    double c;
    double s0;
    double dtpv;    /* dt/pv(i)                  */
    double influx;  /* sum_j min(v_ij, 0)*f(s_j) */
    double outflux; /* sum_j max(v_ij, 0)        */
    int cell;
    const Opm::IncompPropertiesInterface* props;
    const PolymerData* polydata;
};

struct ParametersCRes
{
    double s0;
    double c0;
    double cmax0;
    double dtpv;    /* dt/pv(i)                  */
    double influx;  /* sum_j min(v_ij, 0)*f(s_j) */
    double influx_polymer;
    double outflux; /* sum_j max(v_ij, 0)        */
    double porosity;
    PolymerSolverData* psdata;
    int cell;
    NonlinearSolverCtrl* ctrl;
    double s;
    const PolymerData* polydata;
};


static struct ParametersSRes get_parameters_s(struct PolymerSolverData *d, int cell);
static struct ParametersCRes get_parameters_c(struct PolymerSolverData *d, int cell,
					      NonlinearSolverCtrl* ctrl);
static double residual_s(double s, void *data);
static double residual_c(double c, void *data);
static double fluxfun_props(double s,
			    double c,
			    int cell,
			    const Opm::IncompPropertiesInterface* props, 
			    const PolymerData* polydata);
static double compute_mc(double c,
			 const Opm::IncompPropertiesInterface* props, 
			 const PolymerData* polydata);

void
destroy_solverdata(struct PolymerSolverData *d)
{
    if (d!=NULL)
    {
        free(d->fractionalflow);
        free(d->mc);
    }
    free(d);
}

struct PolymerSolverData *
init_solverdata(struct UnstructuredGrid *grid,
		const Opm::IncompPropertiesInterface* props,
		const PolymerData* polydata,
		const double *darcyflux,
                const double *porevolume,
                const double *porosity,
		const double *source,
                const double dt,
		double *saturation,
		double *concentration,
		double *cmax)
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
        d->porosity   = porosity;
        d->source     = source;
        d->dt         = dt;

        d->saturation     = saturation;
	d->concentration  = concentration;
	d->cmax           = cmax;
        d->fractionalflow = (double*) malloc(grid->number_of_cells *
					     sizeof *d->fractionalflow);
        d->mc             = (double*) malloc(grid->number_of_cells *
					     sizeof *d->mc);
        if (d->fractionalflow == NULL || d->mc == NULL)
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
    struct ParametersCRes   prm = get_parameters_c(d, cell, ctrl);

    d->concentration[cell] = find_zero(residual_c, &prm, ctrl);
    d->cmax[cell] = std::max(d->cmax[cell], d->concentration[cell]);
    d->saturation[cell] = prm.s;
    d->fractionalflow[cell] = fluxfun_props(d->saturation[cell], d->concentration[cell], cell, d->props, d->polydata);
    d->mc[cell] = compute_mc(d->concentration[cell], d->props, d->polydata);
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
    struct ParametersSRes *p = (struct ParametersSRes*) data;
    double c = p->c;
    return s - p->s0 +  p->dtpv*(p->outflux*fluxfun_props(s, c, p->cell, p->props, p->polydata) + p->influx);
}

static double
residual_c(double c, void *data)
{
    struct ParametersCRes *p = (struct ParametersCRes*) data;
    int cell = p->cell;
    struct ParametersSRes prm_s = get_parameters_s(p->psdata, cell);
    prm_s.c = c;
    double s = find_zero(residual_s, &prm_s, p->ctrl);
    p->s = s;
    double ff = fluxfun_props(s, c, p->cell, prm_s.props, p->polydata);
    double mc = compute_mc(c, prm_s.props, p->polydata);
    double dps = p->polydata->dps;
    double s0 = p->s0;
    double c0 = p->c0;
    double rhor = p->polydata->rhor;
    double porosity = p->porosity;
    double cmax0 = p->cmax0;
    double ads0 = p->polydata->adsorbtion(std::max(c0, cmax0));
    double ads = p->polydata->adsorbtion(std::max(c, cmax0));
    return  (s - dps)*c - (s0 - dps)*c0
	+ rhor*((1.0 - porosity)/porosity)*(ads - ads0)
	+ p->dtpv*(p->outflux*ff*mc + p->influx_polymer);
}

static struct ParametersSRes
get_parameters_s(struct PolymerSolverData *d, int cell)
{
    int i;
    struct UnstructuredGrid *g  = d->grid;
    struct ParametersSRes        p;
    double flux;
    int f, other;

    p.c       = d->concentration[cell];
    p.s0      = d->saturation[cell];
    p.influx  = d->source[cell] >  0 ? -d->source[cell] : 0.0;
    p.outflux = d->source[cell] <= 0 ? -d->source[cell] : 0.0;
    p.dtpv    = d->dt/d->porevolume[cell];
    p.cell    = cell;
    p.props   = d->props;

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


static struct ParametersCRes
get_parameters_c(struct PolymerSolverData *d, int cell, NonlinearSolverCtrl* ctrl)
{
    int i;
    struct UnstructuredGrid *g  = d->grid;
    double flux;
    int f, other;

    ParametersCRes p;
    p.c0      = d->concentration[cell];
    p.cmax0   = d->cmax[cell];
    p.s0      = d->saturation[cell];
    p.dtpv    = d->dt/d->porevolume[cell];
    p.influx  = d->source[cell] >  0 ? -d->source[cell] : 0.0;
    p.influx_polymer  = d->source[cell] >  0 ? -d->source[cell] : 0.0; // TODO. Wrong if nonzero, mult by mc.
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
		p.influx_polymer += flux*d->fractionalflow[other]*d->mc[other];
            }
            else {
                p.outflux += flux;
            }
        }
    }
    p.outflux = d->source[cell] <= 0 ? -d->source[cell] : 0.0;
    p.porosity = d->porosity[cell];
    p.psdata = d;
    p.cell = cell;
    p.ctrl = ctrl;
    p.s = -1e100;
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

static double compute_mc(double c,
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
    return c/(inv_mu_w_eff*mu_p_eff);
}


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
