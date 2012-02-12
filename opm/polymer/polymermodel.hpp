/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */
/* Copyright 2012 (c) SINTEF */

#ifndef POLYMER_HPP_INCLUDED
#define POLYMER_HPP_INCLUDED

#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/polymer/TransportModelPolymer.hpp> // For PolymerData.

struct UnstructuredGrid;
namespace Opm
{
    class IncompPropertiesInterface;
}


struct PolymerSolverData  {
    struct UnstructuredGrid *grid;
    const Opm::IncompPropertiesInterface* props;
    const Opm::PolymerData* polydata;
    const double            *darcyflux;   /* one flux per face  in cdata::grid*/
    const double            *porevolume;  /* one volume per cell */
    const double            *porosity;
    const double            *source;      /* one source per cell */
    double                   dt;
    double                   inflow_c;
    double                  *saturation;      /* one per cell */
    double                  *concentration;   /* one per cell */
    double                  *cmax;
    double                  *fractionalflow;  /* one per cell */
    double                  *mc;
};


void
polymer_solvecell(void *data, int cell);

void
destroy_solverdata(struct PolymerSolverData *d);

struct PolymerSolverData *
init_solverdata(struct UnstructuredGrid *grid,
		const Opm::IncompPropertiesInterface* props,
		const Opm::PolymerData* polydata,
		const double *darcyflux,
                const double *porevolume,
		const double *porosity,
		const double *source,
                const double dt,
		const double inflow_c,
		double *saturation,
		double *concentration,
		double *cmax);


#endif /* POLYMER_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
