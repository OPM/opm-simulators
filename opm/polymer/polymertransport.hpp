/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */
/* Copyright 2012 (c) SINTEF */

#ifndef POLYMERTRANSPORT_HPP_INCLUDED
#define POLYMERTRANSPORT_HPP_INCLUDED

namespace Opm
{
    class IncompPropertiesInterface;
}

struct UnstructuredGrid;
struct PolymerData;

void polymertransport(
    const double *porevolume,
    const double *porosity,
    const double *source,
    double dt,
    struct UnstructuredGrid *grid,
    const Opm::IncompPropertiesInterface* props,
    const PolymerData* polydata,
    const double *darcyflux,
    double *saturation,
    double *concentration,
    double *cmax);


#endif /* POLYMERTRANSPORT_HPP_INCLUDED */
