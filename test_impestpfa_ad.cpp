#include "ImpesTPFAAD.hpp"

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>

int
main(int argc, char* argv[])
{
    const Opm::parameter::ParameterGroup param(argc, argv);
    const Opm::GridManager               gm(3, 3);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   props(param, 2, nc);

    typedef Opm::ImpesTPFAAD<Opm::BlackoilPropertiesInterface> PSolver;
    PSolver ps(*g, props);
}
