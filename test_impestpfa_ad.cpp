/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2013 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include "ImpesTPFAAD.hpp"

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>

int
main(int argc, char* argv[])
{
    const Opm::parameter::ParameterGroup param(argc, argv, false);
    const Opm::GridManager               gm(3, 3);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   props(param, 2, nc);

    typedef Opm::ImpesTPFAAD<Opm::BlackoilPropertiesInterface> PSolver;
    PSolver ps(*g, props);

    return 0;
}
