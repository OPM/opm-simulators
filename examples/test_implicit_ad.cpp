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

#include <config.h>

#include <opm/autodiff/FullyImplicitBlackoilSolver.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilPropsAd.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/initState.hpp>

#include <algorithm>

int
main(int argc, char* argv[])
{
    const Opm::parameter::ParameterGroup param(argc, argv, false);
    const Opm::GridManager               gm(20, 1);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   props0(param, 2, nc);
    const Opm::BlackoilPropsAd           props(props0);

    typedef Opm::FullyImplicitBlackoilSolver<Opm::DerivedGeology> BOSolver;

    double grav[] = { 1.0, 0.0 };
    Opm::DerivedGeology geo(*g, props, grav);

    BOSolver solver(*g, props, geo);

    Opm::BlackoilState state;
    initStateBasic(*g, props0, param, 0.0, state);
    initBlackoilSurfvol(*g, props0, state);

    solver.step(1.0, state);

#if 0
    Opm::WellState well_state;
    well_state.init(wells, state);

    ps.solve(1.0, state, well_state);

    std::copy(state.pressure().begin(), state.pressure().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::copy(well_state.bhp().begin(), well_state.bhp().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
#endif

    return 0;
}
