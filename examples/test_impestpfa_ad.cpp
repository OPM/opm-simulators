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

#define HACK_INCOMPRESSIBLE_GRAVITY 0

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/ImpesTPFAAD.hpp>
#include <opm/autodiff/BlackoilPropsAd.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/core/simulator/initState.hpp>

#include <opm/core/wells.h>
// #include <opm/core/WellsManager.hpp>

#include <algorithm>




int
main(int argc, char* argv[])
{
    const Opm::parameter::ParameterGroup param(argc, argv, false);
    const Opm::GridManager               gm(5, 5);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   oldprops(param, 2, nc);
    const Opm::BlackoilPropsAd           props(oldprops);

    typedef AutoDiff::ForwardBlock<double>      ADB;

    Wells* wells = create_wells(2, 2, 5);
    const double inj_frac[] = { 1.0, 0.0 };
    const double prod_frac[] = { 0.0, 0.0 };
    const int num_inj = 3;
    const int inj_cells[num_inj] = { 0, 1, 2 };
    const int num_prod = 2;
    const int prod_cells[num_prod] = { 20, 21 };
    const double WI[3] = { 1e-12, 1e-12, 1e-12 };
    bool ok = add_well(INJECTOR, 0.0, num_inj, inj_frac, inj_cells, WI, "Inj", wells);
    ok = ok && add_well(PRODUCER, 0.0, num_prod, prod_frac, prod_cells, WI, "Prod", wells);
    ok = ok && append_well_controls(BHP, 500.0*Opm::unit::barsa, 0, 0, wells);
    // ok = ok && append_well_controls(BHP, 200.0*Opm::unit::barsa, 0, 1, wells);
    double oildistr[2] = { 0.0, 1.0 };
    ok = ok && append_well_controls(SURFACE_RATE, 1e-3, oildistr, 1, wells);
    if (!ok) {
        THROW("Something went wrong with well init.");
    }
    set_current_control(0, 0, wells);
    set_current_control(1, 0, wells);

    double grav[] = { /*1.0*/ 0.0, 0.0 };
    Opm::DerivedGeology geo(*g, props, grav);
    Opm::LinearSolverFactory linsolver(param);
    Opm::ImpesTPFAAD ps(*g, props, geo, *wells, linsolver);

    Opm::BlackoilState state;
    initStateBasic(*g, oldprops, param, 0.0, state);
    initBlackoilSurfvol(*g, oldprops, state);
    Opm::WellState well_state;
    well_state.init(wells, state);

    ps.solve(1.0, state, well_state);

    std::cout << "Cell pressure:" << std::endl;
    std::copy(state.pressure().begin(), state.pressure().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Face flux:" << std::endl;
    std::copy(state.faceflux().begin(), state.faceflux().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Well bhp pressure:" << std::endl;
    std::copy(well_state.bhp().begin(), well_state.bhp().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    return 0;
}
