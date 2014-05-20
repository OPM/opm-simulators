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
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>

#include <opm/core/grid.h>
#include <opm/core/wells.h>

#include <opm/core/grid/GridManager.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/initState.hpp>

#include <memory>
#include <algorithm>
#include <iostream>


namespace {
    std::shared_ptr<Wells>
    createWellConfig()
    {
        std::shared_ptr<Wells> wells(create_wells(2, 2, 2),
                                       destroy_wells);

        const double inj_frac[] = { 1.0, 0.0 };
        const double prod_frac[] = { 0.0, 0.0 };
        const int num_inj = 1;
        const int inj_cells[num_inj] = { 0 };
        const int num_prod = 1;
        const int prod_cells[num_prod] = { 19 };
        const double WI[3] = { 1e-12, 1e-12, 1e-12 };
        bool ok = add_well(INJECTOR, 0.0, num_inj, inj_frac, inj_cells, WI, "Inj", wells.get());
        ok = ok && add_well(PRODUCER, 0.0, num_prod, prod_frac, prod_cells, WI, "Prod", wells.get());
        ok = ok && append_well_controls(BHP, 500.0*Opm::unit::barsa, 0, 0, wells.get());
        // ok = ok && append_well_controls(BHP, 200.0*Opm::unit::barsa, 0, 1, wells);
        double oildistr[2] = { 0.0, 1.0 };
        ok = ok && append_well_controls(SURFACE_RATE, 1e-3, oildistr, 1, wells.get());
        if (!ok) {
            OPM_THROW(std::runtime_error, "Something went wrong with well init.");
        }
        set_current_control(0, 0, wells.get());
        set_current_control(1, 0, wells.get());

        return wells;
    }

    template <class Ostream, typename T, class A>
    Ostream&
    operator<<(Ostream& os, const std::vector<T,A>& v)
    {
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));

        return os;
    }
}

int
main(int argc, char* argv[])
try
{
    const Opm::parameter::ParameterGroup param(argc, argv, false);
    const Opm::GridManager               gm(20, 1);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   props0(param, 2, nc);
    const Opm::BlackoilPropsAd           props(props0);

    std::shared_ptr<Wells> wells = createWellConfig();

    double grav[] = { 0.0, 0.0 };
    Opm::DerivedGeology geo(*g, props, grav);

    Opm::LinearSolverFactory linsolver(param);
    Opm::NewtonIterationBlackoilSimple fis_solver(linsolver);

    Opm::FullyImplicitBlackoilSolver<UnstructuredGrid> solver(param, *g, props, geo, 0, *wells, fis_solver);

    Opm::BlackoilState state;
    initStateBasic(*g, props0, param, 0.0, state);
    initBlackoilSurfvol(*g, props0, state);

    Opm::WellStateFullyImplicitBlackoil well_state;
    well_state.init(wells.get(), state);

    solver.step(1.0, state, well_state);

    std::cout << state.pressure() << '\n'
              << well_state.bhp() << '\n';

    return 0;
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

