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

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/core/simulator/initState.hpp>

#include <opm/core/wells.h>
// #include <opm/core/WellsManager.hpp>

#include <algorithm>

namespace {
    template <class Geology, class Vector>
    class DerivedGeology {
    public:
        typedef Vector V;

        DerivedGeology(const UnstructuredGrid& grid,
                       const Geology&          geo)
            : pvol_ (grid.number_of_cells)
            , trans_(grid.number_of_faces)
        {
            // Pore volume
            const typename Vector::Index nc = grid.number_of_cells;
            std::transform(grid.cell_volumes, grid.cell_volumes + nc,
                           geo.porosity(), pvol_.data(),
                           std::multiplies<double>());

            // Transmissibility
            Vector htrans(grid.cell_facepos[nc]);
            UnstructuredGrid* ug = const_cast<UnstructuredGrid*>(& grid);
            tpfa_htrans_compute(ug, geo.permeability(), htrans.data());
            tpfa_trans_compute (ug, htrans.data()     , trans_.data());
        }

        const Vector& poreVolume()       const { return pvol_ ; }
        const Vector& transmissibility() const { return trans_; }

    private:
        Vector pvol_ ;
        Vector trans_;
    };
}

int
main(int argc, char* argv[])
{
    const Opm::parameter::ParameterGroup param(argc, argv, false);
    const Opm::GridManager               gm(3, 3);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   props(param, 2, nc);

    typedef AutoDiff::ForwardBlock<double>      ADB;
    typedef Opm::BlackoilPropertiesInterface    Geology;
    typedef DerivedGeology<Geology, ADB::V>     GeoProps;
    typedef Opm::BlackoilPropertiesInterface    BOFluid;
    typedef Opm::ImpesTPFAAD<BOFluid, GeoProps> PSolver;

    Wells* wells = create_wells(2, 2, 2);
    const double inj_frac[] = { 1.0, 0.0 };
    const double prod_frac[] = { 0.0, 0.0 };
    const int inj_cell = 0;
    const int prod_cell = g->number_of_cells - 1;
    const double WI = 1e-8;
    bool ok = add_well(INJECTOR, 0.0, 1, inj_frac, &inj_cell, &WI, "Inj", wells);
    ok = ok && add_well(PRODUCER, 0.0, 1, prod_frac, &prod_cell, &WI, "Prod", wells);
    ok = ok && append_well_controls(BHP, 500.0*Opm::unit::barsa, 0, 0, wells);
    ok = ok && append_well_controls(BHP, 200.0*Opm::unit::barsa, 0, 1, wells);
    if (!ok) {
        THROW("Something went wrong with well init.");
    }

    GeoProps geo(*g, props);
    Opm::LinearSolverFactory linsolver(param);
    PSolver  ps (*g, props, geo, *wells, linsolver);

    Opm::BlackoilState state;
    initStateBasic(*g, props, param, 0.0, state);
    initBlackoilSurfvol(*g, props, state);
    Opm::WellState well_state;
    well_state.init(wells, state);

    ps.solve(1.0, state, well_state);

    std::copy(state.pressure().begin(), state.pressure().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::copy(well_state.bhp().begin(), well_state.bhp().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    return 0;
}
